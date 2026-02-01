import csv
import hashlib
import logging
import re
import shlex
import shutil
import subprocess
import time
from pathlib import Path
from typing import Optional, Union

from fastq_dl.constants import ENA_FAILED, RUN_MERGERS_SUFFIX, SRA_FAILED
from fastq_dl.exceptions import ValidationError

PathLike = Union[str, Path]


def execute(
    cmd: Union[str, list],
    directory: str = str(Path.cwd()),
    capture_stdout: bool = False,
    stdout_file: str = None,
    stderr_file: str = None,
    max_attempts: int = 1,
    is_sra: bool = False,
    sleep: int = 10,
) -> str:
    """Execute a command using subprocess.

    Args:
        cmd (Union[str, list]): A command to execute. Can be a string (will be split with shlex)
            or a list of arguments (preferred for security).
        directory (str, optional): Set the working directory for command. Defaults to str(Path.cwd()).
        capture_stdout (bool, optional): Capture and return the STDOUT of a command. Defaults to False.
        stdout_file (str, optional): File to write STDOUT to. Defaults to None.
        stderr_file (str, optional): File to write STDERR to. Defaults to None.
        max_attempts (int, optional): Maximum times to attempt command execution. Defaults to 1.
        is_sra (bool, optional): The command is from SRA. Defaults to False.
        sleep (int): Minimum amount of time to sleep before retry

    Returns:
        str: Exit code, accepted error message, or STDOUT of command.
    """
    # Convert string commands to list for subprocess
    if isinstance(cmd, str):
        cmd_list = shlex.split(cmd)
    else:
        cmd_list = cmd

    attempt = 0
    while attempt < max_attempts:
        attempt += 1
        logging.debug(f"Executing command: {cmd_list}")
        logging.debug(f"Working directory: {directory}")

        try:
            result = subprocess.run(
                cmd_list,
                cwd=directory,
                capture_output=True,
                text=True,
                check=True,
            )

            logging.debug(f"STDOUT: {result.stdout}")
            logging.debug(f"STDERR: {result.stderr}")

            # Write stdout to file if specified
            if stdout_file and result.stdout:
                with open(stdout_file, "w") as f:
                    f.write(result.stdout)

            # Write stderr to file if specified
            if stderr_file and result.stderr:
                with open(stderr_file, "w") as f:
                    f.write(result.stderr)

            if capture_stdout:
                return result.stdout
            else:
                return result.returncode

        except subprocess.CalledProcessError as e:
            logging.error(f'"{cmd}" return exit code {e.returncode}')
            logging.debug(f"STDOUT: {e.stdout}")
            logging.debug(f"STDERR: {e.stderr}")

            # Write stderr to file even on failure if specified
            if stderr_file and e.stderr:
                with open(stderr_file, "w") as f:
                    f.write(e.stderr)

            if is_sra and e.returncode == 3:
                # The FASTQ isn't on SRA for some reason, try to download from ENA
                error_msg = e.stderr.split("\n")[0] if e.stderr else "Unknown error"
                logging.error(error_msg)
                return SRA_FAILED

            if attempt < max_attempts:
                logging.error(f"Retry execution ({attempt} of {max_attempts})")
                time.sleep(sleep)
            else:
                if is_sra:
                    return SRA_FAILED
                else:
                    return ENA_FAILED


def md5sum(fastq: PathLike) -> Optional[str]:
    """Calculate the MD5 checksum of a file.

    Source: https://stackoverflow.com/a/3431838/5299417

    Args:
        fastq (str): Input FASTQ to calculate MD5 checksum for.

    Returns:
        str: Calculated MD5 checksum.
    """
    fastq = Path(fastq)
    megabyte = 1_048_576
    buffer_size = 10 * megabyte
    if fastq.exists():
        hash_md5 = hashlib.md5()
        with open(fastq, "rb") as fp:
            for chunk in iter(lambda: fp.read(buffer_size), b""):
                hash_md5.update(chunk)

        return hash_md5.hexdigest()
    else:
        return None


def merge_runs(runs: list, output: str) -> None:
    """Merge runs from an experiment or sample.

    Args:
        runs (list): A list of FASTQs to merge.
        output (str): The final merged FASTQ.

    Raises:
        FileNotFoundError: If any of the input files do not exist.
    """
    paths = [Path(r) for r in runs]

    # Validate all files exist before starting
    missing = [str(p) for p in paths if not p.exists()]
    if missing:
        raise FileNotFoundError(f"Missing files for merge: {missing}")

    if len(paths) > 1:
        # concatenate the files in runs into output
        with open(output, "wb") as wfd:
            for p in paths:
                with open(p, "rb") as fd:
                    shutil.copyfileobj(fd, wfd)
        # Only delete source files after successful merge completion
        for p in paths:
            p.unlink()
    else:
        paths[0].rename(output)


def write_tsv(data: Union[list, dict], output: str) -> None:
    """Write a TSV file.

    Args:
        data: Data to be written to TSV. Can be either:
            - list[dict]: List of row dictionaries (for run-info.tsv)
            - dict[str, dict]: Dictionary of accession -> {r1, r2} mappings (for run-mergers.tsv)
        output (str): File to write the TSV to.
    """
    with open(output, "w") as fh:
        if output.endswith(RUN_MERGERS_SUFFIX):
            writer = csv.DictWriter(
                fh, fieldnames=["accession", "r1", "r2"], delimiter="\t"
            )
            writer.writeheader()
            for accession, vals in data.items():
                writer.writerow(
                    {
                        "accession": accession,
                        "r1": ";".join(vals["r1"]),
                        "r2": ";".join(vals["r2"]),
                    }
                )
        else:
            # Handle empty data case
            if not data:
                return
            writer = csv.DictWriter(fh, fieldnames=data[0].keys(), delimiter="\t")
            writer.writeheader()
            for row in data:
                writer.writerow(row)


def validate_query(query: str) -> str:
    """
    Check that query is an accepted accession type and return the accession type. Current
    accepted types are:

        Projects - PRJEB, PRJNA, PRJDA
        Studies - ERP, DRP, SRP
        BioSamples - SAMD, SAME, SAMN
        Samples - ERS, DRS, SRS
        Experiments - ERX, DRX, SRX
        Runs - ERR, DRR, SRR

    Parameters:
        query (str): A string containing an accession.

    Returns:
        str: A string containing the query for ENA search.

    https://ena-docs.readthedocs.io/en/latest/submit/general-guide/accessions.html
    """
    if re.match(r"^PRJ[EDN][A-Z][0-9]+$|^[EDS]RP[0-9]{6,}$", query):
        # Is a project or study accession
        return f"(study_accession={query} OR secondary_study_accession={query})"
    elif re.match(r"^SAM[EDN][A-Z]?[0-9]+$|^[EDS]RS[0-9]{6,}$", query):
        # Is a sample or biosample accession
        return f"(sample_accession={query} OR secondary_sample_accession={query})"
    elif re.match(r"^[EDS]RX[0-9]{6,}$", query):
        # Is an experiment accession
        return f"experiment_accession={query}"
    elif re.match(r"^[EDS]RR[0-9]{6,}$", query):
        # Is a run accession
        return f"run_accession={query}"
    else:
        raise ValidationError(
            f"{query} is not a Study, Sample, Experiment, or Run accession. "
            "See https://ena-docs.readthedocs.io/en/latest/submit/general-guide/accessions.html for valid options"
        )
