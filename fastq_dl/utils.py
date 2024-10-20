import csv
import hashlib
import logging
import re
import shutil
import sys
import time
from pathlib import Path
from typing import Optional, Union

from executor import ExternalCommand, ExternalCommandFailed

from fastq_dl.constants import ENA_FAILED, SRA_FAILED

PathLike = Union[str, Path]


def execute(
    cmd: str,
    directory: str = str(Path.cwd()),
    capture_stdout: bool = False,
    stdout_file: str = None,
    stderr_file: str = None,
    max_attempts: int = 1,
    is_sra: bool = False,
    sleep: int = 10,
) -> str:
    """A simple wrapper around executor.

    Args:
        cmd (str): A command to execute.
        directory (str, optional): Set the working directory for command. Defaults to str(Path.cwd()).
        capture_stdout (bool, optional): Capture and return the STDOUT of a command. Defaults to False.
        stdout_file (str, optional): File to write STDOUT to. Defaults to None.
        stderr_file (str, optional): File to write STDERR to. Defaults to None.
        max_attempts (int, optional): Maximum times to attempt command execution. Defaults to 1.
        is_sra (bool, optional): The command is from SRA. Defaults to False.
        sleep (int): Minimum amount of time to sleep before retry

    Raises:
        error: An unexpected error occurred.

    Returns:
        str: Exit code, accepted error message, or STDOUT of command.
    """
    attempt = 0
    while attempt < max_attempts:
        attempt += 1
        command = ExternalCommand(
            cmd,
            directory=directory,
            capture=True,
            capture_stderr=True,
            stdout_file=stdout_file,
            stderr_file=stderr_file,
        )
        try:
            command.start()
            logging.debug(command.decoded_stdout)
            logging.debug(command.decoded_stderr)

            if capture_stdout:
                return command.decoded_stdout
            else:
                return command.returncode
        except ExternalCommandFailed:
            logging.error(f'"{cmd}" return exit code {command.returncode}')

            if is_sra and command.returncode == 3:
                # The FASTQ isn't on SRA for some reason, try to download from ENA
                error_msg = command.decoded_stderr.split("\n")[0]
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
    """
    if len(runs) > 1:
        # concatenate the files in runs into output
        with open(output, "wb") as wfd:
            for p in map(Path, runs):
                with open(p, "rb") as fd:
                    shutil.copyfileobj(fd, wfd)
                p.unlink()
    else:
        Path(runs[0]).rename(output)


def write_tsv(data: dict, output: str) -> None:
    """Write a TSV file.

    Args:
        data (dict): Data to be written to TSV.
        output (str): File to write the TSV to.
    """
    with open(output, "w") as fh:
        if output.endswith("-run-mergers.tsv"):
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
        logging.error(
            f"{query} is not a Study, Sample, Experiment, or Run accession. See https://ena-docs.readthedocs.io/en/latest/submit/general-guide/accessions.html for valid options"
        )
        sys.exit(1)
