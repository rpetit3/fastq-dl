import csv
import hashlib
import logging
import re
import shlex
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Union

from fastq_dl.constants import ENA_FAILED, SRA_DOWNLOAD_FAILED, SRA_FAILED
from fastq_dl.exceptions import ValidationError

PathLike = Union[str, Path]


def execute(
    cmd: Union[str, list],
    directory: str = str(Path.cwd()),
    capture_stdout: bool = False,
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

            if capture_stdout:
                return result.stdout
            else:
                return result.returncode

        except subprocess.CalledProcessError as e:
            logging.error(f'"{cmd}" return exit code {e.returncode}')
            logging.debug(f"STDOUT: {e.stdout}")
            logging.debug(f"STDERR: {e.stderr}")

            # sracha's Error::NotFound produces stderr "not found: {accession}"
            # (see rnabioco/sracha-rs error.rs). Unlike sra-tools (exit code 3),
            # sracha exits 1 for all errors, so we parse stderr to distinguish
            # "not found" (immediate ENA fallback) from transient failures (retry).
            if is_sra and e.stderr and "not found:" in e.stderr.lower():
                error_msg = e.stderr.split("\n")[0] if e.stderr else "Unknown error"
                logging.error(error_msg)
                return SRA_FAILED

            if attempt < max_attempts:
                logging.error(f"Retry execution ({attempt} of {max_attempts})")
                time.sleep(sleep)
            else:
                if is_sra:
                    return SRA_DOWNLOAD_FAILED
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


def _all_fieldnames(rows: Iterable[Mapping[str, Any]]) -> List[str]:
    """
    Collect the union of all keys across a collection of dictionaries.

    This function preserves first-seen order of keys rather than sorting them,
    which helps maintain stable and human-readable column ordering based on
    input data appearance.

    Args:
        rows (Iterable[Mapping[str, Any]]):
            An iterable of dictionary-like objects representing rows.

    Returns:
        List[str]:
            A list of unique field names (keys) appearing across all rows,
            in order of first occurrence.
    """
    fieldnames: List[str] = []
    seen = set()
    for r in rows:
        for k in r.keys():
            if k not in seen:
                seen.add(k)
                fieldnames.append(k)
    return fieldnames


def write_tsv(
    data: Union[Dict[str, Any], List[Dict[str, Any]]], output: str, na_value: str = ""
) -> None:
    """
    Write heterogeneous dictionary data to a TSV file.

    This function supports two input formats:

    1. A dictionary of the form:
         {
             accession: {"r1": [...], "r2": [...]},
             ...
         }
       When `output` ends with "-run-mergers.tsv", the file will contain
       columns: accession, r1, r2. Lists are joined with ";" and missing
       values are replaced with `na_value`.

    2. A list of dictionaries:
         [
             {"col1": val1, "col2": val2, ...},
             ...
         ]
       In this case, the union of all keys across rows is computed and used
       as the header. Missing keys in any row are filled with `na_value`.

    Args:
        data (Union[Dict[str, Any], List[Dict[str, Any]]]):
            Either:
                - A dictionary mapping identifiers to nested dictionaries
                  (used for run-mergers output), or
                - A list of row dictionaries for general TSV writing.
        output (str):
            Path to the output TSV file.
        na_value (str, optional):
            Value used to fill missing fields. Defaults to "" (empty string).

    Returns:
        None

    Raises:
        TypeError:
            If the provided data structure does not match expected formats.
    """
    with open(output, "w", newline="") as fh:
        if output.endswith("-run-mergers.tsv"):
            writer = csv.DictWriter(
                fh,
                fieldnames=["accession", "r1", "r2"],
                delimiter="\t",
                extrasaction="ignore",
            )
            writer.writeheader()
            # here data is expected to be dict: accession -> {"r1":[...], "r2":[...]}
            for accession, vals in data.items():  # type: ignore[union-attr]
                r1_list = (vals or {}).get("r1", [])
                r2_list = (vals or {}).get("r2", [])
                writer.writerow(
                    {
                        "accession": accession,
                        "r1": ";".join(r1_list) if r1_list else na_value,
                        "r2": ";".join(r2_list) if r2_list else na_value,
                    }
                )
        else:
            # here data is expected to be a list of dict rows
            rows: List[Dict[str, Any]] = list(data)  # type: ignore[arg-type]
            fieldnames = _all_fieldnames(rows)

            writer = csv.DictWriter(
                fh,
                fieldnames=fieldnames,
                delimiter="\t",
                extrasaction="ignore",  # optional safety; unioned fieldnames should already include everything
            )
            writer.writeheader()

            for row in rows:
                # fill missing keys with what was provided by na_value
                out_row = {k: row.get(k, na_value) for k in fieldnames}
                writer.writerow(out_row)


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
