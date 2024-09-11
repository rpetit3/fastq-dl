import csv
import hashlib
import logging
import re
import shutil
import sys
from pathlib import Path
from typing import Optional

from fastq_dl import PathLike


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
