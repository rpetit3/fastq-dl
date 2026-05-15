import logging
import time
from pathlib import Path
from typing import Union

from pysradb import SRAweb

from fastq_dl.constants import (
    PE_R1_SUFFIX,
    PE_R1_SUFFIX_UNCOMPRESSED,
    PE_R2_SUFFIX,
    PE_R2_SUFFIX_UNCOMPRESSED,
    SE_SUFFIX,
    SE_SUFFIX_UNCOMPRESSED,
    SRA_FAILED,
)
from fastq_dl.utils import execute


def get_sra_metadata(query: str, max_attempts: int = 3, sleep: int = 10) -> list:
    """Fetch metadata from SRA.

    Args:
        query (str): The accession to search for.
        max_attempts (int): Maximum number of query attempts. Defaults to 3.
        sleep (int): Seconds to wait between retry attempts. Defaults to 10.

    Returns:
        list: Records associated with the accession.
    """
    for attempt in range(1, max_attempts + 1):
        try:
            db = SRAweb()
            df = db.search_sra(
                query,
                detailed=True,
                sample_attribute=True,
                expand_sample_attributes=True,
            )
            if df is None:
                return [False, []]
            return [True, df.to_dict(orient="records")]
        except Exception as e:
            logging.warning(
                f"pysradb query failed (Attempt {attempt} of {max_attempts}): {e}"
            )
            if attempt < max_attempts:
                time.sleep(sleep)
    return [False, []]


def sra_download(
    accession: str,
    outdir: str,
    cpus: int = 1,
    max_attempts: int = 10,
    force: bool = False,
    no_strict: bool = False,
    sleep: int = 10,
    sra_lite: bool = False,
    compress: bool = True,
    gzip_level: int = 1,
) -> Union[dict, str]:
    """Download FASTQs from SRA using sracha.

    Args:
        accession: The accession to download associated FASTQs.
        outdir: Directory to write FASTQs to.
        cpus: Number of CPUs for download and compression. Defaults to 1.
        max_attempts: Maximum number of download attempts. Defaults to 10.
        force: Force overwrite of existing files.
        no_strict: Downgrade integrity failures to warnings.
        sleep: Seconds to sleep between retry attempts. Defaults to 10.
        sra_lite: If True, prefer SRA Lite downloads (simplified quality scores).
        compress: If True, gzip compress output. Defaults to True.
        gzip_level: Gzip compression level (1-9). Defaults to 1.

    Returns:
        A dictionary of the FASTQs and their paired status, or SRA_FAILED on error.
        The dict contains:
        - r1 (str): Path to R1 FASTQ file
        - r2 (str): Path to R2 FASTQ file (empty string if single-end)
        - single_end (bool): True if single-end, False if paired-end
        - orphan (str | None): Path to orphan reads file if present
    """
    outdir = Path(outdir)
    fastqs = {"r1": "", "r2": "", "single_end": True, "orphan": None}
    if compress:
        se = outdir / f"{accession}{SE_SUFFIX}"
        pe1 = outdir / f"{accession}{PE_R1_SUFFIX}"
        pe2 = outdir / f"{accession}{PE_R2_SUFFIX}"
    else:
        se = outdir / f"{accession}{SE_SUFFIX_UNCOMPRESSED}"
        pe1 = outdir / f"{accession}{PE_R1_SUFFIX_UNCOMPRESSED}"
        pe2 = outdir / f"{accession}{PE_R2_SUFFIX_UNCOMPRESSED}"

    if force:
        for f in [se, pe1, pe2]:
            if f.exists():
                f.unlink()
                logging.warning(f"Overwriting existing file: {f}")

    if not se.exists() and not (pe1.exists() and pe2.exists()):
        outdir.mkdir(parents=True, exist_ok=True)

        sracha_cmd = [
            "sracha",
            "get",
            accession,
            "-O",
            str(outdir),
            "--threads",
            str(cpus),
            "--connections",
            str(cpus),
            "--no-progress",
            "-y",
        ]

        if sra_lite:
            sracha_cmd.extend(["--format", "sralite"])
            logging.debug("Setting preference to SRA Lite")
        else:
            logging.debug("Setting preference to SRA Normalized")

        if force:
            sracha_cmd.append("--force")

        if no_strict:
            sracha_cmd.append("--no-strict")

        if not compress:
            sracha_cmd.append("--no-gzip")
        elif gzip_level != 1:
            sracha_cmd.extend(["--gzip-level", str(gzip_level)])

        outcome = execute(
            sracha_cmd,
            max_attempts=max_attempts,
            directory=str(outdir),
            is_sra=True,
            sleep=sleep,
        )

        if outcome == SRA_FAILED:
            return outcome

        logging.info(f"Downloaded FASTQs for {accession}")
    else:
        if se.exists():
            logging.info(f"Skipping re-download of existing file: {se}")
        elif pe1.exists() and pe2.exists():
            logging.info(f"Skipping re-download of existing file: {pe1}")
            logging.info(f"Skipping re-download of existing file: {pe2}")

    if pe2.exists():
        fastqs["r1"] = str(pe1)
        fastqs["r2"] = str(pe2)
        fastqs["single_end"] = False
        if se.exists():
            fastqs["orphan"] = str(se)
    else:
        fastqs["r1"] = str(se)
        fastqs["single_end"] = True

    return fastqs
