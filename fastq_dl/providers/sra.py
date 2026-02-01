import glob
import logging
import os
from pathlib import Path
from typing import Union

from pysradb import SRAweb

from fastq_dl.constants import PE_R1_SUFFIX, PE_R2_SUFFIX, SE_SUFFIX, SRA_FAILED
from fastq_dl.utils import execute


def get_sra_metadata(query: str) -> list:
    """Fetch metadata from SRA.

    Args:
        query (str): The accession to search for.

    Returns:
        list: Records associated with the accession.
    """
    #
    db = SRAweb()
    df = db.search_sra(
        query, detailed=True, sample_attribute=True, expand_sample_attributes=True
    )
    if df is None:
        return [False, []]
    return [True, df.to_dict(orient="records")]


def sra_download(
    accession: str,
    outdir: str,
    cpus: int = 1,
    max_attempts: int = 10,
    force: bool = False,
    ignore_md5: bool = False,
    sleep: int = 10,
    sra_lite: bool = False,
) -> Union[dict, str]:
    """Download FASTQs from SRA using fasterq-dump.

    Args:
        accession (str): The accession to download associated FASTQs.
        outdir (str): Directory to write FASTQs to.
        cpus (int, optional): Number of CPUs to use. Defaults to 1.
        max_attempts (int, optional): Maximum number of download attempts. Defaults to 10.
        force (bool, optional): Force overwrite of existing files
        ignore_md5 (bool, optional): Ignore MD5 checksums for downloaded files
        sleep (int, default = 10): Amount of seconds to sleep in between attempts
        sra_lite (bool, optional): If True, prefer SRA Lite downloads

    Returns:
        Union[dict, str]: A dictionary of the FASTQs and their paired status, or SRA_FAILED on error.
            The dict contains:
            - r1 (str): Path to R1 FASTQ file
            - r2 (str): Path to R2 FASTQ file (empty string if single-end)
            - single_end (bool): True if single-end, False if paired-end
            - orphan (str | None): Path to orphan reads file if present (paired-end with unpaired reads)
    """
    outdir = Path(outdir)
    fastqs = {"r1": "", "r2": "", "single_end": True, "orphan": None}
    se = outdir / f"{accession}{SE_SUFFIX}"
    pe1 = outdir / f"{accession}{PE_R1_SUFFIX}"
    pe2 = outdir / f"{accession}{PE_R2_SUFFIX}"

    # remove existing files if force is selected.
    if force:
        for f in [se, pe1, pe2]:
            if f.exists():
                f.unlink()
                logging.warning(f"Overwriting existing file: {f}")

    if not se.exists() and not (pe1.exists() and pe2.exists()):
        outdir.mkdir(parents=True, exist_ok=True)

        # Use argument lists instead of string commands to prevent command injection
        vdb_config_cmd = [
            "vdb-config",
            "--simplified-quality-scores",
            "yes" if sra_lite else "no",
        ]
        if sra_lite:
            logging.debug("Setting preference to SRA Lite")
        else:
            logging.debug("Setting preference to SRA Normalized")

        execute(vdb_config_cmd)

        prefetch_cmd = [
            "prefetch",
            accession,
            "--max-size",
            "10T",
            "-o",
            f"{accession}.sra",
            "-f",
            "yes" if force else "no",
            "--verify",
            "no" if ignore_md5 else "yes",
        ]

        outcome = execute(
            prefetch_cmd,
            max_attempts=max_attempts,
            directory=str(outdir),
            is_sra=True,
            sleep=sleep,
        )

        if outcome == SRA_FAILED:
            return outcome

        fasterq_dump_cmd = [
            "fasterq-dump",
            accession,
            "--split-3",
            "--mem",
            "1G",
            "--threads",
            str(cpus),
        ]
        if force:
            fasterq_dump_cmd.append("-f")
        # no need to check MD5 of downloaded fastq as it tests checksums as it reads
        # ref: https://github.com/ncbi/sra-tools/issues/285#issuecomment-586365769

        outcome = execute(
            fasterq_dump_cmd,
            max_attempts=max_attempts,
            directory=str(outdir),
            is_sra=True,
            sleep=sleep,
        )

        # Get list of download fastq files
        fastq_files = ' '.join([os.path.basename(f) for f in glob.glob(f"{str(outdir)}/{accession}*.fastq")])

        if outcome == SRA_FAILED:
            return outcome
        else:
            # pigz with glob pattern needs shell, but accession is validated
            execute(
                f"pigz --force -p {cpus} -n {fastq_files}",
                directory=str(outdir),
            )
            (outdir / f"{accession}.sra").unlink(missing_ok=True)
            logging.info(f"Downloaded FASTQs for {accession}")
    else:
        if se.exists():
            logging.debug(f"Skipping re-download of existing file: {se}")
        elif pe1.exists() and pe2.exists():
            logging.debug(f"Skipping re-download of existing file: {pe1}")
            logging.debug(f"Skipping re-download of existing file: {pe2}")

    if pe2.exists():
        # Paired end
        fastqs["r1"] = str(pe1)
        fastqs["r2"] = str(pe2)
        fastqs["single_end"] = False
        if se.exists():
            # Orphan reads file exists (unpaired reads from paired-end data)
            fastqs["orphan"] = str(se)
    else:
        fastqs["r1"] = str(se)
        fastqs["single_end"] = True

    return fastqs
