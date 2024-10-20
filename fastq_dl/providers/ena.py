import logging
import sys
from pathlib import Path

import requests

from fastq_dl.constants import ENA_FAILED, ENA_URL
from fastq_dl.utils import execute, md5sum


def get_ena_metadata(query: str) -> list:
    """Fetch metadata from ENA.
    https://docs.google.com/document/d/1CwoY84MuZ3SdKYocqssumghBF88PWxUZ/edit#heading=h.ag0eqy2wfin5

    Args:
        query (str): The query to search for.

    Returns:
        list: Records associated with the accession.
    """
    url = f'{ENA_URL}&query="{query}"&fields=all'
    headers = {"Content-type": "application/x-www-form-urlencoded"}
    r = requests.get(url, headers=headers)
    if r.status_code == requests.codes.ok:
        data = []
        col_names = None
        for line in r.text.split("\n"):
            cols = line.rstrip().split("\t")
            if line:
                if col_names:
                    data.append(dict(zip(col_names, cols)))
                else:
                    col_names = cols
        if data:
            return [True, data]
        else:
            return [
                False,
                [r.status_code, "Query was successful, but received an empty response"],
            ]
    else:
        return [False, [r.status_code, r.text]]


def ena_download(
    run: dict,
    outdir: str,
    max_attempts: int = 10,
    force: bool = False,
    ignore_md5: bool = False,
    sleep: int = 10,
) -> dict:
    """Download FASTQs from ENA FTP using wget.

    Args:
        run (dict): Dictionary of run info to download associated FASTQs.
        outdir (str): Directory to write FASTQs to.
        max_attempts (int, optional): Maximum number of download attempts. Defaults to 10.
        force: (bool, optional): Whether to overwrite existing files if the MD5's do not match
        ignore_md5 (bool, optional): Ignore MD5 checksums for downloaded files
        sleep (int): Minimum amount of time to sleep before retry

    Returns:
        dict: A dictionary of the FASTQs and their paired status.
    """
    fastqs = {"r1": "", "r2": "", "single_end": True}
    ftp = run["fastq_ftp"]
    if not ftp:
        return ENA_FAILED

    ftp = ftp.split(";")
    md5 = run["fastq_md5"].split(";")
    for i in range(len(ftp)):
        is_r2 = False
        # If run is paired only include *_1.fastq and *_2.fastq, rarely a
        # run can have 3 files.
        # Example:ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR114/007/ERR1143237
        if run["library_layout"] == "PAIRED":
            if ftp[i].endswith("_2.fastq.gz"):
                # Example: ERR1143237_2.fastq.gz
                is_r2 = True
            elif ftp[i].endswith("_1.fastq.gz"):
                # Example: ERR1143237_1.fastq.gz
                pass
            else:
                # Example: ERR1143237.fastq.gz
                # Not a part of the paired end read, so skip this file. Or,
                # its the only fastq file, and its not a paired
                obs_fq = Path(ftp[i]).name
                exp_fq = f'{run["run_accession"]}.fastq.gz'
                if len(ftp) != 1 and obs_fq != exp_fq:
                    continue

        # Download Run
        if md5[i]:
            fastq = download_ena_fastq(
                ftp[i],
                outdir,
                md5[i],
                max_attempts=max_attempts,
                force=force,
                ignore_md5=ignore_md5,
                sleep=sleep,
            )
            if fastq == ENA_FAILED:
                return ENA_FAILED

            if is_r2:
                fastqs["r2"] = fastq
                fastqs["single_end"] = False
            else:
                fastqs["r1"] = fastq

    return fastqs


def download_ena_fastq(
    ftp: str,
    outdir: str,
    md5: str,
    max_attempts: int = 10,
    force: bool = False,
    ignore_md5: bool = False,
    sleep: int = 10,
) -> str:
    """Download FASTQs from ENA using FTP.

    Args:
        ftp (str): The FTP address of the FASTQ file.
        outdir (str): Directory to download the FASTQ to.
        md5 (str): Expected MD5 checksum of the FASTQ.
        max_attempts (int, optional): Maximum number of download attempts. Defaults to 10.
        force: (bool, optional): Whether to overwrite existing files if the MD5's do not match
        ignore_md5 (bool, optional): Ignore MD5 checksums for downloaded files
        sleep (int): Minimum amount of time to sleep before retry

    Returns:
        str: Path to the downloaded FASTQ.
    """
    success = False
    attempt = 0
    outdir = Path(outdir)
    fastq = outdir / Path(ftp).name
    download_fastq = True

    if fastq.exists() and force:
        logging.warning(f"Overwriting existing file: {fastq}")
        fastq.unlink()
    elif fastq.exists() and not force:
        if ignore_md5:
            logging.warning(f"Skipping re-download of existing file: {fastq}")
            download_fastq = False
        else:
            logging.debug(f"Checking the MD5 of the existing file {fastq}...")
            fastq_md5 = md5sum(fastq)
            if fastq_md5 == md5:
                logging.info(f"MD5s match, skipping re-download of {fastq}")
                download_fastq = False
            else:
                logging.warning(f"MD5s do not match, re-downloading {fastq}")
                fastq.unlink()

    if download_fastq:
        outdir.mkdir(parents=True, exist_ok=True)

        while not success:
            logging.info(f"\t\t{fastq} FTP download attempt {attempt + 1}")
            outcome = execute(
                f"wget --quiet -O {fastq} ftp://{ftp}",
                max_attempts=max_attempts,
                sleep=sleep,
            )
            if outcome == ENA_FAILED:
                return ENA_FAILED

            if ignore_md5:
                logging.debug(f"--ignore used, skipping MD5 check for {fastq}")
                success = True
            else:
                fastq_md5 = md5sum(fastq)
                if fastq_md5 != md5:
                    logging.warning(
                        f"MD5 checksums do not match, attempting re-download of {fastq}"
                    )
                    attempt += 1
                    if fastq.exists():
                        fastq.unlink()
                    if attempt > max_attempts:
                        logging.error(
                            f"Download failed after {max_attempts} attempts. "
                            "Please try again later or manually from SRA/ENA."
                        )
                        sys.exit(1)
                else:
                    logging.info(f"Successfully downloaded {fastq}")
                    success = True

    return str(fastq)
