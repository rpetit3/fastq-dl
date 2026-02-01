import logging
from pathlib import Path
from typing import Literal, Union

import requests

from fastq_dl.constants import ENA_FAILED, ENA_URL, PE_R1_SUFFIX, PE_R2_SUFFIX
from fastq_dl.exceptions import DownloadError
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
    r = requests.get(url, headers=headers, timeout=30)
    if r.status_code == requests.codes.ok:
        data = []
        col_names = None
        for line in r.text.split("\n"):
            cols = line.split("\t")
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
    protocol: Literal['ftp', 'https'] = 'ftp'
) -> Union[dict, str]:
    """Download FASTQs from ENA FTP using wget.

    Args:
        run (dict): Dictionary of run info to download associated FASTQs.
        outdir (str): Directory to write FASTQs to.
        max_attempts (int, optional): Maximum number of download attempts. Defaults to 10.
        force: (bool, optional): Whether to overwrite existing files if the MD5's do not match
        ignore_md5 (bool, optional): Ignore MD5 checksums for downloaded files
        sleep (int): Minimum amount of time to sleep before retry
        protocol (str): Protocol for download (ftp or https)

    Returns:
        Union[dict, str]: A dictionary of the FASTQs and their paired status, or ENA_FAILED on error.
            The dict contains:
            - r1 (str): Path to R1 FASTQ file
            - r2 (str): Path to R2 FASTQ file (empty string if single-end)
            - single_end (bool): True if single-end, False if paired-end
            - orphan (str | None): Always None for ENA (orphan reads not supported)
    """
    fastqs = {"r1": "", "r2": "", "single_end": True, "orphan": None}
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
            if ftp[i].endswith(PE_R2_SUFFIX):
                # Example: ERR1143237_2.fastq.gz
                is_r2 = True
            elif ftp[i].endswith(PE_R1_SUFFIX):
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
                protocol=protocol
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
    protocol: Literal['ftp', 'https'] = 'ftp'
) -> Union[str, str]:
    """Download FASTQs from ENA using FTP or HTTPS.

    Args:
        ftp (str): The FTP address of the FASTQ file.
        outdir (str): Directory to download the FASTQ to.
        md5 (str): Expected MD5 checksum of the FASTQ.
        max_attempts (int, optional): Maximum number of download attempts. Defaults to 10.
        force: (bool, optional): Whether to overwrite existing files if the MD5's do not match
        ignore_md5 (bool, optional): Ignore MD5 checksums for downloaded files
        sleep (int): Minimum amount of time to sleep before retry
        protocol (str): Protocol for download (ftp or https)

    Returns:
        str: Path to the downloaded FASTQ, or ENA_FAILED on error.

    Raises:
        DownloadError: When download fails after max_attempts due to MD5 mismatch.
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
            logging.info(f"{fastq} {protocol.upper()} download attempt {attempt + 1}")
            outcome = execute(
                ["wget", "--quiet", "-O", str(fastq), f"{protocol}://{ftp}"],
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
                        raise DownloadError(
                            f"Download of {fastq} failed after {max_attempts} attempts "
                            "due to MD5 checksum mismatch. Please try again later or "
                            "download manually from SRA/ENA.",
                            accession=Path(ftp).stem,
                            provider="ENA",
                        )
                else:
                    logging.info(f"Successfully downloaded {fastq}")
                    success = True

    return str(fastq)
