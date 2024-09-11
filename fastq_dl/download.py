import logging
import sys
import time
from pathlib import Path

import requests
from executor import ExternalCommand, ExternalCommandFailed
from pysradb import SRAweb

from fastq_dl import ENA, ENA_FAILED, ENA_URL, SRA, SRA_FAILED
from fastq_dl.utils import md5sum


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


def get_run_info(
    accession: str,
    query: str,
    provider: str,
    only_provider: bool,
    max_attempts: int = 10,
    sleep: int = 10,
) -> tuple:
    """Retrieve a list of samples available from ENA.

    The first attempt will be against ENA, and if that fails, SRA will be queried. This should
    capture those samples not yet synced between ENA and SRA.

    Args:
        accession (str): The accession to search for.
        query (str): A formatted query for ENA searches.
        provider (str): Limit queries only to the specified provider (requires only_provider be true)
        only_provider (bool): If true, limit queries to the specified provider
        max_attempts (int, optional): Maximum number of download attempts
        sleep (int): Minimum amount of time to sleep before retry

    Returns:
        tuple: Records associated with the accession.
    """
    # Attempt for when "--only-provider" used, others to allow multiple attempts on fallback
    attempt = 1
    sra_attempt = 1
    ena_attempt = 1
    while True:
        if only_provider:
            logging.debug(
                f"Querying for metadata (Attempt {attempt} of {max_attempts})"
            )
            logging.debug(f"--only-provider supplied, limiting queries to {provider}")
            if provider.lower() == "ena":
                success, ena_data = get_ena_metadata(query)
                if success:
                    return ENA, ena_data
                elif attempt >= max_attempts:
                    logging.error("There was an issue querying ENA, exiting...")
                    logging.error(f"STATUS: {ena_data[0]}")
                    logging.error(f"TEXT: {ena_data[1]}")
                    sys.exit(1)
            else:
                success, sra_data = get_sra_metadata(accession)
                if success:
                    return SRA, sra_data
                elif attempt >= max_attempts:
                    logging.error("There was an issue querying SRA, exiting...")
                    sys.exit(1)
            attempt += 1
            logging.warning(
                f"Querying {provider.lower()} was unsuccessful, retrying after ({sleep} seconds)"
            )
            time.sleep(sleep)
        else:
            if ena_attempt < max_attempts:
                logging.debug(
                    f"Querying ENA for metadata (Attempt {ena_attempt} of {max_attempts})"
                )
            else:
                logging.debug(
                    f"Querying SRA for metadata (Attempt {sra_attempt} of {max_attempts})"
                )

            success, ena_data = get_ena_metadata(query)
            if success:
                return ENA, ena_data
            elif ena_attempt >= max_attempts:
                if ena_attempt == max_attempts:
                    ena_attempt += 1
                    logging.debug("Failed to get metadata from ENA. Trying SRA...")
                success, sra_data = get_sra_metadata(accession)
                if success:
                    return SRA, sra_data
                elif sra_attempt >= max_attempts:
                    logging.error("There was an issue querying ENA and SRA, exiting...")
                    logging.error(f"STATUS: {ena_data[0]}")
                    logging.error(f"TEXT: {ena_data[1]}")
                    sys.exit(1)
                else:
                    sra_attempt += 1
                    logging.warning(
                        f"Querying SRA was unsuccessful, retrying after ({sleep} seconds)"
                    )
                    time.sleep(sleep)
            else:
                ena_attempt += 1
                logging.warning(
                    f"Querying ENA was unsuccessful, retrying after ({sleep} seconds)"
                )
                time.sleep(sleep)


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


def sra_download(
    accession: str,
    outdir: str,
    cpus: int = 1,
    max_attempts: int = 10,
    force: bool = False,
    ignore_md5: bool = False,
    sleep: int = 10,
    sra_lite: bool = False,
) -> dict:
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
        dict: A dictionary of the FASTQs and their paired status.
    """
    outdir = Path(outdir)
    fastqs = {"r1": "", "r2": "", "single_end": True}
    se = outdir / f"{accession}.fastq.gz"
    pe1 = outdir / f"{accession}_1.fastq.gz"
    pe2 = outdir / f"{accession}_2.fastq.gz"

    # remove existing files if force is selected.
    if force:
        for f in [se, pe1, pe2]:
            if f.exists():
                f.unlink()
                logging.warning(f"Overwriting existing file: {f}")

    if not se.exists() and not (pe1.exists() and pe2.exists()):
        outdir.mkdir(parents=True, exist_ok=True)

        vdb_config_cmd = "vdb-config --simplified-quality-scores "
        if sra_lite:
            # Prefer SRA Lite
            logging.debug("Setting preference to SRA Lite")
            vdb_config_cmd += "yes"
        else:
            # Prefer SRA Normalized
            logging.debug("Setting preference to SRA Normalized")
            vdb_config_cmd += "no"

        execute(vdb_config_cmd)

        prefetch_cmd = f"prefetch {accession} --max-size 10T -o {accession}.sra"
        prefetch_cmd += " -f yes" if force else " -f no"
        prefetch_cmd += " --verify no" if ignore_md5 else " --verify yes"

        outcome = execute(
            prefetch_cmd,
            max_attempts=max_attempts,
            directory=str(outdir),
            is_sra=True,
            sleep=sleep,
        )

        if outcome == SRA_FAILED:
            return outcome

        fasterq_dump_cmd = (
            f"fasterq-dump {accession} --split-3 --mem 1G --threads {cpus}"
        )
        fasterq_dump_cmd += " -f" if force else ""
        # no need to check MD5 of downloaded fastq as it tests checksums as it reads
        # ref: https://github.com/ncbi/sra-tools/issues/285#issuecomment-586365769

        outcome = execute(
            fasterq_dump_cmd,
            max_attempts=max_attempts,
            directory=str(outdir),
            is_sra=True,
            sleep=sleep,
        )

        if outcome == SRA_FAILED:
            return outcome
        else:
            execute(
                f"pigz --force -p {cpus} -n {accession}*.fastq", directory=str(outdir)
            )
            (outdir / f"{accession}.sra").unlink()
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
        if se.exists():
            fastqs["single_end"] = str(se)
        else:
            fastqs["single_end"] = False
    else:
        fastqs["r1"] = str(se)

    return fastqs


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
        accession (dict): Dictionary of run info to download associated FASTQs.
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
