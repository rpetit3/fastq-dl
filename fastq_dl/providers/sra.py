import logging
from pathlib import Path

from pysradb import SRAweb

from fastq_dl.constants import SRA_FAILED
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
