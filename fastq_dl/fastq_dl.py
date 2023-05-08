#! /usr/bin/env python3
import csv
import hashlib
import logging
import re
import sys
import time
from pathlib import Path

import requests
import rich
import rich_click as click
from executor import ExternalCommand, ExternalCommandFailed
from pysradb import SRAweb
from rich.logging import RichHandler

click.rich_click.USE_RICH_MARKUP = True
click.rich_click.OPTION_GROUPS = {
    "fastq-dl": [
        {"name": "Required Options", "options": ["--accession"]},
        {
            "name": "Additional Options",
            "options": [
                "--provider",
                "--group-by-experiment",
                "--group-by-sample",
                "--outdir",
                "--prefix",
                "--cpus",
                "--max-attempts",
                "--only-provider",
                "--silent",
                "--version",
                "--verbose",
                "--help",
            ],
        },
    ]
}

PROGRAM = "fastq-dl"
VERSION = "2.0.1"
ENA_FAILED = "ENA_NOT_FOUND"
SRA_FAILED = "SRA_NOT_FOUND"
SRA = "SRA"
ENA = "ENA"
ENA_URL = "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&format=tsv"


def execute(
    cmd: str,
    directory: str = str(Path.cwd()),
    capture_stdout: bool = False,
    stdout_file: str = None,
    stderr_file: str = None,
    max_attempts: int = 1,
    is_sra: bool = False,
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

    Raises:
        error: An unexpected error occurred.

    Returns:
        str: Exit code, accepted error message, or STDOUT of command.
    """
    attempt = 0
    while attempt < max_attempts:
        attempt += 1
        try:
            command = ExternalCommand(
                cmd,
                directory=directory,
                capture=True,
                capture_stderr=True,
                stdout_file=stdout_file,
                stderr_file=stderr_file,
            )

            command.start()
            logging.debug(command.decoded_stdout)
            logging.debug(command.decoded_stderr)

            if capture_stdout:
                return command.decoded_stdout
            else:
                return command.returncode
        except ExternalCommandFailed as error:
            logging.error(f'"{cmd}" return exit code {command.returncode}')

            if is_sra and command.returncode == 3:
                # The FASTQ isn't on SRA for some reason, try to download from ENA
                error_msg = command.decoded_stderr.split("\n")[0]
                logging.error(error_msg)
                return SRA_FAILED

            if attempt < max_attempts:
                logging.error(f"Retry execution ({attempt} of {max_attempts})")
                time.sleep(10)
            else:
                raise error


def sra_download(
    accession: str, outdir: str, cpus: int = 1, max_attempts: int = 10
) -> dict:
    """Download FASTQs from SRA using fasterq-dump.

    Args:
        accession (str): The accession to download associated FASTQs.
        outdir (str): Directory to write FASTQs to.
        cpus (int, optional): Number of CPUs to use. Defaults to 1.
        max_attempts (int, optional): Maximum number of download attempts. Defaults to 10.

    Returns:
        dict: A dictionary of the FASTQs and their paired status.
    """
    fastqs = {"r1": "", "r2": "", "single_end": True}
    se = f"{outdir}/{accession}.fastq.gz"
    pe = f"{outdir}/{accession}_2.fastq.gz"

    if not Path(se).exists() and not Path(pe).exists():
        Path(outdir).mkdir(parents=True, exist_ok=True)
        outcome = execute(
            f"fasterq-dump {accession} --split-files --threads {cpus}",
            max_attempts=max_attempts,
            directory=outdir,
            is_sra=True,
        )
        if outcome == SRA_FAILED:
            return outcome
        else:
            execute(f"pigz --force -p {cpus} -n {accession}*.fastq", directory=outdir)

    if Path(f"{outdir}/{accession}_2.fastq.gz").exists():
        # Paired end
        fastqs["r1"] = f"{outdir}/{accession}_1.fastq.gz"
        fastqs["r2"] = f"{outdir}/{accession}_2.fastq.gz"
        fastqs["single_end"] = False
    else:
        fastqs["r1"] = f"{outdir}/{accession}.fastq.gz"

    return fastqs


def ena_download(run: str, outdir: str, max_attempts: int = 10) -> dict:
    """Download FASTQs from ENA FTP using wget.

    Args:
        accession (str): The accession to download associated FASTQs.
        outdir (str): Directory to write FASTQs to.
        max_attempts (int, optional): Maximum number of download attempts. Defaults to 10.

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
            )

            if is_r2:
                fastqs["r2"] = fastq
                fastqs["single_end"] = False
            else:
                fastqs["r1"] = fastq

    return fastqs


def md5sum(fastq: str) -> str:
    """Calculate the MD5 checksum of a file.

    Source: https://stackoverflow.com/a/3431838/5299417

    Args:
        fastq (str): Input FASTQ to calculate MD5 checksum for.

    Returns:
        str: Calculated MD5 checksum.
    """
    MB = 1_048_576
    BUFFER_SIZE = 10 * MB
    if Path(fastq).exists():
        hash_md5 = hashlib.md5()
        with open(fastq, "rb") as fp:
            for chunk in iter(lambda: fp.read(BUFFER_SIZE), b""):
                hash_md5.update(chunk)

        return hash_md5.hexdigest()
    else:
        return None


def download_ena_fastq(
    ftp: str,
    outdir: str,
    md5: str,
    max_attempts: int = 10,
) -> str:
    """Download FASTQs from ENA using FTP.

    Args:
        ftp (str): The FTP address of the FASTQ file.
        outdir (str): Directory to download the FASTQ to.
        md5 (str): Expected MD5 checksum of the FASTQ.
        max_attempts (int, optional): Maximum number of download attempts. Defaults to 10.

    Returns:
        str: Path to the downloaded FASTQ.
    """
    success = False
    attempt = 0
    fastq = f"{outdir}/{Path(ftp).name}"

    if not Path(fastq).exists():
        Path(outdir).mkdir(parents=True, exist_ok=True)

        while not success:
            logging.info(f"\t\t{Path(ftp).name} FTP download attempt {attempt + 1}")
            execute(f"wget --quiet -O {fastq} ftp://{ftp}", max_attempts=max_attempts)

            fastq_md5 = md5sum(fastq)
            if fastq_md5 != md5:
                logging.debug(f"MD5s, Observed: {fastq_md5}, Expected: {md5}")
                attempt += 1
                if Path(fastq).exists():
                    Path(fastq).unlink()
                if attempt > max_attempts:
                    logging.error(
                        f"Download failed after {max_attempts} attempts. "
                        "Please try again later or manually from SRA/ENA."
                    )
                    sys.exit(1)

                # Hiccup? Wait a bit before trying again.
                time.sleep(10)
            else:
                success = True

    return fastq


def merge_runs(runs: list, output: str) -> None:
    """Merge runs from an experiment or sample.

    Args:
        runs (list): A list of FASTQs to merge.
        output (str): The final merged FASTQ.
    """
    if len(runs) > 1:
        run_fqs = " ".join(runs)
        execute(f"cat {run_fqs} > {output}")
        for p in runs:
            Path(p).unlink()
    else:
        Path(runs[0]).rename(output)


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
        return [True, data]
    else:
        return [False, [r.status_code, r.text]]


def get_run_info(accession: str, query: str) -> tuple:
    """Retrieve a list of samples available from ENA.

    The first attempt will be against ENA, and if that fails, SRA will be queried. This should
    capture those samples not yet synced between ENA and SRA.

    Args:
        accession (str): The accession to search for.
        query (str): A formatted query for ENA searches.

    Returns:
        tuple: Records associated with the accession.
    """

    logging.debug("Querying ENA for metadata...")
    success, ena_data = get_ena_metadata(query)
    if success:
        return ENA, ena_data
    else:
        logging.debug("Failed to get metadata from ENA. Trying SRA...")
        success, sra_data = get_sra_metadata(accession)
        if not success:
            logging.error("There was an issue querying ENA and SRA, exiting...")
            logging.error(f"STATUS: {ena_data[0]}")
            logging.error(f"TEXT: {ena_data[1]}")
            sys.exit(1)
        else:
            return SRA, sra_data


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


@click.command()
@click.version_option(VERSION, "--version", "-V")
@click.option(
    "-a",
    "--accession",
    required=True,
    help="ENA/SRA accession to query. (Study, Sample, Experiment, Run accession)",
)
@click.option(
    "--provider",
    default="ena",
    show_default=True,
    help="Specify which provider (ENA or SRA) to use.",
    type=click.Choice(
        ["ena", "sra"],
        case_sensitive=False,
    ),
)
@click.option(
    "--group-by-experiment", is_flag=True, help="Group Runs by experiment accession."
)
@click.option("--group-by-sample", is_flag=True, help="Group Runs by sample accession.")
@click.option(
    "--outdir",
    "-o",
    default="./",
    show_default=True,
    help="Directory to output downloads to.",
)
@click.option(
    "--prefix",
    default="fastq",
    show_default=True,
    help="Prefix to use for naming log files.",
)
@click.option(
    "--max-attempts",
    "-m",
    default=10,
    show_default=True,
    help="Maximum number of download attempts.",
)
@click.option(
    "-F",
    "--only-provider",
    is_flag=True,
    help="Only attempt download from specified provider.",
)
@click.option(
    "--cpus",
    default=1,
    show_default=True,
    help="Total cpus used for downloading from SRA.",
)
@click.option("--silent", is_flag=True, help="Only critical errors will be printed.")
@click.option("--verbose", "-v", is_flag=True, help="Print debug related text.")
@click.help_option("--help", "-h")
def fastqdl(
    accession,
    provider,
    group_by_experiment,
    group_by_sample,
    outdir,
    prefix,
    max_attempts,
    only_provider,
    cpus,
    silent,
    verbose,
):
    """Download FASTQ files from ENA or SRA."""
    # Setup logs
    logging.basicConfig(
        format="%(asctime)s:%(name)s:%(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            RichHandler(rich_tracebacks=True, console=rich.console.Console(stderr=True))
        ],
    )

    logging.getLogger().setLevel(
        logging.ERROR if silent else logging.DEBUG if verbose else logging.INFO
    )
    # Start Download Process
    query = validate_query(accession)
    data_from, ena_data = get_run_info(accession, query)

    logging.info(f"Query: {accession}")
    logging.info(f"Archive: {provider}")
    logging.info(f"Total Runs To Download: {len(ena_data)}")
    downloaded = {}
    runs = {} if group_by_experiment or group_by_sample else None
    outdir = Path.cwd() if outdir == "./" else f"{outdir}"

    for i, run_info in enumerate(ena_data):
        run_acc = run_info["run_accession"]
        if run_acc not in downloaded:
            downloaded[run_acc] = True
        else:
            logging.warning(f"Duplicate run {run_acc} found, skipping re-download...")
            continue
        logging.info(f"\tWorking on run {run_acc}...")
        fastqs = None
        if provider == "ena" and data_from == ENA:
            fastqs = ena_download(
                run_info,
                outdir,
                max_attempts=max_attempts,
            )

            if fastqs == ENA_FAILED:
                if only_provider:
                    logging.error(f"\tNo fastqs found in ENA for {run_acc}")
                    ena_data[i]["error"] = ENA_FAILED
                    fastqs = None
                else:
                    # Retry download from SRA
                    logging.info(f"\t{run_acc} not found on ENA, retrying from SRA")

                    fastqs = sra_download(
                        run_acc,
                        outdir,
                        cpus=cpus,
                        max_attempts=max_attempts,
                    )
                    if fastqs == SRA_FAILED:
                        logging.error(f"\t{run_acc} not found on SRA")
                        ena_data[i]["error"] = f"{ENA_FAILED}&{SRA_FAILED}"
                        fastqs = None

        else:
            fastqs = sra_download(
                run_acc,
                outdir,
                cpus=cpus,
                max_attempts=max_attempts,
            )
            if fastqs == SRA_FAILED:
                if only_provider or data_from == SRA:
                    logging.error(f"\t{run_acc} not found on SRA or ENA")
                    ena_data[i]["error"] = SRA_FAILED
                    fastqs = None
                else:
                    # Retry download from ENA
                    logging.info(f"\t{run_acc} not found on SRA, retrying from ENA")
                    fastqs = ena_download(
                        run_info,
                        outdir,
                        max_attempts=max_attempts,
                    )
                    if fastqs == ENA_FAILED:
                        logging.error(f"\tNo fastqs found in ENA for {run_acc}")
                        ena_data[i]["error"] = f"{SRA_FAILED}&{ENA_FAILED}"
                        fastqs = None

        # Add the download results
        if fastqs:
            if group_by_experiment or group_by_sample:
                name = run_info["sample_accession"]
                if group_by_experiment:
                    name = run_info["experiment_accession"]

                if name not in runs:
                    runs[name] = {"r1": [], "r2": []}

                if fastqs["single_end"]:
                    runs[name]["r1"].append(fastqs["r1"])
                else:
                    runs[name]["r1"].append(fastqs["r1"])
                    runs[name]["r2"].append(fastqs["r2"])

    # If applicable, merge runs
    if runs:
        for name, vals in runs.items():
            if len(vals["r1"]) and len(vals["r2"]):
                # Not all runs labeled as paired are actually paired.
                if len(vals["r1"]) == len(vals["r2"]):
                    logging.info(f"\tMerging paired end runs to {name}...")
                    merge_runs(vals["r1"], f"{outdir}/{name}_R1.fastq.gz")
                    merge_runs(vals["r2"], f"{outdir}/{name}_R2.fastq.gz")
                else:
                    logging.info("\tMerging single end runs to experiment...")
                    merge_runs(vals["r1"], f"{outdir}/{name}.fastq.gz")
            else:
                logging.info("\tMerging single end runs to experiment...")
                merge_runs(vals["r1"], f"{outdir}/{name}.fastq.gz")
        write_tsv(runs, f"{outdir}/{prefix}-run-mergers.tsv")
    write_tsv(ena_data, f"{outdir}/{prefix}-run-info.tsv")


def main():
    fastqdl()


if __name__ == "__main__":
    main()
