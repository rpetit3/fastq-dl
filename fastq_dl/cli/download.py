#! /usr/bin/env python3
import logging
import sys
from pathlib import Path

import rich
import rich_click as click
from rich.logging import RichHandler

import fastq_dl
from fastq_dl.constants import ENA, ENA_FAILED, SRA, SRA_FAILED
from fastq_dl.providers.ena import ena_download
from fastq_dl.providers.generic import get_run_info
from fastq_dl.providers.sra import sra_download
from fastq_dl.utils import merge_runs, validate_query, write_tsv

click.rich_click.USE_RICH_MARKUP = True
click.rich_click.OPTION_GROUPS = {
    "fastq-dl": [
        {
            "name": "Required Options",
            "options": ["--accession"],
        },
        {
            "name": "Download Options",
            "options": [
                "--provider",
                "--group-by-experiment",
                "--group-by-sample",
                "--max-attempts",
                "--sra-lite",
                "--only-provider",
                "--only-download-metadata",
                "--ignore",
            ],
        },
        {
            "name": "Additional Options",
            "options": [
                "--outdir",
                "--prefix",
                "--cpus",
                "--force",
                "--silent",
                "--sleep",
                "--version",
                "--verbose",
                "--help",
            ],
        },
    ]
}


@click.command()
@click.version_option(fastq_dl.__version__, "--version", "-V")
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
    "--group-by-experiment",
    is_flag=True,
    help="Group Runs by experiment accession.",
)
@click.option(
    "--group-by-sample",
    is_flag=True,
    help="Group Runs by sample accession.",
)
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
    "--sleep",
    "-s",
    default=10,
    show_default=True,
    help="Minimum amount of time to sleep between retries (API query and download)",
)
@click.option(
    "-F",
    "--force",
    is_flag=True,
    help="Overwrite existing files.",
)
@click.option(
    "-I",
    "--ignore",
    "ignore_md5",
    is_flag=True,
    help="Ignore MD5 checksums for downloaded files.",
)
@click.option(
    "--sra-lite",
    is_flag=True,
    help="Set preference to SRA Lite",
)
@click.option(
    "--only-provider",
    is_flag=True,
    help="Only attempt download from specified provider.",
)
@click.option(
    "--only-download-metadata",
    is_flag=True,
    help="Skip FASTQ downloads, and retrieve only the metadata.",
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
    sleep,
    force,
    ignore_md5,
    sra_lite,
    only_provider,
    only_download_metadata,
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
    data_from, ena_data = get_run_info(
        accession,
        query,
        provider,
        only_provider,
        max_attempts=max_attempts,
        sleep=sleep,
    )

    logging.info(f"Query: {accession}")
    logging.info(f"Archive: {provider}")
    if only_download_metadata:
        logging.info(f"Total Runs Found: {len(ena_data)}")
        logging.debug("--only-download-metadata used, skipping FASTQ downloads")
    else:
        logging.info(f"Total Runs To Download: {len(ena_data)}")
    downloaded = {}
    runs = {} if group_by_experiment or group_by_sample else None
    outdir = Path.cwd() if outdir == "./" else f"{outdir}"

    if only_download_metadata:
        Path(outdir).mkdir(parents=True, exist_ok=True)
        logging.info(f"Writing metadata to {outdir}/{prefix}-run-info.tsv")
        write_tsv(ena_data, f"{outdir}/{prefix}-run-info.tsv")
    else:
        for i, run_info in enumerate(ena_data):
            run_acc = run_info["run_accession"]
            if run_acc not in downloaded:
                downloaded[run_acc] = True
            else:
                logging.warning(
                    f"Duplicate run {run_acc} found, skipping re-download..."
                )
                continue
            logging.info(f"\tWorking on run {run_acc}...")
            if provider.lower() == "ena" and data_from == ENA:
                fastqs = ena_download(
                    run_info,
                    outdir,
                    max_attempts=max_attempts,
                    force=force,
                    ignore_md5=ignore_md5,
                    sleep=sleep,
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
                            sleep=sleep,
                            sra_lite=sra_lite,
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
                    sleep=sleep,
                    sra_lite=sra_lite,
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
                            force=force,
                            ignore_md5=ignore_md5,
                            sleep=sleep,
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
            logging.info(
                f"Writing merged run info to {outdir}/{prefix}-run-mergers.tsv"
            )
            write_tsv(runs, f"{outdir}/{prefix}-run-mergers.tsv")
        logging.info(f"Writing metadata to {outdir}/{prefix}-run-info.tsv")
        write_tsv(ena_data, f"{outdir}/{prefix}-run-info.tsv")


def main():
    if len(sys.argv) == 1:
        fastqdl(["--help"])
    else:
        fastqdl()


if __name__ == "__main__":
    main()
