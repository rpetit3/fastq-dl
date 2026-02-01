#! /usr/bin/env python3
import logging
import sys
from pathlib import Path

import rich
import rich_click as click
from rich.logging import RichHandler

import fastq_dl
from fastq_dl.constants import (
    ENA,
    ENA_FAILED,
    MERGED_R1_SUFFIX,
    MERGED_R2_SUFFIX,
    RUN_INFO_SUFFIX,
    RUN_MERGERS_SUFFIX,
    SE_SUFFIX,
    SRA,
    SRA_FAILED,
)
from fastq_dl.exceptions import DownloadError, FastqDLError, ProviderError, ValidationError
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
                "--protocol"
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
    "--protocol",
    default="ftp",
    show_default=True,
    help="Protocol to use for ENA downloads.",
    type=click.Choice(['ftp', 'https'], case_sensitive=False)
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
    protocol
):
    """Download FASTQ files from ENA or SRA."""
    # Setup logs only if no handlers are already configured (allows library usage)
    root_logger = logging.getLogger()
    if not root_logger.handlers:
        logging.basicConfig(
            format="%(asctime)s:%(name)s:%(levelname)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
            handlers=[
                RichHandler(rich_tracebacks=True, console=rich.console.Console(stderr=True))
            ],
        )

    root_logger.setLevel(
        logging.ERROR if silent else logging.DEBUG if verbose else logging.INFO
    )

    try:
        _run_download(
            accession=accession,
            provider=provider,
            group_by_experiment=group_by_experiment,
            group_by_sample=group_by_sample,
            outdir=outdir,
            prefix=prefix,
            max_attempts=max_attempts,
            sleep=sleep,
            force=force,
            ignore_md5=ignore_md5,
            sra_lite=sra_lite,
            only_provider=only_provider,
            only_download_metadata=only_download_metadata,
            cpus=cpus,
            protocol=protocol,
        )
    except ValidationError as e:
        logging.error(str(e))
        sys.exit(1)
    except ProviderError as e:
        logging.error(f"Provider error ({e.provider}): {e}")
        sys.exit(1)
    except DownloadError as e:
        logging.error(f"Download error: {e}")
        sys.exit(1)
    except FastqDLError as e:
        logging.error(f"Error: {e}")
        sys.exit(1)


def _run_download(
    accession: str,
    provider: str,
    group_by_experiment: bool,
    group_by_sample: bool,
    outdir: str,
    prefix: str,
    max_attempts: int,
    sleep: int,
    force: bool,
    ignore_md5: bool,
    sra_lite: bool,
    only_provider: bool,
    only_download_metadata: bool,
    cpus: int,
    protocol: str = "ftp",
) -> None:
    """Internal function that performs the actual download logic.

    This function contains the core download logic and may raise FastqDLError
    subclasses which are caught by the CLI wrapper.

    Raises:
        ValidationError: If the accession is invalid.
        ProviderError: If there's an error querying ENA/SRA.
        DownloadError: If there's an error downloading files.
    """
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
    outdir = Path(outdir).resolve() if outdir else Path.cwd()

    if only_download_metadata:
        outdir.mkdir(parents=True, exist_ok=True)
        run_info_path = outdir / f"{prefix}{RUN_INFO_SUFFIX}"
        logging.info(f"Writing metadata to {run_info_path}")
        write_tsv(ena_data, str(run_info_path))
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
            logging.info(f"Working on run {run_acc}...")

            # Determine primary provider and whether fallback is allowed
            if provider.lower() == "ena" and data_from == ENA:
                primary = "ena"
                no_fallback = only_provider
            else:
                primary = "sra"
                no_fallback = only_provider or data_from == SRA

            fastqs, error = _download_with_fallback(
                run_info=run_info,
                run_acc=run_acc,
                outdir=str(outdir),
                primary_provider=primary,
                only_provider=no_fallback,
                max_attempts=max_attempts,
                force=force,
                ignore_md5=ignore_md5,
                sleep=sleep,
                cpus=cpus,
                sra_lite=sra_lite,
                protocol=protocol,
            )

            if error:
                ena_data[i]["error"] = error

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
                # Check if we have valid paired-end data (both r1 and r2 with equal counts)
                if vals["r1"] and vals["r2"] and len(vals["r1"]) == len(vals["r2"]):
                    logging.info(f"Merging paired end runs to {name}...")
                    merge_runs(vals["r1"], str(outdir / f"{name}{MERGED_R1_SUFFIX}"))
                    merge_runs(vals["r2"], str(outdir / f"{name}{MERGED_R2_SUFFIX}"))
                else:
                    # Single-end merge (covers no r2, empty r2, or mismatched counts)
                    logging.info("Merging single end runs to experiment...")
                    merge_runs(vals["r1"], str(outdir / f"{name}{SE_SUFFIX}"))
            run_mergers_path = outdir / f"{prefix}{RUN_MERGERS_SUFFIX}"
            logging.info(f"Writing merged run info to {run_mergers_path}")
            write_tsv(runs, str(run_mergers_path))
        run_info_path = outdir / f"{prefix}{RUN_INFO_SUFFIX}"
        logging.info(f"Writing metadata to {run_info_path}")
        write_tsv(ena_data, str(run_info_path))


def _download_with_fallback(
    run_info: dict,
    run_acc: str,
    outdir: str,
    primary_provider: str,
    only_provider: bool,
    max_attempts: int,
    force: bool,
    ignore_md5: bool,
    sleep: int,
    cpus: int,
    sra_lite: bool,
    protocol: str = "ftp",
) -> tuple:
    """Download FASTQs with fallback to alternate provider.

    Args:
        run_info: Dictionary containing run metadata (needed for ENA download)
        run_acc: Run accession string
        outdir: Output directory
        primary_provider: 'ena' or 'sra' - which provider to try first
        only_provider: If True, do not fallback to alternate provider
        max_attempts: Maximum download attempts
        force: Overwrite existing files
        ignore_md5: Skip MD5 verification
        sleep: Sleep time between retries
        cpus: Number of CPUs for SRA download
        sra_lite: Use SRA Lite format
        protocol: Protocol for ENA downloads (ftp or https)

    Returns:
        tuple: (fastqs dict or None, error string or None)
    """
    if primary_provider.lower() == "ena":
        primary_download = lambda: ena_download(
            run_info,
            outdir,
            max_attempts=max_attempts,
            force=force,
            ignore_md5=ignore_md5,
            sleep=sleep,
            protocol=protocol,
        )
        primary_failed = ENA_FAILED
        fallback_download = lambda: sra_download(
            run_acc,
            outdir,
            cpus=cpus,
            max_attempts=max_attempts,
            sleep=sleep,
            sra_lite=sra_lite,
        )
        fallback_failed = SRA_FAILED
        primary_name, fallback_name = "ENA", "SRA"
    else:
        primary_download = lambda: sra_download(
            run_acc,
            outdir,
            cpus=cpus,
            max_attempts=max_attempts,
            sleep=sleep,
            sra_lite=sra_lite,
        )
        primary_failed = SRA_FAILED
        fallback_download = lambda: ena_download(
            run_info,
            outdir,
            max_attempts=max_attempts,
            force=force,
            ignore_md5=ignore_md5,
            sleep=sleep,
            protocol=protocol,
        )
        fallback_failed = ENA_FAILED
        primary_name, fallback_name = "SRA", "ENA"

    # Try primary provider
    fastqs = primary_download()

    if fastqs == primary_failed:
        if only_provider:
            logging.error(f"{run_acc} not found on {primary_name}")
            return None, primary_failed

        # Fallback to alternate provider
        logging.info(f"{run_acc} not found on {primary_name}, retrying from {fallback_name}")
        fastqs = fallback_download()

        if fastqs == fallback_failed:
            logging.error(f"{run_acc} not found on {fallback_name}")
            return None, f"{primary_failed}&{fallback_failed}"

    return fastqs, None


def main():
    if len(sys.argv) == 1:
        fastqdl(["--help"])
    else:
        fastqdl()


if __name__ == "__main__":
    main()
