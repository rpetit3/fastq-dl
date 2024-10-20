import logging
import sys
import time

from fastq_dl.constants import ENA, SRA
from fastq_dl.providers.ena import get_ena_metadata
from fastq_dl.providers.sra import get_sra_metadata


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
