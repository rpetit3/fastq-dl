import logging
import time

from fastq_dl.constants import ENA, SRA
from fastq_dl.exceptions import EmptyResultError, ProviderError
from fastq_dl.providers.ena import get_ena_metadata
from fastq_dl.providers.sra import get_sra_metadata


def _is_ena_empty_response(ena_data: list) -> bool:
    """Check if ENA returned HTTP 200 with no data rows."""
    return isinstance(ena_data, list) and len(ena_data) >= 2 and ena_data[0] == 200


def get_run_info(
    accession: str,
    query: str,
    provider: str,
    only_provider: bool,
    max_attempts: int = 3,
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

    # Configure primary/secondary provider order for fallback path
    def _fetch_ena():
        return get_ena_metadata(query)

    def _fetch_sra():
        return get_sra_metadata(accession)

    def _is_sra_empty(data):
        return not data

    if provider.lower() == "sra":
        primary_name, secondary_name = SRA, ENA
        primary_fetch, secondary_fetch = _fetch_sra, _fetch_ena
        primary_is_empty = _is_sra_empty
        secondary_is_empty = _is_ena_empty_response
    else:
        primary_name, secondary_name = ENA, SRA
        primary_fetch, secondary_fetch = _fetch_ena, _fetch_sra
        primary_is_empty = _is_ena_empty_response
        secondary_is_empty = _is_sra_empty

    primary_attempt = 1
    secondary_attempt = 1
    primary_exhausted = False
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
                elif _is_ena_empty_response(ena_data):
                    raise EmptyResultError(
                        f"No results from ENA for '{accession}'.",
                        provider="ENA",
                    )
                elif attempt >= max_attempts:
                    raise ProviderError(
                        f"Failed to query ENA after {max_attempts} attempts. "
                        f"STATUS: {ena_data[0]}, TEXT: {ena_data[1]}",
                        provider="ENA",
                        status_code=ena_data[0]
                        if isinstance(ena_data[0], int)
                        else None,
                    )
            else:
                success, sra_data = get_sra_metadata(accession)
                if success:
                    return SRA, sra_data
                elif not sra_data:
                    raise EmptyResultError(
                        f"No results from SRA for '{accession}'.",
                        provider="SRA",
                    )
                elif attempt >= max_attempts:
                    raise ProviderError(
                        f"Failed to query SRA after {max_attempts} attempts.",
                        provider="SRA",
                    )
            attempt += 1
            logging.warning(
                f"Querying {provider.lower()} was unsuccessful, retrying after ({sleep} seconds)"
            )
            time.sleep(sleep)
        else:
            if not primary_exhausted:
                logging.debug(
                    f"Querying {primary_name} for metadata (Attempt {primary_attempt} of {max_attempts})"
                )
                success, primary_data = primary_fetch()
                if success:
                    return primary_name, primary_data
                elif primary_is_empty(primary_data):
                    logging.info(
                        f"Accession not found in {primary_name}. Trying {secondary_name}..."
                    )
                    primary_exhausted = True
                elif primary_attempt >= max_attempts:
                    logging.debug(
                        f"Failed to get metadata from {primary_name}. Trying {secondary_name}..."
                    )
                    primary_exhausted = True
                else:
                    primary_attempt += 1
                    logging.warning(
                        f"Querying {primary_name} was unsuccessful, retrying after ({sleep} seconds)"
                    )
                    time.sleep(sleep)
                    continue

            logging.debug(
                f"Querying {secondary_name} for metadata (Attempt {secondary_attempt} of {max_attempts})"
            )
            success, secondary_data = secondary_fetch()
            if success:
                return secondary_name, secondary_data
            elif secondary_is_empty(secondary_data):
                raise EmptyResultError(
                    f"No results from {primary_name} or {secondary_name} for '{accession}'.",
                    provider=f"{primary_name}+{secondary_name}",
                )
            elif secondary_attempt >= max_attempts:
                raise ProviderError(
                    f"Failed to query both {primary_name} and {secondary_name} "
                    f"after {max_attempts} attempts each.",
                    provider=f"{primary_name}+{secondary_name}",
                )
            else:
                secondary_attempt += 1
                logging.warning(
                    f"Querying {secondary_name} was unsuccessful, retrying after ({sleep} seconds)"
                )
                time.sleep(sleep)
