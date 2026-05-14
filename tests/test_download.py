"""Legacy integration tests for metadata retrieval.

Note: These tests make real API calls and are marked as integration tests.
They are kept for backwards compatibility but new integration tests should
be added to test_integration.py.
"""

import pytest

from fastq_dl.providers.ena import get_ena_metadata
from fastq_dl.providers.sra import get_sra_metadata


@pytest.mark.integration
@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_get_sra_metadata_success():
    accession = "SRR2838701"
    success, metadata = get_sra_metadata(accession)
    assert success
    assert len(metadata) == 1
    metadata = metadata[0]

    assert metadata["run_accession"] == accession
    assert metadata["study_accession"] == "SRP065366"


@pytest.mark.integration
@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_get_sra_metadata_failure():
    accession = "SRR0000000"
    success, metadata = get_sra_metadata(accession)
    assert not success
    assert not metadata


@pytest.mark.integration
def test_get_ena_metadata_success():
    accession = "SRR2838701"
    success, metadata = get_ena_metadata(f"run_accession={accession}")
    assert success
    assert len(metadata) == 1
    metadata = metadata[0]

    assert metadata["run_accession"] == accession
    assert metadata["sample_accession"] == "SAMN04215065"


@pytest.mark.integration
def test_get_ena_metadata_failure():
    accession = "SRR0000000"
    success, metadata = get_ena_metadata(f"run_accession={accession}")
    assert not success
    assert metadata[1] == "Query was successful, but received an empty response"
