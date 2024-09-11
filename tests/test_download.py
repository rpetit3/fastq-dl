import pytest

from fastq_dl.download import get_ena_metadata, get_sra_metadata


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_get_sra_metadata_success():
    accession = "SRR2838701"
    success, metadata = get_sra_metadata(accession)
    assert success
    assert len(metadata) == 1
    metadata = metadata[0]

    assert metadata["run_accession"] == accession
    assert metadata["study_accession"] == "SRP065366"


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_get_sra_metadata_failure():
    accession = "SRR0000000"
    success, metadata = get_sra_metadata(accession)
    assert not success
    assert not metadata


def test_get_ena_metadata_success():
    accession = "SRR2838701"
    success, metadata = get_ena_metadata(f"run_accession={accession}")
    assert success
    assert len(metadata) == 1
    metadata = metadata[0]

    assert metadata["run_accession"] == accession
    assert metadata["sample_accession"] == "SAMN04215065"


def test_get_ena_metadata_failure():
    accession = "SRR0000000"
    success, metadata = get_ena_metadata(f"run_accession={accession}")
    assert not success
    assert metadata[1] == "Query was successful, but received an empty response"
