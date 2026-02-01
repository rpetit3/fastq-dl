"""Shared fixtures for fastq-dl tests."""

import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch


# ============================================================================
# Path Fixtures
# ============================================================================


@pytest.fixture
def tmp_outdir(tmp_path):
    """Create a temporary output directory."""
    outdir = tmp_path / "outdir"
    outdir.mkdir()
    return outdir


@pytest.fixture
def sample_fastq_content():
    """Sample FASTQ content for testing."""
    return b"@read1\nACGT\n+\n1234\n@read2\nTGCA\n+\n4321\n"


@pytest.fixture
def mock_fastq_files(tmp_path, sample_fastq_content):
    """Create mock FASTQ files for testing."""
    r1 = tmp_path / "SRR123456_1.fastq.gz"
    r2 = tmp_path / "SRR123456_2.fastq.gz"

    r1.write_bytes(sample_fastq_content)
    r2.write_bytes(sample_fastq_content)

    return {"r1": str(r1), "r2": str(r2), "single_end": False, "orphan": None}


# ============================================================================
# Mock Metadata Fixtures
# ============================================================================


@pytest.fixture
def sample_ena_metadata():
    """Sample ENA metadata response."""
    return [
        {
            "run_accession": "SRR2838701",
            "experiment_accession": "SRX1390608",
            "sample_accession": "SAMN04215065",
            "study_accession": "SRP065366",
            "library_layout": "PAIRED",
            "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/SRR283/001/SRR2838701/SRR2838701_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR283/001/SRR2838701/SRR2838701_2.fastq.gz",
            "fastq_md5": "abc123;def456",
        }
    ]


@pytest.fixture
def sample_sra_metadata():
    """Sample SRA metadata response."""
    return [
        {
            "run_accession": "SRR2838701",
            "experiment_accession": "SRX1390608",
            "sample_accession": "SAMN04215065",
            "study_accession": "SRP065366",
        }
    ]


@pytest.fixture
def single_end_ena_metadata():
    """Sample single-end ENA metadata."""
    return [
        {
            "run_accession": "SRR123456",
            "experiment_accession": "SRX123456",
            "sample_accession": "SAMN123456",
            "study_accession": "SRP123456",
            "library_layout": "SINGLE",
            "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/SRR123/SRR123456.fastq.gz",
            "fastq_md5": "abc123",
        }
    ]


# ============================================================================
# Mock External Command Fixtures
# ============================================================================


@pytest.fixture
def mock_execute_success():
    """Mock successful command execution."""
    with patch("fastq_dl.utils.execute") as mock:
        mock.return_value = 0
        yield mock


@pytest.fixture
def mock_execute_failure():
    """Mock failed command execution."""
    with patch("fastq_dl.utils.execute") as mock:
        mock.return_value = "ENA_NOT_FOUND"
        yield mock


@pytest.fixture
def mock_requests_get():
    """Mock requests.get for ENA API calls."""
    with patch("requests.get") as mock:
        yield mock


@pytest.fixture
def mock_sraweb():
    """Mock pysradb SRAweb class."""
    with patch("fastq_dl.providers.sra.SRAweb") as mock:
        yield mock


# ============================================================================
# Integration Test Markers
# ============================================================================


def pytest_configure(config):
    """Configure custom markers."""
    config.addinivalue_line("markers", "integration: mark test as integration test")
    config.addinivalue_line("markers", "slow: mark test as slow running")
