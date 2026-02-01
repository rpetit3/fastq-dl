"""Integration tests that make real API calls.

These tests are marked with @pytest.mark.integration and are skipped
in normal test runs. They run during weekly scheduled CI.

Run these tests with: pytest -m integration
"""

import pytest
from pathlib import Path


@pytest.mark.integration
class TestENAIntegration:
    """Integration tests for ENA API."""

    def test_ena_metadata_query_run(self):
        """Test real ENA metadata query for a run accession."""
        from fastq_dl.providers.ena import get_ena_metadata

        success, data = get_ena_metadata("run_accession=SRR2838701")

        assert success is True
        assert len(data) >= 1
        assert data[0]["run_accession"] == "SRR2838701"

    def test_ena_metadata_query_experiment(self):
        """Test ENA experiment accession query."""
        from fastq_dl.providers.ena import get_ena_metadata

        success, data = get_ena_metadata("experiment_accession=SRX1390608")

        assert success is True
        assert len(data) >= 1

    def test_ena_metadata_query_study(self):
        """Test ENA study accession query."""
        from fastq_dl.providers.ena import get_ena_metadata

        success, data = get_ena_metadata(
            "(study_accession=SRP065366 OR secondary_study_accession=SRP065366)"
        )

        assert success is True
        assert len(data) >= 1

    def test_ena_metadata_query_sample(self):
        """Test ENA sample accession query."""
        from fastq_dl.providers.ena import get_ena_metadata

        success, data = get_ena_metadata(
            "(sample_accession=SAMN04215065 OR secondary_sample_accession=SAMN04215065)"
        )

        assert success is True
        assert len(data) >= 1

    def test_ena_metadata_invalid_accession(self):
        """Test ENA query with invalid accession returns no results."""
        from fastq_dl.providers.ena import get_ena_metadata

        success, data = get_ena_metadata("run_accession=SRR0000000")

        assert success is False


@pytest.mark.integration
class TestSRAIntegration:
    """Integration tests for SRA API."""

    @pytest.mark.filterwarnings("ignore::DeprecationWarning")
    def test_sra_metadata_query_run(self):
        """Test real SRA metadata query for a run accession."""
        from fastq_dl.providers.sra import get_sra_metadata

        success, data = get_sra_metadata("SRR2838701")

        assert success is True
        assert len(data) >= 1
        assert data[0]["run_accession"] == "SRR2838701"

    @pytest.mark.filterwarnings("ignore::DeprecationWarning")
    def test_sra_metadata_query_experiment(self):
        """Test SRA experiment accession query."""
        from fastq_dl.providers.sra import get_sra_metadata

        success, data = get_sra_metadata("SRX1390608")

        assert success is True
        assert len(data) >= 1

    @pytest.mark.filterwarnings("ignore::DeprecationWarning")
    def test_sra_metadata_invalid_accession(self):
        """Test SRA query with invalid accession returns no results."""
        from fastq_dl.providers.sra import get_sra_metadata

        success, data = get_sra_metadata("SRR0000000")

        assert success is False


@pytest.mark.integration
class TestValidateQuery:
    """Integration tests for query validation with real accessions."""

    def test_validate_project_accession(self):
        """Test validation of project accession."""
        from fastq_dl.utils import validate_query

        result = validate_query("PRJNA123456")
        assert "study_accession" in result
        assert "PRJNA123456" in result

    def test_validate_experiment_accession(self):
        """Test validation of experiment accession."""
        from fastq_dl.utils import validate_query

        result = validate_query("SRX1390608")
        assert result == "experiment_accession=SRX1390608"

    def test_validate_run_accession(self):
        """Test validation of run accession."""
        from fastq_dl.utils import validate_query

        result = validate_query("SRR2838701")
        assert result == "run_accession=SRR2838701"


@pytest.mark.integration
class TestGetRunInfo:
    """Integration tests for get_run_info function."""

    def test_get_run_info_ena(self):
        """Test get_run_info with ENA provider."""
        from fastq_dl.providers.generic import get_run_info
        from fastq_dl.constants import ENA

        source, data = get_run_info(
            "SRR2838701",
            "run_accession=SRR2838701",
            "ena",
            only_provider=True,
            max_attempts=3,
            sleep=5,
        )

        assert source == ENA
        assert len(data) >= 1

    @pytest.mark.filterwarnings("ignore::DeprecationWarning")
    def test_get_run_info_sra(self):
        """Test get_run_info with SRA provider."""
        from fastq_dl.providers.generic import get_run_info
        from fastq_dl.constants import SRA

        source, data = get_run_info(
            "SRR2838701",
            "run_accession=SRR2838701",
            "sra",
            only_provider=True,
            max_attempts=3,
            sleep=5,
        )

        assert source == SRA
        assert len(data) >= 1
