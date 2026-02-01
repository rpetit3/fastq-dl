"""Tests for generic provider functionality."""

import pytest
from unittest.mock import patch, MagicMock
import time

from fastq_dl.exceptions import ProviderError
from fastq_dl.providers.generic import get_run_info
from fastq_dl.constants import ENA, SRA


class TestGetRunInfo:
    """Tests for get_run_info function."""

    @patch("fastq_dl.providers.generic.get_ena_metadata")
    def test_ena_success(self, mock_ena, sample_ena_metadata):
        """Test successful ENA query."""
        mock_ena.return_value = (True, sample_ena_metadata)

        source, data = get_run_info(
            "SRR123456",
            "run_accession=SRR123456",
            "ena",
            only_provider=False,
        )

        assert source == ENA
        assert data == sample_ena_metadata

    @patch("fastq_dl.providers.generic.get_sra_metadata")
    def test_sra_only_provider(self, mock_sra, sample_sra_metadata):
        """Test SRA query with --only-provider."""
        mock_sra.return_value = (True, sample_sra_metadata)

        source, data = get_run_info(
            "SRR123456",
            "run_accession=SRR123456",
            "sra",
            only_provider=True,
        )

        assert source == SRA
        assert data == sample_sra_metadata

    @patch("fastq_dl.providers.generic.get_ena_metadata")
    def test_ena_only_provider(self, mock_ena, sample_ena_metadata):
        """Test ENA query with --only-provider."""
        mock_ena.return_value = (True, sample_ena_metadata)

        source, data = get_run_info(
            "SRR123456",
            "run_accession=SRR123456",
            "ena",
            only_provider=True,
        )

        assert source == ENA
        assert data == sample_ena_metadata

    @patch("fastq_dl.providers.generic.get_ena_metadata")
    @patch("fastq_dl.providers.generic.get_sra_metadata")
    def test_fallback_to_sra(self, mock_sra, mock_ena, sample_sra_metadata):
        """Test fallback to SRA when ENA fails."""
        mock_ena.return_value = (False, [500, "Error"])
        mock_sra.return_value = (True, sample_sra_metadata)

        source, data = get_run_info(
            "SRR123456",
            "run_accession=SRR123456",
            "ena",
            only_provider=False,
            max_attempts=1,
            sleep=0,
        )

        assert source == SRA

    @patch("fastq_dl.providers.generic.get_ena_metadata")
    @patch("time.sleep")
    def test_retry_logic_ena(self, mock_sleep, mock_ena, sample_ena_metadata):
        """Test retry logic on ENA failure."""
        # Fail twice, then succeed
        mock_ena.side_effect = [
            (False, [500, "Error"]),
            (False, [500, "Error"]),
            (True, sample_ena_metadata),
        ]

        source, data = get_run_info(
            "SRR123456",
            "run_accession=SRR123456",
            "ena",
            only_provider=True,
            max_attempts=3,
            sleep=1,
        )

        assert source == ENA
        assert mock_sleep.call_count == 2

    @patch("fastq_dl.providers.generic.get_sra_metadata")
    @patch("time.sleep")
    def test_retry_logic_sra(self, mock_sleep, mock_sra, sample_sra_metadata):
        """Test retry logic on SRA failure."""
        # Fail twice, then succeed
        mock_sra.side_effect = [
            (False, []),
            (False, []),
            (True, sample_sra_metadata),
        ]

        source, data = get_run_info(
            "SRR123456",
            "run_accession=SRR123456",
            "sra",
            only_provider=True,
            max_attempts=3,
            sleep=1,
        )

        assert source == SRA
        assert mock_sleep.call_count == 2

    @patch("fastq_dl.providers.generic.get_ena_metadata")
    def test_ena_only_provider_failure_raises(self, mock_ena):
        """Test that ENA-only failure raises ProviderError."""
        mock_ena.return_value = (False, [500, "Error"])

        with pytest.raises(ProviderError) as exc_info:
            get_run_info(
                "INVALID",
                "run_accession=INVALID",
                "ena",
                only_provider=True,
                max_attempts=1,
                sleep=0,
            )
        assert exc_info.value.provider == "ENA"

    @patch("fastq_dl.providers.generic.get_sra_metadata")
    def test_sra_only_provider_failure_raises(self, mock_sra):
        """Test that SRA-only failure raises ProviderError."""
        mock_sra.return_value = (False, [])

        with pytest.raises(ProviderError) as exc_info:
            get_run_info(
                "INVALID",
                "run_accession=INVALID",
                "sra",
                only_provider=True,
                max_attempts=1,
                sleep=0,
            )
        assert exc_info.value.provider == "SRA"

    @patch("fastq_dl.providers.generic.get_ena_metadata")
    @patch("fastq_dl.providers.generic.get_sra_metadata")
    def test_both_fail_raises(self, mock_sra, mock_ena):
        """Test that failure of both providers raises ProviderError."""
        mock_ena.return_value = (False, [500, "Error"])
        mock_sra.return_value = (False, [])

        with pytest.raises(ProviderError) as exc_info:
            get_run_info(
                "INVALID",
                "run_accession=INVALID",
                "ena",
                only_provider=False,
                max_attempts=1,
                sleep=0,
            )
        assert exc_info.value.provider == "ENA+SRA"

    @patch("fastq_dl.providers.generic.get_ena_metadata")
    @patch("fastq_dl.providers.generic.get_sra_metadata")
    @patch("time.sleep")
    def test_fallback_with_retry(self, mock_sleep, mock_sra, mock_ena, sample_sra_metadata):
        """Test fallback from ENA to SRA with retry on SRA."""
        # ENA fails max_attempts times, then SRA fails once, then succeeds
        mock_ena.return_value = (False, [500, "Error"])
        mock_sra.side_effect = [
            (False, []),
            (True, sample_sra_metadata),
        ]

        source, data = get_run_info(
            "SRR123456",
            "run_accession=SRR123456",
            "ena",
            only_provider=False,
            max_attempts=2,
            sleep=1,
        )

        assert source == SRA
        assert data == sample_sra_metadata
