"""Tests for ENA provider functionality."""

import pytest
from unittest.mock import patch, MagicMock
from pathlib import Path
import responses

from fastq_dl.exceptions import DownloadError
from fastq_dl.providers.ena import get_ena_metadata, ena_download, download_ena_fastq
from fastq_dl.constants import ENA_FAILED, ENA_URL


class TestGetEnaMetadata:
    """Tests for get_ena_metadata function."""

    @responses.activate
    def test_successful_query(self):
        """Test successful metadata retrieval from ENA."""
        mock_response = "run_accession\tsample_accession\nSRR123456\tSAMN123456"
        responses.add(
            responses.GET,
            ENA_URL,
            body=mock_response,
            status=200,
            match_querystring=False,
        )

        success, data = get_ena_metadata("run_accession=SRR123456")

        assert success is True
        assert len(data) == 1
        assert data[0]["run_accession"] == "SRR123456"
        assert data[0]["sample_accession"] == "SAMN123456"

    @responses.activate
    def test_multiple_results(self):
        """Test handling of multiple results from ENA."""
        mock_response = "run_accession\tsample_accession\nSRR123456\tSAMN123456\nSRR123457\tSAMN123457"
        responses.add(
            responses.GET,
            ENA_URL,
            body=mock_response,
            status=200,
            match_querystring=False,
        )

        success, data = get_ena_metadata("study_accession=SRP123456")

        assert success is True
        assert len(data) == 2

    @responses.activate
    def test_empty_response(self):
        """Test handling of empty response from ENA."""
        responses.add(
            responses.GET,
            ENA_URL,
            body="",
            status=200,
            match_querystring=False,
        )

        success, data = get_ena_metadata("run_accession=INVALID")

        assert success is False
        assert "empty response" in data[1].lower()

    @responses.activate
    def test_http_error(self):
        """Test handling of HTTP errors."""
        responses.add(
            responses.GET,
            ENA_URL,
            body="Server Error",
            status=500,
            match_querystring=False,
        )

        success, data = get_ena_metadata("run_accession=SRR123456")

        assert success is False
        assert data[0] == 500


class TestEnaDownload:
    """Tests for ena_download function."""

    def test_no_ftp_url(self, tmp_outdir):
        """Test handling when no FTP URL is available."""
        metadata = {
            "run_accession": "SRR123456",
            "fastq_ftp": "",
            "fastq_md5": "",
            "library_layout": "PAIRED",
        }

        result = ena_download(metadata, str(tmp_outdir))

        assert result == ENA_FAILED

    @patch("fastq_dl.providers.ena.download_ena_fastq")
    def test_paired_end_download(self, mock_download, sample_ena_metadata, tmp_outdir):
        """Test paired-end FASTQ download."""
        mock_download.return_value = str(tmp_outdir / "test.fastq.gz")

        result = ena_download(sample_ena_metadata[0], str(tmp_outdir))

        assert mock_download.call_count == 2  # R1 and R2
        assert "r1" in result
        assert "r2" in result
        assert result["single_end"] is False
        assert result["orphan"] is None  # ENA doesn't support orphan reads

    @patch("fastq_dl.providers.ena.download_ena_fastq")
    def test_single_end_download(self, mock_download, single_end_ena_metadata, tmp_outdir):
        """Test single-end FASTQ download."""
        mock_download.return_value = str(tmp_outdir / "test.fastq.gz")

        result = ena_download(single_end_ena_metadata[0], str(tmp_outdir))

        assert mock_download.call_count == 1
        assert result["single_end"] is True
        assert result["orphan"] is None

    @patch("fastq_dl.providers.ena.download_ena_fastq")
    def test_download_failure(self, mock_download, sample_ena_metadata, tmp_outdir):
        """Test handling of download failure."""
        mock_download.return_value = ENA_FAILED

        result = ena_download(sample_ena_metadata[0], str(tmp_outdir))

        assert result == ENA_FAILED


class TestDownloadEnaFastq:
    """Tests for download_ena_fastq function."""

    @patch("fastq_dl.providers.ena.execute")
    @patch("fastq_dl.providers.ena.md5sum")
    def test_successful_download(self, mock_md5sum, mock_execute, tmp_outdir):
        """Test successful FASTQ download."""
        mock_execute.return_value = 0
        mock_md5sum.return_value = "abc123"

        result = download_ena_fastq(
            "ftp.example.com/test.fastq.gz",
            str(tmp_outdir),
            "abc123",
        )

        assert "test.fastq.gz" in result
        assert mock_execute.called

    @patch("fastq_dl.providers.ena.md5sum")
    def test_skip_existing_matching_md5(self, mock_md5sum, tmp_outdir):
        """Test skipping re-download when MD5 matches."""
        # Create existing file
        existing_file = tmp_outdir / "test.fastq.gz"
        existing_file.write_bytes(b"content")
        mock_md5sum.return_value = "abc123"

        result = download_ena_fastq(
            "ftp.example.com/test.fastq.gz",
            str(tmp_outdir),
            "abc123",
        )

        assert result == str(existing_file)

    @patch("fastq_dl.providers.ena.execute")
    @patch("fastq_dl.providers.ena.md5sum")
    def test_force_overwrite(self, mock_md5sum, mock_execute, tmp_outdir):
        """Test --force flag overwrites existing files."""
        existing_file = tmp_outdir / "test.fastq.gz"
        existing_file.write_bytes(b"old content")
        mock_execute.return_value = 0
        mock_md5sum.return_value = "abc123"

        result = download_ena_fastq(
            "ftp.example.com/test.fastq.gz",
            str(tmp_outdir),
            "abc123",
            force=True,
        )

        # execute should have been called for download
        assert mock_execute.called

    def test_ignore_md5_flag(self, tmp_outdir):
        """Test --ignore flag skips MD5 verification."""
        existing_file = tmp_outdir / "test.fastq.gz"
        existing_file.write_bytes(b"content")

        result = download_ena_fastq(
            "ftp.example.com/test.fastq.gz",
            str(tmp_outdir),
            "wrong_md5",
            ignore_md5=True,
        )

        assert result == str(existing_file)

    @patch("fastq_dl.providers.ena.execute")
    def test_download_failure_returns_failed(self, mock_execute, tmp_outdir):
        """Test that download failure returns ENA_FAILED."""
        mock_execute.return_value = ENA_FAILED

        result = download_ena_fastq(
            "ftp.example.com/test.fastq.gz",
            str(tmp_outdir),
            "abc123",
        )

        assert result == ENA_FAILED

    @patch("fastq_dl.providers.ena.execute")
    @patch("fastq_dl.providers.ena.md5sum")
    def test_md5_mismatch_raises_after_max_attempts(self, mock_md5sum, mock_execute, tmp_outdir):
        """Test that repeated MD5 mismatch raises DownloadError."""
        mock_execute.return_value = 0
        # MD5 never matches
        mock_md5sum.return_value = "wrong_md5"

        with pytest.raises(DownloadError) as exc_info:
            download_ena_fastq(
                "ftp.example.com/test.fastq.gz",
                str(tmp_outdir),
                "correct_md5",
                max_attempts=2,
            )

        assert exc_info.value.provider == "ENA"
        assert "MD5 checksum mismatch" in str(exc_info.value)
