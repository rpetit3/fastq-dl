"""Tests for SRA provider functionality."""

import pytest
from unittest.mock import patch, MagicMock
from pathlib import Path
import pandas as pd

from fastq_dl.providers.sra import get_sra_metadata, sra_download
from fastq_dl.constants import SRA_FAILED


class TestGetSraMetadata:
    """Tests for get_sra_metadata function."""

    @patch("fastq_dl.providers.sra.SRAweb")
    def test_successful_query(self, mock_sraweb_class):
        """Test successful metadata retrieval from SRA."""
        mock_db = MagicMock()
        mock_sraweb_class.return_value = mock_db
        mock_db.search_sra.return_value = pd.DataFrame(
            {
                "run_accession": ["SRR123456"],
                "sample_accession": ["SAMN123456"],
            }
        )

        success, data = get_sra_metadata("SRR123456")

        assert success is True
        assert len(data) == 1
        assert data[0]["run_accession"] == "SRR123456"

    @patch("fastq_dl.providers.sra.SRAweb")
    def test_multiple_results(self, mock_sraweb_class):
        """Test handling of multiple results from SRA."""
        mock_db = MagicMock()
        mock_sraweb_class.return_value = mock_db
        mock_db.search_sra.return_value = pd.DataFrame(
            {
                "run_accession": ["SRR123456", "SRR123457"],
                "sample_accession": ["SAMN123456", "SAMN123457"],
            }
        )

        success, data = get_sra_metadata("SRX123456")

        assert success is True
        assert len(data) == 2

    @patch("fastq_dl.providers.sra.SRAweb")
    def test_no_results(self, mock_sraweb_class):
        """Test handling when no results found."""
        mock_db = MagicMock()
        mock_sraweb_class.return_value = mock_db
        mock_db.search_sra.return_value = None

        success, data = get_sra_metadata("INVALID")

        assert success is False
        assert data == []


class TestSraDownload:
    """Tests for sra_download function."""

    @patch("fastq_dl.providers.sra.execute")
    def test_successful_paired_download(self, mock_execute, tmp_outdir):
        """Test successful paired-end download from SRA."""
        mock_execute.return_value = 0

        # Create mock output files that would be created by fasterq-dump
        pe1 = tmp_outdir / "SRR123456_1.fastq.gz"
        pe2 = tmp_outdir / "SRR123456_2.fastq.gz"
        pe1.write_bytes(b"read1")
        pe2.write_bytes(b"read2")

        result = sra_download("SRR123456", str(tmp_outdir))

        assert result["r1"] == str(pe1)
        assert result["r2"] == str(pe2)
        assert result["single_end"] is False
        assert result["orphan"] is None

    @patch("fastq_dl.providers.sra.execute")
    def test_successful_single_download(self, mock_execute, tmp_outdir):
        """Test successful single-end download from SRA."""
        mock_execute.return_value = 0

        # Create mock output file
        se = tmp_outdir / "SRR123456.fastq.gz"
        se.write_bytes(b"reads")

        result = sra_download("SRR123456", str(tmp_outdir))

        assert result["r1"] == str(se)
        assert result["single_end"] is True
        assert result["orphan"] is None

    @patch("fastq_dl.providers.sra.execute")
    def test_paired_with_orphan_reads(self, mock_execute, tmp_outdir):
        """Test paired-end download with orphan reads file."""
        mock_execute.return_value = 0

        # Create mock output files including orphan reads
        pe1 = tmp_outdir / "SRR123456_1.fastq.gz"
        pe2 = tmp_outdir / "SRR123456_2.fastq.gz"
        se = tmp_outdir / "SRR123456.fastq.gz"  # Orphan reads
        pe1.write_bytes(b"read1")
        pe2.write_bytes(b"read2")
        se.write_bytes(b"orphan_reads")

        result = sra_download("SRR123456", str(tmp_outdir))

        assert result["r1"] == str(pe1)
        assert result["r2"] == str(pe2)
        assert result["single_end"] is False
        assert result["orphan"] == str(se)

    @patch("fastq_dl.providers.sra.execute")
    def test_prefetch_failure(self, mock_execute, tmp_outdir):
        """Test handling of prefetch command failure."""
        mock_execute.return_value = SRA_FAILED

        result = sra_download("SRR123456", str(tmp_outdir))

        assert result == SRA_FAILED

    @patch("fastq_dl.providers.sra.execute")
    def test_sra_lite_preference(self, mock_execute, tmp_outdir):
        """Test --sra-lite flag sets correct preference."""
        se = tmp_outdir / "SRR123456.fastq.gz"
        sra_file = tmp_outdir / "SRR123456.sra"

        def create_output_files(*args, **kwargs):
            se.write_bytes(b"reads")
            sra_file.write_bytes(b"sra")
            return 0

        mock_execute.side_effect = create_output_files

        sra_download("SRR123456", str(tmp_outdir), sra_lite=True)

        # Check that vdb-config was called with 'yes'
        calls = [c[0][0] for c in mock_execute.call_args_list]  # Get first positional arg of each call
        assert any(
            isinstance(c, list) and "vdb-config" in c and "yes" in c for c in calls
        )

    @patch("fastq_dl.providers.sra.execute")
    def test_sra_normalized_preference(self, mock_execute, tmp_outdir):
        """Test default SRA Normalized preference."""
        se = tmp_outdir / "SRR123456.fastq.gz"
        sra_file = tmp_outdir / "SRR123456.sra"

        def create_output_files(*args, **kwargs):
            se.write_bytes(b"reads")
            sra_file.write_bytes(b"sra")
            return 0

        mock_execute.side_effect = create_output_files

        sra_download("SRR123456", str(tmp_outdir), sra_lite=False)

        # Check that vdb-config was called with 'no'
        calls = [c[0][0] for c in mock_execute.call_args_list]  # Get first positional arg of each call
        assert any(
            isinstance(c, list) and "vdb-config" in c and "no" in c for c in calls
        )

    def test_skip_existing_paired_files(self, tmp_outdir):
        """Test skipping download when paired files already exist."""
        pe1 = tmp_outdir / "SRR123456_1.fastq.gz"
        pe2 = tmp_outdir / "SRR123456_2.fastq.gz"
        pe1.write_bytes(b"read1")
        pe2.write_bytes(b"read2")

        with patch("fastq_dl.providers.sra.execute") as mock_execute:
            result = sra_download("SRR123456", str(tmp_outdir))

        # execute should not have been called for download commands
        assert result["r1"] == str(pe1)
        assert result["r2"] == str(pe2)

    def test_skip_existing_single_file(self, tmp_outdir):
        """Test skipping download when single-end file already exists."""
        se = tmp_outdir / "SRR123456.fastq.gz"
        se.write_bytes(b"reads")

        with patch("fastq_dl.providers.sra.execute") as mock_execute:
            result = sra_download("SRR123456", str(tmp_outdir))

        assert result["r1"] == str(se)
        assert result["single_end"] is True

    @patch("fastq_dl.providers.sra.execute")
    def test_force_removes_existing_files(self, mock_execute, tmp_outdir):
        """Test --force flag removes existing files before download."""
        mock_execute.return_value = 0

        pe1 = tmp_outdir / "SRR123456_1.fastq.gz"
        pe2 = tmp_outdir / "SRR123456_2.fastq.gz"
        sra_file = tmp_outdir / "SRR123456.sra"
        pe1.write_bytes(b"old")
        pe2.write_bytes(b"old")
        sra_file.write_bytes(b"sra")  # Create .sra file that gets unlinked after download

        # Create output files that would be created by fasterq-dump
        def create_output_files(*args, **kwargs):
            pe1.write_bytes(b"new")
            pe2.write_bytes(b"new")
            return 0

        mock_execute.side_effect = create_output_files

        sra_download("SRR123456", str(tmp_outdir), force=True)

        # execute should have been called since force=True
        assert mock_execute.called

    @patch("fastq_dl.providers.sra.execute")
    def test_fasterq_dump_failure(self, mock_execute, tmp_outdir):
        """Test handling of fasterq-dump command failure."""
        # First call (vdb-config) succeeds, second call (prefetch) succeeds,
        # third call (fasterq-dump) fails
        mock_execute.side_effect = [0, 0, SRA_FAILED]

        result = sra_download("SRR123456", str(tmp_outdir))

        assert result == SRA_FAILED

    @patch("fastq_dl.providers.sra.execute")
    def test_cpus_parameter(self, mock_execute, tmp_outdir):
        """Test cpus parameter is passed to fasterq-dump."""
        se = tmp_outdir / "SRR123456.fastq.gz"
        sra_file = tmp_outdir / "SRR123456.sra"

        def create_output_files(*args, **kwargs):
            se.write_bytes(b"reads")
            sra_file.write_bytes(b"sra")
            return 0

        mock_execute.side_effect = create_output_files

        sra_download("SRR123456", str(tmp_outdir), cpus=4)

        # Get first positional arg of each call
        calls = [c[0][0] for c in mock_execute.call_args_list]
        # Check for fasterq-dump command with --threads and 4
        assert any(
            isinstance(c, list) and "fasterq-dump" in c and "--threads" in c and "4" in c
            for c in calls
        )

    @patch("fastq_dl.providers.sra.execute")
    def test_ignore_md5_skips_verification(self, mock_execute, tmp_outdir):
        """Test ignore_md5 flag adds --verify no to prefetch command."""
        se = tmp_outdir / "SRR123456.fastq.gz"
        sra_file = tmp_outdir / "SRR123456.sra"

        def create_output_files(*args, **kwargs):
            se.write_bytes(b"reads")
            sra_file.write_bytes(b"sra")
            return 0

        mock_execute.side_effect = create_output_files

        sra_download("SRR123456", str(tmp_outdir), ignore_md5=True)

        # Get all commands passed to execute
        calls = [c[0][0] for c in mock_execute.call_args_list]

        # Find the prefetch command and verify --verify no is present
        prefetch_calls = [c for c in calls if isinstance(c, list) and "prefetch" in c]
        assert len(prefetch_calls) == 1
        # Check that --verify is followed by no
        prefetch_cmd = prefetch_calls[0]
        verify_idx = prefetch_cmd.index("--verify")
        assert prefetch_cmd[verify_idx + 1] == "no"

    @patch("fastq_dl.providers.sra.execute")
    def test_md5_verification_enabled_by_default(self, mock_execute, tmp_outdir):
        """Test MD5 verification is enabled by default (--verify yes)."""
        se = tmp_outdir / "SRR123456.fastq.gz"
        sra_file = tmp_outdir / "SRR123456.sra"

        def create_output_files(*args, **kwargs):
            se.write_bytes(b"reads")
            sra_file.write_bytes(b"sra")
            return 0

        mock_execute.side_effect = create_output_files

        sra_download("SRR123456", str(tmp_outdir), ignore_md5=False)

        calls = [c[0][0] for c in mock_execute.call_args_list]
        prefetch_calls = [c for c in calls if isinstance(c, list) and "prefetch" in c]
        assert len(prefetch_calls) == 1
        # Check that --verify is followed by yes
        prefetch_cmd = prefetch_calls[0]
        verify_idx = prefetch_cmd.index("--verify")
        assert prefetch_cmd[verify_idx + 1] == "yes"
