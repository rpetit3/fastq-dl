"""Tests for SRA provider functionality."""

from unittest.mock import MagicMock, patch

import pandas as pd

from fastq_dl.constants import SRA_FAILED
from fastq_dl.providers.sra import get_sra_metadata, sra_download


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

        pe1 = tmp_outdir / "SRR123456_1.fastq.gz"
        pe2 = tmp_outdir / "SRR123456_2.fastq.gz"
        se = tmp_outdir / "SRR123456.fastq.gz"
        pe1.write_bytes(b"read1")
        pe2.write_bytes(b"read2")
        se.write_bytes(b"orphan_reads")

        result = sra_download("SRR123456", str(tmp_outdir))

        assert result["r1"] == str(pe1)
        assert result["r2"] == str(pe2)
        assert result["single_end"] is False
        assert result["orphan"] == str(se)

    @patch("fastq_dl.providers.sra.execute")
    def test_download_failure(self, mock_execute, tmp_outdir):
        """Test handling of sracha command failure."""
        mock_execute.return_value = SRA_FAILED

        result = sra_download("SRR123456", str(tmp_outdir))

        assert result == SRA_FAILED

    @patch("fastq_dl.providers.sra.execute")
    def test_sra_lite_preference(self, mock_execute, tmp_outdir):
        """Test --sra-lite flag passes --format sralite to sracha."""
        se = tmp_outdir / "SRR123456.fastq.gz"

        def create_output_files(*args, **kwargs):
            se.write_bytes(b"reads")
            return 0

        mock_execute.side_effect = create_output_files

        sra_download("SRR123456", str(tmp_outdir), sra_lite=True)

        cmd = mock_execute.call_args_list[0][0][0]
        assert "sracha" in cmd
        assert "--format" in cmd
        assert "sralite" in cmd

    @patch("fastq_dl.providers.sra.execute")
    def test_sra_normalized_preference(self, mock_execute, tmp_outdir):
        """Test default SRA Normalized preference (no --format sralite)."""
        se = tmp_outdir / "SRR123456.fastq.gz"

        def create_output_files(*args, **kwargs):
            se.write_bytes(b"reads")
            return 0

        mock_execute.side_effect = create_output_files

        sra_download("SRR123456", str(tmp_outdir), sra_lite=False)

        cmd = mock_execute.call_args_list[0][0][0]
        assert "sracha" in cmd
        assert "sralite" not in cmd

    def test_skip_existing_paired_files(self, tmp_outdir):
        """Test skipping download when paired files already exist."""
        pe1 = tmp_outdir / "SRR123456_1.fastq.gz"
        pe2 = tmp_outdir / "SRR123456_2.fastq.gz"
        pe1.write_bytes(b"read1")
        pe2.write_bytes(b"read2")

        with patch("fastq_dl.providers.sra.execute"):
            result = sra_download("SRR123456", str(tmp_outdir))

        assert result["r1"] == str(pe1)
        assert result["r2"] == str(pe2)

    def test_skip_existing_single_file(self, tmp_outdir):
        """Test skipping download when single-end file already exists."""
        se = tmp_outdir / "SRR123456.fastq.gz"
        se.write_bytes(b"reads")

        with patch("fastq_dl.providers.sra.execute"):
            result = sra_download("SRR123456", str(tmp_outdir))

        assert result["r1"] == str(se)
        assert result["single_end"] is True

    @patch("fastq_dl.providers.sra.execute")
    def test_force_removes_existing_files(self, mock_execute, tmp_outdir):
        """Test --force flag removes existing files before download."""
        mock_execute.return_value = 0

        pe1 = tmp_outdir / "SRR123456_1.fastq.gz"
        pe2 = tmp_outdir / "SRR123456_2.fastq.gz"
        pe1.write_bytes(b"old")
        pe2.write_bytes(b"old")

        def create_output_files(*args, **kwargs):
            pe1.write_bytes(b"new")
            pe2.write_bytes(b"new")
            return 0

        mock_execute.side_effect = create_output_files

        sra_download("SRR123456", str(tmp_outdir), force=True)

        assert mock_execute.called
        cmd = mock_execute.call_args_list[0][0][0]
        assert "--force" in cmd

    @patch("fastq_dl.providers.sra.execute")
    def test_cpus_parameter(self, mock_execute, tmp_outdir):
        """Test cpus parameter is passed to sracha as --threads."""
        se = tmp_outdir / "SRR123456.fastq.gz"

        def create_output_files(*args, **kwargs):
            se.write_bytes(b"reads")
            return 0

        mock_execute.side_effect = create_output_files

        sra_download("SRR123456", str(tmp_outdir), cpus=4)

        cmd = mock_execute.call_args_list[0][0][0]
        threads_idx = cmd.index("--threads")
        assert cmd[threads_idx + 1] == "4"

    @patch("fastq_dl.providers.sra.execute")
    def test_connections_parameter(self, mock_execute, tmp_outdir):
        """Test connections parameter is passed to sracha as --connections."""
        se = tmp_outdir / "SRR123456.fastq.gz"

        def create_output_files(*args, **kwargs):
            se.write_bytes(b"reads")
            return 0

        mock_execute.side_effect = create_output_files

        sra_download("SRR123456", str(tmp_outdir), cpus=4, connections=12)

        cmd = mock_execute.call_args_list[0][0][0]
        threads_idx = cmd.index("--threads")
        assert cmd[threads_idx + 1] == "4"
        connections_idx = cmd.index("--connections")
        assert cmd[connections_idx + 1] == "12"

    @patch("fastq_dl.providers.sra.execute")
    def test_no_strict_passes_flag(self, mock_execute, tmp_outdir):
        """Test no_strict flag adds --no-strict to sracha command."""
        se = tmp_outdir / "SRR123456.fastq.gz"

        def create_output_files(*args, **kwargs):
            se.write_bytes(b"reads")
            return 0

        mock_execute.side_effect = create_output_files

        sra_download("SRR123456", str(tmp_outdir), no_strict=True)

        cmd = mock_execute.call_args_list[0][0][0]
        assert "--no-strict" in cmd

    @patch("fastq_dl.providers.sra.execute")
    def test_no_strict_not_passed_by_default(self, mock_execute, tmp_outdir):
        """Test --no-strict is not passed when no_strict=False."""
        se = tmp_outdir / "SRR123456.fastq.gz"

        def create_output_files(*args, **kwargs):
            se.write_bytes(b"reads")
            return 0

        mock_execute.side_effect = create_output_files

        sra_download("SRR123456", str(tmp_outdir), no_strict=False)

        cmd = mock_execute.call_args_list[0][0][0]
        assert "--no-strict" not in cmd

    @patch("fastq_dl.providers.sra.execute")
    def test_skip_compress_produces_uncompressed_files(self, mock_execute, tmp_outdir):
        """Test compress=False passes --no-gzip and uses .fastq suffixes."""
        se = tmp_outdir / "SRR123456.fastq"

        def create_output_files(*args, **kwargs):
            se.write_bytes(b"reads")
            return 0

        mock_execute.side_effect = create_output_files

        result = sra_download("SRR123456", str(tmp_outdir), compress=False)

        assert result["r1"] == str(se)
        assert result["single_end"] is True
        assert str(result["r1"]).endswith(".fastq")
        cmd = mock_execute.call_args_list[0][0][0]
        assert "--no-gzip" in cmd

    @patch("fastq_dl.providers.sra.execute")
    def test_skip_compress_paired_end(self, mock_execute, tmp_outdir):
        """Test compress=False with paired-end data uses .fastq suffixes."""
        mock_execute.return_value = 0

        pe1 = tmp_outdir / "SRR123456_1.fastq"
        pe2 = tmp_outdir / "SRR123456_2.fastq"
        pe1.write_bytes(b"read1")
        pe2.write_bytes(b"read2")

        result = sra_download("SRR123456", str(tmp_outdir), compress=False)

        assert result["r1"] == str(pe1)
        assert result["r2"] == str(pe2)
        assert result["single_end"] is False
        assert str(result["r1"]).endswith("_1.fastq")
        assert str(result["r2"]).endswith("_2.fastq")

    @patch("fastq_dl.providers.sra.execute")
    def test_gzip_level_default(self, mock_execute, tmp_outdir):
        """Test default gzip level (1) does not add --gzip-level flag."""
        se = tmp_outdir / "SRR123456.fastq.gz"

        def create_output_files(*args, **kwargs):
            se.write_bytes(b"reads")
            return 0

        mock_execute.side_effect = create_output_files

        sra_download("SRR123456", str(tmp_outdir), gzip_level=1)

        cmd = mock_execute.call_args_list[0][0][0]
        assert "--gzip-level" not in cmd

    @patch("fastq_dl.providers.sra.execute")
    def test_gzip_level_custom(self, mock_execute, tmp_outdir):
        """Test custom gzip level is passed to sracha."""
        se = tmp_outdir / "SRR123456.fastq.gz"

        def create_output_files(*args, **kwargs):
            se.write_bytes(b"reads")
            return 0

        mock_execute.side_effect = create_output_files

        sra_download("SRR123456", str(tmp_outdir), gzip_level=6)

        cmd = mock_execute.call_args_list[0][0][0]
        gzip_idx = cmd.index("--gzip-level")
        assert cmd[gzip_idx + 1] == "6"

    @patch("fastq_dl.providers.sra.execute")
    def test_sracha_command_structure(self, mock_execute, tmp_outdir):
        """Test the full sracha get command is constructed correctly."""
        se = tmp_outdir / "SRR123456.fastq.gz"

        def create_output_files(*args, **kwargs):
            se.write_bytes(b"reads")
            return 0

        mock_execute.side_effect = create_output_files

        sra_download("SRR123456", str(tmp_outdir), cpus=2)

        cmd = mock_execute.call_args_list[0][0][0]
        assert cmd[0] == "sracha"
        assert cmd[1] == "get"
        assert cmd[2] == "SRR123456"
        assert "-O" in cmd
        assert "--no-progress" in cmd
        assert "-y" in cmd
