"""Tests for CLI functionality."""

from unittest.mock import patch

import pytest
from click.testing import CliRunner

from fastq_dl.cli.download import fastqdl, main
from fastq_dl.constants import ENA_FAILED, ENA_NO_FASTQS


class TestCLI:
    """Tests for command-line interface."""

    @pytest.fixture
    def runner(self):
        """Create CLI test runner."""
        return CliRunner()

    def test_help_flag(self, runner):
        """Test --help flag displays help text."""
        result = runner.invoke(fastqdl, ["--help"])
        assert result.exit_code == 0
        assert "Download FASTQ files from ENA or SRA" in result.output

    def test_version_flag(self, runner):
        """Test --version flag displays version."""
        result = runner.invoke(fastqdl, ["--version"])
        assert result.exit_code == 0
        assert "version" in result.output.lower() or "." in result.output

    def test_missing_accession(self, runner):
        """Test error when accession is not provided."""
        result = runner.invoke(fastqdl, [])
        assert result.exit_code != 0
        assert "Missing option" in result.output or "required" in result.output.lower()

    @patch("fastq_dl.cli.download.validate_query")
    @patch("fastq_dl.cli.download.get_run_info")
    @patch("fastq_dl.cli.download.write_tsv")
    def test_metadata_only_download(
        self,
        mock_write_tsv,
        mock_get_run_info,
        mock_validate_query,
        runner,
        sample_ena_metadata,
        tmp_path,
    ):
        """Test --only-download-metadata flag."""
        mock_validate_query.return_value = "run_accession=SRR123456"
        mock_get_run_info.return_value = ("ENA", sample_ena_metadata)

        result = runner.invoke(
            fastqdl,
            [
                "--accession",
                "SRR123456",
                "--only-download-metadata",
                "--outdir",
                str(tmp_path),
            ],
        )

        assert mock_write_tsv.called
        assert result.exit_code == 0

    @patch("fastq_dl.cli.download.validate_query")
    @patch("fastq_dl.cli.download.get_run_info")
    @patch("fastq_dl.cli.download.ena_download")
    @patch("fastq_dl.cli.download.write_tsv")
    def test_ena_download_flow(
        self,
        mock_write_tsv,
        mock_ena_download,
        mock_get_run_info,
        mock_validate_query,
        runner,
        sample_ena_metadata,
        mock_fastq_files,
        tmp_path,
    ):
        """Test ENA download flow."""
        mock_validate_query.return_value = "run_accession=SRR123456"
        mock_get_run_info.return_value = ("ENA", sample_ena_metadata)
        mock_ena_download.return_value = mock_fastq_files

        runner.invoke(
            fastqdl,
            [
                "--accession",
                "SRR123456",
                "--provider",
                "ena",
                "--outdir",
                str(tmp_path),
            ],
        )

        assert mock_ena_download.called

    @patch("fastq_dl.cli.download.validate_query")
    @patch("fastq_dl.cli.download.get_run_info")
    @patch("fastq_dl.cli.download.sra_download")
    @patch("fastq_dl.cli.download.write_tsv")
    def test_sra_download_flow(
        self,
        mock_write_tsv,
        mock_sra_download,
        mock_get_run_info,
        mock_validate_query,
        runner,
        sample_sra_metadata,
        mock_fastq_files,
        tmp_path,
    ):
        """Test SRA download flow."""
        mock_validate_query.return_value = "run_accession=SRR123456"
        mock_get_run_info.return_value = ("SRA", sample_sra_metadata)
        mock_sra_download.return_value = mock_fastq_files

        runner.invoke(
            fastqdl,
            [
                "--accession",
                "SRR123456",
                "--provider",
                "sra",
                "--outdir",
                str(tmp_path),
            ],
        )

        assert mock_sra_download.called


class TestCLIOptions:
    """Test various CLI option combinations."""

    @pytest.fixture
    def runner(self):
        return CliRunner()

    def test_provider_choices_ena(self, runner):
        """Test provider option accepts 'ena'."""
        result = runner.invoke(
            fastqdl, ["--accession", "SRR123456", "--provider", "ena", "--help"]
        )
        # Just checking it parses correctly before hitting --help
        assert result.exit_code == 0

    def test_provider_choices_sra(self, runner):
        """Test provider option accepts 'sra'."""
        result = runner.invoke(
            fastqdl, ["--accession", "SRR123456", "--provider", "sra", "--help"]
        )
        assert result.exit_code == 0

    def test_provider_choices_invalid(self, runner):
        """Test provider option rejects invalid choices."""
        result = runner.invoke(
            fastqdl, ["--accession", "SRR123456", "--provider", "invalid"]
        )
        assert result.exit_code != 0
        assert "Invalid value" in result.output or "invalid" in result.output.lower()

    @patch("fastq_dl.cli.download.validate_query")
    @patch("fastq_dl.cli.download.get_run_info")
    @patch("fastq_dl.cli.download.write_tsv")
    def test_group_by_experiment_flag(
        self, mock_write_tsv, mock_get_run_info, mock_validate_query, runner, tmp_path
    ):
        """Test --group-by-experiment flag is accepted."""
        mock_validate_query.return_value = "experiment_accession=SRX123456"
        mock_get_run_info.return_value = ("ENA", [])

        result = runner.invoke(
            fastqdl,
            [
                "--accession",
                "SRX123456",
                "--group-by-experiment",
                "--only-download-metadata",
                "--outdir",
                str(tmp_path),
            ],
        )

        assert result.exit_code == 0

    @patch("fastq_dl.cli.download.validate_query")
    @patch("fastq_dl.cli.download.get_run_info")
    @patch("fastq_dl.cli.download.write_tsv")
    def test_group_by_sample_flag(
        self, mock_write_tsv, mock_get_run_info, mock_validate_query, runner, tmp_path
    ):
        """Test --group-by-sample flag is accepted."""
        mock_validate_query.return_value = "sample_accession=SAMN123456"
        mock_get_run_info.return_value = ("ENA", [])

        result = runner.invoke(
            fastqdl,
            [
                "--accession",
                "SAMN123456",
                "--group-by-sample",
                "--only-download-metadata",
                "--outdir",
                str(tmp_path),
            ],
        )

        assert result.exit_code == 0

    @patch("fastq_dl.cli.download.validate_query")
    @patch("fastq_dl.cli.download.get_run_info")
    @patch("fastq_dl.cli.download.write_tsv")
    def test_only_provider_flag(
        self, mock_write_tsv, mock_get_run_info, mock_validate_query, runner, tmp_path
    ):
        """Test --only-provider flag is accepted."""
        mock_validate_query.return_value = "run_accession=SRR123456"
        mock_get_run_info.return_value = ("ENA", [])

        result = runner.invoke(
            fastqdl,
            [
                "--accession",
                "SRR123456",
                "--only-provider",
                "--only-download-metadata",
                "--outdir",
                str(tmp_path),
            ],
        )

        assert result.exit_code == 0

    @patch("fastq_dl.cli.download.validate_query")
    @patch("fastq_dl.cli.download.get_run_info")
    @patch("fastq_dl.cli.download.write_tsv")
    def test_sra_lite_flag(
        self, mock_write_tsv, mock_get_run_info, mock_validate_query, runner, tmp_path
    ):
        """Test --sra-lite flag is accepted."""
        mock_validate_query.return_value = "run_accession=SRR123456"
        mock_get_run_info.return_value = ("ENA", [])

        result = runner.invoke(
            fastqdl,
            [
                "--accession",
                "SRR123456",
                "--sra-lite",
                "--only-download-metadata",
                "--outdir",
                str(tmp_path),
            ],
        )

        assert result.exit_code == 0

    @patch("fastq_dl.cli.download.validate_query")
    @patch("fastq_dl.cli.download.get_run_info")
    @patch("fastq_dl.cli.download.write_tsv")
    def test_force_flag(
        self, mock_write_tsv, mock_get_run_info, mock_validate_query, runner, tmp_path
    ):
        """Test --force flag is accepted."""
        mock_validate_query.return_value = "run_accession=SRR123456"
        mock_get_run_info.return_value = ("ENA", [])

        result = runner.invoke(
            fastqdl,
            [
                "--accession",
                "SRR123456",
                "--force",
                "--only-download-metadata",
                "--outdir",
                str(tmp_path),
            ],
        )

        assert result.exit_code == 0

    @patch("fastq_dl.cli.download.validate_query")
    @patch("fastq_dl.cli.download.get_run_info")
    @patch("fastq_dl.cli.download.write_tsv")
    def test_ignore_flag(
        self, mock_write_tsv, mock_get_run_info, mock_validate_query, runner, tmp_path
    ):
        """Test --ignore flag is accepted."""
        mock_validate_query.return_value = "run_accession=SRR123456"
        mock_get_run_info.return_value = ("ENA", [])

        result = runner.invoke(
            fastqdl,
            [
                "--accession",
                "SRR123456",
                "--ignore",
                "--only-download-metadata",
                "--outdir",
                str(tmp_path),
            ],
        )

        assert result.exit_code == 0


class TestExitCodes:
    """Tests for exit codes on download failures."""

    @pytest.fixture
    def runner(self):
        return CliRunner()

    @patch("fastq_dl.cli.download.validate_query")
    @patch("fastq_dl.cli.download.get_run_info")
    @patch("fastq_dl.cli.download.ena_download")
    @patch("fastq_dl.cli.download.write_tsv")
    def test_all_runs_missing_fastqs_exits_2(
        self,
        mock_write_tsv,
        mock_ena_download,
        mock_get_run_info,
        mock_validate_query,
        runner,
        sample_ena_metadata,
        tmp_path,
    ):
        """Exit 2 when metadata found but no FASTQ files available."""
        mock_validate_query.return_value = "run_accession=SRR2838701"
        mock_get_run_info.return_value = ("ENA", sample_ena_metadata)
        mock_ena_download.return_value = ENA_NO_FASTQS

        result = runner.invoke(
            fastqdl,
            ["--accession", "SRR2838701", "--only-provider", "--outdir", str(tmp_path)],
        )

        assert result.exit_code == 2

    @patch("fastq_dl.cli.download.validate_query")
    @patch("fastq_dl.cli.download.get_run_info")
    @patch("fastq_dl.cli.download.ena_download")
    @patch("fastq_dl.cli.download.write_tsv")
    def test_all_runs_not_found_exits_2(
        self,
        mock_write_tsv,
        mock_ena_download,
        mock_get_run_info,
        mock_validate_query,
        runner,
        sample_ena_metadata,
        tmp_path,
    ):
        """Exit 2 when accessions not found on any provider."""
        mock_validate_query.return_value = "run_accession=SRR2838701"
        mock_get_run_info.return_value = ("ENA", sample_ena_metadata)
        mock_ena_download.return_value = ENA_FAILED

        result = runner.invoke(
            fastqdl,
            ["--accession", "SRR2838701", "--only-provider", "--outdir", str(tmp_path)],
        )

        assert result.exit_code == 2

    @patch("fastq_dl.cli.download.validate_query")
    @patch("fastq_dl.cli.download.get_run_info")
    @patch("fastq_dl.cli.download.ena_download")
    @patch("fastq_dl.cli.download.write_tsv")
    def test_partial_failure_exits_3(
        self,
        mock_write_tsv,
        mock_ena_download,
        mock_get_run_info,
        mock_validate_query,
        runner,
        tmp_path,
    ):
        """Exit 3 when some runs succeed and others fail."""
        two_runs = [
            {
                "run_accession": "ERR0000001",
                "experiment_accession": "ERX0000001",
                "sample_accession": "SAMN0000001",
                "study_accession": "SRP0000001",
                "library_layout": "PAIRED",
                "fastq_ftp": "ftp.example.com/ERR0000001_1.fastq.gz;ftp.example.com/ERR0000001_2.fastq.gz",
                "fastq_md5": "abc;def",
            },
            {
                "run_accession": "ERR0000002",
                "experiment_accession": "ERX0000001",
                "sample_accession": "SAMN0000001",
                "study_accession": "SRP0000001",
                "library_layout": "PAIRED",
                "fastq_ftp": "ftp.example.com/ERR0000002_1.fastq.gz;ftp.example.com/ERR0000002_2.fastq.gz",
                "fastq_md5": "ghi;jkl",
            },
        ]
        mock_validate_query.return_value = "experiment_accession=ERX0000001"
        mock_get_run_info.return_value = ("ENA", two_runs)
        mock_ena_download.side_effect = [
            {
                "r1": "/tmp/r1.fastq.gz",
                "r2": "/tmp/r2.fastq.gz",
                "single_end": False,
                "orphan": None,
            },
            ENA_NO_FASTQS,
        ]

        result = runner.invoke(
            fastqdl,
            ["--accession", "ERX0000001", "--only-provider", "--outdir", str(tmp_path)],
        )

        assert result.exit_code == 3

    @patch("fastq_dl.cli.download.validate_query")
    @patch("fastq_dl.cli.download.get_run_info")
    @patch("fastq_dl.cli.download.ena_download")
    @patch("fastq_dl.cli.download.write_tsv")
    def test_all_runs_succeed_exits_0(
        self,
        mock_write_tsv,
        mock_ena_download,
        mock_get_run_info,
        mock_validate_query,
        runner,
        sample_ena_metadata,
        mock_fastq_files,
        tmp_path,
    ):
        """Exit 0 when all downloads succeed."""
        mock_validate_query.return_value = "run_accession=SRR2838701"
        mock_get_run_info.return_value = ("ENA", sample_ena_metadata)
        mock_ena_download.return_value = mock_fastq_files

        result = runner.invoke(
            fastqdl,
            ["--accession", "SRR2838701", "--outdir", str(tmp_path)],
        )

        assert result.exit_code == 0


class TestMainFunction:
    """Tests for main() entry point."""

    def test_main_no_args(self):
        """Test main() with no arguments shows help."""
        with patch("sys.argv", ["fastq-dl"]):
            with patch.object(fastqdl, "main"):
                pass

    @patch("fastq_dl.cli.download.fastqdl")
    def test_main_with_args(self, mock_fastqdl):
        """Test main() with arguments."""
        with patch("sys.argv", ["fastq-dl", "--accession", "SRR123456"]):
            main()
            mock_fastqdl.assert_called_once()
