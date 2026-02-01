import pytest

from fastq_dl.exceptions import ValidationError
from fastq_dl.utils import md5sum, merge_runs, validate_query


@pytest.fixture
def test_file(tmp_path):
    # Create a temporary file with known content
    test_file = tmp_path / "test.fastq"
    test_content = b"@read1\nACGT\n+\n1234\n"
    with open(test_file, "wb") as f:
        f.write(test_content)
    return test_file


@pytest.fixture
def test_files(tmp_path):
    # Create temporary files with known content
    file1 = tmp_path / "file1.fastq"
    file2 = tmp_path / "file2.fastq"
    file1_content = b"@read1\nACGT\n+\n1234\n"
    file2_content = b"@read2\nTGCA\n+\n4321\n"
    with open(file1, "wb") as f:
        f.write(file1_content)
    with open(file2, "wb") as f:
        f.write(file2_content)
    return [file1, file2]


def test_md5sum_valid_file(test_file):
    # Expected MD5 checksum for the known content
    expected_md5 = "428f145dbcbe924a05f49547d29f19fc"
    # Calculate the MD5 checksum using the md5sum function
    calculated_md5 = md5sum(test_file)
    # Compare the expected checksum with the calculated checksum
    assert calculated_md5 == expected_md5


def test_md5sum_nonexistent_file():
    # Test the case where the file does not exist
    assert md5sum("nonexistent.fastq") is None


def test_merge_runs_multiple_files(test_files, tmp_path):
    # Output file path
    output_file = str(tmp_path / "merged.fastq")
    # Merge the files
    merge_runs(test_files, output_file)
    # Expected content of the merged file
    expected_content = b"@read1\nACGT\n+\n1234\n@read2\nTGCA\n+\n4321\n"
    # Verify the content of the merged file
    with open(output_file, "rb") as f:
        assert f.read() == expected_content


def test_merge_runs_single_file(test_files, tmp_path):
    # Output file path
    output_file = tmp_path / "merged.fastq"
    # Use only one file for merging
    merge_runs([test_files[0]], output_file)
    # Expected content of the merged file
    expected_content = b"@read1\nACGT\n+\n1234\n"
    # Verify the content of the merged file
    with open(output_file, "rb") as f:
        assert f.read() == expected_content


def test_validate_query_project_study():
    assert (
        validate_query("PRJNA123456")
        == "(study_accession=PRJNA123456 OR secondary_study_accession=PRJNA123456)"
    )
    assert (
        validate_query("ERP123456")
        == "(study_accession=ERP123456 OR secondary_study_accession=ERP123456)"
    )


def test_validate_query_sample_biosample():
    assert (
        validate_query("SAMN12345678")
        == "(sample_accession=SAMN12345678 OR secondary_sample_accession=SAMN12345678)"
    )
    assert (
        validate_query("SRS123456")
        == "(sample_accession=SRS123456 OR secondary_sample_accession=SRS123456)"
    )


def test_validate_query_experiment():
    assert validate_query("SRX123456") == "experiment_accession=SRX123456"


def test_validate_query_run():
    assert validate_query("SRR123456") == "run_accession=SRR123456"
    assert validate_query("ERR123456") == "run_accession=ERR123456"
    assert validate_query("DRR123456") == "run_accession=DRR123456"


def test_validate_query_invalid():
    with pytest.raises(ValidationError) as exc_info:
        validate_query("INVALID123")
    assert "is not a Study, Sample, Experiment, or Run accession" in str(exc_info.value)


# ============================================================================
# Tests for execute() function
# ============================================================================

from unittest.mock import patch, MagicMock
import csv
import subprocess

from fastq_dl.utils import execute, write_tsv
from fastq_dl.constants import ENA_FAILED, SRA_FAILED


class TestExecute:
    """Tests for execute function."""

    @patch("fastq_dl.utils.subprocess.run")
    def test_successful_execution(self, mock_run):
        """Test successful command execution."""
        mock_run.return_value = MagicMock(
            returncode=0,
            stdout="output",
            stderr="",
        )

        result = execute("echo test")

        assert result == 0

    @patch("fastq_dl.utils.subprocess.run")
    def test_capture_stdout(self, mock_run):
        """Test capturing stdout."""
        mock_run.return_value = MagicMock(
            returncode=0,
            stdout="captured output",
            stderr="",
        )

        result = execute("echo test", capture_stdout=True)

        assert result == "captured output"

    @patch("fastq_dl.utils.subprocess.run")
    def test_command_failure_returns_ena_failed(self, mock_run):
        """Test command failure returns ENA_FAILED by default."""
        mock_run.side_effect = subprocess.CalledProcessError(
            returncode=1,
            cmd="failing command",
            stderr="error",
        )

        with patch("time.sleep"):
            result = execute("failing command", max_attempts=1, sleep=0)

        assert result == ENA_FAILED

    @patch("fastq_dl.utils.subprocess.run")
    def test_sra_exit_code_3_returns_sra_failed(self, mock_run):
        """Test SRA-specific exit code 3 handling."""
        mock_run.side_effect = subprocess.CalledProcessError(
            returncode=3,
            cmd="prefetch SRR123456",
            stderr="Item not found\nMore info",
        )

        result = execute("prefetch SRR123456", is_sra=True)

        assert result == SRA_FAILED

    @patch("fastq_dl.utils.subprocess.run")
    def test_sra_failure_returns_sra_failed(self, mock_run):
        """Test SRA command failure returns SRA_FAILED."""
        mock_run.side_effect = subprocess.CalledProcessError(
            returncode=1,
            cmd="sra command",
            stderr="error",
        )

        with patch("time.sleep"):
            result = execute("sra command", is_sra=True, max_attempts=1, sleep=0)

        assert result == SRA_FAILED

    @patch("fastq_dl.utils.subprocess.run")
    @patch("time.sleep")
    def test_retry_on_failure(self, mock_sleep, mock_run):
        """Test retry logic on command failure."""
        mock_run.side_effect = subprocess.CalledProcessError(
            returncode=1,
            cmd="failing command",
            stderr="error",
        )

        result = execute("failing command", max_attempts=3, sleep=1)

        # Should have slept between retries
        assert mock_sleep.call_count == 2
        assert result == ENA_FAILED

    @patch("fastq_dl.utils.subprocess.run")
    def test_directory_parameter(self, mock_run):
        """Test directory parameter is passed correctly."""
        mock_run.return_value = MagicMock(
            returncode=0,
            stdout="",
            stderr="",
        )

        execute("ls", directory="/tmp")

        mock_run.assert_called_once()
        call_kwargs = mock_run.call_args[1]
        assert call_kwargs["cwd"] == "/tmp"

    @patch("fastq_dl.utils.subprocess.run")
    def test_stdout_file_written(self, mock_run, tmp_path):
        """Test stdout is written to file when specified."""
        mock_run.return_value = MagicMock(
            returncode=0,
            stdout="captured stdout content",
            stderr="",
        )

        stdout_file = tmp_path / "stdout.txt"
        execute("echo test", stdout_file=str(stdout_file))

        assert stdout_file.exists()
        assert stdout_file.read_text() == "captured stdout content"

    @patch("fastq_dl.utils.subprocess.run")
    def test_stderr_file_written(self, mock_run, tmp_path):
        """Test stderr is written to file when specified."""
        mock_run.return_value = MagicMock(
            returncode=0,
            stdout="",
            stderr="captured stderr content",
        )

        stderr_file = tmp_path / "stderr.txt"
        execute("echo test", stderr_file=str(stderr_file))

        assert stderr_file.exists()
        assert stderr_file.read_text() == "captured stderr content"

    @patch("fastq_dl.utils.subprocess.run")
    def test_stderr_file_written_on_failure(self, mock_run, tmp_path):
        """Test stderr is written to file even when command fails."""
        mock_run.side_effect = subprocess.CalledProcessError(
            returncode=1,
            cmd="failing command",
            stderr="error output",
        )

        stderr_file = tmp_path / "stderr.txt"
        with patch("time.sleep"):
            execute("failing command", stderr_file=str(stderr_file), max_attempts=1)

        assert stderr_file.exists()
        assert stderr_file.read_text() == "error output"

    @patch("fastq_dl.utils.subprocess.run")
    def test_stdout_file_not_written_when_empty(self, mock_run, tmp_path):
        """Test stdout file is not written when stdout is empty."""
        mock_run.return_value = MagicMock(
            returncode=0,
            stdout="",
            stderr="",
        )

        stdout_file = tmp_path / "stdout.txt"
        execute("echo test", stdout_file=str(stdout_file))

        # File should not be created when stdout is empty
        assert not stdout_file.exists()


# ============================================================================
# Tests for write_tsv() function
# ============================================================================


class TestWriteTsv:
    """Tests for write_tsv function."""

    def test_write_run_info(self, tmp_path):
        """Test writing run info TSV."""
        data = [
            {"run_accession": "SRR123456", "sample_accession": "SAMN123456"},
            {"run_accession": "SRR789012", "sample_accession": "SAMN789012"},
        ]
        output = tmp_path / "test-run-info.tsv"

        write_tsv(data, str(output))

        assert output.exists()
        with open(output) as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)
            assert len(rows) == 2
            assert rows[0]["run_accession"] == "SRR123456"
            assert rows[1]["run_accession"] == "SRR789012"

    def test_write_run_mergers(self, tmp_path):
        """Test writing run mergers TSV."""
        data = {
            "SRX123456": {
                "r1": ["file1_R1.fastq.gz", "file2_R1.fastq.gz"],
                "r2": ["file1_R2.fastq.gz", "file2_R2.fastq.gz"],
            }
        }
        output = tmp_path / "test-run-mergers.tsv"

        write_tsv(data, str(output))

        assert output.exists()
        with open(output) as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)
            assert len(rows) == 1
            assert rows[0]["accession"] == "SRX123456"
            assert ";" in rows[0]["r1"]
            assert "file1_R1.fastq.gz" in rows[0]["r1"]
            assert "file2_R1.fastq.gz" in rows[0]["r1"]

    def test_write_run_mergers_multiple_accessions(self, tmp_path):
        """Test writing run mergers TSV with multiple accessions."""
        data = {
            "SRX123456": {
                "r1": ["file1_R1.fastq.gz"],
                "r2": ["file1_R2.fastq.gz"],
            },
            "SRX123457": {
                "r1": ["file2_R1.fastq.gz"],
                "r2": ["file2_R2.fastq.gz"],
            },
        }
        output = tmp_path / "test-run-mergers.tsv"

        write_tsv(data, str(output))

        assert output.exists()
        with open(output) as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)
            assert len(rows) == 2

    def test_write_empty_run_info(self, tmp_path):
        """Test writing empty run info TSV."""
        # This would fail with IndexError if data is empty
        # Testing that the function handles at least one row
        data = [{"run_accession": "SRR123456"}]
        output = tmp_path / "test-run-info.tsv"

        write_tsv(data, str(output))

        assert output.exists()

    def test_write_tsv_empty_data(self, tmp_path):
        """Test writing TSV with empty data list does not crash."""
        data = []
        output = tmp_path / "test-run-info.tsv"

        # Should not raise IndexError
        write_tsv(data, str(output))

        # File should exist but be empty
        assert output.exists()
        assert output.read_text() == ""
