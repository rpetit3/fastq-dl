import pytest

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


def test_validate_query_invalid(caplog):
    with pytest.raises(SystemExit):
        validate_query("INVALID123")
    assert "is not a Study, Sample, Experiment, or Run accession" in caplog.text
