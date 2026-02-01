# CLAUDE.md - fastq-dl Development Guide

## Project Overview

**fastq-dl** is a Python CLI tool for downloading FASTQ sequencing files from the European Nucleotide Archive (ENA) or Sequence Read Archive (SRA). It accepts various accession types (Project, Study, BioSample, Sample, Experiment, Run) and handles provider fallback, retry logic, and optional run merging.

## Quick Reference

```bash
# Install dependencies
poetry install

# Run tests (unit only, excludes integration)
just test-cov

# Run integration tests (makes real API calls)
just test-integration

# Format code
just fmt

# Lint code
just lint

# Full check (format + lint)
just check

# Build package
just build
```

## Project Structure

```
fastq-dl/
├── fastq_dl/                    # Main package
│   ├── cli/
│   │   └── download.py          # CLI entry point (click-based)
│   ├── providers/
│   │   ├── generic.py           # Provider coordination and fallback
│   │   ├── ena.py               # ENA FTP downloads via wget
│   │   └── sra.py               # SRA downloads via sra-tools
│   ├── constants.py             # Shared constants (URLs, suffixes)
│   ├── exceptions.py            # Custom exception hierarchy
│   └── utils.py                 # Utilities (execute, md5sum, merge_runs)
├── tests/                       # Test suite
│   ├── conftest.py              # Shared fixtures
│   ├── test_cli.py              # CLI tests
│   ├── test_providers_*.py      # Provider-specific tests
│   ├── test_utils.py            # Utility function tests
│   └── test_integration.py      # Real API integration tests
├── pyproject.toml               # Poetry config, dependencies, pytest settings
├── justfile                     # Task runner commands
└── environment.yml              # Conda environment (for external tools)
```

## Architecture

### Core Flow
1. **CLI** (`cli/download.py`) → parses args, sets up logging, orchestrates workflow
2. **Validation** (`utils.validate_query`) → validates accession format via regex
3. **Metadata Query** (`providers/generic.get_run_info`) → queries ENA first, falls back to SRA
4. **Download** (`_download_with_fallback`) → attempts primary provider, falls back if needed
5. **Post-processing** → optional run merging, writes TSV metadata files

### Provider Pattern
- **ENA**: Uses HTTP API for metadata, wget for FTP downloads, MD5 validation
- **SRA**: Uses pysradb for metadata, prefetch + fasterq-dump for downloads, pigz for compression
- Automatic fallback: ENA → SRA (or vice versa based on `--provider` flag)

### Exception Hierarchy
```
FastqDLError (base)
├── ValidationError     # Invalid accession format
├── ProviderError       # API/query failures (has provider, status_code)
└── DownloadError       # File download failures (has accession, provider)
```

## Key Modules

### cli/download.py
- `fastqdl()` - Click command entry point
- `_run_download()` - Core download logic
- `_download_with_fallback()` - Provider failover handling

### providers/ena.py
- `get_ena_metadata(query)` - Queries ENA Data Warehouse API
- `ena_download(run, outdir, ...)` - Orchestrates ENA download
- `download_ena_fastq(ftp, outdir, md5, ...)` - Downloads single FASTQ via wget

### providers/sra.py
- `get_sra_metadata(query)` - Queries SRA via pysradb
- `sra_download(accession, outdir, ...)` - Uses prefetch + fasterq-dump + pigz

### utils.py
- `execute(cmd, ...)` - Subprocess wrapper with retries
- `md5sum(fastq)` - Calculates MD5 in 10MB chunks
- `merge_runs(runs, output)` - Concatenates FASTQs
- `validate_query(query)` - Regex validation for accession types
- `write_tsv(data, output)` - Writes run-info and run-mergers TSV files

## External Tool Dependencies

The following CLI tools are required at runtime:
- `wget` - FTP downloads (ENA)
- `prefetch` - SRA file download (from sra-tools)
- `fasterq-dump` - SRA to FASTQ conversion (from sra-tools)
- `vdb-config` - SRA configuration (from sra-tools)
- `pigz` - Parallel gzip compression

## Testing

### Test Categories
- **Unit tests**: Fast, mocked, run by default
- **Integration tests**: Real API calls, marked with `@pytest.mark.integration`

### Running Tests
```bash
# Unit tests with coverage (requires 70% minimum)
just test-cov

# Integration tests only
just test-integration

# Specific test file
just test tests/test_utils.py

# Specific test function
just test tests/test_utils.py::TestValidateQuery::test_valid_run_accession
```

### Writing Tests
- Place fixtures in `tests/conftest.py`
- Use `@pytest.mark.integration` for tests making real API calls
- Mock external calls with `responses` library (HTTP) or `unittest.mock` (subprocess)
- Key fixtures: `tmp_outdir`, `sample_ena_metadata`, `mock_execute_success`

## Code Style

- **Formatter**: Black
- **Import sorting**: isort (profile: black)
- **Linter**: flake8
- **Type hints**: Used in function signatures
- **Docstrings**: NumPy style with Args/Returns/Raises

Run formatting:
```bash
just fmt        # Apply formatting
just check-fmt  # Check without changing
```

## Common Development Tasks

### Adding a new CLI option
1. Add `@click.option()` decorator in `cli/download.py:fastqdl()`
2. Add to appropriate option group in `click.rich_click.OPTION_GROUPS`
3. Pass through to `_run_download()` and relevant functions

### Adding a new provider
1. Create `providers/newprovider.py` with `get_*_metadata()` and `*_download()` functions
2. Add constants to `constants.py`
3. Integrate into `providers/generic.py` for fallback logic
4. Update `_download_with_fallback()` in `cli/download.py`

### Modifying accession validation
- Edit regex patterns in `utils.validate_query()`
- See ENA accession docs: https://ena-docs.readthedocs.io/en/latest/submit/general-guide/accessions.html

## Important Patterns

### Command Execution
Always use list form for subprocess commands (security):
```python
# Good
execute(["wget", "--quiet", "-O", str(fastq), f"ftp://{ftp}"])

# Avoid (string splitting with shlex)
execute(f"wget --quiet -O {fastq} ftp://{ftp}")
```

### Return Value Pattern
Provider functions return sentinel values on failure:
```python
from fastq_dl.constants import ENA_FAILED, SRA_FAILED

result = ena_download(...)
if result == ENA_FAILED:
    # Handle failure
```

### Retry Pattern
```python
outcome = execute(
    cmd,
    max_attempts=max_attempts,
    sleep=sleep,
    is_sra=True,  # Enables SRA-specific error handling
)
```

## CI/CD

### GitHub Actions Workflows
- **fastq-dl.yml**: Runs on push/PR, tests Python 3.10-3.13, uploads coverage
- **release.yml**: Triggered by version tags (v*.*.*), publishes to PyPI

### Coverage Requirements
- Minimum: 70%
- Uploaded to Codecov (Python 3.10 only)
- Excludes: `__init__.py`, `if __name__ == "__main__":`

## Version Management

```bash
# Check current version
poetry version -s

# Bump version
poetry version patch  # or minor, major

# Tag release
just tag  # Prints commands to run
```

## Debugging Tips

- Use `--verbose` flag for debug logging
- Check `--only-download-metadata` to test metadata queries without downloading
- Use `--only-provider` to isolate provider-specific issues
- SRA exit code 3 = accession not found (triggers fallback)

## Known Considerations

1. **ENA API**: Returns TSV format, queries use `fields=all`
2. **SRA Lite**: Compressed quality scores, enabled via `--sra-lite` flag
3. **Paired-end detection**: ENA uses `library_layout` field; SRA uses `fasterq-dump --split-3`
4. **Orphan reads**: SRA can produce orphan reads file from paired-end data
5. **MD5 validation**: Enabled by default, skip with `--ignore` flag
