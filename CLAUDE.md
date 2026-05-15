# CLAUDE.md - fastq-dl Development Guide

## Project Overview

**fastq-dl** is a Python CLI tool for downloading FASTQ sequencing files from the European Nucleotide Archive (ENA) or Sequence Read Archive (SRA). It accepts various accession types (Project, Study, BioSample, Sample, Experiment, Run) and handles provider fallback, retry logic, and optional run merging.

- **Version**: 4.0.0
- **License**: MIT
- **Python**: >=3.10, <3.14
- **Repository**: https://github.com/rpetit3/fastq-dl

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
│   ├── __init__.py              # Version from importlib.metadata
│   ├── cli/
│   │   ├── __init__.py
│   │   └── download.py          # CLI entry point (click-based)
│   ├── providers/
│   │   ├── __init__.py
│   │   ├── generic.py           # Provider coordination and fallback
│   │   ├── ena.py               # ENA downloads via wget (FTP or HTTPS)
│   │   └── sra.py               # SRA downloads via sracha
│   ├── constants.py             # Shared constants (URLs, suffixes, sentinels)
│   ├── exceptions.py            # Custom exception hierarchy (8 types)
│   └── utils.py                 # Utilities (execute, md5sum, merge_runs, write_tsv)
├── tests/                       # Test suite
│   ├── conftest.py              # Shared fixtures
│   ├── test_cli.py              # CLI option/argument tests
│   ├── test_download.py         # Download logic tests
│   ├── test_integration.py      # Real API integration tests
│   ├── test_providers_ena.py    # ENA provider tests
│   ├── test_providers_generic.py # Generic provider tests
│   ├── test_providers_sra.py    # SRA provider tests
│   └── test_utils.py            # Utility function tests
├── test/                        # Test data files (sample FASTQs, TSVs)
├── .claude/
│   └── skills/
│       └── update-catalog/
│           ├── skill.md             # Skill definition for /update-catalog
│           └── scripts/
│               └── update_catalog.py # Regenerates catalog.json and llms.txt
├── .github/
│   ├── workflows/
│   │   ├── fastq-dl.yml         # CI: unit tests, integration tests, failure notification
│   │   └── release.yml          # CD: publish to PyPI on version tags
│   └── FUNDING.yml
├── pyproject.toml               # Poetry config, dependencies, pytest/ruff/coverage settings
├── justfile                     # Task runner commands
├── poetry.toml                  # Poetry local config
├── poetry.lock                  # Locked dependencies
├── environment.yml              # Conda environment (for external tools)
├── catalog.json                 # Machine-readable project metadata
├── llms.txt                     # AI-discovery document
├── citation.cff                 # Citation metadata
├── CHANGELOG.md                 # Release history
├── README.md                    # User documentation
├── LICENSE                      # MIT license
├── codecov.yml                  # Codecov configuration
├── .pre-commit-config.yaml      # Pre-commit hooks
└── .gitignore
```

## Architecture

### Core Flow
1. **CLI** (`cli/download.py`) → parses args via Click, sets up logging, orchestrates workflow
2. **Validation** (`utils.validate_query`) → validates accession format via regex
3. **Metadata Query** (`providers/generic.get_run_info`) → queries primary provider, falls back to secondary
4. **Download** (`cli/download.py:_download_with_fallback`) → attempts primary provider, falls back if needed
5. **Post-processing** → optional run merging, writes TSV metadata files

### Provider Pattern
- **ENA**: Uses HTTP API for metadata, wget for FTP/HTTPS downloads, MD5 validation
- **SRA**: Uses pysradb for metadata, sracha for download + FASTQ conversion + compression
- Automatic fallback between providers (configurable with `--provider` and `--only-provider`)

### Entry Point
- `pyproject.toml` defines: `fastq-dl = "fastq_dl.cli.download:main"`
- `main()` calls `fastqdl(["--help"])` when no args given, otherwise `fastqdl()`
- `fastqdl()` is the Click command that wraps `_run_download()` with exception handling

### Exception Hierarchy
```
FastqDLError (base)
├── ValidationError          # Invalid accession format
├── ProviderError            # API/query failures (attrs: provider, status_code)
├── EmptyResultError         # Provider returned HTTP 200 with no data (attrs: provider)
├── DownloadError            # File download failures (attrs: accession, provider)
├── AccessionNotFoundError   # Run not found on any provider (attrs: failed_runs)
├── MissingFastqsError       # FASTQs not yet available/synced (attrs: failed_runs)
└── PartialDownloadError     # Some runs succeeded, others failed (attrs: failed_runs, successful_runs)
```

Exit codes: `1` = validation/provider/download error, `2` = empty result/not found/missing FASTQs, `3` = partial download

## Key Modules

### cli/download.py
- `fastqdl()` — Click command entry point with all CLI options
- `_run_download()` — Core download logic: validate, query metadata, download each run, merge if grouped, write TSV output
- `_download_with_fallback()` — Tries primary provider, falls back to alternate on failure. Returns `(fastqs_dict | None, error_str | None)`
- `main()` — Entry point wrapper; shows help when no args given

### providers/generic.py
- `get_run_info(accession, query, provider, only_provider, ...)` — Provider coordination with dual attempt counters (`primary_attempt`, `secondary_attempt`). Tries primary provider up to `max_attempts`, then falls back to secondary.
- `_is_ena_empty_response(ena_data)` — Detects ENA HTTP 200 with no data rows

### providers/ena.py
- `get_ena_metadata(query)` — Queries ENA Data Warehouse API, returns `[success: bool, data]`
- `ena_download(run, outdir, ...)` — Orchestrates ENA download for a single run
- `download_ena_fastq(ftp, outdir, md5, ...)` — Downloads a single FASTQ via wget (FTP or HTTPS protocol)

### providers/sra.py
- `get_sra_metadata(query)` — Queries SRA via pysradb, returns `[success: bool, data]`
- `sra_download(accession, outdir, ...)` — Uses sracha get for download, conversion, and compression

### utils.py
- `execute(cmd, ...)` — Subprocess wrapper with retries and SRA-specific error handling (exit code 3 = not found)
- `md5sum(fastq)` — Calculates MD5 in 10MB chunks
- `merge_runs(runs, output)` — Concatenates FASTQs, validates all files exist first
- `validate_query(query)` — Regex validation for accession types (Project, Study, BioSample, Sample, Experiment, Run)
- `write_tsv(data, output, na_value="")` — Writes TSV files; handles both run-mergers format (dict of `{accession: {r1, r2}}`) and run-info format (list of dicts with dynamic fieldname union via `_all_fieldnames()`)
- `_all_fieldnames(rows)` — Collects union of all keys across rows preserving first-seen order

### constants.py
```python
# Provider names and sentinel values
ENA = "ENA"
ENA_FAILED = "ENA_NOT_FOUND"
ENA_NO_FASTQS = "ENA_NO_FASTQS"
SRA = "SRA"
SRA_FAILED = "SRA_NOT_FOUND"
SRA_DOWNLOAD_FAILED = "SRA_DOWNLOAD_FAILED"

# File suffixes
RUN_INFO_SUFFIX = "-run-info.tsv"
RUN_MERGERS_SUFFIX = "-run-mergers.tsv"

# FASTQ suffixes (compressed)
PE_R1_SUFFIX = "_1.fastq.gz"
PE_R2_SUFFIX = "_2.fastq.gz"
SE_SUFFIX = ".fastq.gz"

# FASTQ suffixes (uncompressed, for --skip-compression)
PE_R1_SUFFIX_UNCOMPRESSED = "_1.fastq"
PE_R2_SUFFIX_UNCOMPRESSED = "_2.fastq"
SE_SUFFIX_UNCOMPRESSED = ".fastq"

# Merged FASTQ suffixes
MERGED_R1_SUFFIX = "_R1.fastq.gz"
MERGED_R2_SUFFIX = "_R2.fastq.gz"
```

## CLI Options

Organized into groups via `click.rich_click.OPTION_GROUPS`:

### Required Options
| Option | Short | Description |
|--------|-------|-------------|
| `--accession` | `-a` | ENA/SRA accession to query |

### Provider Options
| Option | Default | Description |
|--------|---------|-------------|
| `--provider` | `ena` | Provider to use (ena or sra) |
| `--protocol` | `ftp` | Protocol for ENA downloads (ftp or https) |
| `--sra-lite` | off | Use SRA Lite (compressed quality scores) |
| `--skip-compression` | off | Skip compression of SRA downloads |
| `--gzip-level` | 1 | Gzip compression level for SRA downloads (1-9) |

### Download Options
| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--max-attempts` | `-m` | 3 | Maximum download attempts |
| `--only-provider` | | off | Only attempt specified provider (no fallback) |
| `--only-download-metadata` | | off | Skip downloads, retrieve metadata only |
| `--group-by-experiment` | | off | Group and merge runs by experiment accession |
| `--group-by-sample` | | off | Group and merge runs by sample accession |
| `--ignore` | `-I` | off | Skip MD5 validation (ENA) or relax integrity checks (SRA) |

### Additional Options
| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--outdir` | `-o` | `./` | Output directory |
| `--prefix` | | `fastq` | Prefix for log/metadata files |
| `--cpus` | | 4 | CPUs for SRA downloads |
| `--force` | `-F` | off | Overwrite existing files |
| `--silent` | | off | Only critical errors |
| `--sleep` | `-s` | 10 | Seconds between retries |
| `--verbose` | `-v` | off | Debug logging |

## External Tool Dependencies

Required at runtime (install via conda/environment.yml):
- `wget` — FTP/HTTPS downloads (ENA)
- `sracha` — SRA download, FASTQ conversion, and compression (from sracha-rs)

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
- Key fixtures: `tmp_outdir`, `sample_ena_metadata`, `mock_execute_success`, `mock_fastq_files`, `sample_sra_metadata`, `single_end_ena_metadata`

## Code Style

- **Formatter/Linter**: Ruff (line-length 88, target py310, selects E/W/F/I, ignores E501)
- **Type hints**: Used in function signatures
- **Docstrings**: NumPy style with Args/Returns/Raises

```bash
just fmt        # Apply formatting and auto-fix lint issues
just check-fmt  # Check formatting without changing
just lint       # Check for lint violations
just check      # Full check (format + lint)
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
Provider download functions return either a dict (success) or a sentinel string (failure):
```python
from fastq_dl.constants import ENA_FAILED, ENA_NO_FASTQS, SRA_FAILED, SRA_DOWNLOAD_FAILED

result = ena_download(...)
if isinstance(result, str) and result in {ENA_FAILED, ENA_NO_FASTQS}:
    # Handle failure

# Success dict structure:
# {"r1": str, "r2": str, "single_end": bool, "orphan": str | None}
```

### Retry Pattern
```python
outcome = execute(
    cmd,
    max_attempts=max_attempts,
    sleep=sleep,
    is_sra=True,  # Enables SRA-specific error handling (exit code 3 = not found)
)
```

## CI/CD

### GitHub Actions Workflows

**fastq-dl.yml** — 3-job pipeline:
1. **unit-tests**: Runs on push/PR to main/master/dev. Tests Python 3.10-3.13, checks formatting, linting, runs unit tests with coverage.
2. **integration-tests**: Runs on push to main/master, weekly schedule (Sundays 2am UTC), or manual trigger. Tests real ENA/SRA downloads with various accession types.
3. **notify-failure**: Auto-creates GitHub issue on scheduled run failures (labeled `ci-failure`, `automated`).

**release.yml** — Triggered by version tags (`v*.*.*`), publishes to PyPI.

### Coverage Requirements
- Minimum: 70%
- Uploaded to Codecov (Python 3.10 only)
- Branch coverage enabled
- Omits: `fastq_dl/__init__.py`

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
3. **Paired-end detection**: ENA uses `library_layout` field; SRA uses sracha `--split split-3` (default)
4. **Orphan reads**: SRA can produce orphan reads file from paired-end data
5. **MD5 validation**: ENA validates by default (skip with `--ignore`); SRA (sracha) always verifies integrity (`--ignore` maps to `--no-strict`)
6. **Protocol choice**: ENA supports both FTP (default) and HTTPS via `--protocol`
7. **Compression**: SRA downloads compressed by default via sracha (skip with `--skip-compression`, tune with `--gzip-level`)

## Catalog and LLM Discovery

- `catalog.json` — Machine-readable project metadata (modules, functions, CLI options, dependencies)
- `llms.txt` — AI-discovery document with architecture overview and module index
- Regenerate both with: `python .claude/skills/update-catalog/scripts/update_catalog.py`
- Or use the `/update-catalog` skill in Claude Code
