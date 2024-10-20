# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### TODO

- consider refactoring more
- consider adding back aspera support (now available from bioconda)

## [3.0.0]

### Added

- `--ignore` to ignore MD5 checksums

### Changed

- (**BREAKING**) the `--force` flag now overwrites existing files, but will check the MD5
  checksums of the downloaded files. Previously, it would force download the file, but
  ignore the MD5 as well. See https://github.com/rpetit3/fastq-dl/issues/21#issuecomment-2190686184
  for more details. Use the newly added `--ignore` flag to ignore MD5 checksums too.

[3.0.0]: https://github.com/rpetit3/fastq-dl/compare/v2.0.4...v3.0.0

## [2.0.4]

### Fixed

- mismatch in `--version` output

[2.0.4]: https://github.com/rpetit3/fastq-dl/compare/v2.0.3...v2.0.4

## [2.0.3]

### Added

- Explicit setting of SRA Lite or SRA Normalized (default) downloads
- `--sra-lite` to prefer SRA Lite downloads from SRA
- print usage on empty options

[2.0.3]: https://github.com/rpetit3/fastq-dl/compare/v2.0.2...v2.0.3

## [2.0.2]

### Added

- `--only-download-metadata` to skip FASTQ downloads, and retrieve only the associated metadata @alienzj
- fallback on SRA for metadata retrieval @gtonkinhill
- include `prefetch`for SRA downloads @gtonkinhill
- added `--force` to overwrite existing downloads @gtonkinhill
- update github actions

[2.0.2]: https://github.com/rpetit3/fastq-dl/compare/v2.0.1...v2.0.2

## [2.0.1]

### Fixed

- Invalid fieldName(s) in ENA query @mbhall88
- gitpod dependency and setup
- dependabot PR

[2.0.1]: https://github.com/rpetit3/fastq-dl/compare/v2.0.0...v2.0.1

## [2.0.0]

### Added

- Shorthand option `-m` for `--max-attempts` and `-v` for `--verbose`
- If metadata is unavailable from ENA, try fetching metadata from SRA using [`pysradb`][pysradb]
- Provided query is validated against accession regular expressions
- Support for Biosample/Sample accessions
- Rich click, logging support, type hints
- Support for Gitpod
- Query must now be passed as `--accession/-a`
- Packaging with `poetry`
- README improvements

### Removed

- Deprecated `--sra_only` parameter is now removed
- `--ftp_only` no longer required without Aspera support
- Non-functioning `--debug` option
- Uneeded logging levels

[2.0.0]: https://github.com/rpetit3/fastq-dl/compare/v1.2.0...v2.0.0
[pysradb]: https://github.com/saketkc/pysradb

## [1.2.0]

### Fixed

- Use accession in pigz command

### Removed

- Support for Aspera downloads
- Manual vdb-config (no longer required in sra-tools>=3.0.1)

[1.2.0]: https://github.com/rpetit3/fastq-dl/compare/v1.1.1...v1.2.0

## [1.1.1]

### Fixed

- Handle experiment accessions with duplicate run accessions
    - Example: https://www.ebi.ac.uk/ena/browser/view/ERX012546?show=reads

[1.1.1]: https://github.com/rpetit3/fastq-dl/compare/v1.1.0...v1.1.1

## [1.1.0]

Big thanks to @mbhall88 for submitting https://github.com/rpetit3/fastq-dl/pull/5 with the following improvements!

## Added
- `-o` shorthand option for `--outdir`
- a flag `-F/--only-provider`, which supercedes `--sra_only`. I left `--sra_only` in there for backwards compatibility and added a deprecation notice in the help description for it.
- Conda environment file
- Dockerfile. Feel free to try out an image from my quay.io repo [here](https://quay.io/repository/mbhall88/fastq-dl?tab=tags).

## Changed
- Move everything into a `main` function to avoid [variable shadowing](https://en.wikipedia.org/wiki/Variable_shadowing)
- Provider is now optional and defaults to `ena`. I've also made the option case insensitive, rather than listing the cased versions of the available providers
- If ENA download fails, try SRA
- Reduced a bunch of `execute` calls which were mostly operations easily dealt with by `pathlib.Path`
- Changed `md5sum` to use `hashlib.md5` instead of executing a subprocess call to `md5sum`
- Moved a bunch of `import` statements out of function bodies to the top of the file

## Removed
- Docstring at top of file as this is a mirror of argparse's help menu and creates needless maintenance

[1.1.0]: https://github.com/rpetit3/fastq-dl/compare/v1.0.6...v1.1.0
