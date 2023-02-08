# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [2.0.0]

### Added

- Shorthand option `-a` for `--max_attempts` and `v` for `--verbose`
- If metadata is unavailable from ENA, try fetching metadata from SRA using [`pysradb`][pysradb]
- Provided query is validated against accession regular expressions
- Support for Biosample/Sample accessions
- Rich click and logging support
- type hints
- Support for Gitpod

### Removed

- Deprecated `--sra_only` parameter is now removed
- `--ftp_only` no longer required without Aspera support

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
