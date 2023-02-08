# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.3.0]

### Added

- Shorthand option `-a` for `--max_attempts` and `v` for `--verbose`
- If metadata is unavailable from ENA, try fetching metadata from SRA using [`pysradb`][pysradb]
- Provided query is validated against accession regular expressions

[1.3.0]: https://github.com/rpetit3/fastq-dl/compare/v1.2.0...v1.3.0
[pysradb]: https://github.com/saketkc/pysradb