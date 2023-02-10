[![GitHub release (latest by date)](https://img.shields.io/github/v/release/rpetit3/fastq-dl)](https://github.com/bactopia/rpetit3/fastq-dl)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/fastq-dl/badges/downloads.svg)](https://anaconda.org/bioconda/fastq-dl)
[![Gitpod ready-to-code](https://img.shields.io/badge/Gitpod-ready--to--code-908a85?logo=gitpod)](https://gitpod.io/#https://github.com/rpetit3/fastq-dl)

# fastq-dl

Download FASTQ files from the [European Nucleotide Archive](https://www.ebi.ac.uk/ena) or the
[Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) repositories.

## Introduction

`fastq-dl` takes an ENA/SRA accession (Study, Sample, Experiment, or Run) and queries ENA (via
[Data Warehouse API](https://www.ebi.ac.uk/ena/browse/search-rest)) to determine the associated
metadata. It then downloads FASTQ files for each Run. For Samples or Experiments with multiple
Runs, users can optionally merge the runs.

## Installation

### Bioconda

`fastq-dl` is available from [Bioconda](https://bioconda.github.io/) and I highly recommend you
go this route to for installation.

```{bash}
conda create -n fastq-dl -c conda-forge -c bioconda fastq-dl
conda activate fastq-dl 
```

## Usage

```{bash}
fastq-dl --help
                                                                                          
 Usage: fastq-dl [OPTIONS]                                                                
                                                                                          
 Download FASTQ files from ENA or SRA.                                                    
                                                                                          
╭─ Required Options ─────────────────────────────────────────────────────────────────────╮
│ *  --accession  -a  TEXT  ENA/SRA accession to query. (Study, Sample, Experiment, Run  │
│                           accession) [required]                                        │
╰────────────────────────────────────────────────────────────────────────────────────────╯
╭─ Additional Options ───────────────────────────────────────────────────────────────────╮
│ --provider                 [ena|sra]  Specify which provider (ENA or SRA) to use.      │
│                                       [default: ena]                                   │
│ --group-by-experiment                 Group Runs by experiment accession.              │
│ --group-by-sample                     Group Runs by sample accession.                  │
│ --outdir               -o  TEXT       Directory to output downloads to. [default: ./]  │
│ --prefix                   TEXT       Prefix to use for naming log files.              │
│                                       [default: fastq]                                 │
│ --cpus                     INTEGER    Total cpus used for downloading from SRA.        │
│                                       [default: 1]                                     │
│ --max-attempts         -m  INTEGER    Maximum number of download attempts.             │
│                                       [default: 10]                                    │
│ --only-provider        -F             Only attempt download from specified provider.   │
│ --silent                              Only critical errors will be printed.            │
│ --version              -V             Show the version and exit.                       │
│ --verbose              -v             Print debug related text.                        │
│ --help                 -h             Show this message and exit.                      │
╰────────────────────────────────────────────────────────────────────────────────────────╯
```

*fastq-dl* requires a single ENA/SRA Study, Sample, Experiment, or Run accession and FASTQs
for all Runs that fall under the given accession will be downloaded. For example, if a Study
accession is given all Runs under that studies umbrella will be downloaded. By default, 
`fastq-dl` will try to download from ENA first, then SRA.

### --accession

The accession you would like to download associated FASTQS for. Currently the following types
of accessions are accepted.

| Accession Type | Prefixes            | Example                                  |
|----------------|---------------------|------------------------------------------|
| BioProject     | PRJEB, PRJNA, PRJDB | PRJEB42779, PRJNA480016, PRJDB14838      |
| Study          | ERP, DRP, SRP       | ERP126685, DRP009283, SRP158268          |
| BioSample      | SAMD, SAME, SAMN    | SAMD00258402, SAMEA7997453, SAMN06479985 |
| Sample         | ERS, DRS, SRS       | ERS5684710, DRS259711, SRS2024210        |
| Experiment     | ERX, DRX, SRX       | ERX5050800, DRX406443, SRX4563689        |
| Run            | ERR, DRR, SRR       | ERR5260405, DRR421224, SRR7706354        |

The accessions are using regular expressions from the [ENA Training Modules - Accession Numbers](https://ena-docs.readthedocs.io/en/latest/submit/general-guide/accessions.html#accession-numbers) section.

### --provider

`fastq-dl` gives you the option to download from ENA or SRA. the `--provider` option will
specify which provider you would like to attempt downloads from first. If a download fails
from the first provider, additional attempts will be made using the other provider.

ENA was selected as the default provider because the FASTQs are available directly without
the need for conversion.

### --only-provider

By default, `fastq-dl` will fallback on a secondary provider to attempt downloads. There
may be cases where you would prefer to disable this feature, and that is exactly the
purpose of `--only-provider`. When provided, if a FASTQ cannot be downloaded from the
original provider, no additional attempts will be made.

### --group-by-experiment & --group-by-sample

There maybe times you might want to group Run accessions based on a Experiment or Sample
accessions. This will merge FASTQs associated with a Run accession based its associated
Experiment accession (`--group-by-experiment`) or Sample accession (`--group-by-sample`).

## Output Files

| Extension          | Description                                                                              |
|--------------------|------------------------------------------------------------------------------------------|
| `-run-info.tsv`    | Tab-delimited file containing metadata for each Run downloaded                           |
| `-run-mergers.tsv` | Tab-delimited file merge information from `--group-by-experiment` or `--group-by-sample` |
| `.fastq.gz`        | FASTQ files downloaded from ENA or SRA                                                   |

## Example Usage

#### Download FASTQs associated with a Study

Sometimes you might be reading a paper and they very kindly provided a Bioproject of all
the samples they sequenced. So, you decide you want to download FASTQs for all the samples
asscociated with the Bioproject. `fastq-dl` can help you with that! 

```{bash}
fastq-dl --accession PRJNA248678 --provider SRA
fastq-dl --accession PRJNA248678
```

The above commands will download the 3 Runs that fall under Study accession [PRJNA248678](https://www.ebi.ac.uk/ena/browser/view/PRJNA248678)
from either SRA (`--provider SRA`) or ENA (without `--provider`).

#### Download FASTQs associated with an Experiment

Let's say instead of the whole Bioproject you just want a single Experiment. You can do
that as well.

```{bash}
fastq-dl --accession SRX477044
```

The above command would download the Run accessions from ENA that fall under Experiment SRX477044.

The relationship of Experiment to Run is a 1-to-many relationship, or there can be many Run accessions
associated with a single Experiment Accession (e.g. re-sequencing the same sample). Although in most
cases, it is a 1-to-1 relationship, you can use `--group-by-experiment` to merge multiple runs
associated with an Experiment accession into a single FASTQ file.

#### Download FASTQs associated with an Sample

Ok, this time you just want a single Sample, or Biosample.

```{bash}
fastq-dl --accession SRS1904245 --provider SRA
```

The above command would download the Run accessions from SRA that fall under Sample SRS1904245.

Similar to Experiment accessions, the relationship of Sample to Run is a 1-to-many relationship,
or there can be many Run accessions associated with a single Sample Accession. Although in most
cases, it is a 1-to-1 relationship, you can use `--group-by-sample` to merge multiple runs
associated with an Sample accession into a single FASTQ file.

_Warning! For some type strains (e.g. S. aureus USA300) a Biosample accession might be associated with
100s or 1000s of Run accessions. These Runs are likely associated with many different conditions and
really should not fall under a single BioSample accession. Please consider this when using
`--group-by-sample`.

#### Download FASTQs associated with a Run

Let's keep it super simple and just download a Run.

```
fastq-dl --accession SRR1178105 --provider SRA
```

The above command would download the Run SRR1178105 from SRA. Run accessions are the end of the
line (1-to-1 relationship), so you will always get the expected Run.

## Alternatives
`fastq-dl`, is a spin-off of [ena-dl](https://github.com/rpetit3/ena-dl), that has been developed for
usage with [Bactopia](https://github.com/bactopia/bactopia). With this in mind, EBI/NCBI and provide
their own tools ([enaBrowserTools](https://github.com/enasequence/enaBrowserTools) and
[SRA Toolkit](https://github.com/ncbi/sra-tools)) that offer more extensive access to their databases.
