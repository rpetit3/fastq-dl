[![Anaconda-Server Badge](https://anaconda.org/bioconda/fastq-dl/badges/downloads.svg)](https://anaconda.org/bioconda/fastq-dl)

# fastq-dl
fastq-dl is a tool for downloading FASTQ files from the [European Nucleotide Archive](https://www.ebi.ac.uk/ena) or the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra).

fastq-dl takes an ENA/SRA accession (Study, Experiment, or Run) and queries ENA (via [Data Warehouse API](https://www.ebi.ac.uk/ena/browse/search-rest)) to determine the associated metadata. It then downloads FASTQ files for each Run. For Samples or Experiments with multiple Runs, users can optionally merge the runs.

### Alternatives
fastq-dl, is a spin-off of [ena-dl](https://github.com/rpetit3/ena-dl), that has been developed for usage with [Bactopia](https://github.com/bactopia/bactopia). With this in mind, EBI/NCBI and provide their own tools ([enaBrowserTools](https://github.com/enasequence/enaBrowserTools) and [SRA Toolkit](https://github.com/ncbi/sra-tools)) that offer more extensive access to their databases.

# Installation
### Bioconda
fastq-dl is available from [Bioconda](https://bioconda.github.io/) and I highly recommend you go this route to for installation.
```
conda create -n fastq-dl -c conda-forge -c bioconda fastq-dl
conda activate fastq-dl 
```

# Usage
*fastq-dl* requires a single ENA/SRA Study, Experiment, or Run accession and FASTQs for all Runs that fall under the given accession will be downloaded. For example, if a Study accession is given all Runs under that studies umbrella will be downloaded. Which archive to download from, either *ENA* or *SRA*, is also required.

### Usage Output
```{bash}
fastq-dl --help

 Usage: fastq-dl [OPTIONS]

 Download FASTQ files from ENA or SRA.

╭─ Options ────────────────────────────────────────────────────────────────────────────────────────────╮
│    --version                             Show the version and exit.                                  │
│ *  --accession                TEXT       ENA/SRA accession to query. (Study, Sample, Experiment, Run │
│                                          accession)                                                  │
│                                          [required]                                                  │
│    --provider                 [ena|sra]  Specify which provider (ENA or SRA) to use. [default: ena]  │
│    --group_by_experiment                 Group Runs by experiment accession.                         │
│    --group_by_sample                     Group Runs by sample accession.                             │
│    --outdir               -o  TEXT       Directory to output downloads to. [default: ./]             │
│    --prefix                   TEXT       Prefix to use for naming log files. [default: fastq]        │
│    --max_attempts         -a  INTEGER    Maximum number of download attempts. [default: 10]          │
│    --only-provider        -F             Only attempt download from specified provider.              │
│    --cpus                     INTEGER    Total cpus used for downloading from SRA. [default: 1]      │
│    --silent                              Only critical errors will be printed.                       │
│    --verbose              -v             Print debug related text.                                   │
│    --debug                               Skip downloads, print what will be downloaded.              │
│    --help                                Show this message and exit.                                 │
╰──────────────────────────────────────────────────────────────────────────────────────────────────────╯
```

### Example Usage
#### Download a Study
```
fastq-dl --accession PRJNA248678 --provider SRA
fastq-dl --accession PRJNA248678
```

The above commands will download 3 runs that fall under Study PRJNA248678 from either SRA (`fastq-dl --accession PRJNA248678 --provider SRA`) or ENA (`fastq-dl --accession PRJNA248678`). The relationship of Study to Run is a 1-to-many relationship, or there can be many Run accessions associated with a single Study Accession. You could use `--group_by_experiment` to group these runs by Experiment accession (or `--group_by_sample` for Sample accession).

#### Download an Experiment
```
fastq-dl --accession SRX477044
```

The above command would download the single run from ENA that falls under Experiment SRX477044. The relationship of Experiment to Run is a 1-to-many relationship, or there can be many Run accessions associated with a single Experiment Accession (e.g. re-sequencing the same sample). Although in most cases, especially for bacterial samples, it is a 1-to-1 relationship. In any case, you can use `--group_by_experiment` to merge multiple runs associated with an Experiment accession into a single FASTQ file (or `--group_by_sample` for Sample accession).

#### Download a Run
```
fastq-dl --accession SRR1178105 --provider SRA
```

The above command would download the Run SRR1178105 from SRA. Run accessions are the end of the line (1-to-1 relationship)
