[![Anaconda-Server Badge](https://anaconda.org/bioconda/fastq-dl/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda) [![Anaconda-Server Badge](https://anaconda.org/bioconda/fastq-dl/badges/downloads.svg)](https://anaconda.org/bioconda/fastq-dl)

# fastq-dl
fastq-dl is a tool for downloading FASTQ files from the [European Nucleotide Archive](https://www.ebi.ac.uk/ena) or the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra).

fastq-dl takes an ENA/SRA accession (Study, Experiment, or Run) and queries ENA (via [Data Warehouse API](https://www.ebi.ac.uk/ena/browse/search-rest)) to determine the associated metadata. It then downloads FASTQ files for each Run. For Samples or Experiments with multiple Runs, users can optionally merge the runs.

### Alternatives
fastq-dl, is a spin-off of [ena-dl](https://github.com/rpetit3/ena-dl), that has been developed for usage with [Bactopia](https://github.com/bactopia/bactopia). With this in mind, EBI/NCBI and provide their own tools ([enaBrowserTools](https://github.com/enasequence/enaBrowserTools) and [SRA Toolkit](https://github.com/ncbi/sra-tools)) that offer more extensive access to their databases.

# Installation
### Bioconda
fastq-dl is available from [Bioconda](https://bioconda.github.io/) and I highly recommend you go this route to for installation.
```
conda install -c conda-forge -c bioconda fastq-dl
```

# Usage
*fastq-dl* requires a single ENA/SRA Study, Experiment, or Run accession and FASTQs for all Runs that fall under the given accession will be downloaded. For example, if a Study accession is given all Runs under that studies umbrella will be downloaded. Which archive to download from, either *ENA* or *SRA*, is also required.

### Usage Output
```
usage: fastq-dl [-h] [--aspera STRING] [--aspera_key STRING]
                [--aspera_speed STRING] [--is_study] [--is_experiment]
                [--is_run] [--group_by_experiment] [--group_by_sample]
                [--outdir OUTPUT_DIR] [--prefix PREFIX] [--max_attempts INT]
                [--cpus INT] [--ftp_only] [--silent] [--verbose] [--debug]
                [--version]
                ACCESSION {sra,SRA,ena,ENA}

fastq-dl - Download FASTQs from ENA or SRA

optional arguments:
  -h, --help            show this help message and exit

Required Options:

  ACCESSION             ENA/SRA accession to query. (Study, Experiment, or Run
                        accession)
  {sra,SRA,ena,ENA}     Specify which provider (ENA or SRA) to use. Accepted
                        Values: ENA SRA

Aspera Connect Options:
  --aspera STRING       Path to the Aspera Connect tool "ascp" (Default:
                        "which ascp")
  --aspera_key STRING   Path to Aspera Connect private key, if not given,
                        guess based on ascp path
  --aspera_speed STRING
                        Speed at which Aspera Connect will download. (Default:
                        100M)

Query Related Options:
  --is_study            Query is a Study.
  --is_experiment       Query is an Experiment.
  --is_run              Query is a Run.
  --group_by_experiment
                        Group Runs by experiment accession.
  --group_by_sample     Group Runs by sample accession.

Helpful Options:
  --outdir OUTPUT_DIR   Directory to output downloads to. (Default: ./)
  --prefix PREFIX       Prefix to use for naming log files (Default: fastq)
  --max_attempts INT    Maximum number of download attempts (Default: 10)
  --cpus INT            Total cpus used for downloading from SRA (Default: 1)
  --ftp_only            FTP only downloads.
  --silent              Only critical errors will be printed.
  --verbose             Print debug related text.
  --debug               Skip downloads, print what will be downloaded.
  --version             show program's version number and exit
```

### Example Usage
#### Download a Study
```
fastq-dl PRJNA248678 SRA
fastq-dl PRJNA248678 ENA
```

The above commands will download 3 runs that fall under Study PRJNA248678 from either SRA (`fastq-dl PRJNA248678 SRA`) or ENA (`fastq-dl PRJNA248678 ENA`). The relationship of Study to Run is a 1-to-many relationship, or there can be many Run accessions associated with a single Study Accession. You could use `--group_by_experiment` to group these runs by Experiment accession (or `--group_by_sample` for Sample accession).

#### Download an Experiment
```
fastq-dl SRX477044 ENA
```

The above command would download the single run from ENA that falls under Experiment SRX477044. The relationship of Experiment to Run is a 1-to-many relationship, or there can be many Run accessions associated with a single Experiment Accession (e.g. resequencing the same sample). Although in most cases, especially for bacterial samples, it is a 1-to-1 relationship. In any case, you can use `--group_by_experiment` to merge multiple runs associated with an Experiment accession into a single FASTQ file (or `--group_by_sample` for Sample accession).

#### Download a Run
```
fastq-dl SRR1178105 SRA
```

The above command would download the Run SRR1178105 from SRA. Run accessions are the end of the line (1-to-1 relationship)

#### Using Aspera Connect
Users can use [Aspera Connect](https://downloads.asperasoft.com/connect2/) to speed up the download of FASTQs from ENA. Installation and setup of Aspera Connect is out of the scope of this documentation, but I can assure you its a simple installation.

```
fastq-dl SRR1178105 ENA --aspera /path/to/ascp
```

The above command will attempt to download SRR1178105 using the `ascp` tool. By default it will try to use the private key file that is included during the Aspera Connect installation. If it is not found you will need to use the `--aspera_key` parameter to specify its path.
