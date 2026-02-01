""" Constants used in the fastq_dl package """

# ENA Related
ENA = "ENA"
ENA_FAILED = "ENA_NOT_FOUND"
ENA_URL = "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&format=tsv"

# SRA Related
SRA = "SRA"
SRA_FAILED = "SRA_NOT_FOUND"

# File suffixes
RUN_INFO_SUFFIX = "-run-info.tsv"
RUN_MERGERS_SUFFIX = "-run-mergers.tsv"

# FASTQ suffixes (SRA/ENA raw downloads)
PE_R1_SUFFIX = "_1.fastq.gz"
PE_R2_SUFFIX = "_2.fastq.gz"
SE_SUFFIX = ".fastq.gz"

# Merged FASTQ suffixes
MERGED_R1_SUFFIX = "_R1.fastq.gz"
MERGED_R2_SUFFIX = "_R2.fastq.gz"
