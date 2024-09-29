from pathlib import Path
from typing import Union

PROGRAM = "fastq-dl"
VERSION = "2.0.4"
ENA_FAILED = "ENA_NOT_FOUND"
SRA_FAILED = "SRA_NOT_FOUND"
SRA = "SRA"
ENA = "ENA"
ENA_URL = "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&format=tsv"
PathLike = Union[str, Path]
