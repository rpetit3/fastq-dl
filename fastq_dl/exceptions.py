"""Custom exceptions for fastq-dl."""


class FastqDLError(Exception):
    """Base exception for fastq-dl."""

    pass


class ValidationError(FastqDLError):
    """Invalid accession or input."""

    pass


class ProviderError(FastqDLError):
    """Error from a data provider (ENA/SRA)."""

    def __init__(self, message: str, provider: str, status_code: int = None):
        self.provider = provider
        self.status_code = status_code
        super().__init__(message)


class DownloadError(FastqDLError):
    """Error during file download."""

    def __init__(self, message: str, accession: str = None, provider: str = None):
        self.accession = accession
        self.provider = provider
        super().__init__(message)
