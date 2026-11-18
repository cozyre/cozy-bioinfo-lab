"""Network I/O layer: UniProt, NCBI BLAST/Entrez and EBI Clustal Omega.

Sessions are injected by the caller (see §6.2); no module here imports
Streamlit or writes to disk.
"""

from __future__ import annotations

from tp53.clients.base import (
    RemoteHTTPError,
    RemoteServiceError,
    RemoteTimeout,
    request_with_retry,
)

__all__ = [
    "RemoteHTTPError",
    "RemoteServiceError",
    "RemoteTimeout",
    "request_with_retry",
]
