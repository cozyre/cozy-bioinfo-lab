"""EBI Clustal Omega REST client (FR-14).

Asynchronous run → status → result flow. A valid contact email is required
by EMBL-EBI's fair-use policy and is sent on submission. The job dispatcher
returns aligned FASTA from the ``out`` result type (see config.py).
"""

from __future__ import annotations

import time
from collections.abc import Callable
from typing import Any

from tp53 import config
from tp53.clients.base import (
    RemoteServiceError,
    RemoteTimeout,
    request_with_retry,
)

_FINISHED = "FINISHED"
_FAILURE_STATUSES = ("FAILED", "ERROR", "NOT_FOUND")


def run_clustalo(
    fasta: str,
    email: str,
    session: Any,
    timeout_s: int,
    *,
    poll_interval: int = config.EBI_POLL_INTERVAL_S,
    sleep: Callable[[float], None] = time.sleep,
    clock: Callable[[], float] = time.monotonic,
    timeout: int = config.REQUEST_TIMEOUT_S,
) -> str:
    """Submit sequences to Clustal Omega and return the aligned FASTA.

    Args:
        fasta: multi-FASTA of query plus selected orthologs.
        email: user-supplied contact email sent to EBI.
        session: injected requests-like session.
        timeout_s: overall wall-clock deadline for the job in seconds.
        poll_interval: seconds between status checks.
        sleep: sleep function, injectable for offline tests.
        clock: monotonic clock, injectable for offline tests.
        timeout: per-request timeout in seconds.

    Returns:
        Aligned FASTA text from the ``out`` result type.

    Raises:
        RemoteServiceError: submission failed, no job ID was returned, or the
            job failed.
        RemoteTimeout: the job did not finish before ``timeout_s``.
    """
    response = request_with_retry(
        session,
        "POST",
        f"{config.EBI_CLUSTALO_URL}/run",
        timeout=timeout,
        data={
            "email": email,
            "sequence": fasta,
            "stype": config.EBI_CLUSTALO_SEQUENCE_TYPE,
        },
    )
    job_id = response.text.strip().splitlines()[0].strip() if response.text.strip() else ""
    if not job_id:
        raise RemoteServiceError("EBI Clustal Omega returned no job ID.")

    status_url = f"{config.EBI_CLUSTALO_URL}/status/{job_id}"
    result_url = (
        f"{config.EBI_CLUSTALO_URL}/result/{job_id}/{config.EBI_CLUSTALO_RESULT_TYPE}"
    )
    start = clock()
    while clock() - start < timeout_s:
        status_response = request_with_retry(
            session, "GET", status_url, timeout=timeout
        )
        status = status_response.text.strip().upper()
        if status == _FINISHED:
            result = request_with_retry(
                session, "GET", result_url, timeout=timeout
            )
            return result.text
        if status in _FAILURE_STATUSES:
            raise RemoteServiceError(
                f"EBI Clustal Omega job {job_id} failed with status {status}."
            )
        sleep(poll_interval)
    raise RemoteTimeout(
        f"EBI Clustal Omega job {job_id} did not finish within {timeout_s} s."
    )
