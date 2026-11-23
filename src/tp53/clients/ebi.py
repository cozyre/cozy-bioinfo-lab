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
    PollOutcome,
    RemoteServiceError,
    poll_until,
    request_with_retry,
)

_FINISHED = "FINISHED"
_FAILURE_STATUSES = ("FAILED", "ERROR", "NOT_FOUND")


def _job_urls(job_id: str) -> tuple[str, str]:
    """Status and result URLs of one Clustal Omega job (§8)."""
    status = config.EBI_CLUSTALO_URL + config.EBI_CLUSTALO_STATUS_PATH.format(
        job_id=job_id
    )
    result = config.EBI_CLUSTALO_URL + config.EBI_CLUSTALO_RESULT_PATH.format(
        job_id=job_id, result_type=config.EBI_CLUSTALO_RESULT_TYPE
    )
    return status, result


def submit_clustalo(
    fasta: str,
    email: str,
    session: Any,
    *,
    timeout: int = config.REQUEST_TIMEOUT_S,
) -> str:
    """Submit one Clustal Omega job and return its job ID (FR-14, §12).

    Separating submission from polling lets the app persist the job ID as soon
    as it is issued, so a page refresh can resume the same alignment.

    Args:
        fasta: multi-FASTA of query plus selected orthologs.
        email: contact email sent to EBI as required by its fair-use policy.
        session: injected requests-like session.
        timeout: per-request timeout in seconds.

    Returns:
        EBI job identifier.

    Raises:
        RemoteServiceError: submission failed or no job ID was returned.
    """
    response = request_with_retry(
        session,
        "POST",
        f"{config.EBI_CLUSTALO_URL}{config.EBI_CLUSTALO_RUN_PATH}",
        timeout=timeout,
        data={
            "email": email,
            "sequence": fasta,
            "stype": config.EBI_CLUSTALO_SEQUENCE_TYPE,
        },
        retry=False,
    )
    job_id = response.text.strip().splitlines()[0].strip() if response.text.strip() else ""
    if not job_id:
        raise RemoteServiceError("EBI Clustal Omega returned no job ID.")
    return job_id


def poll_clustalo_once(
    job_id: str,
    session: Any,
    *,
    timeout: int = config.REQUEST_TIMEOUT_S,
) -> PollOutcome:
    """Check one Clustal Omega job and fetch its result when finished.

    Args:
        job_id: identifier returned by :func:`submit_clustalo`.
        session: injected requests-like session.
        timeout: per-request timeout in seconds.

    Returns:
        Pending, complete (with aligned FASTA), or failed outcome.
    """
    status_url, result_url = _job_urls(job_id)
    status_response = request_with_retry(session, "GET", status_url, timeout=timeout)
    status = status_response.text.strip().upper()
    if status == _FINISHED:
        result = request_with_retry(session, "GET", result_url, timeout=timeout)
        return PollOutcome.complete(result.text)
    if status in _FAILURE_STATUSES:
        return PollOutcome.failed(
            f"EBI Clustal Omega job {job_id} failed with status {status}."
        )
    return PollOutcome.pending()


def resume_clustalo(
    job_id: str,
    session: Any,
    timeout_s: int,
    *,
    poll_interval: int = config.EBI_POLL_INTERVAL_S,
    sleep: Callable[[float], None] = time.sleep,
    clock: Callable[[], float] = time.monotonic,
    timeout: int = config.REQUEST_TIMEOUT_S,
) -> str:
    """Poll a previously submitted Clustal Omega job until it finishes (§12)."""
    return poll_until(
        lambda: poll_clustalo_once(job_id, session, timeout=timeout),
        label=f"EBI Clustal Omega job {job_id}",
        timeout_s=timeout_s,
        interval_s=poll_interval,
        sleep=sleep,
        clock=clock,
    )


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
    job_id = submit_clustalo(fasta, email, session, timeout=timeout)
    return resume_clustalo(
        job_id,
        session,
        timeout_s,
        poll_interval=poll_interval,
        sleep=sleep,
        clock=clock,
        timeout=timeout,
    )
