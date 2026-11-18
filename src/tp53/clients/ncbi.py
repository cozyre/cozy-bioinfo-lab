"""NCBI Entrez efetch and remote BLAST URL API client (FR-06, FR-07, FR-12).

NCBI requires a contact email and a tool name on every request; both are sent
here. Polling obeys the current NCBI usage guidelines: at most one poll per
RID per minute (verified 2026-11-18), with an overall wall-clock deadline.
"""

from __future__ import annotations

import re
import time
from collections.abc import Callable
from typing import Any

from tp53 import config
from tp53.clients.base import (
    RemoteServiceError,
    RemoteTimeout,
    request_with_retry,
)
from tp53.config import BlastParams
from tp53.fasta_io import parse_fasta
from tp53.models import SequenceRecord, ValidationError

_RID_RE = re.compile(r"\bRID\s*=\s*([A-Za-z0-9_-]+)")
_STATUS_RE = re.compile(r"\bStatus\s*=\s*([A-Za-z]+)")
_ERROR_DIV_RE = re.compile(r'<div class="error msInf">(.*?)</div>', re.DOTALL)
_ERROR_PARA_RE = re.compile(r'<p class="error">(.*?)</p>', re.DOTALL)
_ERROR_ID_RE = re.compile(r"(Message ID#.*)", re.DOTALL)


def _ncbi_error_message(text: str) -> str:
    """Extract an error message from an NCBI HTML response, if present."""
    for pattern in (_ERROR_DIV_RE, _ERROR_PARA_RE, _ERROR_ID_RE):
        match = pattern.search(text)
        if match:
            return re.sub(r"\s+", " ", match.group(1)).strip()
    return ""


def _identification(email: str, api_key: str | None) -> dict[str, str]:
    """Parameters NCBI uses to identify and rate-limit the caller."""
    params = {"email": email, "tool": config.NCBI_TOOL_NAME}
    if api_key:
        params["api_key"] = api_key
    return params


def fetch_protein(
    accession: str,
    session: Any,
    email: str,
    *,
    api_key: str | None = None,
    timeout: int = config.REQUEST_TIMEOUT_S,
) -> SequenceRecord:
    """Fetch a full-length protein sequence from Entrez efetch (FR-12).

    Args:
        accession: NCBI protein accession, e.g. ``"NP_035770.2"``.
        session: injected requests-like session.
        email: user-supplied contact email sent to NCBI.
        api_key: optional NCBI API key (raises rate limits).
        timeout: per-request timeout in seconds.

    Returns:
        The single parsed FASTA record from efetch.

    Raises:
        RemoteServiceError: efetch returned no FASTA record.
    """
    params = {
        "db": "protein",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text",
        **_identification(email, api_key),
    }
    response = request_with_retry(
        session, "GET", config.NCBI_EFETCH_URL, timeout=timeout, params=params
    )
    try:
        records = parse_fasta(response.text)
    except ValidationError as exc:
        raise RemoteServiceError(
            f"NCBI efetch returned no FASTA record for {accession}: {exc}"
        ) from exc
    if not records:
        raise RemoteServiceError(f"NCBI efetch returned no sequence for {accession}.")
    record = records[0]
    organism_match = re.search(r"\[([^\[\]]+)\]\s*$", record.description)
    if organism_match:
        record.organism = organism_match.group(1)
        record.description = record.description[: organism_match.start()].strip()
    return record


def submit_blastp(
    seq: str,
    organism: str,
    params: BlastParams,
    email: str,
    session: Any,
    *,
    api_key: str | None = None,
    timeout: int = config.REQUEST_TIMEOUT_S,
) -> str:
    """Submit one organism-restricted blastp job and return its RID (FR-07).

    Args:
        seq: query protein sequence.
        organism: target species used as an ``[Organism]`` Entrez filter.
        params: BLAST program/database/E-value/hit-list settings.
        email: user-supplied contact email sent to NCBI.
        session: injected requests-like session.
        api_key: optional NCBI API key.
        timeout: per-request timeout in seconds.

    Returns:
        NCBI request identifier (RID) used by :func:`poll_blast`.

    Raises:
        RemoteServiceError: NCBI did not return an RID.
    """
    data = {
        "CMD": "Put",
        "PROGRAM": params.program,
        "DATABASE": params.database,
        "QUERY": seq,
        "ENTREZ_QUERY": f'"{organism}"[Organism]',
        "EXPECT": str(params.expect),
        "HITLIST_SIZE": str(params.hitlist_size),
        "FORMAT_TYPE": config.BLAST_RESULT_FORMAT,
        **_identification(email, api_key),
    }
    response = request_with_retry(
        session, "POST", config.NCBI_BLAST_URL, timeout=timeout, data=data
    )
    text = response.text
    match = _RID_RE.search(text)
    if not match:
        detail = _ncbi_error_message(text) or "no request identifier returned"
        raise RemoteServiceError(f"NCBI BLAST submission failed: {detail}.")
    return match.group(1)


def submit_blastp_staggered(
    jobs: list[tuple[str, str]],
    params: BlastParams,
    email: str,
    session: Any,
    *,
    stagger_delay_s: int = config.BLAST_STAGGER_DELAY_S,
    sleep: Callable[[float], None] = time.sleep,
    timeout: int = config.REQUEST_TIMEOUT_S,
    api_key: str | None = None,
) -> dict[str, str]:
    """Submit one blastp job per species, spacing submissions out (§9.1).

    Args:
        jobs: ``(species, query_sequence)`` pairs, one per target species.
        params: BLAST parameters shared by all jobs.
        email: user-supplied contact email sent to NCBI.
        session: injected requests-like session.
        stagger_delay_s: pause between submissions (NCBI etiquette).
        sleep: sleep function, injectable for offline tests.
        timeout: per-request timeout in seconds.
        api_key: optional NCBI API key.

    Returns:
        Mapping of species to its RID, in submission order.
    """
    rids: dict[str, str] = {}
    for index, (species, sequence) in enumerate(jobs):
        if index:
            sleep(stagger_delay_s)
        rids[species] = submit_blastp(
            sequence,
            species,
            params,
            email,
            session,
            api_key=api_key,
            timeout=timeout,
        )
    return rids


def _status(text: str) -> str:
    """Status of a BLAST job response: WAITING/RUNNING/READY/FAILED/UNKNOWN."""
    match = _STATUS_RE.search(text)
    return match.group(1).upper() if match else "UNKNOWN"


def _is_blast_xml(text: str) -> bool:
    """Whether the response is the final XML result rather than a status page."""
    return "<BlastOutput" in text


def poll_blast(
    rid: str,
    session: Any,
    timeout_s: int,
    *,
    email: str = "",
    api_key: str | None = None,
    poll_interval: int = config.BLAST_POLL_INTERVAL_S,
    sleep: Callable[[float], None] = time.sleep,
    clock: Callable[[], float] = time.monotonic,
    timeout: int = config.REQUEST_TIMEOUT_S,
) -> str:
    """Poll a BLAST job until it is ready and return the XML text (FR-07).

    One RID is polled at most once per ``poll_interval`` seconds (60 s by
    default, per NCBI's guidelines). Tests inject ``sleep`` and ``clock`` so
    no real waiting is needed.

    Args:
        rid: RID returned by :func:`submit_blastp`.
        session: injected requests-like session.
        timeout_s: overall wall-clock deadline for this job in seconds.
        email: contact email sent on every poll.
        api_key: optional NCBI API key.
        poll_interval: minimum seconds between two polls of this RID.
        sleep: sleep function, injectable for offline tests.
        clock: monotonic clock, injectable for offline tests.
        timeout: per-request timeout in seconds.

    Returns:
        Raw BLAST XML text.

    Raises:
        RemoteTimeout: the job did not finish before ``timeout_s``.
        RemoteServiceError: the job failed or the response was unrecognised.
    """
    params = {
        "CMD": "Get",
        "RID": rid,
        "FORMAT_TYPE": config.BLAST_RESULT_FORMAT,
        **_identification(email, api_key),
    }
    start = clock()
    while clock() - start < timeout_s:
        response = request_with_retry(
            session, "GET", config.NCBI_BLAST_URL, timeout=timeout, params=params
        )
        text = response.text
        status = _status(text)
        if status == "FAILED":
            detail = _ncbi_error_message(text) or "BLAST job failed"
            raise RemoteServiceError(f"NCBI BLAST job {rid} failed: {detail}.")
        if _is_blast_xml(text):
            return text
        if status not in {"WAITING", "RUNNING", "READY", "UNKNOWN"}:
            raise RemoteServiceError(
                f"NCBI BLAST job {rid} returned an unexpected status: {status}."
            )
        sleep(poll_interval)
    raise RemoteTimeout(f"NCBI BLAST job {rid} did not finish within {timeout_s} s.")
