"""UniProt REST client: query sequence, organism, gene and features (FR-06).

The contact email is not needed by UniProt; NCBI/EBI clients require it.
All URLs and the list of interesting feature types live in :mod:`tp53.config`.
"""

from __future__ import annotations

from typing import Any
from urllib.parse import quote

from tp53 import config
from tp53.clients.base import RemoteServiceError, request_with_retry
from tp53.models import DomainFeature, QueryInfo, SequenceRecord


def _feature_name(feature: dict[str, Any]) -> str:
    """Human-readable name of one UniProt feature (type as fallback)."""
    description = str(feature.get("description") or "").strip()
    return description or str(feature.get("type") or "")


def _parse_features(payload: dict[str, Any]) -> list[DomainFeature]:
    """Select the annotated features relevant to conservation (§9.9).

    Only exact end points are supported; features with fuzzy or unknown
    boundaries are skipped because they cannot be mapped to reference
    positions reliably.
    """
    features: list[DomainFeature] = []
    for raw in payload.get("features", []):
        feature_type = raw.get("type")
        if feature_type not in config.UNIPROT_FEATURE_TYPES:
            continue
        location = raw.get("location") or {}
        start_point = location.get("start") or {}
        end_point = location.get("end") or {}
        if any(
            point.get("modifier", config.UNIPROT_EXACT_MODIFIER)
            != config.UNIPROT_EXACT_MODIFIER
            for point in (start_point, end_point)
        ):
            continue
        start = start_point.get("value")
        end = end_point.get("value")
        if not isinstance(start, int) or not isinstance(end, int):
            continue
        if start < 1 or end < start:
            continue
        features.append(
            DomainFeature(
                type=str(feature_type),
                name=_feature_name(raw),
                start=start,
                end=end,
            )
        )
    return features


def _description(payload: dict[str, Any]) -> str:
    """Protein name plus gene name, since models have no separate gene field."""
    protein_description = payload.get("proteinDescription") or {}
    recommended = protein_description.get("recommendedName") or {}
    submitted = protein_description.get("submissionNames") or []
    name = recommended.get("fullName", {}).get("value") or ""
    if not name and submitted:
        name = (submitted[0].get("fullName") or {}).get("value") or ""
    gene = ""
    for gene_payload in payload.get("genes", []):
        gene = (gene_payload.get("geneName") or {}).get("value") or ""
        if gene:
            break
    if gene and gene not in name:
        return f"{name} (gene: {gene})".strip() if name else f"Gene: {gene}"
    return name


def fetch_query_info(
    accession: str,
    session: Any,
    *,
    timeout: int = config.REQUEST_TIMEOUT_S,
) -> QueryInfo:
    """Fetch a UniProt entry as a :class:`QueryInfo` (FR-06).

    Args:
        accession: UniProt accession, e.g. ``"P04637"``.
        session: injected requests-like session.
        timeout: per-request timeout in seconds.

    Returns:
        Query sequence with organism, gene and selected feature annotations.

    Raises:
        RemoteServiceError: the response is not a usable UniProt JSON entry.
    """
    url = f"{config.UNIPROT_BASE_URL}/{quote(accession.strip(), safe='')}.json"
    response = request_with_retry(session, "GET", url, timeout=timeout)
    try:
        payload = response.json()
    except ValueError as exc:
        raise RemoteServiceError(f"UniProt returned invalid JSON for {accession}.") from exc
    if not isinstance(payload, dict):
        raise RemoteServiceError(f"UniProt returned no entry for {accession}.")

    sequence = str((payload.get("sequence") or {}).get("value") or "").strip()
    if not sequence:
        raise RemoteServiceError(f"UniProt returned no sequence for {accession}.")

    record = SequenceRecord(
        id=str(payload.get("primaryAccession") or accession),
        accession=str(payload.get("primaryAccession") or accession),
        organism=str((payload.get("organism") or {}).get("scientificName") or ""),
        description=_description(payload),
        sequence=sequence,
    )
    return QueryInfo(record=record, features=_parse_features(payload))
