"""Deterministic best-hit selection per species (§9.2, FR-11)."""

from __future__ import annotations

import warnings
from dataclasses import dataclass

from tp53 import config
from tp53.models import BlastHit, BlastHsp, BlastResult, SelectedHit

_PREFIX_RANK = {"NP_": 0, "XP_": 1}
_OTHER_ACCESSION_RANK = 2


@dataclass(frozen=True)
class _Candidate:
    """A hit that passed the filters, with the data needed to rank it (§9.2)."""

    hit: BlastHit
    hsp: BlastHsp
    coverage: float


def _accession_rank(accession: str) -> int:
    """Prefer RefSeq NP_ accessions, then XP_, then all others (§9.2)."""
    for prefix, rank in _PREFIX_RANK.items():
        if accession.upper().startswith(prefix):
            return rank
    return _OTHER_ACCESSION_RANK


def _rank_key(candidate: _Candidate) -> tuple:
    """Sort key implementing the four-step ranking of §9.2."""
    return (
        -candidate.hsp.bit_score,
        -candidate.hsp.percent_identity,
        _accession_rank(candidate.hit.accession),
        candidate.hit.accession.upper(),
    )


def _reason(ranked: list[_Candidate]) -> str:
    """Explain the choice, including the tie-break that decided it, if any."""
    winner = ranked[0]
    parts = [
        f"bit score {winner.hsp.bit_score:.1f}",
        f"identity {winner.hsp.percent_identity:.1f}%",
        f"query coverage {winner.coverage:.1%}",
    ]
    if len(ranked) > 1:
        runner_up = ranked[1]
        if winner.hsp.bit_score == runner_up.hsp.bit_score:
            if winner.hsp.percent_identity != runner_up.hsp.percent_identity:
                parts.append("bit-score tie broken by percent identity")
            elif _accession_rank(winner.hit.accession) != _accession_rank(
                runner_up.hit.accession
            ):
                parts.append("bit-score and identity tie broken by NP_ over XP_")
            else:
                parts.append("full tie broken by accession alphabetically")
    return "; ".join(parts)


def select_best_hit(
    result: BlastResult,
    params: config.SelectionParams | None = None,
) -> SelectedHit | None:
    """Pick the best ortholog hit for one species using §9.2.

    Filters hits by E-value and query coverage, then ranks by bit score,
    percent identity, NP_ over XP_, and accession alphabetically. The result
    is deterministic regardless of the order of hits in the BLAST XML.

    Args:
        result: parsed blastp result for one organism-restricted job.
        params: selection thresholds; defaults to :class:`SelectionParams`.

    Returns:
        Selected hit with its ranked candidates and a human-readable reason,
        or ``None`` (with a ``UserWarning``) when no hit qualifies.
    """
    params = params or config.SelectionParams()
    ranked: list[_Candidate] = []
    for hit in result.hits:
        if not hit.hsps:
            continue
        hsp = hit.best_hsp
        if hsp.evalue > params.evalue_cutoff:
            continue
        coverage = (
            (hsp.query_end - hsp.query_start + 1) / result.query_len
            if result.query_len
            else 0.0
        )
        if coverage < params.min_query_coverage:
            continue
        ranked.append(_Candidate(hit=hit, hsp=hsp, coverage=coverage))

    ranked.sort(key=_rank_key)
    if not ranked:
        species = result.species or result.query_id or "query"
        warnings.warn(
            f"No BLAST hit for {species} passed the E-value "
            f"({params.evalue_cutoff:g}) and query-coverage "
            f"({params.min_query_coverage:.0%}) filters.",
            UserWarning,
            stacklevel=2,
        )
        return None

    winner = ranked[0]
    return SelectedHit(
        species=result.species,
        hit=winner.hit,
        query_coverage=winner.coverage,
        percent_identity=winner.hsp.percent_identity,
        reason=_reason(ranked),
        candidates=[candidate.hit for candidate in ranked],
    )
