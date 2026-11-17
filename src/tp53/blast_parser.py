"""Parse NCBI BLAST XML into project dataclasses (§7.2).

Wraps ``Bio.Blast.NCBIXML`` so the rest of the code only sees
:class:`tp53.models.BlastResult`.
"""

from __future__ import annotations

import io
import re

from Bio.Blast import NCBIXML

from tp53.models import BlastHit, BlastHsp, BlastResult

#: Organism tag inside a RefSeq hit title, e.g. ``... [Mus musculus]``.
_ORGANISM_RE = re.compile(r"\[([^\[\]]+)\]\s*$")


def _counts(value) -> tuple[int, int]:
    """Normalise an NCBIXML counter (identity/positives/gaps).

    Current Biopython returns a ``(count, total)`` tuple; older versions
    returned the string ``"count/total (pct%)"``. Both are accepted.
    """
    if value is None:
        return 0, 0
    if isinstance(value, (int, float)):
        return int(value), int(value)
    if isinstance(value, tuple):
        count, total = value
        return int(count or 0), int(total or 0)
    match = re.match(r"\s*(\d+)\s*/\s*(\d+)", str(value))
    if match:
        return int(match.group(1)), int(match.group(2))
    return 0, 0


def _parse_accession(hit_id: str) -> str:
    """Extract the accession from a BLAST ``hit_id`` like ``ref|NP_123.4|``."""
    parts = hit_id.split("|")
    if len(parts) >= 2:
        return parts[1]
    return hit_id


def _parse_organism(title: str) -> str:
    """Extract the organism from a hit title ending in ``[Organism]``."""
    match = _ORGANISM_RE.search(title)
    return match.group(1) if match else ""


def parse_blast_xml(xml_text: str, species: str = "") -> BlastResult:
    """Parse BLAST XML text into a :class:`BlastResult`.

    Args:
        xml_text: raw XML from an NCBI BLASTp job (blastp XML format).
        species: the organism filter that was used for this job, if any.

    Returns:
        Parsed result with all hits and their HSPs, in BLAST ranking order.

    Raises:
        ValueError: if the XML holds no BLAST result.
    """
    blast_records = list(NCBIXML.parse(io.StringIO(xml_text)))
    if not blast_records:
        raise ValueError("No BLAST result found in the XML.")
    record = blast_records[0]

    hits: list[BlastHit] = []
    for alignment in record.alignments:
        hsps = [
            BlastHsp(
                bit_score=float(hsp.bits),
                evalue=float(hsp.expect),
                identity=_counts(hsp.identities)[0],
                positives=_counts(hsp.positives)[0],
                gaps=_counts(hsp.gaps)[0],
                align_len=int(hsp.align_length),
                query_start=int(hsp.query_start),
                query_end=int(hsp.query_end),
                hit_start=int(hsp.sbjct_start),
                hit_end=int(hsp.sbjct_end),
                qseq=hsp.query,
                hseq=hsp.sbjct,
                midline=hsp.match,
            )
            for hsp in alignment.hsps
        ]
        hits.append(
            BlastHit(
                accession=_parse_accession(alignment.hit_id),
                title=alignment.title,
                organism=_parse_organism(alignment.title),
                length=int(alignment.length),
                hsps=hsps,
            )
        )

    return BlastResult(
        query_id=record.query,
        query_len=int(record.query_letters),
        database=record.database,
        program=record.application,
        hits=hits,
        species=species,
    )
