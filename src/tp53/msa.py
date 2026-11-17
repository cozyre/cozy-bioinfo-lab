"""MSA parsing, validation and coordinate mapping (§9.4, §9.7, FR-16, FR-17)."""

from __future__ import annotations

import io

import pandas as pd
from Bio import SeqIO

from tp53 import config
from tp53.constants import (
    ALIGNMENT_STAT_COLUMNS,
    ERR_BAD_CHARS,
    ERR_DEGAP_MISMATCH,
    ERR_EMPTY_ALIGNMENT,
    ERR_REF_ROW_MISSING,
    ERR_ROW_COUNT,
    ERR_ROW_LENGTH,
    WARN_ID_REMAPPED,
)
from tp53.models import Alignment, AlignmentValidationError


def _normalise_id(identifier: str) -> str:
    """Case- and whitespace-insensitive form of a sequence ID."""
    return identifier.strip().upper()


def _find_reference_row(ids: list[str], ref_id: str) -> int:
    """Index of the reference (query) row, tolerant of minor ID rewrites.

    Matches in this order: exact, case/whitespace-insensitive, one ID being
    a prefix of the other (Clustal Omega can truncate long headers).
    """
    ref_norm = _normalise_id(ref_id)
    for index, identifier in enumerate(ids):
        if identifier == ref_id:
            return index
    for index, identifier in enumerate(ids):
        id_norm = _normalise_id(identifier)
        if id_norm == ref_norm:
            return index
    for index, identifier in enumerate(ids):
        id_norm = _normalise_id(identifier)
        if id_norm.startswith(ref_norm) or ref_norm.startswith(id_norm):
            return index
    raise AlignmentValidationError(ERR_REF_ROW_MISSING.format(ref_id=ref_id))


def _build_col_to_ref(rows: list[str], ref_index: int) -> list[int | None]:
    """Map each alignment column to a 1-based reference position (§9.7).

    Every non-gap reference character increments the counter; columns where
    the reference row has a gap map to ``None``. UniProt feature coordinates
    are looked up through this map, never by alignment column.
    """
    col_to_ref: list[int | None] = []
    ref_position = 0
    for character in rows[ref_index]:
        if character != config.GAP_CHAR:
            ref_position += 1
            col_to_ref.append(ref_position)
        else:
            col_to_ref.append(None)
    return col_to_ref


def parse_alignment(text: str, ref_id: str, expected_ids: list[str] | None = None) -> Alignment:
    """Parse an aligned multi-FASTA into an :class:`Alignment`.

    Sequences are upper-cased and whitespace-stripped. All-gap columns are
    removed (§9.4) before the column-to-reference map is built.

    Args:
        text: aligned FASTA text (Clustal Omega output or an uploaded MSA).
        ref_id: ID of the reference (query) row for the coordinate map.
        expected_ids: IDs of the input records, in input order. When the
            service rewrote the IDs, rows are remapped by order (§9.4).

    Returns:
        Alignment with ``col_to_ref`` built from the reference row.
        All-gap columns are NOT dropped here; call
        :func:`trim_all_gap_columns` after validation (it reports the count).

    Raises:
        AlignmentValidationError: if the reference row cannot be located.
    """
    ids: list[str] = []
    rows: list[str] = []
    for parsed in SeqIO.parse(io.StringIO(text), "fasta"):
        ids.append(parsed.id)
        rows.append(str(parsed.seq).upper().replace(" ", ""))
    if not rows:
        raise AlignmentValidationError(ERR_EMPTY_ALIGNMENT)

    ids, ref_index = _remap_ids(ids, ref_id, expected_ids)
    length = len(rows[0]) if rows else 0
    alignment = Alignment(ids=ids, rows=rows, length=length, ref_index=ref_index)
    alignment.col_to_ref = _build_col_to_ref(alignment.rows, alignment.ref_index)
    return alignment


def _remap_ids(ids: list[str], ref_id: str, expected_ids: list[str] | None) -> tuple[list[str], int]:
    """Resolve row identity: as parsed, or remapped by order (§9.4)."""
    if expected_ids and not all(
        _normalise_id(a) == _normalise_id(b)
        for a, b in zip(ids, expected_ids, strict=False)
    ):
        ids = list(expected_ids)
        ref_index = expected_ids.index(ref_id) if ref_id in expected_ids else None
        if ref_index is not None:
            return ids, ref_index
    return ids, _find_reference_row(ids, ref_id)


def validate_alignment(
    alignment: Alignment,
    expected_ids: list[str],
    expected_sequences: list[str] | None = None,
) -> list[str]:
    """Validate an alignment against the input records (§9.4, FR-16).

    Hard failures (raise :class:`AlignmentValidationError`, naming the
    offending sequence): wrong row count, unequal row lengths, characters
    outside the allowed alphabet, and de-gapped rows not reproducing their
    input sequence.

    Args:
        alignment: parsed alignment.
        expected_ids: IDs of the input records, in input order.
        expected_sequences: input sequences, in the same order. Given, the
            §9.4 hard consistency check runs on every row.

    Returns:
        Warnings, e.g. that IDs were remapped by order.
    """
    warnings: list[str] = []

    if len(alignment.rows) != len(expected_ids):
        raise AlignmentValidationError(
            ERR_ROW_COUNT.format(rows=len(alignment.rows), expected=len(expected_ids))
        )

    for index, row in enumerate(alignment.rows):
        identifier = alignment.ids[index]
        if len(row) != alignment.length:
            raise AlignmentValidationError(
                ERR_ROW_LENGTH.format(
                    id=identifier, length=len(row), expected=alignment.length
                )
            )
        invalid = sorted(set(row) - set(config.PROTEIN_ALPHABET) - {config.GAP_CHAR})
        if invalid:
            raise AlignmentValidationError(
                ERR_BAD_CHARS.format(id=identifier, chars=" ".join(invalid))
            )
        if expected_sequences is not None:
            expected = expected_sequences[index].upper().replace("-", "")
            if alignment.ungapped_row(index) != expected:
                raise AlignmentValidationError(
                    ERR_DEGAP_MISMATCH.format(id=identifier)
                )

    if expected_ids and not all(
        _normalise_id(a) == _normalise_id(b)
        for a, b in zip(alignment.ids, expected_ids, strict=False)
    ):
        warnings.append(WARN_ID_REMAPPED)
    return warnings


def trim_all_gap_columns(alignment: Alignment) -> tuple[Alignment, int]:
    """Drop columns that are gaps in every row (§9.4, FR-16).

    Args:
        alignment: parsed alignment.

    Returns:
        Trimmed alignment and the number of removed columns.
    """
    keep = [
        j
        for j in range(alignment.length)
        if any(row[j] != config.GAP_CHAR for row in alignment.rows)
    ]
    removed = alignment.length - len(keep)
    trimmed = Alignment(
        ids=alignment.ids,
        rows=["".join(row[j] for j in keep) for row in alignment.rows],
        length=len(keep),
        ref_index=alignment.ref_index,
    )
    trimmed.col_to_ref = _build_col_to_ref(trimmed.rows, trimmed.ref_index)
    return trimmed, removed


def alignment_statistics(alignment: Alignment) -> pd.DataFrame:
    """Per-sequence alignment statistics (FR-16).

    Returns a table with columns ``id``, ``length`` (aligned), ``ungapped_length``
    and ``gap_pct``.
    """
    stats = []
    for identifier, row in zip(alignment.ids, alignment.rows, strict=True):
        ungapped = row.replace(config.GAP_CHAR, "")
        stats.append(
            {
                ALIGNMENT_STAT_COLUMNS[0]: identifier,
                ALIGNMENT_STAT_COLUMNS[1]: len(row),
                ALIGNMENT_STAT_COLUMNS[2]: len(ungapped),
                ALIGNMENT_STAT_COLUMNS[3]: round(
                    100.0 * (len(row) - len(ungapped)) / len(row), 1
                )
                if row
                else 0.0,
            }
        )
    return pd.DataFrame(stats, columns=ALIGNMENT_STAT_COLUMNS)
