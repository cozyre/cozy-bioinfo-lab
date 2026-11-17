"""FASTA parsing, validation and writing (§7.2)."""

from __future__ import annotations

import io
import textwrap

from Bio import SeqIO

from tp53 import config
from tp53.constants import (
    ERR_EMPTY_FASTA,
    ERR_EMPTY_SEQUENCE,
    ERR_INVALID_CHARS,
    ERR_LENGTH,
    ERR_NOT_FASTA,
    ERR_ONE_RECORD,
)
from tp53.models import SequenceRecord, ValidationError


def parse_fasta(text: str) -> list[SequenceRecord]:
    """Parse FASTA text into records.

    Sequences are upper-cased and stripped of whitespace. Gaps are kept so
    aligned FASTA can pass through the same parser.

    Args:
        text: FASTA text, possibly multi-line and multi-record.

    Returns:
        One :class:`SequenceRecord` per FASTA record, in file order.

    Raises:
        ValidationError: if the text does not look like FASTA at all.
    """
    text = text.strip()
    if not text:
        raise ValidationError(ERR_EMPTY_FASTA)
    if not text.startswith(">"):
        raise ValidationError(ERR_NOT_FASTA)

    records: list[SequenceRecord] = []
    for parsed in SeqIO.parse(io.StringIO(text), "fasta"):
        header = parsed.description.split(maxsplit=1)
        records.append(
            SequenceRecord(
                id=parsed.id,
                accession=parsed.id,
                organism="",
                description=header[1] if len(header) > 1 else "",
                sequence=str(parsed.seq).upper(),
            )
        )
    return records


def parse_single_fasta(text: str) -> SequenceRecord:
    """Parse FASTA text that must hold exactly one record (FR-03).

    Raises:
        ValidationError: if the text is not FASTA or holds more than one
            record.
    """
    records = parse_fasta(text)
    if len(records) != 1:
        raise ValidationError(ERR_ONE_RECORD.format(count=len(records)))
    return records[0]


def validate_protein(
    record: SequenceRecord,
    min_len: int = config.FASTA_MIN_LEN,
    max_len: int = config.FASTA_MAX_LEN,
) -> None:
    """Validate one protein record (FR-03).

    Checks amino-acid alphabet and length bounds (default 20-5000 aa).

    Raises:
        ValidationError: with a message naming the offending record.
    """
    if not record.sequence:
        raise ValidationError(ERR_EMPTY_SEQUENCE.format(id=record.id))
    bad = sorted(set(record.sequence) - set(config.PROTEIN_ALPHABET))
    if bad:
        raise ValidationError(
            ERR_INVALID_CHARS.format(id=record.id, chars=" ".join(bad))
        )
    if not min_len <= len(record.sequence) <= max_len:
        raise ValidationError(
            ERR_LENGTH.format(
                id=record.id,
                length=len(record.sequence),
                min_len=min_len,
                max_len=max_len,
            )
        )


def to_fasta(records: list[SequenceRecord], header_separator: str = config.HEADER_SEPARATOR) -> str:
    """Write records as FASTA text with ``Species|Accession`` headers (FR-13).

    Args:
        records: records in the order they should appear in the MSA input.
        header_separator: separator between species name and accession.

    Returns:
        FASTA text, wrapped at ``config.FASTA_LINE_WIDTH`` characters.
    """
    lines: list[str] = []
    for record in records:
        header = record.description or record.id
        if record.organism:
            header = f"{record.organism}{header_separator}{record.accession}"
        else:
            header = f"{record.accession}"
        lines.append(f">{header}")
        lines.extend(
            textwrap.wrap(record.sequence, width=config.FASTA_LINE_WIDTH)
        )
    return "\n".join(lines) + "\n"
