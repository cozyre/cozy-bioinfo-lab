"""Tests for fasta_io (FR-03)."""

from __future__ import annotations

import pytest

from tp53.fasta_io import parse_fasta, parse_single_fasta, to_fasta, validate_protein
from tp53.models import SequenceRecord, ValidationError


def _record(sequence: str, id: str = "test|X1", organism: str = "") -> SequenceRecord:
    return SequenceRecord(
        id=id,
        accession=id.split("|")[-1],
        organism=organism,
        description="test protein",
        sequence=sequence,
    )


class TestParseFasta:
    def test_single_record(self) -> None:
        records = parse_fasta(">P04637 TP53 human\nMEEPQSDPSVE\nPP")
        assert len(records) == 1
        record = records[0]
        assert record.id == "P04637"
        assert record.description == "TP53 human"
        assert record.sequence == "MEEPQSDPSVEPP"

    def test_multi_record_preserves_order(self) -> None:
        records = parse_fasta(">a\nMEEP\n>b\nQSDP\n>c\nVEPP\n")
        assert [record.id for record in records] == ["a", "b", "c"]
        assert [record.sequence for record in records] == ["MEEP", "QSDP", "VEPP"]

    def test_uppercase_and_whitespace(self) -> None:
        records = parse_fasta(">x\nmeep\n sdps\n")
        assert records[0].sequence == "MEEPSDPS"

    def test_gaps_are_kept_for_aligned_input(self) -> None:
        records = parse_fasta(">x aligned\nAC--GT\n")
        assert records[0].sequence == "AC--GT"

    def test_not_fasta_raises(self) -> None:
        with pytest.raises(ValidationError, match="FASTA"):
            parse_fasta("MEEPQSDPSVEPPLSQETFS\n")

    def test_empty_raises(self) -> None:
        with pytest.raises(ValidationError, match="no FASTA records"):
            parse_fasta("   \n")

    def test_single_record_constraint(self) -> None:
        text = ">a\nMEEPQSDPSVEPP\n>b\nQSDPSVEPPLSQ\n"
        with pytest.raises(ValidationError, match="exactly one record.*2"):
            parse_single_fasta(text)

    def test_single_record_ok(self) -> None:
        record = parse_single_fasta(">a\nMEEPQSDPSVEPP\n")
        assert record.id == "a"


class TestValidateProtein:
    def test_valid(self) -> None:
        validate_protein(_record("M" * 393))

    def test_min_length_boundary(self) -> None:
        validate_protein(_record("M" * 20))

    def test_too_short_raises(self) -> None:
        with pytest.raises(ValidationError, match="length 19"):
            validate_protein(_record("M" * 19))

    def test_too_long_raises(self) -> None:
        with pytest.raises(ValidationError, match="5001"):
            validate_protein(_record("M" * 5001))

    def test_invalid_character_named(self) -> None:
        with pytest.raises(ValidationError, match="\\*"):
            validate_protein(_record("MEEPQSDPSVEPP*Q"))

    def test_empty_sequence_raises(self) -> None:
        with pytest.raises(ValidationError, match="empty"):
            validate_protein(_record(""))

    def test_ambiguity_codes_allowed(self) -> None:
        # X B Z J U O are legal amino-acid codes (config.PROTEIN_ALPHABET).
        validate_protein(_record("M" * 13 + "XBZJUO" + "M" * 19))


class TestToFasta:
    def test_species_header_format(self) -> None:
        text = to_fasta([_record("MEEPQSDPSVEPP", organism="Homo sapiens")])
        assert text.startswith(">Homo sapiens|X1\n")

    def test_no_organism_falls_back_to_accession(self) -> None:
        text = to_fasta([_record("MEEPQSDPSVEPP")])
        assert text.startswith(">X1\n")

    def test_line_wrapping(self) -> None:
        text = to_fasta([_record("M" * 393)])
        body = text.splitlines()[1:]
        assert max(len(line) for line in body) <= 60
        assert "".join(body) == "M" * 393

    def test_round_trip(self) -> None:
        record = _record("MEEPQSDPSVEPP", organism="Homo sapiens")
        parsed = parse_fasta(to_fasta([record]))[0]
        assert parsed.sequence == record.sequence
        # Biopython splits FASTA IDs at the first space; the full species
        # header lands in the description. msa.parse_alignment re-matches
        # such truncated IDs with its prefix fallback (§9.4).
        assert parsed.id == "Homo"
        assert parsed.description.startswith("sapiens|X1")
