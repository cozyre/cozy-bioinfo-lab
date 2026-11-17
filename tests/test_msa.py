"""Tests for msa: parsing, validation, trimming, coordinate map (§9.4, §9.7)."""

from __future__ import annotations

import pandas as pd
import pytest

from tp53.msa import (
    alignment_statistics,
    parse_alignment,
    trim_all_gap_columns,
    validate_alignment,
)
from tp53.models import Alignment, AlignmentValidationError

# 0-indexed columns:   0  1  2  3  4  5
# Query (ref) row:     A  C  -  G  T  -
RAW_ALN = ">Query|P04637\nAC-GT-\n>Mouse|P02340\nAC-GA-\n>Rat|P63059\nACGGT-\n"
REF_ID = "Query|P04637"
EXPECTED_IDS = ["Query|P04637", "Mouse|P02340", "Rat|P63059"]
EXPECTED_SEQS = ["ACGT", "ACGA", "ACGGT"]


class TestParseAlignment:
    def test_basic_fields(self) -> None:
        alignment = parse_alignment(RAW_ALN, REF_ID)
        assert alignment.ids == EXPECTED_IDS
        assert alignment.length == 6
        assert alignment.ref_index == 0

    def test_col_to_ref_gaps_map_to_none(self) -> None:
        alignment = parse_alignment(RAW_ALN, REF_ID)
        assert alignment.col_to_ref == [1, 2, None, 3, 4, None]

    def test_col_to_ref_reference_not_first_row(self) -> None:
        alignment = parse_alignment(RAW_ALN, "Rat|P63059")
        assert alignment.ref_index == 2
        assert alignment.col_to_ref == [1, 2, 3, 4, 5, None]

    def test_truncated_reference_id_found(self) -> None:
        # Clustal Omega may truncate long headers; prefix match must win.
        alignment = parse_alignment(RAW_ALN, "Query")
        assert alignment.ref_index == 0

    def test_unknown_reference_raises(self) -> None:
        with pytest.raises(AlignmentValidationError, match="not found"):
            parse_alignment(RAW_ALN, "Missing")

    def test_id_remap_by_order(self) -> None:
        # Service rewrote the IDs: rows are remapped to the expected IDs.
        rewritten = ">seq1\nAC-GT-\n>seq2\nAC-GA-\n>seq3\nACGGT-\n"
        alignment = parse_alignment(rewritten, "Query|P04637", expected_ids=EXPECTED_IDS)
        assert alignment.ids == EXPECTED_IDS
        assert alignment.ref_index == 0
        assert alignment.rows[0] == "AC-GT-"

    def test_empty_raises(self) -> None:
        with pytest.raises(AlignmentValidationError, match="no FASTA records"):
            parse_alignment("", "Query")

    def test_parse_is_case_insensitive(self) -> None:
        alignment = parse_alignment(RAW_ALN.lower(), REF_ID.lower())
        assert alignment.rows[0] == "AC-GT-"


class TestTrimAllGapColumns:
    def test_trim_count_and_columns(self) -> None:
        alignment = parse_alignment(RAW_ALN, REF_ID)
        trimmed, removed = trim_all_gap_columns(alignment)
        assert removed == 1
        assert trimmed.length == 5
        assert trimmed.rows == ["AC-GT", "AC-GA", "ACGGT"]

    def test_trim_rebuilds_col_to_ref(self) -> None:
        alignment = parse_alignment(RAW_ALN, REF_ID)
        trimmed, _ = trim_all_gap_columns(alignment)
        assert trimmed.col_to_ref == [1, 2, None, 3, 4]

    def test_no_all_gap_columns(self) -> None:
        alignment = parse_alignment(">a\nACGT\n>b\nACGA\n", "a")
        trimmed, removed = trim_all_gap_columns(alignment)
        assert removed == 0
        assert trimmed.rows == ["ACGT", "ACGA"]

    def test_entirely_gap_column_set(self) -> None:
        alignment = parse_alignment(">a\nA-N\n>b\nA-N\n", "a")
        trimmed, removed = trim_all_gap_columns(alignment)
        assert removed == 1
        assert trimmed.rows == ["AN", "AN"]


class TestValidateAlignment:
    def test_valid_returns_no_warnings(self) -> None:
        alignment = parse_alignment(RAW_ALN, REF_ID)
        assert validate_alignment(alignment, EXPECTED_IDS, EXPECTED_SEQS) == []

    def test_row_count_mismatch(self) -> None:
        alignment = parse_alignment(">a\nACGT\n", "a")
        with pytest.raises(AlignmentValidationError, match="rows"):
            validate_alignment(alignment, EXPECTED_IDS)

    def test_row_length_mismatch_names_sequence(self) -> None:
        broken = ">Query|P04637\nAC-GT-\n>Mouse|P02340\nAC-GA\n>Rat|P63059\nACGGT-\n"
        alignment = parse_alignment(broken, REF_ID)
        with pytest.raises(AlignmentValidationError, match="Mouse\\|P02340.*length 5"):
            validate_alignment(alignment, EXPECTED_IDS)

    def test_invalid_characters_named(self) -> None:
        broken = ">Query|P04637\nAC*GT-\n>Mouse|P02340\nAC-GA-\n>Rat|P63059\nACGGT-\n"
        alignment = parse_alignment(broken, REF_ID)
        with pytest.raises(AlignmentValidationError, match="Query\\|P04637.*\\*"):
            validate_alignment(alignment, EXPECTED_IDS)

    def test_degap_mismatch_names_sequence(self) -> None:
        # The Mouse row aligns to a slightly different input sequence.
        alignment = parse_alignment(RAW_ALN, REF_ID)
        wrong_seqs = ["ACGT", "ACGT", "ACGGT"]
        with pytest.raises(AlignmentValidationError, match="Mouse\\|P02340"):
            validate_alignment(alignment, EXPECTED_IDS, wrong_seqs)

    def test_degap_check_accepts_originals(self) -> None:
        alignment = parse_alignment(RAW_ALN, REF_ID)
        assert validate_alignment(alignment, EXPECTED_IDS, EXPECTED_SEQS) == []

    def test_id_mismatch_warns(self) -> None:
        rewritten = ">seq1\nAC-GT-\n>seq2\nAC-GA-\n>seq3\nACGGT-\n"
        alignment = parse_alignment(rewritten, "Query|P04637", expected_ids=EXPECTED_IDS)
        # parse_alignment already remapped, so validation is clean here;
        # a mismatched alignment passed directly must produce a warning.
        stale = Alignment(
            ids=["seq1", "seq2", "seq3"],
            rows=["AC-GT-", "AC-GA-", "ACGGT-"],
            length=6,
            ref_index=0,
            col_to_ref=[],
        )
        warnings = validate_alignment(stale, EXPECTED_IDS)
        assert warnings and "order" in warnings[0]

    def test_case_insensitive_id_check(self) -> None:
        alignment = parse_alignment(RAW_ALN, REF_ID)
        alignment.ids = [identifier.lower() for identifier in alignment.ids]
        assert validate_alignment(alignment, EXPECTED_IDS) == []


class TestAlignmentStatistics:
    def test_columns_and_values(self) -> None:
        alignment = parse_alignment(RAW_ALN, REF_ID)
        stats = alignment_statistics(alignment)
        assert list(stats.columns) == ["id", "length", "ungapped_length", "gap_pct"]
        first = stats.iloc[0]
        assert first["id"] == "Query|P04637"
        assert first["length"] == 6
        assert first["ungapped_length"] == 4
        assert first["gap_pct"] == pytest.approx(100.0 * 2 / 6, abs=0.1)

    def test_returns_dataframe(self) -> None:
        assert isinstance(alignment_statistics(parse_alignment(RAW_ALN, REF_ID)), pd.DataFrame)
