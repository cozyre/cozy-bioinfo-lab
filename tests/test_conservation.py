"""Conservation tests: hand-computed values plus fixture sanity checks (§13)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from tp53 import config
from tp53.conservation import (
    SUMMARY_FULLY_CONSERVED_COUNT,
    SUMMARY_LONGEST_REGION,
    SUMMARY_MEAN_IDENTITY,
    SUMMARY_REGION_COUNT,
    classify_columns,
    column_identity,
    compute_conservation,
    find_regions,
    pairwise_identity,
    summarise_features,
)
from tp53.constants import FEATURE_SUMMARY_COLUMNS, OUTSIDE_FEATURES_LABEL
from tp53.fasta_io import parse_fasta
from tp53.models import Alignment, DomainFeature
from tp53.msa import parse_alignment, trim_all_gap_columns, validate_alignment

DATA_DIR = Path(__file__).parent / "data"


def make_alignment(*rows: str) -> Alignment:
    """Build a small alignment through the real parser (so col_to_ref is set)."""
    identifiers = [f"seq{i}" for i in range(len(rows))]
    text = "".join(f">{name}\n{row}\n" for name, row in zip(identifiers, rows, strict=True))
    aln = parse_alignment(text, identifiers[0], identifiers)
    assert validate_alignment(
        aln, identifiers, [row.replace("-", "") for row in rows]
    ) == []
    return aln


class TestColumnIdentity:
    def test_hand_computed_small_alignment(self) -> None:
        aln = make_alignment(
            "ACDEFGHIKL",
            "ACDEFGHIKL",
            "ACDEFGHIKL",
            "ACDEFGHIKX",  # ambiguous code is excluded from the count
            "ACDEF-HIKL",  # gap stays in the denominator
        )
        identity, gaps = column_identity(aln)
        assert identity == pytest.approx(
            [100, 100, 100, 100, 100, 80, 100, 100, 100, 80]
        )
        assert gaps == pytest.approx([0, 0, 0, 0, 0, 20, 0, 0, 0, 0])
        assert classify_columns(identity, gaps) == [
            "fully conserved",
            "fully conserved",
            "fully conserved",
            "fully conserved",
            "fully conserved",
            "high",
            "fully conserved",
            "fully conserved",
            "fully conserved",
            "high",
        ]

    def test_gap_column_never_reaches_100(self) -> None:
        aln = make_alignment(
            "AAAAA", "A----", "A----", "A----", "A----"
        )
        identity, gaps = column_identity(aln)
        assert identity[1:] == pytest.approx([20, 20, 20, 20])
        assert gaps[1:] == pytest.approx([80, 80, 80, 80])
        assert classify_columns(identity, gaps)[1:] == ["gapped"] * 4

    def test_identical_sequences(self) -> None:
        row = "ACDEFGHIKL"
        aln = make_alignment(*([row] * 5))
        identity, gaps = column_identity(aln)
        assert identity == [100.0] * len(row)
        assert gaps == [0.0] * len(row)
        result = compute_conservation(aln)
        assert result.summary[SUMMARY_MEAN_IDENTITY] == 100.0
        assert result.summary[SUMMARY_FULLY_CONSERVED_COUNT] == len(row)
        assert result.summary[SUMMARY_REGION_COUNT] == 1
        assert result.summary[SUMMARY_LONGEST_REGION] == len(row)
        assert (result.pairwise.to_numpy() == 100.0).all()

    def test_empty_alignment(self) -> None:
        assert column_identity(Alignment([], [], 0, 0)) == ([], [])


class TestRegions:
    def test_min_length_and_threshold(self) -> None:
        aln = make_alignment(
            "ACDEFGHIKL",
            "ACDEFGHIKL",
            "ACDEFGHIKL",
            "ACDEFGHIKX",
            "ACDEF-HIKL",
        )
        identity, gaps = column_identity(aln)
        default_regions = find_regions(aln, identity, gaps)
        assert [(r.aln_start, r.aln_end) for r in default_regions] == [(1, 5)]
        assert default_regions[0].mean_identity == 100.0

        relaxed = config.ConservationParams(
            region_identity_threshold=80.0, min_region_length=3
        )
        regions = find_regions(aln, identity, gaps, relaxed)
        assert [(r.aln_start, r.aln_end) for r in regions] == [(1, 5), (7, 10)]
        assert all(region.length >= 3 for region in regions)

    def test_reference_coordinates_skip_reference_gaps(self) -> None:
        aln = make_alignment(
            "AC-DEFGH",  # reference gaps at alignment column 3
            "ACXDEFGH",
            "ACXDEFGH",
            "ACXDEFGH",
            "ACXDEFGH",
        )
        assert aln.col_to_ref == [1, 2, None, 3, 4, 5, 6, 7]
        identity, gaps = column_identity(aln)
        regions = find_regions(aln, identity, gaps)
        assert len(regions) == 1
        region = regions[0]
        assert (region.aln_start, region.aln_end) == (4, 8)
        assert (region.ref_start, region.ref_end) == (3, 7)

    def test_no_region_below_threshold(self) -> None:
        aln = make_alignment("AAAAA", "AAAAT", "AAATA", "AATAA", "ATAAA")
        identity, gaps = column_identity(aln)
        assert find_regions(aln, identity, gaps) == []


class TestPairwiseIdentity:
    def test_hand_computed_matrix(self) -> None:
        aln = make_alignment("ACGT", "ACGA", "AC-T")
        matrix = pairwise_identity(aln)
        assert matrix.loc["seq0", "seq1"] == pytest.approx(75.0)  # 3/4
        assert matrix.loc["seq0", "seq2"] == pytest.approx(100.0)  # 3/3
        assert matrix.loc["seq1", "seq2"] == pytest.approx(200 / 3)  # 2/3
        assert matrix.loc["seq2", "seq0"] == matrix.loc["seq0", "seq2"]
        assert (matrix.to_numpy().diagonal() == 100.0).all()


class TestFeatureSummary:
    def test_reference_numbering_and_outside_row(self) -> None:
        aln = make_alignment(
            "ACDEFGHIKL",
            "ACDEFGHIKL",
            "ACDEFGHIKL",
            "ACDEFGHIKL",
            "ACDEFGHIKL",
        )
        identity, _ = column_identity(aln)
        features = [
            DomainFeature("Region", "First three", 1, 3),
            DomainFeature("DNA binding", "Core", 5, 7),
        ]
        summary = summarise_features(aln, identity, features)
        assert list(summary.columns) == list(FEATURE_SUMMARY_COLUMNS)
        assert summary.loc[0, "columns"] == 3
        assert summary.loc[0, "mean_identity"] == 100.0
        assert summary.loc[1, "start"] == 5 and summary.loc[1, "end"] == 7
        assert summary.loc[2, "feature"] == OUTSIDE_FEATURES_LABEL
        assert summary.loc[2, "columns"] == 4  # positions 4, 8, 9, 10

    def test_empty_features(self) -> None:
        aln = make_alignment("AAAAA", "AAAAA")
        identity, _ = column_identity(aln)
        summary = summarise_features(aln, identity, [])
        assert len(summary) == 1
        assert summary.loc[0, "feature"] == OUTSIDE_FEATURES_LABEL
        assert summary.loc[0, "columns"] == 5


@pytest.fixture(scope="module")
def fixture_alignment() -> Alignment:
    """The recorded five-species Clustal Omega alignment."""
    inputs = parse_fasta((DATA_DIR / "msa_input.fasta").read_text())
    ids = [record.id for record in inputs]
    aln = parse_alignment(
        (DATA_DIR / "msa_tp53_five_species.fasta").read_text(), ids[0], ids
    )
    assert validate_alignment(
        aln, ids, [record.sequence for record in inputs]
    ) == []
    trimmed, removed = trim_all_gap_columns(aln)
    assert removed == 0
    return trimmed


class TestRecordedAlignment:
    """AC-3, AC-4 and the scientific sanity check on the real MSA."""

    def test_identity_shape_and_bounds(self, fixture_alignment: Alignment) -> None:
        result = compute_conservation(fixture_alignment)
        assert len(result.identity) == fixture_alignment.length
        assert len(result.gap_fraction) == fixture_alignment.length
        assert len(result.classes) == fixture_alignment.length
        assert all(0.0 <= value <= 100.0 for value in result.identity)
        for percent, gaps, label in zip(
            result.identity, result.gap_fraction, result.classes, strict=True
        ):
            if label == "fully conserved":
                assert percent == 100.0 and gaps == 0.0

    def test_regions_respect_rules(self, fixture_alignment: Alignment) -> None:
        result = compute_conservation(fixture_alignment)
        assert result.regions
        for region in result.regions:
            assert region.length >= config.ConservationParams().min_region_length
            assert all(
                identity >= config.ConservationParams().region_identity_threshold
                for identity in result.identity[
                    region.aln_start - 1 : region.aln_end
                ]
            )
            assert 1 <= region.ref_start <= region.ref_end <= 393

    def test_feature_overlay_uses_reference_coordinates(
        self, fixture_alignment: Alignment
    ) -> None:
        payload = json.loads((DATA_DIR / "uniprot_P04637.json").read_text())
        dna_binding = next(
            feature
            for feature in payload["features"]
            if feature["type"] == "DNA binding"
        )
        feature = DomainFeature(
            type="DNA binding",
            name="DNA binding",
            start=dna_binding["location"]["start"]["value"],
            end=dna_binding["location"]["end"]["value"],
        )
        result = compute_conservation(fixture_alignment, [feature])
        row = result.feature_summary.iloc[0]
        assert (row["start"], row["end"]) == (102, 292)
        assert result.feature_summary.iloc[-1]["feature"] == OUTSIDE_FEATURES_LABEL
        # Every column assigned to the feature maps to a reference position
        # inside the UniProt coordinate range.
        columns = [
            index
            for index, position in enumerate(fixture_alignment.col_to_ref)
            if position is not None and 102 <= position <= 292
        ]
        assert row["columns"] == len(columns)
        assert row["mean_identity"] > 90.0

        # Scientific sanity check (§13): the DNA-binding core is more conserved
        # than the N-terminal transactivation region.
        n_terminal = compute_conservation(
            fixture_alignment, [DomainFeature("Region", "N-term", 1, 83)]
        ).feature_summary.iloc[0]
        assert row["mean_identity"] > n_terminal["mean_identity"] + 10
