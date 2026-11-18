"""Deterministic best-hit selection tests (§9.2, FR-11)."""

from __future__ import annotations

from pathlib import Path

import pytest

from tp53 import config
from tp53.blast_parser import parse_blast_xml
from tp53.models import BlastHit, BlastHsp, BlastResult
from tp53.selection import select_best_hit

DATA_DIR = Path(__file__).parent / "data"


def make_hit(
    accession: str,
    *,
    bit_score: float = 100.0,
    identity: int = 300,
    align_len: int = 393,
    evalue: float = 1e-30,
    query_start: int = 1,
    query_end: int = 393,
    organism: str = "Mus musculus",
) -> BlastHit:
    """Build a hit with one HSP whose numbers are easy to reason about."""
    hsp = BlastHsp(
        bit_score=bit_score,
        evalue=evalue,
        identity=identity,
        positives=identity + 20,
        gaps=0,
        align_len=align_len,
        query_start=query_start,
        query_end=query_end,
        hit_start=1,
        hit_end=align_len,
        qseq="",
        hseq="",
        midline="",
    )
    return BlastHit(
        accession=accession,
        title=f"p53 [{organism}]",
        organism=organism,
        length=align_len,
        hsps=[hsp],
    )


def make_result(*hits: BlastHit, species: str = "Mus musculus") -> BlastResult:
    """Wrap hits in a result for the default 393-residue query."""
    return BlastResult(
        query_id="P04637",
        query_len=393,
        database="refseq_protein",
        program="BLASTP",
        hits=list(hits),
        species=species,
    )


class TestFilters:
    def test_evalue_filter(self) -> None:
        poor = make_hit("NP_000000.1", evalue=1e-5)
        good = make_hit("NP_111111.1", bit_score=90, evalue=1e-50)
        selected = select_best_hit(make_result(poor, good))
        assert selected is not None
        assert selected.hit.accession == "NP_111111.1"
        assert [hit.accession for hit in selected.candidates] == ["NP_111111.1"]

    def test_coverage_filter(self) -> None:
        fragment = make_hit("NP_000000.1", bit_score=500, query_end=200)
        good = make_hit("NP_111111.1", bit_score=100)
        selected = select_best_hit(make_result(fragment, good))
        assert selected is not None
        assert selected.hit.accession == "NP_111111.1"

    def test_custom_thresholds(self) -> None:
        params = config.SelectionParams(evalue_cutoff=1e-3, min_query_coverage=0.5)
        fragment = make_hit("NP_000000.1", bit_score=500, query_end=250)
        selected = select_best_hit(make_result(fragment), params)
        assert selected is not None
        assert selected.query_coverage == pytest.approx(250 / 393)

    def test_no_hits_warns_and_returns_none(self) -> None:
        with pytest.warns(UserWarning, match="No BLAST hit for Mus musculus"):
            assert select_best_hit(make_result()) is None

    def test_all_hits_filtered_warns(self) -> None:
        result = make_result(make_hit("NP_000000.1", evalue=1.0))
        with pytest.warns(UserWarning, match="passed"):
            assert select_best_hit(result) is None

    def test_hit_without_hsps_is_skipped(self) -> None:
        empty = BlastHit("NP_000000.1", "empty", "Mus musculus", 393)
        with pytest.warns(UserWarning):
            assert select_best_hit(make_result(empty)) is None


class TestRanking:
    def test_bit_score_wins_over_input_order(self) -> None:
        first = make_hit("NP_000001.1", bit_score=100)
        second = make_hit("NP_000002.1", bit_score=200)
        selected = select_best_hit(make_result(first, second))
        assert selected is not None
        assert selected.hit.accession == "NP_000002.1"

    def test_identity_tie_break(self) -> None:
        lower = make_hit("NP_000001.1", identity=250)
        higher = make_hit("NP_000002.1", identity=300)
        selected = select_best_hit(make_result(lower, higher))
        assert selected is not None
        assert selected.hit.accession == "NP_000002.1"
        assert "percent identity" in selected.reason

    def test_np_preferred_over_xp_on_full_tie(self) -> None:
        xp = make_hit("XP_000001.1")
        np_hit = make_hit("NP_999999.1")
        selected = select_best_hit(make_result(xp, np_hit))
        assert selected is not None
        assert selected.hit.accession == "NP_999999.1"
        assert "NP_ over XP_" in selected.reason

    def test_accession_alphabetical_on_full_tie(self) -> None:
        later = make_hit("NP_999999.1")
        earlier = make_hit("NP_000001.1")
        selected = select_best_hit(make_result(later, earlier))
        assert selected is not None
        assert selected.hit.accession == "NP_000001.1"
        assert "accession alphabetically" in selected.reason

    def test_order_independent(self) -> None:
        a = make_hit("NP_000001.1", bit_score=100, identity=300)
        b = make_hit("NP_000002.1", bit_score=100, identity=300)
        c = make_hit("NP_000003.1", bit_score=100, identity=300)
        first = select_best_hit(make_result(c, a, b))
        second = select_best_hit(make_result(b, c, a))
        assert first is not None and second is not None
        assert first.hit.accession == second.hit.accession == "NP_000001.1"

    def test_reason_contains_score_identity_coverage(self) -> None:
        selected = select_best_hit(make_result(make_hit("NP_000001.1")))
        assert selected is not None
        assert "bit score 100.0" in selected.reason
        assert "identity 76.3%" in selected.reason
        assert "query coverage 100.0%" in selected.reason


@pytest.fixture(scope="module")
def per_species_results() -> dict[str, BlastResult]:
    """Split the recorded combined XML into the per-species jobs of §9.1."""
    combined = parse_blast_xml(
        (DATA_DIR / "blast_tp53_five_species.xml").read_text()
    )
    return {
        species: BlastResult(
            query_id=combined.query_id,
            query_len=combined.query_len,
            database=combined.database,
            program=combined.program,
            hits=[hit for hit in combined.hits if hit.organism == species],
            species=species,
        )
        for species in config.SPECIES
    }


class TestRecordedFixture:
    """AC-1: one viable candidate per species in the recorded BLAST XML."""

    def test_one_hit_per_species(
        self, per_species_results: dict[str, BlastResult]
    ) -> None:
        expected = {
            "Homo sapiens": "NP_001394193.1",
            "Pan troglodytes": "XP_001172077.2",
            "Mus musculus": "NP_035770.2",  # NP-over-XP tie break
            "Rattus norvegicus": "NP_112251.2",
            "Bos taurus": "NP_776626.1",
        }
        for species in config.SPECIES:
            selected = select_best_hit(per_species_results[species])
            assert selected is not None, species
            assert selected.species == species
            assert selected.hit.accession == expected[species]
            assert selected.hit.organism == species

    def test_candidates_are_ranked(
        self, per_species_results: dict[str, BlastResult]
    ) -> None:
        selected = select_best_hit(per_species_results["Homo sapiens"])
        assert selected is not None
        scores = [hit.best_hsp.bit_score for hit in selected.candidates]
        assert scores == sorted(scores, reverse=True)
