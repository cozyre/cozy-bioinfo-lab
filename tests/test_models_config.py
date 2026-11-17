"""Tests for models and config defaults (§16 decisions, §7.1)."""

from __future__ import annotations

import pytest

from tp53 import config
from tp53.models import BlastHit, BlastHsp, SequenceRecord


class TestConfigDefaults:
    """The §16 defaults PLAN.md locks in must hold."""

    def test_blast_database(self) -> None:
        assert config.BlastParams().database == "refseq_protein"

    def test_region_threshold_is_100_percent(self) -> None:
        assert config.ConservationParams().region_identity_threshold == 100.0

    def test_min_region_length(self) -> None:
        assert config.ConservationParams().min_region_length == 5

    def test_species_list(self) -> None:
        assert config.SPECIES == (
            "Homo sapiens",
            "Pan troglodytes",
            "Mus musculus",
            "Rattus norvegicus",
            "Bos taurus",
        )

    def test_blast_cutoffs(self) -> None:
        params = config.BlastParams()
        assert params.expect == 1e-10
        assert params.hitlist_size == 10
        selection = config.SelectionParams()
        assert selection.min_query_coverage == pytest.approx(0.70)

    def test_fasta_length_bounds(self) -> None:
        assert config.FASTA_MIN_LEN == 20
        assert config.FASTA_MAX_LEN == 5000


class TestModels:
    def test_best_hsp_is_highest_bit_score(self) -> None:
        def hsp(bits: float) -> BlastHsp:
            return BlastHsp(
                bit_score=bits, evalue=1e-20, identity=90, positives=95, gaps=0,
                align_len=393, query_start=1, query_end=393, hit_start=1,
                hit_end=393, qseq="A" * 393, hseq="A" * 393, midline="A" * 393,
            )

        hit = BlastHit(
            accession="NP_000537.3", title="P53 [Homo sapiens]",
            organism="Homo sapiens", length=393,
            hsps=[hsp(100.0), hsp(890.5), hsp(50.0)],
        )
        assert hit.best_hsp.bit_score == pytest.approx(890.5)

    def test_hsp_percent_identity(self) -> None:
        hsp = BlastHsp(
            bit_score=800, evalue=0, identity=350, positives=360, gaps=3,
            align_len=393, query_start=1, query_end=393, hit_start=1,
            hit_end=393, qseq="", hseq="", midline="",
        )
        assert hsp.percent_identity == pytest.approx(100.0 * 350 / 393)

    def test_ungapped(self) -> None:
        record = SequenceRecord(
            id="x", accession="x", organism="", description="", sequence="AC--GT"
        )
        assert record.ungapped() == "ACGT"
