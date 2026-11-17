"""Tests for blast_parser against a synthetic XML and the recorded fixture."""

from __future__ import annotations

from pathlib import Path

import pytest

from tp53.blast_parser import parse_blast_xml
from tp53.models import BlastResult

DATA_DIR = Path(__file__).parent / "data"
FIXTURE = DATA_DIR / "blast_tp53_five_species.xml"

# A minimal blastp XML in the NCBI XML schema: one hit with one HSP, whose
# values the tests assert exactly. Note: Bio.Blast.NCBIXML requires the text
# to start with the XML declaration and tolerates no root namespace.
MINI_XML = """<?xml version="1.0"?>
<BlastOutput>
  <BlastOutput_program>blastp</BlastOutput_program>
  <BlastOutput_version>BLASTP 2.16.0+</BlastOutput_version>
  <BlastOutput_db>refseq_protein</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>P04637.4 RecName: Full=Cellular tumor antigen p53</BlastOutput_query-def>
  <BlastOutput_query-len>393</BlastOutput_query-len>
  <BlastOutput_param>
    <Parameters>
      <Parameters_matrix>BLOSUM62</Parameters_matrix>
      <Parameters_expect>1e-10</Parameters_expect>
      <Parameters_gap-open>11</Parameters_gap-open>
      <Parameters_gap-extend>1</Parameters_gap-extend>
      <Parameters_filter>F</Parameters_filter>
    </Parameters>
  </BlastOutput_param>
  <BlastOutput_iterations>
  <Iteration>
    <Iteration_iter-num>1</Iteration_iter-num>
    <Iteration_query-def>P04637.4 RecName: Full=Cellular tumor antigen p53</Iteration_query-def>
    <Iteration_query-len>393</Iteration_query-len>
    <Iteration_hits>
      <Hit>
        <Hit_num>1</Hit_num>
        <Hit_id>ref|NP_038672.1|</Hit_id>
        <Hit_def>cellular tumor antigen p53 [Mus musculus]</Hit_def>
        <Hit_accession>NP_038672</Hit_accession>
        <Hit_len>389</Hit_len>
        <Hsp>
          <Hsp_num>1</Hsp_num>
          <Hsp_bit-score>612.123</Hsp_bit-score>
          <Hsp_score>1573</Hsp_score>
          <Hsp_evalue>0.0</Hsp_evalue>
          <Hsp_query-from>1</Hsp_query-from>
          <Hsp_query-to>385</Hsp_query-to>
          <Hsp_hit-from>3</Hsp_hit-from>
          <Hsp_hit-to>387</Hsp_hit-to>
          <Hsp_query-frame>0</Hsp_query-frame>
          <Hsp_identity>300</Hsp_identity>
          <Hsp_positive>340</Hsp_positive>
          <Hsp_gaps>4</Hsp_gaps>
          <Hsp_align-len>385</Hsp_align-len>
          <Hsp_qseq>MEEPQSDPSVEPPLSQETFSDLWKLLPE</Hsp_qseq>
          <Hsp_hseq>MEEPQSDLSIEPPLSQETFSDLWKLLPE</Hsp_hseq>
          <Hsp_midline>MEEPQSD+ +EPPLSQETFSDLWKLLPE</Hsp_midline>
        </Hsp>
      </Hit>
    </Iteration_hits>
    <Iteration_stat>
      <Statistics>
        <Statistics_db-num>12345</Statistics_db-num>
        <Statistics_db-len>45678901</Statistics_db-len>
        <Statistics_hsp-len>0</Statistics_hsp-len>
        <Statistics_eff-space>12345678901</Statistics_eff-space>
        <Statistics_kappa>0.041</Statistics_kappa>
        <Statistics_lambda>0.267</Statistics_lambda>
        <Statistics_entropy>0.14</Statistics_entropy>
      </Statistics>
    </Iteration_stat>
  </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""


class TestMiniXml:
    """Structural parsing, verified against hand-written values."""

    @pytest.fixture()
    def result(self) -> BlastResult:
        return parse_blast_xml(MINI_XML, species="Mus musculus")

    def test_header(self, result: BlastResult) -> None:
        assert result.program == "BLASTP"
        assert result.database == "refseq_protein"
        assert result.query_len == 393
        assert result.species == "Mus musculus"

    def test_hit_fields(self, result: BlastResult) -> None:
        assert len(result.hits) == 1
        hit = result.hits[0]
        assert hit.accession == "NP_038672.1"
        assert hit.organism == "Mus musculus"
        assert hit.length == 389

    def test_hsp_fields(self, result: BlastResult) -> None:
        hsp = result.hits[0].hsps[0]
        assert hsp.bit_score == pytest.approx(612.123)
        assert hsp.evalue == pytest.approx(0.0)
        assert hsp.identity == 300
        assert hsp.positives == 340
        assert hsp.gaps == 4
        assert hsp.align_len == 385
        assert (hsp.query_start, hsp.query_end) == (1, 385)
        assert (hsp.hit_start, hsp.hit_end) == (3, 387)
        assert hsp.qseq.startswith("MEEPQSDPS")
        assert hsp.hseq.startswith("MEEPQSDLS")
        assert hsp.percent_identity == pytest.approx(100.0 * 300 / 385)

    def test_no_result_raises(self) -> None:
        # NCBIXML raises its own error family for text that parses but holds
        # no BLAST record (ValueError / ExpatError); our parser adds a
        # ValueError when the XML is valid but empty of records.
        from xml.parsers.expat import ExpatError

        with pytest.raises((ValueError, ExpatError)):
            parse_blast_xml("<?xml version='1.0'?><foo/>")


class TestRecordedFixture:
    """Assertions against the real recorded NCBI BLAST XML.

    Expected values are cross-checked against the tabular output of the same
    job (recorded alongside the XML).
    """

    @pytest.fixture(scope="class")
    def result(self) -> BlastResult:
        if not FIXTURE.exists():  # pragma: no cover - recorded fixture guard
            pytest.skip("recorded BLAST XML fixture not present")
        return parse_blast_xml(FIXTURE.read_text(), species="five species")

    def test_header(self, result: BlastResult) -> None:
        assert result.program == "blastp"
        assert result.database == "refseq_protein"
        assert result.query_len == 393
        assert result.species == "five species"

    def test_hits_present(self, result: BlastResult) -> None:
        assert 1 <= len(result.hits) <= 10
        organisms = {hit.organism for hit in result.hits if hit.organism}
        assert {"Homo sapiens", "Pan troglodytes", "Mus musculus",
                "Rattus norvegicus", "Bos taurus"} <= organisms

    def test_top_hit_is_human_refseq(self, result: BlastResult) -> None:
        hit = result.hits[0]
        assert hit.accession.startswith("NP_")
        assert hit.organism == "Homo sapiens"
        assert hit.hsps[0].identity > 380
        assert hit.hsps[0].evalue <= 1e-10

    def test_hsps_have_ranges(self, result: BlastResult) -> None:
        for hit in result.hits:
            for hsp in hit.hsps:
                assert 1 <= hsp.query_start <= hsp.query_end <= 393
                assert hsp.align_len == hsp.query_end - hsp.query_start + 1 + hsp.gaps
                assert len(hsp.qseq) == hsp.align_len
                assert len(hsp.hseq) == hsp.align_len
                assert len(hsp.midline) == hsp.align_len
