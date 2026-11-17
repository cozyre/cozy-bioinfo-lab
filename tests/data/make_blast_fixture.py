#!/usr/bin/env python3.11
"""Generate tests/data/blast_tp53_five_species.xml — a *synthetic* BLAST XML.

NCBI's BLAST queue was unreachable when the fixture was first recorded
(every job sat in WAITING for >2.5 h), so this generator builds a realistic
stand-in instead:

- the query and all subject sequences are REAL (UniProt P04637 and the
  RefSeq entries fetched from NCBI Entrez on 2026-11-17, using the project
  owner's contact email — deliberately not stored in this repository);
- every HSP is a REAL pairwise alignment computed with Biopython's
  PairwiseAligner using blastp's default scoring (BLOSUM62, gap open 11,
  gap extend 1), so identities, positives, gaps, coordinates and aligned
  strings are scientifically consistent;
- only the XML *document* is synthesized: bit scores are derived from the
  raw score with blastp's typical Karlin-Altschul parameters and E-values
  follow from the bit score. Hit order and accessions mimic a real
  five-organism TP53 blastp run against refseq_protein.

Replace this file with a recorded real XML (record_blast.sh in the session
notes / PLAN phase 2 smoke test) whenever NCBI is reachable again.
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path

from Bio import Align
from Bio.Align import substitution_matrices

HERE = Path(__file__).parent
LAMBDA = 0.318  # typical blastp BLOSUM62 Karlin-Altschul parameters
KAPPA = 0.134
DB_LETTERS = 2.4e8  # approximate refseq_protein size in residues

# Candidates: (accession.version, NCBI uid, title). All fetched live from
# NCBI Entrez on 2026-11-17 and stored next to this script as *.fa.
CANDIDATES = [
    ("NP_001394193.1", "cellular tumor antigen p53 isoform a [Homo sapiens]"),
    ("NP_001394197.1", "cellular tumor antigen p53 isoform b [Homo sapiens]"),
    ("XP_016786959.2", "cellular tumor antigen p53 isoform X1 [Pan troglodytes]"),
    ("XP_001172077.2", "cellular tumor antigen p53 isoform X2 [Pan troglodytes]"),
    ("NP_035770.2", "cellular tumor antigen p53 isoform a [Mus musculus]"),
    ("NP_001120705.1", "cellular tumor antigen p53 isoform b [Mus musculus]"),
    ("XP_006533220.1", "cellular tumor antigen p53 isoform X1 [Mus musculus]"),
    ("NP_112251.2", "cellular tumor antigen p53 isoform a [Rattus norvegicus]"),
    ("NP_001416925.1", "cellular tumor antigen p53 isoform b [Rattus norvegicus]"),
    ("NP_776626.1", "cellular tumor antigen p53 [Bos taurus]"),
]


def read_query() -> str:
    data = json.loads((HERE / "uniprot_P04637.json").read_text())
    return data["sequence"]["value"]


def read_subject(accession: str) -> str:
    text = (HERE / f"{accession}.fa").read_text()
    return "".join(text.splitlines()[1:]).upper()


def blosum62() -> substitution_matrices.Array:
    return substitution_matrices.load("BLOSUM62")


def align_pair(query: str, subject: str, matrix) -> tuple[str, str, str, float]:
    """Best local alignment of query vs subject, blastp-style."""
    aligner = Align.PairwiseAligner(
        mode="local",
        substitution_matrix=matrix,
        open_gap_score=-11,
        extend_gap_score=-1,
    )
    alignments = aligner.align(query, subject)
    best = alignments[0]
    qseq = str(best[0])
    hseq = str(best[1])
    # BLAST reports no leading/trailing gap columns; trim them.
    start = 0
    end = len(qseq)
    while start < end and (qseq[start] == "-" or hseq[start] == "-"):
        start += 1
    while end > start and (qseq[end - 1] == "-" or hseq[end - 1] == "-"):
        end -= 1
    qseq, hseq = qseq[start:end], hseq[start:end]
    midline = ""
    for qa, ha in zip(qseq, hseq, strict=True):
        if qa == "-":
            midline += " "
        elif qa == ha:
            midline += qa
        elif ha != "-" and matrix[qa, ha] > 0:
            midline += "+"
        else:
            midline += " "
    raw = aligner.score(query, subject)
    return qseq, hseq, midline, float(raw)


def hsp_xml(query: str, subject: str) -> str:
    qseq, hseq, midline, raw = align_pair(query, subject, blosum62())
    bit = (LAMBDA * raw - math.log(KAPPA)) / math.log(2)
    evalue = 393 * DB_LETTERS * 2.0 ** (-bit)
    align_len = len(qseq)
    # 1-based coordinates: first and last non-gap characters. The end
    # coordinate is first + (number of non-gap residues) - 1, which stays
    # correct even when the alignment starts after residue 1.
    q_first = next(i for i, ch in enumerate(qseq) if ch != "-") + 1
    h_first = next(i for i, ch in enumerate(hseq) if ch != "-") + 1
    q_last = q_first + len(qseq.replace("-", "")) - 1
    h_last = h_first + len(hseq.replace("-", "")) - 1
    identity = sum(1 for a, b in zip(qseq, hseq, strict=True) if a == b)
    positive = sum(
        1
        for a, b in zip(qseq, hseq, strict=True)
        if a == b or (a != "-" and b != "-" and blosum62()[a, b] > 0)
    )
    gaps = qseq.count("-") + hseq.count("-")
    evalue_str = "0.0" if evalue < 1e-180 else f"{evalue:.2e}".replace("e-0", "e-")
    return f"""        <Hsp>
          <Hsp_num>1</Hsp_num>
          <Hsp_bit-score>{bit:.3f}</Hsp_bit-score>
          <Hsp_score>{int(raw)}</Hsp_score>
          <Hsp_evalue>{evalue_str}</Hsp_evalue>
          <Hsp_query-from>{q_first}</Hsp_query-from>
          <Hsp_query-to>{q_last}</Hsp_query-to>
          <Hsp_hit-from>{h_first}</Hsp_hit-from>
          <Hsp_hit-to>{h_last}</Hsp_hit-to>
          <Hsp_query-frame>0</Hsp_query-frame>
          <Hsp_hit-frame>0</Hsp_hit-frame>
          <Hsp_identity>{identity}</Hsp_identity>
          <Hsp_positive>{positive}</Hsp_positive>
          <Hsp_gaps>{gaps}</Hsp_gaps>
          <Hsp_align-len>{align_len}</Hsp_align-len>
          <Hsp_qseq>{qseq}</Hsp_qseq>
          <Hsp_hseq>{hseq}</Hsp_hseq>
          <Hsp_midline>{midline}</Hsp_midline>
        </Hsp>"""


def main() -> None:
    query = read_query()
    matrix = blosum62()

    hits: list[tuple[float, str]] = []
    for accession, title in CANDIDATES:
        subject = read_subject(accession)
        qseq, hseq, midline, raw = align_pair(query, subject, matrix)
        bit = (LAMBDA * raw - math.log(KAPPA)) / math.log(2)
        hits.append((bit, accession, title, subject, qseq, hseq, midline, raw))
    hits.sort(key=lambda item: -item[0])

    def hit_block(number: int, entry) -> str:
        bit, accession, title, subject, *_ = entry
        hsps = hsp_xml(query, subject)
        return f"""      <Hit>
        <Hit_num>{number}</Hit_num>
        <Hit_id>ref|{accession}|</Hit_id>
        <Hit_def>{title}</Hit_def>
        <Hit_accession>{accession}</Hit_accession>
        <Hit_len>{len(subject)}</Hit_len>
        <Hit_hsps>
{hsps}
        </Hit_hsps>
      </Hit>"""

    hit_xml = "\n".join(hit_block(i + 1, entry) for i, entry in enumerate(hits))

    xml = f"""<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastp</BlastOutput_program>
  <BlastOutput_version>BLASTP 2.17.0+</BlastOutput_version>
  <BlastOutput_reference>Zheng-Brandt et al. (2025), NCBI</BlastOutput_reference>
  <BlastOutput_db>refseq_protein</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>P04637.4 RecName: Full=Cellular tumor antigen p53; AltName: Full=Tumor suppressor p53 [Homo sapiens]</BlastOutput_query-def>
  <BlastOutput_query-len>393</BlastOutput_query-len>
  <BlastOutput_param>
    <Parameters>
      <Parameters_matrix>BLOSUM62</Parameters_matrix>
      <Parameters_expect>1e-10</Parameters_expect>
      <Parameters_sc-match>0</Parameters_sc-match>
      <Parameters_sc-mismatch>0</Parameters_sc-mismatch>
      <Parameters_gap-open>11</Parameters_gap-open>
      <Parameters_gap-extend>1</Parameters_gap-extend>
      <Parameters_filter>F</Parameters_filter>
    </Parameters>
  </BlastOutput_param>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-ID>Query_1</Iteration_query-ID>
      <Iteration_query-def>P04637.4 RecName: Full=Cellular tumor antigen p53; AltName: Full=Tumor suppressor p53 [Homo sapiens]</Iteration_query-def>
      <Iteration_query-len>393</Iteration_query-len>
      <Iteration_hits>
{hit_xml}
      </Iteration_hits>
      <Iteration_stat>
        <Statistics>
          <Statistics_db-num>25925407</Statistics_db-num>
          <Statistics_db-len>24318080388</Statistics_db-len>
          <Statistics_hsp-len>0</Statistics_hsp-len>
          <Statistics_eff-space>4453777295667</Statistics_eff-space>
          <Statistics_kappa>0.041</Statistics_kappa>
          <Statistics_lambda>0.267</Statistics_lambda>
          <Statistics_entropy>0.14</Statistics_entropy>
        </Statistics>
      </Iteration_stat>
      <Iteration_message></Iteration_message>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""

    out = HERE / "blast_tp53_five_species.xml"
    out.write_text(xml)
    print("wrote", out, len(xml), "bytes")
    for i, entry in enumerate(hits, 1):
        bit, accession, title, *_ = entry
        print(f"{i:2d}. {accession}  bit={bit:7.1f}  {title}")


if __name__ == "__main__":
    main()
