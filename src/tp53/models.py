"""Data models (§7.1) and exceptions.

All dataclasses from design doc section 7.1, verbatim. Dataclasses are the
only data structures passed between layers.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from typing import Any, Literal

import pandas as pd

from tp53 import config


class ValidationError(Exception):
    """Input failed validation (FR-03); message names the problem."""


class AlignmentValidationError(ValidationError):
    """MSA failed validation (§9.4); message names the offending sequence."""


@dataclass
class SequenceRecord:
    """One protein sequence with its provenance."""

    id: str
    accession: str
    organism: str
    description: str
    sequence: str

    def ungapped(self) -> str:
        """Sequence with all gap characters removed (§9.4)."""
        return self.sequence.replace("-", "")


@dataclass
class DomainFeature:
    """One UniProt feature in 1-based reference coordinates (§9.9)."""

    type: str
    name: str
    start: int
    end: int


@dataclass
class QueryInfo:
    """The query sequence plus its UniProt feature annotations (FR-06)."""

    record: SequenceRecord
    features: list[DomainFeature]


@dataclass
class BlastHsp:
    """One high-scoring segment pair of a BLAST hit (FR-09)."""

    bit_score: float
    evalue: float
    identity: int
    positives: int
    gaps: int
    align_len: int
    query_start: int
    query_end: int
    hit_start: int
    hit_end: int
    qseq: str
    hseq: str
    midline: str

    @property
    def percent_identity(self) -> float:
        """Identity of this HSP as a percentage."""
        return 100.0 * self.identity / self.align_len if self.align_len else 0.0

    @property
    def query_cover(self) -> float:
        """Aligned query length / HSP alignment length, as a fraction."""
        span = self.query_end - self.query_start + 1
        return span / self.align_len if self.align_len else 0.0


@dataclass
class BlastHit:
    """One subject sequence of a BLAST result (FR-09)."""

    accession: str
    title: str
    organism: str
    length: int
    hsps: list[BlastHsp] = field(default_factory=list)

    @property
    def best_hsp(self) -> BlastHsp:
        """HSP with the highest bit score."""
        return max(self.hsps, key=lambda h: h.bit_score)


@dataclass
class BlastResult:
    """Parsed BLAST XML for one per-species job (§7.1)."""

    query_id: str
    query_len: int
    database: str
    program: str
    hits: list[BlastHit] = field(default_factory=list)
    species: str = ""


@dataclass
class SelectedHit:
    """The best hit for one species, with the selection rationale (§9.2)."""

    species: str
    hit: BlastHit
    query_coverage: float
    percent_identity: float
    reason: str
    candidates: list[BlastHit] = field(default_factory=list)


@dataclass
class Alignment:
    """A parsed multiple sequence alignment (§9.4, §9.7)."""

    ids: list[str]
    rows: list[str]
    length: int
    ref_index: int
    col_to_ref: list[int | None] = field(default_factory=list)

    def ungapped_row(self, index: int) -> str:
        """Row with all gap characters removed (§9.4)."""
        return self.rows[index].replace("-", "")


@dataclass
class ConservedRegion:
    """A run of conserved alignment columns (§9.6)."""

    aln_start: int
    aln_end: int
    ref_start: int
    ref_end: int
    length: int
    mean_identity: float


@dataclass
class ConservationResult:
    """Per-column conservation plus regions, pairwise and feature summary."""

    identity: list[float]
    gap_fraction: list[float]
    classes: list[str]
    regions: list[ConservedRegion]
    pairwise: pd.DataFrame
    feature_summary: pd.DataFrame
    summary: dict


@dataclass
class PipelineResult:
    """Everything one run produced, rendered by the UI (§7.1, §10)."""

    query: QueryInfo | None
    blast_results: dict[str, BlastResult]
    selected: list[SelectedHit]
    input_records: list[SequenceRecord]
    alignment: Alignment | None
    conservation: ConservationResult | None
    run_metadata: dict = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    #: Raw artefacts kept in memory for FR-32 downloads.
    raw_blast_xml: dict[str, str] = field(default_factory=dict)
    raw_msa_fasta: str = ""


# --- Phase 3 shared contract (§6.1, §7.2) --------------------------------------

#: Progress events emitted as ``(stage, detail)``. Stages are stable strings
#: used by the UI to update its step indicator; ``detail`` is human-readable.
ProgressCallback = Callable[[str, str], None]


@dataclass(frozen=True)
class PipelineInputs:
    """All inputs of one pipeline run, including offline overrides.

    Exactly one of ``query_accession`` and ``query_fasta`` must be supplied,
    unless ``query_info`` is given (the bundled example snapshot). Uploaded
    BLAST XML / MSA FASTA replace the corresponding remote stages; full-length
    sequences keyed by accession replace Entrez fetches in example mode.
    """

    query_accession: str | None = None
    accession_source: Literal["uniprot", "ncbi"] = "uniprot"
    query_fasta: str | None = None
    query_organism: str = ""
    species: tuple[str, ...] = config.SPECIES
    email: str = ""
    api_key: str | None = None
    blast_params: config.BlastParams = field(default_factory=config.BlastParams)
    selection_params: config.SelectionParams = field(
        default_factory=config.SelectionParams
    )
    conservation_params: config.ConservationParams = field(
        default_factory=config.ConservationParams
    )
    blast_xml_by_species: Mapping[str, str] | None = None
    msa_fasta: str | None = None
    full_sequences: Mapping[str, str] | None = None
    query_info: QueryInfo | None = None
    #: RIDs saved by the UI on a previous page run; skip re-submission (§12).
    existing_rids: Mapping[str, str] | None = None
    #: EBI Clustal Omega job ID saved on a previous page run (§12).
    msa_job_id: str | None = None
    #: Smoothing window chosen in advanced settings (FR-05, reproducibility).
    smoothing_window: int = config.SMOOTHING_WINDOW_DEFAULT


@dataclass
class PipelineClients:
    """Network operations injected into the pipeline (§6.2).

    The app binds real clients with their HTTP sessions, contact email and
    parameters; offline tests substitute fakes. The pipeline never imports
    Streamlit or constructs a session itself.
    """

    fetch_uniprot: Callable[[str], QueryInfo]
    fetch_ncbi_query_info: Callable[[str], QueryInfo]
    fetch_protein: Callable[[str], SequenceRecord]
    submit_blastp_staggered: Callable[[list[tuple[str, str]]], dict[str, str]]
    #: ``(rid, species) -> raw BLAST XML``.
    poll_blast: Callable[[str, str], str]
    run_clustalo: Callable[[str], str]
    #: ``(sequence, species) -> RID``. When supplied, the pipeline submits and
    #: saves each RID immediately, keeping the stagger between species (§9.1).
    submit_blastp: Callable[[str, str], str] | None = None
    #: ``(rid) -> PollOutcome-like (state, value, error)``. Enables one
    #: round-robin polling loop for all RIDs instead of one blocking call each.
    poll_blast_once: Callable[[str], Any] | None = None
    #: ``(fasta) -> EBI job ID``. Lets the app persist the job ID as soon as it
    #: is issued so a page refresh can resume polling (§12).
    submit_clustalo: Callable[[str], str] | None = None
    #: ``(job_id) -> aligned FASTA``. Resumes a saved EBI job.
    resume_clustalo: Callable[[str], str] | None = None
