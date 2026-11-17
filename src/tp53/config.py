"""Central configuration for the TP53 Conservation Analyzer.

All defaults, thresholds, URLs and species lists live here (design doc §5,
"Maintainability"). No magic numbers anywhere else in the code base.
"""

from __future__ import annotations

from dataclasses import dataclass

# --- Species -----------------------------------------------------------------

#: The five case species of the study, in research-question order (§1).
SPECIES: tuple[str, ...] = (
    "Homo sapiens",
    "Pan troglodytes",
    "Mus musculus",
    "Rattus norvegicus",
    "Bos taurus",
)

#: Default query accession (human TP53, canonical isoform, 393 aa, A2).
DEFAULT_QUERY_ACCESSION: str = "P04637"

# --- Validation (FR-03) ------------------------------------------------------

FASTA_MIN_LEN: int = 20
FASTA_MAX_LEN: int = 5000

#: The 20 standard amino acids plus the IUPAC ambiguity codes accepted
#: anywhere in the app (§9.4). Gap handling is separate (GAP_CHAR).
PROTEIN_ALPHABET: str = "ACDEFGHIKLMNPQRSTVWYXBZJUO"

#: Alignment gap character used by Clustal Omega and this app.
GAP_CHAR: str = "-"

#: Line width when writing FASTA text.
FASTA_LINE_WIDTH: int = 60

#: Separator used in MSA headers: ``Species_name|Accession`` (FR-13).
HEADER_SEPARATOR: str = "|"

# --- URLs (§8). Never hard-code these in client logic. ------------------------

UNIPROT_BASE_URL: str = "https://rest.uniprot.org/uniprotkb"
NCBI_EFETCH_URL: str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
NCBI_BLAST_URL: str = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
EBI_CLUSTALO_URL: str = "https://www.ebi.ac.uk/Tools/services/rest/clustalo"

# --- HTTP etiquette (§5 Reliability, Service etiquette) ----------------------

#: Per-request timeout in seconds for every HTTP call.
REQUEST_TIMEOUT_S: int = 30

#: Retry policy: attempts and exponential backoff base in seconds.
HTTP_MAX_RETRIES: int = 3
HTTP_BACKOFF_FACTOR_S: float = 2.0

#: NCBI BLAST polling: never poll more often than this (check the current
#: NCBI usage guidelines again when implementing the clients in phase 2).
BLAST_POLL_INTERVAL_S: int = 60

#: Overall ceiling on a single remote BLAST job.
BLAST_MAX_WAIT_S: int = 15 * 60

#: Delay between staggered per-species BLAST submissions (§9.1).
BLAST_STAGGER_DELAY_S: int = 10

#: EBI Clustal Omega polling interval.
EBI_POLL_INTERVAL_S: int = 3

#: Overall ceiling on a single Clustal Omega job.
EBI_MAX_WAIT_S: int = 10 * 60

#: BLAST hit list size (§8, small hit list).
BLAST_HITLIST_SIZE: int = 10


@dataclass(frozen=True)
class BlastParams:
    """Parameters of one remote BLASTp job (§9.1)."""

    database: str = "refseq_protein"
    expect: float = 1e-10
    hitlist_size: int = BLAST_HITLIST_SIZE
    program: str = "blastp"


#: BLAST databases offered in the advanced settings (FR-05).
BLAST_DATABASES: tuple[str, ...] = ("refseq_protein", "nr", "swissprot")


@dataclass(frozen=True)
class SelectionParams:
    """Best-hit selection rules (§9.2)."""

    evalue_cutoff: float = 1e-10
    #: Minimum aligned-query-length / query-length, as a fraction.
    min_query_coverage: float = 0.70


@dataclass(frozen=True)
class ConservationParams:
    """Conservation thresholds (§9.5, §9.6)."""

    #: A conserved region needs every column at or above this identity.
    region_identity_threshold: float = 100.0
    #: Minimum number of alignment columns in a conserved region.
    min_region_length: int = 5
    #: Class boundaries in percent identity (§9.5). With five sequences the
    #: possible values are 20/40/60/80/100 when no gaps are present.
    high_threshold: float = 80.0
    moderate_threshold: float = 60.0


#: Advanced-settings default for the smoothing window (FR-05).
SMOOTHING_WINDOW_DEFAULT: int = 5

#: Upload size limit in bytes (§5 Security) — 5 MB.
MAX_UPLOAD_BYTES: int = 5 * 1024 * 1024
