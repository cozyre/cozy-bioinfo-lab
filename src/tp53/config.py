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

#: The 20 standard amino acids (counted in column identity, §9.5).
STANDARD_AMINO_ACIDS: str = "ACDEFGHIKLMNPQRSTVWY"

#: The 20 standard amino acids plus the IUPAC ambiguity codes accepted
#: anywhere in the app (§9.4). Gap handling is separate (GAP_CHAR).
PROTEIN_ALPHABET: str = STANDARD_AMINO_ACIDS + "XBZJUO"

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

#: Resource paths appended to :data:`EBI_CLUSTALO_URL` (§8). Clients never
#: assemble endpoints from literals of their own.
EBI_CLUSTALO_RUN_PATH: str = "/run"
EBI_CLUSTALO_STATUS_PATH: str = "/status/{job_id}"
EBI_CLUSTALO_RESULT_PATH: str = "/result/{job_id}/{result_type}"

# --- HTTP etiquette (§5 Reliability, Service etiquette) ----------------------

#: Per-request timeout in seconds for every HTTP call.
REQUEST_TIMEOUT_S: int = 30

#: Retry policy: attempts and exponential backoff base in seconds.
HTTP_MAX_RETRIES: int = 3
HTTP_BACKOFF_FACTOR_S: float = 2.0

#: Only these HTTP statuses are worth retrying; other 4xx errors are I/O
#: errors in the request itself and are surfaced immediately.
HTTP_RETRYABLE_STATUSES: tuple[int, ...] = (429, 500, 502, 503, 504)

#: Tool name sent to NCBI/EBI so they can identify traffic from this app
#: (NCBI BLAST URL API usage guidelines: email + tool are required).
NCBI_TOOL_NAME: str = "tp53_conservation_analyzer"

#: Only features with exact boundaries can be mapped to reference positions.
UNIPROT_EXACT_MODIFIER: str = "EXACT"

#: UniProt feature types shown in the app and included in the feature summary
#: (FR-06, §9.9). "DNA binding" is the spelling in the 2026 UniProt JSON;
#: "DNA-binding" is accepted for forward compatibility.
UNIPROT_FEATURE_TYPES: tuple[str, ...] = (
    "Domain",
    "Region",
    "DNA binding",
    "DNA-binding",
    "Motif",
    "Zinc finger",
)

#: NCBI BLAST polling: never poll one RID more often than this. Current NCBI
#: guidelines (verified 2026-11-18): at most one request every 10 s per client
#: and at most one poll per RID per minute; send email and tool on every call.
BLAST_POLL_INTERVAL_S: int = 60

#: Entrez efetch formats: GenBank text carries the feature table needed for
#: FR-06, while FASTA is used for full-length hit sequences (FR-12).
NCBI_EFETCH_GB_RETTYPE: str = "gb"
NCBI_EFETCH_TEXT_RETMODE: str = "text"

#: Entrez GenBank feature types mapped onto the canonical names used by
#: UniProt (§9.9). This is a vocabulary mapping, not a domain annotation.
NCBI_FEATURE_TYPES: dict[str, str] = {
    "Domain": "Domain",
    "Region": "Region",
    "Motif": "Motif",
    "ZnF_C2H2": "Zinc finger",
    "ZnF_C4": "Zinc finger",
    "DNA binding": "DNA binding",
}

#: NCBI BLAST result format used by the clients (the parser expects XML).
BLAST_RESULT_FORMAT: str = "XML"

#: Overall ceiling on a single remote BLAST job.
BLAST_MAX_WAIT_S: int = 15 * 60

#: Delay between staggered per-species BLAST submissions (§9.1).
BLAST_STAGGER_DELAY_S: int = 10

#: EBI Clustal Omega REST parameters (verified against current job-dispatcher
#: documentation: email is required, ``out`` returns the aligned FASTA).
EBI_CLUSTALO_SEQUENCE_TYPE: str = "protein"
EBI_CLUSTALO_RESULT_TYPE: str = "out"

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

# --- Caching and bundled example snapshot (FR-31, FR-33, A8) -------------------

#: Time-to-live for cached client calls and parsed results, in seconds.
CACHE_TTL_S: int = 6 * 60 * 60

#: The example snapshot is bundled under ``app/`` so it ships with the app.
EXAMPLE_DATA_DIRNAME: str = "example_data"
EXAMPLE_BLAST_XML_FILENAME: str = "blast_tp53_five_species.xml"
EXAMPLE_MSA_FASTA_FILENAME: str = "msa_tp53_five_species.fasta"
EXAMPLE_UNIPROT_JSON_FILENAME: str = "uniprot_P04637.json"
EXAMPLE_SEQUENCE_FILENAMES: tuple[str, ...] = (
    "NP_001120705.1.fa",
    "NP_001394193.1.fa",
    "NP_001394197.1.fa",
    "NP_001416925.1.fa",
    "NP_035770.2.fa",
    "NP_112251.2.fa",
    "NP_776626.1.fa",
    "XP_001172077.2.fa",
    "XP_006533220.1.fa",
    "XP_016786959.2.fa",
)
