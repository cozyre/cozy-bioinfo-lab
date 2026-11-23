"""User-facing strings and colours.

Everything the UI shows lives here (A1: English UI) so it can be localised
later, plus the colour-blind-safe palette used by ``viz`` (§11).

This module is pure data (no imports, no framework dependency), so domain
modules may import labels and messages from it. This is a narrow, deliberate
extension of the §6.2 import list: it keeps a single source of truth for every
displayed string without leaking Streamlit into the core.
"""

from __future__ import annotations

# --- App chrome ---------------------------------------------------------------

APP_TITLE: str = "TP53 Conservation Analyzer"
APP_SUBTITLE: str = (
    "BLASTp ortholog search, multiple sequence alignment and residue "
    "conservation for TP53 across five species."
)

# --- Colour-blind-safe palette (Okabe & Ito) ----------------------------------

COLOR_PALETTE: dict[str, str] = {
    "orange": "#E69F00",
    "sky_blue": "#56B4E9",
    "bluish_green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "reddish_purple": "#CC79A7",
    "black": "#000000",
}

#: Colour per conservation class, ordered from most to least conserved.
CLASS_COLORS: dict[str, str] = {
    "fully conserved": COLOR_PALETTE["bluish_green"],
    "high": COLOR_PALETTE["sky_blue"],
    "moderate": COLOR_PALETTE["yellow"],
    "low": COLOR_PALETTE["orange"],
    "gapped": COLOR_PALETTE["vermillion"],
}

#: Conservation class labels in display order (§9.5).
CONSERVATION_CLASSES: tuple[str, ...] = (
    "fully conserved",
    "high",
    "moderate",
    "low",
    "gapped",
)

#: Feature-summary table columns (FR-22, §9.9).
FEATURE_SUMMARY_COLUMNS: tuple[str, ...] = (
    "feature",
    "type",
    "start",
    "end",
    "columns",
    "mean_identity",
    "fully_conserved_pct",
)

#: Label of the final row that aggregates positions outside all features.
OUTSIDE_FEATURES_LABEL: str = "Outside all annotated features"

# --- Validation error messages (FR-03, FR-16) ----------------------------------

ERR_NOT_FASTA: str = "Input is not a FASTA file: no '>' header line found."
ERR_EMPTY_FASTA: str = "Input contains no FASTA records."
ERR_ONE_RECORD: str = (
    "FASTA must contain exactly one record, found {count}."
)
ERR_INVALID_CHARS: str = (
    "Sequence '{id}' contains characters that are not valid amino-acid codes: "
    "{chars}."
)
ERR_LENGTH: str = (
    "Sequence '{id}' has length {length}, outside the allowed range "
    "{min_len}-{max_len} residues."
)
ERR_EMPTY_SEQUENCE: str = "Sequence '{id}' is empty."

# --- Alignment validation messages ---------------------------------------------

ERR_ROW_LENGTH: str = (
    "Alignment row '{id}' has length {length}, expected {expected}."
)
ERR_ROW_COUNT: str = (
    "Alignment has {rows} rows, expected {expected} input sequences."
)
ERR_BAD_CHARS: str = (
    "Alignment row '{id}' contains invalid character(s): {chars}."
)
ERR_DEGAP_MISMATCH: str = (
    "De-gapped alignment row '{id}' does not reproduce the original input "
    "sequence (hard consistency check, §9.4)."
)
ERR_REF_ROW_MISSING: str = "Reference sequence '{ref_id}' not found in the alignment."
ERR_EMPTY_ALIGNMENT: str = "Alignment contains no FASTA records."
WARN_ID_REMAPPED: str = "Alignment IDs did not match the input headers; mapped rows by order."

# --- Alignment statistics table (FR-16) -----------------------------------------

ALIGNMENT_STAT_COLUMNS: tuple[str, ...] = (
    "id",
    "length",
    "ungapped_length",
    "gap_pct",
)

# --- Phase 3 core (Agent A) ---------------------------------------------------

# Progress details. The stage names themselves are owned by pipeline.py.
PROGRESS_VALIDATE_ACCESSION: str = "validating {accession} from {source}"
PROGRESS_VALIDATE_FASTA: str = "validating FASTA query {identifier}"
PROGRESS_VALIDATE_SNAPSHOT: str = "using the supplied query snapshot"
PROGRESS_RETRIEVE_DONE: str = "retrieved query {accession}"
PROGRESS_BLAST_UPLOAD: str = "uploaded XML: {species}"
PROGRESS_BLAST_QUERY: str = "query organism: {species} (no BLAST job)"
PROGRESS_BLAST_RESUME: str = "resuming: {species} ({rid})"
PROGRESS_BLAST_SUBMITTED: str = "submitted: {species} ({rid})"
PROGRESS_BLAST_POLLING: str = "polling: {species}"
PROGRESS_BLAST_READY: str = "ready: {species}"
PROGRESS_SELECT_SELECTED: str = "selected: {species} ({accession})"
PROGRESS_SELECT_MISSING: str = "no hit: {species}"
PROGRESS_FETCH_SEQUENCE: str = "fetched: {accession}"
PROGRESS_FETCH_SUPPLIED: str = "supplied sequence: {accession}"
PROGRESS_MSA_SERVICE: str = "aligning {count} sequences"
PROGRESS_MSA_UPLOAD: str = "using the uploaded MSA ({count} sequences)"
PROGRESS_CONSERVATION_DONE: str = "conservation computed for {count} columns"
PROGRESS_DONE: str = "complete: {count} sequences"

# Pipeline errors and warnings.
ERR_NO_QUERY_SOURCE: str = (
    "Supply exactly one query source: query_info, query_accession or query_fasta."
)
ERR_EMPTY_ACCESSION: str = "Query accession is empty."
ERR_UNKNOWN_ACCESSION_SOURCE: str = (
    "Unknown accession source '{source}'; expected 'uniprot' or 'ncbi'."
)
ERR_NO_ORTHOLOGS: str = (
    "No ortholog remained after BLAST and hit selection; at least one "
    "non-query sequence is required."
)
ERR_NO_SPECIES: str = "No target species were selected."
ERR_NO_QUERY_RECORD: str = "The query has no sequence record."
SELECT_REASON_QUERY: str = (
    "Query sequence used for its own organism; no BLAST job submitted."
)
WARN_BLAST_FAILED: str = "BLAST failed for {species}: {error}. Species skipped."
WARN_BLAST_TIMEOUT: str = "BLAST job for {species} timed out; species skipped."
WARN_NO_HIT: str = "No qualifying BLAST hit for {species}; species skipped."
WARN_SEQUENCE_FAILED: str = (
    "Could not fetch the full-length sequence for {accession}: {error}. "
    "Hit skipped."
)
WARN_TRIMMED_COLUMNS: str = "Removed {count} all-gap columns from the alignment."

# Plot titles, trace names and labels (V1-V5).
FIG_TITLE_NO_DATA: str = "No data available for this view"
FIG_TITLE_CONSERVATION: str = "Per-column percent identity"
FIG_TITLE_ALIGNMENT: str = "Multiple sequence alignment"
FIG_TITLE_BLAST: str = "Selected hits per species"
FIG_TITLE_FEATURES: str = "Mean identity per annotated feature"
FIG_TITLE_HEATMAP: str = "Pairwise percent identity"
FIG_TRACE_IDENTITY: str = "Percent identity"
FIG_TRACE_THRESHOLD: str = "Region threshold"
FIG_TRACE_REGION: str = "Conserved region {number}"
FIG_TRACE_FEATURE: str = "Feature band"
FIG_TRACE_FEATURE_LABELS: str = "Feature labels"
FIG_TRACE_MEAN_IDENTITY: str = "Mean identity"
FIG_TRACE_FULLY_CONSERVED: str = "Fully conserved share"
FIG_AXIS_REFERENCE: str = "Reference position"
FIG_AXIS_ALIGNMENT: str = "Alignment column"
FIG_AXIS_PERCENT_IDENTITY: str = "Percent identity (%)"
FIG_AXIS_SPECIES: str = "Species"
FIG_AXIS_BIT_SCORE: str = "Bit score"
FIG_AXIS_FEATURE: str = "Feature"
FIG_AXIS_FEATURES: str = "Features"
FIG_AXIS_SEQUENCE: str = "Sequence"
FIG_COLORBAR_IDENTITY: str = "Identity (%)"
FIG_NO_FEATURES_NOTICE: str = (
    "No UniProt features are available for this query (FASTA input)."
)
FIG_HOVER_IDENTITY: str = (
    "Reference position: %{customdata[0]}<br>"
    "Alignment column: %{customdata[1]}<br>"
    "Residue: %{customdata[2]}<br>"
    "Identity: %{customdata[3]:.1f}%<br>"
    "Gap fraction: %{customdata[4]:.1f}%<extra></extra>"
)
FIG_HOVER_REGION: str = (
    "Conserved region %{customdata[0]}<br>Reference %{customdata[1]}-%{customdata[2]}"
    "<br>Mean identity: %{customdata[3]:.1f}%<extra></extra>"
)

# Alignment viewer (V2): CSS class stems, row labels and block heading.
ALN_ROOT_CLASS: str = "tp53-alignment"
ALN_BLOCK_HEADER: str = "Columns {start}-{end} | reference {ref_start}-{ref_end}"
ALN_LABEL_RULER_TENS: str = "Position"
ALN_LABEL_CONSENSUS: str = "Consensus"
ALN_LABEL_IDENTITY: str = "Identity"
ALN_CONSENSUS_NONE: str = "."
ALN_LEGEND_TITLE: str = "Residue properties"
ALN_LEGEND_IDENTITY: str = "Conservation classes"
ALN_PROPERTY_CLASSES: dict[str, str] = {
    "hydrophobic": "aa-hydrophobic",
    "polar": "aa-polar",
    "acidic": "aa-acidic",
    "basic": "aa-basic",
    "special": "aa-special",
    "ambiguous": "aa-ambiguous",
    "gap": "aa-gap",
}
ALN_PROPERTY_RESIDUES: dict[str, str] = {
    "hydrophobic": "AVILMFWY",
    "polar": "STNQ",
    "acidic": "DE",
    "basic": "KRH",
    "special": "CGP",
    "ambiguous": "XBZJUO",
    "gap": "-",
}
ALN_PROPERTY_COLORS: dict[str, str] = {
    "hydrophobic": COLOR_PALETTE["blue"],
    "polar": COLOR_PALETTE["sky_blue"],
    "acidic": COLOR_PALETTE["vermillion"],
    "basic": COLOR_PALETTE["reddish_purple"],
    "special": COLOR_PALETTE["bluish_green"],
    "ambiguous": COLOR_PALETTE["orange"],
    "gap": COLOR_PALETTE["black"],
}
ALN_CLASS_SYMBOLS: dict[str, str] = {
    "fully conserved": "|",
    "high": ":",
    "moderate": ".",
    "low": ",",
    "gapped": "-",
}
ALN_CLASS_CSS: dict[str, str] = {
    label: f"class-{label.replace(' ', '-')}" for label in CLASS_COLORS
}

# The query row in V3 is not a BLAST result (§9.1); it has no bit score.
FIG_QUERY_NO_BLAST: str = "query (no BLAST)"

# --- Phase 4 hardening: progress, error paths and all UI labels (A1, §12) ------

PROGRESS_MSA_SUBMITTED: str = "submitted MSA job {job_id}"
PROGRESS_MSA_RESUME: str = "resuming MSA job {job_id}"
WARN_MSA_FAILED: str = (
    "Clustal Omega failed: {error}. Upload an MSA FASTA to continue."
)
ERR_XML_NO_HITS: str = "BLAST XML for {species} contains no hits."

TAB_OVERVIEW: str = "Overview"
TAB_BLAST: str = "BLAST"
TAB_ALIGNMENT: str = "Alignment"
TAB_CONSERVATION: str = "Conservation"
TAB_DOWNLOADS: str = "Downloads"
TABS: tuple[str, ...] = (
    TAB_OVERVIEW,
    TAB_BLAST,
    TAB_ALIGNMENT,
    TAB_CONSERVATION,
    TAB_DOWNLOADS,
)

QUERY_MODE_ACCESSION: str = "Accession"
QUERY_MODE_FASTA: str = "FASTA sequence"
QUERY_MODES: tuple[str, ...] = (QUERY_MODE_ACCESSION, QUERY_MODE_FASTA)

SOURCE_REMOTE: str = "Remote services"
SOURCE_UPLOAD: str = "Uploaded files"
SOURCE_MODES: tuple[str, ...] = (SOURCE_REMOTE, SOURCE_UPLOAD)

UI_X_AXIS_REFERENCE: str = "Reference position"
UI_X_AXIS_COLUMN: str = "Alignment column"
UI_X_AXIS_CHOICES: tuple[str, ...] = (UI_X_AXIS_REFERENCE, UI_X_AXIS_COLUMN)

ERR_EMAIL_REQUIRED: str = "Enter a contact email before running remote services."
ERR_EMAIL_FORMAT: str = "Enter a valid contact email address."
ERR_ACCESSION_REQUIRED: str = "Enter an accession or switch to the FASTA query mode."
ERR_ACCESSION_FORMAT: str = "Enter a valid UniProt or NCBI protein accession."
ERR_FASTA_REQUIRED: str = "Paste a FASTA sequence or upload a .fa file."
ERR_ORGANISM_REQUIRED: str = "Enter the query organism for a FASTA query."
ERR_SPECIES_REQUIRED: str = "Select at least one target species."
ERR_XML_UPLOAD_REQUIRED: str = "Upload one BLAST XML file for every selected species."
ERR_UPLOAD_TOO_LARGE: str = "File '{name}' exceeds the {limit} upload limit."
ERR_UPLOAD_NOT_TEXT: str = "File '{name}' is not UTF-8 text."
ERR_XML_PARSE: str = "Could not parse BLAST XML for {species}: {detail}"
FEATURE_SKIP_NOTICE: str = (
    "Feature overlays and the feature summary are skipped because this query "
    "has no UniProt feature annotations."
)
IDLE_NOTICE: str = (
    "Select an accession or FASTA query in the sidebar, then choose "
    "Run analysis. Load example works fully offline."
)
IDENTITY_NOTE: str = (
    "Identity counts gaps in the denominator: a column with a gap cannot "
    "reach 100%."
)
RUNNING_NOTICE: str = "Analysis running. Step progress is shown above."
NO_RESULT_NOTICE: str = "No result is available. Run the analysis or load the example."
ANALYSIS_WARNINGS: str = "Analysis finished with warnings."
ANALYSIS_COMPLETE: str = "Analysis complete."
ANALYSIS_FAILED: str = "Analysis failed: {error}"
FAILED_ACTION: str = (
    "Upload the missing BLAST XML and/or MSA FASTA in the sidebar, or load the "
    "bundled example snapshot."
)
STAGE_ACTION: str = "The '{stage}' stage failed: {error}. {hint}"
STAGE_LABELS: dict[str, str] = {
    "validate": "input validation",
    "retrieve": "query retrieval",
    "blast": "BLAST",
    "select": "hit selection",
    "msa": "multiple sequence alignment",
    "conservation": "conservation analysis",
}
UPLOAD_HINT_BLAST: str = "Upload the saved BLAST XML for each species to continue."
UPLOAD_HINT_MSA: str = "Upload the saved MSA FASTA to continue."
UPLOAD_HINT_GENERIC: str = "Load the bundled example snapshot to continue."
EXAMPLE_ACTION: str = "Load example works fully offline."
NO_HITS_NOTICE: str = (
    "No ortholog hits were selected. Try uploading BLAST XML or load the example."
)
NO_ALIGNMENT_NOTICE: str = "No alignment is available for this run."
NO_CONSERVATION_NOTICE: str = "No conservation results are available for this run."
NO_REGIONS_NOTICE: str = "No conserved regions meet the current threshold."
NO_ARTIFACTS_NOTICE: str = "No raw artefacts are available for this run."
PROGRESS_LINE: str = "{stage}: {detail}"
STATUS_STEP: str = "Step {step} of {total}: {detail}"
ALL_STAGES_COMPLETE: str = "All pipeline stages completed."
CLASSES_PREFIX: str = "Classes: "
STATUS_VALIDATING: str = "Step {step} of {total}: validating input"
STATUS_COMPLETE: str = "Analysis complete"
UPLOAD_LIMIT_MB: str = "{count} MB"
UPLOAD_LIMIT_BYTES: str = "{count} bytes"

SIDEBAR_QUERY: str = "Query"
LABEL_QUERY_INPUT: str = "Query input"
LABEL_ACCESSION: str = "Protein accession"
HELP_ACCESSION: str = "UniProt (P04637) or NCBI protein (NP_035770.2) accession."
LABEL_FASTA: str = "Query FASTA"
LABEL_FASTA_UPLOAD: str = "Or upload query FASTA"
LABEL_QUERY_ORGANISM: str = "Query organism"
VALUE_ORGANISM_FROM_ACCESSION: str = "Derived from accession"
SIDEBAR_CONTACT: str = "Contact"
LABEL_EMAIL: str = "Contact email"
HELP_EMAIL: str = "Required by NCBI and EBI; never stored in the repository."
SIDEBAR_TARGETS: str = "Targets and source"
LABEL_SPECIES: str = "Target species"
LABEL_SOURCE: str = "Data source"
LABEL_USE_XML: str = "Use uploaded BLAST XML (one file per species)"
LABEL_USE_MSA: str = "Use uploaded MSA FASTA"
LABEL_MSA_UPLOAD: str = "MSA FASTA"
ADVANCED_SETTINGS: str = "Advanced settings"
LABEL_BLAST_DATABASE: str = "BLAST database"
LABEL_EVALUE: str = "E-value cutoff"
LABEL_COVERAGE: str = "Minimum query coverage (%)"
LABEL_THRESHOLD: str = "Conservation threshold (%)"
LABEL_MIN_REGION: str = "Minimum region length"
LABEL_SMOOTHING_ADVANCED: str = "Smoothing window (Conservation tab)"
BUTTON_RUN: str = "Run analysis"
BUTTON_EXAMPLE: str = "Load example"
BUTTON_RESET: str = "Reset"
BUTTON_RESUME: str = "Resume analysis"

LABEL_RUN_METADATA: str = "Run metadata"
LABEL_CONSERVATION_SUMMARY: str = "Conservation summary"
METRIC_MEAN_IDENTITY: str = "Mean identity"
METRIC_FULLY_CONSERVED: str = "Fully conserved positions"
METRIC_FULLY_SHARE: str = "Fully conserved share"
METRIC_LONGEST_REGION: str = "Longest conserved region"
METRIC_REGIONS: str = "Conserved regions"
METRIC_COLUMNS: str = "Alignment columns"
UNIT_COLUMNS: str = "{count} columns"
LABEL_SELECTED_HIT: str = "{species} — {accession}"
EXPANDER_RAW_ALIGNMENT: str = "Raw hit alignment — {species} / {accession}"
LABEL_HSP: str = "HSP {index}: {query_start}-{query_end}"
LABEL_SMOOTHING_SLIDER: str = "Smoothing window"
LABEL_X_AXIS: str = "X axis"
LABEL_REGIONS: str = "Conserved regions"
LABEL_FEATURE_COMPARISON: str = "Feature comparison"
LABEL_PAIRWISE: str = "Pairwise identity"
DOWNLOAD_BLAST_XML: str = "BLAST XML — {species}"
DOWNLOAD_MSA: str = "MSA FASTA"
DOWNLOAD_CSV: str = "Identity CSV"
FILE_BLAST_XML: str = "blast_{slug}.xml"
FILE_MSA: str = "tp53_alignment.fasta"
FILE_IDENTITY: str = "tp53_identity.csv"
MIME_XML: str = "application/xml"
MIME_TEXT: str = "text/plain"
MIME_CSV: str = "text/csv"
