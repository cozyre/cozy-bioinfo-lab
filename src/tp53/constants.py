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
