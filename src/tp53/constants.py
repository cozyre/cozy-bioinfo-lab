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
