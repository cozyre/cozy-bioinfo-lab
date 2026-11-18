"""Residue conservation: identity, classes, regions and feature summary.

Implements §9.5–§9.9. Gaps stay in the identity denominator, so a column
containing a gap can never reach 100% identity. Feature coordinates are looked
up through :attr:`Alignment.col_to_ref` (reference numbering), never by
alignment column.
"""

from __future__ import annotations

from collections import Counter

import pandas as pd

from tp53 import config
from tp53.constants import (
    CONSERVATION_CLASSES,
    FEATURE_SUMMARY_COLUMNS,
    OUTSIDE_FEATURES_LABEL,
)
from tp53.models import (
    Alignment,
    ConservationResult,
    ConservedRegion,
    DomainFeature,
)

#: The five class labels, from the single source in :mod:`tp53.constants`.
FULLY_CONSERVED, HIGH, MODERATE, LOW, GAPPED = CONSERVATION_CLASSES

#: Keys of :attr:`ConservationResult.summary` (consumed by the UI in phase 3).
SUMMARY_MEAN_IDENTITY = "mean_identity"
SUMMARY_FULLY_CONSERVED_COUNT = "fully_conserved_count"
SUMMARY_FULLY_CONSERVED_PCT = "fully_conserved_pct"
SUMMARY_REGION_COUNT = "region_count"
SUMMARY_LONGEST_REGION = "longest_region"
SUMMARY_ALIGNMENT_LENGTH = "alignment_length"


def column_identity(aln: Alignment) -> tuple[list[float], list[float]]:
    """Percent identity and gap fraction of every alignment column (§9.5).

    Ambiguous codes (XBZJUO) are ignored in the residue count, but all
    sequences remain in the denominator. An empty alignment yields empty
    lists.

    Args:
        aln: parsed and validated alignment.

    Returns:
        ``(identity, gap_fraction)`` lists, both in percent and one entry per
        alignment column.
    """
    identity: list[float] = []
    gap_fraction: list[float] = []
    total = len(aln.rows)
    if not total:
        return identity, gap_fraction

    standard = set(config.STANDARD_AMINO_ACIDS)
    for column_index in range(aln.length):
        column = [row[column_index] for row in aln.rows]
        counts = Counter(char for char in column if char in standard)
        top_count = max(counts.values(), default=0)
        gaps = sum(char == config.GAP_CHAR for char in column)
        identity.append(100.0 * top_count / total)
        gap_fraction.append(100.0 * gaps / total)
    return identity, gap_fraction


def classify_columns(
    identity: list[float],
    gap_fraction: list[float],
    params: config.ConservationParams | None = None,
) -> list[str]:
    """Classify each column as fully conserved/high/moderate/low/gapped (§9.5).

    Identity takes priority: gap columns with a clear majority residue are
    still reported as high/moderate, while a gapped column with no conserved
    majority is labelled ``gapped``.

    Args:
        identity: percent identity per column.
        gap_fraction: percent gaps per column.
        params: thresholds; defaults to :class:`ConservationParams`.

    Returns:
        Class label per column, in the same order.
    """
    params = params or config.ConservationParams()
    classes: list[str] = []
    for percent, gaps in zip(identity, gap_fraction, strict=True):
        if percent >= 100.0 and gaps == 0:
            classes.append(FULLY_CONSERVED)
        elif percent >= params.high_threshold:
            classes.append(HIGH)
        elif percent >= params.moderate_threshold:
            classes.append(MODERATE)
        elif gaps > 0:
            classes.append(GAPPED)
        else:
            classes.append(LOW)
    return classes


def _make_region(
    aln: Alignment,
    identity: list[float],
    start: int,
    end: int,
    params: config.ConservationParams,
) -> ConservedRegion | None:
    """Build one region from an inclusive run of alignment column indices."""
    if end - start + 1 < params.min_region_length:
        return None
    reference_positions = [
        position
        for position in aln.col_to_ref[start : end + 1]
        if position is not None
    ]
    return ConservedRegion(
        aln_start=start + 1,
        aln_end=end + 1,
        ref_start=reference_positions[0] if reference_positions else 0,
        ref_end=reference_positions[-1] if reference_positions else 0,
        length=end - start + 1,
        mean_identity=sum(identity[start : end + 1]) / (end - start + 1),
    )


def find_regions(
    aln: Alignment,
    identity: list[float],
    gap_fraction: list[float],
    params: config.ConservationParams | None = None,
) -> list[ConservedRegion]:
    """Detect runs of conserved columns (§9.6, FR-20).

    A column qualifies when identity is at or above the configured threshold
    and the gap fraction is zero. Runs shorter than the minimum length are
    discarded. Coordinates are 1-based inclusive in both alignment and
    reference numbering; reference positions come from the column map, using
    the first and last mapped column of the run.

    Args:
        aln: alignment providing the reference column map.
        identity: percent identity per column.
        gap_fraction: percent gaps per column.
        params: threshold and minimum length; defaults apply otherwise.

    Returns:
        Conserved regions in alignment order.
    """
    params = params or config.ConservationParams()
    regions: list[ConservedRegion] = []
    run_start: int | None = None

    for index, (percent, gaps) in enumerate(
        zip(identity, gap_fraction, strict=True)
    ):
        qualifies = percent >= params.region_identity_threshold and gaps == 0
        if qualifies and run_start is None:
            run_start = index
        elif not qualifies and run_start is not None:
            region = _make_region(aln, identity, run_start, index - 1, params)
            if region is not None:
                regions.append(region)
            run_start = None

    if run_start is not None:
        region = _make_region(aln, identity, run_start, len(identity) - 1, params)
        if region is not None:
            regions.append(region)
    return regions


def pairwise_identity(aln: Alignment) -> pd.DataFrame:
    """Pairwise percent identity matrix over shared-residue columns (§9.8).

    Only columns where both rows have a non-gap character contribute;
    ambiguous residues count as residues. Rows and columns are labelled with
    alignment IDs, so the table can be rendered directly as a heatmap.

    Args:
        aln: parsed alignment.

    Returns:
        Symmetric DataFrame with IDs as index and columns.
    """
    matrix = pd.DataFrame(0.0, index=aln.ids, columns=aln.ids, dtype=float)
    for first in range(len(aln.rows)):
        for second in range(first, len(aln.rows)):
            pairs = [
                (a, b)
                for a, b in zip(
                    aln.rows[first], aln.rows[second], strict=True
                )
                if a != config.GAP_CHAR and b != config.GAP_CHAR
            ]
            if pairs:
                percent = 100.0 * sum(a == b for a, b in pairs) / len(pairs)
            else:
                percent = 0.0
            matrix.iloc[first, second] = percent
            matrix.iloc[second, first] = percent
    return matrix


def _add_feature_row(
    rows: dict[str, list],
    label: str,
    feature_type: str,
    start: int | str,
    end: int | str,
    columns: list[int],
    identity: list[float],
) -> None:
    """Append one feature (or the outside-features row) to the summary."""
    values = [identity[index] for index in columns]
    rows["feature"].append(label)
    rows["type"].append(feature_type)
    rows["start"].append(start)
    rows["end"].append(end)
    rows["columns"].append(len(columns))
    rows["mean_identity"].append(sum(values) / len(values) if values else 0.0)
    rows["fully_conserved_pct"].append(
        100.0 * sum(value >= 100.0 for value in values) / len(values)
        if values
        else 0.0
    )


def summarise_features(
    aln: Alignment,
    identity: list[float],
    features: list[DomainFeature],
) -> pd.DataFrame:
    """Mean identity per UniProt feature and outside all features (§9.9, FR-22).

    Feature start/end are 1-based reference positions resolved through the
    column map. Overlapping features may share columns; the final row uses
    only columns not covered by any feature. Columns where the reference has
    a gap have no reference position and are excluded from both.

    Args:
        aln: alignment with a reference column map.
        identity: percent identity per alignment column.
        features: annotated query features in reference coordinates.

    Returns:
        DataFrame with :data:`FEATURE_SUMMARY_COLUMNS`, plus one final
        outside-features row.
    """
    rows: dict[str, list] = {name: [] for name in FEATURE_SUMMARY_COLUMNS}
    covered: set[int] = set()
    for feature in features:
        indices = [
            index
            for index, position in enumerate(aln.col_to_ref)
            if position is not None and feature.start <= position <= feature.end
        ]
        covered.update(indices)
        _add_feature_row(
            rows,
            feature.name,
            feature.type,
            feature.start,
            feature.end,
            indices,
            identity,
        )
    outside = [
        index
        for index, position in enumerate(aln.col_to_ref)
        if position is not None and index not in covered
    ]
    _add_feature_row(rows, OUTSIDE_FEATURES_LABEL, "", "", "", outside, identity)
    return pd.DataFrame(rows, columns=list(FEATURE_SUMMARY_COLUMNS))


def summarise_conservation(
    aln: Alignment,
    identity: list[float],
    classes: list[str],
    regions: list[ConservedRegion],
) -> dict:
    """Aggregate the summary metrics shown in the Overview tab (FR-23).

    Args:
        aln: alignment (for its length).
        identity: percent identity per column.
        classes: class label per column.
        regions: detected conserved regions.

    Returns:
        Dict keyed by the ``SUMMARY_*`` constants defined in this module.
    """
    fully_conserved = sum(1 for label in classes if label == FULLY_CONSERVED)
    total = len(identity)
    return {
        SUMMARY_MEAN_IDENTITY: sum(identity) / total if total else 0.0,
        SUMMARY_FULLY_CONSERVED_COUNT: fully_conserved,
        SUMMARY_FULLY_CONSERVED_PCT: 100.0 * fully_conserved / total if total else 0.0,
        SUMMARY_REGION_COUNT: len(regions),
        SUMMARY_LONGEST_REGION: max((region.length for region in regions), default=0),
        SUMMARY_ALIGNMENT_LENGTH: aln.length,
    }


def compute_conservation(
    aln: Alignment,
    features: list[DomainFeature] | None = None,
    params: config.ConservationParams | None = None,
) -> ConservationResult:
    """Run all conservation steps and assemble a :class:`ConservationResult`.

    Args:
        aln: parsed and validated alignment.
        features: query feature annotations; an empty list means no overlays
            (e.g. a bare FASTA query).
        params: conservation thresholds; defaults apply otherwise.

    Returns:
        Full per-column result with regions, pairwise matrix, feature table
        and summary metrics.
    """
    identity, gap_fraction = column_identity(aln)
    classes = classify_columns(identity, gap_fraction, params)
    regions = find_regions(aln, identity, gap_fraction, params)
    pairwise = pairwise_identity(aln)
    feature_summary = summarise_features(aln, identity, features or [])
    summary = summarise_conservation(aln, identity, classes, regions)
    return ConservationResult(
        identity=identity,
        gap_fraction=gap_fraction,
        classes=classes,
        regions=regions,
        pairwise=pairwise,
        feature_summary=feature_summary,
        summary=summary,
    )
