"""Framework-free visualisations: Plotly figures (V1, V3-V5) and alignment HTML (V2).

Every function takes domain result objects from :mod:`tp53.models` and returns
a :class:`plotly.graph_objects.Figure` or a standalone HTML string. Nothing in
this module imports Streamlit or performs I/O, so the figures can be rendered
and tested offline.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import fields
from html import escape
from typing import Any

import plotly.graph_objects as go

from tp53 import config, constants
from tp53.conservation import classify_columns
from tp53.models import Alignment, PipelineResult, SelectedHit

#: Structural constants (not tunable defaults; those live in config.py).
_Y_MIN = 0.0
_Y_MAX = 100.0
_BAND_Y_MIN = 0.0
_BAND_Y_MAX = 1.0
_BAND_DOMAIN = (0.0, 0.12)
_PLOT_DOMAIN = (0.18, 1.0)
_TENS_STEP = 10
_FULL_WINDOW = 1
_HEATMAP_COLORSCALE = (
    (0.0, constants.COLOR_PALETTE["vermillion"]),
    (0.5, constants.COLOR_PALETTE["yellow"]),
    (1.0, constants.COLOR_PALETTE["bluish_green"]),
)
_COLUMN_AXIS_NAMES = frozenset({"alignment", "column", "aln"})


def _empty_figure(title: str, message: str | None = None) -> go.Figure:
    """Return a valid figure with a message instead of raising on empty data."""
    fig = go.Figure()
    fig.add_annotation(
        text=message or constants.FIG_TITLE_NO_DATA,
        showarrow=False,
        xref="paper",
        yref="paper",
        x=0.5,
        y=0.5,
    )
    fig.update_layout(
        title=title,
        xaxis={"visible": False},
        yaxis={"visible": False},
    )
    return fig


def _conservation_params(result: PipelineResult) -> config.ConservationParams:
    """Read applied thresholds from run metadata, falling back to config."""
    raw: Any = result.run_metadata.get("conservation_params")
    if isinstance(raw, config.ConservationParams):
        return raw
    if isinstance(raw, dict):
        allowed = {field.name for field in fields(config.ConservationParams)}
        kwargs = {key: value for key, value in raw.items() if key in allowed}
        try:
            return config.ConservationParams(**kwargs)
        except TypeError:
            pass
    return config.ConservationParams()


def _smooth(values: list[float], window: int) -> list[float]:
    """ centred moving average; ``window <= 1`` leaves the values untouched."""
    if window <= _FULL_WINDOW or not values:
        return list(values)
    half = window // 2
    smoothed: list[float] = []
    for index, value in enumerate(values):
        start = max(0, index - half)
        end = min(len(values), index + half + _FULL_WINDOW)
        chunk = values[start:end]
        smoothed.append(sum(chunk) / len(chunk) if chunk else value)
    return smoothed


def _reference_positions(aln: Alignment) -> list[int | None]:
    """1-based reference position of each column (§9.7), with a safe fallback."""
    if aln.col_to_ref and len(aln.col_to_ref) == aln.length:
        return list(aln.col_to_ref)
    if not aln.rows or not 0 <= aln.ref_index < len(aln.rows):
        return [None] * aln.length
    positions: list[int | None] = []
    reference = 0
    for char in aln.rows[aln.ref_index]:
        if char != config.GAP_CHAR:
            reference += 1
            positions.append(reference)
        else:
            positions.append(None)
    return positions


def _span_x(
    positions: list[int | None], start_col: int, end_col: int, *, reference: bool
) -> tuple[float, float] | None:
    """X-range of an inclusive column span in the selected coordinate system."""
    if not reference:
        return float(start_col), float(end_col)
    mapped = [
        positions[index]
        for index in range(start_col - 1, min(end_col, len(positions)))
        if positions[index] is not None
    ]
    if not mapped:
        return None
    return float(mapped[0]), float(mapped[-1])


def _column_of_reference(positions: list[int | None], position: int) -> int | None:
    """First 1-based alignment column holding a reference position."""
    for index, mapped in enumerate(positions):
        if mapped == position:
            return index + 1
    return None


def conservation_profile(
    result: PipelineResult,
    *,
    smoothing_window: int = config.SMOOTHING_WINDOW_DEFAULT,
    x_axis: str = "reference",
) -> go.Figure:
    """Build the V1 conservation profile with regions, features and an x toggle.

    The applied region threshold is taken from
    ``result.run_metadata["conservation_params"]`` when available, otherwise
    from :class:`config.ConservationParams`. Smoothing is a centred moving
    average over ``smoothing_window`` columns; the hover data always reports
    the unsmoothed identity.

    Args:
        result: completed pipeline result.
        smoothing_window: width of the moving average (1 disables smoothing).
        x_axis: initial axis, ``"reference"`` (default) or ``"alignment"``.

    Returns:
        Plotly figure with a 0-100 identity axis and an optional feature track.
    """
    aln = result.alignment
    conservation = result.conservation
    if aln is None or conservation is None or not conservation.identity or not aln.rows:
        return _empty_figure(constants.FIG_TITLE_CONSERVATION)

    params = _conservation_params(result)
    positions = _reference_positions(aln)
    use_reference = x_axis not in _COLUMN_AXIS_NAMES
    ref_x: list[float | None] = [float(p) if p is not None else None for p in positions]
    col_x: list[float] = [float(index) for index in range(1, aln.length + 1)]
    current_x = ref_x if use_reference else col_x
    raw_identity = conservation.identity
    identity = _smooth(raw_identity, smoothing_window)
    residue_row = aln.rows[aln.ref_index] if 0 <= aln.ref_index < len(aln.rows) else ""

    fig = go.Figure()
    x_by_reference: list[list[float | None]] = []
    x_by_column: list[list[float]] = []

    # Conserved-region shading, kept behind the identity line.
    for number, region in enumerate(conservation.regions, start=1):
        reference_span = _span_x(
            positions, region.aln_start, region.aln_end, reference=True
        )
        column_span = _span_x(
            positions, region.aln_start, region.aln_end, reference=False
        )
        if reference_span is None or column_span is None:
            continue
        start, end = (reference_span if use_reference else column_span)
        hover = [
            number,
            region.ref_start,
            region.ref_end,
            region.mean_identity,
        ]
        fig.add_trace(
            go.Scatter(
                x=[start, end, end, start],
                y=[_Y_MIN, _Y_MAX, _Y_MAX, _Y_MIN],
                fill="toself",
                fillcolor=constants.CLASS_COLORS["fully conserved"],
                line={"width": 0},
                opacity=0.15,
                mode="lines",
                name=constants.FIG_TRACE_REGION.format(number=number),
                legendgroup="regions",
                showlegend=number == 1,
                hovertemplate=constants.FIG_HOVER_REGION,
                customdata=[hover] * 4,
            )
        )
        x_by_reference.append(
            [reference_span[0], reference_span[1], reference_span[1], reference_span[0]]
        )
        x_by_column.append(
            [column_span[0], column_span[1], column_span[1], column_span[0]]
        )

    # Threshold line for conserved-region detection (§9.6).
    mapped = [position for position in positions if position is not None]
    ref_edges = [mapped[0], mapped[-1]] if mapped else [_Y_MIN, float(aln.length)]
    col_edges = [1.0, float(aln.length)]
    edges = ref_edges if use_reference else col_edges
    fig.add_trace(
        go.Scatter(
            x=edges,
            y=[params.region_identity_threshold] * 2,
            mode="lines",
            line={"dash": "dash", "color": constants.COLOR_PALETTE["orange"]},
            name=constants.FIG_TRACE_THRESHOLD,
            hovertemplate=(
                f"{constants.FIG_TRACE_THRESHOLD}: "
                f"{params.region_identity_threshold:g}%<extra></extra>"
            ),
        )
    )
    x_by_reference.append([float(edge) for edge in ref_edges])
    x_by_column.append([float(edge) for edge in col_edges])

    # Identity line with full, unsmoothed hover information.
    customdata = [
        [
            positions[index],
            index + 1,
            residue_row[index] if index < len(residue_row) else "",
            raw_identity[index],
            conservation.gap_fraction[index]
            if index < len(conservation.gap_fraction)
            else 0.0,
        ]
        for index in range(len(identity))
    ]
    fig.add_trace(
        go.Scatter(
            x=list(current_x),
            y=identity,
            mode="lines+markers",
            line={"color": constants.COLOR_PALETTE["blue"], "width": 2},
            marker={"size": 4},
            name=constants.FIG_TRACE_IDENTITY,
            customdata=customdata,
            hovertemplate=constants.FIG_HOVER_IDENTITY,
        )
    )
    x_by_reference.append(ref_x)
    x_by_column.append(col_x)

    # Feature bands live in their own narrow track below the identity plot.
    features = result.query.features if result.query is not None else []
    feature_tracks = 0
    for feature in features:
        reference_span = _span_x(
            positions, feature.start, feature.end, reference=True
        )
        if reference_span is None:
            continue
        column_start = _column_of_reference(positions, feature.start)
        column_end = _column_of_reference(positions, feature.end)
        if column_start is None or column_end is None:
            continue
        column_span = (float(column_start), float(column_end))
        start, end = reference_span if use_reference else column_span
        feature_tracks += 1
        fig.add_trace(
            go.Scatter(
                x=[start, end, end, start],
                y=[_BAND_Y_MIN, _BAND_Y_MIN, _BAND_Y_MAX, _BAND_Y_MAX],
                fill="toself",
                fillcolor=constants.COLOR_PALETTE["sky_blue"],
                line={"width": 0},
                opacity=0.55,
                mode="lines",
                yaxis="y2",
                name=constants.FIG_TRACE_FEATURE,
                legendgroup="features",
                showlegend=feature_tracks == 1,
                hovertemplate=(
                    f"{escape(str(feature.name))}<br>"
                    f"{escape(str(feature.type))} "
                    f"{feature.start}-{feature.end}<extra></extra>"
                ),
            )
        )
        x_by_reference.append(
            [reference_span[0], reference_span[1], reference_span[1], reference_span[0]]
        )
        x_by_column.append(
            [column_span[0], column_span[1], column_span[1], column_span[0]]
        )

    if feature_tracks:
        label_centres_ref: list[float] = []
        label_centres_col: list[float] = []
        label_text: list[str] = []
        for feature in features:
            reference_span = _span_x(
                positions, feature.start, feature.end, reference=True
            )
            if reference_span is None:
                continue
            column_start = _column_of_reference(positions, feature.start)
            column_end = _column_of_reference(positions, feature.end)
            if column_start is None or column_end is None:
                continue
            label_centres_ref.append(
                (reference_span[0] + reference_span[1]) / 2
            )
            label_centres_col.append((column_start + column_end) / 2)
            label_text.append(str(feature.name))
        fig.add_trace(
            go.Scatter(
                x=label_centres_ref if use_reference else label_centres_col,
                y=[(_BAND_Y_MIN + _BAND_Y_MAX) / 2] * len(label_text),
                text=label_text,
                textposition="middle center",
                mode="text",
                yaxis="y2",
                name=constants.FIG_TRACE_FEATURE_LABELS,
                showlegend=False,
                hoverinfo="skip",
            )
        )
        x_by_reference.append(label_centres_ref)
        x_by_column.append(label_centres_col)

    layout: dict[str, Any] = {
        "title": constants.FIG_TITLE_CONSERVATION,
        "xaxis": {
            "title": constants.FIG_AXIS_REFERENCE
            if use_reference
            else constants.FIG_AXIS_ALIGNMENT
        },
        "yaxis": {
            "title": constants.FIG_AXIS_PERCENT_IDENTITY,
            "range": [_Y_MIN, _Y_MAX],
            "domain": list(_PLOT_DOMAIN) if feature_tracks else [0.0, 1.0],
        },
        "margin": {"b": 70, "t": 80},
        "hovermode": "closest",
        "updatemenus": [
            {
                "type": "buttons",
                "direction": "right",
                "x": 0.0,
                "y": 1.12,
                "xanchor": "left",
                "yanchor": "bottom",
                "buttons": [
                    {
                        "label": constants.FIG_AXIS_REFERENCE,
                        "method": "update",
                        "args": [
                            {"x": x_by_reference},
                            {"xaxis.title.text": constants.FIG_AXIS_REFERENCE},
                        ],
                    },
                    {
                        "label": constants.FIG_AXIS_ALIGNMENT,
                        "method": "update",
                        "args": [
                            {"x": x_by_column},
                            {"xaxis.title.text": constants.FIG_AXIS_ALIGNMENT},
                        ],
                    },
                ],
            }
        ],
    }
    if feature_tracks:
        layout["yaxis2"] = {
            "title": constants.FIG_AXIS_FEATURES,
            "range": [_BAND_Y_MIN, _BAND_Y_MAX],
            "domain": list(_BAND_DOMAIN),
            "anchor": "x",
            "showticklabels": False,
        }
    else:
        fig.add_annotation(
            text=constants.FIG_NO_FEATURES_NOTICE,
            showarrow=False,
            xref="paper",
            yref="paper",
            x=0.5,
            y=-0.14,
        )
    fig.update_layout(**layout)
    return fig


def _gap_fraction(aln: Alignment) -> list[float]:
    """Percent of gap characters in each column (used for V2 class shading)."""
    total_rows = len(aln.rows)
    if not total_rows:
        return []
    return [
        100.0
        * sum(row[index] == config.GAP_CHAR for row in aln.rows)
        / total_rows
        for index in range(aln.length)
    ]


def _consensus_char(column: list[str]) -> str:
    """Majority standard residue in a column; ties break alphabetically."""
    counts = Counter(char for char in column if char in config.STANDARD_AMINO_ACIDS)
    if not counts:
        return constants.ALN_CONSENSUS_NONE
    top = max(counts.values())
    return min(char for char, count in counts.items() if count == top)


def _property_class(residue: str) -> str:
    """CSS class for one residue, by amino-acid property."""
    for group, residues in constants.ALN_PROPERTY_RESIDUES.items():
        if residue in residues:
            return constants.ALN_PROPERTY_CLASSES[group]
    return constants.ALN_PROPERTY_CLASSES["special"]


def _escaped_row(label: str, content: str) -> str:
    """One labelled monospace row; ``content`` already contains safe spans."""
    return (
        f'<div class="tp53-row">'
        f'<span class="tp53-label">{escape(label, quote=True)}</span>'
        f"{content}</div>"
    )


def alignment_html(aln: Alignment, identity: list[float], block: int = 60) -> str:
    """Render V2 as a standalone HTML string in monospace 60-column blocks.

    Sequence text, IDs and labels are HTML-escaped. Residues are coloured by
    amino-acid property; the consensus row shows the majority standard residue
    per column; the identity row is shaded by conservation class and carries
    the numeric identity in its ``title`` attribute. The ruler always uses
    1-based reference numbering via :attr:`Alignment.col_to_ref`.

    Args:
        aln: parsed (and ideally validated) alignment.
        identity: percent identity per alignment column.
        block: number of columns per HTML block (60 by default).

    Returns:
        HTML fragment, or an empty string when the alignment has no columns.
    """
    if aln is None or not aln.rows or aln.length == 0:
        return ""
    block_size = max(_FULL_WINDOW, int(block))
    positions = _reference_positions(aln)
    gaps = _gap_fraction(aln)
    params = config.ConservationParams()
    classes = classify_columns(list(identity), gaps, params)

    css = [
        "<style>",
        ".tp53-alignment {font-family: monospace; font-size: 12px;}",
        ".tp53-block {margin-bottom: 1.2em;}",
        ".tp53-block-header {font-weight: bold; margin-bottom: 0.3em;}",
        '.tp53-row {white-space: pre; line-height: 1.25;}',
        ".tp53-label {display: inline-block; width: 24ch; color: "
        f"{constants.COLOR_PALETTE['black']};}}",
    ]
    for group, css_class in constants.ALN_PROPERTY_CLASSES.items():
        css.append(
            f".{css_class} {{color: {constants.ALN_PROPERTY_COLORS[group]}; "
            "font-weight: bold;}"
        )
    for label, color in constants.CLASS_COLORS.items():
        css.append(
            f".{constants.ALN_CLASS_CSS[label]} {{background-color: {color}; "
            f"color: {constants.COLOR_PALETTE['black']};}}"
        )
    css.append("</style>")

    parts = [
        f'<div class="{constants.ALN_ROOT_CLASS}">',
        *css,
        f"<h3>{escape(constants.FIG_TITLE_ALIGNMENT, quote=True)}</h3>",
        f"<p><strong>{escape(constants.ALN_LEGEND_TITLE, quote=True)}:</strong> "
        + ", ".join(
            f"{escape(group, quote=True)} "
            f"({escape(constants.ALN_PROPERTY_COLORS[group], quote=True)})"
            for group in constants.ALN_PROPERTY_CLASSES
        )
        + "</p>",
        f"<p><strong>{escape(constants.ALN_LEGEND_IDENTITY, quote=True)}:</strong> "
        + ", ".join(
            f"{escape(label, quote=True)} "
            f"{constants.ALN_CLASS_SYMBOLS[label]}"
            for label in constants.CONSERVATION_CLASSES
        )
        + "</p>",
    ]

    for start in range(0, aln.length, block_size):
        end = min(start + block_size, aln.length)
        columns = list(range(start, end))
        refs = [positions[index] for index in columns]
        mapped = [position for position in refs if position is not None]
        header = constants.ALN_BLOCK_HEADER.format(
            start=start + 1,
            end=end,
            ref_start=mapped[0] if mapped else "n/a",
            ref_end=mapped[-1] if mapped else "n/a",
        )
        parts.append('<div class="tp53-block">')
        parts.append(
            f'<div class="tp53-block-header">{escape(header, quote=True)}</div>'
        )

        ticks = [" "] * len(columns)
        units = [" "] * len(columns)
        for offset, position in enumerate(refs):
            if position is None:
                continue
            units[offset] = str(position % _TENS_STEP)
            if position % _TENS_STEP == 0:
                ticks[offset] = str((position // _TENS_STEP) % _TENS_STEP)
        parts.append(_escaped_row(constants.ALN_LABEL_RULER_TENS, escape("".join(ticks))))
        parts.append(_escaped_row(constants.ALN_LABEL_RULER_TENS, escape("".join(units))))

        for identifier, row in zip(aln.ids, aln.rows, strict=False):
            residues = "".join(
                f'<span class="{_property_class(char)}">{escape(char, quote=True)}</span>'
                for char in row[start:end]
            )
            parts.append(_escaped_row(str(identifier), residues))

        consensus = "".join(
            f'<span class="{_property_class(char)}">{escape(char, quote=True)}</span>'
            for char in (
                _consensus_char([row[index] for row in aln.rows]) for index in columns
            )
        )
        parts.append(_escaped_row(constants.ALN_LABEL_CONSENSUS, consensus))

        identity_spans = []
        for index in columns:
            label = classes[index] if index < len(classes) else constants.CONSERVATION_CLASSES[-1]
            symbol = constants.ALN_CLASS_SYMBOLS[label]
            value = identity[index] if index < len(identity) else 0.0
            identity_spans.append(
                f'<span class="{constants.ALN_CLASS_CSS[label]}" '
                f'title="{value:.1f}% identity">{escape(symbol, quote=True)}</span>'
            )
        parts.append(
            _escaped_row(constants.ALN_LABEL_IDENTITY, "".join(identity_spans))
        )
        parts.append("</div>")

    parts.append("</div>")
    return "\n".join(parts) + "\n"


def blast_overview(selected: list[SelectedHit]) -> go.Figure:
    """Build the V3 per-species bar chart of identity and bit score.

    Species keep the order in which they were supplied, so the figure is
    deterministic. Values are written directly next to every bar.

    Args:
        selected: selected hits in pipeline order.

    Returns:
        Plotly figure with two horizontal bar traces.
    """
    if not selected:
        return _empty_figure(constants.FIG_TITLE_BLAST)
    # The query stands in for its own organism and has no BLAST bit score
    # (§9.1). Label it explicitly instead of drawing a misleading zero bar.
    is_query = [hit.reason == constants.SELECT_REASON_QUERY for hit in selected]
    species = [
        f"{hit.species} ({constants.FIG_QUERY_NO_BLAST})" if own else hit.species
        for hit, own in zip(selected, is_query, strict=True)
    ]
    identities = [hit.percent_identity for hit in selected]
    bit_scores = [
        None if own else hit.hit.best_hsp.bit_score
        for hit, own in zip(selected, is_query, strict=True)
    ]
    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            y=species,
            x=identities,
            orientation="h",
            name=constants.FIG_TRACE_IDENTITY,
            marker_color=constants.COLOR_PALETTE["blue"],
            text=[f"{value:.1f}%" for value in identities],
            textposition="outside",
            hovertemplate="%{y}<br>Identity: %{x:.1f}%<extra></extra>",
        )
    )
    fig.add_trace(
        go.Bar(
            y=species,
            x=bit_scores,
            orientation="h",
            xaxis="x2",
            name=constants.FIG_AXIS_BIT_SCORE,
            marker_color=constants.COLOR_PALETTE["bluish_green"],
            text=[None if value is None else f"{value:.1f}" for value in bit_scores],
            textposition="outside",
            hovertemplate="%{y}<br>Bit score: %{x:.1f}<extra></extra>",
        )
    )
    fig.update_layout(
        title=constants.FIG_TITLE_BLAST,
        xaxis={
            "title": constants.FIG_AXIS_PERCENT_IDENTITY,
            "range": [_Y_MIN, _Y_MAX],
            "domain": [0.0, 0.46],
        },
        xaxis2={"title": constants.FIG_AXIS_BIT_SCORE, "domain": [0.54, 1.0]},
        yaxis={
            "title": constants.FIG_AXIS_SPECIES,
            "categoryorder": "array",
            "categoryarray": species,
        },
        margin={"b": 60, "t": 60},
    )
    return fig


def feature_bars(result: PipelineResult) -> go.Figure:
    """Build the V4 feature comparison chart from the feature-summary table.

    The second bar series is the share of fully conserved columns. When the
    query has no UniProt features (bare FASTA), a notice is shown instead of
    pretending that feature bars are meaningful.

    Args:
        result: completed pipeline result.

    Returns:
        Plotly figure with mean-identity and fully-conserved bars.
    """
    conservation = result.conservation
    if conservation is None or conservation.feature_summary is None:
        return _empty_figure(constants.FIG_TITLE_FEATURES)
    table = conservation.feature_summary
    if table is None or len(table) == 0:
        return _empty_figure(constants.FIG_TITLE_FEATURES)

    labels = [str(value) for value in table["feature"]]
    means = [float(value) for value in table["mean_identity"]]
    conserved = [float(value) for value in table["fully_conserved_pct"]]
    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=labels,
            y=means,
            name=constants.FIG_TRACE_MEAN_IDENTITY,
            marker_color=constants.COLOR_PALETTE["blue"],
            text=[f"{value:.1f}%" for value in means],
            textposition="outside",
            hovertemplate="%{x}<br>Mean identity: %{y:.1f}%<extra></extra>",
        )
    )
    fig.add_trace(
        go.Bar(
            x=labels,
            y=conserved,
            name=constants.FIG_TRACE_FULLY_CONSERVED,
            marker_color=constants.COLOR_PALETTE["bluish_green"],
            text=[f"{value:.1f}%" for value in conserved],
            textposition="outside",
            hovertemplate="%{x}<br>Fully conserved: %{y:.1f}%<extra></extra>",
        )
    )
    if result.query is None or not result.query.features:
        fig.add_annotation(
            text=constants.FIG_NO_FEATURES_NOTICE,
            showarrow=False,
            xref="paper",
            yref="paper",
            x=0.5,
            y=1.12,
        )
    fig.update_layout(
        title=constants.FIG_TITLE_FEATURES,
        xaxis={"title": constants.FIG_AXIS_FEATURE},
        yaxis={
            "title": constants.FIG_AXIS_PERCENT_IDENTITY,
            "range": [_Y_MIN, _Y_MAX],
        },
        barmode="group",
        margin={"b": 100, "t": 100},
    )
    return fig


def identity_heatmap(result: PipelineResult) -> go.Figure:
    """Build the V5 symmetric pairwise-identity heatmap with numeric labels.

    Args:
        result: completed pipeline result.

    Returns:
        Plotly heatmap figure with a colour legend and printed cell values.
    """
    conservation = result.conservation
    if (
        conservation is None
        or conservation.pairwise is None
        or conservation.pairwise.empty
    ):
        return _empty_figure(constants.FIG_TITLE_HEATMAP)
    matrix = conservation.pairwise
    labels = [str(label) for label in matrix.columns]
    values = matrix.to_numpy()
    fig = go.Figure(
        go.Heatmap(
            z=values,
            x=labels,
            y=[str(label) for label in matrix.index],
            zmin=_Y_MIN,
            zmax=_Y_MAX,
            colorscale=list(_HEATMAP_COLORSCALE),
            showscale=True,
            colorbar={"title": constants.FIG_COLORBAR_IDENTITY},
            text=[[f"{float(value):.1f}" for value in row] for row in values],
            texttemplate="%{text}",
            hovertemplate=(
                "%{y} / %{x}<br>Identity: %{z:.1f}%<extra></extra>"
            ),
        )
    )
    fig.update_layout(
        title=constants.FIG_TITLE_HEATMAP,
        xaxis={"title": constants.FIG_AXIS_SEQUENCE, "tickangle": -45, "side": "bottom"},
        yaxis={"title": constants.FIG_AXIS_SEQUENCE, "autorange": "reversed"},
        margin={"b": 120, "l": 140, "t": 60},
    )
    return fig
