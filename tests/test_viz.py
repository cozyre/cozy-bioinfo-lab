"""Offline rendering tests for all five visualisations (§11, V1-V5)."""

from __future__ import annotations

from math import ceil

import plotly.graph_objects as go
import pytest

from tp53 import constants
from tp53.models import Alignment, PipelineResult
from tp53.msa import parse_alignment
from tp53.viz import (
    alignment_html,
    blast_overview,
    conservation_profile,
    feature_bars,
    identity_heatmap,
)

from tests.test_pipeline import QUERY_INFO, FakeClients, make_inputs, run


def fixture_result(**overrides: object) -> PipelineResult:
    """Build a real result through the fake-client pipeline."""
    result, _events = run(make_inputs(**overrides), FakeClients())
    return result


def empty_result() -> PipelineResult:
    """Minimal result with every optional analysis field missing."""
    return PipelineResult(
        query=None,
        blast_results={},
        selected=[],
        input_records=[],
        alignment=None,
        conservation=None,
    )


def test_conservation_profile_v1() -> None:
    """V1: identity axis, regions, threshold, features and x-axis toggle."""
    result = fixture_result()
    assert result.conservation is not None
    fig = conservation_profile(result)

    names = [trace.name for trace in fig.data]
    assert constants.FIG_TRACE_IDENTITY in names
    assert constants.FIG_TRACE_THRESHOLD in names
    assert any(
        name and name.startswith("Conserved region") for name in names
    )
    assert constants.FIG_TRACE_FEATURE in names
    assert fig.layout.yaxis.range == (0.0, 100.0)
    assert fig.layout.xaxis.title.text == constants.FIG_AXIS_REFERENCE
    assert len(fig.layout.updatemenus[0].buttons) == 2
    assert {
        button.label for button in fig.layout.updatemenus[0].buttons
    } == {constants.FIG_AXIS_REFERENCE, constants.FIG_AXIS_ALIGNMENT}

    identity_trace = next(
        trace for trace in fig.data if trace.name == constants.FIG_TRACE_IDENTITY
    )
    assert len(identity_trace.x) == len(result.conservation.identity)
    assert len(identity_trace.customdata[0]) == 5
    for label in (
        "Reference position",
        "Residue",
        "Identity",
        "Gap fraction",
    ):
        assert label in (identity_trace.hovertemplate or "")

    # Applied thresholds come from run metadata.
    assert result.run_metadata["conservation_params"]["region_identity_threshold"] == 100.0
    threshold = next(
        trace for trace in fig.data if trace.name == constants.FIG_TRACE_THRESHOLD
    )
    assert threshold.y[0] == 100.0

    alternate = conservation_profile(result, smoothing_window=3, x_axis="alignment")
    assert alternate.layout.xaxis.title.text == constants.FIG_AXIS_ALIGNMENT
    alternate_identity = next(
        trace
        for trace in alternate.data
        if trace.name == constants.FIG_TRACE_IDENTITY
    )
    assert alternate_identity.x[0] == 1.0
    assert any(
        smoothed != raw
        for smoothed, raw in zip(
            alternate_identity.y, identity_trace.y, strict=True
        )
    )


def test_alignment_html_v2() -> None:
    """V2: 60-column blocks, consensus, reference ruler and property colours."""
    result = fixture_result()
    assert result.alignment is not None
    assert result.conservation is not None
    html_text = alignment_html(
        result.alignment, result.conservation.identity, block=60
    )

    expected_blocks = ceil(result.alignment.length / 60)
    assert html_text.count('<div class="tp53-block">') == expected_blocks
    assert constants.ALN_LABEL_CONSENSUS in html_text
    assert constants.ALN_LABEL_IDENTITY in html_text
    assert constants.ALN_LABEL_RULER_TENS in html_text
    first_refs = [
        position
        for position in result.alignment.col_to_ref[:60]
        if position is not None
    ]
    assert constants.ALN_BLOCK_HEADER.format(
        start=1, end=60, ref_start=first_refs[0], ref_end=first_refs[-1]
    ) in html_text
    assert any(
        css_class in html_text
        for css_class in constants.ALN_PROPERTY_CLASSES.values()
    )
    assert f".{constants.ALN_CLASS_CSS['fully conserved']}" in html_text
    assert "100.0% identity" in html_text
    assert "<pre" not in html_text  # labelled rows stay in fixed-width divs


def test_alignment_html_escapes_sequence_and_ids() -> None:
    """All sequence and ID text is HTML-escaped before embedding."""
    identifier = 'seq<&">'
    aln = parse_alignment(
        f">{identifier}\nACDEFGHIKLMNPQRSTVWY\n", identifier, [identifier]
    )
    html_text = alignment_html(aln, [100.0] * aln.length)
    assert "seq&lt;&amp;&quot;&gt;" in html_text
    assert identifier not in html_text
    assert "&lt;" in html_text and "&amp;" in html_text and "&quot;" in html_text


def test_blast_overview_v3() -> None:
    """V3: deterministic species order and direct value labels."""
    result = fixture_result()
    fig = blast_overview(result.selected)
    assert isinstance(fig, go.Figure)
    assert len(fig.data) == 2
    assert [trace.name for trace in fig.data] == [
        constants.FIG_TRACE_IDENTITY,
        constants.FIG_AXIS_BIT_SCORE,
    ]
    for trace in fig.data:
        assert len(trace.y) == len(result.selected)
    # Every selected hit has an identity value; the query row has no bit score.
    assert all(fig.data[0].text)
    query_rows = [
        index
        for index, hit in enumerate(result.selected)
        if hit.reason == constants.SELECT_REASON_QUERY
    ]
    for index, _hit in enumerate(result.selected):
        if index in query_rows:
            assert fig.data[1].x[index] is None
            assert fig.data[1].text[index] is None
        else:
            assert fig.data[1].x[index] is not None
            assert fig.data[1].text[index] is not None
    expected_species = [
        f"{hit.species} ({constants.FIG_QUERY_NO_BLAST})"
        if hit.reason == constants.SELECT_REASON_QUERY
        else hit.species
        for hit in result.selected
    ]
    assert list(fig.data[0].y) == expected_species
    assert fig.layout.xaxis.range == (0.0, 100.0)


def test_feature_bars_v4() -> None:
    """V4: mean identity and fully conserved share from the summary table."""
    result = fixture_result()
    assert result.conservation is not None
    table = result.conservation.feature_summary
    fig = feature_bars(result)

    assert [trace.name for trace in fig.data] == [
        constants.FIG_TRACE_MEAN_IDENTITY,
        constants.FIG_TRACE_FULLY_CONSERVED,
    ]
    assert list(fig.data[0].x) == list(table["feature"])
    assert fig.data[0].y == pytest.approx(list(table["mean_identity"]))
    assert fig.data[1].y == pytest.approx(list(table["fully_conserved_pct"]))
    assert constants.OUTSIDE_FEATURES_LABEL in list(table["feature"])
    assert all(trace.text for trace in fig.data)


def test_identity_heatmap_v5() -> None:
    """V5: symmetric annotated heatmap with a colour legend."""
    result = fixture_result()
    assert result.conservation is not None
    fig = identity_heatmap(result)
    assert isinstance(fig, go.Figure)
    assert len(fig.data) == 1
    heatmap = fig.data[0]
    assert heatmap.type == "heatmap"
    matrix = result.conservation.pairwise
    assert list(heatmap.x) == list(matrix.columns)
    assert list(heatmap.y) == list(matrix.index)
    assert heatmap.texttemplate == "%{text}"
    assert heatmap.showscale is True
    assert heatmap.colorbar.title.text == constants.FIG_COLORBAR_IDENTITY
    assert len(heatmap.text) == len(matrix)
    assert heatmap.text[0][0] == f"{matrix.iloc[0, 0]:.1f}"
    assert all(
        heatmap.z[index][other] == heatmap.z[other][index]
        for index in range(len(matrix))
        for other in range(len(matrix))
    )


def test_empty_result_renders_without_errors() -> None:
    """Every renderer degrades gracefully when no analysis has run."""
    result = empty_result()
    figures = [
        conservation_profile(result),
        blast_overview(result.selected),
        feature_bars(result),
        identity_heatmap(result),
    ]
    for fig in figures:
        assert isinstance(fig, go.Figure)
        assert fig.layout.annotations[0].text == constants.FIG_TITLE_NO_DATA
    assert alignment_html(Alignment([], [], 0, 0), []) == ""


def test_fasta_result_omits_feature_overlays() -> None:
    """FASTA queries are rendered with a notice instead of feature bands."""
    fasta = f">query\n{QUERY_INFO.record.sequence}\n"
    result = fixture_result(
        query_accession=None,
        query_fasta=fasta,
        query_organism="Homo sapiens",
    )
    assert result.query is not None
    assert result.query.features == []

    profile = conservation_profile(result)
    assert any(
        annotation.text == constants.FIG_NO_FEATURES_NOTICE
        for annotation in profile.layout.annotations
    )
    assert not any(
        trace.name == constants.FIG_TRACE_FEATURE for trace in profile.data
    )
    bars = feature_bars(result)
    assert any(
        annotation.text == constants.FIG_NO_FEATURES_NOTICE
        for annotation in bars.layout.annotations
    )
