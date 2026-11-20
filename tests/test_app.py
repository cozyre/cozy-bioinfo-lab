"""Offline UI tests driven by ``streamlit.testing.v1.AppTest``.

The pipeline and plotting functions are stubbed: these tests assert that the
app wires widgets, run states, tabs and downloads correctly without calling a
remote service.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pandas as pd
import plotly.graph_objects as go
import pytest
from streamlit.testing.v1 import AppTest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


def _install_import_stubs() -> None:
    """Provide temporary modules only when the parallel core work has not landed."""
    if importlib.util.find_spec("tp53.pipeline") is None:
        module = ModuleType("tp53.pipeline")

        def run_pipeline(inputs: Any, clients: Any, progress_cb: Any) -> Any:
            raise RuntimeError("tp53.pipeline is not implemented yet")

        module.run_pipeline = run_pipeline  # type: ignore[attr-defined]
        sys.modules["tp53.pipeline"] = module

    if importlib.util.find_spec("tp53.viz") is None:
        module = ModuleType("tp53.viz")

        def figure(*args: Any, **kwargs: Any) -> go.Figure:
            return go.Figure()

        module.conservation_profile = figure  # type: ignore[attr-defined]
        module.blast_overview = figure  # type: ignore[attr-defined]
        module.feature_bars = figure  # type: ignore[attr-defined]
        module.identity_heatmap = figure  # type: ignore[attr-defined]
        module.alignment_html = (  # type: ignore[attr-defined]
            lambda aln, identity, block=60: "<div>alignment viewer</div>"
        )
        sys.modules["tp53.viz"] = module


_install_import_stubs()

import app.pages as pages  # noqa: E402  (stubs must be installed first)
from tp53 import constants  # noqa: E402
from tp53.constants import FEATURE_SUMMARY_COLUMNS  # noqa: E402
from tp53.conservation import (  # noqa: E402
    SUMMARY_ALIGNMENT_LENGTH,
    SUMMARY_FULLY_CONSERVED_COUNT,
    SUMMARY_FULLY_CONSERVED_PCT,
    SUMMARY_LONGEST_REGION,
    SUMMARY_MEAN_IDENTITY,
    SUMMARY_REGION_COUNT,
)
from tp53.models import (  # noqa: E402
    Alignment,
    BlastHit,
    BlastHsp,
    BlastResult,
    ConservationResult,
    ConservedRegion,
    DomainFeature,
    PipelineResult,
    QueryInfo,
    SelectedHit,
    SequenceRecord,
)


SEQUENCE = "MEEPQSDPSVEPPLSQETFSDLWKLLPENN"


def _sample_result() -> PipelineResult:
    """Small, fully offline result with all artefacts the UI can render."""
    query_record = SequenceRecord(
        id="P04637",
        accession="P04637",
        organism="Homo sapiens",
        description="Cellular tumor antigen p53",
        sequence=SEQUENCE,
    )
    query = QueryInfo(
        record=query_record,
        features=[DomainFeature("Domain", "DNA-binding core", 5, 15)],
    )
    hsp = BlastHsp(
        bit_score=120.0,
        evalue=1e-40,
        identity=len(SEQUENCE),
        positives=len(SEQUENCE),
        gaps=0,
        align_len=len(SEQUENCE),
        query_start=1,
        query_end=len(SEQUENCE),
        hit_start=1,
        hit_end=len(SEQUENCE),
        qseq=SEQUENCE,
        hseq=SEQUENCE,
        midline="|" * len(SEQUENCE),
    )
    hit = BlastHit(
        accession="NP_000546.2",
        title="cellular tumor antigen p53 [Homo sapiens]",
        organism="Homo sapiens",
        length=len(SEQUENCE),
        hsps=[hsp],
    )
    selected = SelectedHit(
        species="Homo sapiens",
        hit=hit,
        query_coverage=1.0,
        percent_identity=100.0,
        reason="highest bit score; NP_ accession preferred",
        candidates=[hit],
    )
    alignment = Alignment(
        ids=[query_record.accession, hit.accession],
        rows=[SEQUENCE, SEQUENCE],
        length=len(SEQUENCE),
        ref_index=0,
        col_to_ref=list(range(1, len(SEQUENCE) + 1)),
    )
    identity = [100.0] * len(SEQUENCE)
    gaps = [0.0] * len(SEQUENCE)
    classes = ["fully conserved"] * len(SEQUENCE)
    region = ConservedRegion(1, len(SEQUENCE), 1, len(SEQUENCE), len(SEQUENCE), 100.0)
    pairwise = pd.DataFrame(
        [[100.0, 100.0], [100.0, 100.0]],
        index=alignment.ids,
        columns=alignment.ids,
    )
    feature_summary = pd.DataFrame(
        [
            {
                "feature": "DNA-binding core",
                "type": "Domain",
                "start": 5,
                "end": 15,
                "columns": 11,
                "mean_identity": 100.0,
                "fully_conserved_pct": 100.0,
            }
        ],
        columns=list(FEATURE_SUMMARY_COLUMNS),
    )
    summary = {
        SUMMARY_MEAN_IDENTITY: 100.0,
        SUMMARY_FULLY_CONSERVED_COUNT: len(SEQUENCE),
        SUMMARY_FULLY_CONSERVED_PCT: 100.0,
        SUMMARY_REGION_COUNT: 1,
        SUMMARY_LONGEST_REGION: len(SEQUENCE),
        SUMMARY_ALIGNMENT_LENGTH: len(SEQUENCE),
    }
    conservation = ConservationResult(
        identity=identity,
        gap_fraction=gaps,
        classes=classes,
        regions=[region],
        pairwise=pairwise,
        feature_summary=feature_summary,
        summary=summary,
    )
    blast = BlastResult(
        query_id="P04637",
        query_len=len(SEQUENCE),
        database="refseq_protein",
        program="blastp",
        hits=[hit],
        species="Homo sapiens",
    )
    return PipelineResult(
        query=query,
        blast_results={"Homo sapiens": blast},
        selected=[selected],
        input_records=[query_record, hit],
        alignment=alignment,
        conservation=conservation,
        run_metadata={"query": "P04637", "database": "refseq_protein"},
        warnings=[],
        raw_blast_xml={"Homo sapiens": "<BlastOutput/>"},
        raw_msa_fasta=f">P04637\n{SEQUENCE}\n",
    )


def _patch_pipeline(monkeypatch: pytest.MonkeyPatch) -> list[tuple[str, str]]:
    """Replace network/plot calls with deterministic offline fakes."""
    events: list[tuple[str, str]] = []

    def fake_run(inputs: Any, clients: Any, progress_cb: Any) -> PipelineResult:
        progress_cb("validate", "input accepted")
        progress_cb("blast", "submitted: Homo sapiens (ABC123XYZ)")
        progress_cb("done", "finished")
        return _sample_result()

    monkeypatch.setattr(pages, "run_pipeline", fake_run)
    for name in (
        "conservation_profile",
        "blast_overview",
        "feature_bars",
        "identity_heatmap",
    ):
        monkeypatch.setattr(pages, name, lambda *args, **kwargs: go.Figure())
    monkeypatch.setattr(
        pages, "alignment_html", lambda *args, **kwargs: "<div>alignment viewer</div>"
    )
    return events


def _app() -> AppTest:
    """Create an AppTest for the repository entry point."""
    return AppTest.from_file(str(ROOT / "streamlit_app.py"), default_timeout=20)


def _sidebar_button(at: AppTest, label: str) -> Any:
    """Find a uniquely labelled sidebar button."""
    matches = [button for button in at.sidebar.button if button.label == label]
    assert len(matches) == 1, f"expected one {label!r} button, found {len(matches)}"
    return matches[0]


def _load_example(at: AppTest) -> None:
    """Click Load example and let the app rerun."""
    _sidebar_button(at, "Load example").click().run()


def _messages(at: AppTest, kind: str) -> list[str]:
    """Collect widget messages from both the sidebar and main area."""
    main = [str(getattr(item, "value", "")) for item in getattr(at, kind)]
    sidebar = [str(getattr(item, "value", "")) for item in getattr(at.sidebar, kind)]
    return main + sidebar


def test_load_example_renders_all_tabs_and_summary_metrics(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """FR-33 example mode reaches every tab with FR-23 metrics filled."""
    _patch_pipeline(monkeypatch)
    at = _app()
    at.run()
    assert not at.exception

    _load_example(at)
    assert not at.exception
    assert [tab.label for tab in at.tabs] == list(pages.TABS)
    assert len(at.metric) >= 5
    assert {metric.label for metric in at.metric} >= {
        "Mean identity",
        "Fully conserved positions",
        "Fully conserved share",
        "Longest conserved region",
        "Conserved regions",
    }


def test_invalid_fasta_blocks_run_with_specific_message(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A non-FASTA query is rejected before the pipeline is called."""
    called: list[Any] = []
    _patch_pipeline(monkeypatch)
    original = pages.run_pipeline

    def counted(inputs: Any, clients: Any, progress_cb: Any) -> PipelineResult:
        called.append(inputs)
        return original(inputs, clients, progress_cb)

    monkeypatch.setattr(pages, "run_pipeline", counted)
    at = _app()
    at.run()
    query_mode = [radio for radio in at.sidebar.radio if radio.label == "Query input"][0]
    query_mode.set_value(pages.QUERY_MODE_FASTA).run()
    email = [field for field in at.sidebar.text_input if field.label == "Contact email"][0]
    email.set_value("user@example.org").run()
    fasta = at.sidebar.text_area[0]
    fasta.set_value("this is not a FASTA record").run()

    _sidebar_button(at, "Run analysis").click().run()

    assert called == []
    assert any(constants.ERR_NOT_FASTA in message for message in _messages(at, "error"))


def test_downloads_are_present_in_example_mode(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """FR-32 artefacts are downloadable from the example snapshot."""
    _patch_pipeline(monkeypatch)
    at = _app()
    at.run()
    _load_example(at)

    labels = [button.label for button in at.download_button]
    assert any(label.startswith("BLAST XML") for label in labels)
    assert "MSA FASTA" in labels
    assert "Identity CSV" in labels


def test_rids_are_saved_and_reused(monkeypatch: pytest.MonkeyPatch) -> None:
    """§12: BLAST RIDs from progress events survive into the next run."""
    seen: list[dict[str, str]] = []

    def fake_run(inputs: Any, clients: Any, progress_cb: Any) -> PipelineResult:
        seen.append(dict(inputs.existing_rids or {}))
        progress_cb("blast", "submitted: Homo sapiens (RID123)")
        return _sample_result()

    monkeypatch.setattr(pages, "run_pipeline", fake_run)
    for name in (
        "conservation_profile",
        "blast_overview",
        "feature_bars",
        "identity_heatmap",
    ):
        monkeypatch.setattr(pages, name, lambda *args, **kwargs: go.Figure())
    monkeypatch.setattr(
        pages, "alignment_html", lambda *args, **kwargs: "<div>alignment viewer</div>"
    )

    at = _app()
    at.run()
    email = [field for field in at.sidebar.text_input if field.label == "Contact email"][0]
    email.set_value("user@example.org").run()

    _sidebar_button(at, "Run analysis").click().run()
    assert at.session_state["rids"] == {"Homo sapiens": "RID123"}

    _sidebar_button(at, "Run analysis").click().run()
    assert seen[0] == {}
    assert seen[1] == {"Homo sapiens": "RID123"}


def test_reset_returns_to_idle(monkeypatch: pytest.MonkeyPatch) -> None:
    """Reset clears the stored result and restores the idle prompt."""
    _patch_pipeline(monkeypatch)
    at = _app()
    at.run()
    _load_example(at)
    assert at.session_state["run_state"] == "success"

    _sidebar_button(at, "Reset").click().run()

    assert at.session_state["run_state"] == "idle"
    assert at.session_state["result"] is None
    assert any("Load example" in message for message in _messages(at, "info"))


def test_example_snapshot_runs_through_real_pipeline_and_viz() -> None:
    """AC-1/AC-6: the bundled snapshot produces a real, renderable result.

    No pipeline or plotting function is stubbed here. Any attempt to contact a
    remote service fails the test, proving example mode is fully offline.
    """
    from app.state import load_example_inputs
    from tp53 import viz
    from tp53.models import PipelineClients
    from tp53.pipeline import run_pipeline

    def forbidden(*args: Any, **kwargs: Any) -> Any:
        raise AssertionError("example mode attempted a remote call")

    clients = PipelineClients(
        forbidden, forbidden, forbidden, forbidden, forbidden, forbidden
    )
    result = run_pipeline(load_example_inputs(), clients, lambda stage, detail: None)

    assert result.run_metadata["mode"] == "example"
    assert [hit.hit.accession for hit in result.selected] == [
        "P04637",
        "XP_001172077.2",
        "NP_035770.2",
        "NP_112251.2",
        "NP_776626.1",
    ]
    assert len(result.input_records) == 5
    assert result.alignment is not None and len(result.alignment.rows) == 5
    assert result.conservation is not None
    assert result.warnings == []
    assert len(result.raw_blast_xml) == 5
    assert result.raw_msa_fasta

    assert viz.conservation_profile(result)
    assert viz.blast_overview(result.selected)
    assert viz.feature_bars(result)
    assert viz.identity_heatmap(result)
    assert "consensus" in viz.alignment_html(
        result.alignment, result.conservation.identity
    ).lower()
