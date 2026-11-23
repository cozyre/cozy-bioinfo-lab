"""Streamlit pages for the TP53 Conservation Analyzer (§10).

Everything shown in the main area is rendered from the single
:class:`~tp53.models.PipelineResult` held in session state. This module and
:mod:`app.state` are the only places allowed to import Streamlit.
"""

from __future__ import annotations

import re
from collections.abc import Callable
from typing import Any

import pandas as pd
import streamlit as st
from streamlit.delta_generator import DeltaGenerator
from streamlit.errors import StreamlitSecretNotFoundError

from tp53 import config, constants, msa as msa_module
from tp53.blast_parser import parse_blast_xml
from tp53.conservation import (
    SUMMARY_ALIGNMENT_LENGTH,
    SUMMARY_FULLY_CONSERVED_COUNT,
    SUMMARY_FULLY_CONSERVED_PCT,
    SUMMARY_LONGEST_REGION,
    SUMMARY_MEAN_IDENTITY,
    SUMMARY_REGION_COUNT,
)
from tp53.fasta_io import parse_fasta, parse_single_fasta, validate_protein
from tp53.models import (
    AlignmentValidationError,
    BlastHit,
    PipelineInputs,
    PipelineResult,
    ValidationError,
)
from tp53.msa import alignment_statistics
from tp53.pipeline import STAGES as PIPELINE_STAGES
from tp53.pipeline import run_pipeline
from tp53.viz import (
    alignment_html,
    blast_overview,
    conservation_profile,
    feature_bars,
    identity_heatmap,
)

from . import state

# --- UI strings (A1): aliases only; every literal lives in constants.py. ------

TAB_OVERVIEW = constants.TAB_OVERVIEW
TAB_BLAST = constants.TAB_BLAST
TAB_ALIGNMENT = constants.TAB_ALIGNMENT
TAB_CONSERVATION = constants.TAB_CONSERVATION
TAB_DOWNLOADS = constants.TAB_DOWNLOADS
TABS = constants.TABS

QUERY_MODE_ACCESSION = constants.QUERY_MODE_ACCESSION
QUERY_MODE_FASTA = constants.QUERY_MODE_FASTA
QUERY_MODES = constants.QUERY_MODES

SOURCE_REMOTE = constants.SOURCE_REMOTE
SOURCE_UPLOAD = constants.SOURCE_UPLOAD
SOURCE_MODES = constants.SOURCE_MODES

X_AXIS_REFERENCE = constants.UI_X_AXIS_REFERENCE
X_AXIS_COLUMN = constants.UI_X_AXIS_COLUMN
X_AXIS_CHOICES = constants.UI_X_AXIS_CHOICES

# --- Parsing patterns and derived layout values. ------------------------------

EMAIL_RE = re.compile(r"^[^@\s]+@[^@\s]+\.[^@\s]+$")
UNIPROT_ACCESSION_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})"
    r"(?:\.\d+)?$",
    re.IGNORECASE,
)
NCBI_ACCESSION_RE = re.compile(
    r"^(?:NP|XP|WP|YP|ZP)_\d+(?:\.\d+)?$", re.IGNORECASE
)
BYTES_PER_MIB = 1024 * 1024  # unit conversion, not a tunable default

# Derive step numbers from the pipeline's stage order instead of restating it.
STAGE_STEP: dict[str, int] = {
    stage: index + 1
    for index, stage in enumerate(s for s in PIPELINE_STAGES if s != "done")
}
STAGE_STEP["done"] = len(STAGE_STEP)
TOTAL_STEPS = config.PROGRESS_TOTAL_STEPS


def _widget_key(name: str) -> str:
    """Return a generation-scoped widget key so Reset restores defaults."""
    return f"{name}_{st.session_state['widget_generation']}"


def _secret(key: str) -> str:
    """Read an optional secret without raising when secrets are absent."""
    try:
        return str(st.secrets.get(key, ""))
    except (KeyError, TypeError, AttributeError, StreamlitSecretNotFoundError):
        return ""


def _format_limit(byte_count: int) -> str:
    """Human-readable upload limit for validation messages."""
    if byte_count % BYTES_PER_MIB == 0:
        return constants.UPLOAD_LIMIT_MB.format(count=byte_count // BYTES_PER_MIB)
    return constants.UPLOAD_LIMIT_BYTES.format(count=byte_count)


def _read_upload(uploaded: Any, label: str) -> str | None:
    """Decode one uploaded file, enforcing :data:`config.MAX_UPLOAD_BYTES`."""
    if uploaded is None:
        return None
    size = getattr(uploaded, "size", None)
    if isinstance(size, int) and size > config.MAX_UPLOAD_BYTES:
        st.sidebar.error(
            constants.ERR_UPLOAD_TOO_LARGE.format(
                name=getattr(uploaded, "name", label),
                limit=_format_limit(config.MAX_UPLOAD_BYTES),
            )
        )
        return None
    try:
        return uploaded.getvalue().decode("utf-8")
    except UnicodeDecodeError:
        st.sidebar.error(
            constants.ERR_UPLOAD_NOT_TEXT.format(name=getattr(uploaded, "name", label))
        )
        return None


def _validate_accession(value: str) -> tuple[str, str] | None:
    """Return ``(accession, source)`` or ``None`` after an inline error."""
    accession = value.strip()
    if not accession:
        st.sidebar.error(constants.ERR_ACCESSION_REQUIRED)
        return None
    if UNIPROT_ACCESSION_RE.fullmatch(accession):
        return accession, "uniprot"
    if NCBI_ACCESSION_RE.fullmatch(accession):
        return accession, "ncbi"
    st.sidebar.error(constants.ERR_ACCESSION_FORMAT)
    return None


def _validate_fasta_query(text: str) -> ValidationError | None:
    """Validate a pasted/uploaded query FASTA with the shared domain rules."""
    try:
        record = parse_single_fasta(text)
        validate_protein(record)
    except ValidationError as exc:
        return exc
    return None


def _prevalidate_msa(text: str) -> str | None:
    """Check an uploaded alignment before any network call (§12).

    Row count, row lengths and the allowed alphabet are checked by the shared
    validator, so its error names the offending sequence. The hard de-gap
    check is repeated in the pipeline once the expected sequences are known.
    """
    try:
        records = parse_fasta(text)
        if not records:
            raise AlignmentValidationError(constants.ERR_EMPTY_ALIGNMENT)
        ids = [record.id for record in records]
        rows = [record.sequence for record in records]
        alignment = msa_module.parse_alignment(text, ref_id=ids[0], expected_ids=ids)
        msa_module.validate_alignment(alignment, ids, rows)
    except ValidationError as exc:
        return str(exc)
    return None


def _make_progress_callback(status: Any) -> Callable[[str, str], None]:
    """Build a progress callback that updates the step indicator.

    Job IDs are saved by the client wrappers in :mod:`app.state` the moment a
    service issues them, so progress text is display-only (§12).
    """

    def callback(stage: str, detail: str) -> None:
        state.add_progress(stage, detail)
        step = STAGE_STEP.get(stage, 1)
        status.update(
            label=constants.STATUS_STEP.format(
                step=step, total=TOTAL_STEPS, detail=detail
            ),
            state="running",
            expanded=True,
        )
        status.write(constants.PROGRESS_LINE.format(stage=stage, detail=detail))

    return callback


def _failure_message(error: Exception) -> str:
    """Build a stage-aware error with the right upload fallback (§12)."""
    events = st.session_state.get("progress", [])
    stage = events[-1][0] if events else "validate"
    if stage in {"blast", "select"}:
        hint = constants.UPLOAD_HINT_BLAST
    elif stage in {"msa", "conservation"}:
        hint = constants.UPLOAD_HINT_MSA
    else:
        hint = constants.UPLOAD_HINT_GENERIC
    return constants.STAGE_ACTION.format(
        stage=constants.STAGE_LABELS.get(stage, stage),
        error=error,
        hint=hint,
    )


def _execute(inputs: PipelineInputs, status_slot: DeltaGenerator) -> None:
    """Run the pipeline once, updating run state and progress in the UI."""
    state.set_running()
    state.save_inputs(inputs)
    try:
        with status_slot.status(
            constants.STATUS_VALIDATING.format(
                step=STAGE_STEP["validate"], total=TOTAL_STEPS
            ),
            expanded=True,
        ) as status:
            clients = state.build_clients(
                inputs.email, inputs.api_key, blast_params=inputs.blast_params
            )
            result = run_pipeline(inputs, clients, _make_progress_callback(status))
            status.update(
                label=constants.STATUS_COMPLETE, state="complete", expanded=False
            )
        warnings = list(result.warnings or [])
        state.set_result(result, partial_success=bool(warnings))
        if warnings:
            st.warning(constants.ANALYSIS_WARNINGS)
        else:
            st.success(constants.ANALYSIS_COMPLETE)
    except Exception as exc:  # surface every remote or parsing failure in the UI
        message = _failure_message(exc)
        state.set_failed(message)
        st.error(message)


def _load_example(status_slot: DeltaGenerator) -> None:
    """Run the bundled offline snapshot (FR-33)."""
    state.clear_rids()
    st.session_state["query_mode"] = QUERY_MODE_ACCESSION
    st.session_state["query_signature"] = (
        config.DEFAULT_QUERY_ACCESSION,
        None,
        tuple(config.SPECIES),
    )
    _execute(state.load_example_inputs(), status_slot)


def _build_inputs(
    query_mode: str,
    accession: str,
    fasta_text: str,
    query_organism: str,
    species: list[str],
    email: str,
    api_key: str | None,
    source_mode: str,
    use_uploaded_xml: bool,
    use_uploaded_msa: bool,
    xml_uploads: dict[str, Any],
    msa_upload: Any,
    blast_database: str,
    evalue_cutoff: float,
    min_coverage_pct: float,
    conservation_threshold: float,
    min_region_length: int,
    smoothing_window: int,
    msa_job_id: str | None,
) -> PipelineInputs | None:
    """Validate sidebar values and assemble :class:`PipelineInputs`."""
    if not EMAIL_RE.fullmatch(email.strip()):
        st.sidebar.error(
            constants.ERR_EMAIL_REQUIRED
            if not email.strip()
            else constants.ERR_EMAIL_FORMAT
        )
        return None
    if not species:
        st.sidebar.error(constants.ERR_SPECIES_REQUIRED)
        return None

    query_accession: str | None = None
    accession_source = "uniprot"
    query_fasta: str | None = None
    organism = query_organism.strip()

    if query_mode == QUERY_MODE_ACCESSION:
        validated = _validate_accession(accession)
        if validated is None:
            return None
        query_accession, accession_source = validated
    else:
        if not fasta_text.strip():
            st.sidebar.error(constants.ERR_FASTA_REQUIRED)
            return None
        problem = _validate_fasta_query(fasta_text)
        if problem is not None:
            st.sidebar.error(str(problem))
            return None
        if not organism:
            st.sidebar.error(constants.ERR_ORGANISM_REQUIRED)
            return None
        query_fasta = fasta_text

    blast_xml_by_species: dict[str, str] | None = None
    if source_mode == SOURCE_UPLOAD and use_uploaded_xml:
        blast_xml_by_species = {}
        for target in species:
            text = _read_upload(
                xml_uploads.get(target),
                constants.DOWNLOAD_BLAST_XML.format(species=target),
            )
            if text is None:
                return None
            try:
                parse_blast_xml(text, target, require_hits=True)
            except Exception as exc:
                st.sidebar.error(
                    constants.ERR_XML_PARSE.format(species=target, detail=exc)
                )
                return None
            blast_xml_by_species[target] = text
        if not blast_xml_by_species:
            st.sidebar.error(constants.ERR_XML_UPLOAD_REQUIRED)
            return None

    msa_fasta: str | None = None
    if source_mode == SOURCE_UPLOAD and use_uploaded_msa:
        msa_fasta = _read_upload(msa_upload, constants.LABEL_MSA_UPLOAD)
        if msa_fasta is None:
            return None
        problem = _prevalidate_msa(msa_fasta)
        if problem is not None:
            st.sidebar.error(problem)
            return None

    query_signature = (query_accession, query_fasta, tuple(species))
    if st.session_state.get("query_signature") != query_signature:
        # A different query/species set must not resume another run's RIDs.
        state.clear_rids()
        st.session_state["query_signature"] = query_signature

    return PipelineInputs(
        query_accession=query_accession,
        accession_source=accession_source,  # type: ignore[arg-type]
        query_fasta=query_fasta,
        query_organism=organism,
        species=tuple(species),
        email=email.strip(),
        api_key=api_key or None,
        blast_params=config.BlastParams(
            database=blast_database, expect=evalue_cutoff
        ),
        selection_params=config.SelectionParams(
            evalue_cutoff=evalue_cutoff,
            min_query_coverage=min_coverage_pct / 100.0,
        ),
        conservation_params=config.ConservationParams(
            region_identity_threshold=conservation_threshold,
            min_region_length=min_region_length,
        ),
        blast_xml_by_species=blast_xml_by_species,
        msa_fasta=msa_fasta,
        existing_rids=dict(st.session_state["rids"]) or None,
        msa_job_id=msa_job_id,
        smoothing_window=smoothing_window,
    )


def _render_sidebar(status_slot: DeltaGenerator) -> None:
    """Render inputs, advanced settings and the run/load/resume/reset buttons."""
    st.sidebar.header(constants.SIDEBAR_QUERY)
    query_mode = st.sidebar.radio(
        constants.LABEL_QUERY_INPUT,
        list(QUERY_MODES),
        key=_widget_key("query_mode"),
        horizontal=True,
    )
    st.session_state["query_mode"] = query_mode

    accession = st.sidebar.text_input(
        constants.LABEL_ACCESSION,
        value=config.DEFAULT_QUERY_ACCESSION,
        key=_widget_key("accession"),
        disabled=query_mode != QUERY_MODE_ACCESSION,
        help=constants.HELP_ACCESSION,
    )
    if query_mode == QUERY_MODE_FASTA:
        fasta_text = st.sidebar.text_area(
            constants.LABEL_FASTA, key=_widget_key("fasta_text"), height=110
        )
        fasta_upload = st.sidebar.file_uploader(
            constants.LABEL_FASTA_UPLOAD,
            type=["fa", "fasta", "faa", "txt"],
            key=_widget_key("fasta_upload"),
        )
        if fasta_upload is not None:
            uploaded_text = _read_upload(fasta_upload, constants.LABEL_FASTA)
            if uploaded_text:
                fasta_text = uploaded_text
        query_organism = st.sidebar.text_input(
            constants.LABEL_QUERY_ORGANISM, key=_widget_key("query_organism")
        )
    else:
        fasta_text = ""
        query_organism = st.sidebar.text_input(
            constants.LABEL_QUERY_ORGANISM,
            value=constants.VALUE_ORGANISM_FROM_ACCESSION,
            key=_widget_key("query_organism_accession"),
            disabled=True,
        )

    st.sidebar.header(constants.SIDEBAR_CONTACT)
    email = st.sidebar.text_input(
        constants.LABEL_EMAIL,
        value=_secret("NCBI_EMAIL"),
        key=_widget_key("email"),
        help=constants.HELP_EMAIL,
    )
    api_key = _secret("NCBI_API_KEY") or None

    st.sidebar.header(constants.SIDEBAR_TARGETS)
    species = st.sidebar.multiselect(
        constants.LABEL_SPECIES,
        list(config.SPECIES),
        default=list(config.SPECIES),
        key=_widget_key("species"),
    )
    source_mode = st.sidebar.radio(
        constants.LABEL_SOURCE,
        list(SOURCE_MODES),
        key=_widget_key("source_mode"),
    )
    use_uploaded_xml = False
    use_uploaded_msa = False
    xml_uploads: dict[str, Any] = {}
    msa_upload: Any = None
    if source_mode == SOURCE_UPLOAD:
        use_uploaded_xml = st.sidebar.checkbox(
            constants.LABEL_USE_XML,
            value=True,
            key=_widget_key("use_xml"),
        )
        if use_uploaded_xml:
            for target in species:
                xml_uploads[target] = st.sidebar.file_uploader(
                    constants.DOWNLOAD_BLAST_XML.format(species=target),
                    type=["xml"],
                    key=_widget_key(f"xml_{_slug(target)}"),
                )
        use_uploaded_msa = st.sidebar.checkbox(
            constants.LABEL_USE_MSA,
            value=True,
            key=_widget_key("use_msa"),
        )
        if use_uploaded_msa:
            msa_upload = st.sidebar.file_uploader(
                constants.LABEL_MSA_UPLOAD,
                type=["fa", "fasta", "faa", "txt"],
                key=_widget_key("msa_upload"),
            )

    with st.sidebar.expander(constants.ADVANCED_SETTINGS, expanded=False):
        blast_database = st.selectbox(
            constants.LABEL_BLAST_DATABASE,
            list(config.BLAST_DATABASES),
            index=0,
            key=_widget_key("blast_database"),
        )
        defaults = config.SelectionParams()
        cons_defaults = config.ConservationParams()
        evalue_cutoff = st.number_input(
            constants.LABEL_EVALUE,
            min_value=0.0,
            value=float(defaults.evalue_cutoff),
            format="%.1e",
            key=_widget_key("evalue"),
        )
        min_coverage_pct = st.number_input(
            constants.LABEL_COVERAGE,
            min_value=0.0,
            max_value=100.0,
            value=defaults.min_query_coverage * 100.0,
            key=_widget_key("coverage"),
        )
        conservation_threshold = st.number_input(
            constants.LABEL_THRESHOLD,
            min_value=0.0,
            max_value=100.0,
            value=cons_defaults.region_identity_threshold,
            key=_widget_key("threshold"),
        )
        min_region_length = st.number_input(
            constants.LABEL_MIN_REGION,
            min_value=1,
            value=cons_defaults.min_region_length,
            key=_widget_key("min_region"),
        )
        smoothing_window = int(
            st.number_input(
                constants.LABEL_SMOOTHING_ADVANCED,
                min_value=1,
                value=config.SMOOTHING_WINDOW_DEFAULT,
                key=_widget_key("smoothing_advanced"),
            )
        )

    def current_inputs() -> PipelineInputs | None:
        """Build inputs from the current widget values, including saved jobs."""
        return _build_inputs(
            query_mode,
            accession,
            fasta_text,
            query_organism,
            species,
            email,
            api_key,
            source_mode,
            use_uploaded_xml,
            use_uploaded_msa,
            xml_uploads,
            msa_upload,
            blast_database,
            evalue_cutoff,
            min_coverage_pct,
            conservation_threshold,
            int(min_region_length),
            smoothing_window,
            state.get_clustalo_job(),
        )

    run_col, example_col, reset_col = st.sidebar.columns(3, gap="small")
    if run_col.button(constants.BUTTON_RUN, key=_widget_key("run"), type="primary"):
        inputs = current_inputs()
        if inputs is not None:
            _execute(inputs, status_slot)
    if example_col.button(constants.BUTTON_EXAMPLE, key=_widget_key("example")):
        _load_example(status_slot)
    if reset_col.button(constants.BUTTON_RESET, key=_widget_key("reset")):
        state.reset_state()
        st.session_state["widget_generation"] += 1
        st.rerun()

    saved_ids = bool(st.session_state.get("rids")) or bool(state.get_clustalo_job())
    resumable = st.session_state.get("run_state") in {state.FAILED, state.RUNNING}
    if saved_ids and resumable:
        if st.sidebar.button(constants.BUTTON_RESUME, key=_widget_key("resume")):
            inputs = current_inputs()
            if inputs is not None:
                _execute(inputs, status_slot)


def _json_safe(value: Any) -> Any:
    """Convert arbitrary metadata into JSON-serialisable primitives."""
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_json_safe(item) for item in value]
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def _class_legend() -> None:
    """Show the conservation-class legend with labels as well as colours."""
    parts = [
        f'<span style="color:{constants.CLASS_COLORS[label]}">{label}</span>'
        for label in constants.CONSERVATION_CLASSES
    ]
    st.markdown(
        constants.CLASSES_PREFIX + " · ".join(parts), unsafe_allow_html=True
    )


def _render_overview(result: PipelineResult) -> None:
    """Render status, run metadata and FR-23 summary metrics."""
    failed = st.session_state.get("run_state") == state.FAILED
    status_label = (
        constants.ANALYSIS_FAILED.format(
            error=st.session_state.get("run_error", "")
        )
        if failed
        else constants.STATUS_COMPLETE
    )
    status_state = "error" if failed else "complete"
    with st.status(status_label, state=status_state, expanded=False) as status:
        events = st.session_state.get("progress", [])
        if events:
            for stage, detail in events:
                status.write(constants.PROGRESS_LINE.format(stage=stage, detail=detail))
        else:
            status.write(constants.ALL_STAGES_COMPLETE)

    if result.run_metadata:
        st.subheader(constants.LABEL_RUN_METADATA)
        st.json(_json_safe(result.run_metadata))

    summary = result.conservation.summary if result.conservation else {}
    mean = summary.get(SUMMARY_MEAN_IDENTITY, 0.0)
    count = summary.get(SUMMARY_FULLY_CONSERVED_COUNT, 0)
    share = summary.get(SUMMARY_FULLY_CONSERVED_PCT, 0.0)
    longest = summary.get(SUMMARY_LONGEST_REGION, 0)
    regions = summary.get(SUMMARY_REGION_COUNT, 0)
    columns = summary.get(SUMMARY_ALIGNMENT_LENGTH, 0)

    st.subheader(constants.LABEL_CONSERVATION_SUMMARY)
    first, second, third = st.columns(3, gap="small")
    with first:
        st.metric(constants.METRIC_MEAN_IDENTITY, f"{float(mean):.1f}%")
    with second:
        st.metric(constants.METRIC_FULLY_CONSERVED, str(count))
    with third:
        st.metric(constants.METRIC_FULLY_SHARE, f"{float(share):.1f}%")
    first, second, third = st.columns(3, gap="small")
    with first:
        st.metric(
            constants.METRIC_LONGEST_REGION,
            constants.UNIT_COLUMNS.format(count=longest),
        )
    with second:
        st.metric(constants.METRIC_REGIONS, str(regions))
    with third:
        st.metric(constants.METRIC_COLUMNS, str(columns))
    st.caption(constants.IDENTITY_NOTE)


def _style_candidates(frame: pd.DataFrame, accession: str) -> Any:
    """Highlight the selected candidate row without relying on colour alone."""

    def highlight(row: pd.Series) -> list[str]:
        selected = str(row.get("accession", "")) == accession
        return [f"background-color: {constants.COLOR_PALETTE['sky_blue']}; font-weight: 600"] * len(row) if selected else [""] * len(row)

    return frame.style.apply(highlight, axis=1)


def _best_hsp(hit: BlastHit) -> Any:
    """Best HSP of a hit, or ``None`` if the hit has no HSPs."""
    return hit.best_hsp if hit.hsps else None


def _render_blast(result: PipelineResult) -> None:
    """Render candidate tables, the overview chart and raw hit alignments."""
    if not result.selected:
        st.info(constants.NO_HITS_NOTICE)
    else:
        figure = blast_overview(result.selected)
        st.plotly_chart(figure, width="stretch", key="blast_overview")

    for selected in result.selected:
        with st.container(border=True):
            st.subheader(
                constants.LABEL_SELECTED_HIT.format(
                    species=selected.species, accession=selected.hit.accession
                )
            )
            st.caption(selected.reason)
            rows: list[dict[str, Any]] = []
            blast = result.blast_results.get(selected.species)
            for hit in selected.candidates or [selected.hit]:
                hsp = _best_hsp(hit)
                span = (
                    hsp.query_end - hsp.query_start + 1 if hsp is not None else 0
                )
                coverage = (
                    span / blast.query_len
                    if blast is not None and blast.query_len
                    else 0.0
                )
                rows.append(
                    {
                        "accession": hit.accession,
                        "title": hit.title,
                        "organism": hit.organism,
                        "bit_score": round(hsp.bit_score, 1) if hsp else 0.0,
                        "evalue": hsp.evalue if hsp else 0.0,
                        "identity_pct": round(hsp.percent_identity, 1) if hsp else 0.0,
                        "coverage_pct": round(coverage * 100.0, 1),
                        "selected": "yes" if hit.accession == selected.hit.accession else "",
                    }
                )
            frame = pd.DataFrame(rows)
            st.dataframe(
                _style_candidates(frame, selected.hit.accession),
                width="stretch",
                hide_index=True,
            )
            with st.expander(
                constants.EXPANDER_RAW_ALIGNMENT.format(
                    species=selected.species, accession=selected.hit.accession
                ),
                expanded=False,
            ):
                for index, hsp in enumerate(selected.hit.hsps, start=1):
                    st.caption(
                        constants.LABEL_HSP.format(
                            index=index,
                            query_start=hsp.query_start,
                            query_end=hsp.query_end,
                        )
                    )
                    st.code(
                        f"Query  {hsp.qseq}\n       {hsp.midline}\nSubject {hsp.hseq}",
                        language=None,
                    )


def _render_alignment(result: PipelineResult) -> None:
    """Render the block-wise alignment viewer and per-sequence statistics."""
    if result.alignment is None:
        st.info(constants.NO_ALIGNMENT_NOTICE)
        return
    alignment = result.alignment
    identity = result.conservation.identity if result.conservation else []
    st.dataframe(
        alignment_statistics(alignment), width="stretch", hide_index=True
    )
    _class_legend()
    html = alignment_html(alignment, identity, block=config.FASTA_LINE_WIDTH)
    st.iframe(html, height=config.ALIGNMENT_VIEWER_HEIGHT_PX)


def _render_regions(result: PipelineResult) -> None:
    """Render the conserved-region table (§9.6, FR-20)."""
    regions = result.conservation.regions if result.conservation else []
    if not regions:
        st.info(constants.NO_REGIONS_NOTICE)
        return
    frame = pd.DataFrame(
        [
            {
                "alignment_start": region.aln_start,
                "alignment_end": region.aln_end,
                "reference_start": region.ref_start,
                "reference_end": region.ref_end,
                "length": region.length,
                "mean_identity": round(region.mean_identity, 1),
            }
            for region in regions
        ]
    )
    st.dataframe(frame, width="stretch", hide_index=True)


def _render_conservation(result: PipelineResult) -> None:
    """Render the profile, regions, feature summary and pairwise heatmap."""
    if result.conservation is None:
        st.info(constants.NO_CONSERVATION_NOTICE)
        return
    _class_legend()
    window_max = max(
        config.SMOOTHING_WINDOW_DEFAULT * config.SMOOTHING_MAX_MULTIPLIER, 1
    )
    smoothing_key = _widget_key("smoothing")
    if smoothing_key not in st.session_state:
        advanced = int(
            st.session_state.get(
                _widget_key("smoothing_advanced"), config.SMOOTHING_WINDOW_DEFAULT
            )
        )
        st.session_state[smoothing_key] = max(1, min(advanced, window_max))
    # The key is seeded from advanced settings above, so no explicit value is
    # passed: Streamlit would otherwise warn about a duplicate default.
    smoothing = st.slider(
        constants.LABEL_SMOOTHING_SLIDER,
        min_value=1,
        max_value=window_max,
        key=smoothing_key,
    )
    axis_choice = st.radio(
        constants.LABEL_X_AXIS,
        list(X_AXIS_CHOICES),
        horizontal=True,
        key=_widget_key("x_axis"),
    )
    x_axis = "reference" if axis_choice == X_AXIS_REFERENCE else "column"
    st.plotly_chart(
        conservation_profile(result, smoothing_window=smoothing, x_axis=x_axis),
        width="stretch",
        key="conservation_profile",
    )

    st.subheader(constants.LABEL_REGIONS)
    _render_regions(result)

    st.subheader(constants.LABEL_FEATURE_COMPARISON)
    query_source = (result.run_metadata or {}).get("query_source")
    if query_source != "fasta":
        st.dataframe(
            result.conservation.feature_summary,
            width="stretch",
            hide_index=True,
        )
        st.plotly_chart(feature_bars(result), width="stretch", key="feature_bars")
    else:
        st.info(constants.FEATURE_SKIP_NOTICE)

    st.subheader(constants.LABEL_PAIRWISE)
    st.plotly_chart(identity_heatmap(result), width="stretch", key="identity_heatmap")


def _slug(value: str) -> str:
    """Filesystem-safe fragment for download file names."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_") or "hit"


def _identity_csv(result: PipelineResult) -> str | None:
    """Build the downloadable per-column identity table (FR-32)."""
    conservation = result.conservation
    if conservation is None:
        return None
    count = len(conservation.identity)
    reference = (
        result.alignment.col_to_ref
        if result.alignment is not None and len(result.alignment.col_to_ref) == count
        else [None] * count
    )
    frame = pd.DataFrame(
        {
            "alignment_column": range(1, count + 1),
            "reference_position": reference,
            "identity_pct": conservation.identity,
            "gap_fraction_pct": conservation.gap_fraction,
            "class": conservation.classes,
        }
    )
    return frame.to_csv(index=False)


def _render_downloads(result: PipelineResult) -> None:
    """Offer raw BLAST XML, MSA FASTA and identity CSV downloads (FR-32)."""
    any_download = False
    for species, xml_text in (result.raw_blast_xml or {}).items():
        any_download = True
        st.download_button(
            constants.DOWNLOAD_BLAST_XML.format(species=species),
            data=xml_text.encode("utf-8"),
            file_name=constants.FILE_BLAST_XML.format(slug=_slug(species)),
            mime=constants.MIME_XML,
            key=_widget_key(f"download_blast_{_slug(species)}"),
        )
    if result.raw_msa_fasta:
        any_download = True
        st.download_button(
            constants.DOWNLOAD_MSA,
            data=result.raw_msa_fasta.encode("utf-8"),
            file_name=constants.FILE_MSA,
            mime=constants.MIME_TEXT,
            key=_widget_key("download_msa"),
        )
    csv_text = _identity_csv(result)
    if csv_text is not None:
        any_download = True
        st.download_button(
            constants.DOWNLOAD_CSV,
            data=csv_text.encode("utf-8"),
            file_name=constants.FILE_IDENTITY,
            mime=constants.MIME_CSV,
            key=_widget_key("download_identity"),
        )
    if not any_download:
        st.info(constants.NO_ARTIFACTS_NOTICE)


def _render_result(result: PipelineResult) -> None:
    """Render the five result tabs from one stored :class:`PipelineResult`."""
    for warning in result.warnings or []:
        st.warning(warning)

    overview, blast, alignment, conservation, downloads = st.tabs(list(TABS))
    with overview:
        _render_overview(result)
    with blast:
        _render_blast(result)
    with alignment:
        _render_alignment(result)
    with conservation:
        _render_conservation(result)
    with downloads:
        _render_downloads(result)


def render_app() -> None:
    """Render the complete page: sidebar, run state and the five tabs."""
    if "widget_generation" not in st.session_state:
        st.session_state["widget_generation"] = 0
    state.init_state()

    st.title(constants.APP_TITLE)
    st.caption(constants.APP_SUBTITLE)
    status_slot = st.container()
    _render_sidebar(status_slot)

    run_state = st.session_state["run_state"]
    result = state.get_result()

    if run_state == state.IDLE:
        st.info(constants.IDLE_NOTICE)
        return
    if run_state == state.RUNNING:
        st.info(constants.RUNNING_NOTICE)
        return
    if run_state == state.FAILED:
        st.error(st.session_state.get("run_error") or constants.NO_RESULT_NOTICE)
        st.info(constants.FAILED_ACTION)
        if result is None:  # keep earlier artefacts visible after a failed retry
            return
    elif result is None:
        st.error(constants.NO_RESULT_NOTICE)
        return

    _render_result(result)
