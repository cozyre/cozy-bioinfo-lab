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

from tp53 import config, constants
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
    BlastHit,
    PipelineInputs,
    PipelineResult,
    ValidationError,
)
from tp53.msa import alignment_statistics
from tp53.pipeline import run_pipeline
from tp53.viz import (
    alignment_html,
    blast_overview,
    conservation_profile,
    feature_bars,
    identity_heatmap,
)

from . import state

# --- UI strings (A1). Reuse constants.ERR_* for validation messages. -----------

TAB_OVERVIEW = "Overview"
TAB_BLAST = "BLAST"
TAB_ALIGNMENT = "Alignment"
TAB_CONSERVATION = "Conservation"
TAB_DOWNLOADS = "Downloads"
TABS: tuple[str, ...] = (
    TAB_OVERVIEW,
    TAB_BLAST,
    TAB_ALIGNMENT,
    TAB_CONSERVATION,
    TAB_DOWNLOADS,
)

QUERY_MODE_ACCESSION = "Accession"
QUERY_MODE_FASTA = "FASTA sequence"
QUERY_MODES: tuple[str, ...] = (QUERY_MODE_ACCESSION, QUERY_MODE_FASTA)

SOURCE_REMOTE = "Remote services"
SOURCE_UPLOAD = "Uploaded files"
SOURCE_MODES: tuple[str, ...] = (SOURCE_REMOTE, SOURCE_UPLOAD)

X_AXIS_REFERENCE = "Reference position"
X_AXIS_COLUMN = "Alignment column"
X_AXIS_CHOICES: tuple[str, ...] = (X_AXIS_REFERENCE, X_AXIS_COLUMN)

ERR_EMAIL_REQUIRED = "Enter a contact email before running remote services."
ERR_EMAIL_FORMAT = "Enter a valid contact email address."
ERR_ACCESSION_REQUIRED = "Enter an accession or switch to the FASTA query mode."
ERR_ACCESSION_FORMAT = "Enter a valid UniProt or NCBI protein accession."
ERR_FASTA_REQUIRED = "Paste a FASTA sequence or upload a .fa file."
ERR_ORGANISM_REQUIRED = "Enter the query organism for a FASTA query."
ERR_SPECIES_REQUIRED = "Select at least one target species."
ERR_XML_UPLOAD_REQUIRED = "Upload one BLAST XML file for every selected species."
ERR_UPLOAD_TOO_LARGE = "File '{name}' exceeds the {limit} upload limit."
ERR_UPLOAD_NOT_TEXT = "File '{name}' is not UTF-8 text."
ERR_XML_PARSE = "Could not parse BLAST XML for {species}: {detail}"
FEATURE_SKIP_NOTICE = (
    "Feature overlays and the feature summary are skipped: a bare FASTA query "
    "has no UniProt feature annotations."
)
IDLE_NOTICE = (
    "Select an accession or FASTA query in the sidebar, then choose "
    "Run analysis. Load example works fully offline."
)
FAILED_ACTION = (
    "Try uploading BLAST XML and/or MSA FASTA in the sidebar, or load the "
    "bundled example snapshot."
)
IDENTITY_NOTE = (
    "Identity counts gaps in the denominator: a column with a gap cannot "
    "reach 100%."
)

# --- Layout and parsing constants (UI only, not domain defaults). -------------

BYTES_PER_MIB = 1024 * 1024
ALIGNMENT_VIEWER_HEIGHT = 620
MAX_SMOOTHING_MULTIPLIER = 4
EMAIL_RE = re.compile(r"^[^@\s]+@[^@\s]+\.[^@\s]+$")
UNIPROT_ACCESSION_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})"
    r"(?:\.\d+)?$",
    re.IGNORECASE,
)
NCBI_ACCESSION_RE = re.compile(
    r"^(?:NP|XP|WP|YP|ZP)_\d+(?:\.\d+)?$", re.IGNORECASE
)
#: NCBI progress details: ``submitted: Mus musculus (ABC123)`` (§12).
RID_DETAIL_RE = re.compile(
    r"^(?:submitted|resuming):\s*(?P<species>.+?)\s*\((?P<rid>[A-Za-z0-9_-]+)\)\s*$"
)
#: Progress stage -> step number of 7 (§10).
STAGE_STEP: dict[str, int] = {
    "validate": 1,
    "retrieve": 2,
    "blast": 3,
    "select": 4,
    "fetch_sequences": 5,
    "msa": 6,
    "conservation": 7,
    "done": 7,
}
TOTAL_STEPS = 7


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
        return f"{byte_count // BYTES_PER_MIB} MB"
    return f"{byte_count} bytes"


def _read_upload(uploaded: Any, label: str) -> str | None:
    """Decode one uploaded file, enforcing :data:`config.MAX_UPLOAD_BYTES`."""
    if uploaded is None:
        return None
    size = getattr(uploaded, "size", None)
    if isinstance(size, int) and size > config.MAX_UPLOAD_BYTES:
        st.sidebar.error(
            ERR_UPLOAD_TOO_LARGE.format(
                name=getattr(uploaded, "name", label),
                limit=_format_limit(config.MAX_UPLOAD_BYTES),
            )
        )
        return None
    try:
        return uploaded.getvalue().decode("utf-8")
    except UnicodeDecodeError:
        st.sidebar.error(ERR_UPLOAD_NOT_TEXT.format(name=getattr(uploaded, "name", label)))
        return None


def _validate_accession(value: str) -> tuple[str, str] | None:
    """Return ``(accession, source)`` or ``None`` after an inline error."""
    accession = value.strip()
    if not accession:
        st.sidebar.error(ERR_ACCESSION_REQUIRED)
        return None
    if UNIPROT_ACCESSION_RE.fullmatch(accession):
        return accession, "uniprot"
    if NCBI_ACCESSION_RE.fullmatch(accession):
        return accession, "ncbi"
    st.sidebar.error(ERR_ACCESSION_FORMAT)
    return None


def _validate_fasta_query(text: str) -> ValidationError | None:
    """Validate a pasted/uploaded query FASTA with the shared domain rules."""
    try:
        record = parse_single_fasta(text)
        validate_protein(record)
    except ValidationError as exc:
        return exc
    return None


def _make_progress_callback(status: Any) -> Callable[[str, str], None]:
    """Build a progress callback that updates the step indicator and RIDs."""

    def callback(stage: str, detail: str) -> None:
        state.add_progress(stage, detail)
        match = RID_DETAIL_RE.match(detail)
        if stage == "blast" and match:
            state.remember_rid(match.group("species"), match.group("rid"))
        step = STAGE_STEP.get(stage, 1)
        status.update(
            label=f"Step {step} of {TOTAL_STEPS}: {detail}",
            state="running",
            expanded=True,
        )
        status.write(f"{stage}: {detail}")

    return callback


def _execute(inputs: PipelineInputs, status_slot: DeltaGenerator) -> None:
    """Run the pipeline once, updating run state and progress in the UI."""
    state.set_running()
    state.save_inputs(inputs)
    try:
        with status_slot.status(
            f"Step 1 of {TOTAL_STEPS}: validating input", expanded=True
        ) as status:
            clients = state.build_clients(
                inputs.email, inputs.api_key, blast_params=inputs.blast_params
            )
            result = run_pipeline(inputs, clients, _make_progress_callback(status))
            status.update(label="Analysis complete", state="complete", expanded=False)
        warnings = list(result.warnings or [])
        state.set_result(result, partial_success=bool(warnings))
        if warnings:
            st.warning("Analysis finished with warnings.")
        else:
            st.success("Analysis complete.")
    except Exception as exc:  # surface every remote or parsing failure in the UI
        state.set_failed(str(exc))
        st.error(f"Analysis failed: {exc}")


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
) -> PipelineInputs | None:
    """Validate sidebar values and assemble :class:`PipelineInputs`."""
    if not EMAIL_RE.fullmatch(email.strip()):
        st.sidebar.error(ERR_EMAIL_REQUIRED if not email.strip() else ERR_EMAIL_FORMAT)
        return None
    if not species:
        st.sidebar.error(ERR_SPECIES_REQUIRED)
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
            st.sidebar.error(ERR_FASTA_REQUIRED)
            return None
        problem = _validate_fasta_query(fasta_text)
        if problem is not None:
            st.sidebar.error(str(problem))
            return None
        if not organism:
            st.sidebar.error(ERR_ORGANISM_REQUIRED)
            return None
        query_fasta = fasta_text

    blast_xml_by_species: dict[str, str] | None = None
    if source_mode == SOURCE_UPLOAD and use_uploaded_xml:
        blast_xml_by_species = {}
        for target in species:
            text = _read_upload(xml_uploads.get(target), f"BLAST XML — {target}")
            if text is None:
                return None
            try:
                parse_blast_xml(text, target)
            except Exception as exc:
                st.sidebar.error(ERR_XML_PARSE.format(species=target, detail=exc))
                return None
            blast_xml_by_species[target] = text
        if not blast_xml_by_species:
            st.sidebar.error(ERR_XML_UPLOAD_REQUIRED)
            return None

    msa_fasta: str | None = None
    if source_mode == SOURCE_UPLOAD and use_uploaded_msa:
        msa_fasta = _read_upload(msa_upload, "MSA FASTA")
        if msa_fasta is None:
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
    )


def _render_sidebar(status_slot: DeltaGenerator) -> None:
    """Render inputs, advanced settings and the three action buttons."""
    st.sidebar.header("Query")
    query_mode = st.sidebar.radio(
        "Query input",
        list(QUERY_MODES),
        key=_widget_key("query_mode"),
        horizontal=True,
    )
    st.session_state["query_mode"] = query_mode

    accession = st.sidebar.text_input(
        "Protein accession",
        value=config.DEFAULT_QUERY_ACCESSION,
        key=_widget_key("accession"),
        disabled=query_mode != QUERY_MODE_ACCESSION,
        help="UniProt (P04637) or NCBI protein (NP_035770.2) accession.",
    )
    if query_mode == QUERY_MODE_FASTA:
        fasta_text = st.sidebar.text_area(
            "Query FASTA", key=_widget_key("fasta_text"), height=110
        )
        fasta_upload = st.sidebar.file_uploader(
            "Or upload query FASTA",
            type=["fa", "fasta", "faa", "txt"],
            key=_widget_key("fasta_upload"),
        )
        if fasta_upload is not None:
            uploaded_text = _read_upload(fasta_upload, "Query FASTA")
            if uploaded_text:
                fasta_text = uploaded_text
        query_organism = st.sidebar.text_input(
            "Query organism", key=_widget_key("query_organism")
        )
    else:
        fasta_text = ""
        query_organism = st.sidebar.text_input(
            "Query organism",
            value="Derived from accession",
            key=_widget_key("query_organism_accession"),
            disabled=True,
        )

    st.sidebar.header("Contact")
    email = st.sidebar.text_input(
        "Contact email",
        value=_secret("NCBI_EMAIL"),
        key=_widget_key("email"),
        help="Required by NCBI and EBI; never stored in the repository.",
    )
    api_key = _secret("NCBI_API_KEY") or None

    st.sidebar.header("Targets and source")
    species = st.sidebar.multiselect(
        "Target species",
        list(config.SPECIES),
        default=list(config.SPECIES),
        key=_widget_key("species"),
    )
    source_mode = st.sidebar.radio(
        "Data source",
        list(SOURCE_MODES),
        key=_widget_key("source_mode"),
    )
    use_uploaded_xml = False
    use_uploaded_msa = False
    xml_uploads: dict[str, Any] = {}
    msa_upload: Any = None
    if source_mode == SOURCE_UPLOAD:
        use_uploaded_xml = st.sidebar.checkbox(
            "Use uploaded BLAST XML (one file per species)",
            value=True,
            key=_widget_key("use_xml"),
        )
        if use_uploaded_xml:
            for target in species:
                xml_uploads[target] = st.sidebar.file_uploader(
                    f"BLAST XML — {target}",
                    type=["xml"],
                    key=_widget_key(f"xml_{_slug(target)}"),
                )
        use_uploaded_msa = st.sidebar.checkbox(
            "Use uploaded MSA FASTA",
            value=True,
            key=_widget_key("use_msa"),
        )
        if use_uploaded_msa:
            msa_upload = st.sidebar.file_uploader(
                "MSA FASTA",
                type=["fa", "fasta", "faa", "txt"],
                key=_widget_key("msa_upload"),
            )

    with st.sidebar.expander("Advanced settings", expanded=False):
        blast_database = st.selectbox(
            "BLAST database",
            list(config.BLAST_DATABASES),
            index=0,
            key=_widget_key("blast_database"),
        )
        defaults = config.SelectionParams()
        cons_defaults = config.ConservationParams()
        evalue_cutoff = st.number_input(
            "E-value cutoff",
            min_value=0.0,
            value=float(defaults.evalue_cutoff),
            format="%.1e",
            key=_widget_key("evalue"),
        )
        min_coverage_pct = st.number_input(
            "Minimum query coverage (%)",
            min_value=0.0,
            max_value=100.0,
            value=defaults.min_query_coverage * 100.0,
            key=_widget_key("coverage"),
        )
        conservation_threshold = st.number_input(
            "Conservation threshold (%)",
            min_value=0.0,
            max_value=100.0,
            value=cons_defaults.region_identity_threshold,
            key=_widget_key("threshold"),
        )
        min_region_length = st.number_input(
            "Minimum region length",
            min_value=1,
            value=cons_defaults.min_region_length,
            key=_widget_key("min_region"),
        )
        st.number_input(
            "Smoothing window (Conservation tab)",
            min_value=1,
            value=config.SMOOTHING_WINDOW_DEFAULT,
            key=_widget_key("smoothing_advanced"),
        )

    run_col, example_col, reset_col = st.sidebar.columns(3, gap="small")
    if run_col.button("Run analysis", key=_widget_key("run"), type="primary"):
        inputs = _build_inputs(
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
        )
        if inputs is not None:
            _execute(inputs, status_slot)
    if example_col.button("Load example", key=_widget_key("example")):
        _load_example(status_slot)
    if reset_col.button("Reset", key=_widget_key("reset")):
        state.reset_state()
        st.session_state["widget_generation"] += 1
        st.rerun()


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
    st.markdown("Classes: " + " · ".join(parts), unsafe_allow_html=True)


def _render_overview(result: PipelineResult) -> None:
    """Render status, run metadata and FR-23 summary metrics."""
    with st.status("Analysis complete", state="complete", expanded=False) as status:
        events = st.session_state.get("progress", [])
        if events:
            for stage, detail in events:
                status.write(f"{stage}: {detail}")
        else:
            status.write("All pipeline stages completed.")

    if result.run_metadata:
        st.subheader("Run metadata")
        st.json(_json_safe(result.run_metadata))

    summary = result.conservation.summary if result.conservation else {}
    mean = summary.get(SUMMARY_MEAN_IDENTITY, 0.0)
    count = summary.get(SUMMARY_FULLY_CONSERVED_COUNT, 0)
    share = summary.get(SUMMARY_FULLY_CONSERVED_PCT, 0.0)
    longest = summary.get(SUMMARY_LONGEST_REGION, 0)
    regions = summary.get(SUMMARY_REGION_COUNT, 0)
    columns = summary.get(SUMMARY_ALIGNMENT_LENGTH, 0)

    st.subheader("Conservation summary")
    first, second, third = st.columns(3, gap="small")
    with first:
        st.metric("Mean identity", f"{float(mean):.1f}%")
    with second:
        st.metric("Fully conserved positions", str(count))
    with third:
        st.metric("Fully conserved share", f"{float(share):.1f}%")
    first, second, third = st.columns(3, gap="small")
    with first:
        st.metric("Longest conserved region", f"{longest} columns")
    with second:
        st.metric("Conserved regions", str(regions))
    with third:
        st.metric("Alignment columns", str(columns))
    st.caption(IDENTITY_NOTE)


def _style_candidates(frame: pd.DataFrame, accession: str) -> Any:
    """Highlight the selected candidate row without relying on colour alone."""

    def highlight(row: pd.Series) -> list[str]:
        selected = str(row.get("accession", "")) == accession
        return [f"background-color: {constants.COLOR_PALETTE['sky_blue']}; font-weight: 600"] * len(row) if selected else [""] * len(row)

    return frame.style.apply(highlight, axis=1)


def _best_hsp(hit: BlastHit) -> Any:
    """Best HSP of a hit, or ``None`` if the hit has no HSPs."""
    return hit.hsps[0] if hit.hsps else None


def _render_blast(result: PipelineResult) -> None:
    """Render candidate tables, the overview chart and raw hit alignments."""
    if not result.selected:
        st.info("No ortholog hits were selected. Try uploading BLAST XML or load the example.")
    else:
        figure = blast_overview(result.selected)
        st.plotly_chart(figure, width="stretch", key="blast_overview")

    for selected in result.selected:
        with st.container(border=True):
            st.subheader(f"{selected.species} — {selected.hit.accession}")
            st.caption(selected.reason)
            rows: list[dict[str, Any]] = []
            for hit in selected.candidates or [selected.hit]:
                hsp = _best_hsp(hit)
                rows.append(
                    {
                        "accession": hit.accession,
                        "title": hit.title,
                        "organism": hit.organism,
                        "bit_score": round(hsp.bit_score, 1) if hsp else 0.0,
                        "evalue": hsp.evalue if hsp else 0.0,
                        "identity_pct": round(hsp.percent_identity, 1) if hsp else 0.0,
                        "coverage_pct": round(selected.query_coverage * 100.0, 1),
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
                f"Raw hit alignment — {selected.species} / {selected.hit.accession}",
                expanded=False,
            ):
                for index, hsp in enumerate(selected.hit.hsps, start=1):
                    st.caption(f"HSP {index}: {hsp.query_start}-{hsp.query_end}")
                    st.code(
                        f"Query  {hsp.qseq}\n       {hsp.midline}\nSubject {hsp.hseq}",
                        language=None,
                    )


def _render_alignment(result: PipelineResult) -> None:
    """Render the block-wise alignment viewer and per-sequence statistics."""
    if result.alignment is None:
        st.info("No alignment is available for this run.")
        return
    alignment = result.alignment
    identity = result.conservation.identity if result.conservation else []
    st.dataframe(
        alignment_statistics(alignment), width="stretch", hide_index=True
    )
    _class_legend()
    html = alignment_html(alignment, identity, block=60)
    st.iframe(html, height=ALIGNMENT_VIEWER_HEIGHT)


def _render_regions(result: PipelineResult) -> None:
    """Render the conserved-region table (§9.6, FR-20)."""
    regions = result.conservation.regions if result.conservation else []
    if not regions:
        st.info("No conserved regions meet the current threshold.")
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
        st.info("No conservation results are available for this run.")
        return
    _class_legend()
    window_max = max(config.SMOOTHING_WINDOW_DEFAULT * MAX_SMOOTHING_MULTIPLIER, 1)
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
        "Smoothing window",
        min_value=1,
        max_value=window_max,
        key=smoothing_key,
    )
    axis_choice = st.radio(
        "X axis",
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

    st.subheader("Conserved regions")
    _render_regions(result)

    st.subheader("Feature comparison")
    features = result.query.features if result.query else []
    if features:
        st.dataframe(
            result.conservation.feature_summary,
            width="stretch",
            hide_index=True,
        )
        st.plotly_chart(feature_bars(result), width="stretch", key="feature_bars")
    else:
        st.info(FEATURE_SKIP_NOTICE)

    st.subheader("Pairwise identity")
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
            f"BLAST XML — {species}",
            data=xml_text.encode("utf-8"),
            file_name=f"blast_{_slug(species)}.xml",
            mime="application/xml",
            key=_widget_key(f"download_blast_{_slug(species)}"),
        )
    if result.raw_msa_fasta:
        any_download = True
        st.download_button(
            "MSA FASTA",
            data=result.raw_msa_fasta.encode("utf-8"),
            file_name="tp53_alignment.fasta",
            mime="text/plain",
            key=_widget_key("download_msa"),
        )
    csv_text = _identity_csv(result)
    if csv_text is not None:
        any_download = True
        st.download_button(
            "Identity CSV",
            data=csv_text.encode("utf-8"),
            file_name="tp53_identity.csv",
            mime="text/csv",
            key=_widget_key("download_identity"),
        )
    if not any_download:
        st.info("No raw artefacts are available for this run.")


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
        st.info(IDLE_NOTICE)
        return
    if run_state == state.RUNNING:
        st.info("Analysis running. Step progress is shown above.")
        return
    if run_state == state.FAILED:
        st.error(st.session_state.get("run_error", "Analysis failed."))
        st.info(FAILED_ACTION)
        return

    for warning in result.warnings if result is not None else []:
        st.warning(warning)
    if result is None:
        st.error("No result is available. Run the analysis or load the example.")
        return

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
