"""Session state, cached example inputs and real client wiring (FR-31, FR-33).

Only this module and :mod:`app.pages` import Streamlit. The pipeline receives
real clients with an injected :class:`requests.Session`; no HTTP object is
created by the framework-free core.
"""

from __future__ import annotations

import json
from functools import partial
from pathlib import Path
from typing import Any

import requests
import streamlit as st

from tp53 import config
from tp53.clients import ebi, ncbi, uniprot
from tp53.fasta_io import parse_fasta
from tp53.models import (
    PipelineClients,
    PipelineInputs,
    PipelineResult,
)

#: Keys managed by :func:`init_state` and cleared by :func:`reset_state`.
STATE_KEYS: tuple[str, ...] = (
    "run_state",
    "result",
    "run_error",
    "progress",
    "rids",
    "inputs",
    "query_mode",
    "query_signature",
)

#: Valid values of ``st.session_state["run_state"]`` (§10).
IDLE = "idle"
RUNNING = "running"
SUCCESS = "success"
PARTIAL = "partial"
FAILED = "failed"


def init_state() -> None:
    """Populate missing session-state keys with their defaults."""
    defaults: dict[str, Any] = {
        "run_state": IDLE,
        "result": None,
        "run_error": "",
        "progress": [],
        "rids": {},
        "inputs": None,
        "query_mode": "accession",
        "query_signature": None,
    }
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value


def reset_state() -> None:
    """Forget results, progress and saved RIDs, returning the app to idle."""
    for key in STATE_KEYS:
        st.session_state.pop(key, None)
    init_state()


def get_result() -> PipelineResult | None:
    """Return the single :class:`PipelineResult` rendered by the UI."""
    init_state()
    result = st.session_state["result"]
    return result if isinstance(result, PipelineResult) else None


def set_running() -> None:
    """Mark a run as in progress and clear the previous error."""
    init_state()
    st.session_state["run_state"] = RUNNING
    st.session_state["run_error"] = ""
    st.session_state["progress"] = []


def set_result(result: PipelineResult, *, partial_success: bool = False) -> None:
    """Store a completed result and its success state."""
    init_state()
    st.session_state["result"] = result
    st.session_state["run_state"] = PARTIAL if partial_success else SUCCESS
    st.session_state["run_error"] = ""


def set_failed(message: str) -> None:
    """Record a run failure while keeping the previous result available."""
    init_state()
    st.session_state["run_state"] = FAILED
    st.session_state["run_error"] = message


def add_progress(stage: str, detail: str) -> None:
    """Append one ``(stage, detail)`` progress event to the session log."""
    init_state()
    st.session_state["progress"].append((stage, detail))


def remember_rid(species: str, rid: str) -> None:
    """Save a BLAST request identifier so a refresh can resume polling (§12)."""
    init_state()
    st.session_state["rids"][species] = rid


def clear_rids() -> None:
    """Discard saved request identifiers, e.g. after Reset or a new query."""
    init_state()
    st.session_state["rids"] = {}


def save_inputs(inputs: PipelineInputs) -> None:
    """Remember the inputs of the last run for the metadata panel."""
    init_state()
    st.session_state["inputs"] = inputs


@st.cache_resource(show_spinner=False)
def _http_session() -> requests.Session:
    """Return the process-wide HTTP session shared by all real clients."""
    return requests.Session()


def _poll_blast_for_species(
    rid: str,
    species: str,
    *,
    session: requests.Session,
    email: str,
    api_key: str | None,
) -> str:
    """Adapt ``(rid, species)`` to the blocking NCBI polling client (§7.2).

    ``species`` is part of the injected client contract and is used only for
    progress/error context; the RID alone identifies the remote BLAST job.
    """
    del species  # the RID is sufficient; the name keeps the contract explicit
    return ncbi.poll_blast(
        rid,
        session,
        config.BLAST_MAX_WAIT_S,
        email=email,
        api_key=api_key,
    )


def build_clients(
    email: str,
    api_key: str | None,
    *,
    blast_params: config.BlastParams | None = None,
) -> PipelineClients:
    """Bind real remote clients to one injected HTTP session.

    ``functools.partial`` supplies session, contact email, API key and BLAST
    parameters so every attribute matches the zero/extra-argument signatures
    in :class:`~tp53.models.PipelineClients`.

    Args:
        email: user-supplied contact email sent to NCBI and EBI.
        api_key: optional NCBI API key from ``st.secrets``.
        blast_params: BLAST database/E-value/hit-list settings; defaults to
            :class:`~tp53.config.BlastParams`.

    Returns:
        Client bundle ready to pass to ``run_pipeline``.
    """
    session = _http_session()
    params = blast_params or config.BlastParams()
    return PipelineClients(
        fetch_uniprot=partial(uniprot.fetch_query_info, session=session),
        fetch_ncbi_query_info=partial(
            ncbi.fetch_ncbi_query_info,
            session=session,
            email=email,
            api_key=api_key,
        ),
        fetch_protein=partial(
            ncbi.fetch_protein,
            session=session,
            email=email,
            api_key=api_key,
        ),
        submit_blastp_staggered=partial(
            ncbi.submit_blastp_staggered,
            params=params,
            email=email,
            session=session,
            api_key=api_key,
        ),
        poll_blast=partial(
            _poll_blast_for_species,
            session=session,
            email=email,
            api_key=api_key,
        ),
        run_clustalo=partial(
            ebi.run_clustalo,
            email=email,
            session=session,
            timeout_s=config.EBI_MAX_WAIT_S,
        ),
    )


@st.cache_data(ttl=config.CACHE_TTL_S, show_spinner=False)
def _read_example_file(path: str) -> str:
    """Read one bundled snapshot file as UTF-8 text (cached, never written)."""
    return Path(path).read_text(encoding="utf-8")


@st.cache_data(ttl=config.CACHE_TTL_S, show_spinner=False)
def load_example_inputs() -> PipelineInputs:
    """Build offline pipeline inputs from the bundled TP53 snapshot (A8).

    The UniProt JSON is parsed with the same pure function as the live client;
    per-accession FASTA files become the ``full_sequences`` lookup that
    replaces Entrez fetches. Fixture BLAST XML and MSA FASTA stand in for the
    remote stages, so this mode never calls a service.

    Returns:
        Fully offline :class:`~tp53.models.PipelineInputs`.
    """
    directory = Path(__file__).parent / config.EXAMPLE_DATA_DIRNAME
    payload = json.loads(
        _read_example_file(str(directory / config.EXAMPLE_UNIPROT_JSON_FILENAME))
    )
    query_info = uniprot.parse_query_info(payload, config.DEFAULT_QUERY_ACCESSION)

    blast_xml = _read_example_file(str(directory / config.EXAMPLE_BLAST_XML_FILENAME))
    msa_fasta = _read_example_file(str(directory / config.EXAMPLE_MSA_FASTA_FILENAME))
    full_sequences: dict[str, str] = {}
    for filename in config.EXAMPLE_SEQUENCE_FILENAMES:
        records = parse_fasta(_read_example_file(str(directory / filename)))
        if records:
            full_sequences[records[0].accession] = records[0].sequence

    return PipelineInputs(
        query_accession=config.DEFAULT_QUERY_ACCESSION,
        accession_source="uniprot",
        query_info=query_info,
        species=config.SPECIES,
        blast_xml_by_species={species: blast_xml for species in config.SPECIES},
        msa_fasta=msa_fasta,
        full_sequences=full_sequences,
    )
