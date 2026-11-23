"""Offline orchestration tests for the phase-3 pipeline (§6.1, §12, AC-1)."""

from __future__ import annotations

import json
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
from xml.etree import ElementTree

import pytest

from tp53 import config, constants
from tp53.clients.base import PollOutcome, RemoteServiceError
from tp53.clients.uniprot import parse_query_info
from tp53.fasta_io import parse_fasta
from tp53.models import (
    PipelineClients,
    PipelineInputs,
    PipelineResult,
    ProgressCallback,
    QueryInfo,
    ValidationError,
)
from tp53 import pipeline as pipeline_module
from tp53.pipeline import STAGES, run_pipeline

DATA_DIR = Path(__file__).parent / "data"
FULL_XML = (DATA_DIR / config.EXAMPLE_BLAST_XML_FILENAME).read_text()
MSA_FASTA = (DATA_DIR / config.EXAMPLE_MSA_FASTA_FILENAME).read_text()
QUERY_INFO = parse_query_info(
    json.loads((DATA_DIR / config.EXAMPLE_UNIPROT_JSON_FILENAME).read_text()),
    config.DEFAULT_QUERY_ACCESSION,
)

#: Expected deterministic best hit for each species in the fixture (§9.2).
EXPECTED_ACCESSION = {
    "Homo sapiens": config.DEFAULT_QUERY_ACCESSION,
    "Pan troglodytes": "XP_001172077.2",
    "Mus musculus": "NP_035770.2",
    "Rattus norvegicus": "NP_112251.2",
    "Bos taurus": "NP_776626.1",
}
#: IDs used by the recorded Clustal Omega MSA, in input order.
MSA_LABEL = {
    "Homo sapiens": "Homo",
    "Pan troglodytes": "Pan",
    "Mus musculus": "Mus",
    "Rattus norvegicus": "Rattus",
    "Bos taurus": "Bos",
}


def xml_for_species(species: str, *, empty: bool = False) -> str:
    """Derive a per-species BLAST XML from the committed five-species fixture."""
    root = ElementTree.fromstring(FULL_XML)
    hits = root.find("./BlastOutput_iterations/Iteration/Iteration_hits")
    assert hits is not None
    for hit in list(hits):
        definition = hit.findtext("Hit_def", "")
        if empty or f"[{species}]" not in definition:
            hits.remove(hit)
    return '<?xml version="1.0"?>\n' + ElementTree.tostring(root, encoding="unicode")


def aligned_fasta_from_input(fasta_text: str) -> str:
    """Re-key the recorded MSA rows onto the IDs of the current input set."""
    fixture_rows = {
        record.id: record.sequence for record in parse_fasta(MSA_FASTA)
    }
    chunks: list[str] = []
    for record in parse_fasta(fasta_text):
        header = f"{record.id} {record.description}".strip()
        species = header.split(config.HEADER_SEPARATOR)[0]
        label = MSA_LABEL[species]
        chunks.append(f">{header}\n{fixture_rows[label]}\n")
    return "".join(chunks)


@dataclass
class FakeClients:
    """Injected clients backed entirely by committed fixtures."""

    fail_poll: frozenset[str] = frozenset()
    empty_hits: frozenset[str] = frozenset()
    fail_submit: frozenset[str] = frozenset()
    xml_by_species: dict[str, str] = field(default_factory=dict)
    network_calls: list[str] = field(default_factory=list)
    submit_calls: list[list[tuple[str, str]]] = field(default_factory=list)
    poll_calls: list[tuple[str, str]] = field(default_factory=list)
    clustalo_calls: list[str] = field(default_factory=list)
    # Phase-4 split contract, enabled only by ``as_contract(split=True)``.
    msa_job_id: str = "EBI_JOB_1"
    submit_blastp_calls: list[tuple[str, str]] = field(default_factory=list)
    rid_by_species: dict[str, str] = field(default_factory=dict)
    species_by_rid: dict[str, str] = field(default_factory=dict)
    one_shot_polls: dict[str, int] = field(default_factory=dict)
    submit_clustalo_calls: list[str] = field(default_factory=list)
    resume_clustalo_calls: list[str] = field(default_factory=list)
    resume_fasta: str = ""

    def as_contract(self, *, split: bool = False) -> PipelineClients:
        """Bind the fakes to the public injected-clients contract.

        The split clients expose per-job submission, one-shot polling and
        resumable alignment; existing tests keep the original blocking path.
        """
        contract: dict[str, Any] = {
            "fetch_uniprot": self.fetch_uniprot,
            "fetch_ncbi_query_info": self.fetch_ncbi_query_info,
            "fetch_protein": self.fetch_protein,
            "submit_blastp_staggered": self.submit_blastp_staggered,
            "poll_blast": self.poll_blast,
            "run_clustalo": self.run_clustalo,
        }
        if split:
            contract.update(
                submit_blastp=self.submit_blastp,
                poll_blast_once=self.poll_blast_once,
                submit_clustalo=self.submit_clustalo,
                resume_clustalo=self.resume_clustalo,
            )
        return PipelineClients(**contract)

    def fetch_uniprot(self, accession: str) -> QueryInfo:
        self.network_calls.append(f"uniprot:{accession}")
        return QUERY_INFO

    def fetch_ncbi_query_info(self, accession: str) -> QueryInfo:
        self.network_calls.append(f"ncbi-query:{accession}")
        return QUERY_INFO

    def fetch_protein(self, accession: str):
        self.network_calls.append(f"protein:{accession}")
        return parse_fasta((DATA_DIR / f"{accession}.fa").read_text())[0]

    def submit_blastp_staggered(self, jobs: list[tuple[str, str]]) -> dict[str, str]:
        self.submit_calls.append(list(jobs))
        for species, _sequence in jobs:
            if species in self.fail_submit:
                raise RemoteServiceError(f"submission failed for {species}")
        return {
            species: f"RID_{species.split()[0].upper()}{index}"
            for index, (species, _sequence) in enumerate(jobs, start=1)
        }

    def poll_blast(self, rid: str, species: str) -> str:
        self.poll_calls.append((rid, species))
        if species in self.fail_poll:
            raise RemoteServiceError(f"polling failed for {species}")
        if species in self.xml_by_species:
            return self.xml_by_species[species]
        return xml_for_species(species, empty=species in self.empty_hits)

    def run_clustalo(self, fasta: str) -> str:
        self.clustalo_calls.append(fasta)
        return aligned_fasta_from_input(fasta)

    def submit_blastp(self, sequence: str, species: str) -> str:
        """Return a stable RID and record the submission in real time."""
        self.submit_blastp_calls.append((species, sequence))
        rid = f"SPLIT_{species.split()[0].upper()}"
        self.rid_by_species[species] = rid
        self.species_by_rid[rid] = species
        return rid

    def poll_blast_once(self, rid: str) -> PollOutcome:
        """Return pending once, then the fixture XML for that species."""
        self.one_shot_polls[rid] = self.one_shot_polls.get(rid, 0) + 1
        species = self.species_by_rid[rid]
        if self.one_shot_polls[rid] < 2:
            return PollOutcome.pending()
        return PollOutcome.complete(xml_for_species(species))

    def submit_clustalo(self, fasta: str) -> str:
        """Save the submitted FASTA so a resumed job can reuse it."""
        self.submit_clustalo_calls.append(fasta)
        return self.msa_job_id

    def resume_clustalo(self, job_id: str) -> str:
        """Return the alignment previously produced for the saved job."""
        self.resume_clustalo_calls.append(job_id)
        if self.resume_fasta:
            return self.resume_fasta
        return aligned_fasta_from_input(self.submit_clustalo_calls[-1])


def make_inputs(**overrides: object) -> PipelineInputs:
    """Default remote-mode inputs, overridable per test."""
    values: dict[str, object] = {
        "query_accession": config.DEFAULT_QUERY_ACCESSION,
        "species": config.SPECIES,
        "email": "test@example.org",
    }
    values.update(overrides)
    return PipelineInputs(**values)  # type: ignore[arg-type]


def run(
    inputs: PipelineInputs,
    clients: FakeClients,
    *,
    split: bool = False,
    sleep: Callable[[float], None] | None = None,
    clock: Callable[[], float] | None = None,
) -> tuple[PipelineResult, list[tuple[str, str]]]:
    """Run the pipeline while collecting stable progress events."""
    events: list[tuple[str, str]] = []

    def progress(stage: str, detail: str) -> None:
        events.append((stage, detail))

    callback: ProgressCallback = progress
    kwargs: dict[str, Any] = {}
    if sleep is not None:
        kwargs["sleep"] = sleep
    if clock is not None:
        kwargs["clock"] = clock
    result = run_pipeline(inputs, clients.as_contract(split=split), callback, **kwargs)
    return result, events


def fake_timers() -> tuple[list[float], Callable[[float], None], Callable[[], float]]:
    """Injected sleep/clock pair that advances without real waiting."""
    waits: list[float] = []
    now = 0.0

    def sleep(seconds: float) -> None:
        nonlocal now
        waits.append(seconds)
        now += seconds

    return waits, sleep, lambda: now


def test_success_with_five_species( ) -> None:
    """AC-1: one labelled ortholog per species from the fixture."""
    clients = FakeClients()
    result, _events = run(make_inputs(), clients)

    assert [hit.hit.accession for hit in result.selected] == [
        EXPECTED_ACCESSION[species] for species in config.SPECIES
    ]
    assert len(result.input_records) == 5
    assert result.input_records[0].id == "Homo sapiens|P04637"
    assert all(
        record.id == f"{species}{config.HEADER_SEPARATOR}{EXPECTED_ACCESSION[species]}"
        for species, record in zip(
            config.SPECIES[1:], result.input_records[1:], strict=True
        )
    )
    assert result.alignment is not None
    assert len(result.alignment.rows) == 5
    assert len({len(row) for row in result.alignment.rows}) == 1
    assert result.conservation is not None
    assert len(result.conservation.identity) == result.alignment.length
    # Only the benign ID-remapping notice is allowed (Clustal truncates headers
    # containing spaces); no stage may have failed.
    assert all("failed" not in warning.lower() for warning in result.warnings)
    assert len(clients.poll_calls) == 4  # the query organism is not BLASTed
    assert len(clients.clustalo_calls) == 1


def test_partial_success_after_per_species_blast_failure() -> None:
    """A failed species is skipped with a warning when others remain."""
    clients = FakeClients(fail_poll=frozenset({"Mus musculus"}))
    result, _events = run(make_inputs(), clients)

    assert len(result.selected) == 4
    assert "Mus musculus" not in result.raw_blast_xml
    assert any(
        "BLAST failed for Mus musculus" in warning for warning in result.warnings
    )
    assert result.alignment is not None
    assert len(result.alignment.rows) == 4


def test_no_hit_species_is_skipped() -> None:
    """A BLAST result without qualifying hits produces a warning, not a crash."""
    clients = FakeClients(empty_hits=frozenset({"Bos taurus"}))
    result, events = run(make_inputs(), clients)

    assert len(result.selected) == 4
    assert any("No qualifying BLAST hit for Bos taurus" in w for w in result.warnings)
    assert ("select", "no hit: Bos taurus") in events
    assert result.conservation is not None


def test_no_ortholog_failure_is_explicit() -> None:
    """Without any usable ortholog the run fails with a clear message."""
    clients = FakeClients(empty_hits=frozenset(config.SPECIES))
    with pytest.raises(ValidationError, match="ortholog"):
        run(make_inputs(), clients)


def test_example_mode_uses_supplied_artefacts_without_network() -> None:
    """Uploaded/snapshot inputs bypass every client call (FR-33, A8)."""
    xml = {species: xml_for_species(species) for species in config.SPECIES}
    sequences = {
        accession: (DATA_DIR / f"{accession}.fa").read_text()
        for accession in EXPECTED_ACCESSION.values()
        if (DATA_DIR / f"{accession}.fa").exists()
    }
    clients = FakeClients(xml_by_species=xml)
    inputs = make_inputs(
        query_accession=None,
        query_info=QUERY_INFO,
        blast_xml_by_species=xml,
        msa_fasta=MSA_FASTA,
        full_sequences=sequences,
    )
    result, _events = run(inputs, clients)

    assert len(result.selected) == 5
    assert result.raw_msa_fasta == MSA_FASTA
    assert result.raw_blast_xml == xml
    assert clients.network_calls == []
    assert clients.submit_calls == []
    assert clients.poll_calls == []
    assert clients.clustalo_calls == []
    assert result.run_metadata["mode"] == "example"


def test_resume_from_existing_rids() -> None:
    """Saved RIDs are polled instead of resubmitted (§12)."""
    rids = {
        species: f"SAVED_{species.split()[0]}" for species in config.SPECIES
    }
    clients = FakeClients()
    inputs = make_inputs(existing_rids=rids)
    result, events = run(inputs, clients)

    assert clients.submit_calls == []
    targets = [species for species in config.SPECIES if species != "Homo sapiens"]
    assert [rid for rid, _species in clients.poll_calls] == [
        rids[species] for species in targets
    ]
    assert result.blast_results and len(result.blast_results) == 4
    for species in targets:
        assert (
            "blast",
            constants.PROGRESS_BLAST_RESUME.format(
                species=species, rid=rids[species]
            ),
        ) in events


def test_progress_stages_and_submission_details() -> None:
    """Stages are stable and RID details allow the UI to save new jobs."""
    clients = FakeClients()
    _result, events = run(make_inputs(), clients)

    assert list(dict.fromkeys(stage for stage, _ in events)) == list(STAGES)
    details = [detail for _stage, detail in events]
    first_target = config.SPECIES[1]
    first_rid = f"RID_{first_target.split()[0].upper()}1"
    assert f"submitted: {first_target} ({first_rid})" in details
    assert f"polling: {first_target}" in details
    assert f"ready: {first_target}" in details
    assert details[-1].startswith("complete: 5 sequences")


def test_run_metadata_and_raw_artefacts() -> None:
    """Metadata records provenance, parameters, versions and raw artefacts."""
    clients = FakeClients()
    result, _events = run(make_inputs(), clients)
    metadata = result.run_metadata

    assert metadata["query"] == config.DEFAULT_QUERY_ACCESSION
    assert metadata["query_source"] == "uniprot"
    assert metadata["mode"] == "remote"
    assert metadata["species"] == list(config.SPECIES)
    assert metadata["blast_params"] == {
        "database": config.BlastParams().database,
        "expect": config.BlastParams().expect,
        "hitlist_size": config.BlastParams().hitlist_size,
        "program": config.BlastParams().program,
    }
    assert metadata["selection_params"]["min_query_coverage"] == pytest.approx(
        config.SelectionParams().min_query_coverage
    )
    assert metadata["conservation_params"]["region_identity_threshold"] == 100.0
    assert metadata["smoothing_window"] == config.SMOOTHING_WINDOW_DEFAULT
    assert metadata["selected_accessions"] == [
        EXPECTED_ACCESSION[species] for species in config.SPECIES
    ]
    assert set(metadata["library_versions"]) == {
        "biopython",
        "requests",
        "pandas",
        "numpy",
        "plotly",
    }
    assert all(metadata["library_versions"].values())
    assert metadata["raw_artifacts"]["blast_xml"] is True
    assert metadata["raw_artifacts"]["msa_fasta"] is True
    assert set(result.raw_blast_xml) == {
        species for species in config.SPECIES if species != "Homo sapiens"
    }
    assert result.raw_msa_fasta == aligned_fasta_from_input(clients.clustalo_calls[0])


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"query_accession": None}, "exactly one query source"),
        (
            {"query_accession": "P04637", "query_fasta": ">x\nACDEFGHIKLMNPQRSTVWY"},
            "exactly one query source",
        ),
        ({"query_accession": "   "}, "empty"),
        (
            {"query_accession": "P04637", "accession_source": "unknown"},
            "Unknown accession source",
        ),
    ],
)
def test_invalid_query_inputs_are_rejected(
    overrides: dict[str, object], message: str
) -> None:
    """Invalid query sources fail before any remote call (FR-03, §12)."""
    clients = FakeClients()
    with pytest.raises(ValidationError, match=message):
        run(make_inputs(**overrides), clients)
    assert clients.network_calls == []


def test_fasta_query_has_no_feature_annotations() -> None:
    """A bare FASTA query skips UniProt feature overlays (§9.7, §12)."""
    fasta = f">{config.DEFAULT_QUERY_ACCESSION}\n{QUERY_INFO.record.sequence}\n"
    clients = FakeClients()
    inputs = make_inputs(
        query_accession=None,
        query_fasta=fasta,
        query_organism="Homo sapiens",
        species=config.SPECIES,
    )
    result, _events = run(inputs, clients)

    assert result.query is not None
    assert result.query.features == []
    assert result.run_metadata["query_source"] == "fasta"
    assert result.conservation is not None
    assert list(result.conservation.feature_summary["feature"]) == [
        constants.OUTSIDE_FEATURES_LABEL
    ]


def test_split_clients_emit_rids_immediately_and_poll_round_robin() -> None:
    """§9.1/§12: RIDs arrive before polling and all jobs share one poll loop."""
    clients = FakeClients()
    waits, sleep, clock = fake_timers()
    species = config.SPECIES[:3]  # query plus two BLAST targets

    result, events = run(
        make_inputs(species=species, smoothing_window=7),
        clients,
        split=True,
        sleep=sleep,
        clock=clock,
    )

    details = [detail for _stage, detail in events]
    submitted = [detail for detail in details if detail.startswith("submitted:")]
    polling = [detail for detail in details if detail.startswith("polling:")]

    assert submitted[0] == constants.PROGRESS_BLAST_SUBMITTED.format(
        species=config.SPECIES[1], rid="SPLIT_PAN"
    )
    assert details.index(submitted[0]) < details.index(polling[0])
    assert polling == [
        constants.PROGRESS_BLAST_POLLING.format(species=config.SPECIES[1]),
        constants.PROGRESS_BLAST_POLLING.format(species=config.SPECIES[2]),
        constants.PROGRESS_BLAST_POLLING.format(species=config.SPECIES[1]),
        constants.PROGRESS_BLAST_POLLING.format(species=config.SPECIES[2]),
    ]
    assert waits == [
        config.BLAST_STAGGER_DELAY_S,
        config.BLAST_MIN_REQUEST_INTERVAL_S,
        config.BLAST_MIN_REQUEST_INTERVAL_S,
        config.BLAST_POLL_INTERVAL_S - config.BLAST_MIN_REQUEST_INTERVAL_S,
        config.BLAST_MIN_REQUEST_INTERVAL_S,
    ]
    assert clients.poll_calls == []
    assert clients.one_shot_polls == {"SPLIT_PAN": 2, "SPLIT_MUS": 2}
    assert result.run_metadata["smoothing_window"] == 7


def test_split_msa_job_is_submitted_then_resumed() -> None:
    """§12: an EBI job ID is emitted at submission and reused on a new run."""
    clients = FakeClients()
    waits, sleep, clock = fake_timers()
    species = config.SPECIES[:3]  # include two BLAST targets so stagger applies

    first, events = run(
        make_inputs(species=species, smoothing_window=3),
        clients,
        split=True,
        sleep=sleep,
        clock=clock,
    )
    details = [detail for _stage, detail in events]
    assert constants.PROGRESS_MSA_SUBMITTED.format(job_id=clients.msa_job_id) in details
    assert first.run_metadata["smoothing_window"] == 3
    assert waits[0] == config.BLAST_STAGGER_DELAY_S
    assert len(clients.submit_clustalo_calls) == 1

    resumed = FakeClients(
        resume_fasta=aligned_fasta_from_input(clients.submit_clustalo_calls[0])
    )
    second_inputs = make_inputs(
        species=species,
        msa_job_id=clients.msa_job_id,
        smoothing_window=3,
    )
    second, second_events = run(
        second_inputs, resumed, split=True, sleep=lambda _: None, clock=lambda: 0.0
    )
    assert resumed.submit_clustalo_calls == []
    assert resumed.resume_clustalo_calls == [clients.msa_job_id]
    assert constants.PROGRESS_MSA_RESUME.format(
        job_id=clients.msa_job_id
    ) in [detail for _stage, detail in second_events]
    assert second.run_metadata["smoothing_window"] == 3


def test_round_robin_timeout_warns_for_each_pending_species(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """An expired BLAST deadline warns per species instead of hanging (§12)."""
    monkeypatch.setattr(config, "BLAST_MAX_WAIT_S", 0)
    clients = FakeClients()
    warnings_out: list[str] = []

    raw_xml, parsed = pipeline_module._collect_blast(
        make_inputs(),
        clients.as_contract(split=True),
        QUERY_INFO,
        lambda stage, detail: None,
        warnings_out,
        sleep=lambda _: None,
        clock=lambda: 0.0,
    )

    assert parsed == {}
    assert set(raw_xml) <= {"Homo sapiens"}
    for species in config.SPECIES:
        if species != "Homo sapiens":
            expected = constants.WARN_BLAST_TIMEOUT.format(species=species)
            assert warnings_out.count(expected) == 1
