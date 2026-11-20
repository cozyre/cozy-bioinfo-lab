"""End-to-end orchestration of one analysis run (§6.1, §7.2, §12).

The pipeline only coordinates domain modules and injected clients. It never
imports Streamlit, constructs HTTP sessions or writes to disk, so it can be
exercised fully offline with fake clients and committed fixtures.
"""

from __future__ import annotations

import warnings
from dataclasses import asdict, replace
from datetime import datetime, timezone
from typing import Any

import Bio
import numpy
import pandas
import plotly
import requests

from tp53 import blast_parser, config, constants, conservation, msa, selection
from tp53.fasta_io import parse_fasta, parse_single_fasta, to_fasta, validate_protein
from tp53.models import (
    Alignment,
    BlastHit,
    BlastHsp,
    BlastResult,
    PipelineClients,
    PipelineInputs,
    PipelineResult,
    ProgressCallback,
    QueryInfo,
    SelectedHit,
    SequenceRecord,
    ValidationError,
)

#: Stable progress stages consumed by the UI's step indicator.
STAGES: tuple[str, ...] = (
    "validate",
    "retrieve",
    "blast",
    "select",
    "fetch_sequences",
    "msa",
    "conservation",
    "done",
)


def _library_versions() -> dict[str, str]:
    """Installed versions of the libraries recorded in run metadata."""
    modules = {
        "biopython": Bio,
        "requests": requests,
        "pandas": pandas,
        "numpy": numpy,
        "plotly": plotly,
    }
    versions: dict[str, str] = {}
    for name, module in modules.items():
        versions[name] = str(getattr(module, "__version__", "unknown"))
    return versions


def _now() -> str:
    """UTC ISO-8601 timestamp."""
    return datetime.now(timezone.utc).isoformat()


def _header_for(accession: str, organism: str) -> str:
    """``Species_name|Accession`` header ID (FR-13)."""
    if not organism:
        return accession
    return f"{organism}{config.HEADER_SEPARATOR}{accession}"


def _header_id(record: SequenceRecord, fallback_species: str = "") -> str:
    """Header ID for a record, falling back to the target species label."""
    return _header_for(record.accession, record.organism or fallback_species)


def _is_query_species(query: QueryInfo, species: str) -> bool:
    """Whether a target species is the query's own organism (§9.1)."""
    organism = query.record.organism.strip().casefold()
    return bool(organism) and organism == species.strip().casefold()


def _query_selected_hit(query: QueryInfo, species: str) -> SelectedHit:
    """Represent the query as its own species' selected sequence (§9.3)."""
    sequence = query.record.sequence
    hsp = BlastHsp(
        bit_score=0.0,
        evalue=0.0,
        identity=len(sequence),
        positives=len(sequence),
        gaps=0,
        align_len=len(sequence),
        query_start=1,
        query_end=len(sequence),
        hit_start=1,
        hit_end=len(sequence),
        qseq=sequence,
        hseq=sequence,
        midline="",
    )
    hit = BlastHit(
        accession=query.record.accession,
        title=query.record.description,
        organism=query.record.organism or species,
        length=len(sequence),
        hsps=[hsp],
    )
    return SelectedHit(
        species=species,
        hit=hit,
        query_coverage=1.0,
        percent_identity=100.0,
        reason=constants.SELECT_REASON_QUERY,
        candidates=[hit],
    )


def _resolve_query(
    inputs: PipelineInputs, clients: PipelineClients, emit: ProgressCallback
) -> QueryInfo:
    """Validate the query source and retrieve one validated query (§7.1, §12)."""
    supplied = sum(
        source is not None
        for source in (inputs.query_info, inputs.query_accession, inputs.query_fasta)
    )
    if supplied != 1:
        raise ValidationError(constants.ERR_NO_QUERY_SOURCE)

    if inputs.query_info is not None:
        emit("validate", constants.PROGRESS_VALIDATE_SNAPSHOT)
        query = inputs.query_info
    elif inputs.query_accession is not None:
        accession = inputs.query_accession.strip()
        if not accession:
            raise ValidationError(constants.ERR_EMPTY_ACCESSION)
        source = inputs.accession_source
        emit(
            "validate",
            constants.PROGRESS_VALIDATE_ACCESSION.format(
                accession=accession, source=source
            ),
        )
        if source == "uniprot":
            query = clients.fetch_uniprot(accession)
        elif source == "ncbi":
            query = clients.fetch_ncbi_query_info(accession)
        else:
            raise ValidationError(
                constants.ERR_UNKNOWN_ACCESSION_SOURCE.format(source=source)
            )
    else:
        record = parse_single_fasta(inputs.query_fasta or "")
        record.organism = inputs.query_organism
        emit(
            "validate",
            constants.PROGRESS_VALIDATE_FASTA.format(identifier=record.id),
        )
        query = QueryInfo(record=record, features=[])

    if query.record is None:
        raise ValidationError(constants.ERR_NO_QUERY_RECORD)
    validate_protein(query.record)
    emit(
        "retrieve",
        constants.PROGRESS_RETRIEVE_DONE.format(accession=query.record.accession),
    )
    return query


def _submit_jobs(
    inputs: PipelineInputs,
    clients: PipelineClients,
    query: QueryInfo,
    emit: ProgressCallback,
    warnings_out: list[str],
) -> dict[str, str]:
    """Submit only the species without uploaded XML or a saved RID."""
    uploaded = inputs.blast_xml_by_species or {}
    existing = inputs.existing_rids or {}
    pending = [
        species
        for species in inputs.species
        if not _is_query_species(query, species)
        and species not in uploaded
        and species not in existing
    ]
    if not pending:
        return {}
    jobs = [(species, query.record.sequence) for species in pending]
    try:
        return dict(clients.submit_blastp_staggered(jobs))
    except Exception as exc:
        # A single submission failure must not discard the other species. The
        # staggered batch is preferred; this per-job retry isolates failures.
        rids: dict[str, str] = {}
        for job in jobs:
            try:
                rids.update(clients.submit_blastp_staggered([job]))
            except Exception as species_exc:
                warnings_out.append(
                    constants.WARN_BLAST_FAILED.format(
                        species=job[0], error=species_exc
                    )
                )
        return rids


def _collect_blast(
    inputs: PipelineInputs,
    clients: PipelineClients,
    query: QueryInfo,
    emit: ProgressCallback,
    warnings_out: list[str],
) -> tuple[dict[str, str], dict[str, BlastResult]]:
    """Obtain and parse one BLAST XML per species, tolerating local failures."""
    uploaded = inputs.blast_xml_by_species or {}
    rids = dict(inputs.existing_rids or {})
    new_rids = _submit_jobs(inputs, clients, query, emit, warnings_out)
    rids.update(new_rids)

    raw_xml: dict[str, str] = {}
    parsed: dict[str, BlastResult] = {}
    for species in inputs.species:
        if _is_query_species(query, species):
            # The query organism is represented by the query itself (§9.1).
            if species in uploaded:
                raw_xml[species] = uploaded[species]
            emit("blast", constants.PROGRESS_BLAST_QUERY.format(species=species))
            continue
        try:
            if species in uploaded:
                emit("blast", constants.PROGRESS_BLAST_UPLOAD.format(species=species))
                xml = uploaded[species]
            elif species in rids:
                rid = rids[species]
                if species in new_rids:
                    emit(
                        "blast",
                        constants.PROGRESS_BLAST_SUBMITTED.format(
                            species=species, rid=rid
                        ),
                    )
                else:
                    emit(
                        "blast",
                        constants.PROGRESS_BLAST_RESUME.format(
                            species=species, rid=rid
                        ),
                    )
                emit("blast", constants.PROGRESS_BLAST_POLLING.format(species=species))
                xml = clients.poll_blast(rid, species)
                emit("blast", constants.PROGRESS_BLAST_READY.format(species=species))
            else:
                raise ValidationError(f"No BLAST result available for {species}.")
            result = blast_parser.parse_blast_xml(xml, species)
        except Exception as exc:
            warnings_out.append(
                constants.WARN_BLAST_FAILED.format(species=species, error=exc)
            )
            continue
        raw_xml[species] = xml
        parsed[species] = result
    return raw_xml, parsed


def _select_hits(
    inputs: PipelineInputs,
    blast_results: dict[str, BlastResult],
    query: QueryInfo,
    emit: ProgressCallback,
    warnings_out: list[str],
) -> list[SelectedHit]:
    """Apply §9.2 per species; a species without a valid hit is skipped."""
    selected: list[SelectedHit] = []
    for species in inputs.species:
        if _is_query_species(query, species):
            selected.append(_query_selected_hit(query, species))
            emit(
                "select",
                constants.PROGRESS_SELECT_SELECTED.format(
                    species=species, accession=query.record.accession
                ),
            )
            continue
        result = blast_results.get(species)
        if result is None:
            continue
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            hit = selection.select_best_hit(result, inputs.selection_params)
        for warning in caught:
            warnings_out.append(str(warning.message))
        if hit is None:
            warnings_out.append(constants.WARN_NO_HIT.format(species=species))
            emit("select", constants.PROGRESS_SELECT_MISSING.format(species=species))
            continue
        selected.append(hit)
        emit(
            "select",
            constants.PROGRESS_SELECT_SELECTED.format(
                species=species, accession=hit.hit.accession
            ),
        )
    return selected


def _sequence_from_supplied(text: str, accession: str, species: str) -> str:
    """Accept either a FASTA string or a bare sequence from a snapshot."""
    if text.lstrip().startswith(">"):
        records = parse_fasta(text)
        if len(records) != 1:
            raise ValidationError(f"Expected one FASTA record for {accession}.")
        return records[0].sequence
    return "".join(text.split())


def _assemble_records(
    inputs: PipelineInputs,
    clients: PipelineClients,
    query: QueryInfo,
    selected: list[SelectedHit],
    emit: ProgressCallback,
    warnings_out: list[str],
) -> tuple[list[SequenceRecord], list[SelectedHit]]:
    """Build the deduplicated query-plus-ortholog input set (§9.3, FR-13)."""
    query_record = replace(query.record, id=_header_id(query.record))
    records = [query_record]
    seen = {query_record.accession}
    kept: list[SelectedHit] = []

    for hit in selected:
        accession = hit.hit.accession
        if accession in seen:
            # The query stands in for its own species; keep it in `selected`
            # for the UI even though it adds no extra input record.
            if _is_query_species(query, hit.species):
                kept.append(hit)
            continue
        seen.add(accession)
        supplied = inputs.full_sequences or {}
        try:
            if accession in supplied:
                sequence = _sequence_from_supplied(
                    supplied[accession], accession, hit.species
                )
                record = SequenceRecord(
                    id=_header_for(accession, hit.species),
                    accession=accession,
                    organism=hit.species,
                    description=hit.hit.title,
                    sequence=sequence,
                )
                emit(
                    "fetch_sequences",
                    constants.PROGRESS_FETCH_SUPPLIED.format(accession=accession),
                )
            else:
                record = clients.fetch_protein(accession)
                emit(
                    "fetch_sequences",
                    constants.PROGRESS_FETCH_SEQUENCE.format(accession=accession),
                )
                record = replace(
                    record,
                    id=_header_id(record, hit.species),
                    organism=record.organism or hit.species,
                )
        except Exception as exc:
            warnings_out.append(
                constants.WARN_SEQUENCE_FAILED.format(accession=accession, error=exc)
            )
            seen.discard(accession)
            continue
        records.append(record)
        kept.append(hit)
    return records, kept


def _build_metadata(
    inputs: PipelineInputs,
    query: QueryInfo,
    started: str,
    raw_xml: dict[str, str],
    raw_msa_fasta: str,
    selected: list[SelectedHit],
    alignment: Alignment | None,
) -> dict[str, Any]:
    """Record provenance, applied settings and artefact availability."""
    has_overrides = bool(
        inputs.query_info
        or inputs.blast_xml_by_species
        or inputs.msa_fasta
        or inputs.full_sequences
    )
    mode = "example" if inputs.query_info else "upload" if has_overrides else "remote"
    if inputs.query_info is not None:
        source = "supplied"
    elif inputs.query_fasta is not None:
        source = "fasta"
    else:
        source = inputs.accession_source
    return {
        "start_time": started,
        "end_time": _now(),
        "started_at": started,
        "finished_at": _now(),
        "query": query.record.accession,
        "query_id": query.record.id,
        "query_accession": query.record.accession,
        "query_organism": query.record.organism,
        "query_source": source,
        "source": source,
        "mode": mode,
        "species": list(inputs.species),
        "blast_params": asdict(inputs.blast_params),
        "selection_params": asdict(inputs.selection_params),
        "conservation_params": asdict(inputs.conservation_params),
        "smoothing_window": config.SMOOTHING_WINDOW_DEFAULT,
        "selected_accessions": [hit.hit.accession for hit in selected],
        "library_versions": _library_versions(),
        "raw_artifacts": {
            "blast_xml": bool(raw_xml),
            "blast_xml_species": sorted(raw_xml),
            "msa_fasta": bool(raw_msa_fasta),
            "aligned_columns": alignment.length if alignment else 0,
        },
        "raw_artefacts": bool(raw_xml) or bool(raw_msa_fasta),
        "has_raw_blast_xml": bool(raw_xml),
        "has_raw_msa_fasta": bool(raw_msa_fasta),
    }


def run_pipeline(
    inputs: PipelineInputs,
    clients: PipelineClients,
    progress_cb: ProgressCallback,
) -> PipelineResult:
    """Run all stages of §6.1 and return a single renderable result.

    Remote calls happen only through ``clients``. A species whose BLAST job,
    hit selection or sequence fetch fails is skipped with a warning; the run
    fails only when the query plus at least one ortholog cannot be assembled.

    Args:
        inputs: query, species, parameters and optional offline overrides.
        clients: injected network operations (real or fake).
        progress_cb: called as ``(stage, detail)``; stages follow :data:`STAGES`.

    Returns:
        Complete result with raw artefacts, metadata and warnings.

    Raises:
        ValidationError: invalid input, or fewer than two usable sequences.
    """
    started = _now()
    warnings_out: list[str] = []

    def emit(stage: str, detail: str) -> None:
        """Emit one progress event without changing the stable stage names."""
        progress_cb(stage, detail)

    query = _resolve_query(inputs, clients, emit)
    if not inputs.species:
        raise ValidationError(constants.ERR_NO_SPECIES)

    raw_xml, blast_results = _collect_blast(inputs, clients, query, emit, warnings_out)
    selected = _select_hits(inputs, blast_results, query, emit, warnings_out)
    if not selected:
        raise ValidationError(constants.ERR_NO_ORTHOLOGS)
    input_records, selected = _assemble_records(
        inputs, clients, query, selected, emit, warnings_out
    )
    if len(input_records) < 2:
        raise ValidationError(constants.ERR_NO_ORTHOLOGS)

    fasta_text = inputs.msa_fasta or to_fasta(input_records)
    if inputs.msa_fasta:
        emit(
            "msa",
            constants.PROGRESS_MSA_UPLOAD.format(count=len(input_records)),
        )
        msa_fasta = inputs.msa_fasta
    else:
        emit("msa", constants.PROGRESS_MSA_SERVICE.format(count=len(input_records)))
        msa_fasta = clients.run_clustalo(fasta_text)

    expected_ids = [record.id for record in input_records]
    expected_sequences = [record.sequence for record in input_records]
    aln = msa.parse_alignment(
        msa_fasta, ref_id=expected_ids[0], expected_ids=expected_ids
    )
    msa_warnings = msa.validate_alignment(aln, expected_ids, expected_sequences)
    warnings_out.extend(msa_warnings)
    aln, removed = msa.trim_all_gap_columns(aln)
    if removed:
        warnings_out.append(constants.WARN_TRIMMED_COLUMNS.format(count=removed))

    features = query.features if query.features else []
    cons = conservation.compute_conservation(
        aln, features=features, params=inputs.conservation_params
    )
    emit(
        "conservation",
        constants.PROGRESS_CONSERVATION_DONE.format(count=aln.length),
    )

    metadata = _build_metadata(
        inputs, query, started, raw_xml, msa_fasta, selected, aln
    )
    emit("done", constants.PROGRESS_DONE.format(count=len(input_records)))
    return PipelineResult(
        query=query,
        blast_results=blast_results,
        selected=selected,
        input_records=input_records,
        alignment=aln,
        conservation=cons,
        run_metadata=metadata,
        warnings=warnings_out,
        raw_blast_xml=raw_xml,
        raw_msa_fasta=msa_fasta,
    )
