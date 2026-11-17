# TP53 Conservation Analyzer: Requirements & Design Document

Version 1.0 (draft) · 17 November 2026 · Audience: the project owner and the coding agents that will implement it

## 1. Purpose and scope

**Background.** A bioinformatics lab analyses disease-related proteins by hand: reading BLAST output, comparing sequences, calculating residue conservation, drawing charts and writing reports as separate steps. This project automates that pipeline in a Streamlit web app, using TP53 in five species as the study case.

**Research question.** Are there parts of TP53 that are conserved across *Homo sapiens*, *Pan troglodytes*, *Mus musculus*, *Rattus norvegicus* and *Bos taurus*, and how can that result help explain the protein's biological function?

**Goals**

- G1: Take a human TP53 query (accession or FASTA), run BLASTp remotely and automatically select the best ortholog hit for each target species.
- G2: Parse the BLAST XML, build a multiple sequence alignment (MSA) and process it.
- G3: Compute per-position percent identity, detect conserved regions and relate them to annotated functional domains.
- G4: Present everything in an interactive UI so the researcher can answer the research question directly from the app.

**In scope:** the pipeline above, interactive visualisation, input validation, error handling, deployment on Streamlit Community Cloud.

**Out of scope:** report generation or export (Tugas 7 is written manually), user accounts, a database, phylogenetic trees, 3D structure, variant pathogenicity prediction, local BLAST or local alignment binaries.

### 1.1 Task coverage

| Task | Requirement | Delivered by |
| --- | --- | --- |
| 1. Implement in Python | Python backend, Streamlit UI | Whole project |
| 2. Modular program structure | Layered package, framework-free core | Section 6 and 7 |
| 3. Read BLAST XML | Parse remote BLAST results (and uploaded XML) | `blast_parser`, FR-07 to FR-12 |
| 4. Process MSA | Build, parse, validate and map the MSA | `msa`, FR-13 to FR-17 |
| 5. Residue conservation | Percent identity per column, regions, domain summary | `conservation`, FR-18 to FR-23 |
| 6. Visualisation | Interactive charts and alignment viewer | `viz`, FR-24 to FR-29 |
| 7. Report | Written manually from app results (not an app feature) | Appendix A helps with this |

## 2. Users and use cases

**Primary user:** a student or researcher who knows basic molecular biology but does not want to run BLAST, alignment and plotting by hand.

| ID | Use case |
| --- | --- |
| UC1 | Run the full pipeline for human TP53 (default P04637) against the five case species and read the conservation results. |
| UC2 | Run the pipeline with a different accession or pasted FASTA. |
| UC3 | Inspect BLAST candidate hits and see why a specific hit was selected per species. |
| UC4 | Inspect the alignment residue by residue. |
| UC5 | See which conserved regions overlap annotated domains such as the DNA-binding domain. |
| UC6 | Load previously produced BLAST XML or MSA FASTA instead of waiting for remote jobs. |

## 3. Decisions, assumptions and constraints

### 3.1 Confirmed decisions

| Topic | Decision |
| --- | --- |
| Data acquisition | The user supplies a UniProt/NCBI accession or a FASTA sequence. For an accession the app fetches the sequence from UniProt/NCBI. BLAST runs remotely at NCBI. |
| Species selection | BLASTp from the query (human TP53); the best hit for each target species is selected automatically. |
| MSA | Remote EBI Clustal Omega REST service (recommended: no binaries needed on Community Cloud), with optional upload of a ready MSA FASTA. |
| Conservation metric | Percent identity per alignment column only. |
| Reporting | No report export. Results are shown in the app. |
| Deployment | Streamlit Community Cloud. |
| Contact email | Each user enters their own NCBI/EBI contact email in the app before running (owner decision, 2026-11-17). Nothing identifying is stored in the repository; an optional `st.secrets` entry may pre-fill the field for the deployer's convenience. |

### 3.2 Assumptions (please correct any that are wrong)

- A1: UI text is in English. All strings live in one constants module so localising to Indonesian later is easy.
- A2: Default query is UniProt P04637 (human TP53, canonical isoform, 393 aa).
- A3: BLAST program is blastp, default database `refseq_protein`, default E-value cutoff 1e-10, both configurable.
- A4: Domain and region annotations come from the UniProt features of the query at runtime. Nothing is hard-coded.
- A5: Python 3.11 with Biopython, requests, pandas, numpy and plotly.
- A6: Single-user sessions, no authentication, no persistent storage.
- A7: Downloading raw artefacts (BLAST XML, MSA FASTA, identity CSV) is allowed as a convenience. It is not a report.
- A8: A bundled example snapshot (BLAST XML, MSA, UniProt JSON for TP53) is included as a fallback if remote services are down or slow. This is my recommendation, not something you requested.

### 3.3 Constraints

- Streamlit Community Cloud: ephemeral filesystem, limited CPU and memory, apps sleep when idle, deployment from a GitHub repository, secrets through `st.secrets`. No reliance on background workers or local disk persistence.
- External services (NCBI BLAST, NCBI Entrez, UniProt, EBI) are slow, rate limited and occasionally unavailable. The design must tolerate this.
- No external binaries: everything must install from `requirements.txt`.

## 4. Functional requirements

Priority: **Must**, **Should**, **Could**.

### 4.1 Input and retrieval

| ID | Requirement | Priority | Task |
| --- | --- | --- | --- |
| FR-01 | Accept a UniProt or NCBI protein accession as the query. Default P04637 is prefilled. | Must | 1 |
| FR-02 | Accept a pasted or uploaded single-record protein FASTA as an alternative query. | Must | 1 |
| FR-03 | Validate input: accession format, exactly one FASTA record, amino-acid alphabet only, length 20 to 5,000 residues. Show a clear message on failure. | Must | 1 |
| FR-04 | Let the user choose target species (multiselect, default the five case species) and set the query organism (automatic for accessions, manual for FASTA). | Must | 1 |
| FR-05 | Advanced settings: BLAST database, E-value cutoff, minimum query coverage, conservation threshold, minimum region length, smoothing window. | Should | 1 |
| FR-06 | Fetch sequence and metadata (organism, gene, feature annotations) from UniProt for UniProt accessions, or from NCBI Entrez for NCBI accessions. | Must | 1 |

### 4.2 BLAST

| ID | Requirement | Priority | Task |
| --- | --- | --- | --- |
| FR-07 | Submit one remote BLASTp job per target species (organism restricted), poll until complete and report progress per species. | Must | 3 |
| FR-08 | Keep the raw BLAST XML of each job in memory for parsing and optional download. | Must | 3 |
| FR-09 | Parse BLAST XML into structured hits: accession, title, subject length, bit score, E-value, identities, positives, gaps, alignment length, query and subject ranges, aligned strings. | Must | 3 |
| FR-10 | Optionally accept an uploaded BLAST XML file instead of running BLAST. | Should | 3 |
| FR-11 | Select the best hit per species by the rules in 9.2. Show all candidates and the reason for the choice. | Must | 3 |
| FR-12 | Fetch the full-length sequence of each selected hit from NCBI Entrez, because BLAST XML only contains the aligned segment. | Must | 3 |

### 4.3 MSA

| ID | Requirement | Priority | Task |
| --- | --- | --- | --- |
| FR-13 | Build a multi-FASTA of the query plus selected sequences, with headers formatted `Species_name\|Accession`. | Must | 4 |
| FR-14 | Submit it to the EBI Clustal Omega REST service, poll until complete and retrieve the FASTA alignment. | Must | 4 |
| FR-15 | Optionally accept an uploaded MSA FASTA instead of calling the service. | Should | 4 |
| FR-16 | Parse and validate the MSA: all rows equal length, IDs map back to input records, valid characters only. Drop columns that are gaps in every row. Report alignment length, gap percentage and ungapped length per sequence. | Must | 4 |
| FR-17 | Build a map from alignment column to reference (query) residue position, with `None` where the reference has a gap. | Must | 4 |

### 4.4 Conservation

| ID | Requirement | Priority | Task |
| --- | --- | --- | --- |
| FR-18 | Compute percent identity and gap fraction for every alignment column (9.5). | Must | 5 |
| FR-19 | Classify columns: fully conserved, high, moderate, low, gapped. | Must | 5 |
| FR-20 | Detect conserved regions (9.6) in alignment and reference coordinates. | Must | 5 |
| FR-21 | Compute a pairwise percent identity matrix between sequences (9.8). | Should | 5 |
| FR-22 | For each UniProt feature (domain, region, motif, DNA-binding), report mean identity and the share of fully conserved positions, and compare with the rest of the protein. | Must | 5 |
| FR-23 | Show summary metrics: mean identity, number and share of fully conserved positions, longest conserved region, number of regions. | Must | 5 |

### 4.5 Visualisation and UI

| ID | Requirement | Priority | Task |
| --- | --- | --- | --- |
| FR-24 | Conservation profile plot with feature bands and conserved-region shading (V1). | Must | 6 |
| FR-25 | Alignment viewer with consensus and identity row (V2). | Must | 6 |
| FR-26 | BLAST overview table and bar chart for selected hits (V3). | Must | 6 |
| FR-27 | Feature comparison bar chart (V4). | Must | 6 |
| FR-28 | Pairwise identity heatmap (V5). | Should | 6 |
| FR-29 | Optional overlay of user-given residue positions, for example known cancer hotspots (V6). | Could | 6 |
| FR-30 | Step-by-step progress indicator, partial results shown as each stage completes, clear error messages, reset button. | Must | 1 |
| FR-31 | Cache results by input so reruns and widget changes do not repeat remote calls. | Must | 1 |
| FR-32 | Download raw artefacts (BLAST XML, MSA FASTA, identity CSV). | Could | 1 |
| FR-33 | Load the bundled example snapshot with one click (A8). | Should | 1 |

## 5. Non-functional requirements

| Area | Requirement |
| --- | --- |
| Performance | Typical end-to-end run (five species) within about 5 minutes, dominated by remote BLAST. All local computation (parsing, conservation, plots) under 2 seconds. The UI never freezes silently: every long step shows progress. |
| Reliability | Timeouts and retries with exponential backoff on every HTTP call. If one species fails, the pipeline continues if at least the query plus one other species remain, and shows a warning. |
| Modularity | Core logic imports no Streamlit. Only the `app/` package touches Streamlit. |
| Testability | Unit tests run fully offline using recorded fixtures. At least 80% coverage of core modules. |
| Reproducibility | A run-metadata panel shows query, BLAST database and parameters, selected accessions, timestamps and library versions. |
| Service etiquette | Send a contact email to NCBI and EBI as required, stagger BLAST submissions and respect the polling intervals in NCBI's current usage guidelines. Agents must check the guidelines when implementing. |
| Security | No secrets in the repository. Sanitise all user input. Never execute uploaded content. Limit upload size. |
| Maintainability | Type hints, docstrings, dataclasses for data models, pinned dependencies, centralised configuration. |
| Portability | No system binaries, so no `packages.txt` is needed on Community Cloud. |

## 6. Architecture

### 6.1 Pipeline

```
Input (accession / FASTA)
  -> Sequence retrieval (UniProt / NCBI)
  -> Remote BLASTp, one job per target species (NCBI)
  -> BLAST XML parsing
  -> Best-hit selection per species
  -> Full-length sequence fetch (NCBI Entrez)
  -> MSA (EBI Clustal Omega, or uploaded MSA)
  -> MSA parsing, validation, coordinate map
  -> Conservation (percent identity, regions, feature summary)
  -> Visualisation in Streamlit
```

### 6.2 Layers

| Layer | Responsibility | May import |
| --- | --- | --- |
| `app/` (UI) | Streamlit pages, widgets, session state, caching wrappers | `pipeline`, `viz` |
| `pipeline` | Orchestrates the stages, returns one `PipelineResult` | domain modules, `clients` |
| Domain modules | Pure functions: parsing, selection, MSA processing, conservation | `models`, `config` |
| `clients` | All network I/O: UniProt, NCBI BLAST and Entrez, EBI | `models`, `config` |
| `viz` | Functions that take result objects and return Plotly figures or HTML strings | `models` |

HTTP sessions are injected into clients so tests can replace them with mocks.

### 6.3 Project structure

```
tp53_analyzer/
  streamlit_app.py          # entry point for Community Cloud
  requirements.txt
  .streamlit/config.toml
  README.md
  src/tp53/
    config.py               # defaults, species list, URLs, thresholds
    constants.py            # UI strings, colour schemes
    models.py               # dataclasses (section 7)
    fasta_io.py             # FASTA parse/validate/write
    clients/
      uniprot.py            # sequence + feature fetch
      ncbi.py               # Entrez efetch, remote BLAST submit/poll
      ebi.py                # Clustal Omega submit/poll/result
    blast_parser.py         # BLAST XML -> BlastResult
    selection.py            # best-hit selection
    msa.py                  # MSA parse, validate, column map
    conservation.py         # identity, classes, regions, feature summary
    viz.py                  # Plotly figures, alignment HTML
    pipeline.py             # orchestrator
  app/
    pages.py                # tab renderers
    state.py                # session_state helpers, cache wrappers
  tests/
    data/                   # recorded BLAST XML, MSA FASTA, UniProt JSON
    test_*.py
```

## 7. Module design

### 7.1 Data models (`models.py`)

| Model | Fields |
| --- | --- |
| `SequenceRecord` | `id`, `accession`, `organism`, `description`, `sequence` |
| `DomainFeature` | `type`, `name`, `start`, `end` (1-based, reference coordinates) |
| `QueryInfo` | `record: SequenceRecord`, `features: list[DomainFeature]` |
| `BlastHsp` | `bit_score`, `evalue`, `identity`, `positives`, `gaps`, `align_len`, `query_start`, `query_end`, `hit_start`, `hit_end`, `qseq`, `hseq`, `midline` |
| `BlastHit` | `accession`, `title`, `organism`, `length`, `hsps: list[BlastHsp]` with property `best_hsp` |
| `BlastResult` | `query_id`, `query_len`, `database`, `program`, `hits: list[BlastHit]`, `species` (the organism filter used) |
| `SelectedHit` | `species`, `hit: BlastHit`, `query_coverage`, `percent_identity`, `reason`, `candidates: list[BlastHit]` |
| `Alignment` | `ids`, `rows` (aligned strings), `length`, `ref_index` (row of the query), `col_to_ref: list[int or None]` |
| `ConservationResult` | `identity: list[float]`, `gap_fraction: list[float]`, `classes: list[str]`, `regions: list[ConservedRegion]`, `pairwise: DataFrame`, `feature_summary: DataFrame`, `summary: dict` |
| `ConservedRegion` | `aln_start`, `aln_end`, `ref_start`, `ref_end`, `length`, `mean_identity` |
| `PipelineResult` | all of the above plus `run_metadata` and `warnings` |

### 7.2 Public interfaces

```python
# fasta_io.py
parse_fasta(text: str) -> list[SequenceRecord]
validate_protein(record: SequenceRecord, min_len=20, max_len=5000) -> None   # raises ValidationError
to_fasta(records: list[SequenceRecord]) -> str

# clients/uniprot.py
fetch_query_info(accession: str, session) -> QueryInfo

# clients/ncbi.py
fetch_protein(accession: str, session, email: str) -> SequenceRecord
submit_blastp(seq: str, organism: str, params: BlastParams, session) -> str      # returns RID
poll_blast(rid: str, session, timeout_s: int) -> str                              # returns XML text

# clients/ebi.py
run_clustalo(fasta: str, email: str, session, timeout_s: int) -> str             # returns aligned FASTA

# blast_parser.py
parse_blast_xml(xml_text: str, species: str) -> BlastResult

# selection.py
select_best_hit(result: BlastResult, params: SelectionParams) -> SelectedHit | None

# msa.py
parse_alignment(text: str, ref_id: str) -> Alignment
validate_alignment(aln: Alignment, expected_ids: list[str]) -> list[str]          # returns warnings

# conservation.py
column_identity(aln: Alignment) -> tuple[list[float], list[float]]
classify_columns(identity, gap_fraction, params) -> list[str]
find_regions(aln, identity, gap_fraction, params) -> list[ConservedRegion]
pairwise_identity(aln: Alignment) -> DataFrame
summarise_features(aln, identity, features) -> DataFrame

# pipeline.py
run_pipeline(inputs: PipelineInputs, clients, progress_cb) -> PipelineResult

# viz.py
conservation_profile(result) -> plotly Figure
alignment_html(aln, identity, block=60) -> str
blast_overview(selected) -> plotly Figure
feature_bars(result) -> plotly Figure
identity_heatmap(result) -> plotly Figure
```

Biopython's `Bio.Blast.NCBIXML` is the recommended parser for BLAST XML, wrapped so the rest of the code only sees the project's own dataclasses.

## 8. External integrations

| Service | Used for | Notes for implementers |
| --- | --- | --- |
| UniProt REST (`rest.uniprot.org/uniprotkb/{accession}`) | Query sequence, organism, gene and feature annotations (JSON) | Parse feature types: Domain, Region, DNA binding, Motif, Zinc finger. |
| NCBI BLAST URL API or Biopython `NCBIWWW.qblast` | Remote BLASTp | Restrict by `entrez_query` such as `"Mus musculus"[Organism]`, request XML output and a small hit list (10). Stagger submissions and poll at the interval NCBI recommends. |
| NCBI Entrez efetch (db=protein) | Full-length subject sequences, NCBI-accession queries | Set contact email and tool name. API key optional via `st.secrets`. |
| EBI Clustal Omega REST (`ebi.ac.uk/Tools/services/rest/clustalo`) | MSA | Asynchronous: run, poll status, fetch the FASTA alignment result. Requires an email parameter. Verify parameter and result-type names against current EBI documentation. |

All endpoints are read from `config.py`, never hard-coded in client logic. The contact email is supplied by the user at runtime (owner decision, 2026-11-17) and is never stored in the repository.

## 9. Algorithms and rules

### 9.1 BLAST strategy

One BLASTp job per target species, with an organism restriction, instead of one large search filtered afterwards. A single global search can miss a species in its top hits, while per-species jobs reliably return that species' orthologs. The query organism (human) is not BLASTed. The query sequence is used directly as the human entry. Jobs are submitted in a staggered fashion and polled in one loop to keep total wall time near the slowest job.

### 9.2 Best-hit selection

For each species:

1. Keep hits with E-value at or below the cutoff.
2. Compute query coverage as aligned query length of the best HSP divided by query length. Keep hits at or above the coverage threshold (default 70%). This removes fragments and unrelated paralogs.
3. Rank by bit score (descending), then percent identity (descending), then prefer accessions beginning `NP_` over `XP_`, then accession alphabetically so the result is deterministic.
4. The top candidate is selected. The `reason` string states score, identity, coverage and the tie-break used if any.
5. If no candidate remains, the species is recorded as not found and a warning is raised.

Known risk: a species may return several isoforms with similar scores. The ranking above favours the full-length one. The UI shows all candidates so the user can verify.

### 9.3 Sequence set assembly

Set = query sequence plus the full-length sequence of each selected hit (FR-12). Duplicate accessions are removed. Header format is `Species_name|Accession`.

### 9.4 MSA validation

- Alignment rows equal in length, and the number of rows equals the number of input sequences.
- IDs match the input headers. If the service rewrites IDs, map them by order.
- Characters limited to the 20 amino acids, `X`, `B`, `Z`, `J`, `U`, `O`, and `-` (case-insensitive, normalised to upper case).
- Remove columns that are `-` in every row. Record how many were removed.
- Removing gaps from each row must reproduce the original input sequence. This is a hard consistency check.

### 9.5 Percent identity per column

For column *j* with *N* sequences:

- Count each standard residue, ignoring gaps and ambiguous codes in the count.
- `identity_j = (highest residue count / N) * 100`
- `gap_fraction_j = gaps in column j / N`

Gaps stay in the denominator, so a column with a gap can never reach 100%. With five sequences the possible values are 20, 40, 60, 80 and 100 (or lower when gaps are present). The UI states this explicitly.

Column classes: **fully conserved** (100%, no gaps), **high** (80%), **moderate** (60%), **low** (40% or lower), **gapped** (gap fraction above 0 and not otherwise conserved). Thresholds come from `config.py`.

### 9.6 Conserved region detection

A conserved region is a run of consecutive alignment columns where identity is at or above the threshold (default 100%) and gap fraction is 0. Runs shorter than the minimum length (default 5) are discarded. For each region, report alignment start and end, reference start and end (from the column map, using the first and last columns that have a reference position), length and mean identity.

### 9.7 Coordinate mapping

Walk the reference row. Each non-gap character increments the reference position counter (1-based). `col_to_ref[j]` holds that position, or `None` where the reference has a gap. UniProt feature coordinates refer to this reference numbering, so features are overlaid by looking up reference positions, never alignment columns.

This mapping is only meaningful when the reference is the UniProt entry the features came from. If the user supplies a bare FASTA, feature overlays are skipped and the UI says so.

### 9.8 Pairwise identity

For each pair of sequences, count columns where both rows have a residue and the residues are equal, divided by the number of columns where both rows have a residue, times 100.

### 9.9 Feature summary

For each feature, take the reference positions start to end, collect the corresponding columns, and compute mean identity and the share of fully conserved columns. Add a final row for all positions outside any listed feature. This table is the main quantitative answer to the research question.

## 10. UI design

Layout: sidebar for inputs and settings, main area with tabs. Everything is rendered from one `PipelineResult` held in `st.session_state`.

**Sidebar**

- Query input: accession text box (default P04637) or FASTA paste/upload (radio switch).
- Query organism (read-only when derived from an accession).
- Target species multiselect.
- Source options: run BLAST and MSA remotely (default) or upload BLAST XML and/or MSA FASTA.
- Advanced settings expander (FR-05).
- Buttons: **Run analysis**, **Load example**, **Reset**.

**Tabs**

| Tab | Content |
| --- | --- |
| Overview | Run status per stage, run-metadata panel, summary metrics (FR-23), warnings |
| BLAST | Per-species candidate table, selected hit highlighted with the selection reason, overview chart (V3), expandable raw hit alignment |
| Alignment | Alignment viewer (V2) with consensus and identity row, alignment statistics |
| Conservation | Conservation profile (V1) with smoothing slider, regions table, feature summary table and chart (V4), pairwise heatmap (V5) |
| Downloads | Raw artefacts (FR-32), if included |

**Run states:** idle, running (a `st.status` block listing stages with per-species BLAST progress), success, partial success (warnings shown at the top), failed (message and suggested action, such as uploading a file instead).

## 11. Visualisation specification

| ID | Visual | Details |
| --- | --- | --- |
| V1 | Conservation profile | X axis: reference position (toggle for alignment column). Y axis: percent identity 0 to 100. Line or bars, optional smoothing window. Domain and feature bands in a track below the plot. Conserved regions shaded. Threshold line. Hover shows position, residue, identity, gap fraction. |
| V2 | Alignment viewer | Monospace blocks of 60 columns, one row per sequence with labels, residues coloured by amino-acid property, consensus row, identity row (shaded by class), position ruler in reference numbering. Rendered as HTML so five sequences of about 400 columns stay fast. |
| V3 | BLAST overview | Bar chart per species of percent identity and bit score for the selected hit, plus the candidates table. |
| V4 | Feature comparison | Bar chart of mean identity per feature versus the rest of the protein, with fully conserved share as a secondary value. |
| V5 | Pairwise heatmap | Symmetric matrix of pairwise percent identity, annotated values. |
| V6 | Position overlay (Could) | User-entered residue numbers drawn as markers on V1. |

Colours use a colour-blind-safe palette defined in `constants.py`.

## 12. Error handling and edge cases

| Situation | Behaviour |
| --- | --- |
| Invalid accession or sequence | Block the run, show the specific validation message. |
| NCBI or EBI timeout or HTTP error | Retry with backoff, then show failure for that stage with the option to upload the missing file. |
| BLAST job returns no hits for a species | Warn, skip that species, continue if at least two sequences remain. |
| Several isoforms with equal scores | Deterministic tie-break (9.2) and all candidates displayed. |
| MSA row length mismatch or sequence mismatch after de-gapping | Fail with a clear message naming the offending sequence. |
| Uploaded BLAST XML malformed or empty | Show parse error, keep previous results untouched. |
| User supplies FASTA without a UniProt accession | Skip feature overlays and feature summary, show a notice. |
| Page refresh or session rerun during a long job | Job IDs are kept in session state so polling can resume. Completed stages come from cache. |
| Remote services unreachable | Offer the bundled example snapshot (FR-33). |

## 13. Testing and acceptance

**Unit tests (offline, fixtures in `tests/data`)**

- FASTA parsing and validation, including malformed input.
- BLAST XML parser against a recorded real XML file: hit count, scores and ranges match known values.
- Selection: coverage filter, tie-breaks, no-hit case.
- MSA: validation failures, all-gap column removal, de-gap reproduces input, column-to-reference map.
- Conservation: hand-computed identity for a small synthetic alignment, identical sequences give 100 everywhere, gap columns capped below 100, region detection boundaries.
- Clients: mocked HTTP responses for success, timeout and error paths.

**Acceptance criteria**

| ID | Criterion |
| --- | --- |
| AC-1 | With default input, the pipeline produces exactly one selected sequence per species, five sequences in total, each labelled with species and accession. |
| AC-2 | The alignment has five equal-length rows and removing gaps reproduces every input sequence. |
| AC-3 | The identity list has one value per alignment column, all between 0 and 100, and fully conserved columns have identity 100 and no gaps. |
| AC-4 | Conserved regions respect the minimum length and every column inside a region meets the threshold. |
| AC-5 | Feature overlays use reference numbering and match UniProt feature coordinates. |
| AC-6 | Every Must requirement is reachable through the UI without errors, using both remote mode and the example snapshot. |
| AC-7 | The app deploys on Streamlit Community Cloud from the repository with only `requirements.txt` and secrets. |

A scientific sanity check, not an automated test: the DNA-binding core of TP53 should show clearly higher mean identity than the N-terminal and C-terminal regions. If the results disagree, check the alignment and coordinate mapping before trusting them.

## 14. Deployment (Streamlit Community Cloud)

- Public or accessible GitHub repository, entry point `streamlit_app.py` at the repository root.
- `requirements.txt` with pinned versions: streamlit, biopython, requests, pandas, numpy, plotly (plus pytest for development).
- Contact email: each user enters their own NCBI/EBI contact email in the app sidebar before running — NCBI and EBI require it for etiquette, and it is never stored in the repository. A deployer may optionally pre-fill it via `st.secrets`; `NCBI_API_KEY` is optional (raises rate limits).
- Use `st.cache_data` with a time-to-live for client calls and parsed results. Do not write results to disk.
- Keep each remote call within sensible timeouts and surface progress, since apps on the free tier can sleep and sessions can drop.
- The example snapshot ships inside the repository so the app is demonstrable even when external services fail.

## 15. Implementation plan for the coding agents

**Ground rules**

- Follow the layering in 6.2. No Streamlit imports outside `app/` and `streamlit_app.py`.
- Type hints and docstrings on all public functions. Dataclasses for models.
- No magic numbers. Defaults live in `config.py`.
- Every module ships with tests in the same milestone.
- Record real BLAST XML, UniProt JSON and an MSA once and commit them as fixtures. Do not call live services in tests.

| Milestone | Scope | Done when |
| --- | --- | --- |
| M0 | Repository scaffold, config, constants, models | Package imports, tests run |
| M1 | `fasta_io`, `blast_parser`, `msa` with fixtures | Parser and validation tests pass |
| M2 | `clients` (UniProt, NCBI, EBI) with mocks | Client tests pass, one manual live smoke test |
| M3 | `selection`, `conservation`, regions, feature summary | AC-2 to AC-5 pass on fixtures |
| M4 | `viz` | All figures render from fixture results |
| M5 | `pipeline` and Streamlit UI, caching, session state | AC-1 and AC-6 pass in remote and example modes |
| M6 | Hardening: error paths, progress, deployment to Community Cloud | AC-7 passes |
| M7 | You: run the case, write Tugas 7 | Report finished |

## 16. Open items to confirm

1. BLAST database: `refseq_protein` (assumed) or `swissprot` or `nr`. Swissprot gives one reviewed entry per species but may not cover every species.
2. Conserved region default threshold: strict 100% (assumed) or 80%. With five sequences, 100% is the cleanest definition.
3. Should gaps stay in the identity denominator (assumed) or be excluded?
4. Keep the bundled example snapshot (A8) and raw downloads (FR-32)?
5. Should the target species list stay fixed to the five case species, or be freely editable?
6. English UI (A1) or Indonesian?
