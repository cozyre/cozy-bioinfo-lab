# Implementation Plan

Agent-facing execution plan for the TP53 Conservation Analyzer. Condenses the milestone plan in §15 of `"Requirements & Design Document.md"` into four phases. The design document remains the source of truth for every algorithm and data model — this file only sequences the work.

**Before starting any phase:** read `AGENTS.md` and the design doc sections listed there (§6, §7, §9, §13). Commit fixtures once and never call live services in tests.

**Defaults for the open items (§16)** — proceed with these unless the project owner overrides: BLAST database `refseq_protein`; region threshold 100%; gaps stay in the identity denominator; bundled example snapshot and raw downloads are kept; species list editable with the five case species as default; UI in English.

---

## Phase 1 — Scaffold and offline core

**Goal:** the package exists, imports cleanly, and the offline logic (FASTA, BLAST XML, MSA) works against committed fixtures.

1. Scaffold the repository exactly as §6.3: root `streamlit_app.py` (stub), `requirements.txt` with pinned versions (`streamlit`, `biopython`, `requests`, `pandas`, `numpy`, `plotly`, `pytest`), `.streamlit/config.toml`, `src/tp53/` package, `app/`, `tests/`.
2. Write `config.py` (defaults, species list, URLs, thresholds, `BlastParams`/`SelectionParams` dataclasses), `constants.py` (UI strings, colour-blind-safe palette), `models.py` (all dataclasses from §7.1 verbatim).
3. Record fixtures in `tests/data/`: one real NCBI BLAST XML (five-species TP53 run), one UniProt JSON for P04637, one Clustal Omega MSA FASTA. Commit them.
4. Implement with tests in the same commit:
   - `fasta_io.py` — `parse_fasta`, `validate_protein` (one record, amino-acid alphabet, 20–5000 aa, raises `ValidationError`), `to_fasta`.
   - `blast_parser.py` — wrap `Bio.Blast.NCBIXML` so callers only see `BlastResult`; test hit count, scores, ranges against fixture values.
   - `msa.py` — `parse_alignment`, `validate_alignment`: equal row lengths, ID mapping (by order if service rewrote IDs), allowed characters (20 aa + `XBZJUO` + `-`, upper-cased), all-gap column removal with count, and the **hard check that de-gapping reproduces the input sequence** (§9.4). Also build `col_to_ref` map (§9.7).

**Done when:** `pytest` passes fully offline; AC-2-style properties hold on fixtures; package imports from a fresh `pip install -r requirements.txt`.

---

## Phase 2 — Clients and analysis

**Goal:** network I/O isolated in `clients/` and the scientific core (selection, conservation) correct on fixtures.

1. Implement `clients/` with an **injected HTTP session** everywhere:
   - `uniprot.py` — `fetch_query_info`: sequence, organism, gene, features (Domain, Region, DNA binding, Motif, Zinc finger) from `rest.uniprot.org`.
   - `ncbi.py` — `fetch_protein` (Entrez efetch, email + tool set), `submit_blastp` (one job per species, organism-restricted `entrez_query`, staggered submissions, XML output, hitlist size 10), `poll_blast` (respect NCBI's current polling guidelines — verify them at implementation time).
   - `ebi.py` — `run_clustalo`: submit → poll status → fetch aligned FASTA; email parameter required.
   - All calls: timeout + retry with exponential backoff; endpoints/credentials from `config.py` and secrets, never hard-coded.
2. Client tests with mocked HTTP sessions covering success, timeout and HTTP-error paths.
3. Implement `selection.py` per §9.2: E-value filter → coverage ≥ threshold → rank by bit score → identity → `NP_` over `XP_` → accession alphabetical. Deterministic; `reason` string; no-hit → `None` + warning. Test coverage filter, all tie-breaks, no-hit case.
4. Implement `conservation.py` per §9.5–§9.9: `column_identity` (gaps in denominator), `classify_columns` (fully/high/moderate/low/gapped), `find_regions` (threshold + min length, both coordinate systems), `pairwise_identity`, `summarise_features` (feature positions via reference numbering, plus an "outside features" row). Hand-computed unit tests on a small synthetic alignment; identical sequences → 100 everywhere; gap columns capped below 100; region boundary cases.
5. One manual live smoke test (real UniProt/NCBI/EBI round-trip) — documented in the PR, never part of the automated suite.

**Done when:** client and analysis tests pass offline; AC-3, AC-4, AC-5 hold on fixtures; the scientific sanity check direction is plausible (DNA-binding core high identity on fixture data).

---

## Phase 3 — Visualisation and app

**Goal:** the full pipeline runs end to end and every Must requirement is reachable from the UI.

1. Implement `viz.py` per §11: `conservation_profile` (V1, feature bands + shaded regions + threshold line, reference/alignment x-axis toggle, smoothing), `alignment_html` (V2, 60-column blocks, consensus + identity rows, colour by amino-acid property, position ruler in reference numbering), `blast_overview` (V3), `feature_bars` (V4), `identity_heatmap` (V5). All take result objects and return Plotly figures or HTML strings; render each from fixture results in tests.
2. Implement `pipeline.py`: `run_pipeline(inputs, clients, progress_cb) -> PipelineResult` orchestrating §6.1. Per-species BLAST failures are non-fatal if query + one species remain (warning); results carry `run_metadata` and `warnings`.
3. Implement the Streamlit app (`streamlit_app.py`, `app/pages.py`, `app/state.py`) per §10: sidebar (query radio accession/FASTA, species multiselect, source options, advanced settings, Run/Load example/Reset), five tabs, run states via `st.status`, everything rendered from one `PipelineResult` in `st.session_state`. **Only these two files may import Streamlit.**
4. Caching with `st.cache_data` (TTL) for client calls and parsed results keyed by input; keep BLAST RIDs in session state so page reruns resume polling (§12). No disk writes.
5. Bundle the example snapshot (A8) in the repo and wire **Load example** (FR-33) plus downloads of BLAST XML / MSA FASTA / identity CSV (FR-32).

**Done when:** AC-1 and AC-6 pass — full run with defaults yields one selected sequence per species (five total), and every Must requirement is reachable in both remote and example modes.

---

## Phase 4 — Hardening and deployment

**Goal:** error paths are complete and the app is live on Streamlit Community Cloud.

1. Cover every case in §12: invalid input blocks with a specific message; service timeouts surface an upload-the-file fallback; malformed uploaded XML keeps previous results; FASTA-only queries skip feature overlays with a notice; MSA validation failures name the offending sequence.
2. Verify performance budget: local computation under 2 s (parse, conservation, plots); end-to-end ≈ 5 min for five species with visible progress throughout.
3. Security pass: sanitise all input, never execute uploaded content, enforce upload size limits, no secrets anywhere in the repo or logs.
4. Run the acceptance checklist §13 (AC-1 through AC-6) against the deployed behaviour; fix gaps.
5. Deploy to Streamlit Community Cloud from the repository with only `requirements.txt` and `st.secrets` (AC-7). Confirm the deployed app works in both remote and example modes, including a cold start from sleep.

**Done when:** AC-1…AC-7 all pass; the live URL is recorded in the README.

*(Non-agent milestone: the project owner runs the final case study and writes the report manually — no app feature required.)*
