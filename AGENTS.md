# AGENTS.md

## Project status and source of truth

- **No code exists yet.** This repository currently contains only the design document: `"Requirements & Design Document.md"` (note the spaces and `&` — quote the path in shell commands). Read it before implementing anything.
- Read first: §6 (architecture/layering), §7 (module design + public interfaces), §9 (algorithms — exact definitions, do not improvise), §13 (testing + acceptance criteria). For sequencing the work, follow `PLAN.md` (4 phases; each module ships with tests in the same commit).
- §16 lists open design decisions that are explicitly **not** confirmed. Do not silently pick one; ask the project owner.

## Stack

- Python 3.11, Streamlit, Biopython, requests, pandas, numpy, plotly. Dependencies pinned in `requirements.txt`.
- Entry point is `streamlit_app.py` **at the repository root** — required for Streamlit Community Cloud deployment (AC-7).

## Hard constraints from the design doc

- **Layering (§6.2):** no Streamlit imports outside `app/` and `streamlit_app.py`. Core modules (`pipeline`, domain modules, `clients`, `viz`) are framework-free and get HTTP sessions injected so tests can mock them.
- **No external binaries, ever.** BLAST, alignment (EBI Clustal Omega REST), everything installs from `requirements.txt`. Do not add `packages.txt`.
- **Tests run fully offline.** Record real BLAST XML, MSA FASTA and UniProt JSON once, commit as fixtures under `tests/data/`. Never call live services in tests; clients are tested with mocked HTTP sessions. Target ≥80% coverage of core modules.
- **Style:** type hints and docstrings on all public functions, dataclasses for models, no magic numbers — all defaults/thresholds/URLs/species live in `config.py`, UI strings in `constants.py`.
- **Secrets:** NCBI/EBI contact emails come from `st.secrets` (`NCBI_EMAIL` required). Never hard-code endpoints, credentials, or domain annotations.
- **Feature annotations are fetched at runtime from UniProt** — nothing domain-related is hard-coded. Feature coordinates are 1-based reference numbering, mapped through the alignment column map (§9.7); features are never matched on alignment columns.
- **Determinism:** best-hit selection (§9.2) must be deterministic (bit score → identity → `NP_` over `XP_` → accession alphabetical). De-gapping each alignment row must reproduce the input sequence exactly — this is a hard consistency check (§9.4).
- **Remote-service etiquette:** every HTTP call needs timeout + retry with backoff; send contact emails to NCBI/EBI; stagger BLAST submissions and respect current NCBI polling intervals (check the guidelines when implementing). Deployment target is ephemeral — no disk persistence, cache with `st.cache_data`.
