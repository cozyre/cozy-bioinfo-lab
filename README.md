# TP53 Conservation Analyzer

A Streamlit web app that automates a protein conservation analysis for **TP53** across five species — *Homo sapiens*, *Pan troglodytes*, *Mus musculus*, *Rattus norvegicus* and *Bos taurus* — and answers one research question from the browser:

> Which parts of TP53 are conserved across these species, and what does that tell us about the protein's function?

It replaces a manual lab workflow (running BLAST, aligning sequences, computing conservation by hand, plotting, writing up) with a single automated pipeline.

> **Status:** 🚧 under construction. This repository currently contains the design document and implementation plan; the app itself is being built phase by phase (see `PLAN.md`).

---

## What the app does

For a protein query (by default human TP53, UniProt **P04637**), the app:

1. **Fetches the sequence and annotations** from UniProt (or NCBI) — including functional domains such as the DNA-binding domain.
2. **Runs remote BLASTp** at NCBI — one search per target species — and automatically picks the best ortholog hit for each species, showing all candidates and why each one was chosen.
3. **Builds a multiple sequence alignment (MSA)** of the query plus the selected sequences using the EBI Clustal Omega REST service.
4. **Computes conservation** — percent identity for every alignment column, detects conserved regions, and compares conservation inside annotated functional domains versus the rest of the protein.
5. **Visualises everything** — interactive Plotly charts, a residue-level alignment viewer, and summary tables.

No installations, no command-line tools, no local BLAST binaries: everything runs through web services, so the app deploys as a plain Python package.

---

## Running the app

### Locally

```bash
pip install -r requirements.txt
streamlit run streamlit_app.py
```

Python 3.11 is recommended. NCBI and EBI ask that automated use identifies itself, so set your contact email in `.streamlit/secrets.toml`:

```toml
NCBI_EMAIL = "you@example.com"
# optional: NCBI_API_KEY = "..."   # raises NCBI rate limits
# optional: EBI_EMAIL = "..."      # defaults to NCBI_EMAIL
```

### On Streamlit Community Cloud

Deploy the repository directly — the entry point is `streamlit_app.py` at the repo root, and the same secrets are configured in the dashboard. No extra system packages are needed.

---

## Using the app

**Sidebar** (all inputs live here):

| Control | What it does |
| --- | --- |
| Query | A UniProt/NCBI accession (default `P04637`, prefilled) **or** a pasted/uploaded protein FASTA (one record, amino acids, 20–5,000 residues) |
| Target species | Multiselect, default the five case species |
| Source options | Run BLAST + MSA remotely (default), **or** upload a previously saved BLAST XML and/or MSA FASTA to skip the remote services |
| Advanced settings | BLAST database, E-value cutoff, query-coverage threshold, conservation threshold, minimum region length, smoothing window |
| Buttons | **Run analysis** · **Load example** · **Reset** |

**Load example** bundles a previously recorded run (BLAST XML, alignment, annotations) so the app is fully usable even when NCBI or EBI is slow or offline.

**Main area** — five tabs:

| Tab | Contents |
| --- | --- |
| **Overview** | Run status per stage, run metadata (parameters, accessions, timestamps, library versions), summary metrics, warnings |
| **BLAST** | Per-species candidate table, the selected hit highlighted with its selection reason, score charts, expandable raw alignments |
| **Alignment** | Residue-by-residue alignment viewer with consensus and identity rows, plus alignment statistics |
| **Conservation** | Conservation profile plot with domain bands and shaded conserved regions, regions table, feature-vs-rest comparison, pairwise identity heatmap |
| **Downloads** | Raw artefacts: BLAST XML, MSA FASTA, identity CSV |

A typical five-species run takes **about 5 minutes**, dominated by remote BLAST; results stream in as each stage completes. Every long step shows progress — the app never appears frozen.

---

## How to read the results

- **Percent identity** is computed per alignment column, with **gaps kept in the denominator**. With five sequences a gapped column can never reach 100%, so possible clean values are 20 / 40 / 60 / 80 / 100%. A column at 100% means every sequence has the same residue and no gaps.
- **Conserved regions** are runs of consecutive columns at or above the conservation threshold (default: strict 100%) that are at least 5 residues long. Both alignment and reference (human TP53) coordinates are reported.
- **Feature summary** is the key table: mean identity and fully-conserved share inside each annotated UniProt feature (domains, motifs, DNA-binding etc.) compared against the rest of the protein. For TP53, the DNA-binding core should stand out clearly — that contrast is the answer to the research question.
- All feature positions use **human TP53 residue numbering**. If you supply a bare FASTA instead of an accession, feature overlays are skipped (there is no annotation source) and the app says so.

---

## Scope and limitations

- Results are **displayed in the app only** — there is no report generator. The downloads are raw artefacts for your own write-up, not a formatted report.
- **Percent identity per column is the only conservation metric** (no phylogenetic trees, no 3D structure, no variant pathogenicity prediction).
- Runs depend on external services (NCBI BLAST/Entrez, UniProt, EBI) that are rate-limited and occasionally slow; the app retries with backoff, tolerates single-species failures, and offers the bundled example as a fallback.
- Single-user sessions, no accounts, no database, nothing persisted between sessions.
- Best-hit selection is automatic and deterministic; isoform ambiguity is handled by ranking rules and every candidate is shown so you can verify the choice.

## Documentation

- `Requirements & Design Document.md` — full specification: requirements, architecture, algorithms, acceptance criteria.
- `PLAN.md` — phased implementation plan for coding agents.
- `AGENTS.md` — working conventions for agents contributing to this repo.

## Development

The core logic lives in `src/tp53/` and is framework-free (no Streamlit imports outside `app/` and `streamlit_app.py`). Tests run fully offline against recorded fixtures in `tests/data/`:

```bash
python3.11 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
pytest
```

`src/tp53/config.py` holds every default, threshold and URL; `src/tp53/constants.py` holds the UI strings and colour palette. Nothing domain-specific is hard-coded in the logic.
