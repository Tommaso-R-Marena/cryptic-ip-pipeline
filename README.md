# Cryptic IP Binding Site Detection Pipeline

A computational pipeline for identifying cryptic inositol phosphate (IP) binding sites in proteins using pocket geometry, electrostatics, and structural comparison.

## Overview

This pipeline integrates multiple bioinformatics tools to characterize potential IP binding pockets in protein structures, discriminating deeply buried cryptic sites from surface-exposed PH-domain pockets.

### Pipeline Components

| Tool | Purpose | Version |
|------|---------|---------|
| **fpocket** | Pocket detection, volume, depth | 3.1.4.2 |
| **FreeSASA** | Solvent-accessible surface area | Python bindings |
| **PDB2PQR** | PDB to PQR conversion (charge assignment) | 3.7.1 |
| **APBS** | Adaptive Poisson-Boltzmann electrostatics | 3.4.1 |
| **BioPython** | Structural alignment and RMSD | Latest |

### Metrics Computed

1. **Pocket detection** — fpocket rank, score, druggability
2. **Pocket geometry** — Volume (Å³), depth (Å, from "Cent. of mass - Alpha Sphere max dist")
3. **Solvent accessibility** — Per-residue sidechain SASA via FreeSASA
4. **Electrostatics** — APBS potential at pocket center (kT/e)
5. **Charge density** — Basic/acidic residue counts at 5/8/10 Å using sidechain terminal atoms
6. **pLDDT filtering** — AlphaFold confidence scoring (≥70 threshold)
7. **RMSD comparison** — Crystal vs AlphaFold binding region (SVD superposition)
8. **Composite scoring** — Weighted multi-metric score with volume and APBS components

## Validation Set

| Structure | PDB | Category | IP Type | Description |
|-----------|-----|----------|---------|-------------|
| ADAR2 crystal | 1ZY7 | Positive | IP6 | Buried IP6 site (known coordinating residues) |
| ADAR2 AlphaFold | AF-P78563 | Positive | IP6 | Predicted structure without IP6 |
| HDAC1 | 5ICN | Positive | IP4 | Buried IP4 site |
| HDAC3 | 4A69 | Positive | IP4 | Buried IP4 site |
| Pds5B | 5HDT | Positive | IP6 | Buried IP6 site |
| PLCδ1 PH | 1MAI | Negative | IP3 | Surface PH domain |
| Btk PH | 1BTK | Negative | IP4 | Surface PH domain |
| DAPP1 PH | 1FAO | Negative | IP4 | Surface PH domain |
| Grp1 PH | 1FGY | Negative | IP4 | Surface PH domain |

## Key Results

### Success Criteria (Document Sections 5-6)

| Criterion | Target | Result | Status |
|-----------|--------|--------|--------|
| ADAR2 IP6 pocket rank | Top 3 | #1 of 41 pockets | **PASS** |
| Electrostatic potential | >5 kT/e | +74.6 kT/e | **PASS** |
| Basic residues near site | ≥6 within 8 Å | 7 basic residues | **PASS** |
| AlphaFold pLDDT | ≥70 average | 89.0 | **PASS** |
| Pocket depth separation | Positives deeper | 20.7 vs 10.5 Å | **PASS** |
| AF vs crystal RMSD | <2 Å | 14.56 Å | **FAIL** (real) |

The large RMSD (14.56 Å) between crystal and AlphaFold binding regions is a genuine finding: AlphaFold, without the IP6 ligand as a structural template, predicts a substantially different conformation for the binding site. This confirms that IP6 binding induces significant conformational change.

## Repository Structure

```
cryptic-ip-pipeline/
├── pipeline/
│   ├── expanded_analysis.py   # Full multi-metric analysis module
│   ├── full_analysis.py       # Original analysis module
│   └── generate_figures.py    # Figure generation
├── data/
│   ├── pdb/                   # Crystal structures (PDB format)
│   ├── alphafold/             # AlphaFold predictions
│   └── fpocket_results/       # fpocket output directories
├── results/
│   ├── expanded_validation_results.json
│   ├── expanded_validation_summary.csv
│   ├── validation_results.json
│   └── validation_summary.csv
├── figures/                   # Publication-quality figures (8 total)
├── tests/
│   ├── test_expanded.py       # 50 tests covering all new features
│   └── test_scoring.py        # 12 original scoring tests
└── README.md
```

## Running the Pipeline (v0.5+ — `crypticip` CLI)

The original validation scripts under `pipeline/` and the validated
results in `results/` are kept as-is for the scientific record. The
production pipeline is now a Python package `crypticip` with a single
CLI of the same name.

### Install

```bash
# Conda (recommended — also installs fpocket / FreeSASA / APBS / pdb2pqr)
conda env create -f environment.yml
conda activate crypticip

# Pip-only (assumes binaries are installed separately)
pip install -e .[dev]
```

### Quick start

```bash
crypticip --version
crypticip check-env                     # verify external tools
crypticip download-validation           # fetch crystal + AlphaFold structures
crypticip validate --config config/validation.yaml
crypticip report --validation

# Proteome screening
crypticip download-proteome --organism yeast --resume
crypticip index-proteome --organism yeast
crypticip screen --organism yeast --mode fast --workers 8 --resume
crypticip report --organism yeast
crypticip screen --organism yeast --mode full --workers 8 --resume
crypticip pymol --organism yeast --top 50
crypticip experimental-plan --organism yeast --top 25
# repeat for --organism human, --organism dictyostelium, then:
crypticip report --all
crypticip experimental-plan --all --top 25
```

### Legacy entry point

The original v3/v4 pipeline still works for re-running the published
validation:

```bash
python pipeline/expanded_analysis.py
```

### Tests

```bash
python -m pytest -q              # ~ 90 tests, no external binaries required
```

See `IMPLEMENTATION_REPORT.md` for the v0.5 refactor summary, `docs/`
for the per-stage documentation, and `IMPLEMENTATION_PLAN.md` for the
checkpoint structure used to build it.

### Run on Google Colab

The `notebooks/` directory contains five Colab-ready notebooks, each
designed to be opened directly from GitHub in a fresh Colab session
and run top-to-bottom. Every notebook starts with a *Run this first —
fresh-Colab setup* section that clones the repo, `pip install -e .`s
it, optionally mounts Drive, and verifies `crypticip --version`,
`crypticip check-env`, and `crypticip list-configs`:

| Notebook | Purpose | Open |
|----------|---------|------|
| `00_colab_quickstart.ipynb` | Clone, install, `check-env`, tests, synthetic smoke workflow | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/cryptic-ip-pipeline/blob/main/notebooks/00_colab_quickstart.ipynb) |
| `01_validation_colab.ipynb` | Run validation on the curated control set (with optional conda install of external binaries) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/cryptic-ip-pipeline/blob/main/notebooks/01_validation_colab.ipynb) |
| `02_yeast_screening_colab.ipynb` | Yeast proteome download / index / screen / PyMOL / experimental plan | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/cryptic-ip-pipeline/blob/main/notebooks/02_yeast_screening_colab.ipynb) |
| `03_results_analysis_colab.ipynb` | Load screening CSVs and plot distributions / tiers / feature heatmap | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/cryptic-ip-pipeline/blob/main/notebooks/03_results_analysis_colab.ipynb) |
| `04_experimental_prioritization_colab.ipynb` | Generate mutagenesis + DSF plans for top candidates | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/cryptic-ip-pipeline/blob/main/notebooks/04_experimental_prioritization_colab.ipynb) |

The external scientific binaries (fpocket / FreeSASA / APBS / PyMOL) are
not available on a vanilla Colab runtime; the notebooks document the
fallback behaviour and a `condacolab + mamba` install path. See
[`docs/colab.md`](docs/colab.md) for the full guide.

## Technical Notes

### Chain-Aware Analysis
For multimeric structures (e.g., 1ZY7 homodimer), the pipeline uses only the first chain to avoid averaging coordinates across monomers.

### Sidechain Atom Distances
Charged residue distances are measured from sidechain terminal atoms (NZ for Lys, NH1 for Arg, NE2 for His) rather than Cα atoms, providing more accurate proximity measurements to the IP molecule.

### APBS DX Parser
The DX file parser handles non-numeric trailing lines ("attribute", "component") that APBS appends after the data block.

## Author

Tommaso Marena, Department of Chemistry, The Catholic University of America

## Acknowledgments

Dr. Gregory Miller, Department of Chemistry, CUA — Project supervision and guidance.
