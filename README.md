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

## Running the Pipeline

### Prerequisites

```bash
# Install system dependencies
apt-get install fpocket apbs

# Install Python dependencies
pip install numpy biopython freesasa pdb2pqr matplotlib pytest
```

### Full Analysis

```bash
python pipeline/expanded_analysis.py
```

This runs all 9 structures through the complete pipeline and outputs:
- `results/expanded_validation_results.json` — Full results with all metrics
- `results/expanded_validation_summary.csv` — Summary table

### Tests

```bash
# Run all 62 tests
pytest tests/ -v

# Run only expanded tests (50 tests)
pytest tests/test_expanded.py -v
```

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
