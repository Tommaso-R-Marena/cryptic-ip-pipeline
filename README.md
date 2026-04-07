# Cryptic Buried Inositol Phosphate Binding Site Discovery Pipeline

**Computational identification of buried, structurally required inositol phosphate (IP) binding sites across proteomes.**

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

## Overview

Most known inositol phosphate (IP) binding proteins use surface-exposed PH, FYVE, or C2 domains for reversible signaling interactions. However, a small but growing number of proteins — beginning with ADAR2 (Macbeth et al., 2005) — require a **completely buried IP cofactor** for structural integrity and catalytic function.

This pipeline systematically identifies buried IP-binding pockets across proteomes using:

| Tool | Purpose |
|------|---------|
| **fpocket** | Alpha-sphere pocket detection |
| **FreeSASA** | Per-residue solvent-accessible surface area |
| **Charge analysis** | Local electrostatic environment profiling |
| **Composite scoring** | Weighted burial score (depth 30%, SASA 35%, electrostatics 20%, basic residues 15%) |

## Key Results

### Phase 1 Validation (9 structures, 5 positive + 4 negative controls)

| Metric | Positive Controls | Negative Controls | Cohen's d | p-value |
|--------|-------------------|-------------------|-----------|---------|
| Composite Score | 0.56 ± 0.09 | 0.42 ± 0.08 | +1.62 | 0.111 |
| Pocket Depth (Å) | 40.1 ± 16.2 | 9.4 ± 5.7 | **+2.52** | **0.016** |
| IP-Site SASA (Å²) | 57.8 ± 35.8 | 33.6 ± 16.7 | +0.87 | 0.393 |
| Net Charge (8 Å) | +2.3 ± 2.9 | +2.2 ± 2.0 | +0.02 | 1.000 |

**Pocket depth is the strongest single discriminator** (Cohen's d = 2.52, Mann-Whitney p = 0.016), cleanly separating buried IP sites (mean 40 Å) from surface PH-domain sites (mean 9 Å).

### Validation Structures

**Positive controls (buried IP):**
- ADAR2 deaminase — IP6 (PDB 1ZY7 + AlphaFold AF-P78563-F1)
- HDAC1 deacetylase — IP4 (PDB 5ICN)
- HDAC3 deacetylase — IP4 (PDB 4A69)
- Pds5B cohesin regulator — IP6 (PDB 5HDT)

**Negative controls (surface PH domains):**
- PLCδ1 (PDB 1MAI), Btk (PDB 1BTK), DAPP1 (PDB 1FAO), Grp1 (PDB 1FGY)

## Installation

```bash
git clone https://github.com/Tommaso-R-Marena/cryptic-ip-pipeline.git
cd cryptic-ip-pipeline
pip install -r requirements.txt
```

### External dependencies

- **fpocket 3.x** — build from source ([Discngine/fpocket](https://github.com/Discngine/fpocket))
- **FreeSASA** — installed via `pip install freesasa`

## Usage

### Run the full validation pipeline

```bash
# 1. Prepare structures (strip HETATM, clean chains)
python pipeline/prepare_structures.py .

# 2. Run fpocket on all structures
for pdb in data/cleaned/*.pdb; do
    fpocket -f "$pdb"
    mv "${pdb%.pdb}_out" data/fpocket_results/
done

# 3. Run full analysis (SASA + charge + scoring)
python pipeline/full_analysis.py

# 4. Generate figures
python pipeline/generate_figures.py
```

### Run tests

```bash
python -m pytest tests/ -v
```

## Repository Structure

```
cryptic-ip-pipeline/
├── README.md
├── LICENSE
├── requirements.txt
├── pipeline/
│   ├── prepare_structures.py    # PDB cleaning and chain extraction
│   ├── full_analysis.py         # Main analysis: fpocket + SASA + scoring
│   ├── generate_figures.py      # Publication-quality figure generation
│   └── analyze_all.py           # Legacy analysis module
├── data/
│   ├── pdb/                     # Raw PDB crystal structures
│   ├── alphafold/               # AlphaFold structural models
│   ├── cleaned/                 # Processed protein-only PDBs
│   ├── fpocket_results/         # fpocket output directories
│   └── sasa_results/            # FreeSASA output
├── results/
│   ├── validation_results.json  # Full results with all metrics
│   └── validation_summary.csv   # Summary table
├── figures/
│   ├── fig1_composite_scores.png
│   ├── fig2_multi_panel.png
│   ├── fig3_score_breakdown.png
│   ├── fig4_adar2_sasa.png
│   ├── fig5_scatter_depth_sasa.png
│   ├── fig6_pipeline.png
│   └── fig7_fpocket_stats.png
├── tests/
│   └── test_scoring.py
├── docs/
│   └── methods.md
└── validation/
    └── controls.yaml
```

## Composite Scoring Function

```
S = 0.30·norm(depth) + 0.35·(1 − norm(SASA)) + 0.20·norm(charge) + 0.15·norm(basic)
```

| Component | Weight | Normalization | Rationale |
|-----------|--------|---------------|-----------|
| Pocket Depth | 30% | Cap at 30 Å | Burial distance from protein surface |
| Inverse SASA | 35% | Cap at 150 Å², inverted | Low SASA = deeply buried residues |
| Charge Density | 20% | Cap at +15 formal charges | IP coordination requires positive charges |
| Basic Residues | 15% | Cap at 8 within 5 Å | Arg/Lys/His direct IP coordination |

## Proteome-Scale Screening (Planned)

| Proteome | AlphaFold Models | [IP6] (µM) | Status |
|----------|------------------|------------|--------|
| *S. cerevisiae* | 6,049 | 15–25 | Planned |
| *H. sapiens* | 23,391 | 15–50 | Planned |
| *D. discoideum* | 12,622 | ~520 | Planned |

*D. discoideum* has 10–30× higher intracellular IP6 than mammalian cells, enabling a co-evolutionary test: if buried-site frequency tracks [IP6], protein architecture co-evolves with IP metabolism.

## References

1. Macbeth MR et al. (2005) Inositol hexakisphosphate is bound in the ADAR2 core and required for RNA editing. *Science* 309:1534–1539.
2. Le Guilloux V et al. (2009) Fpocket: An open source platform for ligand pocket detection. *BMC Bioinformatics* 10:168.
3. Jumper J et al. (2021) Highly accurate protein structure prediction with AlphaFold. *Nature* 596:583–589.
4. Blind RD (2020) Structural analyses of inositol phosphate second messengers. *Adv Biol Regul* 75:100667.
5. Dick RA et al. (2018) Inositol phosphates are assembly co-factors for HIV-1. *Nature* 560:509–512.
6. Yuan L et al. (2022) IP6 binds an allosteric site on Uba6. *Nat Commun* 13:4871.

## Author

**Tommaso Marena**  
Department of Chemistry, The Catholic University of America  
Advisor: Dr. Gregory Miller

## License

MIT License — see [LICENSE](LICENSE) for details.
