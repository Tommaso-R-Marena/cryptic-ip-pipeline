# Methods Documentation

## Structure Acquisition

- Crystal structures downloaded from RCSB PDB (https://www.rcsb.org)
- AlphaFold models from EBI AlphaFold Database v6 (https://alphafold.ebi.ac.uk)
- Structures cleaned by stripping all HETATM records, water molecules, and non-protein atoms

## Pocket Detection (fpocket 3.1)

- Alpha-sphere based cavity detection
- Default parameters (min alpha sphere radius 3.0 Å, max 6.0 Å)
- Pockets ranked by druggability score
- Pocket-lining residues extracted from `_atm.pdb` output files

## Solvent Accessibility (FreeSASA)

- Lee–Richards algorithm with 1.4 Å probe radius
- Per-residue SASA computed for all atoms
- IP-binding residue SASA extracted by matching known residue numbers from crystal structures

## Charge Analysis

- Formal charges assigned: ARG +1, LYS +1, HIS +0.5, ASP −1, GLU −1
- Net charge computed within 8 Å radius of pocket center
- Basic residue count (ARG, LYS, HIS) within 5 Å

## Composite Scoring

Score = 0.30·norm(depth) + 0.35·(1−norm(SASA)) + 0.20·norm(charge) + 0.15·norm(basic)

Normalization bounds: depth/30Å, SASA/150Å², charge/15, basic/8

## Statistical Analysis

- Cohen's d effect size for group comparisons
- Mann-Whitney U test (non-parametric, appropriate for small n)
- All statistics computed with SciPy 1.x
