---
title: "Cryptic structural inositol-phosphate binding sites in protein structure"
author: "Tommaso Marena"
date: 2026
bibliography: references.bib
---

# Abstract

Buried structural inositol phosphates (IP3/IP4/IP5/IP6) are essential
cofactors of an unexpectedly diverse set of proteins, including ADAR2
[@macbeth2005], cohesin (Pds5B), class I HDACs, and the HIV-1 capsid
[@dick2018; @mallery2018]. Unlike surface PH-domain IP binders, these
sites are deeply buried and electropositive, often hidden in solved
structures because the IP molecule was omitted during crystallisation.

We present a reproducible, modular pipeline (`crypticip`) that integrates
fpocket [@leguilloux2009], FreeSASA [@mitternacht2016], and
APBS/pdb2pqr [@jurrus2018] to scan AlphaFold proteomes [@varadi2024]
for pockets that match a small set of biophysical criteria of the
ADAR2/IP6 site: depth > 15 Å, sidechain SASA < 5 Å², electrostatic
potential > 5 kT/e, ≥ 4 basic residues within 5 Å, volume 300–800 Å³,
mean pocket pLDDT ≥ 70.

# Methods

See `docs/methods.md`.

# Validation

(See `docs/validation_plan.md`.)

The published v3/v4 analyses in `results/expanded_validation_results.json`
report:

- ADAR2 IP6 pocket: rank 1 of 41.
- Mean pocket depth in positives vs negatives: 20.7 Å vs 10.5 Å.
- APBS potential at ADAR2/IP6 centre: +74.6 kT/e.
- Basic residues within 8 Å: 7.
- Mean pLDDT of ADAR2 AlphaFold model: 89.0.

# Screening (in progress)

(Plumbing in place; per-organism summaries land in
`results/screening/<organism>/`.)

# Limitations

- 4 positives + 4 negatives is too few to validate weighting at high
  confidence; the bootstrap AUC CI is wide.
- AlphaFold predicts the *apo* fold of IP-coordinating proteins;
  ADAR2's binding region RMSD between crystal (1ZY7) and AF
  (AF-P78563-F1) is large (~14.5 Å in the published v3 results)
  because IP6 binding induces conformational change. This is real, not
  a pipeline failure.
- APBS is computationally expensive; the screening "fast" mode skips
  it for the first pass.

# References
