# Methods

(Supersedes [`methods_v1.md`](methods_v1.md); the original document is kept verbatim alongside.)

## Structure acquisition

- Crystal structures: RCSB PDB (https://www.rcsb.org).
- AlphaFold models: AlphaFold DB v4 (https://alphafold.ebi.ac.uk).
- Cleaning:
  - **AlphaFold** structures: all HETATM records are removed
    (`crypticip.pdb_io.clean_for_alphafold`).
  - **Crystals**: waters, ions and buffer molecules are dropped but
    IP-family ligands (IHP, I3P, I4P, IP6, …) are preserved so the IP
    site can be located.
- Per-chain analysis: ``chain`` argument is honoured throughout so
  multimers don't smear coordinates across monomers.

## Pocket detection (fpocket)

- Binary: fpocket >= 3.1 (`tools.fpocket.binary` in config).
- The ``*_info.txt`` file is parsed in `crypticip.fpocket.parse_info_text`;
  the parser is robust to ordering changes, unit annotations, and to
  keys that share substrings with each other (e.g. "Score" vs "Local
  hydrophobic density score").
- Per-pocket residue lists and centroids come from
  `pocket{N}_atm.pdb`.
- Pocket depth = "Cent. of mass − Alpha Sphere max dist".

## Solvent accessibility (FreeSASA)

- Python bindings preferred; binary fallback supported via PATH.
- Lee–Richards, probe radius 1.4 Å (configurable).
- Per-residue total / sidechain / mainchain SASA;
  `crypticip.sasa.mean_sidechain_sasa(result, residues)` reports the
  mean sidechain SASA over a list of residues.
- If FreeSASA is missing the pipeline records `sasa_status="missing"`
  and continues; the SASA filter is treated as informational.

## Electrostatics (pdb2pqr + APBS)

- pdb2pqr at pH 7.0 with AMBER force field.
- APBS `mg-auto` with finer focusing.
- DX grid parser (`crypticip.electrostatics.parse_dx`) handles trailing
  `attribute` / `component` lines APBS appends after the data block.
- Potential reported at the alpha-sphere centroid of each pocket.
- If APBS is unavailable the pipeline falls back to a Coulomb-like
  estimate over formal side-chain charges within 15 Å of the centre and
  tags the result `apbs_status="fallback"`. Such pockets cannot be
  promoted to Tier 1.

## Charge / residue analysis

- Formal charges: ARG +1, LYS +1, HIS +0.5, ASP −1, GLU −1.
- Distances are measured from sidechain terminal atoms (NZ, CZ/NH, NE2,
  CG/OD, CD/OE, …) rather than from Cα — see
  `crypticip.residues.SIDECHAIN_TERMINAL`.
- Multi-shell summaries at 5/8/10 Å are produced for every pocket and
  surface in the per-protein JSON output.

## Composite scoring

Default (configurable in `config/scoring.yaml`):

```
score = 0.25·norm(depth) + 0.25·(1−norm(SASA)) + 0.20·norm(elec)
      + 0.20·norm(basic) + 0.10·volume_fit − plddt_penalty
```

Normalisations: depth/30 Å, SASA/150 Å², elec/30 kT/e, basic/8, volume_fit
is 1 inside [300, 800] Å³ and decays linearly outside. The pLDDT penalty
(default 0.2) applies if mean pLDDT < 70.

## Tiering

`crypticip.scoring.tier`:
- **Tier 1** — all filter flags pass AND score ≥ 0.65 AND APBS status = ok.
- **Tier 2** — all-but-one filter pass AND score ≥ 0.55.
- **Tier 3** — score ≥ 0.40.
- **Reject** — otherwise.

## Statistical evaluation

- ROC AUC (`crypticip.statistics.roc_auc`).
- Bootstrap AUC 95 % CI.
- Permutation test for AUC > 0.5.
- Cohen's d + Mann–Whitney U for positive vs negative score
  distributions (small-N caveats noted).
- Weight sensitivity: Dirichlet random sweep over the five normalised
  components, AUC and pos–neg separation reported per draw.
