# Scoring

The composite score lives in `crypticip.scoring.composite_score`.

## Feature vector

`FeatureVector(depth, sasa, elec, basic_count, volume, plddt, apbs_status)`:

| field         | unit  | source                                                |
|---------------|-------|-------------------------------------------------------|
| `depth`       | Å     | fpocket `Cent. of mass − Alpha Sphere max dist`       |
| `sasa`        | Å²    | mean sidechain SASA of pocket-lining residues         |
| `elec`        | kT/e  | APBS potential at pocket centre (Coulomb fallback)    |
| `basic_count` | int   | basic residues (R/K/H) within 5 Å of centre           |
| `volume`      | Å³    | fpocket monte-carlo / real volume                     |
| `plddt`       | 0-100 | mean per-residue B-factor of AlphaFold structures     |
| `apbs_status` | enum  | `ok` / `fallback` / `failed` / `missing`              |

## Normalisations

```yaml
norms:
  depth_A:     30.0    # depth >= 30 saturates the component
  sasa_A2:     150.0   # SASA == 0 → inv_sasa = 1.0
  elec_kT_e:   30.0    # potential >= +30 saturates the elec component
  basic_count: 8.0
  volume_lo:   300.0   # below: scales linearly from 0 → 1
  volume_hi:   800.0   # above: linear decay until 2x hi
```

## Weights

`config/scoring.yaml` ships four named schemes:

- `default` — 25/25/20/20/10 across depth / inv_sasa / elec / basic / volume_fit.
- `depth_heavy` — emphasises depth (40%).
- `electrostatics_heavy` — emphasises APBS (35%).
- `no_apbs` — zero APBS weight (for re-scoring when APBS is unavailable).

## Tiering

See `docs/methods.md`. Tier 1 requires all six filters to pass and
`apbs_status == "ok"`, so a pocket whose only signal is the Coulomb
fallback cannot reach Tier 1 — it is at best Tier 2.

## Sensitivity / overfitting

`crypticip.statistics.weight_sensitivity` runs a Dirichlet sweep over
the five weights and reports the AUC distribution. Use
`bootstrap_auc(scores, labels)` for a 95 % CI and `permutation_p(...)`
to test AUC > 0.5. With only 4 positives + 4 negatives the bounds are
wide; treat any AUC ≥ 0.9 as encouraging but not definitive.
