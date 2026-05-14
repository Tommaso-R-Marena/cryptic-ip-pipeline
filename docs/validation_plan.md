# Validation plan (Phase 1)

## Goal

Before running on any whole proteome we must show, with the curated
validation set, that the scoring scheme separates **buried IP sites**
(positive controls) from **surface PH-domain IP sites** (negative
controls), and in particular that the ADAR2/IP6 site is in the top 3
pockets of the ADAR2 crystal structure.

## Inputs

`config/validation.yaml` declares two groups:

- **Positive controls**: ADAR2 (1ZY7 + AF-P78563-F1), Pds5B (5HDT),
  HDAC1 (5ICN), HDAC3 (4A69).
- **Negative controls**: PLCδ1 PH (1MAI), Btk PH (1BTK), DAPP1 PH
  (1FAO), Grp1 PH (1FGY).

## Gate

`validation_gate:` in `config/validation.yaml`:

- `adar2_rank_ok` — best ADAR2 IP6 pocket ranked ≤ 3.
- `separation_ok` — mean positive composite score − mean negative ≥ 0.05.
- `depth_separation_ok` — mean positive pocket depth − mean negative ≥ 5 Å.
- `positives_passing_ok` — ≥ 2 positives pass both depth and basic flags.

All four must pass for `gate.overall_pass = True`.

## Outputs

`crypticip validate` writes (under `results/reports/validation/`):

- `validation_results.json` — full per-structure dump.
- `validation_summary.csv` — one row per structure with the key fields.
- `validation_gate_report.md` — pass/fail summary, env, AUC, bootstrap CI.
- `validation_report.md` / `.html` — rendered report.
- `scores.png` — composite-score bar chart by structure.

## Pre-existing scientific record

The original repository already shipped, in `results/`:

- `validation_results.json` and `validation_summary.csv` from the
  initial 9-structure run with real fpocket/FreeSASA/APBS values.
- `expanded_validation_results.json` and `expanded_validation_summary.csv`
  from the v3 expanded analysis.

The new pipeline **does not overwrite** these files. Re-running
`crypticip validate` writes into `results/reports/validation/` (a
different directory) so the historical record remains the canonical
scientific result. Any files that would clash are timestamp-backed-up.

## Failure handling

If the gate fails:

1. Inspect the `validation_results.json` for each failing positive.
2. Look at `flags` for each structure — which criterion is missing?
3. If multiple positives fail because APBS is missing, run a full APBS
   sweep and re-run validate; APBS-failed pockets cannot reach Tier 1.
4. Try `--config config/scoring.yaml` with `no_apbs` or `depth_heavy`
   weighting; record the result, but **do not** promote a weighting
   that was hand-tuned to the gate without running the leave-one-out
   AUC report (`crypticip.statistics.leave_one_out_auc`).
