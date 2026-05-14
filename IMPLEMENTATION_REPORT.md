# Implementation Report — v0.5 production refactor

## Summary

The repository now ships, alongside the original v3/v4 scripts and the
historical scientific record, a production Python package `crypticip`
with a single CLI of the same name. The new package gives:

- a packaged, importable, unit-tested library covering every pipeline
  stage (preprocessing, fpocket, FreeSASA, APBS, residue features,
  composite scoring, validation, screening, reporting, experimental
  planning),
- a single CLI `crypticip` with subcommands `check-env`,
  `download-validation`, `validate`, `download-proteome`,
  `index-proteome`, `screen`, `report`, `pymol`, `experimental-plan`,
  `list-configs`,
- mocks / fallbacks for every external binary so the pipeline runs
  (with status flags) when fpocket / FreeSASA / APBS / pdb2pqr / PyMOL
  are absent,
- proteome ingestion + screening orchestration plumbed end-to-end,
  including a parallel process pool and resume-safe per-protein JSON
  output,
- a smoke-test CI workflow that exercises the screening orchestration
  with a synthetic mini-proteome (no external binaries needed),
- docs covering methods, scoring, validation, ingestion, screening,
  troubleshooting, reproducibility, and the public API,
- a manuscript skeleton + BibTeX with the core references requested.

The pre-existing v3/v4 scripts and their committed results
(`results/validation_results.json`, `results/expanded_validation_*.{json,csv}`,
all figures, the poster) are **untouched**. The new CLI writes into
`results/reports/validation/` and `results/screening/<organism>/`, so
both records co-exist.

## Files changed

80 files changed, ~5400 insertions, ~70 deletions. Highlights:

```
crypticip/                      new package (15 modules)
  __init__.py
  paths.py logging_utils.py config.py external_tools.py
  pdb_io.py structures.py alphafold.py
  fpocket.py sasa.py electrostatics.py residues.py scoring.py
  qc.py statistics.py
  validation.py screening.py reporting.py pymol.py experimental.py
  conservation.py literature.py
  cli.py

config/                         new
  default.yaml validation.yaml yeast.yaml human.yaml dictyostelium.yaml scoring.yaml

scripts/                        thin CLI wrappers
  download_validation_structures.py download_alphafold_proteome.py
  build_proteome_index.py run_validation.py run_screen.py
  generate_reports.py make_pymol_sessions.py

tests/                          new modules + fixtures
  conftest.py
  fixtures/{tiny.pdb,fpocket_info.txt,potential.dx}
  test_pdb_io.py test_structures.py test_alphafold.py
  test_fpocket.py test_sasa.py test_electrostatics.py
  test_residues.py test_scoring_v2.py test_statistics.py
  test_config.py test_cli.py test_qc_and_screening.py
  test_pymol_experimental.py

docs/                           new (methods.md replaces v1 — preserved as methods_v1.md)
  methods.md validation_plan.md scoring.md
  proteome_ingestion.md screening.md experimental_followup.md
  troubleshooting.md reproducibility.md api.md

manuscript/                     new
  paper.md supplement.md references.bib

.github/workflows/              tests.yml + validation-smoke.yml (replaces ci.yml)

Top-level                       pyproject.toml, environment.yml, Dockerfile,
                                docker-compose.yml, Makefile, CITATION.cff,
                                IMPLEMENTATION_PLAN.md, IMPLEMENTATION_REPORT.md,
                                updated .gitignore + README
```

## Tests run / results

```
$ python -m pytest -q
90 passed, 25 skipped in 0.36s
```

- 90 unit tests pass with **no** external bioinformatics binaries
  installed (fpocket, FreeSASA, APBS, pdb2pqr, PyMOL absent).
- 25 of the pre-existing `tests/test_expanded.py` tests skip
  gracefully when the corresponding raw data files are not yet
  downloaded — this is the original behaviour.
- A live end-to-end smoke test ran `index-proteome → screen →
  report → experimental-plan` against a synthetic 3-protein
  proteome (made from `tests/fixtures/tiny.pdb`); all four
  subcommands returned 0 and wrote the expected output files. The
  screening produced 0 candidates because fpocket is missing from
  the sandbox; that's the documented degraded mode.

## Scientific judgement / limitations

- **Do not** treat the new CLI's `validate` output as scientifically
  authoritative until it has been re-run with real fpocket, FreeSASA,
  pdb2pqr, and APBS installed. The historical authoritative result
  lives in `results/validation_results.json` and
  `results/expanded_validation_results.json` (v3 run with all four
  tools). The CLI deliberately writes into a separate directory.
- Mass screening of yeast / human / dictyostelium is **plumbed but
  not executed** here. Running on the full human proteome (~23k
  models) takes ~24 h on 8 cores in "full" mode and is intentionally
  left to the user.
- With only 4 positive + 4 negative validation structures, weight
  optimisation is brittle. The `weight_sensitivity` / bootstrap /
  permutation helpers report wide confidence intervals; treat any
  high AUC as encouraging only, not conclusive.
- APBS-failed pockets cannot be promoted to Tier 1 even if all other
  filters pass, by construction in `crypticip.scoring.tier`. If a
  candidate looks compelling in the per-residue dump but is APBS-missing,
  it should be re-run with APBS available before any wet-lab investment.
- The conservation module is a placeholder; full ConSurf/Rate4Site
  integration would be a follow-up.

## Branch / commit info

- Working branch: **`feat/production-pipeline`** (created from main).
- Initial commit: pending (will be created once the user approves the
  diff). Run:
  ```bash
  git add -A
  git commit -m "feat: v0.5 production refactor — crypticip package + CLI"
  git push -u origin feat/production-pipeline
  gh pr create
  ```

## Exact commands for the user

```bash
# 1. install (one of)
conda env create -f environment.yml && conda activate crypticip
# or
pip install -e .[dev]

# 2. environment check
crypticip check-env

# 3. tests
python -m pytest -q

# 4. validation (requires real binaries; writes into results/reports/validation/)
crypticip download-validation
crypticip validate --config config/validation.yaml
crypticip report --validation

# 5. proteome screening (yeast)
crypticip download-proteome --organism yeast --resume
crypticip index-proteome --organism yeast
crypticip screen --organism yeast --mode fast --workers 8 --resume
crypticip report --organism yeast
crypticip screen --organism yeast --mode full --workers 8 --resume
crypticip report --organism yeast

# 6. PyMOL + experimental plan
crypticip pymol --organism yeast --top 50
crypticip experimental-plan --organism yeast --top 25

# 7. repeat for human / dictyostelium, then
crypticip report --all
crypticip experimental-plan --all --top 25
```
