# Reproducibility

## What is captured

Every CLI command stamps run metadata into its output (under
`metadata` in the validation report, and under `summary.json` for
screening):

- crypticip version
- argv used to launch the command
- host name, Python version, platform
- start time (UTC)
- git commit (resolved with `git rev-parse HEAD`)
- config hash (`Config.hash()` — first 12 hex chars of a SHA-256 over
  the merged YAML)
- per-tool versions from `crypticip.external_tools.env_report`

## Deterministic seeds

Statistical helpers in `crypticip.statistics` (`weight_sensitivity`,
`bootstrap_auc`, `permutation_p`) take a `seed=` argument; the CLI does
not currently expose them but downstream notebooks should pass the
same seed across runs.

## Environment

Two complementary install routes:

- `conda env create -f environment.yml && conda activate crypticip`
  pulls fpocket / freesasa / APBS / pdb2pqr from bioconda.
- `pip install -e .[dev]` installs only the Python deps. External
  binaries must be on PATH for the relevant stages to run; otherwise
  the pipeline degrades to documented fallbacks.

A `Dockerfile` builds an image with the conda environment baked in.

## CI

- `.github/workflows/tests.yml` runs unit tests on 3.10/3.11/3.12.
- `.github/workflows/validation-smoke.yml` runs an in-repo synthetic
  proteome through `screen_proteome` to ensure the orchestration path
  doesn't regress.

CI does not run with external bioinformatics binaries installed; that
would couple every push to bioconda's availability. End-to-end runs
with real binaries belong on a server / cluster.

## What is *not* reproducible

- AlphaFold model versions move forward (v4 today). Pin the URL in
  `config/<organism>.yaml` if you need long-term reproducibility.
- The pre-existing scientific record in `results/validation_results.json`
  was produced by the original `pipeline/expanded_analysis.py`; the new
  CLI writes into `results/reports/validation/` so the two are kept
  side-by-side.
