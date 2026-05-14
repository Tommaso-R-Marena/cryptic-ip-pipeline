# Notebooks

Exploratory analyses live here. They are intentionally thin — the heavy
lifting belongs in `crypticip/` so it is testable. Start a notebook by
calling `crypticip.config.load_config(...)` and `crypticip.validation.validate_all(...)`.

## Google Colab notebooks

The five `*_colab.ipynb` notebooks below are **designed to be opened
directly from GitHub in a fresh Google Colab session and run
top-to-bottom**. Every notebook starts with the same *Run this first —
fresh-Colab setup* section that:

- exposes `REPO_URL`, `BRANCH`, `PROJECT_DIR`, `MOUNT_DRIVE`,
  `DRIVE_ROOT`, `DRIVE_RESULTS` as configurable variables,
- clones (or updates) the repo into `/content/` and `pip install -e .`,
- optionally mounts Drive and redirects `results/` onto it,
- verifies `crypticip --version`, `crypticip check-env`, and
  `crypticip list-configs`.

Re-running the bootstrap is idempotent. See `docs/colab.md` for full
caveats about storage and external binaries.

| File | Purpose | Open |
|------|---------|------|
| `00_colab_quickstart.ipynb` | Bootstrap, `check-env`, pytest, synthetic mini-proteome smoke workflow. | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/cryptic-ip-pipeline/blob/main/notebooks/00_colab_quickstart.ipynb) |
| `01_validation_colab.ipynb` | Optional condacolab install of binaries; `download-validation`; `validate --config config/validation.yaml`. | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/cryptic-ip-pipeline/blob/main/notebooks/01_validation_colab.ipynb) |
| `02_yeast_screening_colab.ipynb` | `download-proteome --resume` (smoke `--limit 50` first; full run commented), index, fast `--limit` screen, PyMOL, experimental plan. | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/cryptic-ip-pipeline/blob/main/notebooks/02_yeast_screening_colab.ipynb) |
| `03_results_analysis_colab.ipynb` | Load screening CSVs (Drive or on-the-fly smoke), rank, plot distributions / tiers / feature heatmap. | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/cryptic-ip-pipeline/blob/main/notebooks/03_results_analysis_colab.ipynb) |
| `04_experimental_prioritization_colab.ipynb` | Generate mutagenesis + DSF plans for top candidates. | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/cryptic-ip-pipeline/blob/main/notebooks/04_experimental_prioritization_colab.ipynb) |

The shared bootstrap logic also lives in
[`colab_bootstrap.py`](colab_bootstrap.py); notebooks **embed** inline
copies of it so they work without it, but `import colab_bootstrap` is
available after cloning if you prefer a one-liner.
