# Running on Google Colab

The five notebooks under `notebooks/` are **designed to be opened
directly from GitHub in a fresh Google Colab runtime and run
top-to-bottom**. Each starts with a *Run this first — fresh-Colab
setup* section that:

- exposes `REPO_URL`, `BRANCH`, `PROJECT_DIR`, `MOUNT_DRIVE`,
  `DRIVE_RESULTS` as configurable variables in the first code cell,
- clones (or updates) the repo into `/content/`,
- `pip install -e .` (idempotent),
- optionally mounts Google Drive and redirects `results/` onto it,
- verifies `crypticip --version`, `crypticip check-env`, and
  `crypticip list-configs`.

Re-running any bootstrap cell is safe: the clone is `fetch + reset`ed
in-place and the pip install is editable. After a `condacolab` kernel
restart you must re-run the bootstrap cells.

The `crypticip` CLI itself runs fine in a stock Colab Python image.
**External scientific binaries** (`fpocket`, `freesasa`, `pdb2pqr`,
`apbs`, `pymol`) do not. Each notebook with a real scientific step has
an opt-in *Install external scientific binaries* section that uses
`condacolab + mamba`. Without it, the CLI runs in clearly-labelled
fallback mode — useful for plumbing checks, **not** for scientific
conclusions.

## Notebooks

| Order | Notebook | What it does | Time | Disk | Open |
|-------|----------|--------------|------|------|------|
| 1 | [`00_colab_quickstart.ipynb`](../notebooks/00_colab_quickstart.ipynb) | Bootstrap, `check-env`, `pytest`, synthetic mini-proteome `index → screen → report → pymol → experimental-plan` | ~5 min | < 100 MB | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/cryptic-ip-pipeline/blob/main/notebooks/00_colab_quickstart.ipynb) |
| 2 | [`01_validation_colab.ipynb`](../notebooks/01_validation_colab.ipynb) | Optional condacolab install of binaries; `download-validation`; `validate --config config/validation.yaml`; ADAR2 diagnostics | ~10-15 min real / ~3 min fallback | ~50 MB | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/cryptic-ip-pipeline/blob/main/notebooks/01_validation_colab.ipynb) |
| 3 | [`02_yeast_screening_colab.ipynb`](../notebooks/02_yeast_screening_colab.ipynb) | Yeast AlphaFold proteome `download-proteome --resume` (smoke `--limit 50` first; full run commented out), index, fast `--limit` screen, PyMOL bundles, experimental plan | 30 min – several h | up to ~6 GB | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/cryptic-ip-pipeline/blob/main/notebooks/02_yeast_screening_colab.ipynb) |
| 4 | [`03_results_analysis_colab.ipynb`](../notebooks/03_results_analysis_colab.ipynb) | Load screening CSVs (from Drive or via on-the-fly smoke screen), rank, plot distributions / tiers / feature heatmap | minutes | < 10 MB | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/cryptic-ip-pipeline/blob/main/notebooks/03_results_analysis_colab.ipynb) |
| 5 | [`04_experimental_prioritization_colab.ipynb`](../notebooks/04_experimental_prioritization_colab.ipynb) | Generate mutagenesis + DSF plans for top candidates (same Drive-or-smoke pattern as 03) | minutes | < 10 MB | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/cryptic-ip-pipeline/blob/main/notebooks/04_experimental_prioritization_colab.ipynb) |

Open each by clicking the **Open in Colab** badge on GitHub, or directly with:

```
https://colab.research.google.com/github/Tommaso-R-Marena/cryptic-ip-pipeline/blob/main/notebooks/00_colab_quickstart.ipynb
```

(replace the filename for each notebook).

## Bootstrap helper

`notebooks/colab_bootstrap.py` is an optional helper that wraps the
clone / install / Drive-mount / CLI-verify steps. Notebooks **embed**
inline copies of these steps so they work without it — but after
cloning you can also do:

```python
import sys; sys.path.insert(0, '/content/cryptic-ip-pipeline/notebooks')
import colab_bootstrap as cb
cb.bootstrap(branch='main', mount_drive_flag=True,
             drive_results='/content/drive/MyDrive/crypticip/results')
```

## External binaries on Colab

The only reliable Colab path for `fpocket` + `apbs` + `pdb2pqr` +
`pymol-open-source` + `freesasa` is **condacolab + mamba**:

```python
!pip install -q condacolab
import condacolab; condacolab.install()       # restarts the kernel
# after restart:
!mamba install -y -c conda-forge -c bioconda fpocket freesasa pdb2pqr apbs pymol-open-source python-freesasa
```

Caveats:
- `condacolab.install()` **restarts the kernel** — you must re-run the
  fresh-Colab bootstrap cells (clone + `pip install -e .`) after the
  restart, because Colab loses cwd and the editable install.
- The install takes 5-10 minutes and ~3 GB of disk.
- PyMOL's `apt-get` package on Colab pulls in a heavy GUI stack;
  prefer `pymol-open-source` from conda-forge.
- For reproducible runs, use the project `Dockerfile` instead of Colab
  (see `docs/reproducibility.md`).

## Without external binaries (fallback mode)

If you skip the condacolab section, the CLI still runs but produces:
- `fpocket_status = missing` → zero detected pockets, empty screening CSV.
- `freesasa missing` → `sasa_status = missing` on every record.
- `apbs / pdb2pqr missing` → `apbs_status = fallback` (Coulomb potential).
- `pymol missing` → `.pml` sessions are still written; PNG renders are
  skipped.

Treat fallback-mode numbers as plumbing diagnostics, **not** scientific
findings.

## Persisting results across runtimes

Colab runtimes are ephemeral. To keep outputs, set `MOUNT_DRIVE=True`
and `DRIVE_RESULTS=/content/drive/MyDrive/...` in the bootstrap cell.
The notebook will then symlink `<project>/results` onto Drive (any
pre-existing local `results/` is backed up to a timestamped sibling,
not deleted). Alternatively, every notebook ends with a *What to
download from Colab* section that zips outputs and triggers a
`google.colab.files.download` call.

## What NOT to do on Colab

- Don't download the **human** AlphaFold proteome on free-tier Colab —
  ~80 GB extracted does not fit reliably.
- Don't run `crypticip screen --mode full` on the whole yeast proteome
  without `--limit` first — the runtime may reclaim itself before it
  finishes.
- Don't commit any proteome data or screening CSVs back to the repo;
  they're large and not reproducible from source.
