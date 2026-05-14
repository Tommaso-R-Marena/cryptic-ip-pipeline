# Running on Google Colab

The `crypticip` CLI runs fine in a stock Google Colab Python environment.
**External scientific binaries** (`fpocket`, `freesasa`, `pdb2pqr`, `apbs`,
`pymol`) do not. The notebooks under `notebooks/` walk through the full
workflow with honest fallback behaviour where binaries are unavailable.

## Notebooks

| Order | Notebook | What it does | Time | Disk | Open |
|-------|----------|--------------|------|------|------|
| 1 | [`00_colab_quickstart.ipynb`](../notebooks/00_colab_quickstart.ipynb) | Clone, install, `check-env`, run pytest, synthetic smoke workflow | ~5 min | < 100 MB | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/cryptic-ip-pipeline/blob/main/notebooks/00_colab_quickstart.ipynb) |
| 2 | [`01_validation_colab.ipynb`](../notebooks/01_validation_colab.ipynb) | Optional condacolab install of external binaries; `download-validation`; `validate`; ADAR2 diagnostics | ~10-15 min real / ~3 min fallback | ~50 MB | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/cryptic-ip-pipeline/blob/main/notebooks/01_validation_colab.ipynb) |
| 3 | [`02_yeast_screening_colab.ipynb`](../notebooks/02_yeast_screening_colab.ipynb) | Yeast AlphaFold proteome download, index, fast/full screen, PyMOL bundles, experimental plan | 30 min – several h | up to ~6 GB | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/cryptic-ip-pipeline/blob/main/notebooks/02_yeast_screening_colab.ipynb) |
| 4 | [`03_results_analysis_colab.ipynb`](../notebooks/03_results_analysis_colab.ipynb) | Load screening CSVs, rank, plot distributions / tiers / feature heatmap | minutes | < 10 MB | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/cryptic-ip-pipeline/blob/main/notebooks/03_results_analysis_colab.ipynb) |
| 5 | [`04_experimental_prioritization_colab.ipynb`](../notebooks/04_experimental_prioritization_colab.ipynb) | Generate mutagenesis + DSF plans for top candidates | minutes | < 10 MB | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/cryptic-ip-pipeline/blob/main/notebooks/04_experimental_prioritization_colab.ipynb) |

Badge links target the `main` branch so they remain valid after PR merge.
Open each by clicking the **Open in Colab** badge above, or directly with:

```
https://colab.research.google.com/github/Tommaso-R-Marena/cryptic-ip-pipeline/blob/main/notebooks/00_colab_quickstart.ipynb
```

(replace the filename for each notebook).

## External binaries on Colab

The only reliable Colab path for `fpocket` + `apbs` + `pdb2pqr` + `pymol-open-source`
+ `freesasa` is **condacolab + mamba**:

```python
!pip install -q condacolab
import condacolab; condacolab.install()       # restarts the kernel
# after restart:
!mamba install -y -c conda-forge -c bioconda fpocket freesasa pdb2pqr apbs pymol-open-source python-freesasa
```

Caveats:
- `condacolab.install()` **restarts the kernel** — you must re-`cd` into
  the repo and re-`pip install -e .` after the restart.
- The install takes 5-10 minutes and ~3 GB of disk.
- PyMOL's `apt-get` package on Colab pulls in a heavy GUI stack; prefer
  `pymol-open-source` from conda-forge.
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

Colab runtimes are ephemeral. To keep outputs:

- Mount Drive: `from google.colab import drive; drive.mount('/content/drive')`,
  then symlink `results/` → a Drive folder before running CLI commands.
- Or zip + download at the end of a run (see `02_yeast_screening_colab.ipynb`,
  section 11).

## What NOT to do on Colab

- Don't download the **human** AlphaFold proteome on free-tier Colab —
  ~80 GB extracted does not fit reliably.
- Don't run `crypticip screen --mode full` on the whole yeast proteome
  without `--limit` first — the runtime may reclaim itself before it
  finishes.
- Don't commit any proteome data or screening CSVs back to the repo;
  they're large and not reproducible from source.
