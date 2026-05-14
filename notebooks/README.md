# Notebooks

Exploratory analyses live here. They are intentionally thin — the heavy
lifting belongs in `crypticip/` so it is testable. Start a notebook by
calling `crypticip.config.load_config(...)` and `crypticip.validation.validate_all(...)`.

## Google Colab notebooks

Step-by-step notebooks that can be opened and run on Colab. They clone
this repo, install the `crypticip` package, and walk through the CLI.
See `docs/colab.md` for the index and caveats.

| File | Purpose | Open |
|------|---------|------|
| `00_colab_quickstart.ipynb` | Clone, install, `check-env`, run tests, synthetic smoke workflow. | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/cryptic-ip-pipeline/blob/main/notebooks/00_colab_quickstart.ipynb) |
| `01_validation_colab.ipynb` | Run validation on the curated control set; install external binaries via condacolab + mamba (optional). | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/cryptic-ip-pipeline/blob/main/notebooks/01_validation_colab.ipynb) |
| `02_yeast_screening_colab.ipynb` | Download + index + screen the yeast AlphaFold proteome; PyMOL + experimental plan. | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/cryptic-ip-pipeline/blob/main/notebooks/02_yeast_screening_colab.ipynb) |
| `03_results_analysis_colab.ipynb` | Load screening CSVs, rank, plot score distribution / tier counts / feature heatmap. | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/cryptic-ip-pipeline/blob/main/notebooks/03_results_analysis_colab.ipynb) |
| `04_experimental_prioritization_colab.ipynb` | Generate mutagenesis + DSF plans for top candidates. | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/cryptic-ip-pipeline/blob/main/notebooks/04_experimental_prioritization_colab.ipynb) |
