"""Lightweight checks for the Colab notebooks.

These tests do NOT execute notebook cells (the heavy ones download
proteomes). They only verify that each notebook:
- is valid JSON,
- is a valid nbformat notebook,
- starts with the expected H1 heading,
- contains the expected install + crypticip CLI calls.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

nbformat = pytest.importorskip("nbformat")

REPO_ROOT = Path(__file__).resolve().parent.parent
NB_DIR = REPO_ROOT / "notebooks"

EXPECTED = {
    "00_colab_quickstart.ipynb": {
        "h1": "00 - Cryptic IP Pipeline: Colab Quickstart",
        "must_contain": [
            "pip install -q -e .",
            "crypticip --version",
            "crypticip check-env",
            "pytest",
            "index-proteome",
            "screen",
        ],
    },
    "01_validation_colab.ipynb": {
        "h1": "01 - Validation on the Curated Control Set",
        "must_contain": [
            "download-validation",
            "crypticip validate",
            "report --validation",
            "condacolab",
        ],
    },
    "02_yeast_screening_colab.ipynb": {
        "h1": "02 - Yeast Proteome Screening on Colab",
        "must_contain": [
            "download-proteome --organism yeast",
            "index-proteome --organism yeast",
            "screen --organism yeast",
            "pymol --organism yeast",
            "experimental-plan --organism yeast",
            "--limit",
        ],
    },
    "03_results_analysis_colab.ipynb": {
        "h1": "03 - Results Analysis",
        "must_contain": [
            "screening_results.csv",
            "score",
            "tier",
        ],
    },
    "04_experimental_prioritization_colab.ipynb": {
        "h1": "04 - Experimental Prioritization",
        "must_contain": [
            "experimental-plan",
            "mutagenesis",
            "dsf",
        ],
    },
}


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_notebook_is_valid_json(name):
    path = NB_DIR / name
    assert path.exists(), f"missing notebook: {path}"
    json.loads(path.read_text())


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_notebook_is_valid_nbformat(name):
    nb = nbformat.read(str(NB_DIR / name), as_version=4)
    nbformat.validate(nb)
    assert nb.cells, f"{name}: no cells"


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_notebook_starts_with_h1(name):
    spec = EXPECTED[name]
    nb = nbformat.read(str(NB_DIR / name), as_version=4)
    first_md = next((c for c in nb.cells if c.cell_type == "markdown"), None)
    assert first_md is not None, f"{name}: no markdown cell"
    src = first_md.source.lstrip()
    assert src.startswith(f"# {spec['h1']}"), (
        f"{name}: expected H1 '# {spec['h1']}', got: {src[:80]!r}"
    )


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_notebook_contains_expected_strings(name):
    spec = EXPECTED[name]
    nb = nbformat.read(str(NB_DIR / name), as_version=4)
    blob = "\n".join(c.source for c in nb.cells)
    missing = [s for s in spec["must_contain"] if s.lower() not in blob.lower()]
    assert not missing, f"{name}: missing required snippets: {missing}"


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_notebook_outputs_are_empty(name):
    """We commit notebooks without execution output to keep diffs clean."""
    nb = nbformat.read(str(NB_DIR / name), as_version=4)
    for i, c in enumerate(nb.cells):
        if c.cell_type == "code":
            assert not c.get("outputs"), (
                f"{name} cell {i}: unexpected output in committed notebook"
            )
            assert c.get("execution_count") in (None, 0), (
                f"{name} cell {i}: unexpected execution_count"
            )
