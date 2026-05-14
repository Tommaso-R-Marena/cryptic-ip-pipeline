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
import os
from pathlib import Path

import pytest

# nbformat is a declared dev/test dependency. Allow opting out of notebook
# validation in environments where it's intentionally unavailable, but by
# default require it so CI fails loudly if it's missing.
if os.environ.get("CRYPTICIP_SKIP_NOTEBOOK_TESTS") == "1":
    nbformat = pytest.importorskip("nbformat")
else:
    import nbformat  # noqa: F401  (hard dependency; declared in pyproject [dev])

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


# ---------------------------------------------------------------------------
# Regression tests for previously reviewed bugs (Codex review on PR #2).
# ---------------------------------------------------------------------------


def _nb_code_cells(name):
    nb = nbformat.read(str(NB_DIR / name), as_version=4)
    return [c for c in nb.cells if c.cell_type == "code"]


def test_notebook_02_drive_symlink_cell_is_robust():
    """The Drive symlink cell must not crash when `results/` already exists.

    Regression test for Codex review: previously the cell did
    `if target.is_symlink() or target.exists(): pass` then unconditionally
    called `target.symlink_to(...)`, which raises FileExistsError.
    """
    cells = _nb_code_cells("02_yeast_screening_colab.ipynb")
    matching = [c for c in cells if "DRIVE_RESULTS" in "".join(c.source)]
    assert matching, "no Drive symlink cell found"
    src = "".join(matching[0].source)

    # Must handle the pre-existing-symlink case (unlink / replace).
    assert "unlink" in src, (
        "Drive symlink cell must unlink stale symlinks before replacing"
    )
    # Must back up an existing non-symlink results directory rather than
    # destructively deleting it.
    assert "local_backup_" in src, (
        "Drive symlink cell must preserve existing results/ via timestamped backup"
    )
    # Must ensure the Drive target exists before symlinking.
    assert "mkdir(parents=True, exist_ok=True)" in src
    # Should NOT call symlink_to unconditionally after a noop branch.
    assert "        pass\n    target.symlink_to" not in src, (
        "found old buggy pattern: pass-then-symlink unconditional"
    )


def test_notebook_04_uses_organism_variable():
    """The experimental-plan call must honor the configured ORGANISM variable.

    Regression test for Codex review: command was hardcoded to `yeast`.
    """
    cells = _nb_code_cells("04_experimental_prioritization_colab.ipynb")
    cmds = [
        "".join(c.source) for c in cells
        if "experimental-plan" in "".join(c.source)
    ]
    assert cmds, "no experimental-plan cell found"
    # At least one cell must use {ORGANISM} substitution.
    parametrized = [c for c in cmds if "{ORGANISM}" in c]
    assert parametrized, (
        "experimental-plan call must use {ORGANISM}, not hardcoded organism. "
        f"Got: {cmds!r}"
    )
    # The hardcoded form `--organism yeast` should not appear in any
    # experimental-plan invocation in this notebook.
    for c in cmds:
        assert "--organism yeast" not in c, (
            f"experimental-plan still hardcodes --organism yeast: {c!r}"
        )


def test_nbformat_declared_in_dev_dependencies():
    """nbformat must be listed in pyproject.toml [project.optional-dependencies] dev.

    Regression test for Codex review: missing dep made notebook tests
    silently skip in CI.
    """
    pyproject = (REPO_ROOT / "pyproject.toml").read_text()
    # Find the dev = [...] line and assert nbformat appears in the value.
    import re
    m = re.search(r"^dev\s*=\s*\[([^\]]*)\]", pyproject, re.MULTILINE)
    assert m, "no `dev = [...]` entry found in pyproject.toml"
    assert "nbformat" in m.group(1), (
        f"nbformat missing from dev deps. Got: {m.group(1)!r}"
    )
