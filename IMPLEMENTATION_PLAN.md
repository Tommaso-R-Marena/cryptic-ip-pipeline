# Implementation Plan — Cryptic IP Pipeline (Production Refactor)

## Status of pre-existing work (preserved)

The original repo already contains a working validation pipeline in
`pipeline/*.py`, real results in `results/`, and figures in `figures/`. These
are **kept unchanged** as scientific record and reference implementations. The
production refactor lives under a new package `crypticip/` and a new CLI
`crypticip`. The old scripts are wrapped or referenced; no real results are
overwritten.

## Goals

1. Refactor the existing one-off scripts into a packaged, importable,
   testable Python library `crypticip` with a single `crypticip` CLI.
2. Add proteome ingestion + screening for yeast / human / dictyostelium.
3. Add proper validation reports, scoring sensitivity, experimental planning.
4. Make every external tool (fpocket / FreeSASA / pdb2pqr / APBS / PyMOL)
   optional with graceful fallbacks and unit-tested mocks.
5. Make the whole thing runnable in CI smoke mode using fixtures.

## Checkpoint structure

Each chunk below is implemented + committed independently.

- **C1 — Plan & scaffold**
  - `IMPLEMENTATION_PLAN.md` (this file)
  - Empty `crypticip/` package skeleton, `pyproject.toml`, `environment.yml`,
    `Dockerfile`, `Makefile`, `.github/workflows/{tests.yml,validation-smoke.yml}`,
    `config/*.yaml`.
- **C2 — Core utilities**
  - `crypticip/{paths,logging_utils,config,external_tools,cli}.py`
  - Wire up `crypticip check-env`.
- **C3 — IO + structure cleaning**
  - `crypticip/{pdb_io,structures,alphafold}.py`
- **C4 — Analysis modules with mocks**
  - `crypticip/{fpocket,sasa,electrostatics,residues,scoring}.py`
- **C5 — Validation + reporting**
  - `crypticip/{validation,reporting,statistics,qc,pymol,experimental}.py`
- **C6 — Proteome download + index + screening**
  - `crypticip/{screening}.py` + scripts under `scripts/`
- **C7 — Tests + CI smoke**
  - `tests/test_*.py` for every module, all runnable offline.
- **C8 — Docs + final report**

## Mock / fallback policy

External binaries (fpocket, freesasa, pdb2pqr, apbs, pymol, gh) may be
absent. The wrappers in `crypticip/external_tools.py` return a status
object `{available, path, version, error}`. Modules read this and either
run the real tool or fall back to a documented heuristic (e.g.,
distance-based charge potential instead of APBS), tagging the result
with `apbs_status: failed` etc. Unit tests use fixtures and never
require a real binary.

## Test policy

Every module ships a `tests/test_<module>.py` exercising the
import-and-parse path with checked-in fixtures. The CI smoke job runs
`pytest -m "not requires_external"` on Linux and only `crypticip
check-env` and a tiny synthetic validation against fixtures. The full
validation against real PDB/AlphaFold structures is documented in the
README but not run in CI (those structures are already in `results/`
from the historical run and remain authoritative).

## What is intentionally *not* done

- We do not re-download or re-run the original 9-structure validation;
  the existing `results/validation_results.json` and figures remain the
  scientific record. The new CLI can regenerate them if the user runs
  `crypticip download-validation && crypticip validate` with the tools
  installed.
- We do not actually screen yeast/human/dictyostelium proteomes in this
  PR. The plumbing is in place and tested with synthetic mini-proteomes;
  the user runs the long jobs themselves.
- We do not invent new scientific results. Anywhere a real tool is
  unavailable, the output records a `*_status: missing|failed` field.

## Final deliverable

- New `crypticip` Python package + `crypticip` CLI on the
  `feat/production-pipeline` branch.
- Unit tests passing under `pytest -q` with no external binaries.
- `IMPLEMENTATION_REPORT.md` summarizing what was built, what was
  tested, and what the user must run.
