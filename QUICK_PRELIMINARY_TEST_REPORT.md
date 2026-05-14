# Quick Preliminary Test Report

Date: 2026-05-14
Branch: `feat/colab-notebooks-and-smoke-tests` (base: `main`, which already includes the production pipeline merged via PR #1)
Python: 3.12, install: `pip install -e .`

This report captures a fast, sandbox-feasible sweep of the `crypticip` CLI
intended only as a smoke test of plumbing. **It does NOT replace real
validation or screening** — those require external binaries (fpocket,
FreeSASA, PDB2PQR, APBS, PyMOL) and large AlphaFold proteome downloads that
are unavailable in this sandbox.

## 1. Environment check

Command:

```
crypticip check-env
```

Result (`/tmp/crypticip_logs/check_env.txt`):

```
External tool check:
  - fpocket  missing  fpocket not in PATH
  - freesasa missing  No module named 'freesasa'
  - pdb2pqr  missing  pdb2pqr not in PATH
  - apbs     missing  apbs not in PATH
  - pymol    missing  pymol not in PATH (rendering disabled)
```

All five external dependencies are absent in the sandbox. The CLI prints the
intended fallback messaging. **Status:** PASS (CLI behavior correct;
degraded analysis expected).

## 2. Test suite

Command:

```
python -m pytest -q
```

Result: **100 passed, 25 skipped in 0.52s.**

The 25 skips are tests marked `requires_external` (fpocket / FreeSASA / APBS /
PyMOL not installed) — this is the intended behavior of the test markers.
**Status:** PASS.

## 3. CLI help

Command (and per-subcommand help captured in `/tmp/crypticip_logs/help_cmds.txt`):

```
crypticip --help
crypticip <sub> --help    # for each of:
  check-env, list-configs, download-validation, validate,
  download-proteome, index-proteome, screen, report, pymol,
  experimental-plan
```

All 10 subcommands print help with no traceback. **Status:** PASS.

## 4. `list-configs`

```
crypticip list-configs
```

Lists six bundled YAML configs (`default`, `validation`, `scoring`, `yeast`,
`human`, `dictyostelium`). **Status:** PASS.

## 5. Synthetic mini-proteome smoke workflow

To avoid the ~3 GB AlphaFold yeast tar, a synthetic 3-file "yeast" proteome
was placed at `/tmp/crypticip_smoke/data/proteomes/yeast/` using copies of
`tests/fixtures/tiny.pdb` renamed as `AF-P00001-F1-model_v4.pdb`,
`AF-P00002-...`, `AF-P00003-...`. A smoke config (`/tmp/crypticip_smoke/smoke_config.yaml`)
redirected `paths.*` to a scratch directory.

### 5a. `index-proteome`

```
crypticip index-proteome --organism yeast --config /tmp/crypticip_smoke/smoke_config.yaml
```

Output:

```
index_csv=/tmp/crypticip_smoke/data/proteomes/yeast/index.csv
qc_spotcheck=/tmp/crypticip_smoke/data/proteomes/yeast/qc_spotcheck.json
n_files=3
```

The `index.csv` contains accession / size / atom count / pLDDT columns for
all three synthetic structures. `qc_spotcheck.json` reports
`status=ok, n_pass=3, n_fail=0`. **Status:** PASS.

### 5b. `screen --mode fast`

```
crypticip screen --organism yeast --mode fast --workers 1 --limit 3 \
    --config /tmp/crypticip_smoke/smoke_config.yaml
```

Output:

```
{
  "csv": "...screening_results.csv",
  "top_csv": "...screening_top.csv",
  "n_rows": 0,
  "n_tier1": 0, "n_tier2": 0, "n_tier3": 0,
  "n_files": 3, "n_processed_this_run": 3, "n_failures_this_run": 0
}
```

`n_rows = 0` is the **expected** degraded-mode behavior: fpocket is not in
PATH, so the screening pipeline finds zero pockets. The plumbing is correct
(3 files iterated, 0 failures, CSVs written). **Status:** PASS (plumbing);
**DEGRADED** (no biological signal without fpocket).

### 5c. `report --organism yeast`

```
crypticip report --organism yeast --config /tmp/crypticip_smoke/smoke_config.yaml
```

Writes `screening_report.md`. **Status:** PASS.

### 5d. `pymol --organism yeast`

```
crypticip pymol --organism yeast --top 5 --config /tmp/crypticip_smoke/smoke_config.yaml
```

Returns `{"status": "ok", "n": 0}` (no candidates to render because of 5b).
**Status:** PASS (plumbing); **DEGRADED** (no PyMOL binary in PATH; even if
candidates existed, no PNG renders would be produced — `.pml` files would be
written, but the smoke set produced no candidates).

### 5e. `experimental-plan --organism yeast`

```
crypticip experimental-plan --organism yeast --top 5 \
    --config /tmp/crypticip_smoke/smoke_config.yaml
```

Writes mutagenesis / DSF / experimental-plan CSVs + Markdown. With
`n_candidates=0`, files are created but empty. **Status:** PASS (plumbing).

## 6. Validation smoke

`crypticip validate` requires the curated validation set (downloaded via
`download-validation` — needs network access to RCSB + AlphaFold) AND
fpocket/FreeSASA/APBS/PDB2PQR to produce meaningful pocket / SASA /
electrostatic features. Both prerequisites are unavailable in this sandbox.

To avoid fabricating validation conclusions, validation was **NOT** executed
here. The unit-test suite already exercises the validation report writer
(`tests/test_*` — included in the 100-pass count above) with synthetic
records.

**Status:** SKIPPED — must be run on a workstation with binaries +
internet access. The Colab notebook `01_validation_colab.ipynb` provides
a reproducible recipe.

## 7. Summary

| Step | Status | Notes |
|------|--------|-------|
| `check-env` | PASS | All 5 external tools correctly reported missing |
| `pytest` | PASS | 100 passed, 25 marker-skipped (no external tools) |
| CLI `--help` (10 subcommands) | PASS | All print usage with no traceback |
| `list-configs` | PASS | 6 bundled YAMLs |
| `index-proteome` (synthetic, 3 files) | PASS | index.csv + qc_spotcheck.json written |
| `screen --mode fast` | PASS plumbing / DEGRADED science | 0 pockets because fpocket missing |
| `report` | PASS | Markdown rendered |
| `pymol` | PASS plumbing / DEGRADED science | 0 candidates → 0 sessions |
| `experimental-plan` | PASS | CSVs + MD written (empty due to upstream) |
| `validate` | SKIPPED | Requires network + external binaries |

## 8. What must be run on a real workstation / server

1. **`crypticip check-env`** must show all five tools available. Install
   options:
   - conda: `mamba install -c conda-forge fpocket freesasa pdb2pqr apbs pymol-open-source`
   - or use the supplied `Dockerfile` (top of repo).
2. **`crypticip download-validation`** — needs internet (RCSB + AlphaFold).
3. **`crypticip validate`** — full real validation with fpocket pockets and
   APBS electrostatics. Expected output gates documented in
   `docs/validation_plan.md`.
4. **`crypticip download-proteome --organism yeast`** — pulls the
   ~3 GB UP000002311 AlphaFold tar (and ~20–80 GB for human). Use
   `--limit` to extract only N files for a small test.
5. **`crypticip screen --organism yeast --mode fast`** — first run with
   `--limit 100` to time it, then scale up. `--resume` makes long runs
   restartable.
6. **`crypticip pymol --render`** — produces PNG renders (requires
   PyMOL with a headless display, e.g. `pymol -cq`).

## 9. Caveats and honest notes

- The synthetic 3-residue PDB has no real pocket geometry; the
  `n_rows=0` / `n_pockets=0` results are diagnostic of plumbing, not of
  the biological method.
- This sandbox cannot validate fpocket-, APBS-, or PyMOL-dependent
  numerical results. **Do not** treat any number in this report as a
  scientific finding.
- The validation gate cannot be assessed here. The published gate
  outcomes are documented in `README.md` / `IMPLEMENTATION_REPORT.md`
  and were produced on a real workstation, not by this smoke run.
