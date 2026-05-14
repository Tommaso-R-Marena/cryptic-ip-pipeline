# Troubleshooting

## `crypticip check-env` says everything is missing

You're running with pure pip and don't have the bioconda tools.
Either:

```bash
conda env create -f environment.yml
conda activate crypticip
```

or install `fpocket`, `freesasa`, `pdb2pqr`, `apbs` manually. The
package and tests run fine without these binaries — you'll just get
`fpocket_status="missing"`, `apbs_status="fallback"`, etc.

## fpocket prints "could not open input file"

The pipeline copies the cleaned PDB into a temp dir before invoking
fpocket so paths with spaces / Unicode don't trip it up. If you still
see this, run `fpocket -h` directly — older versions are picky about
the working directory.

## APBS exits with code 1, no DX file

Almost always pdb2pqr produced an empty `.pqr` (often because the input
PDB has incomplete residues). Re-clean the structure with
`preprocess_structure(..., chain="A")` and try again. If APBS itself is
the problem, the run falls back to Coulomb and the result is tagged
`apbs_status="failed"`.

## FreeSASA python bindings install fails on macOS

Use `pip install --no-binary :all: freesasa` or install via bioconda.

## Slow proteome download

The AlphaFold tar for human is ~30 GB. Use `--resume` (default) and
prefer wired Ethernet; the streamed download writes its progress every
5 s.

## Tests don't find `crypticip`

`pip install -e .` first, or run from the repo root so the implicit
`sys.path` insertion in `tests/conftest.py` picks the package up.

## "PyMOL not in PATH" warning

PyMOL is optional. The `.pml` files are still produced; open them by
hand in PyMOL or install `pymol-open-source` from bioconda and re-run
`crypticip pymol --organism <…> --top 50 --render`.
