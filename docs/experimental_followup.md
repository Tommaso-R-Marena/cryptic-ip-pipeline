# Experimental follow-up

`crypticip experimental-plan --organism yeast --top 25` reads the
screening top list and emits, under
`results/experimental/<organism>/`:

- `<organism>_top25_experimental_plan.csv` — one row per candidate
  with composite score, tier, basic-residue counts (5/8 Å),
  depth, SASA, APBS status, expected IP species, DSF / MS priorities,
  pLDDT, review flag (e.g. "APBS-missing — manual review"), and the
  PyMOL `.pml` path.
- `<organism>_top25_mutagenesis.csv` — proposed mutations for basic
  pocket residues (K→A/E, R→A/E, H→A, plus W→A and Y→F where applicable).
- `<organism>_top25_dsf.csv` — DSF priority list (high / medium / low).
- `<organism>_top25_experimental.md` — same data rendered for humans.

## Heuristics

- **Expected IP species** — IP6 if volume ≥ 600 Å³ and ≥ 6 basic
  residues within 5 Å; IP4/IP5 for 400–600 Å³; otherwise IP3/IP4.
- **DSF priority** — Tier 1 → high, Tier 2 → medium, otherwise low.
- **MS priority** — high if ≥ 6 basic residues AND depth ≥ 20 Å, else
  medium if composite ≥ 0.55, else low.
- **Review flags** — APBS-missing pockets, low pLDDT (< 70).

These are starting points for a wet-lab triage, not statements of fact.
Always cross-check the PyMOL session and the literature before
investing reagent time.
