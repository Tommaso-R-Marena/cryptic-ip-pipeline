# Screening

Per-protein pipeline:

1. **Clean** (`crypticip.structures.preprocess_structure`) — strip HETATMs,
   keep chain A, dump `<accession>_clean.pdb` under
   `results/screening/<organism>/_cleaned/`.
2. **fpocket** — run, parse `*_info.txt`, attach residue list and
   centroid from `pocket{N}_atm.pdb`.
3. **FreeSASA** (full mode) — per-residue SASA.
4. **APBS** (full mode) — pdb2pqr → APBS → potential at each pocket centre.
5. **Residue features** — basic/acidic counts, net charge, hydropathy,
   H-bond donors/acceptors at 5 / 8 / 10 Å shells.
6. **Score** — composite from `crypticip.scoring`, tier by
   `crypticip.scoring.tier`.

## Modes

- `--mode fast`  — fpocket + SASA + residue features (no APBS).
  Suitable for the first pass on 6 000–23 000 proteins.
- `--mode full`  — adds APBS; ~10–60 s/protein, run on the Tier-2+
  shortlist or on a server.

## CLI

```bash
crypticip screen --organism yeast --mode fast --workers 8 --resume
crypticip screen --organism yeast --mode full --workers 8 --resume \
                 --limit 500
```

`--resume` (default) skips proteins whose per-protein JSON already
exists. `--force` re-runs everything. Failures are isolated:
`--workers > 1` uses a process pool so a crash in one worker doesn't
take down the others; the failed protein's record carries an `error`
field.

## Outputs

Under `results/screening/<organism>/`:

- `per_protein/<accession>.json` — full per-protein result, JSON.
- `screening_results.csv` — flat row-per-pocket (top-K per protein).
- `screening_top.csv` — top 500 pockets globally by composite score.
- `summary.json` — run metadata, tier counts, env snapshot.

The CSV columns are: accession, mean_plddt, rank, depth, sasa, elec,
basic, volume, composite, tier, apbs_status.
