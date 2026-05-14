# Proteome ingestion

## Sources

| organism       | UniProt    | AlphaFold tar                                              |
|----------------|------------|------------------------------------------------------------|
| yeast          | UP000002311 | UP000002311_559292_YEAST_v4.tar  (~6 049 entries)         |
| human          | UP000005640 | UP000005640_9606_HUMAN_v4.tar    (~23 391 entries)        |
| dictyostelium  | UP000002195 | UP000002195_44689_DICDI_v4.tar   (~12 622 entries)        |

URLs and expected counts live in `config/{yeast,human,dictyostelium}.yaml`.

## Download

```bash
crypticip download-proteome --organism yeast --resume
```

`--resume` uses an HTTP `Range` header so partial tar files are picked
up. The default destination is `data/proteomes/<organism>/`. Add
`--no-extract` to keep the tar without extracting, or `--verify-only`
to just check the local tar exists and is non-zero.

## Indexing

```bash
crypticip index-proteome --organism yeast
```

Produces:

- `data/proteomes/yeast/index.csv` — one row per AF file: accession,
  size, atom count, residue count, mean/median pLDDT, fraction
  pLDDT≥70, fraction pLDDT<50, status (`ok`/`empty`/`truncated`).
- `data/proteomes/yeast/qc_spotcheck.json` — random 25-file spot-check.

## QC failure modes

The pipeline detects (and tags) the following before screening:

- `empty`         — zero-byte file (download truncation).
- `truncated`     — parsed atoms present but no CA residues (malformed).
- `unreadable`    — file couldn't be opened or parsed.
- `ok`            — file has CA atoms and parses cleanly.

`crypticip screen` skips non-`ok` files automatically.
