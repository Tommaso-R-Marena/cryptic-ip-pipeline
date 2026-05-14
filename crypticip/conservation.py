"""Conservation analysis placeholder.

Full sequence-conservation analysis requires either downloading
pre-computed scores (ConSurf-DB, AlphaFold conservation tracks) or
running multiple-sequence-alignment + Rate4Site, neither of which we
want to make a hard dependency. This module exposes a thin façade so the
screening output gains a ``conservation_status`` field; if the user
provides a CSV at ``data/cache/conservation/<organism>.csv`` we will
read it.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass
class ConservationStore:
    by_residue: dict[tuple[str, int], float]
    source: str

    def score(self, chain: str, resseq: int) -> float | None:
        return self.by_residue.get((chain, resseq))


def load_organism_conservation(organism: str, *, cache_dir: Path) -> ConservationStore | None:
    p = cache_dir / "conservation" / f"{organism}.csv"
    if not p.exists():
        return None
    out: dict[tuple[str, int], float] = {}
    with p.open() as fh:
        for row in csv.DictReader(fh):
            try:
                out[(row["chain"], int(row["resseq"]))] = float(row["score"])
            except (KeyError, ValueError):
                continue
    return ConservationStore(by_residue=out, source=str(p))


def conservation_for_pocket(store: ConservationStore | None,
                            residues: Iterable[tuple[str, int, str]]) -> float | None:
    """Mean conservation score across pocket residues, if data is loaded."""
    if store is None:
        return None
    vals = [v for v in (store.score(c, r) for c, r, _ in residues) if v is not None]
    if not vals:
        return None
    return float(sum(vals) / len(vals))
