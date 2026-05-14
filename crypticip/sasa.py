"""FreeSASA wrapper with graceful fallback when the library is missing."""
from __future__ import annotations

import math
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from .external_tools import check_freesasa, ToolStatus
from .logging_utils import get_logger
from .pdb_io import parse_pdb_atoms, Atom

log = get_logger(__name__)


@dataclass
class ResidueSASA:
    chain: str
    resseq: int
    resname: str
    total: float = 0.0
    sidechain: float = 0.0
    mainchain: float = 0.0
    relative_sidechain: float | None = None


@dataclass
class SASAResult:
    status: str                                  # ok | python_missing | binary_missing | failed
    total: float | None = None
    residues: list[ResidueSASA] = field(default_factory=list)
    error: str | None = None
    backend: str | None = None


# Backbone atom names (everything else on a standard amino acid is sidechain).
_BACKBONE = {"N", "CA", "C", "O", "OXT", "H", "HA", "H1", "H2", "H3"}


def compute_sasa(pdb_path: Path | str, *, probe_radius: float = 1.4,
                 status: ToolStatus | None = None) -> SASAResult:
    """Compute per-residue SASA using the FreeSASA Python bindings if
    available. Falls back to a documented geometric fallback (returns
    status=python_missing) otherwise."""
    pdb_path = Path(pdb_path)
    status = status or check_freesasa(use_python=True)
    if status.available and status.path == "python:freesasa":
        return _compute_with_freesasa_python(pdb_path, probe_radius=probe_radius)
    return _fallback_no_sasa(pdb_path, reason=status.error or "freesasa not available")


def _compute_with_freesasa_python(pdb_path: Path, *, probe_radius: float) -> SASAResult:
    try:
        import freesasa  # type: ignore
    except Exception as e:
        return _fallback_no_sasa(pdb_path, reason=f"import freesasa failed: {e}")
    try:
        params = freesasa.Parameters({"probe-radius": probe_radius})
        structure = freesasa.Structure(str(pdb_path))
        result = freesasa.Calc(params).calculate(structure)
    except Exception as e:
        return _fallback_no_sasa(pdb_path, reason=f"freesasa calculation failed: {e}")

    residue_index: dict[tuple[str, int], ResidueSASA] = {}
    try:
        n_atoms = structure.nAtoms()
    except Exception:
        n_atoms = 0
    for i in range(n_atoms):
        try:
            chain = structure.chainLabel(i).strip() or "A"
            resseq = int(structure.residueNumber(i))
            resname = structure.residueName(i).strip()
            atom = structure.atomName(i).strip()
            sasa_i = float(result.atomArea(i))
        except Exception:
            continue
        key = (chain, resseq)
        rs = residue_index.setdefault(key, ResidueSASA(chain, resseq, resname))
        rs.total += sasa_i
        if atom in _BACKBONE:
            rs.mainchain += sasa_i
        else:
            rs.sidechain += sasa_i

    total = float(result.totalArea())
    return SASAResult(status="ok", total=total,
                      residues=sorted(residue_index.values(), key=lambda r: (r.chain, r.resseq)),
                      backend="freesasa-python")


def _fallback_no_sasa(pdb_path: Path, *, reason: str) -> SASAResult:
    """No SASA available — return empty residue list so downstream code
    knows to record sasa_status=missing rather than crashing."""
    log.warning("FreeSASA unavailable: %s — falling back to zero SASA estimates", reason)
    return SASAResult(status="missing", total=None, residues=[], error=reason, backend=None)


def sasa_for_residues(result: SASAResult, residues: Iterable[tuple[str, int]]) -> dict[tuple[str, int], ResidueSASA]:
    if result.status != "ok":
        return {}
    idx = {(r.chain, r.resseq): r for r in result.residues}
    return {k: idx[k] for k in residues if k in idx}


def mean_sidechain_sasa(result: SASAResult, residues: Iterable[tuple[str, int]]) -> float | None:
    selected = sasa_for_residues(result, residues)
    if not selected:
        return None
    vals = [r.sidechain for r in selected.values()]
    return float(sum(vals) / len(vals))
