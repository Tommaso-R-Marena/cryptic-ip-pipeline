"""APBS / PDB2PQR wrapper. DX parser is robust to trailing 'attribute' /
'component' lines APBS appends after the data block.

If pdb2pqr / APBS aren't installed, a distance-based Coulomb-like heuristic
is used and the result is tagged ``apbs_status='fallback'``.
"""
from __future__ import annotations

import math
import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from .external_tools import check_apbs, check_pdb2pqr, ToolStatus
from .logging_utils import get_logger
from .pdb_io import parse_pdb_atoms, Atom
from .residues import sidechain_terminal_atom, FORMAL_CHARGE

log = get_logger(__name__)


@dataclass
class DXGrid:
    origin: tuple[float, float, float]
    delta: tuple[float, float, float]            # diag of delta matrix
    shape: tuple[int, int, int]
    data: list[float] = field(default_factory=list)

    def potential_at(self, x: float, y: float, z: float) -> float:
        ox, oy, oz = self.origin
        dx, dy, dz = self.delta
        nx, ny, nz = self.shape
        i = int(round((x - ox) / dx)) if dx else 0
        j = int(round((y - oy) / dy)) if dy else 0
        k = int(round((z - oz) / dz)) if dz else 0
        i = max(0, min(nx - 1, i))
        j = max(0, min(ny - 1, j))
        k = max(0, min(nz - 1, k))
        return float(self.data[(i * ny + j) * nz + k])


_NUM = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")


def parse_dx(path: Path | str) -> DXGrid:
    """Parse an OpenDX scalar grid produced by APBS."""
    path = Path(path)
    origin = (0.0, 0.0, 0.0)
    deltas: list[tuple[float, float, float]] = []
    shape = (0, 0, 0)
    data: list[float] = []
    in_data = False
    with path.open() as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("object 1"):
                m = re.search(r"counts\s+(\d+)\s+(\d+)\s+(\d+)", line)
                if m:
                    shape = (int(m.group(1)), int(m.group(2)), int(m.group(3)))
                continue
            if line.startswith("origin"):
                nums = _NUM.findall(line)
                if len(nums) >= 3:
                    origin = (float(nums[0]), float(nums[1]), float(nums[2]))
                continue
            if line.startswith("delta"):
                nums = _NUM.findall(line)
                if len(nums) >= 3:
                    deltas.append((float(nums[0]), float(nums[1]), float(nums[2])))
                continue
            if line.startswith("object 3 class array"):
                in_data = True
                continue
            if line.startswith("object") or line.startswith("attribute") or line.startswith("component"):
                in_data = False
                continue
            if in_data:
                for tok in line.split():
                    try:
                        data.append(float(tok))
                    except ValueError:
                        in_data = False
                        break

    diag = (deltas[0][0] if deltas else 1.0,
            deltas[1][1] if len(deltas) > 1 else 1.0,
            deltas[2][2] if len(deltas) > 2 else 1.0)
    if not shape[0] and data:
        n = round(len(data) ** (1 / 3))
        shape = (n, n, n)
    return DXGrid(origin=origin, delta=diag, shape=shape, data=data)


def potential_at_center(dx_path: Path | str, center: tuple[float, float, float]) -> float:
    grid = parse_dx(dx_path)
    return grid.potential_at(*center)


def run_apbs(pdb_path: Path | str, center: tuple[float, float, float], *,
             ph: float = 7.0, force_field: str = "AMBER",
             apbs_status: ToolStatus | None = None,
             pdb2pqr_status: ToolStatus | None = None,
             keep_files: bool = False) -> dict:
    """Run pdb2pqr + apbs and return potential at ``center``.

    Result has either status=ok with `potential_kT_e` or status=fallback
    with a Coulomb-like estimate from formal-charge residue distances.
    """
    pdb_path = Path(pdb_path)
    apbs_status = apbs_status or check_apbs()
    pdb2pqr_status = pdb2pqr_status or check_pdb2pqr()
    if not (apbs_status.available and pdb2pqr_status.available):
        pot = _coulomb_fallback(pdb_path, center)
        return {"status": "fallback", "potential_kT_e": pot, "backend": "coulomb",
                "error": f"apbs={apbs_status.available} pdb2pqr={pdb2pqr_status.available}"}

    workdir = Path(tempfile.mkdtemp(prefix="crypticip_apbs_"))
    try:
        local = workdir / pdb_path.name
        shutil.copy(pdb_path, local)
        pqr = workdir / (pdb_path.stem + ".pqr")
        rc = subprocess.run([pdb2pqr_status.path or "pdb2pqr", "--ff", force_field,
                             f"--with-ph={ph}", str(local), str(pqr)],
                            capture_output=True, text=True, timeout=600).returncode
        if rc != 0 or not pqr.exists():
            pot = _coulomb_fallback(pdb_path, center)
            return {"status": "failed", "potential_kT_e": pot, "backend": "coulomb",
                    "error": "pdb2pqr failed"}

        # Minimal APBS input file
        apbsin = workdir / "apbs.in"
        apbsin.write_text(_apbs_input_template(pqr))
        rc = subprocess.run([apbs_status.path or "apbs", str(apbsin)],
                            capture_output=True, text=True, cwd=workdir, timeout=1800).returncode
        dx = workdir / "pot.dx"
        if rc != 0 or not dx.exists():
            pot = _coulomb_fallback(pdb_path, center)
            return {"status": "failed", "potential_kT_e": pot, "backend": "coulomb",
                    "error": "apbs failed"}

        pot = potential_at_center(dx, center)
        return {"status": "ok", "potential_kT_e": pot, "backend": "apbs",
                "dx_path": str(dx) if keep_files else None}
    finally:
        if not keep_files:
            shutil.rmtree(workdir, ignore_errors=True)


def _apbs_input_template(pqr: Path) -> str:
    return f"""read
  mol pqr {pqr.name}
end
elec name pot
  mg-auto
  dime 65 65 65
  cglen 60.0 60.0 60.0
  fglen 50.0 50.0 50.0
  cgcent mol 1
  fgcent mol 1
  mol 1
  lpbe
  bcfl sdh
  pdie 2.0
  sdie 78.54
  srfm smol
  chgm spl2
  sdens 10.0
  srad 1.4
  swin 0.3
  temp 298.15
  calcenergy no
  calcforce no
  write pot dx pot
end
print elecEnergy pot end
quit
"""


def _coulomb_fallback(pdb_path: Path, center: tuple[float, float, float],
                      *, dielectric: float = 4.0,
                      radius_A: float = 10.0) -> float:
    """Rough Coulomb-like potential at ``center`` in kT/e units, using
    formal charges on the side chain terminal atoms within ``radius_A``.
    Intentionally crude — meant only to keep the pipeline going when APBS
    is unavailable. Result is tagged ``apbs_status='fallback'`` upstream.
    """
    atoms = parse_pdb_atoms(pdb_path)
    by_residue: dict[tuple[str, int, str], list[Atom]] = {}
    for a in atoms:
        if a.record != "ATOM":
            continue
        by_residue.setdefault((a.chain, a.resseq, a.resname), []).append(a)
    k_kT_per_e2_A = 7.0   # ~kT/e units for kc * e^2 in water at 298 K, dielectric-scaled
    cx, cy, cz = center
    pot = 0.0
    for (chain, resseq, resname), residue_atoms in by_residue.items():
        charge = FORMAL_CHARGE.get(resname.upper(), 0.0)
        if charge == 0:
            continue
        term = sidechain_terminal_atom(resname, residue_atoms)
        if term is None:
            continue
        dx = term.x - cx; dy = term.y - cy; dz = term.z - cz
        d = math.sqrt(dx * dx + dy * dy + dz * dz)
        if d < 1e-3 or d > 15.0:
            continue
        pot += k_kT_per_e2_A * charge / (dielectric * d)
    return float(pot)
