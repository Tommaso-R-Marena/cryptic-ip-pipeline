"""fpocket wrapper + parser. Works on real fpocket output and on fixtures."""
from __future__ import annotations

import math
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Iterable

from .external_tools import check_fpocket, ToolStatus
from .logging_utils import get_logger
from .pdb_io import parse_pdb_atoms, Atom

log = get_logger(__name__)


@dataclass
class Pocket:
    rank: int
    score: float = 0.0
    druggability: float = 0.0
    volume: float = 0.0
    depth: float = 0.0
    mean_alpha_sphere_radius: float = 0.0
    mean_local_hydrophobic_density: float = 0.0
    charge_score: float = 0.0
    polarity_score: float = 0.0
    n_alpha_spheres: int = 0
    residues: list[tuple[str, int, str]] = field(default_factory=list)   # (chain, resseq, resname)
    center: tuple[float, float, float] = (0.0, 0.0, 0.0)
    atm_pdb: str | None = None

    def to_dict(self) -> dict:
        d = asdict(self)
        d["residues"] = [list(t) for t in self.residues]
        return d


# Ordered list of (regex matching the *full* normalised key, attribute name).
# Order matters: most-specific patterns first so e.g. "Local hydrophobic
# density score" doesn't collide with the bare "Score" key.
_KEY_RULES: list[tuple[str, str]] = [
    (r"^drug\s+score$", "druggability"),
    (r"^druggability\s+score$", "druggability"),
    (r"^real\s+volume", "volume"),
    (r"^pocket\s+volume", "volume"),
    (r"^cent\.\s*of\s*mass\s*-\s*alpha\s*sphere\s*max\s*dist", "depth"),
    (r"^mean\s+local\s+hydrophobic\s+density", "mean_local_hydrophobic_density"),
    (r"^local\s+hydrophobic\s+density\s+score", "mean_local_hydrophobic_density"),
    (r"^mean\s+alpha[\-\s]sphere\s+radius", "mean_alpha_sphere_radius"),
    (r"^number\s+of\s+alpha\s+spheres", "n_alpha_spheres"),
    (r"^charge\s+score", "charge_score"),
    (r"^polarity\s+score", "polarity_score"),
    (r"^score$", "score"),
]


def parse_info_text(text: str) -> list[Pocket]:
    """Parse the contents of an ``*_info.txt`` file."""
    pockets: list[Pocket] = []
    current: Pocket | None = None
    for raw in text.splitlines():
        line = raw.strip()
        m = re.match(r"Pocket\s+(\d+)\s*:", line)
        if m:
            if current is not None:
                pockets.append(current)
            current = Pocket(rank=int(m.group(1)))
            continue
        if current is None or ":" not in line:
            continue
        key, val = line.split(":", 1)
        # strip parenthesised units (we keep them in the key for matching)
        norm_key = key.strip().lower()
        target = None
        for pattern, attr in _KEY_RULES:
            if re.match(pattern, norm_key):
                target = attr
                break
        if target is None:
            continue
        try:
            fval = float(val.split()[0])
        except (ValueError, IndexError):
            continue
        if target == "n_alpha_spheres":
            setattr(current, target, int(fval))
        else:
            setattr(current, target, fval)
    if current is not None:
        pockets.append(current)
    return pockets


def parse_info_file(path: Path) -> list[Pocket]:
    return parse_info_text(Path(path).read_text())


def _residue_set(atoms: Iterable[Atom]) -> list[tuple[str, int, str]]:
    seen: dict[tuple[str, int], str] = {}
    for a in atoms:
        if a.record == "ATOM":
            seen.setdefault((a.chain, a.resseq), a.resname)
    return [(c, n, rn) for (c, n), rn in seen.items()]


def _centroid(atoms: list[Atom]) -> tuple[float, float, float]:
    if not atoms:
        return (0.0, 0.0, 0.0)
    n = len(atoms)
    return (sum(a.x for a in atoms) / n, sum(a.y for a in atoms) / n,
            sum(a.z for a in atoms) / n)


def _alpha_sphere_centroid(out_pdb: Path) -> tuple[float, float, float] | None:
    """Compute alpha-sphere centroid for each pocket from the
    ``*_out.pdb`` file (alpha spheres are HETATM POL/STP). Returns global centroid;
    per-pocket centers come from _residue_set centroid in the wrapper."""
    if not out_pdb.exists():
        return None
    xs, ys, zs = [], [], []
    for a in parse_pdb_atoms(out_pdb):
        if a.record == "HETATM" and a.name in ("POL", "STP", "APOL"):
            xs.append(a.x); ys.append(a.y); zs.append(a.z)
    if not xs:
        return None
    return (sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs))


def parse_pocket_atm(atm_pdb: Path) -> tuple[list[tuple[str, int, str]], tuple[float, float, float]]:
    """Read a ``pocket{N}_atm.pdb`` file: return (residues, centroid)."""
    atoms = parse_pdb_atoms(atm_pdb)
    residues = _residue_set(atoms)
    return residues, _centroid([a for a in atoms if a.record == "ATOM"])


def enrich_pockets_with_residues(pockets: list[Pocket], pockets_dir: Path) -> list[Pocket]:
    """Add residues + center from per-pocket ``pocket{rank}_atm.pdb`` files."""
    for pk in pockets:
        candidate = pockets_dir / f"pocket{pk.rank}_atm.pdb"
        if candidate.exists():
            residues, centroid = parse_pocket_atm(candidate)
            pk.residues = residues
            pk.center = centroid
            pk.atm_pdb = str(candidate)
    return pockets


def parse_fpocket_output(out_dir: Path, basename: str) -> list[Pocket]:
    """Top-level parser: find ``<basename>_info.txt`` + enrich from pocket atm files."""
    info = out_dir / f"{basename}_info.txt"
    if not info.exists():
        # Some fpocket versions name it without the basename
        infos = list(out_dir.glob("*_info.txt"))
        if not infos:
            return []
        info = infos[0]
    pockets = parse_info_file(info)
    pockets_dir = out_dir / "pockets"
    if pockets_dir.exists():
        pockets = enrich_pockets_with_residues(pockets, pockets_dir)
    return pockets


def run_fpocket(pdb_path: Path | str, *, work_dir: Path | None = None,
                keep_output: bool = True, status: ToolStatus | None = None) -> tuple[list[Pocket], dict]:
    """Run fpocket on ``pdb_path``. Returns ``(pockets, run_meta)``.

    If fpocket is not available, ``pockets`` is ``[]`` and
    ``run_meta['status']`` is ``"missing"``.
    """
    pdb_path = Path(pdb_path)
    status = status or check_fpocket()
    if not status.available:
        return [], {"status": "missing", "error": status.error, "out_dir": None}

    work = work_dir or Path(tempfile.mkdtemp(prefix="crypticip_fpocket_"))
    work.mkdir(parents=True, exist_ok=True)
    local_pdb = work / pdb_path.name
    shutil.copy(pdb_path, local_pdb)

    try:
        proc = subprocess.run(
            [status.path or "fpocket", "-f", str(local_pdb)],
            capture_output=True, text=True, cwd=work, timeout=600,
        )
    except subprocess.TimeoutExpired as e:
        return [], {"status": "timeout", "error": str(e), "out_dir": None}

    basename = local_pdb.stem
    out_dir = work / f"{basename}_out"
    if not out_dir.exists():
        return [], {"status": "failed", "error": proc.stderr[-500:] or "no output",
                    "out_dir": None, "returncode": proc.returncode}

    pockets = parse_fpocket_output(out_dir, basename)
    return pockets, {
        "status": "ok",
        "out_dir": str(out_dir),
        "n_pockets": len(pockets),
        "returncode": proc.returncode,
        "fpocket_version": status.version,
        "keep_output": keep_output,
    }
