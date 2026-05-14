"""Proteome QC and run metadata."""
from __future__ import annotations

import json
import os
import random
import socket
import subprocess
import sys
import time
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Iterable

from . import __version__
from .alphafold import accession_from_filename
from .logging_utils import get_logger
from .pdb_io import parse_pdb_atoms
from .structures import per_residue_ca_plddt, plddt_summary

log = get_logger(__name__)


def _git_commit(root: Path) -> str | None:
    try:
        out = subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"],
                             capture_output=True, text=True, timeout=5)
        return out.stdout.strip() or None
    except Exception:
        return None


def run_metadata(cfg, *, args: list[str], extra: dict | None = None) -> dict:
    """Common run-metadata block used by every CLI command."""
    from .paths import REPO_ROOT
    meta = {
        "tool": "crypticip",
        "version": __version__,
        "argv": args,
        "host": socket.gethostname(),
        "python": sys.version.split()[0],
        "platform": sys.platform,
        "started_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "git_commit": _git_commit(REPO_ROOT),
        "config_hash": cfg.hash() if hasattr(cfg, "hash") else None,
    }
    if extra:
        meta.update(extra)
    return meta


@dataclass
class FileQC:
    accession: str | None
    path: str
    size_bytes: int
    n_atoms: int
    n_residues: int
    mean_plddt: float | None
    median_plddt: float | None
    fraction_plddt_high: float | None
    fraction_plddt_low: float | None
    status: str = "ok"        # ok | empty | truncated | unreadable

    def to_row(self) -> dict:
        return asdict(self)


def qc_one_file(p: Path) -> FileQC:
    acc = accession_from_filename(p)
    try:
        size = p.stat().st_size
    except OSError:
        return FileQC(acc, str(p), 0, 0, 0, None, None, None, None, status="unreadable")
    if size == 0:
        return FileQC(acc, str(p), 0, 0, 0, None, None, None, None, status="empty")
    try:
        atoms = parse_pdb_atoms(p)
    except Exception:
        return FileQC(acc, str(p), size, 0, 0, None, None, None, None, status="unreadable")
    ca_plddt = per_residue_ca_plddt([a for a in atoms if a.record == "ATOM"])
    stats = plddt_summary(ca_plddt)
    n_res = len(ca_plddt)
    status = "ok" if (n_res > 0 and any(a.record == "ATOM" for a in atoms)) else "truncated"
    return FileQC(acc, str(p), size, len(atoms), n_res,
                  stats["mean"], stats["median"],
                  stats["fraction_high"], stats["fraction_low"], status=status)


def spot_check_proteome(proteome_dir: Path, *, n: int = 25,
                        seed: int = 0) -> dict:
    files = sorted(proteome_dir.glob("AF-*.pdb"))
    if not files:
        return {"status": "empty", "n_files": 0}
    rng = random.Random(seed)
    sample = rng.sample(files, min(n, len(files)))
    rows = [qc_one_file(p).to_row() for p in sample]
    ok = sum(1 for r in rows if r["status"] == "ok")
    return {
        "status": "ok",
        "n_files_total": len(files),
        "n_sampled": len(rows),
        "n_pass": ok,
        "n_fail": len(rows) - ok,
        "sample": rows,
    }
