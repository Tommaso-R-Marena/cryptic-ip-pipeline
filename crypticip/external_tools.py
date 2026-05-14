"""Wrappers around external binaries (fpocket / freesasa / pdb2pqr / apbs / pymol).

Each ``ToolStatus`` records whether the tool is available, where, and at what
version, so downstream modules can run a real call or fall back to a documented
heuristic without crashing.
"""
from __future__ import annotations

import dataclasses as _dc
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Iterable

from .logging_utils import get_logger

log = get_logger(__name__)


@_dc.dataclass
class ToolStatus:
    name: str
    available: bool
    path: str | None = None
    version: str | None = None
    error: str | None = None

    def to_dict(self) -> dict:
        return _dc.asdict(self)


def _capture(cmd: Iterable[str], timeout: float = 10.0) -> tuple[int, str, str]:
    try:
        proc = subprocess.run(list(cmd), capture_output=True, text=True, timeout=timeout)
        return proc.returncode, proc.stdout or "", proc.stderr or ""
    except FileNotFoundError as e:
        return 127, "", str(e)
    except subprocess.TimeoutExpired as e:
        return 124, "", f"timeout: {e}"
    except OSError as e:
        return 126, "", str(e)


def _parse_version(text: str) -> str | None:
    m = re.search(r"(\d+\.\d+(?:\.\d+)?)", text or "")
    return m.group(1) if m else None


def which(name: str) -> str | None:
    p = shutil.which(name)
    return p


def check_fpocket(binary: str = "fpocket") -> ToolStatus:
    p = which(binary)
    if not p:
        return ToolStatus("fpocket", False, error=f"{binary} not in PATH")
    rc, out, err = _capture([p, "-h"])
    text = (out + err)
    return ToolStatus("fpocket", True, path=p, version=_parse_version(text))


def check_freesasa(use_python: bool = True) -> ToolStatus:
    if use_python:
        try:
            import freesasa  # type: ignore
            return ToolStatus("freesasa", True, path="python:freesasa",
                              version=getattr(freesasa, "__version__", "unknown"))
        except ImportError as e:
            py_err = str(e)
        # fall through to binary
    else:
        py_err = "python freesasa disabled by config"
    p = which("freesasa")
    if p:
        rc, out, err = _capture([p, "--version"])
        return ToolStatus("freesasa", True, path=p, version=_parse_version(out + err))
    return ToolStatus("freesasa", False, error=py_err)


def check_pdb2pqr(binary: str = "pdb2pqr") -> ToolStatus:
    p = which(binary) or which("pdb2pqr30") or which("pdb2pqr.py")
    if not p:
        return ToolStatus("pdb2pqr", False, error="pdb2pqr not in PATH")
    rc, out, err = _capture([p, "--version"])
    return ToolStatus("pdb2pqr", True, path=p, version=_parse_version(out + err))


def check_apbs(binary: str = "apbs") -> ToolStatus:
    p = which(binary)
    if not p:
        return ToolStatus("apbs", False, error="apbs not in PATH")
    rc, out, err = _capture([p, "--version"])
    return ToolStatus("apbs", True, path=p, version=_parse_version(out + err))


def check_pymol(binary: str = "pymol") -> ToolStatus:
    p = which(binary)
    if not p:
        return ToolStatus("pymol", False, error="pymol not in PATH (rendering disabled)")
    rc, out, err = _capture([p, "-c", "-q", "-d", "print('ok')"], timeout=20.0)
    return ToolStatus("pymol", True, path=p, version=_parse_version(out + err))


def env_report(cfg) -> dict[str, ToolStatus]:
    """Return a {name: ToolStatus} mapping for all required externals."""
    tools = cfg.get("tools", {}) if isinstance(cfg, dict) else {}
    return {
        "fpocket":  check_fpocket((tools.get("fpocket") or {}).get("binary", "fpocket")),
        "freesasa": check_freesasa((tools.get("freesasa") or {}).get("use_python", True)),
        "pdb2pqr":  check_pdb2pqr((tools.get("pdb2pqr") or {}).get("binary", "pdb2pqr")),
        "apbs":     check_apbs((tools.get("apbs") or {}).get("binary", "apbs")),
        "pymol":    check_pymol((tools.get("pymol") or {}).get("binary", "pymol")),
    }


def format_env_report(report: dict[str, ToolStatus]) -> str:
    lines = ["External tool check:"]
    for name, st in report.items():
        mark = "OK" if st.available else "missing"
        extra = f"v{st.version}" if st.version else (st.error or "")
        path = f" [{st.path}]" if st.path else ""
        lines.append(f"  - {name:8s} {mark:8s} {extra}{path}")
    return "\n".join(lines)
