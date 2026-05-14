"""PyMOL .pml session generation. The pml files are valid whether or not
PyMOL is installed; the optional render step uses PyMOL via subprocess.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Iterable

from .external_tools import check_pymol
from .logging_utils import get_logger

log = get_logger(__name__)


PML_TEMPLATE = """\
# PyMOL session for {accession} (cryptic-ip-pipeline)
load {pdb_path}, prot
hide everything
show cartoon
color grey80, prot
{plddt_block}
# pocket residues
select pocket_res, resi {pocket_resis}
show sticks, pocket_res
color cyan, pocket_res

# basic / acidic colouring
select basics, resn LYS+ARG+HIS and pocket_res
select acidics, resn ASP+GLU and pocket_res
color marine, basics
color salmon, acidics

# centre pseudoatom
pseudoatom pocket_center, pos=[{cx}, {cy}, {cz}]
show spheres, pocket_center
color yellow, pocket_center

# labels
label pocket_res and name CA, "%s%s" % (resn, resi)
{apbs_block}
orient pocket_res
zoom pocket_res, 5
bg_color white
set ray_shadows, 0
"""

PLDDT_PYMOL_SCRIPT = "spectrum b, blue_white_red, prot, minimum=50, maximum=90\n"


def write_pml(out_path: Path, *, pdb_path: Path, accession: str,
              pocket_residues: list[int],
              center: tuple[float, float, float],
              is_alphafold: bool = False,
              apbs_dx: Path | None = None) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pocket_resis = "+".join(str(r) for r in pocket_residues) or "1"
    plddt = PLDDT_PYMOL_SCRIPT if is_alphafold else ""
    apbs_block = ""
    if apbs_dx and apbs_dx.exists():
        apbs_block = (
            f"\n# APBS electrostatic surface\n"
            f"load {apbs_dx}, pot\n"
            "isosurface pos_iso, pot, 1.0\n"
            "isosurface neg_iso, pot, -1.0\n"
            "color marine, pos_iso\n"
            "color salmon, neg_iso\n"
            "set transparency, 0.4\n"
        )
    out_path.write_text(PML_TEMPLATE.format(
        accession=accession,
        pdb_path=pdb_path,
        pocket_resis=pocket_resis,
        cx=center[0], cy=center[1], cz=center[2],
        plddt_block=plddt,
        apbs_block=apbs_block,
    ))
    return out_path


def render_pml(pml_path: Path, *, image_path: Path | None = None,
               timeout: float = 120.0) -> dict:
    st = check_pymol()
    if not st.available:
        return {"status": "missing", "error": st.error}
    cmd = [st.path or "pymol", "-cq", str(pml_path)]
    if image_path:
        cmd += ["-d", f"png {image_path}"]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return {"status": "timeout"}
    return {"status": "ok" if proc.returncode == 0 else "failed",
            "stderr": proc.stderr[-500:], "stdout": proc.stdout[-500:]}


def generate_pymol_bundle(organism: str, *, paths, top_n: int = 50) -> dict:
    """Walk per-protein screening output and emit one .pml per top pocket."""
    import json
    from .paths import ProjectPaths
    paths: ProjectPaths
    out_dir = paths.organism_report_dir(organism) / "pymol"
    out_dir.mkdir(parents=True, exist_ok=True)
    top_csv = paths.organism_screening_dir(organism) / "screening_top.csv"
    if not top_csv.exists():
        return {"status": "no_results", "out_dir": str(out_dir), "n": 0}

    import csv
    seen = []
    with top_csv.open() as fh:
        for row in csv.DictReader(fh):
            seen.append(row)
            if len(seen) >= top_n:
                break

    n = 0
    for row in seen:
        acc = row["accession"]
        per_protein = paths.organism_screening_dir(organism) / "per_protein" / f"{acc}.json"
        if not per_protein.exists():
            continue
        data = json.loads(per_protein.read_text())
        cleaned_pdb = paths.organism_screening_dir(organism) / "_cleaned" / f"{acc}_clean.pdb"
        if not cleaned_pdb.exists():
            continue
        for pk in data.get("pockets", []):
            if str(pk["rank"]) != row["rank"]:
                continue
            resis: list[int] = []
            atm_file = pk.get("atm_pdb")
            # We don't always have the atm file path — fall back to nothing
            # and let the user open the .pml and select interactively.
            pml_path = out_dir / f"{acc}_p{pk['rank']}.pml"
            write_pml(pml_path,
                      pdb_path=cleaned_pdb,
                      accession=acc,
                      pocket_residues=resis,
                      center=tuple(pk.get("center", [0, 0, 0])),
                      is_alphafold=True)
            n += 1
            break
    return {"status": "ok", "out_dir": str(out_dir), "n": n}
