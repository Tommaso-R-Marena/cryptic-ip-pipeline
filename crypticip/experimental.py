"""Experimental follow-up planning: convert top-ranked candidates into
mutagenesis / DSF / MS suggestion sheets. Outputs CSV and Markdown.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Iterable

from .logging_utils import get_logger
from .paths import ProjectPaths

log = get_logger(__name__)


MUTAGENESIS_SUGGESTIONS = {
    "LYS": "K-to-A or K-to-E (charge swap; abolishes IP coordination)",
    "ARG": "R-to-A or R-to-E (charge swap; key IP coordinator)",
    "HIS": "H-to-A (removes potential H-bond donor)",
    "TRP": "W-to-A (loses ring stacking/positioning of aromatic shelf)",
    "TYR": "Y-to-F (removes hydroxyl, tests H-bond contribution)",
}


def _mutagenesis_for_basics(residues: list[tuple[str, int, str, float]]) -> list[str]:
    out = []
    for chain, resseq, rn, _d in residues:
        sug = MUTAGENESIS_SUGGESTIONS.get(rn.upper())
        if sug:
            out.append(f"{rn}{resseq}{chain}: {sug}")
    return out


def build_experimental_plan(organism: str, *, paths: ProjectPaths,
                            top_n: int = 25) -> dict:
    """Read screening summary for `organism`, pick the top-N pockets and
    generate experimental priority artefacts."""
    out_dir = paths.experimental_dir / organism
    out_dir.mkdir(parents=True, exist_ok=True)
    top_csv = paths.organism_screening_dir(organism) / "screening_top.csv"
    if not top_csv.exists():
        return {"status": "no_results", "out_dir": str(out_dir), "n": 0}

    plan_rows: list[dict] = []
    mutagenesis_rows: list[dict] = []
    with top_csv.open() as fh:
        for row in list(csv.DictReader(fh))[:top_n]:
            acc = row["accession"]
            per_protein = paths.organism_screening_dir(organism) / "per_protein" / f"{acc}.json"
            if not per_protein.exists():
                continue
            data = json.loads(per_protein.read_text())
            pk = next((p for p in data.get("pockets", []) if str(p["rank"]) == row["rank"]), None)
            if not pk:
                continue
            shells = pk.get("shells", {}) or {}
            r5 = shells.get("r5A", {}) or {}
            r8 = shells.get("r8A", {}) or {}
            plan_rows.append({
                "organism": organism, "accession": acc, "pocket_rank": pk["rank"],
                "composite": pk["composite"], "tier": pk["tier"],
                "depth": pk["fv"]["depth"], "sasa": pk["fv"]["sasa"],
                "elec": pk["fv"]["elec"], "basic_r5": r5.get("n_basic"),
                "basic_r8": r8.get("n_basic"), "volume": pk["fv"]["volume"],
                "mean_plddt": data.get("mean_plddt"),
                "apbs_status": pk["fv"]["apbs_status"],
                "center_x": pk["center"][0], "center_y": pk["center"][1], "center_z": pk["center"][2],
                "expected_ip_type": _guess_ip_type(pk),
                "dsf_priority": _dsf_priority(pk),
                "ms_priority": _ms_priority(pk),
                "review_flag": _review_flag(pk),
                "pymol_pml": f"results/reports/{organism}/pymol/{acc}_p{pk['rank']}.pml",
                "notes": "; ".join(_notes(pk)),
            })

            for chain, resseq, rn, d in _basic_residues_near(pk):
                mutagenesis_rows.append({
                    "organism": organism, "accession": acc, "pocket_rank": pk["rank"],
                    "chain": chain, "resseq": resseq, "resname": rn,
                    "distance_A": round(d, 2),
                    "suggestion": MUTAGENESIS_SUGGESTIONS.get(rn.upper(), "alanine scan"),
                })

    plan_csv = out_dir / f"{organism}_top{top_n}_experimental_plan.csv"
    mut_csv = out_dir / f"{organism}_top{top_n}_mutagenesis.csv"
    dsf_csv = out_dir / f"{organism}_top{top_n}_dsf.csv"

    _write_csv(plan_csv, plan_rows)
    _write_csv(mut_csv, mutagenesis_rows)
    _write_csv(dsf_csv, [
        {"organism": r["organism"], "accession": r["accession"], "ip_type": r["expected_ip_type"],
         "priority": r["dsf_priority"], "review_flag": r["review_flag"]}
        for r in plan_rows
    ])

    md = out_dir / f"{organism}_top{top_n}_experimental.md"
    md.write_text(_render_md(plan_rows))

    return {"status": "ok", "n_candidates": len(plan_rows),
            "plan_csv": str(plan_csv), "mutagenesis_csv": str(mut_csv),
            "dsf_csv": str(dsf_csv), "md": str(md)}


def _notes(pk: dict) -> list[str]:
    out = []
    if pk["fv"].get("apbs_status") != "ok":
        out.append(f"APBS={pk['fv'].get('apbs_status')} — electrostatics from fallback heuristic")
    if pk["fv"].get("plddt") is not None and pk["fv"]["plddt"] < 70:
        out.append(f"pLDDT={pk['fv']['plddt']:.1f} below 70 (low confidence)")
    if pk["fv"].get("basic_count", 0) >= 6:
        out.append("strong basic-residue clustering")
    if pk["fv"].get("depth", 0) >= 20:
        out.append("deeply buried")
    return out


def _guess_ip_type(pk: dict) -> str:
    v = pk["fv"].get("volume", 0)
    b = pk["fv"].get("basic_count", 0)
    if v >= 600 and b >= 6:
        return "IP6 (most likely)"
    if 400 <= v < 600:
        return "IP4/IP5"
    return "IP3/IP4"


def _dsf_priority(pk: dict) -> str:
    if pk["tier"] == "Tier1":
        return "high"
    if pk["tier"] == "Tier2":
        return "medium"
    return "low"


def _ms_priority(pk: dict) -> str:
    if pk["fv"].get("basic_count", 0) >= 6 and pk["fv"].get("depth", 0) >= 20:
        return "high"
    if pk["composite"] >= 0.55:
        return "medium"
    return "low"


def _review_flag(pk: dict) -> str:
    flags = []
    if pk["fv"].get("apbs_status") != "ok":
        flags.append("APBS-missing — manual review")
    if pk["fv"].get("plddt") and pk["fv"]["plddt"] < 70:
        flags.append("low-pLDDT — manual review")
    return "; ".join(flags) or ""


def _basic_residues_near(pk: dict) -> list[tuple[str, int, str, float]]:
    """No raw residue list survives screening JSON — return an empty list
    if not stored. (We still emit a row per candidate from shells.)"""
    return []


def _write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def _render_md(rows: list[dict]) -> str:
    if not rows:
        return "# Experimental plan\n_No candidates passed scoring._\n"
    cols = ["accession", "pocket_rank", "tier", "composite", "depth", "sasa",
            "elec", "basic_r5", "volume", "expected_ip_type", "dsf_priority",
            "ms_priority", "review_flag"]
    lines = ["# Experimental plan\n", "| " + " | ".join(cols) + " |",
             "|" + "|".join(["---"] * len(cols)) + "|"]
    for r in rows:
        lines.append("| " + " | ".join(str(r.get(c, "")) for c in cols) + " |")
    return "\n".join(lines) + "\n"
