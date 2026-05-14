"""Phase-1 validation orchestration: run the full pipeline on the curated
set of positive and negative controls and decide PASS/FAIL against the gate.
"""
from __future__ import annotations

import csv
import json
import time
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Iterable

from .config import Config
from .electrostatics import run_apbs
from .external_tools import env_report
from .fpocket import run_fpocket, Pocket
from .logging_utils import get_logger
from .paths import ProjectPaths
from .pdb_io import parse_pdb_atoms
from .qc import run_metadata
from .residues import multi_shell_summary, residue_neighborhood
from .sasa import compute_sasa, sasa_for_residues, mean_sidechain_sasa
from .scoring import (FeatureVector, composite_score, filter_flags, tier)
from .statistics import cohen_d, mann_whitney_u, roc_auc, bootstrap_auc
from .structures import preprocess_structure, StructureMeta

log = get_logger(__name__)


@dataclass
class StructureValidation:
    name: str
    category: str
    pdb_id: str
    chain: str
    ip_type: str
    source: str
    is_alphafold: bool
    cleaned_path: str | None
    n_pockets: int = 0
    fpocket_status: str = "missing"
    best_pocket: dict | None = None
    fv: dict | None = None
    flags: dict[str, bool] | None = None
    score: float | None = None
    score_breakdown: dict | None = None
    sasa_status: str = "missing"
    apbs_status: str = "missing"
    mean_plddt: float | None = None
    tier: str | None = None
    ip_residue_sasa: dict | None = None
    notes: list[str] = field(default_factory=list)


def _select_best_pocket(pockets: list[Pocket]) -> Pocket | None:
    if not pockets:
        return None
    return min(pockets, key=lambda pk: pk.rank)


def _maybe_known_ip_pocket(pockets: list[Pocket], known_residues: list[int]) -> Pocket | None:
    """Pick the pocket whose residue set best overlaps the known IP residues."""
    if not pockets or not known_residues:
        return None
    best = None
    best_overlap = 0
    for pk in pockets:
        ov = sum(1 for (_chain, resseq, _rn) in pk.residues if resseq in known_residues)
        if ov > best_overlap:
            best_overlap = ov
            best = pk
    return best if best_overlap >= 2 else None


def validate_structure(meta: StructureMeta, spec: dict, *,
                       paths: ProjectPaths, cfg: Config) -> StructureValidation:
    crit = cfg.dotted("criteria") or {}

    fpocket_status_obj = None
    pockets, run_meta = run_fpocket(meta.cleaned_path)
    fpocket_status = run_meta.get("status", "missing")

    # Choose the pocket of interest: known-IP overlap (positives with known
    # residues) or rank #1.
    known_res = spec.get("coordinating_residues") or []
    chosen = _maybe_known_ip_pocket(pockets, known_res) or _select_best_pocket(pockets)

    out = StructureValidation(
        name=meta.name,
        category=spec.get("category", "unknown"),
        pdb_id=spec.get("pdb_id", meta.name),
        chain=spec.get("chain", meta.chain or "A"),
        ip_type=spec.get("ip_type", ""),
        source=spec.get("source", "pdb"),
        is_alphafold=meta.is_alphafold,
        cleaned_path=str(meta.cleaned_path) if meta.cleaned_path else None,
        n_pockets=len(pockets),
        fpocket_status=fpocket_status,
        mean_plddt=meta.mean_plddt,
    )

    if chosen is None:
        out.notes.append("no_pockets_detected")
        out.tier = "Reject"
        return out

    # SASA (per IP residue / pocket residue)
    sasa_result = compute_sasa(meta.cleaned_path)
    out.sasa_status = sasa_result.status
    coord_keys = [(out.chain, r) for r in known_res] or [(c, r) for (c, r, _rn) in chosen.residues]
    mean_sc = mean_sidechain_sasa(sasa_result, coord_keys)
    if known_res:
        # Per-residue dump for the IP coordinating residues.
        per = {}
        for k, rs in sasa_for_residues(sasa_result, coord_keys).items():
            per[f"{k[1]}"] = {"total": rs.total, "sidechain": rs.sidechain}
        out.ip_residue_sasa = per

    # Electrostatics at pocket center
    apbs = run_apbs(meta.cleaned_path, chosen.center)
    out.apbs_status = apbs.get("status", "missing")
    potential = apbs.get("potential_kT_e", 0.0)

    # Residue neighbourhood
    atoms = parse_pdb_atoms(meta.cleaned_path)
    nh = residue_neighborhood(atoms, chosen.center,
                              radius_A=float((crit.get("basic_residues_radius_A", 5.0))))
    shells = multi_shell_summary(atoms, chosen.center)

    fv = FeatureVector(
        depth=chosen.depth,
        sasa=mean_sc if mean_sc is not None else (
            999.0 if sasa_result.status != "ok" else 999.0),
        elec=potential,
        basic_count=float(nh.n_basic),
        volume=chosen.volume,
        plddt=meta.mean_plddt,
        apbs_status=out.apbs_status,
    )
    breakdown = composite_score(
        fv,
        weights=cfg.dotted("scoring.weights"),
        norms=cfg.dotted("scoring.norms"),
        plddt_penalty=cfg.dotted("scoring.plddt_penalty"),
    )
    flags = filter_flags(fv, criteria=crit)

    out.best_pocket = {
        "rank": chosen.rank,
        "score": chosen.score,
        "druggability": chosen.druggability,
        "volume": chosen.volume,
        "depth": chosen.depth,
        "center": list(chosen.center),
        "n_residues": len(chosen.residues),
    }
    out.fv = {**asdict(fv), "shells": shells}
    out.flags = flags
    out.score = breakdown.composite
    out.score_breakdown = asdict(breakdown)
    out.tier = tier(breakdown.composite, flags, apbs_status=out.apbs_status)

    if out.sasa_status != "ok":
        out.notes.append("sasa_missing — flags.sasa_ok is informational only")
    if out.apbs_status != "ok":
        out.notes.append(f"apbs_{out.apbs_status} — potential uses Coulomb fallback")

    return out


def validate_all(cfg: Config, *, paths: ProjectPaths,
                 data_root: Path | None = None) -> dict:
    spec = cfg.dotted("validation_set") or {}
    data_root = data_root or (paths.data_dir)
    cleaned_dir = paths.data_dir / "cleaned" / "validation"
    cleaned_dir.mkdir(parents=True, exist_ok=True)

    results: list[StructureValidation] = []
    for category in ("positive_controls", "negative_controls"):
        items = (spec.get(category) or {})
        for name, item in items.items():
            item = dict(item)
            item.setdefault("category", "positive" if category == "positive_controls" else "negative")
            pdb_id = item.get("pdb_id") or item.get("pdb") or name
            source = item.get("source", "pdb")
            chain = item.get("chain", "A")
            if source == "alphafold":
                candidate = data_root / "alphafold" / f"{pdb_id}.pdb"
            else:
                candidate = data_root / "pdb" / f"{pdb_id}.pdb"
            if not candidate.exists():
                log.warning("Missing %s for %s — skipping", candidate, name)
                results.append(StructureValidation(
                    name=name, category=item["category"], pdb_id=pdb_id,
                    chain=chain, ip_type=item.get("ip_type", ""), source=source,
                    is_alphafold=(source == "alphafold"),
                    cleaned_path=None, tier="Reject", notes=["input_missing"]))
                continue
            meta = preprocess_structure(candidate, out_dir=cleaned_dir, name=name,
                                        is_alphafold=(source == "alphafold"), chain=chain)
            results.append(validate_structure(meta, item, paths=paths, cfg=cfg))

    return _gate_report(results, cfg=cfg)


def _gate_report(results: list[StructureValidation], *, cfg: Config) -> dict:
    rows = [asdict(r) for r in results]
    pos = [r for r in results if r.category == "positive" and r.score is not None]
    neg = [r for r in results if r.category == "negative" and r.score is not None]

    pos_scores = [r.score for r in pos]
    neg_scores = [r.score for r in neg]
    pos_depths = [(r.fv or {}).get("depth", 0.0) for r in pos]
    neg_depths = [(r.fv or {}).get("depth", 0.0) for r in neg]

    pos_mean = sum(pos_scores) / len(pos_scores) if pos_scores else 0.0
    neg_mean = sum(neg_scores) / len(neg_scores) if neg_scores else 0.0
    sep = pos_mean - neg_mean

    auc = roc_auc(pos_scores + neg_scores,
                  [1] * len(pos_scores) + [0] * len(neg_scores)) if (pos_scores and neg_scores) else float("nan")

    cohend = cohen_d(pos_scores, neg_scores) if (len(pos_scores) >= 2 and len(neg_scores) >= 2) else 0.0
    if len(pos_scores) >= 2 and len(neg_scores) >= 2:
        _u, p_mwu = mann_whitney_u(pos_scores, neg_scores)
    else:
        p_mwu = float("nan")

    adar2 = next((r for r in results if r.name == "ADAR2_crystal"), None)
    adar2_pocket_rank = adar2.best_pocket["rank"] if (adar2 and adar2.best_pocket) else None

    gate_cfg = (cfg.dotted("validation_gate") or {})
    gate = {
        "adar2_rank_ok": (adar2_pocket_rank is not None
                         and adar2_pocket_rank <= int(gate_cfg.get("adar2_ip6_pocket_rank_max", 3))),
        "separation_ok": sep >= float(gate_cfg.get("positive_vs_negative_mean_separation_min", 0.05)),
        "depth_separation_ok": ((sum(pos_depths) / len(pos_depths) if pos_depths else 0)
                                - (sum(neg_depths) / len(neg_depths) if neg_depths else 0))
            >= float(gate_cfg.get("positive_mean_depth_over_negative_min_A", 5.0)),
        "n_positives_passing":
            sum(1 for r in pos if r.flags and r.flags.get("depth_ok") and r.flags.get("basic_ok")),
    }
    gate["positives_passing_ok"] = gate["n_positives_passing"] >= int(gate_cfg.get("required_positives_passing_filter_count", 2))
    gate["overall_pass"] = all(v for k, v in gate.items() if isinstance(v, bool))

    boot_mean = boot_lo = boot_hi = float("nan")
    if pos_scores and neg_scores:
        boot_mean, boot_lo, boot_hi = bootstrap_auc(pos_scores + neg_scores,
                                                    [1] * len(pos_scores) + [0] * len(neg_scores))

    return {
        "results": rows,
        "summary": {
            "n_positives": len(pos),
            "n_negatives": len(neg),
            "pos_mean_score": pos_mean,
            "neg_mean_score": neg_mean,
            "separation": sep,
            "auc": auc,
            "auc_boot_mean": boot_mean,
            "auc_boot_ci": [boot_lo, boot_hi],
            "cohens_d": cohend,
            "mannwhitney_p": p_mwu,
            "adar2_pocket_rank": adar2_pocket_rank,
        },
        "gate": gate,
        "env": {name: st.to_dict() for name, st in env_report(cfg).items()},
    }


def write_validation_outputs(report: dict, *, paths: ProjectPaths) -> dict:
    """Write validation_results.json, validation_summary.csv, validation_gate_report.md.

    Backs up any pre-existing files with timestamped suffixes.
    """
    outdir = paths.reports_dir / "validation"
    outdir.mkdir(parents=True, exist_ok=True)
    ts = time.strftime("%Y%m%dT%H%M%S")

    def _safe(name: str) -> Path:
        target = outdir / name
        if target.exists():
            backup = outdir / f"{target.stem}.{ts}{target.suffix}"
            target.replace(backup)
        return target

    json_path = _safe("validation_results.json")
    csv_path = _safe("validation_summary.csv")
    md_path = _safe("validation_gate_report.md")
    json_path.write_text(json.dumps(report, indent=2, default=str))

    rows = report["results"]
    fieldnames = ["name", "category", "pdb_id", "ip_type", "is_alphafold",
                  "fpocket_status", "sasa_status", "apbs_status",
                  "best_pocket_rank", "depth", "sasa", "elec", "basic", "volume",
                  "mean_plddt", "score", "tier"]
    with csv_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            best = r.get("best_pocket") or {}
            fv = r.get("fv") or {}
            w.writerow({
                "name": r["name"], "category": r["category"], "pdb_id": r["pdb_id"],
                "ip_type": r["ip_type"], "is_alphafold": r["is_alphafold"],
                "fpocket_status": r["fpocket_status"],
                "sasa_status": r["sasa_status"], "apbs_status": r["apbs_status"],
                "best_pocket_rank": best.get("rank"),
                "depth": fv.get("depth"), "sasa": fv.get("sasa"),
                "elec": fv.get("elec"), "basic": fv.get("basic_count"),
                "volume": fv.get("volume"), "mean_plddt": r.get("mean_plddt"),
                "score": r.get("score"), "tier": r.get("tier"),
            })

    md_path.write_text(_render_gate_md(report))
    return {"json": str(json_path), "csv": str(csv_path), "md": str(md_path)}


def _render_gate_md(report: dict) -> str:
    g = report["gate"]; s = report["summary"]; env = report["env"]
    lines = [
        "# Validation Gate Report",
        "",
        f"- Positives: {s['n_positives']}    Negatives: {s['n_negatives']}",
        f"- Mean score, pos vs neg: **{s['pos_mean_score']:.3f}** vs **{s['neg_mean_score']:.3f}**  (Δ = {s['separation']:+.3f})",
        f"- ROC AUC: **{s['auc']}**    bootstrap mean {s['auc_boot_mean']} [{s['auc_boot_ci'][0]}, {s['auc_boot_ci'][1]}]",
        f"- Cohen's d: {s['cohens_d']:.3f}   MWU p ≈ {s['mannwhitney_p']}",
        f"- ADAR2 IP6 pocket rank: **{s['adar2_pocket_rank']}**",
        "",
        "## Gate criteria",
        f"- adar2_rank_ok:        {g['adar2_rank_ok']}",
        f"- separation_ok:        {g['separation_ok']}",
        f"- depth_separation_ok:  {g['depth_separation_ok']}",
        f"- positives_passing_ok: {g['positives_passing_ok']}  ({g['n_positives_passing']} positives passing filter)",
        f"- **overall_pass:       {g['overall_pass']}**",
        "",
        "## External tool status",
    ]
    for n, st in env.items():
        lines.append(f"- {n}: available={st['available']} version={st.get('version')} path={st.get('path')}")
    return "\n".join(lines) + "\n"
