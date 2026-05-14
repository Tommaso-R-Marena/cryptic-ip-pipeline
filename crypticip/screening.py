"""Proteome screening pipeline (yeast / human / dictyostelium / arbitrary)."""
from __future__ import annotations

import csv
import json
import os
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable

from .alphafold import accession_from_filename, iter_proteome_files
from .config import Config
from .electrostatics import run_apbs
from .external_tools import env_report, check_apbs, check_pdb2pqr, check_fpocket, check_freesasa
from .fpocket import run_fpocket
from .logging_utils import get_logger
from .paths import ProjectPaths
from .pdb_io import parse_pdb_atoms
from .qc import qc_one_file
from .residues import residue_neighborhood, multi_shell_summary
from .sasa import compute_sasa, mean_sidechain_sasa
from .scoring import FeatureVector, composite_score, filter_flags, tier
from .structures import preprocess_structure

log = get_logger(__name__)


@dataclass
class ProteinResult:
    accession: str
    file: str
    mean_plddt: float | None
    n_residues: int
    n_pockets: int
    fpocket_status: str
    sasa_status: str
    apbs_status: str
    pockets: list[dict]
    error: str | None = None
    runtime_s: float | None = None


def _process_protein(args: tuple) -> dict:
    """Worker: returns a JSON-serialisable dict; never raises."""
    (pdb_path_str, cleaned_dir_str, top_k, mode,
     fpocket_avail, sasa_avail, apbs_avail) = args
    pdb_path = Path(pdb_path_str)
    cleaned_dir = Path(cleaned_dir_str)
    started = time.monotonic()
    accession = accession_from_filename(pdb_path) or pdb_path.stem

    try:
        meta = preprocess_structure(pdb_path, out_dir=cleaned_dir,
                                    name=accession, is_alphafold=True)
    except Exception as e:
        return asdict(ProteinResult(
            accession=accession, file=str(pdb_path), mean_plddt=None,
            n_residues=0, n_pockets=0, fpocket_status="failed", sasa_status="missing",
            apbs_status="missing", pockets=[],
            error=f"preprocess: {e}", runtime_s=time.monotonic() - started))

    if not fpocket_avail:
        return asdict(ProteinResult(
            accession=accession, file=str(pdb_path), mean_plddt=meta.mean_plddt,
            n_residues=meta.n_residues, n_pockets=0,
            fpocket_status="missing", sasa_status="skipped", apbs_status="skipped",
            pockets=[], error=None, runtime_s=time.monotonic() - started))

    pockets, run_meta = run_fpocket(meta.cleaned_path)
    if run_meta.get("status") != "ok":
        return asdict(ProteinResult(
            accession=accession, file=str(pdb_path), mean_plddt=meta.mean_plddt,
            n_residues=meta.n_residues, n_pockets=0,
            fpocket_status=run_meta.get("status", "missing"),
            sasa_status="skipped", apbs_status="skipped", pockets=[],
            error=run_meta.get("error"), runtime_s=time.monotonic() - started))

    sasa = compute_sasa(meta.cleaned_path) if sasa_avail else None
    sasa_status = sasa.status if sasa else "skipped"

    atoms = parse_pdb_atoms(meta.cleaned_path)
    candidate_pockets: list[dict] = []
    apbs_used = False
    for pk in pockets[:top_k]:
        nh = residue_neighborhood(atoms, pk.center, radius_A=5.0)
        shells = multi_shell_summary(atoms, pk.center)
        pocket_residue_keys = [(c, r) for (c, r, _rn) in pk.residues]
        sc_sasa = mean_sidechain_sasa(sasa, pocket_residue_keys) if sasa else None

        if mode == "full" and apbs_avail:
            apbs = run_apbs(meta.cleaned_path, pk.center)
            apbs_used = True
        else:
            apbs = run_apbs(meta.cleaned_path, pk.center,
                            apbs_status=type("S", (), {"available": False, "path": None, "error": "fast mode"})(),
                            pdb2pqr_status=type("S", (), {"available": False, "path": None, "error": "fast mode"})())

        fv = FeatureVector(
            depth=pk.depth,
            sasa=sc_sasa if sc_sasa is not None else 999.0,
            elec=apbs.get("potential_kT_e", 0.0),
            basic_count=float(nh.n_basic),
            volume=pk.volume,
            plddt=meta.mean_plddt,
            apbs_status=apbs.get("status", "missing"),
        )
        breakdown = composite_score(fv)
        flags = filter_flags(fv)
        candidate_pockets.append({
            "rank": pk.rank,
            "score": pk.score,
            "druggability": pk.druggability,
            "volume": pk.volume,
            "depth": pk.depth,
            "center": list(pk.center),
            "n_residues": len(pk.residues),
            "fv": asdict(fv),
            "shells": shells,
            "flags": flags,
            "composite": breakdown.composite,
            "tier": tier(breakdown.composite, flags, apbs_status=apbs.get("status", "missing")),
            "apbs_backend": apbs.get("backend"),
        })

    return asdict(ProteinResult(
        accession=accession, file=str(pdb_path), mean_plddt=meta.mean_plddt,
        n_residues=meta.n_residues, n_pockets=len(pockets),
        fpocket_status="ok",
        sasa_status=sasa_status,
        apbs_status="ok" if apbs_used else "skipped",
        pockets=candidate_pockets, error=None,
        runtime_s=time.monotonic() - started))


def screen_proteome(organism: str, *, cfg: Config, paths: ProjectPaths,
                    mode: str = "fast", workers: int = 4,
                    limit: int | None = None, resume: bool = True,
                    force: bool = False) -> dict:
    """Run screening over an organism's proteome dir. Resumable, parallel."""
    proteome_dir = paths.organism_proteome_dir(organism)
    out_dir = paths.organism_screening_dir(organism)
    out_dir.mkdir(parents=True, exist_ok=True)
    cleaned_dir = out_dir / "_cleaned"
    per_protein_dir = out_dir / "per_protein"
    per_protein_dir.mkdir(parents=True, exist_ok=True)
    cleaned_dir.mkdir(parents=True, exist_ok=True)

    files = list(iter_proteome_files(proteome_dir))
    if limit is not None:
        files = files[:limit]
    if not files:
        log.warning("No AlphaFold files in %s", proteome_dir)

    fpocket_avail = check_fpocket().available
    sasa_avail = check_freesasa().available
    apbs_avail = check_apbs().available and check_pdb2pqr().available

    pending: list[Path] = []
    for p in files:
        acc = accession_from_filename(p) or p.stem
        target = per_protein_dir / f"{acc}.json"
        if target.exists() and resume and not force:
            continue
        pending.append(p)

    log.info("Screening %s: %d files (%d pending, %d already done)",
             organism, len(files), len(pending), len(files) - len(pending))

    top_k = int(cfg.dotted("screening.per_protein_top_k", 5))
    args_iter = [(str(p), str(cleaned_dir), top_k, mode,
                  fpocket_avail, sasa_avail, apbs_avail) for p in pending]

    done = 0
    failures = 0
    if workers > 1 and pending:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(_process_protein, a): a for a in args_iter}
            for fut in as_completed(futures):
                rec = fut.result()
                _write_protein(rec, per_protein_dir)
                done += 1
                if rec.get("error"):
                    failures += 1
    else:
        for a in args_iter:
            rec = _process_protein(a)
            _write_protein(rec, per_protein_dir)
            done += 1
            if rec.get("error"):
                failures += 1

    aggregate = _aggregate(per_protein_dir, out_dir)
    aggregate["env"] = {n: st.to_dict() for n, st in env_report(cfg).items()}
    aggregate["n_files"] = len(files)
    aggregate["n_processed_this_run"] = done
    aggregate["n_failures_this_run"] = failures
    summary_path = out_dir / "summary.json"
    summary_path.write_text(json.dumps(aggregate, indent=2, default=str))
    return aggregate


def _write_protein(rec: dict, per_protein_dir: Path) -> None:
    acc = rec.get("accession") or "unknown"
    (per_protein_dir / f"{acc}.json").write_text(json.dumps(rec, indent=2, default=str))


def _aggregate(per_protein_dir: Path, out_dir: Path) -> dict:
    rows = []
    for jf in sorted(per_protein_dir.glob("*.json")):
        try:
            data = json.loads(jf.read_text())
        except Exception:
            continue
        for pk in data.get("pockets") or []:
            rows.append({
                "accession": data.get("accession"),
                "mean_plddt": data.get("mean_plddt"),
                "rank": pk["rank"],
                "depth": pk["fv"]["depth"],
                "sasa": pk["fv"]["sasa"],
                "elec": pk["fv"]["elec"],
                "basic": pk["fv"]["basic_count"],
                "volume": pk["fv"]["volume"],
                "composite": pk["composite"],
                "tier": pk["tier"],
                "apbs_status": pk["fv"]["apbs_status"],
            })

    csv_path = out_dir / "screening_results.csv"
    if rows:
        with csv_path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
    else:
        csv_path.write_text("accession,mean_plddt,rank,depth,sasa,elec,basic,volume,composite,tier,apbs_status\n")

    rows.sort(key=lambda r: r["composite"], reverse=True)
    top_csv = out_dir / "screening_top.csv"
    with top_csv.open("w", newline="") as fh:
        if rows:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows[:500])

    return {
        "csv": str(csv_path),
        "top_csv": str(top_csv),
        "n_rows": len(rows),
        "n_tier1": sum(1 for r in rows if r["tier"] == "Tier1"),
        "n_tier2": sum(1 for r in rows if r["tier"] == "Tier2"),
        "n_tier3": sum(1 for r in rows if r["tier"] == "Tier3"),
    }
