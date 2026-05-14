"""``crypticip`` command-line interface."""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Iterable

from . import __version__
from .config import load_config, list_named_configs
from .external_tools import env_report, format_env_report
from .logging_utils import configure_logging, get_logger
from .paths import ProjectPaths, REPO_ROOT
from .qc import run_metadata, spot_check_proteome

log = get_logger(__name__)


# -- argparse plumbing -------------------------------------------------------

def _add_common(p: argparse.ArgumentParser) -> None:
    p.add_argument("--config", action="append", default=None,
                   help="Optional YAML config (may be passed multiple times)")
    p.add_argument("--log-level", default="INFO")
    p.add_argument("--log-json", action="store_true")


def _load_cfg(args) -> "Config":
    return load_config(*(args.config or []),
                       organism=getattr(args, "organism", None))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="crypticip",
                                     description="Cryptic inositol-phosphate binding site detection")
    parser.add_argument("--version", action="version", version=f"crypticip {__version__}")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("check-env", help="Check external tools + python deps"); _add_common(p)

    p = sub.add_parser("list-configs", help="List bundled configs"); _add_common(p)

    p = sub.add_parser("download-validation", help="Download crystal + AlphaFold structures for validation")
    _add_common(p)

    p = sub.add_parser("validate", help="Run phase-1 validation against controls")
    _add_common(p)

    p = sub.add_parser("download-proteome", help="Download an AlphaFold proteome tar")
    p.add_argument("--organism", required=True, choices=["yeast", "human", "dictyostelium"])
    p.add_argument("--outdir", default=None)
    p.add_argument("--resume", action="store_true", default=True)
    p.add_argument("--force", action="store_true")
    p.add_argument("--no-extract", action="store_true")
    p.add_argument("--verify-only", action="store_true")
    p.add_argument("--limit", type=int, default=None,
                   help="Extract at most N files (smoke testing)")
    _add_common(p)

    p = sub.add_parser("index-proteome", help="QC and build proteome index CSV")
    p.add_argument("--organism", required=True)
    _add_common(p)

    p = sub.add_parser("screen", help="Run cryptic-IP screening on a proteome")
    p.add_argument("--organism", required=True)
    p.add_argument("--mode", choices=["fast", "full"], default="fast")
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--limit", type=int, default=None)
    p.add_argument("--resume", action="store_true", default=True)
    p.add_argument("--force", action="store_true")
    _add_common(p)

    p = sub.add_parser("report", help="Generate Markdown/HTML reports")
    p.add_argument("--organism", default=None)
    p.add_argument("--all", action="store_true",
                   help="Combine all organisms")
    p.add_argument("--validation", action="store_true",
                   help="Render the validation report")
    _add_common(p)

    p = sub.add_parser("pymol", help="Generate PyMOL .pml sessions for top hits")
    p.add_argument("--organism", required=True)
    p.add_argument("--top", type=int, default=50)
    p.add_argument("--render", action="store_true")
    _add_common(p)

    p = sub.add_parser("experimental-plan", help="Build experimental follow-up sheets")
    p.add_argument("--organism", default=None)
    p.add_argument("--all", action="store_true")
    p.add_argument("--top", type=int, default=25)
    _add_common(p)

    return parser


# -- command implementations -------------------------------------------------

def cmd_check_env(args) -> int:
    cfg = _load_cfg(args)
    report = env_report(cfg)
    print(format_env_report(report))
    missing = [n for n, st in report.items() if not st.available]
    if missing:
        print("\nNOTE: missing tools degrade analysis to fallback heuristics:")
        print("  - fpocket missing  -> screening returns zero pockets")
        print("  - freesasa missing -> SASA filter cannot be applied (records sasa_status='missing')")
        print("  - apbs/pdb2pqr missing -> electrostatic potential uses Coulomb fallback (apbs_status='fallback')")
        print("  - pymol missing    -> session generation works but no PNG render")
    return 0


def cmd_list_configs(args) -> int:
    for p in list_named_configs():
        print(p.relative_to(REPO_ROOT))
    return 0


def cmd_download_validation(args) -> int:
    from .alphafold import download_alphafold_pdb
    cfg = _load_cfg(args)
    paths = ProjectPaths.from_config(cfg).ensure()
    spec = cfg.dotted("validation_set") or {}
    pdb_dir = paths.data_dir / "pdb"; pdb_dir.mkdir(parents=True, exist_ok=True)
    af_dir = paths.data_dir / "alphafold"; af_dir.mkdir(parents=True, exist_ok=True)

    n_ok = 0; n_fail = 0
    import requests
    for category in ("positive_controls", "negative_controls"):
        for name, item in (spec.get(category) or {}).items():
            src = item.get("source", "pdb")
            pdb_id = item.get("pdb_id") or item.get("pdb") or name
            if src == "alphafold":
                accession = pdb_id
                # AF accession looks like AF-P78563-F1 — model URL needs full base
                res = download_alphafold_pdb(accession + "-model_v4", af_dir, force=False)
                target = af_dir / f"{accession}.pdb"
                if res.ok and res.path.exists():
                    res.path.rename(target)
                    n_ok += 1
                    print(f"OK  AF  {accession} -> {target}")
                else:
                    n_fail += 1
                    print(f"FAIL AF {accession}: {res.error}")
            else:
                target = pdb_dir / f"{pdb_id}.pdb"
                if target.exists():
                    n_ok += 1; continue
                try:
                    r = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=60)
                    r.raise_for_status()
                    target.write_bytes(r.content)
                    n_ok += 1
                    print(f"OK  PDB {pdb_id} -> {target}")
                except Exception as e:
                    n_fail += 1
                    print(f"FAIL PDB {pdb_id}: {e}")
    print(f"\nDownloaded: ok={n_ok} fail={n_fail}")
    return 0 if n_fail == 0 else 1


def cmd_validate(args) -> int:
    from .validation import validate_all, write_validation_outputs
    cfg = _load_cfg(args)
    paths = ProjectPaths.from_config(cfg).ensure()
    meta = run_metadata(cfg, args=sys.argv)
    print(f"crypticip validate (config_hash={meta['config_hash']})")
    report = validate_all(cfg, paths=paths)
    report["metadata"] = meta
    out = write_validation_outputs(report, paths=paths)
    print(json.dumps({"summary": report["summary"], "gate": report["gate"],
                      "outputs": out}, indent=2, default=str))
    return 0 if report["gate"]["overall_pass"] else 2


def cmd_download_proteome(args) -> int:
    cfg = _load_cfg(args)
    paths = ProjectPaths.from_config(cfg).ensure()
    p = cfg.dotted("proteome") or {}
    if not p:
        log.error("No proteome config for --organism=%s", args.organism)
        return 2
    out_dir = Path(args.outdir) if args.outdir else paths.organism_proteome_dir(args.organism)
    out_dir.mkdir(parents=True, exist_ok=True)
    tar = out_dir / Path(p["alphafold_tar_url"]).name

    if args.verify_only:
        if tar.exists():
            print(f"OK   {tar} ({tar.stat().st_size} bytes)")
            return 0
        print(f"MISSING {tar}")
        return 1

    from .alphafold import stream_download, extract_proteome_tar
    if not tar.exists() or args.force:
        res = stream_download(p["alphafold_tar_url"], tar, resume=args.resume)
        if not res.ok:
            log.error("download failed: %s", res.error)
            return 1
        print(f"downloaded -> {tar} ({res.bytes} bytes)")
    else:
        print(f"already present: {tar}")

    if args.no_extract:
        return 0
    report = extract_proteome_tar(tar, out_dir, limit=args.limit)
    print(json.dumps(report, indent=2))
    return 0


def cmd_index_proteome(args) -> int:
    import csv
    cfg = _load_cfg(args)
    paths = ProjectPaths.from_config(cfg).ensure()
    proteome_dir = paths.organism_proteome_dir(args.organism)
    files = sorted(proteome_dir.glob("AF-*.pdb"))
    if not files:
        log.error("No AlphaFold .pdb files under %s", proteome_dir)
        return 1
    index_csv = proteome_dir / "index.csv"
    qc_path = proteome_dir / "qc_spotcheck.json"
    from .qc import qc_one_file
    rows = []
    for fp in files:
        row = qc_one_file(fp).to_row()
        rows.append(row)
    with index_csv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    qc_path.write_text(json.dumps(spot_check_proteome(proteome_dir), indent=2))
    print(f"index_csv={index_csv}\nqc_spotcheck={qc_path}\nn_files={len(rows)}")
    return 0


def cmd_screen(args) -> int:
    from .screening import screen_proteome
    cfg = _load_cfg(args)
    paths = ProjectPaths.from_config(cfg).ensure()
    summary = screen_proteome(args.organism, cfg=cfg, paths=paths,
                              mode=args.mode, workers=args.workers,
                              limit=args.limit, resume=args.resume, force=args.force)
    print(json.dumps({k: v for k, v in summary.items() if k != "env"}, indent=2, default=str))
    return 0


def cmd_report(args) -> int:
    from .reporting import (write_validation_report, write_screening_report,
                            write_all_organism_report)
    cfg = _load_cfg(args)
    paths = ProjectPaths.from_config(cfg).ensure()
    if args.validation:
        report_path = paths.reports_dir / "validation" / "validation_results.json"
        if not report_path.exists():
            log.error("Run `crypticip validate` first — %s not found", report_path)
            return 1
        report = json.loads(report_path.read_text())
        out = write_validation_report(report, paths=paths)
        print(json.dumps(out, indent=2))
        return 0
    if args.all:
        organisms = ["yeast", "human", "dictyostelium"]
        out = write_all_organism_report(organisms, paths=paths)
        print(json.dumps(out, indent=2))
        return 0
    if not args.organism:
        log.error("Pass --organism, --all, or --validation")
        return 2
    summary_path = paths.organism_screening_dir(args.organism) / "summary.json"
    if not summary_path.exists():
        log.error("No screening output at %s", summary_path)
        return 1
    summary = json.loads(summary_path.read_text())
    out = write_screening_report(args.organism, summary, paths=paths)
    print(json.dumps(out, indent=2))
    return 0


def cmd_pymol(args) -> int:
    from .pymol import generate_pymol_bundle
    cfg = _load_cfg(args)
    paths = ProjectPaths.from_config(cfg).ensure()
    out = generate_pymol_bundle(args.organism, paths=paths, top_n=args.top)
    print(json.dumps(out, indent=2))
    return 0


def cmd_experimental_plan(args) -> int:
    from .experimental import build_experimental_plan
    cfg = _load_cfg(args)
    paths = ProjectPaths.from_config(cfg).ensure()
    targets = (["yeast", "human", "dictyostelium"] if args.all
               else ([args.organism] if args.organism else []))
    if not targets:
        log.error("Pass --organism or --all")
        return 2
    out: dict = {}
    for org in targets:
        out[org] = build_experimental_plan(org, paths=paths, top_n=args.top)
    print(json.dumps(out, indent=2))
    return 0


COMMANDS = {
    "check-env": cmd_check_env,
    "list-configs": cmd_list_configs,
    "download-validation": cmd_download_validation,
    "validate": cmd_validate,
    "download-proteome": cmd_download_proteome,
    "index-proteome": cmd_index_proteome,
    "screen": cmd_screen,
    "report": cmd_report,
    "pymol": cmd_pymol,
    "experimental-plan": cmd_experimental_plan,
}


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    configure_logging(args.log_level, args.log_json)
    handler = COMMANDS[args.cmd]
    try:
        return handler(args)
    except KeyboardInterrupt:
        print("interrupted", file=sys.stderr)
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
