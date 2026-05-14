"""Markdown / HTML report generation. Matplotlib is imported lazily so it
isn't a hard dependency for the package import path used in tests."""
from __future__ import annotations

import csv
import json
import time
from pathlib import Path
from typing import Iterable

from .logging_utils import get_logger
from .paths import ProjectPaths

log = get_logger(__name__)


def _md_table(rows: list[dict], cols: list[str]) -> str:
    if not rows:
        return "_(no rows)_\n"
    out = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for r in rows:
        out.append("| " + " | ".join(_fmt(r.get(c)) for c in cols) + " |")
    return "\n".join(out) + "\n"


def _fmt(v) -> str:
    if v is None:
        return ""
    if isinstance(v, float):
        return f"{v:.3f}" if abs(v) < 1e4 else f"{v:.3g}"
    return str(v)


def write_validation_report(report: dict, *, paths: ProjectPaths) -> dict:
    out_dir = paths.reports_dir / "validation"
    out_dir.mkdir(parents=True, exist_ok=True)
    md = out_dir / "validation_report.md"
    html = out_dir / "validation_report.html"

    sections = []
    sections.append("# Cryptic IP Pipeline — Validation Report\n")
    sections.append(f"_Generated {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}_\n")
    summary = report["summary"]; gate = report["gate"]
    sections.append("## Summary\n")
    sections.append(_md_table([summary], list(summary.keys())))
    sections.append("## Gate\n")
    sections.append(_md_table([gate], list(gate.keys())))
    sections.append("## Per-structure results\n")
    cols = ["name", "category", "ip_type", "is_alphafold", "fpocket_status",
            "sasa_status", "apbs_status", "score", "tier"]
    rows = []
    for r in report["results"]:
        rows.append({k: r.get(k) for k in cols})
    sections.append(_md_table(rows, cols))

    md.write_text("\n".join(sections))
    html.write_text("<html><body><pre>\n" + "\n".join(sections) + "\n</pre></body></html>")

    figs = _make_validation_figures(report, out_dir)
    return {"md": str(md), "html": str(html), "figures": figs}


def _make_validation_figures(report: dict, out_dir: Path) -> list[str]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        log.warning("matplotlib missing — skipping figure generation")
        return []

    fig_paths: list[str] = []
    rows = report["results"]
    pos = [r for r in rows if r.get("category") == "positive" and r.get("score") is not None]
    neg = [r for r in rows if r.get("category") == "negative" and r.get("score") is not None]
    if pos or neg:
        fig, ax = plt.subplots(figsize=(6, 4))
        if pos:
            ax.bar([r["name"] for r in pos], [r["score"] for r in pos], color="tab:green", label="positive")
        if neg:
            ax.bar([r["name"] for r in neg], [r["score"] for r in neg], color="tab:red", label="negative")
        ax.set_ylabel("composite score")
        ax.set_ylim(0, 1)
        ax.set_xticklabels([r["name"] for r in pos + neg], rotation=45, ha="right")
        ax.legend()
        ax.set_title("Validation composite scores")
        fig.tight_layout()
        p = out_dir / "scores.png"
        fig.savefig(p, dpi=120)
        plt.close(fig)
        fig_paths.append(str(p))

    return fig_paths


def write_screening_report(organism: str, summary: dict, *,
                           paths: ProjectPaths) -> dict:
    out_dir = paths.organism_report_dir(organism)
    out_dir.mkdir(parents=True, exist_ok=True)
    md = out_dir / "screening_report.md"
    top_csv = Path(summary.get("top_csv", "")) if summary.get("top_csv") else None
    rows = []
    if top_csv and top_csv.exists():
        with top_csv.open() as fh:
            rows = list(csv.DictReader(fh))[:50]

    cols = ["accession", "rank", "depth", "sasa", "elec", "basic", "volume",
            "composite", "tier", "mean_plddt", "apbs_status"]
    md.write_text(
        f"# Screening report — {organism}\n\n"
        f"- Rows: {summary.get('n_rows')}    Tier1: {summary.get('n_tier1')}    "
        f"Tier2: {summary.get('n_tier2')}    Tier3: {summary.get('n_tier3')}\n\n"
        "## Top 50 pockets\n\n"
        + _md_table(rows, cols)
    )
    return {"md": str(md)}


def write_all_organism_report(organisms: Iterable[str], *,
                              paths: ProjectPaths) -> dict:
    out = paths.reports_dir / "all_organisms.md"
    sections = ["# Combined screening report\n"]
    for org in organisms:
        sum_path = paths.organism_screening_dir(org) / "summary.json"
        if not sum_path.exists():
            sections.append(f"## {org}\n_No results found_\n")
            continue
        s = json.loads(sum_path.read_text())
        sections.append(f"## {org}\n- rows={s.get('n_rows')} Tier1={s.get('n_tier1')} "
                        f"Tier2={s.get('n_tier2')} Tier3={s.get('n_tier3')}\n")
    out.write_text("\n".join(sections))
    return {"md": str(out)}
