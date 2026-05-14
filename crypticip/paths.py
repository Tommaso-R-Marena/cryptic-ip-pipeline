"""Project path resolution.

All directories the pipeline writes to are derived from ``ProjectPaths``,
which prefers config values but falls back to a stable layout rooted at
the repository top.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


def _repo_root() -> Path:
    here = Path(__file__).resolve()
    for cand in (here.parent.parent, *here.parents):
        if (cand / "pyproject.toml").exists() or (cand / ".git").exists():
            return cand
    return Path.cwd()


REPO_ROOT = _repo_root()


@dataclass
class ProjectPaths:
    root: Path = REPO_ROOT
    data_dir: Path = REPO_ROOT / "data"
    proteomes_dir: Path = REPO_ROOT / "data" / "proteomes"
    cache_dir: Path = REPO_ROOT / "data" / "cache"
    results_dir: Path = REPO_ROOT / "results"
    reports_dir: Path = REPO_ROOT / "results" / "reports"
    screening_dir: Path = REPO_ROOT / "results" / "screening"
    experimental_dir: Path = REPO_ROOT / "results" / "experimental"
    figures_dir: Path = REPO_ROOT / "figures"

    @classmethod
    def from_config(cls, cfg) -> "ProjectPaths":
        paths_cfg = (cfg.get("paths") if hasattr(cfg, "get") else cfg.paths) or {}
        root = REPO_ROOT

        def _abs(value: str | Path) -> Path:
            p = Path(value)
            return p if p.is_absolute() else (root / p)

        defaults = cls()
        return cls(
            root=root,
            data_dir=_abs(paths_cfg.get("data_dir", defaults.data_dir)),
            proteomes_dir=_abs(paths_cfg.get("proteomes_dir", defaults.proteomes_dir)),
            cache_dir=_abs(paths_cfg.get("cache_dir", defaults.cache_dir)),
            results_dir=_abs(paths_cfg.get("results_dir", defaults.results_dir)),
            reports_dir=_abs(paths_cfg.get("reports_dir", defaults.reports_dir)),
            screening_dir=_abs(paths_cfg.get("screening_dir", defaults.screening_dir)),
            experimental_dir=_abs(paths_cfg.get("experimental_dir", defaults.experimental_dir)),
            figures_dir=defaults.figures_dir,
        )

    def ensure(self) -> "ProjectPaths":
        for p in (self.data_dir, self.proteomes_dir, self.cache_dir, self.results_dir,
                  self.reports_dir, self.screening_dir, self.experimental_dir):
            p.mkdir(parents=True, exist_ok=True)
        return self

    def organism_proteome_dir(self, organism: str) -> Path:
        return self.proteomes_dir / organism

    def organism_screening_dir(self, organism: str) -> Path:
        return self.screening_dir / organism

    def organism_report_dir(self, organism: str) -> Path:
        return self.reports_dir / organism
