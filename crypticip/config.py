"""Configuration loading and merging.

Layering: ``default.yaml`` -> organism/validation/scoring config -> user
overrides. Everything goes through :func:`load_config` so the CLI and the
library see the same merged dictionary.
"""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

import yaml

from .paths import REPO_ROOT


DEFAULT_CONFIG_PATH = REPO_ROOT / "config" / "default.yaml"


class Config(dict):
    """``dict`` subclass with attribute access and dotted-path get."""

    def __getattr__(self, name: str) -> Any:
        try:
            return self[name]
        except KeyError as e:
            raise AttributeError(name) from e

    def __setattr__(self, name: str, value: Any) -> None:
        self[name] = value

    def dotted(self, path: str, default: Any = None) -> Any:
        node: Any = self
        for part in path.split("."):
            if isinstance(node, dict) and part in node:
                node = node[part]
            else:
                return default
        return node

    def hash(self) -> str:
        blob = json.dumps(self, sort_keys=True, default=str).encode()
        return hashlib.sha256(blob).hexdigest()[:12]


def _deep_merge(base: dict, override: dict) -> dict:
    out = copy.deepcopy(base)
    for k, v in (override or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = v
    return out


def _load_yaml(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")
    with path.open() as f:
        data = yaml.safe_load(f) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Top-level YAML in {path} must be a mapping")
    return data


def load_config(*configs: str | Path | dict | None,
                organism: str | None = None) -> Config:
    """Load and merge the default config with any number of overrides.

    Each positional arg is one of:
      * a path to a YAML file
      * a dict already in memory
      * ``None`` (skipped)
    If ``organism`` is given and ``config/<organism>.yaml`` exists it is
    merged in after the defaults.
    """
    merged: dict = _load_yaml(DEFAULT_CONFIG_PATH)
    if organism:
        org_path = REPO_ROOT / "config" / f"{organism}.yaml"
        if org_path.exists():
            merged = _deep_merge(merged, _load_yaml(org_path))

    for cfg in configs:
        if cfg is None:
            continue
        if isinstance(cfg, dict):
            merged = _deep_merge(merged, cfg)
        else:
            merged = _deep_merge(merged, _load_yaml(Path(cfg)))

    return Config(merged)


def list_named_configs() -> list[Path]:
    return sorted((REPO_ROOT / "config").glob("*.yaml"))
