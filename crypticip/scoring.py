"""Composite scoring for candidate pockets.

The default scheme matches what the existing analysis used, expressed
through a clean weight + normalization config so it can be swept by the
sensitivity tools in :mod:`crypticip.statistics`.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable


@dataclass
class FeatureVector:
    depth: float = 0.0
    sasa: float = 0.0
    elec: float = 0.0
    basic_count: float = 0.0
    volume: float = 0.0
    plddt: float | None = None
    apbs_status: str = "ok"      # ok | fallback | failed | missing


@dataclass
class ScoreBreakdown:
    composite: float
    n_depth: float
    n_inv_sasa: float
    n_elec: float
    n_basic: float
    n_volume_fit: float
    plddt_penalty: float
    weights: dict[str, float]


def _clamp01(v: float) -> float:
    if math.isnan(v):
        return 0.0
    return max(0.0, min(1.0, v))


def _norm_depth(d: float, norm: float) -> float:
    return _clamp01(d / norm)


def _norm_inv_sasa(s: float, norm: float) -> float:
    return _clamp01(1.0 - s / norm)


def _norm_elec(e: float, norm: float) -> float:
    return _clamp01(max(e, 0.0) / norm)


def _norm_basic(b: float, norm: float) -> float:
    return _clamp01(b / norm)


def _norm_volume_fit(v: float, lo: float, hi: float) -> float:
    if v <= 0:
        return 0.0
    if lo <= v <= hi:
        return 1.0
    if v < lo:
        return _clamp01(v / lo)
    # over-volume falls off linearly until 2x the upper bound
    return _clamp01(max(0.0, 1.0 - (v - hi) / max(hi, 1.0)))


def composite_score(fv: FeatureVector, *, weights: dict[str, float] | None = None,
                    norms: dict[str, float] | None = None,
                    plddt_penalty: dict | None = None) -> ScoreBreakdown:
    """Compute the composite score and a per-component breakdown."""
    w = dict({"depth": 0.25, "inv_sasa": 0.25, "elec": 0.20,
              "basic": 0.20, "volume_fit": 0.10}, **(weights or {}))
    norms = dict({"depth_A": 30.0, "sasa_A2": 150.0, "elec_kT_e": 30.0,
                  "basic_count": 8.0, "volume_lo": 300.0,
                  "volume_hi": 800.0}, **(norms or {}))

    n_depth = _norm_depth(fv.depth, norms["depth_A"])
    n_inv_sasa = _norm_inv_sasa(fv.sasa, norms["sasa_A2"])
    n_elec = _norm_elec(fv.elec, norms["elec_kT_e"])
    n_basic = _norm_basic(fv.basic_count, norms["basic_count"])
    n_vfit = _norm_volume_fit(fv.volume, norms["volume_lo"], norms["volume_hi"])

    raw = (w["depth"] * n_depth + w["inv_sasa"] * n_inv_sasa
           + w["elec"] * n_elec + w["basic"] * n_basic
           + w["volume_fit"] * n_vfit)

    penalty = 0.0
    pp = plddt_penalty or {"enabled": True, "threshold": 70.0, "penalty": 0.2}
    if pp.get("enabled") and fv.plddt is not None and fv.plddt < pp["threshold"]:
        penalty = float(pp["penalty"])

    composite = max(0.0, min(1.0, raw - penalty))
    return ScoreBreakdown(
        composite=composite,
        n_depth=n_depth, n_inv_sasa=n_inv_sasa, n_elec=n_elec,
        n_basic=n_basic, n_volume_fit=n_vfit, plddt_penalty=penalty,
        weights=w,
    )


def filter_flags(fv: FeatureVector, *, criteria: dict | None = None) -> dict[str, bool]:
    """Per-criterion pass/fail flags. Always returns ``ok`` booleans."""
    c = dict({
        "depth_min_A": 15.0,
        "sasa_max_A2": 5.0,
        "potential_min_kT_e": 5.0,
        "basic_residues_min": 4,
        "volume_min_A3": 300.0,
        "volume_max_A3": 800.0,
        "plddt_min": 70.0,
    }, **(criteria or {}))
    return {
        "depth_ok": fv.depth >= c["depth_min_A"],
        "sasa_ok": fv.sasa <= c["sasa_max_A2"],
        "potential_ok": fv.elec >= c["potential_min_kT_e"],
        "basic_ok": fv.basic_count >= c["basic_residues_min"],
        "volume_ok": c["volume_min_A3"] <= fv.volume <= c["volume_max_A3"],
        "plddt_ok": fv.plddt is None or fv.plddt >= c["plddt_min"],
    }


def tier(score: float, flags: dict[str, bool], apbs_status: str = "ok") -> str:
    """Promote pockets to Tier 1/2/3 or reject. APBS-missing pockets cannot
    enter Tier 1 unless the user later promotes them manually."""
    passes = sum(1 for v in flags.values() if v)
    total = len(flags)
    if passes == total and score >= 0.65 and apbs_status == "ok":
        return "Tier1"
    if passes >= total - 1 and score >= 0.55:
        return "Tier2"
    if score >= 0.40:
        return "Tier3"
    return "Reject"
