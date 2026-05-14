"""Lightweight statistics used by validation + scoring sensitivity.

Pure NumPy / stdlib so we don't carry a SciPy import at module load.
"""
from __future__ import annotations

import math
import random
from dataclasses import dataclass
from typing import Iterable, Sequence


def cohen_d(a: Sequence[float], b: Sequence[float]) -> float:
    if len(a) < 2 or len(b) < 2:
        return 0.0
    ma = sum(a) / len(a); mb = sum(b) / len(b)
    va = sum((x - ma) ** 2 for x in a) / (len(a) - 1)
    vb = sum((x - mb) ** 2 for x in b) / (len(b) - 1)
    pooled = math.sqrt(((len(a) - 1) * va + (len(b) - 1) * vb) / (len(a) + len(b) - 2))
    if pooled == 0:
        return 0.0
    return (ma - mb) / pooled


def mann_whitney_u(a: Sequence[float], b: Sequence[float]) -> tuple[float, float]:
    """Two-sided MWU. Returns (U, approx_p)."""
    combined = sorted(((v, 0) for v in a), key=lambda x: x[0])
    combined += [(v, 1) for v in b]
    combined.sort(key=lambda x: x[0])
    ranks = {}
    for i, (v, _) in enumerate(combined, 1):
        ranks.setdefault(v, []).append(i)
    rank_a = 0.0
    for v in a:
        rs = ranks[v]
        rank_a += sum(rs) / len(rs)
    n1, n2 = len(a), len(b)
    u = rank_a - n1 * (n1 + 1) / 2
    u_other = n1 * n2 - u
    u_stat = min(u, u_other)
    mean_u = n1 * n2 / 2
    sd_u = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
    if sd_u == 0:
        return u_stat, 1.0
    z = (u_stat - mean_u) / sd_u
    p = 2 * (1 - _phi(abs(z)))
    return u_stat, max(min(p, 1.0), 0.0)


def _phi(x: float) -> float:
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))


def roc_auc(scores: Sequence[float], labels: Sequence[int]) -> float:
    """Compute ROC AUC. labels are 1=positive, 0=negative."""
    pairs = sorted(zip(scores, labels), key=lambda x: -x[0])
    pos = sum(labels)
    neg = len(labels) - pos
    if pos == 0 or neg == 0:
        return float("nan")
    tp = fp = 0
    auc = 0.0
    prev_fp = 0
    for _, y in pairs:
        if y:
            tp += 1
        else:
            fp += 1
            auc += tp
    return auc / (pos * neg)


@dataclass
class WeightSensitivityRow:
    weights: dict[str, float]
    auc: float
    pos_mean: float
    neg_mean: float
    separation: float


def weight_sensitivity(
    feature_table: list[dict],   # each row: {"label":1/0, "depth":..., "inv_sasa":..., "elec":..., "basic":..., "volume_fit":..., "plddt_penalty":...}
    *, n_samples: int = 200, seed: int = 0,
) -> list[WeightSensitivityRow]:
    """Dirichlet-random weight sweep over the 5 normalised components."""
    rng = random.Random(seed)
    out: list[WeightSensitivityRow] = []
    keys = ("depth", "inv_sasa", "elec", "basic", "volume_fit")
    for _ in range(n_samples):
        gammas = [rng.gammavariate(1.0, 1.0) for _ in keys]
        s = sum(gammas)
        w = {k: g / s for k, g in zip(keys, gammas)}
        scores = []
        labels = []
        for row in feature_table:
            sc = (w["depth"] * row["depth"] + w["inv_sasa"] * row["inv_sasa"]
                  + w["elec"] * row["elec"] + w["basic"] * row["basic"]
                  + w["volume_fit"] * row["volume_fit"]) - row.get("plddt_penalty", 0.0)
            scores.append(sc)
            labels.append(row["label"])
        auc = roc_auc(scores, labels)
        pos = [s for s, y in zip(scores, labels) if y]
        neg = [s for s, y in zip(scores, labels) if not y]
        pos_mean = sum(pos) / len(pos) if pos else 0.0
        neg_mean = sum(neg) / len(neg) if neg else 0.0
        out.append(WeightSensitivityRow(weights=w, auc=auc,
                                        pos_mean=pos_mean, neg_mean=neg_mean,
                                        separation=pos_mean - neg_mean))
    return out


def bootstrap_auc(scores: Sequence[float], labels: Sequence[int], *,
                  n_boot: int = 1000, seed: int = 0) -> tuple[float, float, float]:
    """Return mean AUC and (lo,hi) 95% percentile CI."""
    rng = random.Random(seed)
    n = len(scores)
    aucs: list[float] = []
    for _ in range(n_boot):
        idx = [rng.randrange(n) for _ in range(n)]
        s = [scores[i] for i in idx]
        l = [labels[i] for i in idx]
        if sum(l) == 0 or sum(l) == n:
            continue
        aucs.append(roc_auc(s, l))
    if not aucs:
        return float("nan"), float("nan"), float("nan")
    aucs.sort()
    lo = aucs[int(0.025 * len(aucs))]
    hi = aucs[int(0.975 * len(aucs))]
    return sum(aucs) / len(aucs), lo, hi


def permutation_p(scores: Sequence[float], labels: Sequence[int], *,
                  n_perm: int = 1000, seed: int = 0) -> float:
    """Permutation test for AUC > 0.5."""
    rng = random.Random(seed)
    baseline = roc_auc(scores, labels)
    if math.isnan(baseline):
        return float("nan")
    labels = list(labels)
    count = 0
    for _ in range(n_perm):
        rng.shuffle(labels)
        if roc_auc(scores, labels) >= baseline:
            count += 1
    return (count + 1) / (n_perm + 1)


def leave_one_out_auc(scores: Sequence[float], labels: Sequence[int]) -> list[float]:
    out: list[float] = []
    for i in range(len(scores)):
        s = list(scores[:i]) + list(scores[i + 1:])
        l = list(labels[:i]) + list(labels[i + 1:])
        if sum(l) == 0 or sum(l) == len(l):
            out.append(float("nan"))
        else:
            out.append(roc_auc(s, l))
    return out
