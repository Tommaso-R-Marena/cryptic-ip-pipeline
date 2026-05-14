from crypticip.statistics import (cohen_d, mann_whitney_u, roc_auc,
                                  weight_sensitivity, bootstrap_auc,
                                  permutation_p, leave_one_out_auc)


def test_cohen_d_basic():
    d = cohen_d([10, 11, 12], [1, 2, 3])
    assert d > 2.0


def test_mann_whitney_separable():
    u, p = mann_whitney_u([10, 11, 12, 13], [1, 2, 3, 4])
    assert p < 0.1


def test_roc_auc_perfect():
    assert roc_auc([0.9, 0.8, 0.7, 0.1, 0.2, 0.3], [1, 1, 1, 0, 0, 0]) == 1.0


def test_roc_auc_random():
    auc = roc_auc([0.5] * 10, [1] * 5 + [0] * 5)
    assert 0.0 <= auc <= 1.0


def test_weight_sensitivity_runs():
    rows = [
        {"label": 1, "depth": 0.8, "inv_sasa": 0.9, "elec": 0.6, "basic": 0.7, "volume_fit": 0.8, "plddt_penalty": 0.0},
        {"label": 1, "depth": 0.7, "inv_sasa": 0.85, "elec": 0.55, "basic": 0.6, "volume_fit": 0.8, "plddt_penalty": 0.0},
        {"label": 0, "depth": 0.2, "inv_sasa": 0.1, "elec": 0.05, "basic": 0.1, "volume_fit": 0.4, "plddt_penalty": 0.0},
        {"label": 0, "depth": 0.25, "inv_sasa": 0.15, "elec": 0.1, "basic": 0.15, "volume_fit": 0.4, "plddt_penalty": 0.0},
    ]
    out = weight_sensitivity(rows, n_samples=10, seed=0)
    assert len(out) == 10
    assert all(0 <= r.auc <= 1 for r in out)


def test_bootstrap_auc():
    mean, lo, hi = bootstrap_auc([0.9, 0.8, 0.7, 0.1, 0.2, 0.3],
                                 [1, 1, 1, 0, 0, 0], n_boot=200, seed=0)
    assert mean > 0.5
    assert lo <= mean <= hi


def test_permutation_p():
    p = permutation_p([0.9, 0.8, 0.7, 0.1, 0.2, 0.3], [1, 1, 1, 0, 0, 0],
                      n_perm=200, seed=0)
    assert p < 0.2


def test_leave_one_out():
    out = leave_one_out_auc([0.9, 0.8, 0.7, 0.1, 0.2, 0.3], [1, 1, 1, 0, 0, 0])
    assert len(out) == 6
