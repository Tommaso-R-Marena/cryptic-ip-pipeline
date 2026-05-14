import random

from crypticip.scoring import (FeatureVector, composite_score, filter_flags, tier)


def test_composite_monotonic_in_depth():
    base = FeatureVector(depth=5, sasa=50, elec=5, basic_count=3, volume=400, plddt=85)
    deep = FeatureVector(depth=25, sasa=50, elec=5, basic_count=3, volume=400, plddt=85)
    assert composite_score(deep).composite > composite_score(base).composite


def test_composite_inverse_sasa():
    a = FeatureVector(depth=20, sasa=100, elec=5, basic_count=3, volume=400, plddt=85)
    b = FeatureVector(depth=20, sasa=10, elec=5, basic_count=3, volume=400, plddt=85)
    assert composite_score(b).composite > composite_score(a).composite


def test_composite_bounded():
    rng = random.Random(0)
    for _ in range(200):
        fv = FeatureVector(
            depth=rng.uniform(0, 50),
            sasa=rng.uniform(0, 200),
            elec=rng.uniform(-5, 40),
            basic_count=rng.randint(0, 12),
            volume=rng.uniform(0, 1500),
            plddt=rng.uniform(40, 95),
        )
        s = composite_score(fv).composite
        assert 0.0 <= s <= 1.0


def test_filter_flags_basic_pass():
    fv = FeatureVector(depth=20, sasa=2, elec=10, basic_count=6, volume=500, plddt=85)
    flags = filter_flags(fv)
    assert all(flags.values())


def test_filter_flags_volume_cap():
    fv = FeatureVector(depth=20, sasa=2, elec=10, basic_count=6, volume=1500, plddt=85)
    flags = filter_flags(fv)
    assert flags["volume_ok"] is False


def test_plddt_penalty_applied():
    high = FeatureVector(depth=20, sasa=2, elec=10, basic_count=6, volume=500, plddt=85)
    low = FeatureVector(depth=20, sasa=2, elec=10, basic_count=6, volume=500, plddt=50)
    assert composite_score(low).composite < composite_score(high).composite


def test_tier_promotion_requires_apbs():
    fv = FeatureVector(depth=25, sasa=1, elec=15, basic_count=7, volume=550, plddt=85,
                       apbs_status="fallback")
    s = composite_score(fv).composite
    flags = filter_flags(fv)
    assert tier(s, flags, apbs_status="fallback") != "Tier1"
    assert tier(s, flags, apbs_status="ok") == "Tier1"
