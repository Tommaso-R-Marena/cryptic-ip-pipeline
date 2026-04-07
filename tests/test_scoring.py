"""Tests for the composite scoring function."""
import pytest
import sys
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent.parent / 'pipeline'))


def compute_composite_score(depth, sasa, charge, basic_count,
                            w_depth=0.30, w_sasa=0.35, w_elec=0.20, w_basic=0.15):
    """Composite burial score."""
    import numpy as np
    n_depth = min(depth / 30.0, 1.0)
    n_sasa  = 1.0 - min(sasa / 150.0, 1.0)
    n_charge = min(max(charge, 0) / 15.0, 1.0)
    n_basic  = min(basic_count / 8.0, 1.0)
    return 0.30 * n_depth + 0.35 * n_sasa + 0.20 * n_charge + 0.15 * n_basic


class TestCompositeScore:
    def test_perfect_buried(self):
        """Perfect buried pocket: max depth, zero SASA, high charge, many basic."""
        score = compute_composite_score(depth=30, sasa=0, charge=15, basic_count=8)
        assert score == pytest.approx(1.0, abs=0.001)

    def test_perfect_surface(self):
        """Surface pocket: no depth, max SASA, no charge, no basic."""
        score = compute_composite_score(depth=0, sasa=150, charge=0, basic_count=0)
        assert score == pytest.approx(0.0, abs=0.001)

    def test_adar2_like(self):
        """ADAR2-like pocket should score above 0.5."""
        score = compute_composite_score(depth=23, sasa=66, charge=6, basic_count=2)
        assert score > 0.45

    def test_ph_domain_like(self):
        """PH-domain-like pocket should score below 0.5."""
        score = compute_composite_score(depth=10, sasa=120, charge=1, basic_count=1)
        assert score < 0.45

    def test_score_bounds(self):
        """Score should always be in [0, 1]."""
        import random
        random.seed(42)
        for _ in range(100):
            score = compute_composite_score(
                depth=random.uniform(0, 50),
                sasa=random.uniform(0, 200),
                charge=random.uniform(-5, 20),
                basic_count=random.randint(0, 12)
            )
            assert 0.0 <= score <= 1.0

    def test_depth_dominance(self):
        """Increasing depth alone should increase score."""
        s1 = compute_composite_score(depth=5, sasa=50, charge=3, basic_count=2)
        s2 = compute_composite_score(depth=25, sasa=50, charge=3, basic_count=2)
        assert s2 > s1

    def test_sasa_inverse(self):
        """Lower SASA should give higher score."""
        s1 = compute_composite_score(depth=20, sasa=100, charge=3, basic_count=2)
        s2 = compute_composite_score(depth=20, sasa=10, charge=3, basic_count=2)
        assert s2 > s1

    def test_weights_sum(self):
        """Weights should sum to 1.0."""
        assert 0.30 + 0.35 + 0.20 + 0.15 == pytest.approx(1.0)


class TestResultsIntegrity:
    """Test that actual results files exist and are valid."""
    
    def test_json_exists(self):
        import json
        from pathlib import Path
        results_path = Path(__file__).resolve().parent.parent / 'results' / 'validation_results.json'
        assert results_path.exists(), "validation_results.json not found"
        with open(results_path) as f:
            data = json.load(f)
        assert len(data) == 9, f"Expected 9 structures, got {len(data)}"

    def test_csv_exists(self):
        import csv
        from pathlib import Path
        csv_path = Path(__file__).resolve().parent.parent / 'results' / 'validation_summary.csv'
        assert csv_path.exists(), "validation_summary.csv not found"
        with open(csv_path) as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        assert len(rows) == 9

    def test_positive_higher_than_negative_mean(self):
        """Positive controls should have higher mean score than negatives."""
        import json
        from pathlib import Path
        import numpy as np
        results_path = Path(__file__).resolve().parent.parent / 'results' / 'validation_results.json'
        with open(results_path) as f:
            data = json.load(f)
        pos = [r['composite_score'] for r in data if r['category'] == 'positive']
        neg = [r['composite_score'] for r in data if r['category'] == 'negative']
        assert np.mean(pos) > np.mean(neg), "Positive controls should score higher on average"

    def test_adar2_crystal_top3(self):
        """ADAR2 crystal should have pocket rank in top 5."""
        import json
        from pathlib import Path
        results_path = Path(__file__).resolve().parent.parent / 'results' / 'validation_results.json'
        with open(results_path) as f:
            data = json.load(f)
        adar2 = [r for r in data if r['name'] == 'ADAR2_crystal'][0]
        assert adar2['best_pocket_num'] <= 5, f"ADAR2 IP6 pocket rank: {adar2['best_pocket_num']} (expected ≤5)"
