#!/usr/bin/env python3
"""
Expanded test suite for cryptic IP binding site pipeline.
Covers: fpocket parsing, SASA, pLDDT, APBS/DX parsing, RMSD, charge analysis,
composite scoring, success criteria, and end-to-end validation.
"""
import os, sys, json, tempfile, shutil
import numpy as np
import pytest

# Add pipeline to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'pipeline'))
from expanded_analysis import (
    parse_fpocket_info, parse_fpocket_pdb, get_pocket_residues, get_pocket_center,
    get_ip_residue_centroid, compute_sasa, extract_plddt, pocket_plddt_check,
    read_dx_at_point, analyze_charge, compute_composite_score,
    compute_binding_region_rmsd, evaluate_success_criteria,
    cohens_d, welch_t, incomplete_beta,
    VALIDATION_SET, ADAR2_IP_RESNUMS, PROJ, DATA
)


# ─── Fixtures ───────────────────────────────────────────────────────────────

@pytest.fixture
def adar2_crystal_pdb():
    return DATA / 'pdb' / '1ZY7.pdb'

@pytest.fixture
def adar2_af_pdb():
    return DATA / 'alphafold' / 'AF-P78563-F1.pdb'

@pytest.fixture
def fpocket_info_file():
    return DATA / 'fpocket_results' / 'ADAR2_crystal_out' / 'ADAR2_crystal_info.txt'

@pytest.fixture
def sample_dx_file():
    """Create a minimal DX file for testing the parser."""
    content = """# OpenDX format
object 1 class gridpositions counts 3 3 3
origin 0.0 0.0 0.0
delta 1.0 0.0 0.0
delta 0.0 1.0 0.0
delta 0.0 0.0 1.0
object 2 class gridconnections counts 3 3 3
object 3 class array type double rank 0 items 27 data follows
1.0 2.0 3.0
4.0 5.0 6.0
7.0 8.0 9.0
10.0 11.0 12.0
13.0 14.0 15.0
16.0 17.0 18.0
19.0 20.0 21.0
22.0 23.0 24.0
25.0 26.0 27.0
attribute "dep" string "positions"
object "regular positions regular connections" class field
component "positions" value 1
component "connections" value 2
component "data" value 3
"""
    tmpdir = tempfile.mkdtemp()
    dx_path = os.path.join(tmpdir, 'test.dx')
    with open(dx_path, 'w') as f:
        f.write(content)
    yield dx_path
    shutil.rmtree(tmpdir)


@pytest.fixture
def results_json():
    path = PROJ / 'results' / 'expanded_validation_results.json'
    if path.exists():
        with open(path) as f:
            return json.load(f)
    return None


# ─── 1. fpocket parsing tests ──────────────────────────────────────────────

class TestFpocketParsing:
    def test_parse_info_returns_pockets(self, fpocket_info_file):
        if not fpocket_info_file.exists():
            pytest.skip("fpocket info file not available")
        pockets = parse_fpocket_info(fpocket_info_file)
        assert len(pockets) > 0, "Should parse at least one pocket"

    def test_parse_info_pocket_has_volume(self, fpocket_info_file):
        if not fpocket_info_file.exists():
            pytest.skip("fpocket info file not available")
        pockets = parse_fpocket_info(fpocket_info_file)
        assert pockets[0]['volume'] > 0, "Pocket 1 should have positive volume"

    def test_parse_info_pocket_has_depth(self, fpocket_info_file):
        if not fpocket_info_file.exists():
            pytest.skip("fpocket info file not available")
        pockets = parse_fpocket_info(fpocket_info_file)
        assert pockets[0]['pocket_depth'] > 0, "Pocket 1 should have positive depth"

    def test_adar2_pocket1_volume_near_expected(self, fpocket_info_file):
        if not fpocket_info_file.exists():
            pytest.skip("fpocket info file not available")
        pockets = parse_fpocket_info(fpocket_info_file)
        # ADAR2 pocket 1 volume is ~972 Å³ per fpocket
        assert 800 < pockets[0]['volume'] < 1200, \
            f"ADAR2 pocket 1 volume should be ~972 Å³, got {pockets[0]['volume']}"

    def test_adar2_pocket1_depth_near_expected(self, fpocket_info_file):
        if not fpocket_info_file.exists():
            pytest.skip("fpocket info file not available")
        pockets = parse_fpocket_info(fpocket_info_file)
        # Expected: Cent. of mass - Alpha Sphere max dist: 16.933
        assert 15 < pockets[0]['pocket_depth'] < 20, \
            f"ADAR2 pocket 1 depth should be ~16.9 Å, got {pockets[0]['pocket_depth']}"

    def test_parse_info_pocket_has_score(self, fpocket_info_file):
        if not fpocket_info_file.exists():
            pytest.skip("fpocket info file not available")
        pockets = parse_fpocket_info(fpocket_info_file)
        assert pockets[0]['score'] > 0, "Pocket 1 should have a positive score"

    def test_parse_info_pocket_has_druggability(self, fpocket_info_file):
        if not fpocket_info_file.exists():
            pytest.skip("fpocket info file not available")
        pockets = parse_fpocket_info(fpocket_info_file)
        assert pockets[0]['druggability'] >= 0, "Druggability should be non-negative"

    def test_pocket_residues_extraction(self):
        outdir = DATA / 'fpocket_results' / 'ADAR2_crystal_out'
        if not outdir.exists():
            pytest.skip("fpocket output not available")
        residues = get_pocket_residues(outdir, 1)
        assert len(residues) > 5, "Pocket 1 should have >5 residues"

    def test_pocket_center_computation(self):
        outdir = DATA / 'fpocket_results' / 'ADAR2_crystal_out'
        if not outdir.exists():
            pytest.skip("fpocket output not available")
        center = get_pocket_center(outdir, 1)
        assert center is not None, "Should compute pocket center"
        assert len(center) == 3, "Center should be 3D"


# ─── 2. IP residue centroid tests ──────────────────────────────────────────

class TestIPCentroid:
    def test_centroid_returns_3d_point(self, adar2_crystal_pdb):
        if not adar2_crystal_pdb.exists():
            pytest.skip("PDB not available")
        centroid = get_ip_residue_centroid(adar2_crystal_pdb, ADAR2_IP_RESNUMS)
        assert centroid is not None
        assert len(centroid) == 3

    def test_centroid_uses_first_chain_only(self, adar2_crystal_pdb):
        """1ZY7 is a homodimer; centroid should use one chain, not average both."""
        if not adar2_crystal_pdb.exists():
            pytest.skip("PDB not available")
        centroid = get_ip_residue_centroid(adar2_crystal_pdb, ADAR2_IP_RESNUMS)
        # If averaging both chains, centroid would be ~(57, 36, 40)
        # Chain A centroid should be ~(49, 22, 24) using sidechain atoms
        assert centroid[0] < 55, \
            f"Centroid x={centroid[0]:.1f} suggests both chains averaged (bug)"

    def test_centroid_returns_none_for_missing_residues(self, adar2_crystal_pdb):
        if not adar2_crystal_pdb.exists():
            pytest.skip("PDB not available")
        result = get_ip_residue_centroid(adar2_crystal_pdb, [9999, 9998, 9997])
        assert result is None


# ─── 3. SASA tests ─────────────────────────────────────────────────────────

class TestSASA:
    def test_sasa_returns_residue_dict(self, adar2_crystal_pdb):
        if not adar2_crystal_pdb.exists():
            pytest.skip("PDB not available")
        sasa = compute_sasa(adar2_crystal_pdb)
        assert len(sasa) > 100, "ADAR2 should have >100 residues"

    def test_sasa_has_sidechain_component(self, adar2_crystal_pdb):
        if not adar2_crystal_pdb.exists():
            pytest.skip("PDB not available")
        sasa = compute_sasa(adar2_crystal_pdb)
        first_key = next(iter(sasa))
        assert 'sidechain_sasa' in sasa[first_key]
        assert 'backbone_sasa' in sasa[first_key]
        assert 'total_sasa' in sasa[first_key]

    def test_sidechain_plus_backbone_equals_total(self, adar2_crystal_pdb):
        if not adar2_crystal_pdb.exists():
            pytest.skip("PDB not available")
        sasa = compute_sasa(adar2_crystal_pdb)
        for resnum, data in list(sasa.items())[:10]:
            total = data['sidechain_sasa'] + data['backbone_sasa']
            assert abs(total - data['total_sasa']) < 0.01, \
                f"Residue {resnum}: sidechain+backbone != total SASA"


# ─── 4. pLDDT tests ────────────────────────────────────────────────────────

class TestPLDDT:
    def test_plddt_extraction(self, adar2_af_pdb):
        if not adar2_af_pdb.exists():
            pytest.skip("AlphaFold PDB not available")
        plddt = extract_plddt(adar2_af_pdb)
        assert len(plddt) > 100, "Should extract pLDDT for >100 residues"

    def test_plddt_values_in_range(self, adar2_af_pdb):
        if not adar2_af_pdb.exists():
            pytest.skip("AlphaFold PDB not available")
        plddt = extract_plddt(adar2_af_pdb)
        for resnum, score in plddt.items():
            assert 0 <= score <= 100, f"pLDDT for residue {resnum} out of range: {score}"

    def test_plddt_check_passes_high_confidence(self):
        scores = {1: 90, 2: 85, 3: 95, 4: 80}
        result = pocket_plddt_check(scores, [1, 2, 3, 4])
        assert result['passes'] == True
        assert result['avg_plddt'] > 70

    def test_plddt_check_fails_low_confidence(self):
        scores = {1: 40, 2: 50, 3: 60, 4: 30}
        result = pocket_plddt_check(scores, [1, 2, 3, 4])
        assert result['passes'] == False
        assert result['avg_plddt'] < 70

    def test_plddt_check_handles_empty_input(self):
        result = pocket_plddt_check({}, [])
        assert result['passes'] is True
        assert result['avg_plddt'] is None


# ─── 5. DX parser tests ────────────────────────────────────────────────────

class TestDXParser:
    def test_read_dx_at_origin(self, sample_dx_file):
        """Value at origin (0,0,0) should be 1.0 (first data point)."""
        result = read_dx_at_point(sample_dx_file, np.array([0.0, 0.0, 0.0]))
        assert result is not None
        assert abs(result - 1.0) < 0.01

    def test_read_dx_at_center(self, sample_dx_file):
        """Value at (1,1,1) should be 14.0 (center of 3x3x3 grid)."""
        result = read_dx_at_point(sample_dx_file, np.array([1.0, 1.0, 1.0]))
        assert result is not None
        assert abs(result - 14.0) < 0.01

    def test_read_dx_interpolation(self, sample_dx_file):
        """Value at (0.5, 0.5, 0.5) should be between 1 and 14."""
        result = read_dx_at_point(sample_dx_file, np.array([0.5, 0.5, 0.5]))
        assert result is not None
        assert 1 < result < 14

    def test_read_dx_outside_grid(self, sample_dx_file):
        """Point outside grid should return None."""
        result = read_dx_at_point(sample_dx_file, np.array([10.0, 10.0, 10.0]))
        assert result is None

    def test_dx_parser_handles_trailing_text(self):
        """DX parser should stop at non-numeric lines like 'component'."""
        content = """object 1 class gridpositions counts 2 2 2
origin 0.0 0.0 0.0
delta 1.0 0.0 0.0
delta 0.0 1.0 0.0
delta 0.0 0.0 1.0
object 2 class gridconnections counts 2 2 2
object 3 class array type double rank 0 items 8 data follows
1.0 2.0 3.0 4.0
5.0 6.0 7.0 8.0
attribute "dep" string "positions"
component "positions" value 1
component "data" value 3
"""
        tmpdir = tempfile.mkdtemp()
        dx_path = os.path.join(tmpdir, 'test.dx')
        with open(dx_path, 'w') as f:
            f.write(content)
        result = read_dx_at_point(dx_path, np.array([0.0, 0.0, 0.0]))
        shutil.rmtree(tmpdir)
        assert result is not None, "Parser should handle trailing 'component' lines"
        assert abs(result - 1.0) < 0.01


# ─── 6. Charge analysis tests ──────────────────────────────────────────────

class TestChargeAnalysis:
    def test_charge_analysis_returns_radii(self, adar2_crystal_pdb):
        if not adar2_crystal_pdb.exists():
            pytest.skip("PDB not available")
        centroid = get_ip_residue_centroid(adar2_crystal_pdb, ADAR2_IP_RESNUMS)
        if centroid is None:
            pytest.skip("Centroid computation failed")
        charge = analyze_charge(adar2_crystal_pdb, centroid)
        assert '5.0A' in charge
        assert '8.0A' in charge
        assert '10.0A' in charge

    def test_charge_adar2_has_basic_residues_8A(self, adar2_crystal_pdb):
        if not adar2_crystal_pdb.exists():
            pytest.skip("PDB not available")
        centroid = get_ip_residue_centroid(adar2_crystal_pdb, ADAR2_IP_RESNUMS)
        if centroid is None:
            pytest.skip("Centroid computation failed")
        charge = analyze_charge(adar2_crystal_pdb, centroid)
        # Should find ≥6 basic residues within 8 Å (sidechain atoms)
        assert charge['8.0A']['basic_count'] >= 6, \
            f"Expected ≥6 basic residues within 8 Å, got {charge['8.0A']['basic_count']}"

    def test_charge_net_positive(self, adar2_crystal_pdb):
        if not adar2_crystal_pdb.exists():
            pytest.skip("PDB not available")
        centroid = get_ip_residue_centroid(adar2_crystal_pdb, ADAR2_IP_RESNUMS)
        if centroid is None:
            pytest.skip("Centroid computation failed")
        charge = analyze_charge(adar2_crystal_pdb, centroid)
        assert charge['8.0A']['net_charge'] > 0, \
            "Net charge near IP6 site should be positive"

    def test_charge_uses_sidechain_atoms(self, adar2_crystal_pdb):
        """Sidechain atom distances should differ from CA-based distances."""
        if not adar2_crystal_pdb.exists():
            pytest.skip("PDB not available")
        centroid = get_ip_residue_centroid(adar2_crystal_pdb, ADAR2_IP_RESNUMS)
        if centroid is None:
            pytest.skip("Centroid computation failed")
        charge = analyze_charge(adar2_crystal_pdb, centroid)
        # With sidechain atoms, basic residues should be found closer
        assert charge['8.0A']['basic_count'] > 0, \
            "Should find basic residues using sidechain atom distances"


# ─── 7. Composite scoring tests ────────────────────────────────────────────

class TestCompositeScore:
    def test_ideal_pocket_scores_high(self):
        """A deep, low-SASA, high-charge, good-volume pocket should score well."""
        result = compute_composite_score(
            depth=20, sasa=10, net_charge_8A=10, basic_5A=8,
            volume=500, apbs_potential=20
        )
        assert result['composite_score'] > 0.7

    def test_poor_pocket_scores_low(self):
        """A shallow, high-SASA, no-charge, bad-volume pocket should score poorly."""
        result = compute_composite_score(
            depth=2, sasa=100, net_charge_8A=0, basic_5A=0,
            volume=50, apbs_potential=-5
        )
        assert result['composite_score'] < 0.2

    def test_volume_in_ideal_range_scores_1(self):
        result = compute_composite_score(
            depth=15, sasa=30, net_charge_8A=5, basic_5A=4,
            volume=500, apbs_potential=5
        )
        assert result['score_volume'] == 1.0

    def test_volume_outside_range_penalized(self):
        result = compute_composite_score(
            depth=15, sasa=30, net_charge_8A=5, basic_5A=4,
            volume=2000, apbs_potential=5
        )
        assert result['score_volume'] < 0.5

    def test_score_components_sum_to_1_weight(self):
        result = compute_composite_score(
            depth=15, sasa=30, net_charge_8A=5, basic_5A=4,
            volume=500, apbs_potential=5
        )
        total_weight = sum(result['weights'].values())
        assert abs(total_weight - 1.0) < 0.01


# ─── 8. RMSD tests ─────────────────────────────────────────────────────────

class TestRMSD:
    def test_rmsd_computation(self, adar2_crystal_pdb, adar2_af_pdb):
        if not adar2_crystal_pdb.exists() or not adar2_af_pdb.exists():
            pytest.skip("PDB files not available")
        result = compute_binding_region_rmsd(
            adar2_crystal_pdb, adar2_af_pdb, ADAR2_IP_RESNUMS
        )
        assert result['rmsd'] is not None
        assert result['n_aligned'] >= 4

    def test_rmsd_is_large_for_adar2(self, adar2_crystal_pdb, adar2_af_pdb):
        """AlphaFold without IP6 predicts different conformation — expect large RMSD."""
        if not adar2_crystal_pdb.exists() or not adar2_af_pdb.exists():
            pytest.skip("PDB files not available")
        result = compute_binding_region_rmsd(
            adar2_crystal_pdb, adar2_af_pdb, ADAR2_IP_RESNUMS
        )
        assert result['rmsd'] > 5.0, \
            f"Expected large RMSD (>5 Å) for ADAR2 crystal vs AF, got {result['rmsd']}"

    def test_rmsd_self_comparison_near_zero(self, adar2_crystal_pdb):
        """Comparing a structure to itself should give ~0 RMSD."""
        if not adar2_crystal_pdb.exists():
            pytest.skip("PDB not available")
        result = compute_binding_region_rmsd(
            adar2_crystal_pdb, adar2_crystal_pdb, ADAR2_IP_RESNUMS
        )
        assert result['rmsd'] is not None
        assert result['rmsd'] < 0.1, \
            f"Self-comparison RMSD should be ~0, got {result['rmsd']}"

    def test_rmsd_reports_per_residue(self, adar2_crystal_pdb, adar2_af_pdb):
        if not adar2_crystal_pdb.exists() or not adar2_af_pdb.exists():
            pytest.skip("PDB files not available")
        result = compute_binding_region_rmsd(
            adar2_crystal_pdb, adar2_af_pdb, ADAR2_IP_RESNUMS
        )
        assert 'per_residue_distances' in result
        assert len(result['per_residue_distances']) > 0


# ─── 9. Statistics helper tests ─────────────────────────────────────────────

class TestStatistics:
    def test_cohens_d_identical_groups(self):
        d = cohens_d([1, 2, 3], [1, 2, 3])
        assert abs(d) < 0.01

    def test_cohens_d_different_groups(self):
        d = cohens_d([10, 11, 12], [1, 2, 3])
        assert d > 3.0  # Very large effect

    def test_welch_t_significant(self):
        t, p = welch_t([10, 11, 12, 13, 14], [1, 2, 3, 4, 5])
        assert p < 0.01

    def test_welch_t_nonsignificant(self):
        t, p = welch_t([5, 5, 5], [5, 5, 5])
        assert p > 0.5

    def test_incomplete_beta_bounds(self):
        assert incomplete_beta(1, 1, 0) == 0
        assert incomplete_beta(1, 1, 1) == 1


# ─── 10. End-to-end validation ──────────────────────────────────────────────

class TestEndToEnd:
    def test_results_json_exists(self):
        path = PROJ / 'results' / 'expanded_validation_results.json'
        assert path.exists(), "Results JSON should exist after analysis"

    def test_results_has_all_structures(self, results_json):
        if results_json is None:
            pytest.skip("Results not available")
        names = [r['name'] for r in results_json['validation_results']]
        assert 'ADAR2_crystal' in names
        assert 'ADAR2_alphafold' in names
        assert len(names) == 9

    def test_adar2_pocket_rank_1(self, results_json):
        if results_json is None:
            pytest.skip("Results not available")
        adar2 = next(r for r in results_json['validation_results'] if r['name'] == 'ADAR2_crystal')
        assert adar2['best_pocket_rank'] == 1

    def test_adar2_apbs_positive(self, results_json):
        if results_json is None:
            pytest.skip("Results not available")
        adar2 = next(r for r in results_json['validation_results'] if r['name'] == 'ADAR2_crystal')
        assert adar2['apbs_potential_kTe'] > 5, \
            f"ADAR2 APBS should be >5 kT/e, got {adar2['apbs_potential_kTe']}"

    def test_alphafold_plddt_passes(self, results_json):
        if results_json is None:
            pytest.skip("Results not available")
        af = next(r for r in results_json['validation_results'] if r['name'] == 'ADAR2_alphafold')
        # JSON serialization may convert numpy bool to string
        assert af['plddt_check']['passes'] in (True, 'True')

    def test_success_criteria_present(self, results_json):
        if results_json is None:
            pytest.skip("Results not available")
        assert len(results_json['success_criteria']) >= 5

    def test_positive_mean_depth_greater_than_negative(self, results_json):
        if results_json is None:
            pytest.skip("Results not available")
        results = results_json['validation_results']
        pos_depths = [r['pocket_depth'] for r in results if r['category'] == 'positive']
        neg_depths = [r['pocket_depth'] for r in results if r['category'] == 'negative']
        assert np.mean(pos_depths) > np.mean(neg_depths), \
            "Positive controls should have deeper pockets on average"


if __name__ == '__main__':
    pytest.main([__file__, '-v', '--tb=short'])
