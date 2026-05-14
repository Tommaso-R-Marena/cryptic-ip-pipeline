from crypticip.electrostatics import parse_dx, potential_at_center, run_apbs


def test_parse_dx_basic(dx_file):
    grid = parse_dx(dx_file)
    assert grid.shape == (3, 3, 3)
    assert grid.origin == (0.0, 0.0, 0.0)
    assert len(grid.data) == 27
    assert grid.data[0] == 1.0
    assert grid.data[-1] == 27.0


def test_parse_dx_trailing_lines(dx_file):
    # Should not include the "attribute"/"component" lines in data.
    grid = parse_dx(dx_file)
    assert len(grid.data) == 27


def test_potential_at_center(dx_file):
    # Centre of a 3x3x3 grid is index (1,1,1) → flat index 13 → value 14.0
    v = potential_at_center(dx_file, (1.0, 1.0, 1.0))
    assert v == 14.0


def test_run_apbs_fallback(tiny_pdb):
    result = run_apbs(tiny_pdb, (12.0, 17.0, 16.0))
    # On a machine without APBS we should get a fallback, not a crash.
    assert result["status"] in ("ok", "fallback", "failed")
    assert "potential_kT_e" in result
    if result["status"] == "fallback":
        assert result["backend"] == "coulomb"
