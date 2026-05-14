from crypticip.sasa import compute_sasa, mean_sidechain_sasa


def test_sasa_handles_missing(tiny_pdb):
    res = compute_sasa(tiny_pdb)
    # In CI without freesasa installed we want status=missing, no exception.
    assert res.status in ("ok", "missing")
    if res.status == "missing":
        assert res.residues == []
        assert mean_sidechain_sasa(res, [("A", 1)]) is None


def test_sasa_returns_residue_list_when_available(tiny_pdb):
    res = compute_sasa(tiny_pdb)
    if res.status == "ok":
        assert res.residues
        m = mean_sidechain_sasa(res, [("A", 2)])
        assert m is None or m >= 0
