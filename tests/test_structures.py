from crypticip.structures import preprocess_structure, per_residue_ca_plddt, plddt_summary
from crypticip.pdb_io import parse_pdb_atoms


def test_preprocess_alphafold(tiny_pdb, tmp_path):
    meta = preprocess_structure(tiny_pdb, out_dir=tmp_path, name="t",
                                is_alphafold=True)
    assert meta.cleaned_path.exists()
    assert meta.n_protein_atoms > 0
    assert meta.n_residues == 4
    assert meta.mean_plddt is not None
    # tiny.pdb CA b-factors: 82, 90, 75, 60 -> mean ≈ 76.75
    assert abs(meta.mean_plddt - 76.75) < 1e-6


def test_preprocess_crystal_keeps_ip(tiny_pdb, tmp_path):
    meta = preprocess_structure(tiny_pdb, out_dir=tmp_path, name="c",
                                is_alphafold=False)
    cleaned = parse_pdb_atoms(meta.cleaned_path)
    het = [a for a in cleaned if a.record == "HETATM"]
    assert any(a.resname == "IHP" for a in het)
    # waters should be removed
    assert not any(a.resname == "HOH" for a in cleaned)


def test_plddt_summary_empty():
    s = plddt_summary([])
    assert s["mean"] is None
    assert s["n"] == 0


def test_plddt_summary_basic():
    s = plddt_summary([90, 80, 70, 50, 30])
    assert s["mean"] == 64.0
    assert s["fraction_high"] == 3/5    # 90,80,70 >=70
    assert s["fraction_low"] == 1/5     # only 30 < 50
