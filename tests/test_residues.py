from crypticip.pdb_io import parse_pdb_atoms
from crypticip.residues import (sidechain_terminal_atom, residue_neighborhood,
                                multi_shell_summary, count_basic_within,
                                residue_to_terminal_distance, FORMAL_CHARGE)


def test_sidechain_terminal_lys(tiny_pdb):
    atoms = parse_pdb_atoms(tiny_pdb)
    lys = [a for a in atoms if a.resname == "LYS"]
    nz = sidechain_terminal_atom("LYS", lys)
    assert nz.name == "NZ"


def test_sidechain_terminal_arg(tiny_pdb):
    atoms = parse_pdb_atoms(tiny_pdb)
    arg = [a for a in atoms if a.resname == "ARG"]
    t = sidechain_terminal_atom("ARG", arg)
    # CZ is preferred, NH1 fallback — both should resolve
    assert t.name in ("CZ", "NH1", "NH2")


def test_terminal_distance(tiny_pdb):
    atoms = parse_pdb_atoms(tiny_pdb)
    lys = [a for a in atoms if a.resname == "LYS"]
    d = residue_to_terminal_distance((12.0, 17.0, 16.0), lys)
    # NZ at (12,17,16) is the terminal — distance should be 0
    assert d == 0.0


def test_count_basic_within_radius(tiny_pdb):
    atoms = parse_pdb_atoms(tiny_pdb)
    # Centre around the IP-equivalent point
    n = count_basic_within(atoms, (12.0, 17.0, 17.5), radius_A=5.0)
    # Both LYS(NZ@17,16) and ARG(CZ/NH@17/18,19) are nearby
    assert n >= 1


def test_neighborhood_charge(tiny_pdb):
    atoms = parse_pdb_atoms(tiny_pdb)
    nh = residue_neighborhood(atoms, (12.0, 16.0, 17.0), radius_A=4.0)
    # LYS contributes +1; ARG +1; ASP -1 may be out of range
    assert nh.net_charge >= 1
    assert nh.n_basic >= 1


def test_multi_shell(tiny_pdb):
    atoms = parse_pdb_atoms(tiny_pdb)
    s = multi_shell_summary(atoms, (12.0, 16.0, 17.0), radii_A=(5.0, 8.0))
    assert "r5A" in s and "r8A" in s
    assert s["r8A"]["n_basic"] >= s["r5A"]["n_basic"]


def test_formal_charges_table():
    assert FORMAL_CHARGE["ARG"] == 1.0
    assert FORMAL_CHARGE["ASP"] == -1.0
    assert FORMAL_CHARGE["HIS"] == 0.5
