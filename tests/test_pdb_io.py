from crypticip.pdb_io import (parse_pdb_atoms, detect_ip_ligands,
                              clean_for_alphafold, clean_preserving_ip,
                              protein_atoms_only, chain_subset)


def test_parse_pdb(tiny_pdb):
    atoms = parse_pdb_atoms(tiny_pdb)
    assert len(atoms) == 36
    assert atoms[0].name == "N"
    assert atoms[0].resname == "ALA"
    assert atoms[0].chain == "A"


def test_detect_ip_ligands(tiny_pdb):
    atoms = parse_pdb_atoms(tiny_pdb)
    ligands = detect_ip_ligands(atoms)
    assert len(ligands) == 1
    assert ligands[0]["resname"] == "IHP"
    assert ligands[0]["ip_type"] == "IP6"


def test_clean_for_alphafold_drops_hetatm(tiny_pdb):
    atoms = parse_pdb_atoms(tiny_pdb)
    cleaned = clean_for_alphafold(atoms)
    assert not any(a.record == "HETATM" for a in cleaned)


def test_clean_preserving_ip_keeps_ip(tiny_pdb):
    atoms = parse_pdb_atoms(tiny_pdb)
    cleaned = clean_preserving_ip(atoms)
    hetatms = [a for a in cleaned if a.record == "HETATM"]
    assert all(a.resname == "IHP" for a in hetatms)
    assert hetatms, "should still have IP ligand"


def test_chain_subset(tiny_pdb):
    atoms = parse_pdb_atoms(tiny_pdb)
    a_only = chain_subset(atoms, "A")
    assert len(a_only) == len(atoms)
