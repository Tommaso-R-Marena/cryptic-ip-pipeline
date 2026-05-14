from pathlib import Path

from crypticip.pymol import write_pml
from crypticip.experimental import build_experimental_plan
from crypticip.config import load_config
from crypticip.paths import ProjectPaths


def test_write_pml(tmp_path):
    out = tmp_path / "x.pml"
    write_pml(out, pdb_path=Path("/tmp/x.pdb"), accession="P00001",
              pocket_residues=[10, 20, 30], center=(1.0, 2.0, 3.0),
              is_alphafold=True)
    content = out.read_text()
    assert "load /tmp/x.pdb, prot" in content
    assert "pseudoatom pocket_center" in content
    assert "spectrum b" in content


def test_experimental_plan_no_results(tmp_path):
    cfg = load_config()
    cfg["paths"] = {
        "data_dir": str(tmp_path / "data"),
        "proteomes_dir": str(tmp_path / "data" / "proteomes"),
        "results_dir": str(tmp_path / "results"),
        "reports_dir": str(tmp_path / "results" / "reports"),
        "screening_dir": str(tmp_path / "results" / "screening"),
        "experimental_dir": str(tmp_path / "results" / "experimental"),
    }
    paths = ProjectPaths.from_config(cfg).ensure()
    result = build_experimental_plan("yeast", paths=paths, top_n=10)
    assert result["status"] == "no_results"
