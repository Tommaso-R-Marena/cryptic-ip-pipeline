import json
from pathlib import Path

from crypticip.config import load_config
from crypticip.paths import ProjectPaths
from crypticip.qc import qc_one_file, spot_check_proteome
from crypticip.screening import screen_proteome


def _make_mini_proteome(dir_: Path, src_pdb: Path, n: int = 3) -> None:
    dir_.mkdir(parents=True, exist_ok=True)
    text = src_pdb.read_text()
    for i in range(n):
        acc = f"P0000{i+1}"
        (dir_ / f"AF-{acc}-F1-model_v4.pdb").write_text(text)


def test_qc_one_file(tiny_pdb):
    r = qc_one_file(tiny_pdb)
    assert r.status in ("ok", "truncated")
    assert r.n_atoms > 0
    assert r.n_residues > 0


def test_qc_spot_check(tmp_path, tiny_pdb):
    p = tmp_path / "prot"
    _make_mini_proteome(p, tiny_pdb, n=5)
    report = spot_check_proteome(p, n=5)
    assert report["status"] == "ok"
    assert report["n_files_total"] == 5
    assert report["n_pass"] >= 4


def test_screen_mini_proteome(tmp_path, tiny_pdb, monkeypatch):
    cfg_path = tmp_path / "cfg.yaml"
    cfg_path.write_text(
        "paths:\n"
        f"  data_dir: {tmp_path}/data\n"
        f"  proteomes_dir: {tmp_path}/data/proteomes\n"
        f"  results_dir: {tmp_path}/results\n"
        f"  reports_dir: {tmp_path}/results/reports\n"
        f"  screening_dir: {tmp_path}/results/screening\n"
        f"  experimental_dir: {tmp_path}/results/experimental\n"
    )
    organism_dir = tmp_path / "data" / "proteomes" / "yeast"
    _make_mini_proteome(organism_dir, tiny_pdb, n=2)

    cfg = load_config(str(cfg_path))
    paths = ProjectPaths.from_config(cfg).ensure()
    summary = screen_proteome("yeast", cfg=cfg, paths=paths, workers=1, limit=2)
    # Without fpocket installed every protein records a missing pocket count,
    # but the pipeline must still complete and emit aggregate output.
    assert summary["n_processed_this_run"] == 2
    out_dir = paths.organism_screening_dir("yeast")
    assert (out_dir / "summary.json").exists()
    assert (out_dir / "screening_results.csv").exists()
