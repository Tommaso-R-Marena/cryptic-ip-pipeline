import json
import subprocess
import sys
from pathlib import Path

from crypticip.cli import build_parser, main


def test_help_runs():
    p = build_parser()
    out = p.format_help()
    assert "check-env" in out
    assert "validate" in out
    assert "screen" in out


def test_check_env_runs(capsys):
    rc = main(["check-env"])
    captured = capsys.readouterr()
    assert rc == 0
    assert "External tool check" in captured.out


def test_list_configs(capsys):
    rc = main(["list-configs"])
    captured = capsys.readouterr()
    assert rc == 0
    assert "default.yaml" in captured.out


def test_screen_with_empty_proteome(tmp_path, monkeypatch, capsys):
    from crypticip.paths import ProjectPaths
    # No proteome files -> should still produce a (mostly empty) summary, no crash.
    monkeypatch.setenv("CRYPTICIP_LOG_JSON", "")
    # Pass a tmp config that points data dir to tmp_path.
    cfg_path = tmp_path / "tcfg.yaml"
    cfg_path.write_text(
        "paths:\n"
        f"  data_dir: {tmp_path}/data\n"
        f"  proteomes_dir: {tmp_path}/data/proteomes\n"
        f"  results_dir: {tmp_path}/results\n"
        f"  reports_dir: {tmp_path}/results/reports\n"
        f"  screening_dir: {tmp_path}/results/screening\n"
        f"  experimental_dir: {tmp_path}/results/experimental\n"
    )
    (tmp_path / "data" / "proteomes" / "yeast").mkdir(parents=True)
    rc = main(["screen", "--organism", "yeast", "--config", str(cfg_path), "--workers", "1"])
    assert rc == 0
