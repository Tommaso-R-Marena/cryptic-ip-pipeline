"""Bootstrap helper for running the cryptic-IP notebooks on Google Colab.

The notebooks under ``notebooks/`` are designed to be opened directly in
Colab from GitHub and run top-to-bottom from a fresh runtime. This module
collects the idempotent setup steps that each notebook performs:

- clone or update the repo at a configurable branch
- pip install the project (editable) plus runtime deps
- optionally mount Google Drive and symlink ``results/`` onto Drive
- best-effort install of external scientific binaries via condacolab
- print versions and ``check-env`` so the user can confirm state

The helper is intentionally dependency-free (stdlib only at import time);
``google.colab`` is only imported inside the Drive-mount path.

Each notebook also embeds an inline copy of these steps so that opening a
notebook in Colab works even before this file has been fetched.
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

DEFAULT_REPO_URL = "https://github.com/Tommaso-R-Marena/cryptic-ip-pipeline.git"
DEFAULT_BRANCH = "main"
DEFAULT_PROJECT_DIR = "/content/cryptic-ip-pipeline"


def _run(cmd, check=True, **kwargs):
    """Run a shell command, echoing it first."""
    print("$", " ".join(cmd) if isinstance(cmd, list) else cmd)
    return subprocess.run(cmd, check=check, **kwargs)


def clone_or_update(
    repo_url: str = DEFAULT_REPO_URL,
    branch: str = DEFAULT_BRANCH,
    project_dir: str = DEFAULT_PROJECT_DIR,
) -> Path:
    """Clone ``repo_url`` at ``branch`` into ``project_dir`` (idempotent).

    If the directory already exists and is a git checkout, fetch + reset to
    ``origin/<branch>`` so re-running this cell never leaves a half-updated
    tree.
    """
    project = Path(project_dir)
    project.parent.mkdir(parents=True, exist_ok=True)
    if (project / ".git").exists():
        _run(["git", "-C", str(project), "fetch", "--quiet", "origin", branch])
        _run(["git", "-C", str(project), "checkout", branch])
        _run(["git", "-C", str(project), "reset", "--hard", f"origin/{branch}"])
    else:
        _run(["git", "clone", "--quiet", "--branch", branch, repo_url, str(project)])
    os.chdir(project)
    print("project_dir:", project)
    return project


def pip_install(project_dir: str = DEFAULT_PROJECT_DIR, dev: bool = False) -> None:
    """``pip install -e .`` (plus optional ``[dev]`` extras) for the repo."""
    spec = f"{project_dir}[dev]" if dev else project_dir
    _run([sys.executable, "-m", "pip", "install", "-q", "-e", spec])
    # nbformat is convenient for self-validating notebook checks.
    _run([sys.executable, "-m", "pip", "install", "-q", "nbformat"])


def mount_drive(
    drive_root: str = "/content/drive",
    drive_results: str | None = None,
    project_dir: str = DEFAULT_PROJECT_DIR,
) -> Path | None:
    """Mount Google Drive and (optionally) symlink ``<project>/results`` onto it.

    ``drive_results`` is the absolute path under ``/content/drive/MyDrive/...``
    where outputs should live. If provided, the helper ensures that path
    exists, then replaces ``<project>/results`` with a symlink pointing
    there. Pre-existing real ``results/`` directories are backed up to a
    timestamped sibling rather than deleted.
    """
    try:
        from google.colab import drive  # type: ignore[import-not-found]
    except ImportError:
        print("not running on Colab; skipping drive.mount()")
        return None
    drive.mount(drive_root)
    if drive_results is None:
        return Path(drive_root)
    target = Path(project_dir) / "results"
    drive_path = Path(drive_results)
    drive_path.mkdir(parents=True, exist_ok=True)
    if target.is_symlink():
        target.unlink()
    elif target.exists():
        import datetime, shutil
        backup = target.with_name(
            f"results.local_backup_{datetime.datetime.now():%Y%m%d_%H%M%S}"
        )
        shutil.move(str(target), backup)
        print("existing results/ backed up to:", backup)
    target.symlink_to(drive_path, target_is_directory=True)
    print("results/ ->", drive_path)
    return drive_path


def verify_cli() -> None:
    """Print ``crypticip --version``, ``check-env`` and ``list-configs``."""
    for cmd in (["crypticip", "--version"],
                ["crypticip", "check-env"],
                ["crypticip", "list-configs"]):
        try:
            _run(cmd, check=False)
        except FileNotFoundError:
            print("crypticip not on PATH yet — did pip install fail?")
            raise


def try_install_external_tools(use_condacolab: bool = False) -> None:
    """Best-effort install of fpocket / freesasa / pdb2pqr / apbs / pymol.

    Vanilla Colab does **not** ship these. The only reliable path is
    ``condacolab + mamba``, which **restarts the kernel** mid-cell — so
    this function is opt-in via ``use_condacolab=True`` and the rest of
    the notebook must re-bootstrap after the restart.
    """
    if not use_condacolab:
        print("skipping external-tool install (use_condacolab=False).")
        print("CLI will run in fallback mode; see crypticip check-env.")
        return
    _run([sys.executable, "-m", "pip", "install", "-q", "condacolab"])
    print("Now run, in a NEW cell:")
    print("    import condacolab; condacolab.install()")
    print("…wait for the kernel to restart, then run the post-restart cell.")


def bootstrap(
    repo_url: str = DEFAULT_REPO_URL,
    branch: str = DEFAULT_BRANCH,
    project_dir: str = DEFAULT_PROJECT_DIR,
    mount_drive_flag: bool = False,
    drive_root: str = "/content/drive",
    drive_results: str | None = None,
    dev: bool = False,
) -> Path:
    """One-shot bootstrap: clone + install + optional Drive mount + verify."""
    project = clone_or_update(repo_url, branch, project_dir)
    pip_install(str(project), dev=dev)
    if mount_drive_flag:
        mount_drive(drive_root, drive_results, str(project))
    verify_cli()
    return project
