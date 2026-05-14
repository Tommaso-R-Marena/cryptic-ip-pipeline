"""AlphaFold DB download helpers (per-protein and full proteomes)."""
from __future__ import annotations

import gzip
import io
import os
import re
import shutil
import tarfile
import time
import urllib.parse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

import requests

from .logging_utils import get_logger

log = get_logger(__name__)

AF_PROTEIN_URL = "https://alphafold.ebi.ac.uk/files/{accession}-model_v4.pdb"


@dataclass
class DownloadResult:
    path: Path
    ok: bool
    bytes: int
    skipped: bool = False
    error: str | None = None


def download_alphafold_pdb(accession: str, out_dir: Path, *, force: bool = False,
                           timeout: float = 60.0) -> DownloadResult:
    """Download AF-<UniProt>-F1-model_v4.pdb from EBI."""
    out_dir.mkdir(parents=True, exist_ok=True)
    target = out_dir / f"{accession}.pdb"
    if target.exists() and not force and target.stat().st_size > 0:
        return DownloadResult(target, ok=True, bytes=target.stat().st_size, skipped=True)
    url = AF_PROTEIN_URL.format(accession=accession)
    try:
        r = requests.get(url, timeout=timeout)
        r.raise_for_status()
    except Exception as e:
        return DownloadResult(target, ok=False, bytes=0, error=str(e))
    target.write_bytes(r.content)
    return DownloadResult(target, ok=True, bytes=len(r.content))


def stream_download(url: str, dest: Path, *, resume: bool = True,
                    chunk: int = 1 << 20, timeout: float = 300.0,
                    progress: bool = True) -> DownloadResult:
    """Streamed download with optional HTTP Range resume.

    If a Range request is issued but the server responds with 200 OK (full
    payload) instead of 206 Partial Content, the existing partial file is
    discarded and the download restarts cleanly to avoid corruption.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    headers: dict[str, str] = {}
    mode = "wb"
    existing = 0
    ranged = False
    if resume and dest.exists():
        existing = dest.stat().st_size
        if existing > 0:
            headers["Range"] = f"bytes={existing}-"
            mode = "ab"
            ranged = True
    try:
        with requests.get(url, stream=True, timeout=timeout, headers=headers) as r:
            if r.status_code == 416:
                return DownloadResult(dest, ok=True, bytes=existing, skipped=True)
            r.raise_for_status()
            done = existing
            if ranged:
                if r.status_code == 206:
                    # Sanity-check Content-Range start matches our offset, if present.
                    cr = r.headers.get("Content-Range", "")
                    if cr:
                        m = re.match(r"bytes\s+(\d+)-", cr)
                        if m and int(m.group(1)) != existing:
                            # Server resumed from a different offset; restart cleanly.
                            mode = "wb"
                            done = 0
                else:
                    # Server ignored Range and is returning the full payload (or other
                    # non-206 success). Discard partial file and restart cleanly so we
                    # don't append the full body onto our existing prefix.
                    mode = "wb"
                    done = 0
            with dest.open(mode) as fh:
                last = time.monotonic()
                for buf in r.iter_content(chunk_size=chunk):
                    if not buf:
                        continue
                    fh.write(buf)
                    done += len(buf)
                    if progress and time.monotonic() - last > 5.0:
                        log.info("download %s: %d MB", dest.name, done >> 20)
                        last = time.monotonic()
    except Exception as e:
        return DownloadResult(dest, ok=False, bytes=dest.stat().st_size if dest.exists() else 0, error=str(e))
    return DownloadResult(dest, ok=True, bytes=dest.stat().st_size)


_ACC_RE = re.compile(r"AF-([A-Z0-9]+)-F\d+-model_v\d+\.pdb(?:\.gz)?$")


def extract_proteome_tar(tarpath: Path, out_dir: Path, *, decompress_gz: bool = True,
                         limit: int | None = None) -> dict:
    """Extract an AlphaFold proteome tar (e.g. UP000002311_*.tar) into ``out_dir``.

    Files are placed directly into ``out_dir`` (no nested subdirs) and named
    by UniProt accession (``AF-<acc>-F1-model_v4.pdb[.gz]``).
    Returns an extraction report.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    n_total = 0
    n_extracted = 0
    n_skipped = 0
    accs: list[str] = []
    with tarfile.open(tarpath, "r:*") as tar:
        for member in tar:
            if limit is not None and n_extracted >= limit:
                break
            n_total += 1
            if not member.isfile():
                continue
            name = Path(member.name).name
            if not name.startswith("AF-") or "model" not in name:
                continue
            target = out_dir / name
            if target.exists() and target.stat().st_size > 0:
                n_skipped += 1
                m = _ACC_RE.search(name)
                if m:
                    accs.append(m.group(1))
                continue
            f = tar.extractfile(member)
            if f is None:
                continue
            data = f.read()
            if not data:
                continue
            if decompress_gz and name.endswith(".gz"):
                try:
                    data = gzip.decompress(data)
                    name = name[:-3]
                    target = out_dir / name
                except OSError:
                    pass
            target.write_bytes(data)
            n_extracted += 1
            m = _ACC_RE.search(member.name)
            if m:
                accs.append(m.group(1))

    return {
        "tar": str(tarpath),
        "out_dir": str(out_dir),
        "n_tar_entries": n_total,
        "n_extracted": n_extracted,
        "n_skipped_existing": n_skipped,
        "n_files": len(list(out_dir.glob("AF-*.pdb"))),
        "n_unique_accessions": len(set(accs)),
    }


def iter_proteome_files(out_dir: Path) -> Iterator[Path]:
    yield from sorted(out_dir.glob("AF-*.pdb"))
    yield from sorted(out_dir.glob("AF-*.pdb.gz"))


def accession_from_filename(p: Path) -> str | None:
    m = _ACC_RE.search(p.name)
    return m.group(1) if m else None
