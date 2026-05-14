import io
import tarfile
from pathlib import Path

import crypticip.alphafold as af
from crypticip.alphafold import (extract_proteome_tar, accession_from_filename,
                                  iter_proteome_files, stream_download)


def _mk_tar(path: Path, names_to_payload: dict) -> Path:
    with tarfile.open(path, "w") as tar:
        for name, payload in names_to_payload.items():
            data = payload.encode() if isinstance(payload, str) else payload
            info = tarfile.TarInfo(name=name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    return path


def test_extract_proteome_tar(tmp_path, tiny_pdb):
    pdb = tiny_pdb.read_text()
    tar = _mk_tar(tmp_path / "tiny.tar", {
        "AF-P00001-F1-model_v4.pdb": pdb,
        "AF-P00002-F1-model_v4.pdb": pdb,
        "AF-P00003-F1-model_v4.pdb": pdb,
    })
    out = tmp_path / "out"
    report = extract_proteome_tar(tar, out)
    assert report["n_extracted"] == 3
    assert report["n_files"] == 3
    assert report["n_unique_accessions"] == 3
    files = list(iter_proteome_files(out))
    assert len(files) == 3
    assert accession_from_filename(files[0]) is not None


class _FakeResp:
    def __init__(self, *, status_code: int, body: bytes, headers: dict | None = None):
        self.status_code = status_code
        self._body = body
        self.headers = headers or {}

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")

    def iter_content(self, chunk_size: int = 1):
        for i in range(0, len(self._body), chunk_size):
            yield self._body[i:i + chunk_size]


def test_stream_download_resume_non_206_overwrites(tmp_path, monkeypatch):
    """If a partial file exists and the server returns 200 (ignoring Range),
    the existing prefix must be discarded — not appended to."""
    dest = tmp_path / "blob.bin"
    dest.write_bytes(b"OLDPARTIAL")  # 10 bytes of stale partial data

    full_body = b"FULLPAYLOAD-1234567890"
    captured: dict = {}

    def fake_get(url, stream=True, timeout=None, headers=None):
        captured["headers"] = headers or {}
        # Server ignored Range and returned full body with 200 OK.
        return _FakeResp(status_code=200, body=full_body,
                         headers={"Content-Length": str(len(full_body))})

    monkeypatch.setattr(af.requests, "get", fake_get)

    res = stream_download("http://example/blob", dest, progress=False, chunk=4)
    assert res.ok
    assert captured["headers"].get("Range") == "bytes=10-"
    # File must equal the full payload — not OLDPARTIAL + full_body.
    assert dest.read_bytes() == full_body
    assert res.bytes == len(full_body)


def test_stream_download_resume_206_appends(tmp_path, monkeypatch):
    """If the server returns 206 Partial Content, append the body to existing partial."""
    dest = tmp_path / "blob.bin"
    existing = b"OLDPARTIAL"
    dest.write_bytes(existing)

    tail = b"-CONTINUED-TAIL"

    def fake_get(url, stream=True, timeout=None, headers=None):
        assert (headers or {}).get("Range") == f"bytes={len(existing)}-"
        return _FakeResp(
            status_code=206, body=tail,
            headers={"Content-Length": str(len(tail)),
                     "Content-Range": f"bytes {len(existing)}-{len(existing) + len(tail) - 1}/*"},
        )

    monkeypatch.setattr(af.requests, "get", fake_get)

    res = stream_download("http://example/blob", dest, progress=False, chunk=4)
    assert res.ok
    assert dest.read_bytes() == existing + tail
    assert res.bytes == len(existing) + len(tail)


def test_stream_download_resume_206_wrong_offset_restarts(tmp_path, monkeypatch):
    """If 206 Content-Range starts at a different offset than requested, restart cleanly."""
    dest = tmp_path / "blob.bin"
    dest.write_bytes(b"OLDPARTIAL")  # 10 bytes

    body = b"BRAND-NEW-FULL-PAYLOAD"

    def fake_get(url, stream=True, timeout=None, headers=None):
        return _FakeResp(
            status_code=206, body=body,
            headers={"Content-Length": str(len(body)),
                     "Content-Range": f"bytes 0-{len(body) - 1}/{len(body)}"},
        )

    monkeypatch.setattr(af.requests, "get", fake_get)
    res = stream_download("http://example/blob", dest, progress=False, chunk=4)
    assert res.ok
    assert dest.read_bytes() == body


def test_stream_download_416_returns_skipped(tmp_path, monkeypatch):
    dest = tmp_path / "blob.bin"
    dest.write_bytes(b"ALREADYCOMPLETE")

    def fake_get(url, stream=True, timeout=None, headers=None):
        return _FakeResp(status_code=416, body=b"")

    monkeypatch.setattr(af.requests, "get", fake_get)
    res = stream_download("http://example/blob", dest, progress=False)
    assert res.ok and res.skipped
    assert res.bytes == len(b"ALREADYCOMPLETE")


def test_stream_download_fresh_200(tmp_path, monkeypatch):
    """No existing file -> no Range header, 200 OK writes full body."""
    dest = tmp_path / "blob.bin"
    body = b"hello-world-payload"

    captured: dict = {}

    def fake_get(url, stream=True, timeout=None, headers=None):
        captured["headers"] = headers or {}
        return _FakeResp(status_code=200, body=body,
                         headers={"Content-Length": str(len(body))})

    monkeypatch.setattr(af.requests, "get", fake_get)
    res = stream_download("http://example/blob", dest, progress=False, chunk=5)
    assert res.ok
    assert "Range" not in captured["headers"]
    assert dest.read_bytes() == body
    assert res.bytes == len(body)
