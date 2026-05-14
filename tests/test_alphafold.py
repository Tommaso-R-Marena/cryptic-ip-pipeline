import io
import tarfile
from pathlib import Path

from crypticip.alphafold import (extract_proteome_tar, accession_from_filename,
                                  iter_proteome_files)


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
