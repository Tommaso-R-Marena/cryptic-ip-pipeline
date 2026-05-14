"""Minimal PDB IO helpers — pure Python, no Biopython hard dependency."""
from __future__ import annotations

import gzip
import io
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

# IP-family residue names (matched in PDB HETATM lines).
IP_LIGAND_RESNAMES = {"IHP", "IHE", "IPT", "I3P", "I4P", "I5P", "I6P",
                      "INS", "I7P", "IP6", "IP4", "IP3"}

# Map PDB resname -> human-readable IP type.
IP_LIGAND_TO_TYPE = {
    "IHP": "IP6", "IP6": "IP6",
    "I4P": "IP4", "IP4": "IP4",
    "I3P": "IP3", "IP3": "IP3",
    "I5P": "IP5",
    "INS": "INS",
}


@dataclass
class Atom:
    record: str         # "ATOM" or "HETATM"
    serial: int
    name: str
    altloc: str
    resname: str
    chain: str
    resseq: int
    icode: str
    x: float
    y: float
    z: float
    occupancy: float
    bfactor: float
    element: str
    raw: str

    @property
    def is_water(self) -> bool:
        return self.resname.strip() in ("HOH", "WAT", "DOD")

    @property
    def is_ion_or_buffer(self) -> bool:
        return self.resname.strip() in {"SO4", "PO4", "CL", "NA", "K", "MG",
                                        "CA", "ZN", "FE", "MN", "CU", "EDO",
                                        "GOL", "PEG", "ACT", "DMS", "TRS",
                                        "HEPES", "BME", "FMT", "IPA", "MES",
                                        "BES", "CIT"}


def _open(path: Path | str) -> io.TextIOBase:
    p = Path(path)
    if str(p).endswith(".gz"):
        return io.TextIOWrapper(gzip.open(p, "rb"), encoding="utf-8", errors="replace")
    return p.open("r", encoding="utf-8", errors="replace")


def parse_pdb_atoms(path: Path | str) -> list[Atom]:
    """Parse all ATOM / HETATM records in a PDB (gz ok). Lenient: bad
    lines are skipped."""
    atoms: list[Atom] = []
    with _open(path) as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            try:
                a = Atom(
                    record=rec,
                    serial=int(line[6:11]) if line[6:11].strip() else 0,
                    name=line[12:16].strip(),
                    altloc=line[16].strip(),
                    resname=line[17:20].strip(),
                    chain=line[21:22].strip() or "A",
                    resseq=int(line[22:26]) if line[22:26].strip() else 0,
                    icode=line[26].strip(),
                    x=float(line[30:38]),
                    y=float(line[38:46]),
                    z=float(line[46:54]),
                    occupancy=float(line[54:60]) if line[54:60].strip() else 1.0,
                    bfactor=float(line[60:66]) if line[60:66].strip() else 0.0,
                    element=line[76:78].strip(),
                    raw=line.rstrip("\n"),
                )
            except (ValueError, IndexError):
                continue
            atoms.append(a)
    return atoms


def iter_chains(atoms: Iterable[Atom]) -> dict[str, list[Atom]]:
    out: dict[str, list[Atom]] = {}
    for a in atoms:
        out.setdefault(a.chain, []).append(a)
    return out


def detect_ip_ligands(atoms: Iterable[Atom]) -> list[dict]:
    """Return one record per IP molecule found in HETATM (chain, resseq, resname, ip_type, n_atoms)."""
    seen: dict[tuple[str, int, str], list[Atom]] = {}
    for a in atoms:
        if a.record == "HETATM" and a.resname.strip() in IP_LIGAND_RESNAMES:
            seen.setdefault((a.chain, a.resseq, a.resname.strip()), []).append(a)
    out = []
    for (chain, resseq, resname), group in seen.items():
        out.append({
            "chain": chain,
            "resseq": resseq,
            "resname": resname,
            "ip_type": IP_LIGAND_TO_TYPE.get(resname, resname),
            "n_atoms": len(group),
        })
    return out


def write_atoms(path: Path | str, atoms: Iterable[Atom]) -> None:
    """Write atoms as PDB; appends an END record."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w") as fh:
        for a in atoms:
            fh.write(a.raw + "\n")
        fh.write("END\n")


def clean_for_alphafold(atoms: Iterable[Atom]) -> list[Atom]:
    """Drop all HETATM records (waters, ligands, ions, buffers). AlphaFold
    structures should never contain ligands; if we are reusing a crystal
    PDB this strip is required for fpocket to behave like it would on AF."""
    return [a for a in atoms if a.record != "HETATM"]


def clean_preserving_ip(atoms: Iterable[Atom]) -> list[Atom]:
    """Drop waters/buffers/ions but KEEP IP-family ligands (for crystals)."""
    out: list[Atom] = []
    for a in atoms:
        if a.is_water:
            continue
        if a.record == "HETATM":
            rn = a.resname.strip()
            if rn in IP_LIGAND_RESNAMES:
                out.append(a)
                continue
            if a.is_ion_or_buffer:
                continue
        out.append(a)
    return out


def protein_atoms_only(atoms: Iterable[Atom]) -> list[Atom]:
    return [a for a in atoms if a.record == "ATOM"]


def chain_subset(atoms: Iterable[Atom], chain: str) -> list[Atom]:
    return [a for a in atoms if a.chain == chain]
