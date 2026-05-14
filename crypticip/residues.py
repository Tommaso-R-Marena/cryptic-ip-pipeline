"""Residue chemistry helpers: charges, sidechain terminal atoms, distance counts."""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Iterable

from .pdb_io import Atom


FORMAL_CHARGE: dict[str, float] = {
    "ARG": +1.0, "LYS": +1.0, "HIS": +0.5,
    "ASP": -1.0, "GLU": -1.0,
}
BASIC_RESIDUES = {"ARG", "LYS", "HIS"}
ACIDIC_RESIDUES = {"ASP", "GLU"}
AROMATIC_RESIDUES = {"PHE", "TYR", "TRP", "HIS"}

# Side-chain terminal atom used for distance / electrostatics measurements.
SIDECHAIN_TERMINAL: dict[str, tuple[str, ...]] = {
    "LYS": ("NZ",),
    "ARG": ("CZ", "NH1", "NH2"),
    "HIS": ("NE2", "ND1"),
    "ASP": ("CG", "OD1", "OD2"),
    "GLU": ("CD", "OE1", "OE2"),
    "ASN": ("CG", "OD1", "ND2"),
    "GLN": ("CD", "OE1", "NE2"),
    "TRP": ("NE1", "CZ2"),
    "TYR": ("OH",),
    "SER": ("OG",),
    "THR": ("OG1",),
    "CYS": ("SG",),
    "MET": ("SD",),
    "PHE": ("CZ",),
    "GLY": ("CA",),
    "ALA": ("CB",),
    "VAL": ("CB",),
    "LEU": ("CG",),
    "ILE": ("CG1",),
    "PRO": ("CG",),
}

# Kyte–Doolittle hydropathy.
HYDROPATHY: dict[str, float] = {
    "ILE": 4.5, "VAL": 4.2, "LEU": 3.8, "PHE": 2.8, "CYS": 2.5, "MET": 1.9,
    "ALA": 1.8, "GLY": -0.4, "THR": -0.7, "SER": -0.8, "TRP": -0.9,
    "TYR": -1.3, "PRO": -1.6, "HIS": -3.2, "GLU": -3.5, "GLN": -3.5,
    "ASP": -3.5, "ASN": -3.5, "LYS": -3.9, "ARG": -4.5,
}

HBOND_DONOR_RESIDUES = {"ARG", "LYS", "ASN", "GLN", "SER", "THR", "TRP", "TYR", "HIS"}
HBOND_ACCEPTOR_RESIDUES = {"ASP", "GLU", "ASN", "GLN", "SER", "THR", "TYR", "HIS"}


@dataclass
class ResidueNeighborhood:
    center: tuple[float, float, float]
    radius_A: float
    n_basic: int = 0
    n_acidic: int = 0
    n_aromatic: int = 0
    net_charge: float = 0.0
    mean_hydropathy: float | None = None
    n_hbond_donors: int = 0
    n_hbond_acceptors: int = 0
    residues: list[tuple[str, int, str, float]] = field(default_factory=list)  # (chain, resseq, resname, dist)


def sidechain_terminal_atom(resname: str, atoms: Iterable[Atom]) -> Atom | None:
    """Return the terminal sidechain atom (NZ for LYS, etc.); falls back to CA."""
    candidates = SIDECHAIN_TERMINAL.get(resname.upper(), ("CA",))
    by_name = {a.name: a for a in atoms if a.record == "ATOM"}
    for name in candidates:
        if name in by_name:
            return by_name[name]
    return by_name.get("CA")


def residue_to_terminal_distance(center: tuple[float, float, float],
                                 residue_atoms: list[Atom]) -> float | None:
    if not residue_atoms:
        return None
    resname = residue_atoms[0].resname
    a = sidechain_terminal_atom(resname, residue_atoms)
    if a is None:
        return None
    dx = a.x - center[0]; dy = a.y - center[1]; dz = a.z - center[2]
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def _group_by_residue(atoms: Iterable[Atom]) -> dict[tuple[str, int, str], list[Atom]]:
    groups: dict[tuple[str, int, str], list[Atom]] = {}
    for a in atoms:
        if a.record != "ATOM":
            continue
        groups.setdefault((a.chain, a.resseq, a.resname), []).append(a)
    return groups


def residue_neighborhood(atoms: Iterable[Atom], center: tuple[float, float, float],
                         radius_A: float = 5.0) -> ResidueNeighborhood:
    """Summarize residues whose sidechain terminal atom is within
    ``radius_A`` of ``center``."""
    nh = ResidueNeighborhood(center=center, radius_A=radius_A)
    hydropathies: list[float] = []
    for (chain, resseq, resname), atoms_list in _group_by_residue(atoms).items():
        d = residue_to_terminal_distance(center, atoms_list)
        if d is None or d > radius_A:
            continue
        nh.residues.append((chain, resseq, resname, d))
        rn = resname.upper()
        nh.net_charge += FORMAL_CHARGE.get(rn, 0.0)
        if rn in BASIC_RESIDUES:
            nh.n_basic += 1
        if rn in ACIDIC_RESIDUES:
            nh.n_acidic += 1
        if rn in AROMATIC_RESIDUES:
            nh.n_aromatic += 1
        if rn in HBOND_DONOR_RESIDUES:
            nh.n_hbond_donors += 1
        if rn in HBOND_ACCEPTOR_RESIDUES:
            nh.n_hbond_acceptors += 1
        if rn in HYDROPATHY:
            hydropathies.append(HYDROPATHY[rn])
    if hydropathies:
        nh.mean_hydropathy = float(sum(hydropathies) / len(hydropathies))
    nh.residues.sort(key=lambda x: x[3])
    return nh


def count_basic_within(atoms: Iterable[Atom], center: tuple[float, float, float],
                       radius_A: float = 5.0) -> int:
    return residue_neighborhood(atoms, center, radius_A).n_basic


def multi_shell_summary(atoms: Iterable[Atom], center: tuple[float, float, float],
                        radii_A: tuple[float, ...] = (5.0, 8.0, 10.0)) -> dict:
    """Return basic / acidic / charge for each radius."""
    out: dict[str, dict] = {}
    for r in radii_A:
        nh = residue_neighborhood(atoms, center, r)
        out[f"r{int(r)}A"] = {
            "n_basic": nh.n_basic,
            "n_acidic": nh.n_acidic,
            "n_aromatic": nh.n_aromatic,
            "net_charge": nh.net_charge,
            "n_hbond_donors": nh.n_hbond_donors,
            "n_hbond_acceptors": nh.n_hbond_acceptors,
            "mean_hydropathy": nh.mean_hydropathy,
        }
    return out
