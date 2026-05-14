"""Structure preprocessing: cleaning, chain selection, pLDDT extraction."""
from __future__ import annotations

import statistics
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from .logging_utils import get_logger
from .pdb_io import (Atom, parse_pdb_atoms, clean_for_alphafold,
                     clean_preserving_ip, chain_subset, detect_ip_ligands,
                     protein_atoms_only, write_atoms)

log = get_logger(__name__)


@dataclass
class StructureMeta:
    name: str
    source_path: Path
    cleaned_path: Path | None = None
    is_alphafold: bool = False
    chain: str | None = None
    n_atoms_raw: int = 0
    n_atoms_clean: int = 0
    n_protein_atoms: int = 0
    n_residues: int = 0
    detected_ip_ligands: list[dict] = field(default_factory=list)
    mean_plddt: float | None = None
    median_plddt: float | None = None
    fraction_plddt_high: float | None = None
    fraction_plddt_low: float | None = None


def _residue_keys(atoms: Iterable[Atom]) -> set[tuple[str, int, str]]:
    return {(a.chain, a.resseq, a.icode) for a in atoms if a.record == "ATOM"}


def per_residue_ca_plddt(atoms: Iterable[Atom]) -> list[float]:
    """Return list of pLDDT values, one per residue (using CA B-factor)."""
    out: list[float] = []
    seen: set[tuple[str, int, str]] = set()
    for a in atoms:
        if a.record != "ATOM" or a.name != "CA":
            continue
        key = (a.chain, a.resseq, a.icode)
        if key in seen:
            continue
        seen.add(key)
        out.append(float(a.bfactor))
    return out


def plddt_summary(values: list[float]) -> dict:
    if not values:
        return {"mean": None, "median": None, "fraction_high": None, "fraction_low": None, "n": 0}
    n = len(values)
    high = sum(1 for v in values if v >= 70.0) / n
    low = sum(1 for v in values if v < 50.0) / n
    return {
        "mean": float(sum(values) / n),
        "median": float(statistics.median(values)),
        "fraction_high": float(high),
        "fraction_low": float(low),
        "n": n,
    }


def preprocess_structure(source: Path | str, *,
                         out_dir: Path,
                         name: str | None = None,
                         is_alphafold: bool = False,
                         chain: str | None = None) -> StructureMeta:
    """Clean a structure, dump a normalized PDB, and return metadata.

    For AlphaFold structures, all HETATM are dropped. For crystals,
    IP-family ligands are preserved (so we can locate the IP site) but
    waters/buffers/ions are dropped.
    """
    source = Path(source)
    name = name or source.stem
    atoms = parse_pdb_atoms(source)
    n_raw = len(atoms)

    if is_alphafold:
        cleaned = clean_for_alphafold(atoms)
    else:
        cleaned = clean_preserving_ip(atoms)

    if chain:
        cleaned = [a for a in cleaned if a.chain == chain]

    protein = protein_atoms_only(cleaned)
    plddt = per_residue_ca_plddt(protein) if is_alphafold else []
    plddt_stats = plddt_summary(plddt)

    out_dir.mkdir(parents=True, exist_ok=True)
    cleaned_path = out_dir / f"{name}_clean.pdb"
    write_atoms(cleaned_path, cleaned)

    ligands = detect_ip_ligands(atoms)

    residues = len(_residue_keys(protein))

    return StructureMeta(
        name=name,
        source_path=source,
        cleaned_path=cleaned_path,
        is_alphafold=is_alphafold,
        chain=chain,
        n_atoms_raw=n_raw,
        n_atoms_clean=len(cleaned),
        n_protein_atoms=len(protein),
        n_residues=residues,
        detected_ip_ligands=ligands,
        mean_plddt=plddt_stats["mean"],
        median_plddt=plddt_stats["median"],
        fraction_plddt_high=plddt_stats["fraction_high"],
        fraction_plddt_low=plddt_stats["fraction_low"],
    )
