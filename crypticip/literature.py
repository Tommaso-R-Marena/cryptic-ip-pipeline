"""Light-touch literature / annotation enrichment via cached UniProt fetches."""
from __future__ import annotations

import json
import time
from dataclasses import dataclass, asdict
from pathlib import Path

import requests

from .logging_utils import get_logger

log = get_logger(__name__)


@dataclass
class UniProtAnnotation:
    accession: str
    primary_name: str | None = None
    gene_name: str | None = None
    organism: str | None = None
    length: int | None = None
    subcellular_locations: list[str] | None = None
    keywords: list[str] | None = None
    function_summary: str | None = None


def fetch_uniprot(accession: str, *, cache_dir: Path, force: bool = False,
                  timeout: float = 15.0) -> UniProtAnnotation | None:
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache = cache_dir / f"{accession}.json"
    if cache.exists() and not force:
        try:
            return UniProtAnnotation(**json.loads(cache.read_text()))
        except Exception:
            pass
    try:
        r = requests.get(f"https://rest.uniprot.org/uniprotkb/{accession}.json", timeout=timeout)
        r.raise_for_status()
    except Exception as e:
        log.debug("uniprot fetch failed %s: %s", accession, e)
        return None
    data = r.json()
    ann = UniProtAnnotation(
        accession=accession,
        primary_name=(data.get("proteinDescription") or {}).get("recommendedName", {}).get("fullName", {}).get("value"),
        gene_name=(data.get("genes") or [{}])[0].get("geneName", {}).get("value"),
        organism=(data.get("organism") or {}).get("scientificName"),
        length=(data.get("sequence") or {}).get("length"),
        subcellular_locations=_locations(data),
        keywords=[kw.get("name") for kw in data.get("keywords", []) if kw.get("name")],
        function_summary=_function(data),
    )
    cache.write_text(json.dumps(asdict(ann), indent=2))
    return ann


def _locations(data: dict) -> list[str]:
    out = []
    for c in data.get("comments", []):
        if c.get("commentType") == "SUBCELLULAR LOCATION":
            for sl in c.get("subcellularLocations", []):
                loc = (sl.get("location") or {}).get("value")
                if loc:
                    out.append(loc)
    return out


def _function(data: dict) -> str | None:
    for c in data.get("comments", []):
        if c.get("commentType") == "FUNCTION":
            for t in c.get("texts", []):
                if t.get("value"):
                    return t["value"]
    return None
