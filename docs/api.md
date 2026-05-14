# API quick reference

## Configuration

```python
from crypticip.config import load_config
cfg = load_config("config/validation.yaml", organism="yeast")
cfg.dotted("scoring.weights.depth")      # 0.25
cfg.hash()                                # "ab12cd34ef56"
```

## Paths

```python
from crypticip.paths import ProjectPaths
paths = ProjectPaths.from_config(cfg).ensure()
paths.organism_proteome_dir("yeast")
```

## Pipeline steps in isolation

```python
from crypticip.structures import preprocess_structure
from crypticip.fpocket    import run_fpocket
from crypticip.sasa       import compute_sasa, mean_sidechain_sasa
from crypticip.electrostatics import run_apbs
from crypticip.residues   import residue_neighborhood
from crypticip.scoring    import FeatureVector, composite_score, filter_flags, tier

meta = preprocess_structure("data/pdb/1ZY7.pdb", out_dir=paths.cache_dir,
                            name="ADAR2", is_alphafold=False, chain="A")
pockets, run_meta = run_fpocket(meta.cleaned_path)
best = pockets[0]
sasa = compute_sasa(meta.cleaned_path)
nh = residue_neighborhood(parse_pdb_atoms(meta.cleaned_path), best.center, 5.0)
elec = run_apbs(meta.cleaned_path, best.center)["potential_kT_e"]

fv = FeatureVector(depth=best.depth, sasa=..., elec=elec,
                   basic_count=nh.n_basic, volume=best.volume,
                   plddt=meta.mean_plddt)
score = composite_score(fv).composite
flags = filter_flags(fv)
the_tier = tier(score, flags, apbs_status="ok")
```

## Validation

```python
from crypticip.validation import validate_all, write_validation_outputs
report = validate_all(cfg, paths=paths)
out = write_validation_outputs(report, paths=paths)
```

## Screening

```python
from crypticip.screening import screen_proteome
summary = screen_proteome("yeast", cfg=cfg, paths=paths,
                          mode="fast", workers=8, limit=200)
```

## Reporting

```python
from crypticip.reporting import (write_validation_report,
                                 write_screening_report,
                                 write_all_organism_report)
```

## Experimental + PyMOL

```python
from crypticip.experimental import build_experimental_plan
from crypticip.pymol        import generate_pymol_bundle, write_pml
```
