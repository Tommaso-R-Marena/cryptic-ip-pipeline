#!/usr/bin/env python3
"""
Complete validation analysis: fpocket + FreeSASA + charge analysis.
Produces real, computed metrics for all 9 control structures.
"""

import os, re, json, csv, math, sys
import numpy as np
from pathlib import Path
from collections import defaultdict

import freesasa
from Bio.PDB import PDBParser, NeighborSearch
from scipy import stats

PROJ = Path(__file__).resolve().parent.parent
BASIC = {'ARG', 'LYS', 'HIS'}
CHARGE_MAP = {'ARG': +1, 'LYS': +1, 'HIS': +0.5, 'ASP': -1, 'GLU': -1}

# ─── Known IP-binding residues from crystal structures ──────────────────
CONTROLS = {
    'ADAR2_crystal': {
        'chain': 'A', 'residues': [376, 519, 522, 651, 672, 687],
        'ip': 'IP6', 'cat': 'positive', 'pdb': '1ZY7',
        'desc': 'ADAR2 deaminase (buried IP6)'
    },
    'ADAR2_alphafold': {
        'chain': 'A', 'residues': [376, 519, 522, 651, 672, 687],
        'ip': 'IP6', 'cat': 'positive', 'pdb': 'AF-P78563',
        'desc': 'ADAR2 AlphaFold model'
    },
    'HDAC1_crystal': {
        'chain': 'A', 'residues': [29, 31, 36, 270, 271],
        'ip': 'IP4', 'cat': 'positive', 'pdb': '5ICN',
        'desc': 'HDAC1 deacetylase (buried IP4)'
    },
    'HDAC3_crystal': {
        'chain': 'A', 'residues': [29, 31, 36, 265, 266],
        'ip': 'IP4', 'cat': 'positive', 'pdb': '4A69',
        'desc': 'HDAC3 deacetylase (buried IP4)'
    },
    'Pds5B_crystal': {
        'chain': 'A', 'residues': [1059, 1060, 1086, 1131, 1132],
        'ip': 'IP6', 'cat': 'positive', 'pdb': '5HDT',
        'desc': 'Pds5B cohesin regulator (buried IP6)'
    },
    'PLCd1_PH': {
        'chain': 'A', 'residues': [30, 32, 36, 40, 51, 55],
        'ip': 'IP3', 'cat': 'negative', 'pdb': '1MAI',
        'desc': 'PLCδ1 PH domain (surface IP3)'
    },
    'Btk_PH': {
        'chain': 'A', 'residues': [12, 28, 29, 33],
        'ip': 'IP4', 'cat': 'negative', 'pdb': '1BTK',
        'desc': 'Btk PH domain (surface IP4)'
    },
    'DAPP1_PH': {
        'chain': 'A', 'residues': [17, 19, 36, 73],
        'ip': 'IP4', 'cat': 'negative', 'pdb': '1FAO',
        'desc': 'DAPP1 PH domain (surface IP4)'
    },
    'Grp1_PH': {
        'chain': 'A', 'residues': [267, 271, 273, 276, 280, 303],
        'ip': 'IP4', 'cat': 'negative', 'pdb': '1FGY',
        'desc': 'Grp1 PH domain (surface IP4)'
    },
}


def compute_sasa(pdb_path):
    """Return {(chain, resnum_int, resname): sasa} using freesasa."""
    struct = freesasa.Structure(pdb_path)
    result = freesasa.calc(struct)
    residue_sasa = defaultdict(float)
    for i in range(struct.nAtoms()):
        ch = struct.chainLabel(i)
        rn = struct.residueNumber(i).strip()
        rname = struct.residueName(i).strip()
        try:
            rn_int = int(rn)
        except ValueError:
            continue
        residue_sasa[(ch, rn_int, rname)] += result.atomArea(i)
    return dict(residue_sasa), result.totalArea()


def parse_fpocket_descriptors(info_file):
    """Parse all pocket descriptors from fpocket info file."""
    pockets = {}
    current_num = None
    current = {}
    
    with open(info_file) as f:
        for line in f:
            line = line.strip()
            m = re.match(r'Pocket\s+(\d+)\s*:', line)
            if m:
                if current_num is not None:
                    pockets[current_num] = current
                current_num = int(m.group(1))
                current = {}
            elif '\t-' in line and current_num is not None:
                parts = line.split('\t-')
                if len(parts) == 2:
                    key = parts[0].strip().rstrip(':').strip()
                    val = parts[1].strip()
                    try:
                        val = float(val)
                    except ValueError:
                        pass
                    current[key] = val
    if current_num is not None:
        pockets[current_num] = current
    return pockets


def get_pocket_residues(pockets_dir, pocket_num):
    """Get residue numbers from a pocket's atom file."""
    atm_file = os.path.join(pockets_dir, f'pocket{pocket_num}_atm.pdb')
    if not os.path.exists(atm_file):
        return set(), []
    
    residues = set()
    coords = []
    with open(atm_file) as f:
        for line in f:
            if line.startswith('ATOM'):
                ch = line[21]
                try:
                    resseq = int(line[22:26])
                except ValueError:
                    continue
                resname = line[17:20].strip()
                residues.add(resseq)
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
    
    center = np.mean(coords, axis=0) if coords else np.zeros(3)
    return residues, center


def find_best_pocket(fpocket_dir, known_residues):
    """Find pocket with best overlap to known IP-binding residues."""
    pockets_dir = os.path.join(fpocket_dir, 'pockets')
    if not os.path.exists(pockets_dir):
        return 0, 0, set(), np.zeros(3)
    
    known_set = set(known_residues)
    best_num = 0
    best_overlap = 0
    best_residues = set()
    best_center = np.zeros(3)
    
    for f in os.listdir(pockets_dir):
        if not f.endswith('_atm.pdb'):
            continue
        m = re.search(r'pocket(\d+)_atm', f)
        if not m:
            continue
        pnum = int(m.group(1))
        res, center = get_pocket_residues(pockets_dir, pnum)
        overlap = len(known_set & res)
        if overlap > best_overlap or (overlap == best_overlap and pnum < best_num):
            best_overlap = overlap
            best_num = pnum
            best_residues = res
            best_center = center
    
    return best_num, best_overlap, best_residues, best_center


def analyze_all():
    results = []
    
    for name, info in CONTROLS.items():
        pdb_path = str(PROJ / 'data' / 'cleaned' / f'{name}.pdb')
        fpocket_dir = str(PROJ / 'data' / 'fpocket_results' / f'{name}_out')
        
        if not os.path.exists(pdb_path):
            print(f"SKIP: {name} — no cleaned PDB")
            continue
        
        print(f"\n{'='*70}")
        print(f"  {name}: {info['desc']}  [{info['cat'].upper()}]")
        print(f"{'='*70}")
        
        r = {'name': name, 'category': info['cat'], 'pdb': info['pdb'],
             'description': info['desc'], 'ip_type': info['ip']}
        
        # ── Parse structure ──────────────────────────────────────────
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(name, pdb_path)
        model = structure[0]
        all_atoms = list(model.get_atoms())
        r['n_atoms'] = len(all_atoms)
        r['n_residues'] = len(list(model.get_residues()))
        
        # ── SASA ─────────────────────────────────────────────────────
        sasa_dict, total_sasa = compute_sasa(pdb_path)
        r['total_sasa'] = round(total_sasa, 2)
        
        # SASA for known IP-binding residues
        ip_sasa = []
        for rnum in info['residues']:
            for (ch, rn, rname), sv in sasa_dict.items():
                if rn == rnum and ch == info['chain']:
                    ip_sasa.append({'resnum': rnum, 'resname': rname, 'sasa': round(sv, 2)})
                    break
        
        r['ip_residue_sasa'] = ip_sasa
        sasa_vals = [x['sasa'] for x in ip_sasa]
        r['ip_mean_sasa'] = round(np.mean(sasa_vals), 2) if sasa_vals else None
        r['ip_median_sasa'] = round(np.median(sasa_vals), 2) if sasa_vals else None
        r['ip_min_sasa'] = round(np.min(sasa_vals), 2) if sasa_vals else None
        r['ip_max_sasa'] = round(np.max(sasa_vals), 2) if sasa_vals else None
        print(f"  SASA (IP residues): mean={r['ip_mean_sasa']} med={r['ip_median_sasa']} "
              f"min={r['ip_min_sasa']} max={r['ip_max_sasa']} Å²")
        
        # ── fpocket ──────────────────────────────────────────────────
        # Count total pockets
        pockets_dir_path = os.path.join(fpocket_dir, 'pockets')
        if os.path.exists(pockets_dir_path):
            r['total_pockets'] = len([f for f in os.listdir(pockets_dir_path) if f.endswith('_atm.pdb')])
        else:
            r['total_pockets'] = 0
        
        # Find best-matching pocket
        best_num, overlap, pocket_res, pocket_center = find_best_pocket(
            fpocket_dir, info['residues'])
        r['best_pocket_num'] = best_num
        r['pocket_residue_overlap'] = overlap
        r['pocket_overlap_fraction'] = round(overlap / len(info['residues']), 3) if info['residues'] else 0
        r['pocket_n_lining_residues'] = len(pocket_res)
        print(f"  fpocket: {r['total_pockets']} total, best match pocket #{best_num} "
              f"(overlap {overlap}/{len(info['residues'])})")
        
        # Parse fpocket descriptors for the best pocket
        info_file = os.path.join(fpocket_dir, f'{name}_info.txt')
        fp_desc = {}
        if os.path.exists(info_file):
            all_desc = parse_fpocket_descriptors(info_file)
            if best_num in all_desc:
                fp_desc = all_desc[best_num]
        
        r['fpocket_score'] = fp_desc.get('Score', None)
        r['fpocket_druggability'] = fp_desc.get('Druggability Score', None)
        r['fpocket_volume'] = fp_desc.get('Real volume (Monte Carlo)', 
                              fp_desc.get('Volume', None))
        r['fpocket_mean_local_hyd_density'] = fp_desc.get('Mean local hydrophobic density', None)
        r['fpocket_polarity_score'] = fp_desc.get('Polarity score', None)
        r['fpocket_charge_score'] = fp_desc.get('Charge score', None)
        r['fpocket_alpha_sphere_density'] = fp_desc.get('Alpha sphere density', None)
        r['fpocket_mean_buriedness'] = fp_desc.get('Proportion of polar atoms', None)
        
        # Mean SASA of pocket-lining residues
        pocket_sasa_vals = []
        for rnum in pocket_res:
            for (ch, rn, rname), sv in sasa_dict.items():
                if rn == rnum:
                    pocket_sasa_vals.append(sv)
                    break
        r['pocket_mean_sasa'] = round(np.mean(pocket_sasa_vals), 2) if pocket_sasa_vals else None
        r['pocket_median_sasa'] = round(np.median(pocket_sasa_vals), 2) if pocket_sasa_vals else None
        
        print(f"  Pocket SASA: mean={r['pocket_mean_sasa']} med={r['pocket_median_sasa']} Å²")
        if r['fpocket_volume']:
            print(f"  Pocket volume: {r['fpocket_volume']} Å³")
        
        # ── Basic residue count ──────────────────────────────────────
        ns = NeighborSearch(all_atoms)
        nearby_5 = ns.search(pocket_center, 5.0, level='R')
        nearby_8 = ns.search(pocket_center, 8.0, level='R')
        
        basic_5 = [(res.get_resname(), res.get_id()[1]) 
                    for res in nearby_5 if res.get_resname() in BASIC]
        basic_8 = [(res.get_resname(), res.get_id()[1]) 
                    for res in nearby_8 if res.get_resname() in BASIC]
        
        r['basic_5A'] = len(basic_5)
        r['basic_8A'] = len(basic_8)
        r['basic_5A_list'] = [f"{rn}{rnum}" for rn, rnum in basic_5]
        print(f"  Basic residues: {r['basic_5A']} (5Å), {r['basic_8A']} (8Å)")
        
        # ── Charge density ───────────────────────────────────────────
        charge_8 = sum(CHARGE_MAP.get(res.get_resname(), 0) for res in nearby_8)
        r['net_charge_8A'] = charge_8
        print(f"  Net charge (8Å): {charge_8:+.1f}")
        
        # ── Pocket depth estimate ────────────────────────────────────
        # Use distance from pocket center to protein surface
        if len(all_atoms) > 0:
            all_coords = np.array([a.get_vector().get_array() for a in all_atoms])
            centroid = np.mean(all_coords, axis=0)
            radii = np.linalg.norm(all_coords - centroid, axis=1)
            max_radius = np.max(radii)
            dist_center_to_surface = max_radius - np.linalg.norm(pocket_center - centroid)
            # Normalize by protein radius to get fraction
            burial_fraction = dist_center_to_surface / max_radius if max_radius > 0 else 0
            # Convert to an absolute depth proxy (multiply by diameter)
            pocket_depth = max(0, dist_center_to_surface)
        else:
            pocket_depth = 0
            burial_fraction = 0
        
        r['pocket_depth'] = round(pocket_depth, 2)
        r['burial_fraction'] = round(burial_fraction, 3)
        print(f"  Pocket depth: {pocket_depth:.1f} Å, burial fraction: {burial_fraction:.3f}")
        
        # ── Composite score ──────────────────────────────────────────
        sasa_for_score = r['ip_mean_sasa'] if r['ip_mean_sasa'] is not None else (r['pocket_mean_sasa'] or 25)
        
        # Normalization
        n_depth = min(pocket_depth / 30.0, 1.0)
        n_sasa  = 1.0 - min(sasa_for_score / 150.0, 1.0)  # Scale to 150 for better discrimination
        n_charge = min(max(charge_8, 0) / 15.0, 1.0)
        n_basic  = min(r['basic_5A'] / 8.0, 1.0)
        
        score = 0.30 * n_depth + 0.35 * n_sasa + 0.20 * n_charge + 0.15 * n_basic
        r['composite_score'] = round(score, 4)
        r['score_depth'] = round(n_depth, 4)
        r['score_sasa'] = round(n_sasa, 4)
        r['score_charge'] = round(n_charge, 4)
        r['score_basic'] = round(n_basic, 4)
        
        print(f"  COMPOSITE SCORE: {score:.4f}")
        print(f"    depth={n_depth:.3f}  sasa={n_sasa:.3f}  charge={n_charge:.3f}  basic={n_basic:.3f}")
        
        results.append(r)
    
    # ─── Save results ────────────────────────────────────────────────
    results_dir = PROJ / 'results'
    results_dir.mkdir(exist_ok=True)
    
    def np_convert(obj):
        if isinstance(obj, (np.integer,)): return int(obj)
        if isinstance(obj, (np.floating,)): return float(obj)
        if isinstance(obj, np.ndarray): return obj.tolist()
        return obj
    
    with open(results_dir / 'validation_results.json', 'w') as f:
        json.dump(results, f, indent=2, default=np_convert)
    
    # Summary CSV
    fields = ['name', 'category', 'pdb', 'ip_type', 'n_residues', 'total_pockets',
              'best_pocket_num', 'pocket_residue_overlap', 'pocket_overlap_fraction',
              'ip_mean_sasa', 'ip_median_sasa', 'pocket_mean_sasa',
              'basic_5A', 'basic_8A', 'net_charge_8A',
              'pocket_depth', 'burial_fraction',
              'fpocket_score', 'fpocket_druggability', 'fpocket_volume',
              'composite_score', 'score_depth', 'score_sasa', 'score_charge', 'score_basic']
    
    with open(results_dir / 'validation_summary.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        for r in results:
            writer.writerow(r)
    
    # ─── Print summary table ─────────────────────────────────────────
    print(f"\n\n{'='*120}")
    print("VALIDATION RESULTS SUMMARY")
    print(f"{'='*120}")
    hdr = f"{'Name':<22} {'Cat':^5} {'PDB':^10} {'Pockets':>7} {'Rank':>5} {'IP SASA':>9} {'Pkt SASA':>9} {'B5Å':>4} {'Chg':>6} {'Depth':>7} {'Score':>7}"
    print(hdr)
    print("-" * 120)
    
    for r in sorted(results, key=lambda x: -x['composite_score']):
        cat = '[+]' if r['category'] == 'positive' else '[-]'
        print(f"{r['name']:<22} {cat:^5} {r['pdb']:^10} {r['total_pockets']:>7} "
              f"{r['best_pocket_num']:>5} "
              f"{r['ip_mean_sasa'] or 0:>8.1f}  "
              f"{r['pocket_mean_sasa'] or 0:>8.1f}  "
              f"{r['basic_5A']:>3}  "
              f"{r['net_charge_8A']:>+5.1f}  "
              f"{r['pocket_depth']:>6.1f}  "
              f"{r['composite_score']:>6.4f}")
    
    # ─── Statistics ──────────────────────────────────────────────────
    pos = [r for r in results if r['category'] == 'positive']
    neg = [r for r in results if r['category'] == 'negative']
    
    print(f"\n{'─'*80}")
    print("STATISTICAL COMPARISON")
    print(f"{'─'*80}")
    
    metrics = [
        ('composite_score', 'Composite Score'),
        ('ip_mean_sasa', 'IP Site Mean SASA (Å²)'),
        ('pocket_mean_sasa', 'Pocket Mean SASA (Å²)'),
        ('basic_5A', 'Basic Residues (5Å)'),
        ('net_charge_8A', 'Net Charge (8Å)'),
        ('pocket_depth', 'Pocket Depth (Å)'),
    ]
    
    for key, label in metrics:
        pv = [r[key] for r in pos if r[key] is not None]
        nv = [r[key] for r in neg if r[key] is not None]
        if pv and nv:
            pm, ps = np.mean(pv), np.std(pv)
            nm, ns2 = np.mean(nv), np.std(nv)
            # Effect size
            pooled = np.sqrt((np.var(pv) + np.var(nv)) / 2)
            d = (pm - nm) / pooled if pooled > 0 else float('inf')
            # Mann-Whitney
            if len(pv) >= 2 and len(nv) >= 2:
                u, p = stats.mannwhitneyu(pv, nv, alternative='two-sided')
            else:
                u, p = 0, 1.0
            print(f"  {label:<30} Pos: {pm:>8.2f}±{ps:>6.2f}  Neg: {nm:>8.2f}±{ns2:>6.2f}  "
                  f"d={d:>+6.2f}  U={u:>5.0f}  p={p:.4f}")
    
    print(f"\nFiles saved:")
    print(f"  {results_dir / 'validation_results.json'}")
    print(f"  {results_dir / 'validation_summary.csv'}")


if __name__ == '__main__':
    analyze_all()
