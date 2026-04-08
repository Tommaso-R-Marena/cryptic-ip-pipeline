#!/usr/bin/env python3
"""
Advanced analyses for cryptic IP binding site pipeline.
Adds: bootstrap/permutation stats, ROC/AUC, B-factor analysis,
H-bond profiling, hydrophobicity scoring, and conservation analysis.
"""
import json, os, re, math
import numpy as np
from pathlib import Path
from collections import defaultdict

PROJ = Path(__file__).resolve().parent.parent
DATA = PROJ / 'data'
RESULTS = PROJ / 'results'
FIGURES = PROJ / 'figures'

# ─── Load existing results ───────────────────────────────────────────────────
def load_results():
    with open(RESULTS / 'expanded_validation_results.json') as f:
        return json.load(f)


# ═══════════════════════════════════════════════════════════════════════════════
# 1. BOOTSTRAP CONFIDENCE INTERVALS + PERMUTATION TEST
# ═══════════════════════════════════════════════════════════════════════════════

def bootstrap_separation(pos_values, neg_values, n_boot=10000, ci=0.95):
    """Bootstrap CI for the difference in means between positive and negative groups."""
    rng = np.random.default_rng(42)
    pos = np.array(pos_values, dtype=float)
    neg = np.array(neg_values, dtype=float)
    observed_diff = np.mean(pos) - np.mean(neg)

    boot_diffs = []
    for _ in range(n_boot):
        boot_pos = rng.choice(pos, size=len(pos), replace=True)
        boot_neg = rng.choice(neg, size=len(neg), replace=True)
        boot_diffs.append(np.mean(boot_pos) - np.mean(boot_neg))

    boot_diffs = np.array(boot_diffs)
    alpha = 1 - ci
    ci_low = np.percentile(boot_diffs, 100 * alpha / 2)
    ci_high = np.percentile(boot_diffs, 100 * (1 - alpha / 2))
    se = np.std(boot_diffs, ddof=1)

    return {
        'observed_diff': round(float(observed_diff), 4),
        'bootstrap_mean': round(float(np.mean(boot_diffs)), 4),
        'bootstrap_se': round(float(se), 4),
        'ci_lower': round(float(ci_low), 4),
        'ci_upper': round(float(ci_high), 4),
        'ci_level': ci,
        'n_bootstrap': n_boot,
    }


def permutation_test(pos_values, neg_values, n_perm=10000):
    """Permutation test for group separation significance."""
    rng = np.random.default_rng(42)
    pos = np.array(pos_values, dtype=float)
    neg = np.array(neg_values, dtype=float)
    all_vals = np.concatenate([pos, neg])
    n_pos = len(pos)
    observed_diff = np.mean(pos) - np.mean(neg)

    count_extreme = 0
    for _ in range(n_perm):
        perm = rng.permutation(all_vals)
        perm_diff = np.mean(perm[:n_pos]) - np.mean(perm[n_pos:])
        if abs(perm_diff) >= abs(observed_diff):
            count_extreme += 1

    p_value = (count_extreme + 1) / (n_perm + 1)  # +1 for continuity correction

    return {
        'observed_diff': round(float(observed_diff), 4),
        'p_value': round(float(p_value), 6),
        'n_permutations': n_perm,
        'count_extreme': count_extreme,
        'significant_0.05': p_value < 0.05,
        'significant_0.01': p_value < 0.01,
    }


def weight_sensitivity_analysis(results, n_trials=5000):
    """Assess how robust the separation is across random scoring weight perturbations."""
    rng = np.random.default_rng(42)
    base_weights = {'depth': 0.25, 'sasa': 0.30, 'charge': 0.15,
                    'basic': 0.10, 'volume': 0.10, 'apbs': 0.10}

    pos_results = [r for r in results if r['category'] == 'positive']
    neg_results = [r for r in results if r['category'] == 'negative']

    # Extract individual score components
    def get_components(r):
        return {
            'depth': r.get('score_depth', 0),
            'sasa': r.get('score_sasa', 0),
            'charge': r.get('score_charge', 0),
            'basic': r.get('score_basic', 0),
            'volume': r.get('score_volume', 0),
            'apbs': r.get('score_apbs', 0),
        }

    pos_components = [get_components(r) for r in pos_results]
    neg_components = [get_components(r) for r in neg_results]

    separation_count = 0
    separations = []

    for _ in range(n_trials):
        # Random weights (Dirichlet distribution ensures they sum to 1)
        w = rng.dirichlet(np.ones(6))
        keys = list(base_weights.keys())
        weights = {k: float(w[i]) for i, k in enumerate(keys)}

        # Recompute composite scores
        pos_scores = []
        for comp in pos_components:
            s = sum(weights[k] * comp[k] for k in keys)
            pos_scores.append(s)

        neg_scores = []
        for comp in neg_components:
            s = sum(weights[k] * comp[k] for k in keys)
            neg_scores.append(s)

        sep = min(pos_scores) - max(neg_scores)
        separations.append(sep)
        if sep > 0:
            separation_count += 1

    return {
        'n_trials': n_trials,
        'fraction_separated': round(separation_count / n_trials, 4),
        'mean_separation': round(float(np.mean(separations)), 4),
        'median_separation': round(float(np.median(separations)), 4),
        'separation_95ci': [round(float(np.percentile(separations, 2.5)), 4),
                            round(float(np.percentile(separations, 97.5)), 4)],
    }


# ═══════════════════════════════════════════════════════════════════════════════
# 2. ROC / AUC ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════════

def compute_roc_auc(results, metric_key, higher_is_positive=True):
    """Compute ROC curve and AUC for a given metric as binary classifier."""
    pairs = []
    for r in results:
        val = r.get(metric_key)
        if val is None:
            continue
        label = 1 if r['category'] == 'positive' else 0
        pairs.append((float(val), label))

    if not pairs:
        return None

    if not higher_is_positive:
        pairs = [(-v, l) for v, l in pairs]

    pairs.sort(key=lambda x: -x[0])  # descending by score

    n_pos = sum(1 for _, l in pairs if l == 1)
    n_neg = sum(1 for _, l in pairs if l == 0)

    if n_pos == 0 or n_neg == 0:
        return None

    tpr_list = [0.0]
    fpr_list = [0.0]
    tp = 0
    fp = 0

    for val, label in pairs:
        if label == 1:
            tp += 1
        else:
            fp += 1
        tpr_list.append(tp / n_pos)
        fpr_list.append(fp / n_neg)

    # AUC by trapezoidal rule
    auc = 0.0
    for i in range(1, len(tpr_list)):
        auc += (fpr_list[i] - fpr_list[i-1]) * (tpr_list[i] + tpr_list[i-1]) / 2

    return {
        'metric': metric_key,
        'auc': round(float(auc), 4),
        'n_positive': n_pos,
        'n_negative': n_neg,
        'tpr': [round(t, 4) for t in tpr_list],
        'fpr': [round(f, 4) for f in fpr_list],
    }


def multi_metric_roc(results):
    """Compute ROC/AUC for all relevant metrics."""
    metrics = [
        ('pocket_depth', True),        # deeper = more likely positive
        ('composite_score', True),      # higher = more likely positive
        ('apbs_potential_kTe', True),   # higher positive potential
        ('pocket_volume', True),        # larger volumes
        ('basic_8A', True),             # more basic residues
        ('net_charge_8A', True),        # more positive charge
    ]

    roc_results = {}
    for metric, higher_pos in metrics:
        roc = compute_roc_auc(results, metric, higher_pos)
        if roc:
            roc_results[metric] = roc
            print(f'  ROC AUC for {metric}: {roc["auc"]:.3f}')

    return roc_results


# ═══════════════════════════════════════════════════════════════════════════════
# 3. B-FACTOR (CRYSTALLOGRAPHIC FLEXIBILITY) ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════════

def analyze_bfactors(pdb_path, pocket_residues=None, ip_residues=None):
    """Analyze B-factors at binding site vs global structure.
    Low B-factors at pocket = rigid pocket. High = flexible."""
    all_bfactors = []
    pocket_bfactors = []
    ip_bfactors = []
    first_chain = None

    with open(pdb_path) as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            chain = line[21]
            if first_chain is None:
                first_chain = chain
            if chain != first_chain:
                continue

            atomname = line[12:16].strip()
            if atomname != 'CA':
                continue

            try:
                resnum = int(line[22:26].strip())
                bfac = float(line[60:66].strip())
            except ValueError:
                continue

            all_bfactors.append(bfac)

            if pocket_residues and resnum in pocket_residues:
                pocket_bfactors.append(bfac)

            if ip_residues and resnum in ip_residues:
                ip_bfactors.append(bfac)

    result = {
        'global_mean_bfactor': round(float(np.mean(all_bfactors)), 2) if all_bfactors else None,
        'global_median_bfactor': round(float(np.median(all_bfactors)), 2) if all_bfactors else None,
        'global_std_bfactor': round(float(np.std(all_bfactors)), 2) if all_bfactors else None,
        'n_residues': len(all_bfactors),
    }

    if pocket_bfactors:
        result['pocket_mean_bfactor'] = round(float(np.mean(pocket_bfactors)), 2)
        result['pocket_median_bfactor'] = round(float(np.median(pocket_bfactors)), 2)
        result['pocket_n_residues'] = len(pocket_bfactors)
        # Z-score: how different is pocket from global?
        if len(all_bfactors) > 1:
            zscore = (np.mean(pocket_bfactors) - np.mean(all_bfactors)) / np.std(all_bfactors)
            result['pocket_bfactor_zscore'] = round(float(zscore), 3)
            result['pocket_less_flexible'] = zscore < 0  # negative z = less flexible than average

    if ip_bfactors:
        result['ip_mean_bfactor'] = round(float(np.mean(ip_bfactors)), 2)
        result['ip_residue_bfactors'] = {str(r): round(float(b), 2)
                                          for r, b in zip(sorted(ip_residues), ip_bfactors)
                                          if len(ip_bfactors) == len(ip_residues)}

    return result


# ═══════════════════════════════════════════════════════════════════════════════
# 4. HYDROGEN BOND DONOR/ACCEPTOR PROFILING
# ═══════════════════════════════════════════════════════════════════════════════

# H-bond capable atoms by residue type
HBOND_DONORS = {
    'ARG': ['NH1', 'NH2', 'NE'], 'LYS': ['NZ'],
    'HIS': ['ND1', 'NE2'], 'ASN': ['ND2'], 'GLN': ['NE2'],
    'SER': ['OG'], 'THR': ['OG1'], 'TYR': ['OH'], 'TRP': ['NE1'],
    'CYS': ['SG'],  # weak
}
HBOND_ACCEPTORS = {
    'ASP': ['OD1', 'OD2'], 'GLU': ['OE1', 'OE2'],
    'ASN': ['OD1'], 'GLN': ['OE1'],
    'SER': ['OG'], 'THR': ['OG1'], 'TYR': ['OH'],
    'HIS': ['ND1', 'NE2'],
}
# Backbone always contributes: N (donor), O (acceptor)


def analyze_hbond_potential(pdb_path, pocket_center, radius=8.0):
    """Count H-bond donors and acceptors within radius of pocket center.
    IP molecules (polyanionic) need many H-bond donors to stabilize binding."""
    first_chain = None
    donors = 0
    acceptors = 0
    backbone_donors = 0
    backbone_acceptors = 0
    donor_residues = []
    acceptor_residues = []

    with open(pdb_path) as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            chain = line[21]
            if first_chain is None:
                first_chain = chain
            if chain != first_chain:
                continue

            try:
                atomname = line[12:16].strip()
                resname = line[17:20].strip()
                resnum = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue

            dist = np.sqrt((x - pocket_center[0])**2 +
                          (y - pocket_center[1])**2 +
                          (z - pocket_center[2])**2)

            if dist > radius:
                continue

            # Check sidechain donors
            if resname in HBOND_DONORS and atomname in HBOND_DONORS[resname]:
                donors += 1
                donor_residues.append(f'{resname}{resnum}')

            # Check sidechain acceptors
            if resname in HBOND_ACCEPTORS and atomname in HBOND_ACCEPTORS[resname]:
                acceptors += 1
                acceptor_residues.append(f'{resname}{resnum}')

            # Backbone
            if atomname == 'N':
                backbone_donors += 1
            elif atomname == 'O':
                backbone_acceptors += 1

    return {
        'sidechain_donors': donors,
        'sidechain_acceptors': acceptors,
        'backbone_donors': backbone_donors,
        'backbone_acceptors': backbone_acceptors,
        'total_donors': donors + backbone_donors,
        'total_acceptors': acceptors + backbone_acceptors,
        'donor_to_acceptor_ratio': round((donors + backbone_donors) /
                                          max(1, acceptors + backbone_acceptors), 3),
        'donor_residues': list(set(donor_residues)),
        'acceptor_residues': list(set(acceptor_residues)),
        'radius': radius,
    }


# ═══════════════════════════════════════════════════════════════════════════════
# 5. POCKET HYDROPHOBICITY ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════════

# Kyte-Doolittle hydropathy index
HYDROPATHY = {
    'ILE': 4.5, 'VAL': 4.2, 'LEU': 3.8, 'PHE': 2.8, 'CYS': 2.5,
    'MET': 1.9, 'ALA': 1.8, 'GLY': -0.4, 'THR': -0.7, 'SER': -0.8,
    'TRP': -0.9, 'TYR': -1.3, 'PRO': -1.6, 'HIS': -3.2, 'GLU': -3.5,
    'GLN': -3.5, 'ASP': -3.5, 'ASN': -3.5, 'LYS': -3.9, 'ARG': -4.5,
}


def analyze_pocket_hydrophobicity(pdb_path, pocket_residues):
    """Compute hydrophobicity profile of pocket-lining residues.
    IP binding pockets should be relatively hydrophilic (negative hydropathy)
    due to the charged nature of inositol phosphates."""
    first_chain = None
    seen = set()
    hydropathy_scores = []
    residue_details = []

    with open(pdb_path) as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            chain = line[21]
            if first_chain is None:
                first_chain = chain
            if chain != first_chain:
                continue

            atomname = line[12:16].strip()
            if atomname != 'CA':
                continue

            resname = line[17:20].strip()
            resnum = int(line[22:26].strip())

            if resnum in pocket_residues and resnum not in seen:
                seen.add(resnum)
                h = HYDROPATHY.get(resname, 0)
                hydropathy_scores.append(h)
                residue_details.append({
                    'resname': resname, 'resnum': resnum, 'hydropathy': h
                })

    if not hydropathy_scores:
        return {'mean_hydropathy': None, 'n_residues': 0}

    scores = np.array(hydropathy_scores)
    n_hydrophilic = sum(1 for s in scores if s < -1.0)
    n_hydrophobic = sum(1 for s in scores if s > 1.0)

    return {
        'mean_hydropathy': round(float(np.mean(scores)), 3),
        'median_hydropathy': round(float(np.median(scores)), 3),
        'std_hydropathy': round(float(np.std(scores)), 3),
        'n_residues': len(scores),
        'n_hydrophilic': n_hydrophilic,
        'n_hydrophobic': n_hydrophobic,
        'fraction_hydrophilic': round(n_hydrophilic / len(scores), 3),
        'fraction_hydrophobic': round(n_hydrophobic / len(scores), 3),
        'is_hydrophilic_pocket': float(np.mean(scores)) < 0,
    }


# ═══════════════════════════════════════════════════════════════════════════════
# 6. SEQUENCE CONSERVATION (ENTROPY) ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════════

def compute_position_entropy(pdb_path, ip_residues):
    """Approximate position-specific conservation by computing amino acid
    composition statistics at IP-coordinating positions.
    For a true conservation analysis, we'd need a MSA — here we characterize
    the amino acid types at known IP-binding positions as a proxy."""
    residue_types = {}
    first_chain = None

    with open(pdb_path) as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            chain = line[21]
            if first_chain is None:
                first_chain = chain
            if chain != first_chain:
                continue

            atomname = line[12:16].strip()
            if atomname != 'CA':
                continue

            resname = line[17:20].strip()
            resnum = int(line[22:26].strip())

            if resnum in ip_residues:
                residue_types[resnum] = resname

    # Characterize the binding motif
    charged_types = {'ARG', 'LYS', 'HIS', 'ASP', 'GLU'}
    basic_types = {'ARG', 'LYS', 'HIS'}
    aromatic_types = {'PHE', 'TYR', 'TRP', 'HIS'}

    composition = {
        'basic': sum(1 for r in residue_types.values() if r in basic_types),
        'acidic': sum(1 for r in residue_types.values() if r in {'ASP', 'GLU'}),
        'aromatic': sum(1 for r in residue_types.values() if r in aromatic_types),
        'polar': sum(1 for r in residue_types.values() if r in {'SER', 'THR', 'ASN', 'GLN'}),
        'hydrophobic': sum(1 for r in residue_types.values() if r in {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'PRO'}),
    }
    total = len(residue_types)

    return {
        'binding_motif_residues': {str(k): v for k, v in sorted(residue_types.items())},
        'n_positions': total,
        'composition': composition,
        'fraction_basic': round(composition['basic'] / max(1, total), 3),
        'fraction_aromatic': round(composition['aromatic'] / max(1, total), 3),
        'motif_signature': f"{composition['basic']}B-{composition['aromatic']}A-{composition['polar']}P",
    }


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN: Run all advanced analyses
# ═══════════════════════════════════════════════════════════════════════════════

# IP-coordinating residues for ADAR2
ADAR2_IP_RESNUMS = {376, 519, 522, 651, 672, 687}

# Structure configurations
VALIDATION_SET = [
    {'name': 'ADAR2_crystal', 'pdb': '1ZY7.pdb', 'category': 'positive',
     'ip_residues': ADAR2_IP_RESNUMS, 'subdir': 'pdb'},
    {'name': 'ADAR2_alphafold', 'pdb': 'AF-P78563-F1.pdb', 'category': 'positive',
     'ip_residues': ADAR2_IP_RESNUMS, 'subdir': 'alphafold'},
    {'name': 'HDAC1', 'pdb': '5ICN.pdb', 'category': 'positive', 'ip_residues': set(), 'subdir': 'pdb'},
    {'name': 'HDAC3', 'pdb': '4A69.pdb', 'category': 'positive', 'ip_residues': set(), 'subdir': 'pdb'},
    {'name': 'Pds5B', 'pdb': '5HDT.pdb', 'category': 'positive', 'ip_residues': set(), 'subdir': 'pdb'},
    {'name': 'PLCd1_PH', 'pdb': '1MAI.pdb', 'category': 'negative', 'ip_residues': set(), 'subdir': 'pdb'},
    {'name': 'Btk_PH', 'pdb': '1BTK.pdb', 'category': 'negative', 'ip_residues': set(), 'subdir': 'pdb'},
    {'name': 'DAPP1_PH', 'pdb': '1FAO.pdb', 'category': 'negative', 'ip_residues': set(), 'subdir': 'pdb'},
    {'name': 'Grp1_PH', 'pdb': '1FGY.pdb', 'category': 'negative', 'ip_residues': set(), 'subdir': 'pdb'},
]


def main():
    print("=" * 70)
    print("ADVANCED ANALYSES FOR CRYPTIC IP BINDING SITE PIPELINE")
    print("=" * 70)

    data = load_results()
    results = data['validation_results']
    pos_results = [r for r in results if r['category'] == 'positive']
    neg_results = [r for r in results if r['category'] == 'negative']

    advanced = {}

    # ── 1. Bootstrap & Permutation for Pocket Depth ──
    print("\n1. BOOTSTRAP & PERMUTATION TESTS")
    pos_depths = [r.get('pocket_depth', 0) for r in pos_results]
    neg_depths = [r.get('pocket_depth', 0) for r in neg_results]
    print(f"  Positive depths: {[f'{d:.1f}' for d in pos_depths]}")
    print(f"  Negative depths: {[f'{d:.1f}' for d in neg_depths]}")

    boot_depth = bootstrap_separation(pos_depths, neg_depths)
    print(f"  Bootstrap depth diff: {boot_depth['observed_diff']:.2f} Å "
          f"(95% CI: [{boot_depth['ci_lower']:.2f}, {boot_depth['ci_upper']:.2f}])")
    advanced['bootstrap_depth'] = boot_depth

    perm_depth = permutation_test(pos_depths, neg_depths)
    print(f"  Permutation test p-value (depth): {perm_depth['p_value']:.4f} "
          f"({'significant' if perm_depth['significant_0.05'] else 'not significant'} at α=0.05)")
    advanced['permutation_depth'] = perm_depth

    # Composite score
    pos_scores = [r.get('composite_score', 0) for r in pos_results]
    neg_scores = [r.get('composite_score', 0) for r in neg_results]
    boot_score = bootstrap_separation(pos_scores, neg_scores)
    perm_score = permutation_test(pos_scores, neg_scores)
    print(f"  Bootstrap score diff: {boot_score['observed_diff']:.4f} "
          f"(95% CI: [{boot_score['ci_lower']:.4f}, {boot_score['ci_upper']:.4f}])")
    print(f"  Permutation test p-value (score): {perm_score['p_value']:.4f}")
    advanced['bootstrap_score'] = boot_score
    advanced['permutation_score'] = perm_score

    # APBS potential
    pos_apbs = [r.get('apbs_potential_kTe', 0) or 0 for r in pos_results]
    neg_apbs = [r.get('apbs_potential_kTe', 0) or 0 for r in neg_results]
    boot_apbs = bootstrap_separation(pos_apbs, neg_apbs)
    perm_apbs = permutation_test(pos_apbs, neg_apbs)
    print(f"  Bootstrap APBS diff: {boot_apbs['observed_diff']:.2f} kT/e "
          f"(95% CI: [{boot_apbs['ci_lower']:.2f}, {boot_apbs['ci_upper']:.2f}])")
    print(f"  Permutation test p-value (APBS): {perm_apbs['p_value']:.4f}")
    advanced['bootstrap_apbs'] = boot_apbs
    advanced['permutation_apbs'] = perm_apbs

    # Weight sensitivity
    print("\n  Weight sensitivity analysis...")
    weight_sens = weight_sensitivity_analysis(results)
    print(f"  Fraction of random weight combos with clean separation: "
          f"{weight_sens['fraction_separated']:.1%}")
    advanced['weight_sensitivity'] = weight_sens

    # ── 2. ROC / AUC ──
    print("\n2. ROC / AUC ANALYSIS")
    roc_results = multi_metric_roc(results)
    advanced['roc_auc'] = roc_results

    # ── 3. B-factor analysis ──
    print("\n3. B-FACTOR FLEXIBILITY ANALYSIS")
    bfactor_results = {}
    for struct in VALIDATION_SET:
        pdb_path = DATA / struct['subdir'] / struct['pdb']
        if not pdb_path.exists():
            continue

        # Get pocket residues from fpocket
        pocket_residues = set()
        fp_outdir = DATA / 'fpocket_results' / f"{struct['name']}_out"
        if fp_outdir.exists():
            pocket1_pdb = fp_outdir / 'pockets' / 'pocket1_atm.pdb'
            if pocket1_pdb.exists():
                with open(pocket1_pdb) as f:
                    for line in f:
                        if line.startswith(('ATOM', 'HETATM')):
                            try:
                                pocket_residues.add(int(line[22:26].strip()))
                            except ValueError:
                                pass

        bfac = analyze_bfactors(pdb_path, pocket_residues, struct['ip_residues'] or None)
        bfactor_results[struct['name']] = bfac

        pocket_str = f"pocket={bfac.get('pocket_mean_bfactor', 'N/A')}"
        global_str = f"global={bfac['global_mean_bfactor']}"
        z_str = f"z={bfac.get('pocket_bfactor_zscore', 'N/A')}"
        print(f"  {struct['name']:20s}: {global_str}, {pocket_str}, {z_str}")

    advanced['bfactor'] = bfactor_results

    # ── 4. H-bond profiling ──
    print("\n4. HYDROGEN BOND PROFILING")
    hbond_results = {}
    for struct in VALIDATION_SET:
        pdb_path = DATA / struct['subdir'] / struct['pdb']
        if not pdb_path.exists():
            continue

        # Use pocket center from previous results
        matching = [r for r in results if r['name'] == struct['name']]
        if not matching:
            continue
        r = matching[0]
        pc = r.get('pocket_center')
        if not pc:
            continue

        hbond = analyze_hbond_potential(pdb_path, np.array(pc))
        hbond_results[struct['name']] = hbond
        print(f"  {struct['name']:20s}: donors={hbond['total_donors']:3d}, "
              f"acceptors={hbond['total_acceptors']:3d}, "
              f"D/A ratio={hbond['donor_to_acceptor_ratio']:.2f}")

    advanced['hbond'] = hbond_results

    # ── 5. Pocket hydrophobicity ──
    print("\n5. POCKET HYDROPHOBICITY")
    hydro_results = {}
    for struct in VALIDATION_SET:
        pdb_path = DATA / struct['subdir'] / struct['pdb']
        if not pdb_path.exists():
            continue

        # Get pocket residues
        pocket_residues = set()
        fp_outdir = DATA / 'fpocket_results' / f"{struct['name']}_out"
        if fp_outdir.exists():
            pocket1_pdb = fp_outdir / 'pockets' / 'pocket1_atm.pdb'
            if pocket1_pdb.exists():
                with open(pocket1_pdb) as f:
                    for line in f:
                        if line.startswith(('ATOM', 'HETATM')):
                            try:
                                pocket_residues.add(int(line[22:26].strip()))
                            except ValueError:
                                pass

        if not pocket_residues:
            continue

        hydro = analyze_pocket_hydrophobicity(pdb_path, pocket_residues)
        hydro_results[struct['name']] = hydro
        print(f"  {struct['name']:20s}: mean hydropathy={hydro['mean_hydropathy']:+.2f}, "
              f"hydrophilic={hydro['fraction_hydrophilic']:.0%}, "
              f"hydrophobic={hydro['fraction_hydrophobic']:.0%}")

    advanced['hydrophobicity'] = hydro_results

    # ── 6. Conservation / motif analysis ──
    print("\n6. BINDING MOTIF ANALYSIS")
    # Only for ADAR2 where we know IP residues
    for struct in VALIDATION_SET:
        if not struct['ip_residues']:
            continue
        pdb_path = DATA / struct['subdir'] / struct['pdb']
        if not pdb_path.exists():
            continue
        motif = compute_position_entropy(pdb_path, struct['ip_residues'])
        print(f"  {struct['name']}: motif={motif['motif_signature']}, "
              f"residues={motif['binding_motif_residues']}")
        advanced[f'motif_{struct["name"]}'] = motif

    # ── Save results ──
    output_path = RESULTS / 'advanced_analysis_results.json'
    with open(output_path, 'w') as f:
        json.dump(advanced, f, indent=2, default=str)
    print(f"\nResults saved to {output_path}")

    # ── Summary statistics ──
    print("\n" + "=" * 70)
    print("SUMMARY OF ADVANCED ANALYSES")
    print("=" * 70)

    best_auc = max(roc_results.items(), key=lambda x: x[1]['auc'])
    print(f"Best discriminator (AUC): {best_auc[0]} = {best_auc[1]['auc']:.3f}")
    print(f"Depth separation: {boot_depth['observed_diff']:.1f} Å "
          f"(95% CI: [{boot_depth['ci_lower']:.1f}, {boot_depth['ci_upper']:.1f}], "
          f"p={perm_depth['p_value']:.4f})")
    print(f"Weight sensitivity: {weight_sens['fraction_separated']:.1%} of random weights separate groups")

    # H-bond comparison
    pos_donors = [hbond_results[n]['total_donors'] for n in hbond_results
                  if any(s['name'] == n and s['category'] == 'positive' for s in VALIDATION_SET)]
    neg_donors = [hbond_results[n]['total_donors'] for n in hbond_results
                  if any(s['name'] == n and s['category'] == 'negative' for s in VALIDATION_SET)]
    if pos_donors and neg_donors:
        print(f"H-bond donors at pocket: positives {np.mean(pos_donors):.0f} vs negatives {np.mean(neg_donors):.0f}")

    # Hydrophobicity comparison
    pos_hydro = [hydro_results[n]['mean_hydropathy'] for n in hydro_results
                 if any(s['name'] == n and s['category'] == 'positive' for s in VALIDATION_SET)]
    neg_hydro = [hydro_results[n]['mean_hydropathy'] for n in hydro_results
                 if any(s['name'] == n and s['category'] == 'negative' for s in VALIDATION_SET)]
    if pos_hydro and neg_hydro:
        print(f"Mean pocket hydropathy: positives {np.mean(pos_hydro):+.2f} vs negatives {np.mean(neg_hydro):+.2f}")

    return advanced


if __name__ == '__main__':
    main()
