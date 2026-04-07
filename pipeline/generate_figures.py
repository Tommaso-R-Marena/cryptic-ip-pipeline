#!/usr/bin/env python3
"""
Generate publication-quality figures from real validation results.
Reads from results/validation_results.json.
"""

import json
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
from pathlib import Path

PROJ = Path(__file__).resolve().parent.parent
FIG_DIR = PROJ / 'figures'
FIG_DIR.mkdir(exist_ok=True)

# ─── Load results ─────────────────────────────────────────────────────────
with open(PROJ / 'results' / 'validation_results.json') as f:
    data = json.load(f)

pos = [r for r in data if r['category'] == 'positive']
neg = [r for r in data if r['category'] == 'negative']

# ─── Color palette ─────────────────────────────────────────────────────────
TEAL      = '#20808D'
DARK_TEAL = '#1B474D'
BURGUNDY  = '#7B1042'
ORANGE    = '#A84B2F'
GOLD      = '#C8A43A'
NEUTRAL   = '#2E2E2E'
MUTED     = '#7A7974'
BG        = '#F7F6F2'
LIGHT_CYAN = '#DCF0F2'

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.linewidth': 0.8,
    'figure.facecolor': BG,
    'axes.facecolor': BG,
})


def get_vals(records, key):
    """Extract values safely, replacing None with 0."""
    return [r.get(key) or 0 for r in records]


# ════════════════════════════════════════════════════════════════════════════
# FIGURE 1: Composite Score Bar Chart — pos vs neg
# ════════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(10, 5))

names = [r['name'].replace('_', '\n') for r in sorted(data, key=lambda x: -x['composite_score'])]
scores = [r['composite_score'] for r in sorted(data, key=lambda x: -x['composite_score'])]
colors = [TEAL if r['category'] == 'positive' else ORANGE 
          for r in sorted(data, key=lambda x: -x['composite_score'])]

bars = ax.bar(range(len(names)), scores, color=colors, width=0.65, edgecolor='white', linewidth=0.5)
ax.set_xticks(range(len(names)))
ax.set_xticklabels(names, fontsize=8.5, rotation=0, ha='center')
ax.set_ylabel('Composite Burial Score', fontsize=11)
ax.set_title('Composite Score: Positive Controls vs. Negative Controls', 
             fontsize=12, fontweight='bold', color=NEUTRAL, loc='left')

# Score labels on bars
for i, (s, bar) in enumerate(zip(scores, bars)):
    ax.text(bar.get_x() + bar.get_width()/2, s + 0.015, f'{s:.3f}', 
            ha='center', va='bottom', fontsize=8.5, fontweight='bold', color=NEUTRAL)

# Threshold lines
pos_mean = np.mean([r['composite_score'] for r in pos])
neg_mean = np.mean([r['composite_score'] for r in neg])
ax.axhline(pos_mean, color=TEAL, linestyle='--', lw=1.2, alpha=0.7,
           label=f'Positive mean: {pos_mean:.3f}')
ax.axhline(neg_mean, color=ORANGE, linestyle='--', lw=1.2, alpha=0.7,
           label=f'Negative mean: {neg_mean:.3f}')

# Legend
pos_patch = mpatches.Patch(color=TEAL, label=f'Positive controls (n={len(pos)})')
neg_patch = mpatches.Patch(color=ORANGE, label=f'Negative controls (n={len(neg)})')
ax.legend(handles=[pos_patch, neg_patch], fontsize=9, frameon=False, loc='upper right')

ax.set_ylim(0, max(scores) * 1.15)
plt.tight_layout()
plt.savefig(FIG_DIR / 'fig1_composite_scores.png', dpi=250, bbox_inches='tight', facecolor=BG)
plt.close()
print("Fig 1 saved: Composite scores")


# ════════════════════════════════════════════════════════════════════════════
# FIGURE 2: Multi-panel comparison (4 metrics)
# ════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(2, 2, figsize=(12, 9))

sorted_data = sorted(data, key=lambda x: -x['composite_score'])
labels = [r['name'].replace('_crystal', '').replace('_alphafold', '\n(AF)').replace('_PH', '\nPH') 
          for r in sorted_data]
colors_all = [TEAL if r['category'] == 'positive' else ORANGE for r in sorted_data]

# Panel A: IP Site Mean SASA
ax = axes[0, 0]
vals = [r.get('ip_mean_sasa') or 0 for r in sorted_data]
ax.bar(range(len(labels)), vals, color=colors_all, width=0.65)
for i, v in enumerate(vals):
    if v > 0:
        ax.text(i, v + 1, f'{v:.1f}', ha='center', va='bottom', fontsize=7.5, color=NEUTRAL)
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, fontsize=7.5)
ax.set_ylabel('Mean SASA (Å²)', fontsize=10)
ax.set_title('A   IP-Site Residue SASA', fontsize=10.5, fontweight='bold', loc='left', color=NEUTRAL)
ax.axhline(50, color=DARK_TEAL, ls='--', lw=1, label='50 Å² threshold')
ax.legend(fontsize=8, frameon=False)

# Panel B: Pocket Depth
ax = axes[0, 1]
vals = [r.get('pocket_depth') or 0 for r in sorted_data]
ax.bar(range(len(labels)), vals, color=colors_all, width=0.65)
for i, v in enumerate(vals):
    if v > 0:
        ax.text(i, v + 0.5, f'{v:.1f}', ha='center', va='bottom', fontsize=7.5, color=NEUTRAL)
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, fontsize=7.5)
ax.set_ylabel('Pocket Depth (Å)', fontsize=10)
ax.set_title('B   Pocket Depth (Distance to Surface)', fontsize=10.5, fontweight='bold', loc='left', color=NEUTRAL)
ax.axhline(15, color=DARK_TEAL, ls='--', lw=1, label='15 Å buried threshold')
ax.legend(fontsize=8, frameon=False)

# Panel C: Net Charge (8Å)
ax = axes[1, 0]
vals = [r.get('net_charge_8A') or 0 for r in sorted_data]
ax.bar(range(len(labels)), vals, color=colors_all, width=0.65)
for i, v in enumerate(vals):
    ax.text(i, v + (0.15 if v >= 0 else -0.3), f'{v:+.1f}', 
            ha='center', va='bottom', fontsize=7.5, color=NEUTRAL)
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, fontsize=7.5)
ax.set_ylabel('Net Formal Charge', fontsize=10)
ax.set_title('C   Local Charge Density (8 Å radius)', fontsize=10.5, fontweight='bold', loc='left', color=NEUTRAL)
ax.axhline(0, color='#BBB', ls='-', lw=0.8)
ax.axhline(4, color=DARK_TEAL, ls='--', lw=1, label='Positive threshold (+4)')
ax.legend(fontsize=8, frameon=False)

# Panel D: Basic Residues (5Å and 8Å)
ax = axes[1, 1]
vals_5 = [r.get('basic_5A') or 0 for r in sorted_data]
vals_8 = [r.get('basic_8A') or 0 for r in sorted_data]
x = np.arange(len(labels))
w = 0.32
bars5 = ax.bar(x - w/2, vals_5, w, color=[TEAL if c == TEAL else ORANGE for c in colors_all], 
               label='Within 5 Å', alpha=0.9)
bars8 = ax.bar(x + w/2, vals_8, w, color=[DARK_TEAL if c == TEAL else '#D4825A' for c in colors_all],
               label='Within 8 Å', alpha=0.7)
for i, (v5, v8) in enumerate(zip(vals_5, vals_8)):
    ax.text(i - w/2, v5 + 0.1, str(v5), ha='center', va='bottom', fontsize=7.5, color=NEUTRAL)
    ax.text(i + w/2, v8 + 0.1, str(v8), ha='center', va='bottom', fontsize=7.5, color=NEUTRAL)
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, fontsize=7.5)
ax.set_ylabel('Count', fontsize=10)
ax.set_title('D   Basic Residues Near Pocket', fontsize=10.5, fontweight='bold', loc='left', color=NEUTRAL)
ax.legend(fontsize=8, frameon=False)

fig.suptitle('Figure 2. Real computed validation metrics across 9 control structures\n'
             '(fpocket + FreeSASA + charge analysis on PDB crystal structures and AlphaFold model)',
             fontsize=9, color=MUTED, y=0.01, style='italic')

plt.tight_layout(pad=1.2)
plt.savefig(FIG_DIR / 'fig2_multi_panel.png', dpi=250, bbox_inches='tight', facecolor=BG)
plt.close()
print("Fig 2 saved: Multi-panel metrics")


# ════════════════════════════════════════════════════════════════════════════
# FIGURE 3: Score Component Breakdown (stacked bar)
# ════════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(10, 5.5))

sorted_data2 = sorted(data, key=lambda x: -x['composite_score'])
names2 = [r['name'].replace('_crystal', '').replace('_alphafold', ' (AF)').replace('_PH', ' PH') 
          for r in sorted_data2]

depth_vals  = [r['score_depth'] * 0.30 for r in sorted_data2]
sasa_vals   = [r['score_sasa'] * 0.35 for r in sorted_data2]
charge_vals = [r['score_charge'] * 0.20 for r in sorted_data2]
basic_vals  = [r['score_basic'] * 0.15 for r in sorted_data2]

x = np.arange(len(names2))
w = 0.55

ax.bar(x, depth_vals, w, label='Pocket Depth (30%)', color=DARK_TEAL)
ax.bar(x, sasa_vals, w, bottom=depth_vals, label='Inverse SASA (35%)', color=TEAL)
ax.bar(x, charge_vals, w, bottom=np.array(depth_vals)+np.array(sasa_vals), 
       label='Electrostatics (20%)', color=ORANGE)
ax.bar(x, basic_vals, w, 
       bottom=np.array(depth_vals)+np.array(sasa_vals)+np.array(charge_vals),
       label='Basic Residues (15%)', color=GOLD)

# Total score labels
for i, r in enumerate(sorted_data2):
    total = r['composite_score']
    ax.text(i, total + 0.01, f'{total:.3f}', ha='center', va='bottom', 
            fontsize=8.5, fontweight='bold', color=NEUTRAL)

# Category markers
for i, r in enumerate(sorted_data2):
    marker = '+' if r['category'] == 'positive' else '−'
    color = TEAL if r['category'] == 'positive' else ORANGE
    ax.text(i, -0.04, marker, ha='center', va='top', fontsize=14, 
            fontweight='bold', color=color)

ax.set_xticks(x)
ax.set_xticklabels(names2, fontsize=8, rotation=20, ha='right')
ax.set_ylabel('Weighted Score Contribution', fontsize=10)
ax.set_title('Score Component Breakdown: Weighted Contributions to Composite Score',
             fontsize=11, fontweight='bold', loc='left', color=NEUTRAL)
ax.legend(fontsize=9, frameon=False, loc='upper right')
ax.set_ylim(-0.08, max([r['composite_score'] for r in data]) * 1.18)

plt.tight_layout()
plt.savefig(FIG_DIR / 'fig3_score_breakdown.png', dpi=250, bbox_inches='tight', facecolor=BG)
plt.close()
print("Fig 3 saved: Score breakdown")


# ════════════════════════════════════════════════════════════════════════════
# FIGURE 4: Per-residue SASA for ADAR2 IP6-binding residues
# ════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

# ADAR2 crystal
adar2_crystal = [r for r in data if r['name'] == 'ADAR2_crystal'][0]
adar2_af = [r for r in data if r['name'] == 'ADAR2_alphafold'][0]

for ax, record, title_suffix in zip(axes, [adar2_crystal, adar2_af], 
                                     ['PDB 1ZY7 (Crystal)', 'AF-P78563 (AlphaFold)']):
    residues = record.get('ip_residue_sasa', [])
    if residues:
        labels_r = [f"{r['resname']}{r['resnum']}" for r in residues]
        sasa_r = [r['sasa'] for r in residues]
        colors_r = [BURGUNDY if s > 80 else TEAL if s < 30 else DARK_TEAL for s in sasa_r]
        
        bars = ax.bar(range(len(labels_r)), sasa_r, color=colors_r, width=0.55)
        for i, (s, b) in enumerate(zip(sasa_r, bars)):
            ax.text(b.get_x() + b.get_width()/2, s + 1.5, f'{s:.1f}', 
                    ha='center', va='bottom', fontsize=8, color=NEUTRAL)
        
        ax.set_xticks(range(len(labels_r)))
        ax.set_xticklabels(labels_r, fontsize=9, rotation=30, ha='right')
        ax.set_ylabel('SASA (Å²)', fontsize=10)
        ax.set_title(f'ADAR2 IP6-Binding Residues — {title_suffix}',
                     fontsize=10, fontweight='bold', color=NEUTRAL)
        ax.axhline(50, color=ORANGE, ls='--', lw=1, label='Surface threshold (50 Å²)')
        ax.legend(fontsize=8, frameon=False)

fig.text(0.5, -0.02,
         'Figure 4. Per-residue SASA for the six ADAR2 IP6-coordinating residues.\n'
         'Crystal structure (left) vs. AlphaFold prediction (right). '
         'SASA computed by FreeSASA (Lee-Richards, probe 1.4 Å).',
         ha='center', fontsize=8.5, color=MUTED, style='italic')

plt.tight_layout(pad=1.0)
plt.savefig(FIG_DIR / 'fig4_adar2_sasa.png', dpi=250, bbox_inches='tight', facecolor=BG)
plt.close()
print("Fig 4 saved: ADAR2 per-residue SASA")


# ════════════════════════════════════════════════════════════════════════════
# FIGURE 5: Scatter plot — Depth vs SASA colored by category
# ════════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(8, 6))

for r in data:
    color = TEAL if r['category'] == 'positive' else ORANGE
    marker = 'o' if r['category'] == 'positive' else 's'
    sasa = r.get('ip_mean_sasa') or 0
    depth = r.get('pocket_depth') or 0
    
    ax.scatter(sasa, depth, c=color, s=120, marker=marker, edgecolors='white', 
               linewidths=1, zorder=3)
    
    # Label
    label = r['name'].replace('_crystal', '').replace('_alphafold', '\n(AF)').replace('_PH', '')
    offset = (5, 5)
    ax.annotate(label, (sasa, depth), textcoords='offset points', xytext=offset,
                fontsize=7.5, color=NEUTRAL)

ax.set_xlabel('IP-Site Mean SASA (Å²)', fontsize=11)
ax.set_ylabel('Pocket Depth (Å)', fontsize=11)
ax.set_title('Pocket Depth vs. SASA: Buried Sites are Deep and Low-SASA',
             fontsize=11, fontweight='bold', color=NEUTRAL, loc='left')

# Add quadrant lines
ax.axvline(50, color='#CCC', ls='--', lw=0.8)
ax.axhline(15, color='#CCC', ls='--', lw=0.8)

# Quadrant labels
ax.text(10, 65, 'BURIED\n(low SASA, deep)', fontsize=9, color=TEAL, 
        fontweight='bold', ha='center', alpha=0.7)
ax.text(100, 5, 'SURFACE\n(high SASA, shallow)', fontsize=9, color=ORANGE, 
        fontweight='bold', ha='center', alpha=0.7)

pos_patch = mpatches.Patch(color=TEAL, label='Positive controls (buried IP)')
neg_patch = mpatches.Patch(color=ORANGE, label='Negative controls (surface PH)')
ax.legend(handles=[pos_patch, neg_patch], fontsize=9, frameon=False, loc='upper right')

plt.tight_layout()
plt.savefig(FIG_DIR / 'fig5_scatter_depth_sasa.png', dpi=250, bbox_inches='tight', facecolor=BG)
plt.close()
print("Fig 5 saved: Depth vs SASA scatter")


# ════════════════════════════════════════════════════════════════════════════
# FIGURE 6: Pipeline schematic
# ════════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(14, 4))
ax.set_xlim(0, 14)
ax.set_ylim(0, 4)
ax.axis('off')

steps = [
    ("AlphaFold /\nPDB Structures", "PDB + AF DB\n(v4/v6 models)", DARK_TEAL),
    ("Structure\nPreparation", "Strip HETATM\nClean chains\nCheck pLDDT", TEAL),
    ("Pocket\nDetection", "fpocket 3.1\n(α-sphere\nclustering)", TEAL),
    ("Solvent\nAccessibility", "FreeSASA\n(Lee-Richards\nprobe 1.4 Å)", TEAL),
    ("Charge\nDensity", "Per-residue\nformal charge\n(R/K/H ≥ 5Å)", TEAL),
    ("Composite\nScoring", "Depth 30%\nSASA 35%\nCharge 20%\nBasic 15%", ORANGE),
    ("Candidate\nRanking", "Top hits\npLDDT ≥ 70\nManual QC", BURGUNDY),
]

box_w, box_h = 1.55, 2.6
gap = 0.38
start_x = 0.25

for i, (title, detail, color) in enumerate(steps):
    x0 = start_x + i * (box_w + gap)
    fancy = FancyBboxPatch((x0, 0.5), box_w, box_h,
                            boxstyle='round,pad=0.08', facecolor=color,
                            edgecolor='none', alpha=0.92)
    ax.add_patch(fancy)
    ax.text(x0 + box_w/2, 0.5 + box_h * 0.74,
            title, ha='center', va='center', fontsize=9, fontweight='bold',
            color='white', multialignment='center')
    ax.text(x0 + box_w/2, 0.5 + box_h * 0.32,
            detail, ha='center', va='center', fontsize=7.5,
            color='#E8E8E8', multialignment='center')
    ax.text(x0 + box_w/2, 0.5 + box_h + 0.14,
            f'Step {i+1}', ha='center', va='bottom', fontsize=7.5,
            color=color, fontweight='bold')
    if i < len(steps) - 1:
        ax_x = x0 + box_w + 0.06
        ax.annotate('', xy=(ax_x + gap - 0.06, 0.5 + box_h/2),
                    xytext=(ax_x, 0.5 + box_h/2),
                    arrowprops=dict(arrowstyle='->', color=NEUTRAL, lw=1.5))

ax.text(7, 0.15, 'Validated on: ADAR2 (1ZY7), HDAC1 (5ICN), HDAC3 (4A69), Pds5B (5HDT) vs. PLCδ1, Btk, DAPP1, Grp1 PH domains',
        ha='center', va='bottom', fontsize=8, color=MUTED, style='italic')

plt.tight_layout(pad=0.2)
plt.savefig(FIG_DIR / 'fig6_pipeline.png', dpi=250, bbox_inches='tight', facecolor=BG)
plt.close()
print("Fig 6 saved: Pipeline schematic")


# ════════════════════════════════════════════════════════════════════════════
# FIGURE 7: fpocket statistics comparison
# ════════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(8, 4.5))

sorted_all = sorted(data, key=lambda x: -x['composite_score'])
names_short = [r['name'].replace('_crystal', '').replace('_alphafold', '\n(AF)').replace('_PH', '\nPH')
               for r in sorted_all]
pocket_counts = [r['total_pockets'] for r in sorted_all]
colors_c = [TEAL if r['category'] == 'positive' else ORANGE for r in sorted_all]

ax.bar(range(len(names_short)), pocket_counts, color=colors_c, width=0.55)
for i, (pc, r) in enumerate(zip(pocket_counts, sorted_all)):
    ax.text(i, pc + 2, f'{pc}', ha='center', va='bottom', fontsize=8, color=NEUTRAL)
    # Annotate best pocket rank
    ax.text(i, -8, f'#{r["best_pocket_num"]}', ha='center', va='top', 
            fontsize=7.5, color=BURGUNDY, fontweight='bold')

ax.set_xticks(range(len(names_short)))
ax.set_xticklabels(names_short, fontsize=8)
ax.set_ylabel('Total Pockets Detected (fpocket)', fontsize=10)
ax.set_title('fpocket Pocket Counts and Best-Match Pocket Rank',
             fontsize=11, fontweight='bold', loc='left', color=NEUTRAL)
ax.text(0.02, 0.02, 'Numbers below bars: rank of pocket matching known IP site',
        transform=ax.transAxes, fontsize=8, color=MUTED, style='italic')
ax.set_ylim(-15, max(pocket_counts) * 1.15)

plt.tight_layout()
plt.savefig(FIG_DIR / 'fig7_fpocket_stats.png', dpi=250, bbox_inches='tight', facecolor=BG)
plt.close()
print("Fig 7 saved: fpocket statistics")

print(f"\nAll figures saved to {FIG_DIR}/")
print("Files:", sorted(os.listdir(FIG_DIR)))
