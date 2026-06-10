"""
Figure 1 — SCEPTIC workflow + HCTL verdict. Publication-ready.
A: Problem, B: 7 layers, C: HCTL verdict cards, D: recommendation.
Clean, vector-quality, no literature audit, no pie chart.
"""
import os, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
matplotlib.rcParams['font.size'] = 9

OUT = '/Users/guoxutao/sc_analysis/figures_crm'; os.makedirs(OUT, exist_ok=True)

C = {
    'dk': '#1F3A5F', 'tx': '#2C3E50', 'lt': '#95A5A6',
    'green': '#27AE60', 'red': '#E74C3C', 'amber': '#F39C12',
    'bg': '#F8F9FA', 'white': '#FFFFFF', 'blue': '#2980B9',
}

fig = plt.figure(figsize=(12.5, 7.5))

# === PANEL A: Problem statement (top-left) ===
axA = fig.add_axes([0.04, 0.52, 0.28, 0.44])
axA.set_xlim(0, 10); axA.set_ylim(0, 10); axA.axis('off')
axA.text(0, 9.8, 'A  Problem', fontsize=11, weight='bold', color=C['dk'])
# Input claim box
axA.add_patch(FancyBboxPatch((0.45, 6.25), 9.1, 3.05, boxstyle='round,pad=0.2',
                              facecolor='#FDEDEC', edgecolor=C['red'], lw=2))
axA.text(5, 8.65, 'Per-sample TF-score\ncorrelation claim', fontsize=9.5, weight='bold',
         color=C['red'], ha='center', linespacing=1.2)
axA.text(5, 7.7, 'from scRNA-seq pseudobulk', fontsize=8.5, color=C['tx'], ha='center')
# Risk factors
risks = ['Global expr.\nconfounding', 'Cell-type\ncomposition', 'Batch /\ncohort effects', 'Severity /\ndisease state']
for i, r in enumerate(risks):
    x = 1.35 + i * 2.08
    axA.add_patch(FancyBboxPatch((x-0.85, 5.4), 1.8, 1.2, boxstyle='round,pad=0.08',
                                  facecolor=C['white'], edgecolor=C['lt'], lw=0.8))
    axA.text(x, 6.0, r, fontsize=6.2, ha='center', va='center', color=C['lt'], linespacing=1.1)
# Question
axA.text(5, 4.7, '↑ Are these confounders driving the correlation?',
         fontsize=8, ha='center', style='italic', color=C['dk'],
         bbox=dict(boxstyle='round', facecolor=C['bg'], edgecolor=C['lt'], pad=0.5))

# === PANEL B: SCEPTIC layers (bottom-left, below Problem) ===
axB = fig.add_axes([0.04, 0.05, 0.28, 0.42])
axB.set_xlim(0, 12); axB.set_ylim(0, 16); axB.axis('off')
axB.text(0.2, 15.7, 'B  SCEPTIC 7-layer validation', fontsize=11, weight='bold', color=C['dk'])

layers = [
    ('L1', 'Orthogonal scoring', 'Raw target-gene mean + regulon-based', C['green']),
    ('L2', 'Cross-assay (supportive)', 'Plasma proteomics / CyTOF readouts', C['amber']),
    ('L3', 'Bulk / context replication', 'Independent bulk RNA-seq or microarray', C['amber']),
    ('L4', 'Random TF-pair null', 'All pairwise TF-score |r| as empirical null', C['red']),
    ('L5', 'Composition adjustment', 'Partial r after regressing cell proportions', C['green']),
    ('L6', 'Leave-one-cohort-out', 'Meta-analysis; LOO sensitivity', C['red']),
    ('L7', 'Label permutation', 'Disease-label shuffling', C['green']),
]
ly = 14.0
for lab, name, desc, color in layers:
    ly -= 1.8
    axB.add_patch(FancyBboxPatch((0.2, ly), 11, 1.3, boxstyle='round,pad=0.08',
                                  facecolor=C['white'], edgecolor=color, lw=1.3))
    axB.add_patch(FancyBboxPatch((0.2, ly), 0.06, 1.3, boxstyle='round,pad=0.01',
                                  facecolor=color, edgecolor='none', alpha=0.7))
    axB.text(0.6, ly+0.85, f'{lab}  {name}', fontsize=8, weight='bold', color=C['dk'])
    axB.text(0.6, ly+0.25, desc, fontsize=7, color=C['tx'])

# Down arrows
for i in range(6):
    y = 14.0 - (i+1)*1.8 + 0.35
    axB.annotate('', xy=(5.7, y-0.15), xytext=(5.7, y+0.05),
                 arrowprops=dict(arrowstyle='->', color=C['lt'], lw=0.8))

# === CENTER: HCTL triad verdict (main wide panel) ===
axC = fig.add_axes([0.36, 0.05, 0.62, 0.91])
axC.set_xlim(0, 16); axC.set_ylim(0, 19); axC.axis('off')
axC.text(0.2, 18.5, 'C  HCTL triad case study', fontsize=12, weight='bold', color=C['dk'])
axC.text(0.2, 17.8, 'STAT4 / NFKB1 / GATA1  |  1,161 patients  |  6 cohorts  |  3 assays',
         fontsize=8.5, color=C['lt'])

verdicts = [
    # pair, L1-L7 status list, verdict text, color
    ('STAT4–NFKB1',
     [('L1', '+0.56 SCP548, +0.54 E-MTAB', C['green']),
      ('L2', '+0.31 Olink, +0.80 CyTOF', C['amber']),
      ('L3', 'r = +0.07 bulk (p = 0.054)', C['red']),
      ('L4', '74th %ile (not top 5%)', C['red']),
      ('L5', 'Partial r = 0.49', C['green']),
      ('L6', 'I² = 46%, cohort-dependent', C['red']),
      ('L7', 'Permutation p < 0.001', C['green'])],
     'Context-dependent\nPBMC covariation\nNot a universal axis', C['amber']),
    ('NFKB1–GATA1',
     [('L1', 'Raw mean r = +0.57 → AUCell r = +0.15', C['red']),
      ('L2', 'N/A (no GATA1 proteins in panel)', C['lt']),
      ('L3', 'Bulk r = -0.07', C['red']),
      ('L4', '13th %ile', C['red'])],
     'Rejected\nMethod-dependent\nErythroid/ambient signal', C['red']),
    ('STAT4–GATA1',
     [('L1', 'Inconsistent across methods', C['red']),
      ('L2', 'N/A (no GATA1 proteins in panel)', C['lt']),
      ('L3', 'Bulk r = -0.37 (erythroid contamination)', C['red']),
      ('L4', '2nd %ile', C['red'])],
     'Rejected\nNo reproducible signal\nin any cohort', C['red']),
]

ypos = 17.0
for pair, layer_items, verdict_text, color in verdicts:
    n_items = len(layer_items)
    card_h = 4.8 if pair == 'STAT4–NFKB1' else 2.6
    axC.add_patch(FancyBboxPatch((0.2, ypos-card_h), 15.3, card_h,
                                  boxstyle='round,pad=0.15',
                                  facecolor=C['white'], edgecolor=color, lw=2.2))

    # Pair name (left)
    axC.text(0.6, ypos-0.4, pair, fontsize=11, weight='bold', color=C['dk'])

    # Layer tick marks (middle)
    lx_base = 3.8
    cols_per_layer = 2
    for j, (layer, detail, lcol) in enumerate(layer_items):
        row = j // cols_per_layer
        col = j % cols_per_layer
        lx = lx_base + col * 4.3
        ly = ypos - 0.85 - row * 0.78
        axC.add_patch(plt.Circle((lx, ly+0.25), 0.13, color=lcol))
        axC.text(lx+0.3, ly+0.35, f'{layer}: ', fontsize=7, weight='bold', color=C['dk'], va='center')
        axC.text(lx+0.3, ly-0.05, detail, fontsize=6.5, color=C['tx'], va='center')

    # Verdict badge (right)
    badge_w = 3.45
    badge_h = card_h - 1.0
    axC.add_patch(FancyBboxPatch((11.95, ypos - card_h/2 - badge_h/2), badge_w, badge_h,
                                  boxstyle='round,pad=0.12',
                                  facecolor=color, edgecolor='none', alpha=0.9))
    axC.text(11.95 + badge_w/2, ypos - card_h/2, verdict_text,
             fontsize=8.0, weight='bold', color='white', ha='center', va='center', linespacing=1.25)

    ypos -= (card_h + 0.4)

# Bottom: recommendation
axC.add_patch(FancyBboxPatch((0.2, 1.0), 15.3, 1.2, boxstyle='round,pad=0.12',
                              facecolor='#FFF3E0', edgecolor=C['amber'], lw=1.5))
axC.text(7.85, 1.6, 'Recommendation: L1 (orthogonal scoring) + L4 (random TF-pair null) as minimum reporting checks before claiming TF co-regulation.',
         fontsize=8, ha='center', va='center', color=C['dk'], style='italic')

fig.savefig(f'{OUT}/Fig1_final.tif', dpi=300, bbox_inches='tight', format='tiff', facecolor='white', edgecolor='none')
fig.savefig(f'{OUT}/Fig1_final.png', dpi=300, bbox_inches='tight', format='png', facecolor='white', edgecolor='none')
plt.close()
print("Fig1_final saved")
