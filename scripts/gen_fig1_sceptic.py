"""
Figure 1 — SCEPTIC framework: 7-layer null-control validation pipeline.
Scientific Reports style.
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUTDIR = '/Users/guoxutao/sc_analysis/figures_crm'
os.makedirs(OUTDIR, exist_ok=True)

C = {
    'pass':    '#1ABC9C',
    'fail':    '#E74C3C',
    'partial': '#F39C12',
    'layer':   '#2C3E50',
    'bg':      '#FFFFFF',
    'text':    '#2C3E50',
    'light':   '#7F8C8D',
    'border':  '#BDC3C7',
    'box_bg':  '#F8F9FA',
    'accent':  '#2980B9',
}

fig, ax = plt.subplots(figsize=(10, 7.5))
ax.set_xlim(0, 20)
ax.set_ylim(0, 15)
ax.axis('off')

# Title
ax.text(10, 14.6, 'Figure 1', fontsize=14, weight='bold', color=C['layer'], ha='center')
ax.text(10, 14.0, 'The SCEPTIC framework: seven validation layers for single-cell TF co-expression claims',
        fontsize=10, color=C['light'], ha='center', style='italic')

# === LEFT: 7 layer boxes ===
layers = [
    ('L1', 'Orthogonal scoring', 'Raw target-gene mean\n+ regulon-based (AUCell)',
     'Cross-method consistency', C['pass']),
    ('L2', 'Cross-assay replication', 'Plasma proteomics (Olink),\nintracellular staining (CyTOF)',
     'Protein-level support', C['pass']),
    ('L3', 'Bulk-tissue replication', 'Independent bulk\nRNA-seq/microarray cohort',
     'Composition-independent', C['partial']),
    ('L4', 'Random TF triplet null', '1,000 random TF triplets\nfrom same regulon database',
     'Specificity vs. chance', C['fail']),
    ('L5', 'Composition-adjusted\ncorrelation', 'Partial correlation\n| cell-type proportions',
     'Cell-type confound test', C['pass']),
    ('L6', 'Leave-one-cohort-out', 'Random-effects meta-analysis,\nI², Cochran Q',
     'Cohort robustness', C['fail']),
    ('L7', 'Label permutation', '10,000 disease-label\npermutations per cohort',
     'Sample-level significance', C['pass']),
]

y_start = 13.0
box_h = 1.25
box_w = 5.8
gap = 0.1

for i, (label, title, desc, test_type, color) in enumerate(layers):
    y = y_start - i * (box_h + gap)
    # Layer number circle
    ax.add_patch(plt.Circle((1.3, y + box_h/2), 0.28, color=C['accent'], zorder=3))
    ax.text(1.3, y + box_h/2, label[1], fontsize=11, weight='bold', color='white', ha='center', va='center')
    # Main box
    ax.add_patch(FancyBboxPatch((1.8, y), box_w, box_h,
                                 boxstyle='round,pad=0.08',
                                 facecolor=C['box_bg'], edgecolor=C['border'], linewidth=1.0, zorder=1))
    # Title
    ax.text(2.0, y + box_h - 0.2, title, fontsize=9.5, weight='bold', color=C['layer'], va='top')
    # Description
    ax.text(2.0, y + box_h - 0.55, desc, fontsize=8, color=C['text'], va='top', linespacing=1.3)
    # Test type label
    ax.text(2.0, y + 0.1, test_type, fontsize=7.5, color=color, va='bottom', style='italic')
    # Color strip on left
    ax.add_patch(FancyBboxPatch((1.8, y), 0.06, box_h,
                                 boxstyle='round,pad=0.01',
                                 facecolor=color, edgecolor='none', alpha=0.7, zorder=2))

# Down arrows between layers
for i in range(6):
    y_top = y_start - i * (box_h + gap) - gap
    ax.annotate('', xy=(1.3, y_top + 0.05), xytext=(1.3, y_top - 0.28),
                arrowprops=dict(arrowstyle='->', color=C['border'], lw=1.2))

# === RIGHT: Verdict panel ===
right_x = 8.5
right_w = 10.8

# Verdict title
ax.add_patch(FancyBboxPatch((right_x, 12.8), right_w, 1.5,
                             boxstyle='round,pad=0.1',
                             facecolor='white', edgecolor=C['accent'], linewidth=1.5))
ax.text(right_x + 0.3, 14.0, 'SCEPTIC assessment: HCTL triad case study',
        fontsize=10, weight='bold', color=C['accent'])
ax.text(right_x + 0.3, 13.45, '3 TF pairs tested across 7 layers. 2/3 pairs rejected; 1/3 conditionally supported.',
        fontsize=8.5, color=C['text'])

verdicts = [
    ('STAT4–NFKB1', 'CONDITIONALLY\nSUPPORTED', 'L1 ✓ L2 ✓ L5 ✓', 'L3 ✗ L4 ✗ L6 ✗', C['partial']),
    ('NFKB1–GATA1', 'REJECTED', '—', 'L1 ✗ L3 ✗ L4 ✗', C['fail']),
    ('STAT4–GATA1', 'REJECTED', '—', 'L1 ✗ L2 ✗ L3 ✗', C['fail']),
]

for i, (pair, verdict, passes, fails, color) in enumerate(verdicts):
    y = 11.8 - i * 1.3
    # Verdict card
    ax.add_patch(FancyBboxPatch((right_x, y - 0.1), right_w, 1.1,
                                 boxstyle='round,pad=0.08',
                                 facecolor='white', edgecolor=color, linewidth=1.5))
    ax.text(right_x + 0.3, y + 0.75, pair, fontsize=11, weight='bold', color=C['layer'])
    # Verdict badge
    badge_bg = color
    ax.add_patch(FancyBboxPatch((right_x + 3.5, y + 0.35), 3.2, 0.45,
                                 boxstyle='round,pad=0.05',
                                 facecolor=badge_bg, edgecolor='none', alpha=0.85))
    ax.text(right_x + 5.1, y + 0.57, verdict, fontsize=8, weight='bold', color='white', ha='center', va='center')
    # Passed / Failed
    ax.text(right_x + 7.0, y + 0.75, passes, fontsize=8, color=C['pass'])
    ax.text(right_x + 7.0, y + 0.35, fails, fontsize=8, color=C['fail'])

# === BOTTOM: key insight ===
ax.text(10, 7.5,
        'Core function: transforms an initially promising 3-pair finding (all p < 0.001 under raw target-gene means) '
        'into a nuanced evidence table\nwhere 2/3 pairs are rejected and the surviving pair is downgraded to '
        '"context-dependent PBMC covariation."',
        fontsize=8.5, color=C['layer'], ha='center', style='italic',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='#FFF3E0', edgecolor=C['partial'], linewidth=0.8))

# === Legend ===
legend_y = 6.8
ax.text(0.5, legend_y, '● Passed', fontsize=8, color=C['pass'])
ax.text(3.0, legend_y, '● Failed', fontsize=8, color=C['fail'])
ax.text(5.5, legend_y, '● Partial / context-dependent', fontsize=8, color=C['partial'])

plt.tight_layout()
fig.savefig(f'{OUTDIR}/Fig1_sceptic.tif', dpi=300, bbox_inches='tight', format='tiff',
            facecolor='white', edgecolor='none')
fig.savefig(f'{OUTDIR}/Fig1_sceptic.png', dpi=300, bbox_inches='tight', format='png',
            facecolor='white', edgecolor='none')
plt.close()
print(f"Saved: {OUTDIR}/Fig1_sceptic.tif/png")
