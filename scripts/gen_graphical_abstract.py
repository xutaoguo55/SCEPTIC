"""
SCEPTIC graphical abstract for Scientific Reports.
Visual summary: 7-layer framework validating TF co-expression claims.
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
    'pass':    '#27AE60',
    'fail':    '#E74C3C',
    'partial': '#F39C12',
    'dark':    '#1F3A5F',
    'text':    '#2C3E50',
    'light':   '#95A5A6',
    'box_bg':  '#F8F9FA',
    'border':  '#BDC3C7',
    'white':   '#FFFFFF',
    'accent':  '#2980B9',
}

fig, ax = plt.subplots(figsize=(7.5, 5.0))
ax.set_xlim(0, 18)
ax.set_ylim(0, 12)
ax.axis('off')

# === LEFT: 7 layers as a stacked pipeline ===
ax.text(0.5, 11.7, 'The SCEPTIC Framework', fontsize=12, weight='bold', color=C['dark'])

layers = [
    ('L1', 'Orthogonal\nscoring', C['pass']),
    ('L2', 'Cross-assay\nreplication', C['pass']),
    ('L3', 'Bulk-tissue\nreplication', C['partial']),
    ('L4', 'Random TF\ntriplet null', C['fail']),
    ('L5', 'Composition-\nadjusted r', C['pass']),
    ('L6', 'Leave-one-\ncohort-out', C['fail']),
    ('L7', 'Label\npermutation', C['pass']),
]

y_top = 11.0
box_h = 1.1
box_w = 1.6
gap = -0.02

for i, (label, name, color) in enumerate(layers):
    y = y_top - i * (box_h + gap)
    # Layer box
    ax.add_patch(FancyBboxPatch((0.5, y), box_w, box_h,
                                 boxstyle='round,pad=0.05',
                                 facecolor=C['box_bg'], edgecolor=color, linewidth=1.5))
    ax.text(0.95, y + box_h - 0.25, label, fontsize=10, weight='bold', color=color)
    ax.text(0.95, y + box_h - 0.65, name, fontsize=7.5, color=C['text'], linespacing=1.2)
    # Color dot
    ax.add_patch(plt.Circle((0.3 + box_w, y + box_h/2), 0.12, color=color))

# Arrows between layers
for i in range(6):
    y_mid = y_top - i * (box_h + gap) - gap/2
    ax.annotate('', xy=(1.3, y_mid - 0.06), xytext=(1.3, y_mid + 0.06),
                arrowprops=dict(arrowstyle='->', color=C['light'], lw=0.8))

# === MIDDLE: Input → SCEPTIC → Output ===
mid_x = 3.2

# Input box
ax.add_patch(FancyBboxPatch((mid_x, 9.5), 4.5, 1.8,
                             boxstyle='round,pad=0.1',
                             facecolor='#EBF5FB', edgecolor=C['accent'], linewidth=1.5))
ax.text(mid_x + 2.25, 11.0, 'Input', fontsize=9, weight='bold', color=C['accent'], ha='center')
ax.text(mid_x + 2.25, 10.4, 'TF co-expression claims\nfrom per-sample aggregate\nscRNA-seq scoring',
        fontsize=8, color=C['text'], ha='center', linespacing=1.4)

# Arrow
ax.annotate('', xy=(mid_x + 2.25, 9.3), xytext=(mid_x + 2.25, 9.6),
            arrowprops=dict(arrowstyle='->', color=C['dark'], lw=2))

# SCEPTIC central box
ax.add_patch(FancyBboxPatch((mid_x, 6.0), 4.5, 3.0,
                             boxstyle='round,pad=0.12',
                             facecolor=C['white'], edgecolor=C['dark'], linewidth=2.5))
ax.text(mid_x + 2.25, 8.6, 'SCEPTIC', fontsize=13, weight='bold', color=C['dark'], ha='center')
ax.text(mid_x + 2.25, 7.9, '7-layer null-control\nvalidation pipeline',
        fontsize=9, color=C['text'], ha='center', linespacing=1.3)
ax.text(mid_x + 2.25, 7.0, 'Tests: method-dependence,\ncomposition confounding,\nspecificity, cohort robustness',
        fontsize=7.5, color=C['light'], ha='center', linespacing=1.3)
ax.text(mid_x + 2.25, 6.3, 'Case study: 1,161 patients\n6 cohorts · 3 assay types',
        fontsize=7.5, color=C['accent'], ha='center', linespacing=1.3)

# Arrow down
ax.annotate('', xy=(mid_x + 2.25, 5.8), xytext=(mid_x + 2.25, 6.1),
            arrowprops=dict(arrowstyle='->', color=C['dark'], lw=2))

# Output
ax.add_patch(FancyBboxPatch((mid_x, 3.8), 4.5, 1.7,
                             boxstyle='round,pad=0.1',
                             facecolor='#E8F8F5', edgecolor=C['pass'], linewidth=1.5))
ax.text(mid_x + 2.25, 5.2, 'Output', fontsize=9, weight='bold', color=C['pass'], ha='center')
ax.text(mid_x + 2.25, 4.6, 'Evidence-matrix verdict\nper TF pair with explicit\nuncertainty quantification',
        fontsize=8, color=C['text'], ha='center', linespacing=1.4)

# === RIGHT: Case study results ===
right_x = 8.5

ax.add_patch(FancyBboxPatch((right_x, 9.3), 8.8, 2.2,
                             boxstyle='round,pad=0.1',
                             facecolor=C['white'], edgecolor=C['border'], linewidth=1.0))
ax.text(right_x + 0.3, 11.2, 'HCTL Triad Case Study', fontsize=9, weight='bold', color=C['dark'])

results = [
    ('STAT4–NFKB1', 'CONDITIONALLY\nSUPPORTED', 'PBMC contexts only;\nfails bulk + random null', C['partial']),
    ('NFKB1–GATA1', 'REJECTED', 'Method-dependent artifact;\nraw-mean scoring false positive', C['fail']),
    ('STAT4–GATA1', 'REJECTED', 'Inconsistent direction\nacross methods and cohorts', C['fail']),
]

for i, (pair, verdict, detail, color) in enumerate(results):
    y = 10.3 - i * 0.7
    # Color bar
    ax.add_patch(FancyBboxPatch((right_x + 0.3, y - 0.1), 0.08, 0.55,
                                 boxstyle='round,pad=0.01',
                                 facecolor=color, edgecolor='none'))
    ax.text(right_x + 0.6, y + 0.35, pair, fontsize=9, weight='bold', color=C['dark'], va='center')
    ax.text(right_x + 2.8, y + 0.35, verdict, fontsize=7.5, weight='bold', color=color, va='center')
    ax.text(right_x + 5.3, y + 0.35, detail, fontsize=7, color=C['text'], va='center', linespacing=1.2)

# Key message
ax.text(right_x + 0.3, 6.5, 'Key message',
        fontsize=9, weight='bold', color=C['dark'])
ax.text(right_x + 0.3, 5.9,
        '2/3 initially significant TF pairs were method-\ndependent artifacts. Raw target-gene mean\n'
        'scoring without orthogonal validation produces\nfalse-positive co-expression claims.\n\n'
        'SCEPTIC layers 1 (orthogonal scoring) and 4\n(random TF null) are recommended as minimum\n'
        'validation standards.',
        fontsize=7.5, color=C['text'], linespacing=1.3)

# Bottom banner
ax.add_patch(FancyBboxPatch((0.3, 0.3), 17.2, 0.6,
                             boxstyle='round,pad=0.05',
                             facecolor=C['dark'], edgecolor='none'))
ax.text(9.1, 0.6, 'A multi-layer null-control framework for single-cell TF co-expression claims  ·  Open-source checklist  ·  Domain-agnostic',
        fontsize=8.5, color='white', ha='center', va='center', weight='bold')

plt.tight_layout()
fig.savefig(f'{OUTDIR}/srep_graphical_abstract.tif', dpi=300, bbox_inches='tight', format='tiff',
            facecolor='white', edgecolor='none')
fig.savefig(f'{OUTDIR}/srep_graphical_abstract.png', dpi=300, bbox_inches='tight', format='png',
            facecolor='white', edgecolor='none')
plt.close()
print(f"Saved: {OUTDIR}/srep_graphical_abstract.tif/png")
