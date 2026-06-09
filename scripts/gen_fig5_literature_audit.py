"""
Figure 5 for SCEPTIC paper — Literature Audit: L4 random TF null on published claims.
"""
import os, numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUTDIR = '/Users/guoxutao/sc_analysis/figures_crm'
os.makedirs(OUTDIR, exist_ok=True)

# Load data
rdf = pd.read_csv('/Users/guoxutao/sc_analysis/results/literature_audit/literature_audit_scp548.csv')
pdf = pd.read_csv('/Users/guoxutao/sc_analysis/results/literature_audit/all_tf_pairs_scp548.csv')

# Get abs r column regardless of how CSV saved it
ar_col = [c for c in rdf.columns if '#' not in c and 'r' in c.lower() and 'tf' not in c.lower()][0]
r_abs = rdf[ar_col].values
p_abs = pdf[['r']].abs().values.flatten()

C = {'pass': '#27AE60', 'marginal': '#F39C12', 'fail': '#E74C3C',
     'null': '#BDC3C7', 'dark': '#1F3A5F', 'text': '#2C3E50'}

fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.5), gridspec_kw={'width_ratios': [1, 1.2]})
ax1, ax2 = axes
np.random.seed(42)

# Panel A: Null histogram + literature test points
ax1.hist(p_abs, bins=45, color=C['null'], edgecolor='white', alpha=0.7)
ax1.axvline(np.median(p_abs), color=C['dark'], ls='--', lw=1.5,
            label='Median ({:.3f})'.format(float(np.median(p_abs))))
ax1.axvline(np.percentile(p_abs, 95), color=C['fail'], ls='--', lw=1.5,
            label='95th ({:.3f})'.format(float(np.percentile(p_abs, 95))))

for i in range(len(rdf)):
    row = rdf.iloc[i]
    verdict = row['Verdict']
    v = float(row[ar_col])
    if 'PASS' in str(verdict):
        c, m, z = C['pass'], '^', 5
    elif 'MARGINAL' in str(verdict):
        c, m, z = C['marginal'], 's', 4
    else:
        c, m, z = C['fail'], 'o', 3
    yj = float(np.random.uniform(2, 10))
    ax1.scatter(v, yj, c=c, marker=m, s=70, zorder=z, edgecolors='white', linewidth=0.5)
    if v > 0.8 or 'HCTL' in str(row['Source']):
        ax1.annotate('{}-{}'.format(row['TF1'], row['TF2']),
                     (v, yj+1.5), fontsize=6.5, ha='center', color=c, rotation=45)

ax1.set_xlabel('|Pearson r| between per-sample TF scores', fontsize=10)
ax1.set_ylabel('# TF pairs', fontsize=10)
ax1.set_title('L4 Null Distribution (496 TF pairs) +\nliterature test pairs', fontsize=11, weight='bold')
ax1.legend(fontsize=8, loc='upper left')
ax1.set_xlim(0, 1.0)

# Panel B: Ranked bar chart
idx = np.argsort(r_abs)
r_sorted = r_abs[idx]
verdicts = rdf['Verdict'].values[idx]
labels = ['{}-{}'.format(rdf.iloc[i]['TF1'], rdf.iloc[i]['TF2']) for i in idx]

colors = []
for v in verdicts:
    if 'PASS' in str(v):
        colors.append(C['pass'])
    elif 'MARGINAL' in str(v):
        colors.append(C['marginal'])
    else:
        colors.append(C['fail'])

ax2.barh(range(len(r_sorted)), r_sorted, color=colors, edgecolor='white', height=0.7)
ax2.axvline(np.percentile(p_abs, 95), color=C['fail'], ls='--', lw=1, alpha=0.5, label='95th %ile')
ax2.axvline(np.percentile(p_abs, 75), color=C['marginal'], ls='--', lw=1, alpha=0.5, label='75th %ile')
ax2.axvline(np.median(p_abs), color=C['dark'], ls=':', lw=1, alpha=0.5, label='Median')

ax2.set_yticks(range(len(r_sorted)))
ax2.set_yticklabels(labels, fontsize=7)
ax2.set_xlabel('|Pearson r|', fontsize=10)
ax2.set_title('Literature TF Pairs Ranked by |r|\n(SCEPTIC L4, SCP548 PBMC, n=65)', fontsize=11, weight='bold')

# Source annotations
for i in range(len(r_sorted)):
    src = str(rdf.iloc[idx[i]]['Source'])[:55]
    ax2.text(r_sorted[i] + 0.01, i, src, fontsize=5, va='center', color=C['text'], alpha=0.6)

from matplotlib.lines import Line2D
legend_el = [
    Line2D([0], [0], marker='s', color='w', markerfacecolor=C['pass'], markersize=10, label='PASS (top 5%)'),
    Line2D([0], [0], marker='s', color='w', markerfacecolor=C['marginal'], markersize=10, label='MARGINAL (top 25%)'),
    Line2D([0], [0], marker='s', color='w', markerfacecolor=C['fail'], markersize=10, label='FAIL'),
]
ax2.legend(handles=legend_el + ax2.get_legend_handles_labels()[0][:3],
           fontsize=7, loc='lower right')

plt.tight_layout()
fig.savefig('{}/Fig5_literature_audit.tif'.format(OUTDIR), dpi=300, bbox_inches='tight',
            format='tiff', facecolor='white', edgecolor='none')
fig.savefig('{}/Fig5_literature_audit.png'.format(OUTDIR), dpi=300, bbox_inches='tight',
            format='png', facecolor='white', edgecolor='none')
plt.close()
print("Saved: Fig5_literature_audit.tif/png")
