"""
Figure 5 — Literature audit (simplified).
A: PBMC null distribution with key pairs
B: 15 pairs ranked by percentile (not |r|)
No source labels in main figure; no "FAIL" label.
"""
import os, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size'] = 9

OUT = '/Users/guoxutao/sc_analysis/figures_crm'; os.makedirs(OUT, exist_ok=True)
C = {'dk':'#1F3A5F','tx':'#2C3E50','lt':'#95A5A6','green':'#27AE60','red':'#E74C3C','amber':'#F39C12'}

# Load audit results
rdf = pd.read_csv('/Users/guoxutao/sc_analysis/results/literature_audit/literature_audit_scp548.csv')
pdf = pd.read_csv('/Users/guoxutao/sc_analysis/results/literature_audit/all_tf_pairs_scp548.csv')

# Get abs_r column
abs_col = [c for c in rdf.columns if c != 'r' and 'r' in c.lower() and 'emp' not in c.lower()][0]
abs_r_vals = rdf[abs_col].values
# Get |r| from all pairs
all_abs = pdf[['r']].abs().values.flatten()

fig, axes = plt.subplots(1, 2, figsize=(12, 5.8), gridspec_kw={'width_ratios': [1.1, 1]})

# === Panel A: Null distribution ===
ax = axes[0]
counts, bins, patches = ax.hist(all_abs, bins=45, color=C['lt'], edgecolor='white', alpha=0.5, zorder=1)
med = np.median(all_abs)
p75 = np.percentile(all_abs, 75)
p95 = np.percentile(all_abs, 95)
ax.axvline(med, color=C['dk'], ls='--', lw=1.5, label=f'Median ({med:.3f})')
ax.axvline(p75, color=C['amber'], ls='--', lw=1.2, label=f'75th ({p75:.3f})')
ax.axvline(p95, color=C['red'], ls='--', lw=1.5, label=f'95th ({p95:.3f})')

# Plot a few key literature pairs
highlight = ['RUNX1-ETS1', 'RELA-NFKB1', 'SPI1-CEBPA', 'BCL6-PRDM1', 'STAT4-NFKB1']
np.random.seed(42)
for _, row in rdf.iterrows():
    pair = f"{row['TF1']}-{row['TF2']}"
    v = row[abs_col]
    if 'PASS' in str(row['Verdict']) or 'MARGINAL' in str(row['Verdict']):
        c, m = C['green'], '^'
    elif pair in highlight:
        c, m = C['red'], 'o'
    else:
        continue  # skip non-highlighted non-pass pairs
    yj = np.random.uniform(2, 8)
    ax.scatter(v, yj, c=c, marker=m, s=70, zorder=4, edgecolors='white', lw=0.5)
    ax.annotate(pair, (v, yj+1.2), fontsize=6.5, ha='center', color=c, rotation=45)

ax.set_xlabel('|Pearson r| between per-sample TF scores', fontsize=10)
ax.set_ylabel('# TF pairs', fontsize=10)
ax.set_title('A  PBMC TF-pair null distribution\n(SCP548, 32 TFs, 496 pairs)', fontsize=11, weight='bold', color=C['dk'], loc='left')
ax.legend(fontsize=7.5, loc='upper left')

# === Panel B: Ranked by percentile ===
ax = axes[1]
# Sort by percentile descending
idx = np.argsort(rdf['Percentile'].values)[::-1]  # descending
pct_vals = rdf['Percentile'].values[idx]
labels_sorted = [f"{rdf.iloc[i]['TF1']}-{rdf.iloc[i]['TF2']}" for i in idx]
verdicts = rdf['Verdict'].values[idx]

colors = []
for v in verdicts:
    if 'PASS' in str(v):
        colors.append(C['green'])
    elif 'MARGINAL' in str(v):
        colors.append(C['amber'])
    else:
        colors.append(C['red'])

bars = ax.barh(range(len(pct_vals)), pct_vals, color=colors, edgecolor='white', height=0.65)
ax.axvline(95, color=C['green'], ls='--', lw=1.2, alpha=0.6, label='95th %ile (PASS)')
ax.axvline(75, color=C['amber'], ls='--', lw=1.2, alpha=0.6, label='75th %ile')
ax.set_yticks(range(len(pct_vals)))
ax.set_yticklabels(labels_sorted, fontsize=7.5)
ax.set_xlabel('Percentile in PBMC null', fontsize=10)
ax.set_title('B  15 literature TF pairs\nranked by PBMC null percentile', fontsize=11, weight='bold', color=C['dk'], loc='left')
ax.set_xlim(0, 105)
ax.legend(fontsize=7.5, loc='lower right')

# Add region labels inside the panel to avoid title overlap.
ax.text(0.95, 0.98, 'Top 5%\n(PASS)', transform=ax.transAxes,
        fontsize=7, ha='center', va='top', color=C['green'], weight='bold',
        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='none', alpha=0.85))
ax.text(0.78, 0.98, 'Top 25%\n(MARGINAL)', transform=ax.transAxes,
        fontsize=7, ha='center', va='top', color=C['amber'], weight='bold',
        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='none', alpha=0.85))
ax.text(0.43, 0.98, '<75th %ile\n(not PBMC-specific)', transform=ax.transAxes,
        fontsize=7, ha='center', va='top', color=C['red'], weight='bold',
        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='none', alpha=0.85))

fig.suptitle('Figure 5. SCEPTIC L4 literature audit (SCP548 PBMC context)',
             fontsize=13, weight='bold', y=0.985)
fig.text(0.72, 0.025,
         'L4 tests PBMC context-specificity of TF-score covariation; it does not invalidate TF biology in other contexts.',
         fontsize=7.5, ha='center', color=C['tx'],
         bbox=dict(boxstyle='round,pad=0.35', facecolor='#FFF3E0', edgecolor=C['amber'], alpha=0.9))
plt.tight_layout(rect=[0, 0.06, 1, 0.95])
fig.savefig(f'{OUT}/Fig5_final.tif', dpi=300, bbox_inches='tight', format='tiff', facecolor='white', edgecolor='none')
fig.savefig(f'{OUT}/Fig5_final.png', dpi=300, bbox_inches='tight', format='png', facecolor='white', edgecolor='none')
plt.close()
print("Fig5_final saved")
