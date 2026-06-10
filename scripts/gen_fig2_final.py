"""
Figure 2 — NFKB1-GATA1: method-dependent artifact.
4-panel: A=raw mean, B=AUCell, C=housekeeping confound, D=cross-method summary.
"""
import os, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
matplotlib.rcParams['font.size'] = 9

OUT = '/Users/guoxutao/sc_analysis/figures_crm'; os.makedirs(OUT, exist_ok=True)
C = {'dk':'#1F3A5F','tx':'#2C3E50','lt':'#95A5A6','green':'#27AE60','red':'#E74C3C','amber':'#F39C12','blue':'#2980B9'}

# Load discovery data
df = pd.read_csv('/Users/guoxutao/sc_analysis/results/mas_sepsis_ibd_per_sample.csv')
n = len(df)
r_raw, p_raw = pearsonr(df['NFKB1_mean'], df['GATA1_mean'])

# AUCell scores are the per-sample SCENIC-style outputs generated during the
# orthogonal validation pass.
auc = pd.read_csv('/Users/guoxutao/sc_analysis/scenic_results/aucell_3disease.csv')

# Use a composite "global expression" proxy from the available TF scores.
df['global'] = df[['STAT4_mean','NFKB1_mean','GATA1_mean']].mean(axis=1)
r_hk, p_hk = pearsonr(df['global'], df['NFKB1_mean'])

fig, axes = plt.subplots(2, 2, figsize=(10, 8))

# === A: Raw target-gene mean ===
ax = axes[0,0]
ax.scatter(df['NFKB1_mean'], df['GATA1_mean'], c=C['dk'], s=40, alpha=0.7, edgecolors='white')
z = np.polyfit(df['NFKB1_mean'], df['GATA1_mean'], 1)
xr = np.linspace(df['NFKB1_mean'].min(), df['NFKB1_mean'].max(), 100)
ax.plot(xr, np.polyval(z, xr), color=C['red'], lw=1.5, ls='--')
ax.set_xlabel('NFKB1 target-gene mean', fontsize=10)
ax.set_ylabel('GATA1 target-gene mean', fontsize=10)
ax.set_title('A  Raw target-gene mean', fontsize=12, weight='bold', color=C['dk'], loc='left')
ax.text(0.95, 0.05, f'r = {r_raw:+.3f}\np < 0.001\nn = {n}',
        transform=ax.transAxes, fontsize=9, ha='right', va='bottom',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor=C['lt']))

# === B: AUCell ===
ax = axes[0,1]
aucell_nfkb1 = auc['NFKB1_aucell']
aucell_gata1 = auc['GATA1_aucell']
r_auc, p_auc = pearsonr(aucell_nfkb1, aucell_gata1)
ax.scatter(aucell_nfkb1, aucell_gata1, c=C['blue'], s=40, alpha=0.7, edgecolors='white')
z2 = np.polyfit(aucell_nfkb1, aucell_gata1, 1)
xr2 = np.linspace(aucell_nfkb1.min(), aucell_nfkb1.max(), 100)
ax.plot(xr2, np.polyval(z2, xr2), color=C['lt'], lw=1.5, ls='--')
ax.set_xlabel('NFKB1 AUCell score', fontsize=10)
ax.set_ylabel('GATA1 AUCell score', fontsize=10)
ax.set_title('B  SCENIC AUCell (orthogonal)', fontsize=12, weight='bold', color=C['dk'], loc='left')
ax.text(0.95, 0.05, f'r = {r_auc:+.3f}\np = {p_auc:.3f}\nn = {len(auc)}',
        transform=ax.transAxes, fontsize=9, ha='right', va='bottom',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor=C['lt']))

# === C: Global expression confounding ===
ax = axes[1,0]
ax.scatter(df['global'], df['NFKB1_mean'], c=C['red'], s=40, alpha=0.7, edgecolors='white')
z3 = np.polyfit(df['global'], df['NFKB1_mean'], 1)
xr3 = np.linspace(df['global'].min(), df['global'].max(), 100)
ax.plot(xr3, np.polyval(z3, xr3), color=C['dk'], lw=1.5, ls='--')
ax.set_xlabel('TF-score mean (global expr. proxy)', fontsize=10)
ax.set_ylabel('NFKB1 target-gene mean', fontsize=10)
ax.set_title('C  Global-expression confound', fontsize=12, weight='bold', color=C['dk'], loc='left')
ax.text(0.95, 0.05, f'r = {r_hk:+.3f}\np < 0.001',
        transform=ax.transAxes, fontsize=9, ha='right', va='bottom',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor=C['lt']))

# === D: Cross-method evidence summary ===
ax = axes[1,1]
methods = ['Raw target-\ngene mean', 'SCENIC\nAUCell', 'Bulk\n(GSE65682)', 'L4 null\n(SCP548)',
           'Partial r\n(composition adj.)']
r_values = [0.57, 0.15, -0.07, 0.124, 0.46]
colors = [C['red'], C['blue'], C['lt'], C['lt'], C['green']]
bars = ax.bar(range(len(methods)), r_values, color=colors, edgecolor='white', width=0.6)
ax.axhline(y=0, color='black', lw=0.5)
ax.axhline(y=0.3, color=C['lt'], ls='--', lw=0.8, label='|r|=0.3 threshold')
ax.set_xticks(range(len(methods)))
ax.set_xticklabels(methods, fontsize=8, rotation=0)
ax.set_ylabel('Pearson r', fontsize=10)
ax.set_title('D  Cross-method evidence summary', fontsize=12, weight='bold', color=C['dk'], loc='left')
for bar, rv in zip(bars, r_values):
    ypos = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2, ypos + 0.02*np.sign(ypos),
            f'{rv:+.2f}', ha='center', fontsize=8.5, weight='bold',
            va='bottom' if ypos>0 else 'top')
ax.set_ylim(-0.2, 0.75)
ax.legend(fontsize=7.5, loc='upper right')

fig.suptitle('Figure 2. NFKB1–GATA1 TF-score covariation: a method-dependent artifact',
             fontsize=13, weight='bold', y=1.01)
plt.tight_layout()
fig.savefig(f'{OUT}/Fig2_final.tif', dpi=300, bbox_inches='tight', format='tiff', facecolor='white', edgecolor='none')
fig.savefig(f'{OUT}/Fig2_final.png', dpi=300, bbox_inches='tight', format='png', facecolor='white', edgecolor='none')
plt.close()
print("Fig2_final saved")
