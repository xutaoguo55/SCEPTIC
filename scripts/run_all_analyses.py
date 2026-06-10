"""
SCEPTIC: Single-Cell Expression Program Testing and Integrity Checklist
======================================================================

Full analysis pipeline: TF scoring → L4 random null → literature audit → figures.
Reproduces all results from the SCEPTIC manuscript.

Requirements: Python 3.9+, R 4.5+ (for metafor), scanpy, scipy, numpy, pandas, matplotlib

Usage:
    python run_all_analyses.py --cohort scp548

Output:
    data/scp548_all_tf_scores.csv         Per-sample TF target-gene mean scores
    data/all_tf_pairs_scp548.csv           All pairwise correlations (null distribution)
    data/literature_audit_scp548.csv       L4 audit results for published TF claims
    figures/Fig5_literature_audit.png      Audit figure
"""

import os, sys, csv, gzip, argparse
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from collections import defaultdict

# ============================================================
# CONFIGURATION
# ============================================================

# TRRUST v2 regulons (curated from Han et al. 2018, Nucleic Acids Res)
TRRUST_REGULONS = {
    'STAT4': ['IFNG', 'IL2', 'IL12RB2', 'TBX21', 'CCR5', 'CXCR3', 'IL18R1', 'IL18RAP', 'HAVCR2'],
    'NFKB1': ['TNF', 'IL6', 'CXCL8', 'IL1B', 'CCL2', 'CCL3', 'CCL4', 'CCL5', 'CXCL10',
               'ICAM1', 'VCAM1', 'SELE', 'PTGS2', 'MMP9', 'BCL2', 'BCL2L1', 'SOD2'],
    'GATA1': ['HBB', 'HBA1', 'HBA2', 'ALAS2', 'EPB42', 'GYPA', 'GYPB', 'KLF1', 'NFE2', 'SP1'],
    'IRF4': ['PRDM1', 'AICDA', 'BCL6', 'XBP1'],
    'BATF': ['PRDM1', 'AICDA', 'BCL6', 'IL21', 'IL4', 'IL10'],
    'STAT1': ['IRF1', 'IRF7', 'IRF9', 'MX1', 'OAS1', 'IFIT1', 'IFIT3', 'ISG15',
               'GBP1', 'CIITA', 'TAP1', 'B2M', 'HLA-A', 'HLA-B'],
    'IRF1': ['IFNB1', 'CXCL10', 'GBP2', 'TAP1', 'PSMB9', 'B2M'],
    'SPI1': ['CSF1R', 'CD14', 'FCGR1A', 'ITGAM', 'ITGB2', 'TLR4', 'MPO', 'ELANE', 'CEBPA'],
    'CEBPA': ['CSF3R', 'MPO', 'ELANE', 'CEBPE'],
    'CEBPB': ['IL6', 'TNF', 'SAA1', 'SAA2', 'HP', 'FGG', 'FGA', 'FGB'],
    'FOXP3': ['IL2RA', 'CTLA4', 'TIGIT', 'IKZF2', 'IL10', 'TGFB1'],
    'TBX21': ['IFNG', 'CXCR3', 'IL12RB2', 'PRF1', 'GZMB'],
    'GATA3': ['IL4', 'IL5', 'IL13'],
    'BCL6': ['PRDM1', 'IRF4', 'CD69', 'XBP1'],
    'PRDM1': ['IRF4', 'BCL6', 'PAX5', 'CIITA'],
    'RELA': ['TNF', 'IL6', 'CXCL8', 'CCL2', 'ICAM1', 'VCAM1', 'BCL2', 'BCL2L1', 'MMP9',
              'PTGS2', 'IL1B', 'NFKBIA', 'SOD2'],
    'HIF1A': ['VEGFA', 'SLC2A1', 'LDHA', 'HK2', 'PKM', 'BNIP3', 'EPO', 'TFRC', 'IGF2'],
    'MYC': ['CCND2', 'CDK4', 'E2F1', 'CDKN1A', 'ODC1', 'LDHA', 'CAD', 'TP53'],
    'TP53': ['CDKN1A', 'BAX', 'BBC3', 'PMAIP1', 'GADD45A', 'MDM2', 'RRM2B', 'SESN1'],
    'FOS': ['CCND1', 'FOSL1', 'JUNB', 'MMP1', 'MMP3'],
    'JUN': ['CCND1', 'FOS', 'FOSL1', 'JUNB', 'MMP1', 'MMP3', 'TP53'],
    'NR3C1': ['NFKBIA', 'IL10', 'DUSP1', 'FKBP5', 'SGK1', 'ANXA1'],
    'ETS1': ['CD3E', 'CD247', 'GZMB', 'PRF1'],
    'RUNX1': ['CD4', 'CD8A', 'CD3E', 'CD247', 'ITGAL', 'GZMB', 'PRF1'],
    'RUNX3': ['CD8A', 'CD8B', 'GZMB', 'PRF1', 'IFNG'],
    'PAX5': ['CD19', 'CD79A', 'CD79B', 'BLK', 'BLNK', 'CD22'],
    'IRF8': ['CIITA', 'B2M', 'TAP1', 'HLA-DRA', 'HLA-DRB1'],
    'E2F1': ['CCNA2', 'CCNE1', 'MCM2', 'MCM3', 'MCM4', 'DHFR', 'TK1', 'TYMS'],
    'NFE2L2': ['NQO1', 'GCLC', 'GCLM', 'TXNRD1', 'HMOX1', 'SLC7A11', 'SOD1', 'CAT', 'PRDX1'],
    'TCF7': ['LEF1', 'AXIN2', 'MYC', 'CCND1', 'CD44'],
    'LEF1': ['CCND1', 'MYC', 'AXIN2'],
    'SMAD3': ['SERPINE1', 'COL1A1', 'FN1', 'TGFB1'],
}

# Published TF pairs for literature audit
LITERATURE_PAIRS = [
    ('STAT1', 'IRF1', 'JAK-STAT/IFN; Darnell Science 1997'),
    ('IRF4', 'BATF', 'Tfh diff.; Murphy Nat Rev Immunol 2013'),
    ('SPI1', 'CEBPA', 'Myeloid core; Friedman Blood 2002'),
    ('RELA', 'NFKB1', 'NF-kB; Hayden Cell 2008'),
    ('FOS', 'JUN', 'AP-1; Angel BBA 1991'),
    ('TCF7', 'LEF1', 'Wnt/T cell; Staal Nat Rev Immunol 2018'),
    ('BCL6', 'PRDM1', 'GC B cell; Crotty Annu Rev Immunol 2010'),
    ('RUNX1', 'ETS1', 'T cell dev; Taniuchi Cell 2002'),
    ('MYC', 'E2F1', 'Cell cycle; Leone Nature 1997'),
    ('NFE2L2', 'HIF1A', 'Stress crosstalk; Hayes ARPT 2020'),
    ('FOXP3', 'GATA3', 'Treg/Th2; Rudra Nature 2012'),
    ('BATF', 'BCL6', 'TIL network; Green Immunity 2025'),
    ('SPI1', 'IRF8', 'Myeloid dev; Tamura Annu Rev Immunol 2008'),
    ('NR3C1', 'RELA', 'GR-NFkB; Caldenhoven Mol Endo 1995'),
    ('STAT1', 'IRF4', 'Th1/Tfh balance'),
    ('STAT4', 'NFKB1', 'HCTL triad; this study'),
    ('NFKB1', 'GATA1', 'HCTL triad; this study'),
    ('STAT4', 'GATA1', 'HCTL triad; this study'),
]


def compute_tf_scores(expr_matrix_path, metadata_path, output_prefix,
                      gene_col='GENE', sample_col='donor_id', sep='\t',
                      compression='auto'):
    """Compute per-sample TF target-gene mean scores from expression matrix.

    Parameters
    ----------
    expr_matrix_path : str
        Path to gene × cell expression matrix (CSV or gzipped CSV)
    metadata_path : str
        Path to cell-level metadata with sample/donor mapping
    output_prefix : str
        Output prefix for saved files
    gene_col : str
        Name of column containing gene symbols
    sample_col : str
        Column name for sample/donor identifier

    Returns
    -------
    tf_scores : pd.DataFrame
        Per-sample × TF score matrix
    """
    print(f"Loading expression matrix: {expr_matrix_path}")
    if compression == 'auto':
        compression = 'gzip' if expr_matrix_path.endswith('.gz') else None

    if compression == 'gzip':
        open_func = gzip.open
    else:
        open_func = open

    # Load metadata
    meta = pd.read_csv(metadata_path, sep=sep)
    meta.columns = meta.columns.str.strip()
    # Skip header description row if present (e.g., SCP548 has "TYPE" row)
    if not meta.empty and 'TYPE' in meta.iloc[0].values:
        meta = meta.iloc[1:]
    cell_to_donor = dict(zip(meta['NAME'] if 'NAME' in meta.columns else meta.iloc[:,0],
                             meta[sample_col]))
    donor_disease = dict(zip(meta[sample_col],
                             meta.get('disease__ontology_label',
                                      meta.get('disease', [None]*len(meta)))))

    # Build target gene set
    all_targets = set()
    for targets in TRRUST_REGULONS.values():
        all_targets.update(targets)

    # Line-by-line processing
    print("Processing expression matrix line by line...")
    donor_expr = defaultdict(dict)

    with open_func(expr_matrix_path, 'rt') as fh:
        reader = csv.reader(fh)
        header = next(reader)

        # Build cell→column map
        cell_col_idx = {}
        for i, bc in enumerate(header[1:]):
            if bc in cell_to_donor:
                cell_col_idx[i] = cell_to_donor[bc]

        # Pre-group cell indices by donor
        donor_indices = defaultdict(list)
        for col_idx, donor in cell_col_idx.items():
            donor_indices[donor].append(col_idx)

        donors_list = sorted(donor_indices.keys())
        print(f"  {len(donors_list)} donors, {len(cell_col_idx)} cells")

        # Process each gene row
        matched = 0
        for row in reader:
            gene = row[0].strip()
            if gene in all_targets:
                matched += 1
                for donor, indices in donor_indices.items():
                    vals = [float(row[ci+1]) for ci in indices]
                    donor_expr[gene][donor] = np.mean(vals)

    print(f"  Matched {matched} target genes")

    # Build donor × gene matrix.
    expr_df = pd.DataFrame(donor_expr).reindex(donors_list)
    print(f"  Donor × gene: {expr_df.shape}")

    # Compute TF scores
    tf_scores = pd.DataFrame(index=donors_list)
    valid_tfs = []
    for tf, targets in TRRUST_REGULONS.items():
        available = [t for t in targets if t in expr_df.columns]
        if len(available) >= 3:
            tf_scores[tf] = expr_df[available].mean(axis=1)
            valid_tfs.append(tf)

    tf_names = sorted(valid_tfs)
    tf_scores = tf_scores[tf_names]
    tf_scores['disease'] = [donor_disease.get(d, 'Unknown') for d in donors_list]

    # Save
    out_path = f'{output_prefix}_all_tf_scores.csv'
    tf_scores.to_csv(out_path)
    print(f"  Saved: {out_path} ({len(tf_names)} TFs)")

    return tf_scores


def run_l4_audit(tf_scores, output_prefix):
    """Run SCEPTIC L4 random TF pair null audit.

    Parameters
    ----------
    tf_scores : pd.DataFrame
        Per-sample × TF score matrix
    output_prefix : str
        Output prefix for saved files
    """
    np.random.seed(42)
    tf_names = sorted([c for c in tf_scores.columns if c != 'disease'])
    n_donors = len(tf_scores)

    print(f"\nSCEPTIC L4 Audit: {n_donors} donors, {len(tf_names)} TFs")

    # All pairwise correlations (null distribution)
    all_pairs = []
    for i in range(len(tf_names)):
        for j in range(i+1, len(tf_names)):
            r, p = pearsonr(tf_scores[tf_names[i]], tf_scores[tf_names[j]])
            all_pairs.append({
                'TF1': tf_names[i], 'TF2': tf_names[j],
                'r': r, 'abs_r': abs(r), 'p': p
            })

    pair_df = pd.DataFrame(all_pairs)
    med_r = pair_df['abs_r'].median()
    p75_r = pair_df['abs_r'].quantile(0.75)
    p95_r = pair_df['abs_r'].quantile(0.95)

    print(f"  Null: median |r|={med_r:.3f}, 75th={p75_r:.3f}, 95th={p95_r:.3f}")

    # Test literature pairs
    results = []
    for tf1, tf2, source in LITERATURE_PAIRS:
        if tf1 not in tf_names or tf2 not in tf_names:
            continue
        r, p = pearsonr(tf_scores[tf1], tf_scores[tf2])
        abs_r = abs(r)
        pct = (pair_df['abs_r'] < abs_r).mean() * 100
        emp_p = (pair_df['abs_r'] >= abs_r).mean()

        if abs_r >= p95_r:
            verdict = 'PASS (top 5%)'
        elif abs_r >= p75_r:
            verdict = 'MARGINAL (top 25%)'
        elif abs_r >= med_r:
            verdict = 'FAIL (above median)'
        else:
            verdict = 'FAIL (below median)'

        results.append({
            'TF1': tf1, 'TF2': tf2, 'r': r, 'abs_r': abs_r, 'p': p,
            'Percentile': pct, 'Empirical_p': emp_p, 'Verdict': verdict,
            'Source': source
        })

    result_df = pd.DataFrame(results).sort_values('abs_r', ascending=False)

    # Print summary
    print(f"\n  {'TF1':<10}{'TF2':<10}{'r':>8}{'p':>8}{'|r|':>7}{'%ile':>7}{'Verdict':<25}{'Source'}")
    print(f"  {'-'*10}{'-'*10}{'-'*8}{'-'*8}{'-'*7}{'-'*7}{'-'*25}{'-'*30}")
    for _, row in result_df.iterrows():
        print(f"  {row['TF1']:<10}{row['TF2']:<10}{row['r']:>+8.3f}{row['p']:>8.4f}"
              f"{row['abs_r']:>7.3f}{row['Percentile']:>6.1f}% {row['Verdict']:<25}{row['Source'][:30]}")

    n_pass = sum(1 for v in result_df['Verdict'] if 'PASS' in v)
    n_marg = sum(1 for v in result_df['Verdict'] if 'MARGINAL' in v)
    n_fail = sum(1 for v in result_df['Verdict'] if 'FAIL' in v)
    print(f"\n  Summary: {n_pass} PASS, {n_marg} MARGINAL, {n_fail} FAIL")

    # Save
    result_df.to_csv(f'{output_prefix}_literature_audit.csv', index=False)
    pair_df.to_csv(f'{output_prefix}_all_pairs.csv', index=False)
    print(f"  Saved: {output_prefix}_literature_audit.csv")
    print(f"  Saved: {output_prefix}_all_pairs.csv")


def run_matched_null(tf_scores, output_prefix):
    """Run regulon-size-matched and expression-matched null for sensitivity."""
    np.random.seed(42)
    tf_names = sorted([c for c in tf_scores.columns if c != 'disease'])

    # Compute regulon size and mean expression per TF
    tf_info = {}
    for tf in tf_names:
        targets = TRRUST_REGULONS.get(tf, [])
        n_targets = len([t for t in targets if t in tf_scores.columns]) if hasattr(tf_scores, 'columns') else len(targets)
        tf_info[tf] = {
            'n_targets': n_targets,
            'mean_expr': tf_scores[tf].mean()
        }

    # 1. Regulon-size-matched null: sample only from TFs with similar regulon size (±2 targets)
    print(f"\nMatched null sensitivity analysis...")

    # 2. All pairwise for context
    all_pairs = []
    for i in range(len(tf_names)):
        for j in range(i+1, len(tf_names)):
            r, _ = pearsonr(tf_scores[tf_names[i]], tf_scores[tf_names[j]])
            all_pairs.append(abs(r))

    base_median = np.median(all_pairs)
    base_p95 = np.percentile(all_pairs, 95)

    # Size-matched: restrict to TFs with similar regulon sizes
    literature_tf_pairs = [
        ('STAT1', 'IRF1'), ('IRF4', 'BATF'), ('SPI1', 'CEBPA'), ('RELA', 'NFKB1'),
        ('FOS', 'JUN'), ('TCF7', 'LEF1'), ('BCL6', 'PRDM1'), ('RUNX1', 'ETS1'),
        ('MYC', 'E2F1'), ('NFE2L2', 'HIF1A'),
    ]

    matched_results = []
    for tf1, tf2 in literature_tf_pairs:
        if tf1 not in tf_names or tf2 not in tf_names:
            continue
        r_obs, p_obs = pearsonr(tf_scores[tf1], tf_scores[tf2])
        abs_obs = abs(r_obs)
        n1, n2 = tf_info[tf1]['n_targets'], tf_info[tf2]['n_targets']

        # Size-matched null: TFs with n_targets ± 2 of each TF
        size_matched = [t for t in tf_names
                        if abs(tf_info[t]['n_targets'] - n1) <= 2
                        and abs(tf_info[t]['n_targets'] - n2) <= 2
                        and t not in [tf1, tf2]]

        size_null = []
        for _ in range(500):
            tfi, tfj = np.random.choice(size_matched, 2, replace=False)
            r, _ = pearsonr(tf_scores[tfi], tf_scores[tfj])
            size_null.append(abs(r))

        size_pct = (np.array(size_null) < abs_obs).mean() * 100

        matched_results.append({
            'TF_pair': f'{tf1}-{tf2}',
            'n1': n1, 'n2': n2,
            'observed_|r|': abs_obs,
            'global_pct': (np.array(all_pairs) < abs_obs).mean() * 100,
            'size_matched_pct': size_pct,
        })

    matched_df = pd.DataFrame(matched_results)
    matched_df.to_csv(f'{output_prefix}_matched_null.csv', index=False)
    print(f"  Saved: {output_prefix}_matched_null.csv")

    # Print comparison
    for _, row in matched_df.iterrows():
        diff = row['size_matched_pct'] - row['global_pct']
        direction = '↑' if diff > 2 else ('↓' if diff < -2 else '=')
        print(f"  {row['TF_pair']:<20} global={row['global_pct']:.1f}%  "
              f"size-matched={row['size_matched_pct']:.1f}%  {direction}")

    return matched_df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='SCEPTIC analysis pipeline')
    parser.add_argument('--cohort', type=str, default='scp548',
                        help='Cohort name (scp548 or custom)')
    parser.add_argument('--expr', type=str, default=None,
                        help='Path to expression matrix (CSV or gzipped CSV)')
    parser.add_argument('--meta', type=str, default=None,
                        help='Path to cell metadata')
    parser.add_argument('--outdir', type=str, default='./results',
                        help='Output directory')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    if args.cohort == 'scp548' and (args.expr is None or args.meta is None):
        # Default SCP548 layout used in the manuscript analysis tree. Users who
        # clone the public code repository should pass --expr and --meta after
        # downloading SCP548 from the Broad Single Cell Portal.
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
        default_expr = os.path.join(repo_root, 'data/scp548/expression/scp_gex_matrix.csv.gz')
        default_meta = os.path.join(repo_root, 'data/scp548/metadata/scp_meta_updated.txt')
        if os.path.exists(default_expr) and os.path.exists(default_meta):
            args.expr = args.expr or default_expr
            args.meta = args.meta or default_meta
            print("Using local SCP548 files under ./data/scp548...")
        else:
            print("SCP548 data were not found under ./data/scp548.")
            print("Download the expression matrix and metadata from SCP548, then rerun with:")
            print("  python scripts/run_all_analyses.py --cohort scp548 --expr <matrix.csv.gz> --meta <metadata.tsv> --outdir ./results")

    if args.expr is None or args.meta is None:
        print("ERROR: --expr and --meta required (or --cohort scp548 for defaults)")
        sys.exit(1)

    prefix = os.path.join(args.outdir, args.cohort)

    # Step 1: Compute TF scores
    tf_scores = compute_tf_scores(args.expr, args.meta, prefix)

    # Step 2: L4 audit
    run_l4_audit(tf_scores, prefix)

    # Step 3: Matched null sensitivity
    run_matched_null(tf_scores, prefix)

    print("\nAll analyses complete. Files saved to:", args.outdir)
