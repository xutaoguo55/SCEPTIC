# SCEPTIC

**Single-Cell Expression Program Testing and Integrity Checklist**

A seven-layer null-control framework for validating single-cell transcription factor co-expression claims.

## Overview

Single-cell transcriptomic studies frequently report transcription factor (TF) co-expression from per-sample aggregate scores. However, such claims are rarely subjected to systematic null-model validation. SCEPTIC provides a standardised, quantitative framework for assessing whether a per-sample TF co-expression claim is distinguishable from background co-variation.

### The Seven Layers

| Layer | Name | What it tests |
|-------|------|---------------|
| L1 | Orthogonal scoring | Method-dependence of raw mean correlation |
| L2 | Cross-assay replication | Transcript-level finding at protein level |
| L3 | Bulk-tissue replication | Composition-driven correlation vs. cell-intrinsic |
| L4 | Random TF null | Specificity vs. random TF pairs |
| L5 | Composition-adjusted | Cell-type proportion confounding |
| L6 | Leave-one-cohort-out | Single-cohort dependence |
| L7 | Label permutation | Sample-level false positive rate |

## Repository Structure

```
SCEPTIC/
├── scripts/
│   ├── run_all_analyses.py     # Main analysis pipeline
│   ├── gen_fig1_sceptic.py     # Figure 1 generation
│   ├── gen_fig5_literature_audit.py  # Figure 5 generation
│   └── gen_graphical_abstract.py    # Graphical abstract
├── data/
│   ├── scp548_all_tf_scores.csv     # TF scores (32 TFs × 65 donors)
│   ├── all_tf_pairs_scp548.csv      # All pairwise correlations
│   ├── literature_audit_scp548.csv  # L4 audit results
│   ├── hctl_3cohort_crossvalidation.csv  # Cross-cohort validation
│   ├── discovery_tf_scores.csv      # Discovery cohort
│   ├── scp548_tf_means.csv          # SCP548 HCTL triad
│   ├── emtab_tf_means.csv           # E-MTAB-10026 HCTL triad
│   └── gse65682_scores.csv          # GSE65682 bulk
├── figures/
│   ├── Fig1_sceptic.png             # Framework overview
│   ├── Fig2_gata1_artifact.png      # NFKB1-GATA1 artifact
│   ├── Fig3_meta_analysis.png       # Meta-analysis
│   ├── Fig4_evidence_matrix.png     # Evidence matrix
│   ├── Fig5_literature_audit.png    # Literature audit
│   ├── SuppFig1_cross_cohort.png    # Cross-cohort validation
│   └── graphical_abstract.png       # Graphical abstract
└── README.md
```

## Quick Start

### Requirements

Python 3.9+ with: `numpy scipy pandas matplotlib scanpy`
R 4.5+ with: `metafor` (for meta-analysis)

### Run the Analysis

```bash
# Run full pipeline on SCP548
python scripts/run_all_analyses.py --cohort scp548

# Run on custom dataset
python scripts/run_all_analyses.py \
    --expr path/to/expression_matrix.csv.gz \
    --meta path/to/cell_metadata.tsv \
    --outdir ./results
```

### Input Format

**Expression matrix** (CSV or gzipped CSV): genes × cells, first column = gene symbols.

**Metadata** (TSV): must contain a cell barcode column matching the expression column names, and a sample/donor ID column.

## Data Sources

Public datasets used in the manuscript:

- **SCP548**: Broad Single Cell Portal, Reyes et al. Nat Med 2020
- **E-MTAB-10026**: EBI ArrayExpress, Stephenson et al. Nat Med 2021
- **GSE65682**: NCBI GEO, Scicluna et al. Lancet Respir Med 2017
- **Pienkos et al.**: Zenodo 15199328 (Olink + CyTOF)

## Citation

Wei X, Zheng H, Jiang H, Huang J, Wei Q, Wei Y, Feng R, Guo X. SCEPTIC: a multi-layer null-control framework for validating single-cell transcription factor co-expression claims. *Submitted*.

## License

MIT

## Contact

Xutao Guo (gxt827@126.com, ORCID: 0000-0001-6191-2204)
