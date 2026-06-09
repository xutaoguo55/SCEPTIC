# SCEPTIC

**Single-Cell Expression Program Testing and Integrity Checklist**

A seven-layer null-control framework for validating single-cell transcription factor co-expression claims.

## Overview

SCEPTIC provides a standardised, quantitative framework for assessing whether per-sample TF co-expression claims are distinguishable from background co-variation in scRNA-seq data.

### Requirements

Python 3.9+ with: `numpy scipy pandas matplotlib scanpy`
R 4.5+ with: `metafor`

### Repository Contents

```
SCEPTIC/
├── scripts/
│   ├── run_all_analyses.py         # Main analysis pipeline
│   ├── gen_fig1_sceptic.py         # Framework figure
│   ├── gen_fig5_literature_audit.py # Audit figure
│   ├── gen_graphical_abstract.py   # Graphical abstract
│   └── run_literature_audit.py     # Literature audit
├── results/                        # Output directory (user-generated)
└── README.md
```

### Usage

```bash
python scripts/run_all_analyses.py --expr <expression_matrix.csv.gz> --meta <metadata.tsv> --outdir ./results
```

### Input Format

- **Expression matrix** (CSV/gzipped CSV): genes × cells, first column = gene symbols
- **Metadata** (TSV): cell barcode column + sample/donor ID column

### Data

Processed per-sample scores and intermediate data used in the manuscript are provided as Supplementary Data files accompanying the manuscript submission. Raw public datasets are available from their respective repositories (SCP548, E-MTAB-10026, GSE65682).

## Citation

Wei X, Zheng H, Jiang H, Huang J, Wei Q, Wei Y, Feng R, Guo X. SCEPTIC: a multi-layer null-control framework for validating single-cell transcription factor co-expression claims. *Under review*.

## License

MIT
