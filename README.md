# SCEPTIC

**Single-Cell Expression Program Testing and Integrity Checklist**

A seven-layer null-control framework for validating single-cell transcription factor co-expression claims.

[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

## Overview

SCEPTIC is not a method for proving TF co-regulation; it is a **reporting and robustness framework** for determining whether per-sample TF-score covariation is method-consistent, context-replicable, and distinguishable from background TF co-variation.

## The Seven Layers

| Layer | Name | What it tests |
|-------|------|---------------|
| L1 | Orthogonal scoring | Method-dependence of raw target-gene mean correlation |
| L2 | Cross-assay | Protein-level pathway support (supportive only) |
| L3 | Bulk-tissue | Cross-platform / bulk-context consistency |
| L4 | Random TF-pair null | Specificity vs. random TF-pair background |
| L5 | Composition adjustment | Cell-type proportion confounding |
| L6 | Leave-one-cohort-out | Single-cohort dependence |
| L7 | Label permutation | Disease-group dependence |

## Quick Start

### Installation

```bash
pip install sceptic
```

Or from source:

```bash
git clone https://github.com/xutaoguo55/SCEPTIC.git
cd SCEPTIC
pip install -e ".[all]"
```

### Docker

```bash
docker build -t sceptic -f docker/Dockerfile .
docker run sceptic --help
```

### Command Line

```bash
# Validate a single TF-pair claim
sceptic validate \
    --expr expression_matrix.csv.gz \
    --meta cell_metadata.tsv \
    --tf1 STAT4 --tf2 NFKB1 \
    --source "HCTL triad case study"

# Run L4 literature audit
sceptic audit \
    --scores tf_scores.csv \
    --pairs literature_pairs.csv \
    --output audit_results.csv
```

### Python API

```python
from sceptic import SCEPTICValidator, validate_claim
import pandas as pd

# Load per-sample TF scores
tf_scores = pd.read_csv("tf_scores.csv", index_col=0)

# Create validator
validator = SCEPTICValidator(tf_scores, random_seed=42)

# Run full validation
report = validator.validate("STAT4", "NFKB1", source="HCTL triad")

print(report.summary)
# → "STAT4-NFKB1: passes L1, fails L4; context-dependent covariation only"
```

## Input Data Format

### Expression Matrix (CSV / gzipped CSV)
Genes × cells, first column = gene symbols:

```
GENE,AAACCTGAGACGCTTT,AAACCTGAGCAATCTC,...
STAT4,0.0,2.3,...
NFKB1,1.5,0.0,...
```

### Metadata (TSV)
Must contain cell barcode (matching expression columns) and sample/donor ID:

```
NAME	donor_id	disease
AAACCTGAGACGCTTT	P1	sepsis
AAACCTGAGCAATCTC	P1	sepsis
```

### Pre-computed TF Scores (for audit mode)
Samples × TFs, with sample IDs as index:

```
sample,STAT4,NFKB1,GATA1,STAT1,...
P1,0.0378,0.0737,0.1240,...
P2,0.1248,0.1516,0.0679,...
```

## Citation

Wei X, Zheng H, Jiang H, Huang J, Wei Q, Wei Y, Feng R, Guo X. SCEPTIC: a multi-layer null-control framework for validating single-cell transcription factor co-expression claims. *Cell Reports Methods* (submitted).

## License

MIT © 2026 Xutao Guo et al.
