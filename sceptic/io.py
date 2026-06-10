"""
SCEPTIC I/O: data loading utilities.
"""
from typing import Dict, List, Optional, Tuple
import csv
import gzip
from collections import defaultdict
import numpy as np
import pandas as pd

from sceptic.core import TRRUST_REGULONS


def load_expression_matrix(
    path: str,
    compression: str = "gzip",
    gene_col: int = 0,
) -> pd.DataFrame:
    """Load a gene × cell expression matrix.

    Parameters
    ----------
    path : str
        CSV or gzipped CSV with genes as rows and cells as columns.
    compression : str
        "gzip" or "none"
    gene_col : int
        Column index for gene names.

    Returns
    -------
    pd.DataFrame with genes as index and cell barcodes as columns.
    """
    open_func = gzip.open if compression == "gzip" else open
    with open_func(path, "rt") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        gene_names = []
        rows = []
        for row in reader:
            gene_names.append(row[gene_col])
            rows.append([float(v) for v in row[gene_col + 1:]])

    barcodes = header[gene_col + 1:]
    expr = pd.DataFrame(rows, index=gene_names, columns=barcodes)
    return expr


def load_metadata(
    path: str,
    sample_col: str = "donor_id",
    sep: str = "\t",
) -> pd.DataFrame:
    """Load cell-level metadata.

    Parameters
    ----------
    path : str
        TSV or CSV metadata file.
    sample_col : str
        Column with sample/donor identifier.
    sep : str
        Delimiter.

    Returns
    -------
    pd.DataFrame with cell barcodes as index.
    """
    meta = pd.read_csv(path, sep=sep)
    meta.columns = meta.columns.str.strip()

    # Skip header description row if it starts with "TYPE"
    if meta.iloc[0, 0] == "TYPE":
        meta = meta.iloc[1:]

    # Use first column as index (cell barcodes)
    meta = meta.set_index(meta.columns[0])
    return meta


def load_bulk_scores(path: str) -> pd.DataFrame:
    """Load bulk tissue TF scores."""
    return pd.read_csv(path, index_col=0)


def compute_tf_scores(
    expr: pd.DataFrame,
    meta: pd.DataFrame,
    sample_col: str = "donor_id",
    regulons: Optional[Dict[str, List[str]]] = None,
    min_targets: int = 3,
) -> pd.DataFrame:
    """Compute per-sample TF target-gene mean scores.

    Parameters
    ----------
    expr : pd.DataFrame
        Gene × cell expression matrix.
    meta : pd.DataFrame
        Cell-level metadata.
    sample_col : str
        Column with sample identifier.
    regulons : dict, optional
        TF → target gene list. Defaults to TRRUST v2.
    min_targets : int
        Minimum number of expressed targets required per TF.

    Returns
    -------
    pd.DataFrame: samples × TFs.
    """
    regulons = regulons or TRRUST_REGULONS

    # Build cell→sample map
    cell_to_sample = {}
    for bc in expr.columns:
        if bc in meta.index:
            cell_to_sample[bc] = meta.loc[bc, sample_col]

    # Aggregate per sample
    sample_groups = defaultdict(list)
    for bc, sample in cell_to_sample.items():
        sample_groups[sample].append(bc)

    samples = sorted(sample_groups.keys())

    # Expression per sample
    sample_expr = pd.DataFrame(index=samples, columns=expr.index)
    for sample, barcodes in sample_groups.items():
        shared = [b for b in barcodes if b in expr.columns]
        if shared:
            sample_expr.loc[sample] = expr[shared].mean(axis=1)

    sample_expr = sample_expr.astype(float)

    # Compute TF scores
    tf_scores = pd.DataFrame(index=samples)
    for tf, targets in regulons.items():
        available = [t for t in targets if t in sample_expr.columns]
        if len(available) >= min_targets:
            tf_scores[tf] = sample_expr[available].mean(axis=1)

    return tf_scores


def load_pairs_csv(path: str) -> List[Tuple[str, str, str]]:
    """Load literature pairs from CSV."""
    df = pd.read_csv(path)
    pairs = []
    for _, row in df.iterrows():
        pairs.append((
            str(row["TF1"]),
            str(row["TF2"]),
            str(row.get("Source", "")),
        ))
    return pairs
