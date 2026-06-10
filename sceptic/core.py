"""
SCEPTIC core: data types and the main validator engine.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Tuple
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr


class LayerVerdict(Enum):
    PASS = "pass"
    FAIL = "fail"
    SUPPORTIVE = "supportive"
    NA = "na"


class SCEPTICVerdict(Enum):
    SUPPORTED = "supported"
    CONTEXT_DEPENDENT = "context_dependent"
    REJECTED = "rejected"
    INSUFFICIENT = "insufficient_evidence"


@dataclass
class LayerResult:
    layer: int
    name: str
    verdict: LayerVerdict
    r: Optional[float] = None
    p: Optional[float] = None
    ci_95: Optional[Tuple[float, float]] = None
    detail: str = ""


@dataclass
class TFClaim:
    tf1: str
    tf2: str
    source: str = ""


@dataclass
class ValidationReport:
    claim: TFClaim
    layers: List[LayerResult] = field(default_factory=list)
    overall_verdict: SCEPTICVerdict = SCEPTICVerdict.INSUFFICIENT
    summary: str = ""


# ============================================================
# TRRUST v2 regulons (curated from Han et al. 2018)
# ============================================================
TRRUST_REGULONS: Dict[str, List[str]] = {
    "STAT4": ["IFNG", "IL2", "IL12RB2", "TBX21", "CCR5", "CXCR3",
              "IL18R1", "IL18RAP", "HAVCR2"],
    "NFKB1": ["TNF", "IL6", "CXCL8", "IL1B", "CCL2", "CCL3", "CCL4",
              "CCL5", "CXCL10", "ICAM1", "VCAM1", "SELE", "PTGS2",
              "MMP9", "BCL2", "BCL2L1", "SOD2"],
    "GATA1": ["HBB", "HBA1", "HBA2", "ALAS2", "EPB42", "GYPA", "GYPB",
              "KLF1", "NFE2", "SP1"],
    "IRF4": ["PRDM1", "AICDA", "BCL6", "XBP1"],
    "BATF": ["PRDM1", "AICDA", "BCL6", "IL21", "IL4", "IL10"],
    "STAT1": ["IRF1", "IRF7", "IRF9", "MX1", "OAS1", "IFIT1", "IFIT3",
              "ISG15", "GBP1", "CIITA", "TAP1", "B2M", "HLA-A", "HLA-B"],
    "IRF1": ["IFNB1", "CXCL10", "GBP2", "TAP1", "PSMB9", "B2M"],
    "SPI1": ["CSF1R", "CD14", "FCGR1A", "ITGAM", "ITGB2", "TLR4",
             "MPO", "ELANE", "CEBPA"],
    "CEBPA": ["CSF3R", "MPO", "ELANE", "CEBPE"],
    "CEBPB": ["IL6", "TNF", "SAA1", "SAA2", "HP", "FGG", "FGA", "FGB"],
    "FOXP3": ["IL2RA", "CTLA4", "TIGIT", "IKZF2", "IL10", "TGFB1"],
    "TBX21": ["IFNG", "CXCR3", "IL12RB2", "PRF1", "GZMB"],
    "GATA3": ["IL4", "IL5", "IL13"],
    "BCL6": ["PRDM1", "IRF4", "CD69", "XBP1"],
    "PRDM1": ["IRF4", "BCL6", "PAX5", "CIITA"],
    "RELA": ["TNF", "IL6", "CXCL8", "CCL2", "ICAM1", "VCAM1", "BCL2",
             "BCL2L1", "MMP9", "PTGS2", "IL1B", "NFKBIA", "SOD2"],
    "HIF1A": ["VEGFA", "SLC2A1", "LDHA", "HK2", "PKM", "BNIP3", "EPO",
              "TFRC", "IGF2"],
    "MYC": ["CCND2", "CDK4", "E2F1", "CDKN1A", "ODC1", "LDHA", "CAD", "TP53"],
    "TP53": ["CDKN1A", "BAX", "BBC3", "PMAIP1", "GADD45A", "MDM2",
             "RRM2B", "SESN1"],
    "FOS": ["CCND1", "FOSL1", "JUNB", "MMP1", "MMP3"],
    "JUN": ["CCND1", "FOS", "FOSL1", "JUNB", "MMP1", "MMP3", "TP53"],
    "NR3C1": ["NFKBIA", "IL10", "DUSP1", "FKBP5", "SGK1", "ANXA1"],
    "ETS1": ["CD3E", "CD247", "GZMB", "PRF1"],
    "RUNX1": ["CD4", "CD8A", "CD3E", "CD247", "ITGAL", "GZMB", "PRF1"],
    "RUNX3": ["CD8A", "CD8B", "GZMB", "PRF1", "IFNG"],
    "PAX5": ["CD19", "CD79A", "CD79B", "BLK", "BLNK", "CD22"],
    "IRF8": ["CIITA", "B2M", "TAP1", "HLA-DRA", "HLA-DRB1"],
    "E2F1": ["CCNA2", "CCNE1", "MCM2", "MCM3", "MCM4", "DHFR", "TK1", "TYMS"],
    "NFE2L2": ["NQO1", "GCLC", "GCLM", "TXNRD1", "HMOX1", "SLC7A11",
               "SOD1", "CAT", "PRDX1"],
    "TCF7": ["LEF1", "AXIN2", "MYC", "CCND1", "CD44"],
    "LEF1": ["CCND1", "MYC", "AXIN2"],
    "SMAD3": ["SERPINE1", "COL1A1", "FN1", "TGFB1"],
}

# Lineage-specific regulons that may be confounded by contamination
ERYTHROID_REGULONS = {"GATA1"}


# ============================================================
# Utility: Fisher z transform and CI
# ============================================================
def _fisher_z(r: float) -> float:
    return np.arctanh(np.clip(r, -0.9999, 0.9999))


def _fisher_z_inv(z: float) -> float:
    return np.tanh(z)


def _pearson_ci(r: float, n: int, alpha: float = 0.05) -> Tuple[float, float]:
    z = _fisher_z(r)
    se = 1.0 / np.sqrt(n - 3)
    z_crit = 1.96  # for 95% CI
    z_lo, z_hi = z - z_crit * se, z + z_crit * se
    return (_fisher_z_inv(z_lo), _fisher_z_inv(z_hi))
