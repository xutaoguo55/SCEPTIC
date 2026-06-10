"""SCEPTIC: Single-Cell Expression Program Testing and Integrity Checklist.

A seven-layer null-control framework for validating single-cell
transcription factor co-expression claims.
"""

__version__ = "1.0.0"
__author__ = "Xutao Guo et al."

from sceptic.core import (
    SCEPTICVerdict,
    TFClaim,
    LayerVerdict,
    LayerResult,
    ValidationReport,
    TRRUST_REGULONS,
)
from sceptic.validator import (
    SCEPTICValidator,
    validate_claim,
)
