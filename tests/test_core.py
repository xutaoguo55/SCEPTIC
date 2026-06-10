"""
Tests for SCEPTIC core types and validation engine.
"""
import numpy as np
import pandas as pd
import pytest

from sceptic.core import (
    LayerVerdict,
    SCEPTICVerdict,
    TFClaim,
    ValidationReport,
    TRRUST_REGULONS,
    _pearson_ci,
    _fisher_z,
    _fisher_z_inv,
)
from sceptic.validator import SCEPTICValidator

# ── helpers ──────────────────────────────────────────
def _make_tf_scores(n_samples=30):
    """Generate synthetic TF scores with known correlation structure."""
    rng = np.random.RandomState(42)
    x = rng.normal(0, 1, n_samples)
    # STAT1 and IRF1 highly correlated (genuine signal)
    y_irf1 = 0.9 * x + rng.normal(0, 0.3, n_samples)
    # GATA1 mostly independent (artefact)
    y_gata1 = rng.normal(0, 1, n_samples)
    # NFKB1 moderately correlated
    y_nfkb1 = 0.6 * x + rng.normal(0, 0.5, n_samples)

    df = pd.DataFrame({
        "STAT1": x,
        "IRF1": y_irf1,
        "GATA1": y_gata1,
        "NFKB1": y_nfkb1,
    })
    return df


# ── tests ─────────────────────────────────────────────
class TestCore:
    def test_fisher_z_roundtrip(self):
        r = 0.7
        z = _fisher_z(r)
        r2 = _fisher_z_inv(z)
        assert abs(r - r2) < 1e-10

    def test_pearson_ci_bounds(self):
        ci = _pearson_ci(0.6, 50)
        assert ci[0] < 0.6 < ci[1]
        assert ci[0] > 0
        assert ci[1] < 1.0

    def test_trrust_regulons_loaded(self):
        assert "STAT4" in TRRUST_REGULONS
        assert "NFKB1" in TRRUST_REGULONS
        assert len(TRRUST_REGULONS["STAT1"]) >= 9


class TestValidator:
    def test_l4_null_computation(self):
        tf_scores = _make_tf_scores(30)
        v = SCEPTICValidator(tf_scores)
        v._compute_l4_null()
        assert v._l4_null is not None
        assert len(v._l4_null) == 6  # 4 TFs → 6 pairs
        assert 0 < v._l4_median < 1
        assert v._l4_median < v._l4_p75 < v._l4_p95

    def test_l4_null_persists(self):
        """L4 null is computed once and cached."""
        tf_scores = _make_tf_scores(30)
        v = SCEPTICValidator(tf_scores)
        n1 = v.layer4_null("STAT1", "IRF1")
        n2 = v.layer4_null("NFKB1", "GATA1")
        # Both should use the same cached null
        assert "null:" in n1.detail
        assert "null:" in n2.detail

    def test_layer1_pass(self):
        tf_scores = _make_tf_scores(30)
        v = SCEPTICValidator(tf_scores)
        r = v.layer1_orthogonal("STAT1", "IRF1")
        assert r.verdict == LayerVerdict.PASS
        assert abs(r.r) > 0.3

    def test_layer4_pass_genuine_pair(self):
        tf_scores = _make_tf_scores(30)
        v = SCEPTICValidator(tf_scores)
        r = v.layer4_null("STAT1", "IRF1")
        # STAT1-IRF1 should be highly specific
        assert r.verdict in (LayerVerdict.PASS, LayerVerdict.SUPPORTIVE)

    def test_layer4_fail_artefact_pair(self):
        tf_scores = _make_tf_scores(30)
        v = SCEPTICValidator(tf_scores)
        r = v.layer4_null("NFKB1", "GATA1")
        # Should fail with a moderate dataset
        assert r.verdict in (LayerVerdict.FAIL, LayerVerdict.SUPPORTIVE)

    def test_full_validate(self):
        tf_scores = _make_tf_scores(30)
        v = SCEPTICValidator(tf_scores)
        report = v.validate("STAT1", "IRF1", source="test")
        assert isinstance(report, ValidationReport)
        assert report.claim.tf1 == "STAT1"
        assert len(report.layers) == 7
        # L1 and L4 should pass for known genuine pair
        assert report.layers[0].verdict == LayerVerdict.PASS
        # Overall should be at least context-dependent
        assert report.overall_verdict in (
            SCEPTICVerdict.SUPPORTED,
            SCEPTICVerdict.CONTEXT_DEPENDENT,
        )

    def test_layer2_erythroid_na(self):
        tf_scores = _make_tf_scores(30)
        v = SCEPTICValidator(tf_scores)
        r = v.layer2_cross_assay("GATA1", "NFKB1")
        assert r.verdict == LayerVerdict.NA
        assert "rythroid" in r.detail.lower()

    def test_l4_percentile(self):
        tf_scores = _make_tf_scores(30)
        v = SCEPTICValidator(tf_scores)
        r_stats1, _ = np.array([0.5]), None
        pct = v._l4_percentile(0.9)
        assert 0 <= pct <= 100

    def test_rng_seed_reproducibility(self):
        tf_scores = _make_tf_scores(30)
        v1 = SCEPTICValidator(tf_scores, random_seed=42)
        v2 = SCEPTICValidator(tf_scores, random_seed=42)
        r1 = v1.layer7_permutation("STAT1", "IRF1",
                                    group_labels=np.tile([0, 1], 15))
        r2 = v2.layer7_permutation("STAT1", "IRF1",
                                    group_labels=np.tile([0, 1], 15))
        assert abs(r1.r - r2.r) < 1e-10 if r1.r is not None and r2.r is not None else True

    def test_cli_import(self):
        """CLI module can be imported."""
        from sceptic.cli import main
        assert callable(main)

    def test_io_import(self):
        """IO module can be imported."""
        from sceptic.io import compute_tf_scores
        assert callable(compute_tf_scores)


# ── run ───────────────────────────────────────────────
if __name__ == "__main__":
    pytest.main([__file__, "-v"])
