"""
SCEPTIC validator: the 7-layer validation engine.
"""
from __future__ import annotations
import warnings
from typing import Dict, List, Optional, Tuple
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from sceptic.core import (
    LayerVerdict,
    SCEPTICVerdict,
    LayerResult,
    TFClaim,
    ValidationReport,
    TRRUST_REGULONS,
    ERYTHROID_REGULONS,
    _fisher_z,
    _pearson_ci,
)

# ============================================================
# SCEPTICValidator class
# ============================================================
class SCEPTICValidator:
    """Seven-layer validation engine for TF co-expression claims."""

    def __init__(
        self,
        tf_scores: pd.DataFrame,
        regulons: Optional[Dict[str, List[str]]] = None,
        random_seed: int = 42,
        n_perms: int = 10000,
    ):
        """
        Parameters
        ----------
        tf_scores : pd.DataFrame
            Per-sample × TF-score matrix. Rows = samples, columns = TFs.
        regulons : dict, optional
            TF → target gene list mapping. Defaults to TRRUST v2.
        random_seed : int
        n_perms : int
            Number of permutations for L7.
        """
        self.tf_scores = tf_scores.copy()
        self.tf_names = sorted([c for c in tf_scores.columns])
        self.regulons = regulons or TRRUST_REGULONS
        self.random_seed = random_seed
        self.n_perms = n_perms
        self.rng = np.random.RandomState(random_seed)

        # Precompute L4 null distribution
        self._l4_null = None
        self._l4_median = None
        self._l4_p75 = None
        self._l4_p95 = None

    # ================================================================
    # L4 null (shared across all claims)
    # ================================================================
    def _compute_l4_null(self):
        """Compute all-vs-all pairwise |r| empirical null."""
        if self._l4_null is not None:
            return
        tfs = self.tf_names
        n = len(tfs)
        pairs = []
        for i in range(n):
            for j in range(i + 1, n):
                r, _ = pearsonr(self.tf_scores[tfs[i]], self.tf_scores[tfs[j]])
                pairs.append(abs(r))
        self._l4_null = np.array(pairs)
        self._l4_median = float(np.median(self._l4_null))
        self._l4_p75 = float(np.percentile(self._l4_null, 75))
        self._l4_p95 = float(np.percentile(self._l4_null, 95))

    def _l4_percentile(self, abs_r: float) -> float:
        self._compute_l4_null()
        return float((self._l4_null < abs_r).mean() * 100)

    def _l4_empirical_p(self, abs_r: float) -> float:
        self._compute_l4_null()
        return float((self._l4_null >= abs_r).mean())

    # ================================================================
    # Individual layers
    # ================================================================
    def layer1_orthogonal(
        self, tf1: str, tf2: str, scores_method2: Optional[pd.DataFrame] = None
    ) -> LayerResult:
        """L1: compare raw target-gene mean with orthogonal scoring method."""
        if tf1 not in self.tf_scores.columns or tf2 not in self.tf_scores.columns:
            return LayerResult(1, "Orthogonal scoring", LayerVerdict.NA,
                               detail=f"{tf1} or {tf2} not found in scores")

        r, p = pearsonr(self.tf_scores[tf1], self.tf_scores[tf2])
        ci = _pearson_ci(r, len(self.tf_scores))
        n = len(self.tf_scores)

        if scores_method2 is not None and tf1 in scores_method2.columns and tf2 in scores_method2.columns:
            r2, p2 = pearsonr(scores_method2[tf1], scores_method2[tf2])
            same_sign = (r * r2) >= 0
            both_sig = (p < 0.05 and p2 < 0.05)
            abs_r_ok = abs(r) > 0.3

            if same_sign and both_sig and abs_r_ok:
                verdict = LayerVerdict.PASS
                detail = (f"r_raw={r:+.3f}, p_raw={p:.4f}; "
                          f"r_orth={r2:+.3f}, p_orth={p2:.4f}; "
                          f"CI95=[{ci[0]:+.3f}, {ci[1]:+.3f}]; n={n}")
            elif not same_sign:
                verdict = LayerVerdict.FAIL
                detail = (f"Direction mismatch: raw r={r:+.3f} vs orth r={r2:+.3f}")
            else:
                verdict = LayerVerdict.FAIL
                detail = (f"Not robust: r_raw={r:+.3f}, r_orth={r2:+.3f}; "
                          f"one or both conditions fail")
        else:
            # Only raw mean available
            if abs(r) > 0.3 and p < 0.05:
                verdict = LayerVerdict.PASS
                detail = f"r={r:+.3f}, p={p:.4f}, |r|>0.3; orthogonal method not provided"
            else:
                verdict = LayerVerdict.FAIL
                detail = f"r={r:+.3f}, p={p:.4f}; fails |r|>0.3 or p<0.05"

        return LayerResult(1, "Orthogonal scoring", verdict, r=r, p=p, ci_95=ci, detail=detail)

    def layer2_cross_assay(
        self, tf1: str, tf2: str, protein_r: Optional[float] = None,
        protein_p: Optional[float] = None, protein_label: str = "",
    ) -> LayerResult:
        """L2: cross-assay protein-level support."""
        if any(r in ERYTHROID_REGULONS for r in [tf1, tf2]):
            return LayerResult(2, "Cross-assay", LayerVerdict.NA,
                               detail=f"Erythroid regulon; protein-level support N/A in PBMC/whole blood")

        if protein_r is None:
            return LayerResult(2, "Cross-assay", LayerVerdict.NA,
                               detail="No protein-level data provided")

        if protein_p is not None and protein_p < 0.05:
            return LayerResult(2, "Cross-assay", LayerVerdict.SUPPORTIVE,
                               r=protein_r, p=protein_p,
                               detail=(f"{protein_label} r={protein_r:+.3f}, p={protein_p:.4f}; "
                                       "protein-level pathway support, not direct TF activity"))
        else:
            return LayerResult(2, "Cross-assay", LayerVerdict.FAIL,
                               r=protein_r, p=protein_p,
                               detail=f"Protein support not significant")

    def layer3_bulk(
        self, tf1: str, tf2: str, bulk_df: Optional[pd.DataFrame] = None,
    ) -> LayerResult:
        """L3: bulk-tissue cross-platform replication."""
        if bulk_df is None:
            return LayerResult(3, "Bulk replication", LayerVerdict.NA,
                               detail="No bulk data provided")

        if tf1 not in bulk_df.columns or tf2 not in bulk_df.columns:
            return LayerResult(3, "Bulk replication", LayerVerdict.NA,
                               detail=f"{tf1} or {tf2} not found in bulk data")

        r, p = pearsonr(bulk_df[tf1], bulk_df[tf2])
        if abs(r) > 0.15 and p < 0.05:
            return LayerResult(3, "Bulk replication", LayerVerdict.PASS, r=r, p=p,
                               detail=f"Bulk r={r:+.3f}, p={p:.4f}")
        else:
            return LayerResult(3, "Bulk replication", LayerVerdict.FAIL, r=r, p=p,
                               detail=f"Bulk r={r:+.3f}, p={p:.4f}; fails cross-platform consistency")

    def layer4_null(self, tf1: str, tf2: str) -> LayerResult:
        """L4: random TF-pair null specificity test."""
        if tf1 not in self.tf_scores.columns or tf2 not in self.tf_scores.columns:
            return LayerResult(4, "Random TF-pair null", LayerVerdict.NA,
                               detail="TF not in score matrix")

        r, p = pearsonr(self.tf_scores[tf1], self.tf_scores[tf2])
        abs_r = abs(r)
        pct = self._l4_percentile(abs_r)
        emp_p = self._l4_empirical_p(abs_r)

        if abs_r >= self._l4_p95:
            verdict = LayerVerdict.PASS
        elif abs_r >= self._l4_p75:
            verdict = LayerVerdict.SUPPORTIVE
        elif abs_r >= self._l4_median:
            verdict = LayerVerdict.FAIL
        else:
            verdict = LayerVerdict.FAIL

        n_tfs = len(self.tf_names)
        n_pairs = n_tfs * (n_tfs - 1) // 2
        detail = (f"|r|={abs_r:.3f}, pct={pct:.1f}%, emp_p={emp_p:.4f}; "
                  f"null: {n_tfs} TFs, {n_pairs} pairs, "
                  f"median={self._l4_median:.3f}, 75th={self._l4_p75:.3f}, "
                  f"95th={self._l4_p95:.3f}")

        return LayerResult(4, "Random TF-pair null", verdict, r=r, p=p, detail=detail)

    def layer5_composition(
        self, tf1: str, tf2: str, cell_type_props: Optional[pd.DataFrame] = None,
    ) -> LayerResult:
        """L5: composition-adjusted partial correlation."""
        if cell_type_props is None:
            return LayerResult(5, "Composition adjustment", LayerVerdict.NA,
                               detail="No cell-type proportion data provided")

        if tf1 not in self.tf_scores.columns or tf2 not in self.tf_scores.columns:
            return LayerResult(5, "Composition adjustment", LayerVerdict.NA,
                               detail="TF not in score matrix")

        common = self.tf_scores.index.intersection(cell_type_props.index)
        if len(common) < 10:
            return LayerResult(5, "Composition adjustment", LayerVerdict.NA,
                               detail=f"Insufficient overlapping samples: {len(common)}")

        from sklearn.linear_model import LinearRegression
        X = cell_type_props.loc[common].values
        y1 = self.tf_scores.loc[common, tf1].values
        y2 = self.tf_scores.loc[common, tf2].values

        resid1 = y1 - LinearRegression().fit(X, y1).predict(X)
        resid2 = y2 - LinearRegression().fit(X, y2).predict(X)
        r_partial, p_partial = pearsonr(resid1, resid2)

        if abs(r_partial) > 0.3 and p_partial < 0.05:
            verdict = LayerVerdict.PASS
        else:
            verdict = LayerVerdict.FAIL

        return LayerResult(5, "Composition adjustment", verdict,
                           r=r_partial, p=p_partial,
                           detail=f"Partial r={r_partial:+.3f}, p={p_partial:.4f}; n={len(common)}")

    def layer6_looco(
        self, tf1: str, tf2: str,
        cohort_labels: Optional[List[str]] = None,
    ) -> LayerResult:
        """L6: leave-one-cohort-out sensitivity."""
        if cohort_labels is None or len(set(cohort_labels)) < 3:
            return LayerResult(6, "Leave-one-cohort-out", LayerVerdict.NA,
                               detail=f"Need ≥3 cohorts; got {len(set(cohort_labels)) if cohort_labels else 0}")

        if tf1 not in self.tf_scores.columns or tf2 not in self.tf_scores.columns:
            return LayerResult(6, "Leave-one-cohort-out", LayerVerdict.NA,
                               detail="TF not in score matrix")

        cohorts = np.array(cohort_labels)
        all_r = []
        for cohort in sorted(set(cohorts)):
            mask = cohorts == cohort
            if mask.sum() >= 5:
                r, _ = pearsonr(
                    self.tf_scores[tf1].values[mask],
                    self.tf_scores[tf2].values[mask],
                )
                all_r.append(r)

        if len(all_r) < 3:
            return LayerResult(6, "Leave-one-cohort-out", LayerVerdict.NA,
                               detail=f"Fewer than 3 cohorts with ≥5 samples")

        # LOO meta-analysis
        loo_effects = []
        for i in range(len(all_r)):
            z_loo = np.mean([_fisher_z(all_r[j]) for j in range(len(all_r)) if j != i])
            loo_effects.append(z_loo)

        z_all = np.mean([_fisher_z(r) for r in all_r])
        z_loo_min = min(loo_effects)

        stable = all(z * z_all > 0 for z in loo_effects)

        if stable:
            verdict = LayerVerdict.PASS
            detail = (f"All LOO estimates same direction; "
                      f"pooled z={z_all:.3f}, min LOO z={z_loo_min:.3f}")
        else:
            verdict = LayerVerdict.FAIL
            detail = (f"Direction reversal in LOO; "
                      f"pooled z={z_all:.3f}, min LOO z={z_loo_min:.3f}")

        return LayerResult(6, "Leave-one-cohort-out", verdict,
                           detail=detail)

    def layer7_permutation(
        self, tf1: str, tf2: str, group_labels: Optional[np.ndarray] = None,
    ) -> LayerResult:
        """L7: disease-label permutation test."""
        if group_labels is None:
            return LayerResult(7, "Label permutation", LayerVerdict.NA,
                               detail="No group labels provided")

        if tf1 not in self.tf_scores.columns or tf2 not in self.tf_scores.columns:
            return LayerResult(7, "Label permutation", LayerVerdict.NA,
                               detail="TF not in score matrix")

        x = self.tf_scores[tf1].values
        y = self.tf_scores[tf2].values
        r_obs, _ = pearsonr(x, y)
        abs_obs = abs(r_obs)

        null_vals = []
        labels = np.array(group_labels)
        for _ in range(self.n_perms):
            perm = self.rng.permutation(labels)
            # Within each permuted group, compute r
            rs = []
            for grp in sorted(set(perm)):
                mask = perm == grp
                if mask.sum() >= 5:
                    rp, _ = pearsonr(x[mask], y[mask])
                    rs.append(abs(rp))
            if rs:
                null_vals.append(np.mean(rs))

        null_vals = np.array(null_vals)
        emp_p = (null_vals >= abs_obs).mean()
        p95 = np.percentile(null_vals, 95)

        if abs_obs >= p95:
            verdict = LayerVerdict.PASS
        else:
            verdict = LayerVerdict.FAIL

        return LayerResult(7, "Label permutation", verdict,
                           detail=(f"|r_obs|={abs_obs:.3f}, 95th null={p95:.3f}, "
                                   f"emp_p={emp_p:.4f}, n_perms={self.n_perms}"))

    # ================================================================
    # Full validation
    # ================================================================
    def validate(
        self,
        tf1: str,
        tf2: str,
        source: str = "",
        scores_method2: Optional[pd.DataFrame] = None,
        protein_r: Optional[float] = None,
        protein_p: Optional[float] = None,
        protein_label: str = "",
        bulk_df: Optional[pd.DataFrame] = None,
        cell_type_props: Optional[pd.DataFrame] = None,
        cohort_labels: Optional[List[str]] = None,
        group_labels: Optional[np.ndarray] = None,
    ) -> ValidationReport:
        """Run all seven SCEPTIC layers on a single TF-pair claim.

        Returns a ValidationReport with per-layer results and overall verdict.
        """
        claim = TFClaim(tf1=tf1, tf2=tf2, source=source)
        layers = []

        # L1
        layers.append(self.layer1_orthogonal(tf1, tf2, scores_method2))

        # L2
        layers.append(self.layer2_cross_assay(tf1, tf2, protein_r, protein_p, protein_label))

        # L3
        layers.append(self.layer3_bulk(tf1, tf2, bulk_df))

        # L4
        layers.append(self.layer4_null(tf1, tf2))

        # L5
        layers.append(self.layer5_composition(tf1, tf2, cell_type_props))

        # L6
        layers.append(self.layer6_looco(tf1, tf2, cohort_labels))

        # L7
        layers.append(self.layer7_permutation(tf1, tf2, group_labels))

        # Overall verdict
        l1_ok = layers[0].verdict == LayerVerdict.PASS
        l4_ok = layers[3].verdict in (LayerVerdict.PASS, LayerVerdict.SUPPORTIVE)
        l2_ok = layers[1].verdict != LayerVerdict.FAIL
        l3_ok = layers[2].verdict != LayerVerdict.FAIL

        if l1_ok and l4_ok:
            overall = SCEPTICVerdict.SUPPORTED
            summary = (f"{tf1}-{tf2}: passes L1 and L4; "
                       "supported in this context")
        elif not l1_ok and not l4_ok:
            overall = SCEPTICVerdict.REJECTED
            summary = (f"{tf1}-{tf2}: fails both L1 and L4; "
                       "rejected as non-specific or method-dependent")
        elif l1_ok and not l4_ok:
            if l2_ok or l3_ok:
                overall = SCEPTICVerdict.CONTEXT_DEPENDENT
                summary = (f"{tf1}-{tf2}: passes L1, fails L4; "
                           "context-dependent covariation only")
            else:
                overall = SCEPTICVerdict.INSUFFICIENT
                summary = (f"{tf1}-{tf2}: passes L1, fails L4; "
                           "no cross-context support")
        else:
            overall = SCEPTICVerdict.REJECTED
            summary = f"{tf1}-{tf2}: fails L1; method-dependent signal"

        return ValidationReport(
            claim=claim,
            layers=layers,
            overall_verdict=overall,
            summary=summary,
        )


def validate_claim(
    tf_scores: pd.DataFrame,
    tf1: str,
    tf2: str,
    **kwargs,
) -> ValidationReport:
    """Convenience function for single-claim validation."""
    validator = SCEPTICValidator(tf_scores)
    return validator.validate(tf1, tf2, **kwargs)
