"""
SCEPTIC CLI: command-line interface for the 7-layer validation framework.
"""
import sys
import os
import json
import argparse
import numpy as np
import pandas as pd

from sceptic.core import TRRUST_REGULONS
from sceptic.validator import SCEPTICValidator, validate_claim
from sceptic.io import load_expression_matrix, load_metadata, load_bulk_scores


def cmd_validate(args):
    """Run SCEPTIC validation on a TF pair."""
    # Load data
    print(f"Loading expression matrix: {args.expr}", file=sys.stderr)
    expr = load_expression_matrix(args.expr, compression=args.compression)

    print(f"Loading metadata: {args.meta}", file=sys.stderr)
    meta = load_metadata(args.meta, sample_col=args.sample_col)

    print("Computing per-sample TF scores...", file=sys.stderr)
    from sceptic.io import compute_tf_scores
    tf_scores = compute_tf_scores(expr, meta, args.sample_col,
                                   min_targets=args.min_targets)

    # Optional orthogonal scoring
    scores_method2 = None
    if args.orthogonal:
        print(f"Loading orthogonal scores: {args.orthogonal}", file=sys.stderr)
        scores_method2 = pd.read_csv(args.orthogonal, index_col=0)

    # Optional bulk data
    bulk_df = None
    if args.bulk:
        bulk_df = load_bulk_scores(args.bulk)

    # Optional cell-type proportions
    cell_type_props = None
    if args.cell_types:
        cell_type_props = pd.read_csv(args.cell_types, index_col=0)

    # Cohort labels
    cohort_labels = None
    if args.cohorts:
        with open(args.cohorts) as f:
            cohort_labels = [line.strip() for line in f]

    # Group labels
    group_labels = None
    if args.groups:
        group_labels = np.array(meta[args.group_col].values)

    # Run validation
    validator = SCEPTICValidator(
        tf_scores,
        n_perms=args.n_perms,
        random_seed=args.seed,
    )

    report = validator.validate(
        tf1=args.tf1,
        tf2=args.tf2,
        source=args.source or "",
        scores_method2=scores_method2,
        protein_r=args.protein_r,
        protein_p=args.protein_p,
        protein_label=args.protein_label or "",
        bulk_df=bulk_df,
        cell_type_props=cell_type_props,
        cohort_labels=cohort_labels,
        group_labels=group_labels,
    )

    # Output
    _print_report(report, args.format)


def cmd_audit(args):
    """Run SCEPTIC L4 literature audit."""
    print(f"Loading TF scores: {args.scores}", file=sys.stderr)
    tf_scores = pd.read_csv(args.scores, index_col=0)
    # Remove non-TF columns
    tf_cols = [c for c in tf_scores.columns if c != "disease"]
    tf_scores = tf_scores[tf_cols]

    # Load test pairs
    pairs = []
    if args.pairs:
        pairs_df = pd.read_csv(args.pairs)
        for _, row in pairs_df.iterrows():
            pairs.append((
                str(row["TF1"]), str(row["TF2"]), str(row.get("Source", ""))
            ))
    elif args.tf1 and args.tf2:
        pairs.append((args.tf1, args.tf2, args.source or ""))

    validator = SCEPTICValidator(tf_scores, n_perms=args.n_perms, random_seed=args.seed)

    print(f"\n{'='*80}")
    print(f"SCEPTIC L4 Literature Audit")
    print(f"N samples = {len(tf_scores)}, N TFs = {len(tf_scores.columns)}")
    print(f"{'='*80}\n")

    results = []
    for tf1, tf2, source in pairs:
        r = validator.layer4_null(tf1, tf2)
        results.append({
            "TF1": tf1, "TF2": tf2, "r": r.r, "|r|": abs(r.r) if r.r else 0,
            "p": r.p, "Percentile": validator._l4_percentile(abs(r.r)) if r.r else 0,
            "Empirical_p": validator._l4_empirical_p(abs(r.r)) if r.r else 0,
            "L4_verdict": r.verdict.value,
            "Source": source,
        })

        verdict_label = {
            "pass": "PASS (top 5%)",
            "supportive": "MARGINAL (top 25%)",
            "fail": "Not PBMC-specific",
            "na": "N/A",
        }[r.verdict.value]

        print(f"  {tf1:<10} {tf2:<10} |r|={abs(r.r):.3f} "
              f"pct={validator._l4_percentile(abs(r.r)):.1f}%  {verdict_label}  {source}")

    if args.output:
        pd.DataFrame(results).to_csv(args.output, index=False)
        print(f"\nResults saved to {args.output}")

    # Summary
    res_df = pd.DataFrame(results)
    n_pairs = len(res_df)
    n_pass = sum(res_df["L4_verdict"] == "pass")
    n_marg = sum(res_df["L4_verdict"] == "supportive")
    n_fail = sum(res_df["L4_verdict"] == "fail")
    print(f"\nSummary: {n_pass}/{n_pairs} PASS, {n_marg}/{n_pairs} MARGINAL, "
          f"{n_fail}/{n_pairs} not PBMC-specific")


def _print_report(report, fmt="text"):
    """Print validation report."""
    if fmt == "json":
        out = {
            "claim": {"tf1": report.claim.tf1, "tf2": report.claim.tf2,
                       "source": report.claim.source},
            "overall_verdict": report.overall_verdict.value,
            "summary": report.summary,
            "layers": [
                {
                    "layer": l.layer,
                    "name": l.name,
                    "verdict": l.verdict.value,
                    "r": l.r,
                    "p": l.p,
                    "detail": l.detail,
                }
                for l in report.layers
            ],
        }
        print(json.dumps(out, indent=2))
        return

    # Text format
    print(f"\n{'='*80}")
    print(f"SCEPTIC Validation Report")
    print(f"{'='*80}")
    print(f"Claim:    {report.claim.tf1}–{report.claim.tf2}")
    if report.claim.source:
        print(f"Source:   {report.claim.source}")
    print(f"\n{'Layer':<6} {'Name':<30} {'Verdict':<12} {'Detail'}")
    print(f"{'-'*6} {'-'*30} {'-'*12} {'-'*40}")
    for l in report.layers:
        v_label = {
            "pass": "✓ PASS",
            "fail": "✗ FAIL",
            "supportive": "~ SUPPORTIVE",
            "na": "— N/A",
        }[l.verdict.value]
        detail = l.detail[:60] if l.detail else ""
        print(f"L{l.layer:<5} {l.name:<30} {v_label:<12} {detail}")

    print(f"\n{'─'*80}")
    verdict_label = {
        "supported": "SUPPORTED",
        "context_dependent": "CONTEXT-DEPENDENT",
        "rejected": "REJECTED",
        "insufficient_evidence": "INSUFFICIENT EVIDENCE",
    }[report.overall_verdict.value]
    print(f"SCEPTIC VERDICT: {verdict_label}")
    print(f"Summary: {report.summary}")


def main():
    parser = argparse.ArgumentParser(
        description="SCEPTIC: Single-Cell Expression Program Testing and Integrity Checklist",
    )
    parser.add_argument("--version", action="version", version="sceptic 1.0.0")

    sub = parser.add_subparsers(dest="command")

    # ---- validate ----
    p_val = sub.add_parser("validate", help="Validate a TF co-expression claim")
    p_val.add_argument("--expr", required=True, help="Expression matrix (CSV or gzipped CSV)")
    p_val.add_argument("--meta", required=True, help="Cell metadata (TSV)")
    p_val.add_argument("--tf1", required=True, help="First TF")
    p_val.add_argument("--tf2", required=True, help="Second TF")
    p_val.add_argument("--source", help="Source description")
    p_val.add_argument("--sample-col", default="donor_id",
                       help="Metadata column for sample/donor ID")
    p_val.add_argument("--compression", default="gzip",
                       help="Compression (gzip/none)")
    p_val.add_argument("--min-targets", type=int, default=3,
                       help="Minimum expressed targets per TF")
    p_val.add_argument("--orthogonal", help="Orthogonal scoring matrix (CSV)")
    p_val.add_argument("--bulk", help="Bulk tissue TF scores (CSV)")
    p_val.add_argument("--cell-types", help="Cell-type proportion matrix (CSV)")
    p_val.add_argument("--cohorts", help="Cohort labels file (one per sample)")
    p_val.add_argument("--groups", help="Group label column in metadata")
    p_val.add_argument("--group-col", default="disease",
                       help="Group label column name")
    p_val.add_argument("--protein-r", type=float, help="Protein-level r")
    p_val.add_argument("--protein-p", type=float, help="Protein-level p")
    p_val.add_argument("--protein-label", help="Protein assay label")
    p_val.add_argument("--n-perms", type=int, default=10000,
                       help="Number of permutations")
    p_val.add_argument("--seed", type=int, default=42, help="Random seed")
    p_val.add_argument("--format", default="text", choices=["text", "json"],
                       help="Output format")
    p_val.set_defaults(func=cmd_validate)

    # ---- audit ----
    p_audit = sub.add_parser("audit", help="Run L4 literature audit")
    p_audit.add_argument("--scores", required=True, help="TF score matrix (CSV)")
    p_audit.add_argument("--pairs", help="CSV with TF1,TF2,Source columns")
    p_audit.add_argument("--tf1", help="Single pair TF1")
    p_audit.add_argument("--tf2", help="Single pair TF2")
    p_audit.add_argument("--source", help="Source description")
    p_audit.add_argument("--output", help="Output CSV path")
    p_audit.add_argument("--n-perms", type=int, default=10000)
    p_audit.add_argument("--seed", type=int, default=42)
    p_audit.set_defaults(func=cmd_audit)

    args = parser.parse_args()
    if args.command is None:
        parser.print_help()
        sys.exit(1)

    args.func(args)


if __name__ == "__main__":
    main()
