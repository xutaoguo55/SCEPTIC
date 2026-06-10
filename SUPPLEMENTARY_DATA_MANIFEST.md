# Supplementary Data Manifest

Processed manuscript data are distributed with the submission package rather
than committed to the public code repository.

Expected archive name:

`SCEPTIC_Supplementary_Data.zip`

Required files:

- `processed_scores/mas_sepsis_ibd_per_sample.csv`
- `processed_scores/scp548_per_sample_tf_means.csv`
- `processed_scores/emtab10026_per_sample_tf_means.csv`
- `processed_scores/gse65682_full_merged.csv`
- `literature_audit/literature_audit_scp548.csv`
- `literature_audit/all_tf_pairs_scp548.csv`
- `literature_audit/scp548_all_tf_scores.csv`
- `literature_audit/hctl_3cohort_crossvalidation.csv`
- `supplementary_tables/table_s3_qc_metrics.csv`
- `supplementary_tables/table_s4_matched_null_sensitivity.csv`
- `supplementary_tables/table_s5_spearman_robustness.csv`
- `supplementary_tables/table_s6_statistical_testing_families.csv`
- `supplementary_tables/table_s7_sceptic_parameters_thresholds.csv`
- `supplementary_tables/machine_readable_sceptic_checklist.csv`
- `supplementary_tables/machine_readable_sceptic_checklist.json`

Notes:

- `literature_audit/literature_audit_scp548.csv` contains 18 rows: the
  15 published TF-pair claims reported in Supplementary Table S2 plus the
  3 HCTL case-study pairs used for internal comparison in the same SCP548
  PBMC null distribution.
- `supplementary_tables/table_s5_spearman_robustness.csv` includes Fisher-z
  95% confidence intervals for the primary Pearson correlations.
- `supplementary_tables/table_s6_statistical_testing_families.csv` defines
  testing families, p-value types, confidence-interval reporting, and
  interpretation limits for each SCEPTIC layer.
- `supplementary_tables/table_s7_sceptic_parameters_thresholds.csv` documents
  thresholds as reporting heuristics rather than universal biological cutoffs.
- `supplementary_tables/machine_readable_sceptic_checklist.csv` and `.json`
  provide reviewer-facing checklist outputs in machine-readable formats.

Raw public data are available from SCP548, E-MTAB-10026, GSE65682, and the
Pienkos et al. OFID data DOI `10.5281/zenodo.15199327`. In-house discovery data are available from the corresponding
author under institutional data access approval.
