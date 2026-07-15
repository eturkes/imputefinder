# imputefinder 0.99.1

- Introduces condition-aware classification of likely MAR and MNAR
  protein-intensity missingness for matrices and `SummarizedExperiment`
  objects.
- Preserves plausible condition-specific on/off proteins with a reproducible,
  audited one-cell rescue while filtering globally absent features, blocks
  without enough MAR evidence, and all-MNAR features.
- Provides manual and automatic condition-specific cutoffs, stored missingness
  profiles, diagnostic plots, seed provenance, and feature-level decision
  records for downstream normalisation and imputation.
