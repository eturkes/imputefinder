# imputefinder 0.99.5 (2026-07-18)

- Adds the validated mandatory sidecar design core: canonical declared model
  matrices, SVD rank and null-space aliases, row-space contrast estimability,
  and block-aware independent-unit accounting.

# imputefinder 0.99.4 (2026-07-18)

- Adds experimental `analyze_missingness()` input-first sidecars with exact
  classic-fit compatibility checks, packed original masks, pre-rescue evidence,
  typed module routing, portable unavailable results, and compact printing.

# imputefinder 0.99.3 (2026-07-18)

- Adds the experimental `missingness_design()` constructor for explicit
  condition, nuisance, block, acquisition, and interaction roles, with strict
  sample alignment and versioned schema validation.

# imputefinder 0.99.2 (2026-07-15)

- Completes release-candidate audits of public claims, edge cases, process
  side effects, named-order invariance, and seed-only output mutation.
- Excludes local check artefacts from source packages and passes the full local
  `R CMD check --as-cran` and repository-controlled `BiocCheck` gates.
- Validates the manual and automatic paths on a deterministic 10,000-feature,
  50-sample workload with predeclared runtime and allocation limits.

# imputefinder 0.99.1 (2026-07-15)

- Introduces condition-aware classification of likely MAR and MNAR
  protein-intensity missingness for matrices and `SummarizedExperiment`
  objects.
- Preserves plausible condition-specific on/off proteins with a reproducible,
  audited one-cell rescue while filtering globally absent features, features
  with a block lacking enough MAR evidence, and all-MNAR features.
- Provides manual and automatic condition-specific cutoffs, stored missingness
  profiles, diagnostic plots, seed provenance, and feature-level decision
  records for downstream normalisation and imputation.
- Pins rescue sampling to a restored local R RNG configuration so a fixed seed
  is independent of caller RNG settings.
- Records numerical warnings from the automatic trend model and rejects the
  affected cutoff instead of silently accepting a suspect fit.
