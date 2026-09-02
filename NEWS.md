# imputefinder 0.99.4 (2026-09-01)

## Documentation and examples

- Adds a worked comparison showing why on/off features need separate condition
  states.
- Documents the synthetic regression results and explains their limitations.
- Adds a benchmark-style cutoff sweep with known MAR and MNAR states to show the
  classifier's state-recall trade-off.
- Adds a package citation and removes a stale reference to deleted development
  validation tooling.

# imputefinder 0.99.3 (2026-09-01)

- Records the maintainer ORCID and regenerates package documentation with the
  current roxygen2 release.
- Removes development-only scientific, cutoff, and performance tooling from
  the default submission branch so it contains package source only.
- Repairs repository documentation links.

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
  recorded one-cell rescue while filtering globally absent features, features
  with a block lacking enough MAR evidence, and all-MNAR features.
- Provides manual and automatic condition-specific cutoffs, stored missingness
  profiles, diagnostic plots, a rescue log, and feature-level decisions for
  downstream normalisation and imputation.
- Pins rescue sampling to a restored local R RNG configuration so a fixed seed
  is independent of global RNG settings.
- Records numerical warnings from the automatic trend model and rejects the
  affected cutoff instead of silently accepting a suspect fit.
