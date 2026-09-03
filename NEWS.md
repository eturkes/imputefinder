# imputefinder 0.99.4 (2026-09-01)

## Documentation and examples

- Adds a worked comparison showing why on/off features need separate condition
  states.
- Adds a cutoff-sensitivity example with known MAR and MNAR states and explains
  the limits of synthetic validation.
- Adds a classification workflow, a cell-level replacement audit, and annotated
  cutoff confusion heatmaps with explicit pre-filter denominators.
- Adds a package citation.

# imputefinder 0.99.3 (2026-09-01)

- Records the maintainer ORCID.
- Repairs repository documentation links.

# imputefinder 0.99.2 (2026-07-15)

- Stabilizes results across input order and caller RNG state.
- Ensures rescue insertions are the only changes made to retained matrix values.
- Tightens manual and automatic cutoff handling for unsupported profiles and
  numerical warnings.

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
