# imputefinder 0.99.6

- Adds `rescue = "none"` to retain on/off patterns without inserting values.
  Fully missing blocks receive an `absent` state and a separate downstream
  group. They do not contribute to intensity fitting or need a cutoff.
- Handles unobserved conditions and globally empty data in no-rescue mode,
  including labelled empty plots and complete feature-drop records.
- Adds original observed counts, missing counts, and means to classifications.
  Result printing flags complete blocks supported only by synthetic rescue.
- Adds `ambiguity = "warn"` as the default for competing automatic transitions.
  Use `"error"` to require a manual cutoff, or `"lowest"` to accept the existing
  selection explicitly. Diagnostics retain all candidate transitions.
- Advances automatic diagnostics to version 3 without changing cutoff values
  under the existing lowest-transition policy.
- Keeps minimum rescue as the default. State, group, result, and summary schemas
  gain additive fields; code that assumes exact field lists must be updated.
- Adds a reproducible package-check command and cross-platform CI configuration.

# imputefinder 0.99.5

- Records a diagnostic warning when multiple credible downward transitions
  exist. The plot caption shows the warning. The automatic method version is
  now 2; the existing lowest-intensity selection rule and cutoff calculations
  are unchanged.
- Uses a fixed sample order for feature means so numerical cancellation cannot
  change classifications when input columns are reordered.
- Documents that rescue in a one-sample condition produces a `complete` block,
  with its artificial value recorded in `seeded` and `seed_log`.

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
