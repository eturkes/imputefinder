# ImputeFinder

ImputeFinder classifies missing protein intensities separately within each
experimental condition. This preserves condition-specific on/off patterns that
can disappear when all samples are pooled. The result contains the filtered
data, the classification of each feature-condition block, and the reason each
feature was retained or dropped.

The package expects **unnormalised log2 protein intensities**. It classifies and
prepares data; it does not normalise or impute them.

Use it when you need separate downstream treatment for low-intensity missing
blocks and better-observed blocks, without pooling experimental conditions.
The labels mean **likely missing-not-at-random (MNAR)** and **likely
missing-at-random (MAR)** under the cutoff rule. They do not establish why a
measurement is missing.

## Installation

ImputeFinder is under development and is not yet available from a Bioconductor
package repository. The current development package requires R 4.6 or later
and a compatible Bioconductor library. Install the GitHub version in R:

```r
install.packages(c("BiocManager", "remotes"))
BiocManager::install("eturkes/imputefinder")
```

## First result

Run this from a writable working folder after installation. It uses seven
synthetic features in two conditions, downloads no data, and saves the figure
shown below. The manual cutoffs are illustrative, not recommendations for your
experiment.

```r
source(system.file("examples", "quickstart.R", package = "imputefinder",
    mustWork = TRUE))
```

```text
<imputefinder_result>
Retained 4 of 7 features; dropped 3.
Analysed 8 samples in 2 conditions; inserted 1 rescue value.
A: retained features include 2 MNAR, 1 MAR, 1 complete, 0 absent; cutoff 12 (manual)
B: retained features include 0 MNAR, 1 MAR, 3 complete, 0 absent; cutoff 12 (manual)
```

![Classification by condition and retention for seven synthetic features. Four are kept; three are dropped. The on_off feature is MNAR after rescue in A and complete in B.](man/figures/quickstart.png)

Read each row across conditions A and B, then inspect its retention decision.
The `on_off` feature stays because B is complete. One logged artificial value
allows classification in A; it does not complete an imputation.

The sparse, high-intensity `sparse_mar` block lacks an observed majority in A,
so the whole feature is dropped. Low intensity alone does not cause exclusion:
`complete_low` stays because nothing is missing.

The script writes these files under `imputefinder-results/`:

| File | Contents |
|---|---|
| `classifications.png` | The labelled state and retention figure, at 180 dpi |
| `classifications.csv` | Original counts, processed means, cutoffs, states, and rescue flags |
| `feature-status.csv` | Every original feature's retention decision and drop reason |
| `retained-intensities.csv` | The four retained rows, including the rescue value |
| `seed-log.csv` | Every inserted value and its location |
| `result.rds` | The full result, for `readRDS()` |

Re-running replaces these example files. The script leaves `x`, `condition`,
`fit`, and `state_plot` in your R session. Read the
[example script](inst/examples/quickstart.R) to change the inputs or figure.

## Classify your data

Replace `x` with your unnormalised log2 intensity matrix. Give its rows feature
identifiers and its columns unique sample names. Align `condition` with the
matrix columns. Use `NA` for missing values; zero remains an observed value.

Choose cutoffs for your experiment rather than copying 12 from this example.
The same public call runs the demonstration:

```r
library(imputefinder)

fit <- classify_missingness(
    x, group = condition, cutoffs = c(A = 12, B = 12),
    rescue = "minimum", seed = 1L
)
fit$feature_status
fit$seed_log
plot_missingness(fit, "A")
```

The diagnostic plot shows the missingness profile and the recorded cutoff.
Its vertical line is the manual value here, not an estimated boundary.
For larger datasets, `cutoffs = NULL` requests automatic selection; unsupported
profiles require a manual decision.

Use `rescue = "none"` to keep unobserved blocks `absent` without inserting values.
Keep those blocks separate from the MAR and MNAR branches. The result still
needs downstream normalisation and imputation choices.

<details>
<summary>Classification rules, rescue, cutoffs, and downstream workflow</summary>

## Classification rules

Each feature-condition block receives one state after the selected rescue policy:

| State | Rule | Downstream group |
|---|---|---|
| `complete` | No values are missing | `MAR_or_complete` |
| `MNAR` | Values are missing and the observed arithmetic mean is below the condition cutoff | `MNAR` |
| `MAR` | Values are missing, the mean is at or above the cutoff, and strictly more than half of the samples are observed | `MAR_or_complete` |
| `insufficient` | Values are missing, the mean is at or above the cutoff, and no strict observed majority exists | Feature is dropped |
| `absent` | No original values are observed and `rescue = "none"` | `absent`; no mechanism is assigned |

A mean equal to the cutoff is on the MAR side. The majority rule applies only
to MAR candidates. A feature is retained only when no condition is
`insufficient` and at least one condition is `complete` or `MAR`.
Globally missing features and features that are MNAR or absent everywhere are
dropped.

## How rescue works

An on/off feature has no observed value from which to calculate its mean in the
absent condition. The default `rescue = "minimum"` handles it before classification:

1. Features missing in every sample are marked `all_missing` and removed.
2. The minimum finite intensity in each condition is calculated once, before
   any artificial values are inserted.
3. For each surviving feature fully missing in a condition, one sample from
   that condition is selected uniformly and receives the condition minimum.
   The remaining cells stay `NA`.
4. Every inserted value is recorded in `seed_log`.

For the worked fixture:

```r
fit$seed_log
#>   feature condition sample old_value inserted_value seed
#> 1  on_off         A     s1        NA              8    1
```

The rescue value is a low-intensity anchor used to enable classification; it is
not a completed imputation. With a fixed dataset and `seed`, the selected sample
does not change when rows or columns are reordered. The default is `seed = 1L`.
ImputeFinder leaves the global random-number state unchanged, does not modify
the input, and preserves all originally observed values. The inserted value
remains in `fit$data` and can affect later normalization and inference.

### Classify without synthetic values

Use `rescue = "none"` to preserve missing cells and record absence separately:

```r
unseeded <- classify_missingness(
    x,
    group = condition,
    cutoffs = c(A = 12, B = 12),
    rescue = "none"
)
unseeded$groups$A$absent
#> [1] "on_off"
unseeded$classifications[
    unseeded$classifications$feature == "on_off",
    c("condition", "state", "original_observed_count",
      "original_missing_count", "original_mean_intensity")
]
```

The absent block has zero observations, an undefined mean, and no cutoff.
It does not enter density estimation or automatic-cutoff evidence.
An entirely unobserved condition receives an empty, labelled diagnostic plot.
Globally missing features are still dropped with `all_missing` records.
Features that are MNAR or absent everywhere receive `MNAR_or_absent_all_conditions`
when at least one block is absent.

Absence is not evidence for either a MAR or MNAR mechanism. Keep the `absent`
group separate until the analysis specifies how to handle unobserved blocks.
This mode prevents synthetic values; it does not establish biological absence.

## Manual and automatic cutoffs

Use manual cutoffs when the boundary is known in advance or fixed by an analysis
plan. Supply a named, finite value for each relevant condition. A partial vector
is allowed; automatic selection fills the omitted conditions. A cutoff outside
the observed feature-mean range is accepted but recorded with a warning.

With `cutoffs = NULL`, each condition containing partially observed blocks is
assessed independently for an automatic cutoff. Complete or fully absent blocks
need no cutoff:

```r
fit_auto <- classify_missingness(
    intensity_matrix,
    group = sample_condition
)
```

For each condition, automatic selection estimates how the proportion of
incomplete features changes with mean intensity. At intensity `x`, the plotted
curve is:

```text
n_missing f_missing(x)
-----------------------------------------------
n_missing f_missing(x) + n_complete f_complete(x)
```

The detector uses the lower, right-hand boundary of the lowest-intensity
credible downward transition. The method version, supporting values, and
warnings are stored in `cutoff_diagnostics`. Multiple credible transitions now
signal a warning by default, as well as appearing in the plot caption.
Use `ambiguity = "error"` to require a manual decision in such cases.
Use `ambiguity = "lowest"` to accept the lowest transition explicitly without
signalling a warning. Stored warnings remain available under every policy.
The `quality$candidate_transitions` table records candidate cutoffs and support.
These candidates are not confidence intervals or calibrated probabilities.
This selection policy does not establish a unique biological boundary.

Strict failures inherit from `imputefinder_cutoff_ambiguous` and
`imputefinder_cutoff_error`. They retain `profile`, `diagnostic`, and `condition`
fields. A named manual cutoff bypasses automatic ambiguity for that condition.

When the profile is too small, flat, or unsupported, the function reports why
it could not estimate a cutoff and asks for a manual value.
Small examples such as the one above should use manual cutoffs.

Plotting uses the stored profile and cutoff:

```r
plot_missingness(fit, "A")
fit$cutoff_diagnostics$A
```

Bottom rug ticks identify rescued feature-condition blocks. The vertical line
is the recorded cutoff.

## What the result contains

`classify_missingness()` returns an `imputefinder_result` with these components:

| Component | Contents |
|---|---|
| `data` | Filtered data with any recorded rescue values; the input representation and metadata are preserved |
| `classifications` | Processed statistics, original observed and missing counts, original means, and final states for each surviving feature-condition block |
| `groups` | Retained `MNAR`, `MAR`, `complete`, `MAR_or_complete`, and `absent` feature names for every condition |
| `feature_status` | One row per original feature with its retention decision and drop reason |
| `cutoffs` | Named condition-specific numeric cutoffs |
| `cutoff_diagnostics` | Whether each cutoff was manual, automatic, or unnecessary, together with method details and warnings |
| `profiles` | Feature statistics and the density grid used for each condition |
| `seed_log` | The feature, condition, selected sample, inserted value, and seed for every rescue, including rescues on features later dropped |
| `groups_by_sample` | Condition vector named and ordered by output samples |
| `call` | The matched call used to create the result |
| `policies` | Resolved rescue and ambiguity policies |
| `assay_index` | Position of the classified experiment assay, or `NULL` for matrix input |

`summary(fit)` reports all state and drop counts. `plot_missingness()` displays
the profile stored in the result, so it does not need the original input.
`summary(fit)$seeded_complete` counts blocks that became complete solely through
rescue. Original provenance uses `original_observed_count`,
`original_missing_count`, and `original_mean_intensity` in `classifications`.

For a `SummarizedExperiment`, select the condition column explicitly and select
the assay whenever more than one is present:

```r
fit_se <- classify_missingness(
    se,
    group_col = "condition",
    assay = "log2_intensity",
    cutoffs = c(A = 12, B = 12)
)
```

## Downstream normalisation and imputation

The intended pipeline is:

1. Run ImputeFinder on unnormalised log2 intensities.
2. Use `fit$data`, which excludes unsupported features and includes rescue
   values on retained features. `seed_log` also records rescues on features that
   were later dropped.
3. Apply the chosen normalisation to that complete returned object.
4. Split columns by `fit$groups_by_sample`.
5. Within each condition, apply an MNAR method only to rows named in
   `fit$groups[[condition]]$MNAR` and a MAR method only to rows named in
   `fit$groups[[condition]]$MAR_or_complete`.
6. Recombine conditions in the original sample order.

With `rescue = "none"`, handle `fit$groups[[condition]]$absent` separately.
The package does not route unobserved blocks to either imputation branch.
The returned object can still contain missing values under either policy.

For example, these indices isolate the two branches without imposing an
imputation package:

```r
condition_name <- "A"
sample_names <- names(fit$groups_by_sample)[
    fit$groups_by_sample == condition_name
]
condition_data <- fit$data[, sample_names, drop = FALSE]

mnar_rows <- fit$groups[[condition_name]]$MNAR
mar_rows <- fit$groups[[condition_name]]$MAR_or_complete

condition_data[mnar_rows, , drop = FALSE]
condition_data[mar_rows, , drop = FALSE]
```

### Select only missing cells

Use `missingness_mask()` to identify replacement targets in the retained data:

```r
mnar_mask <- missingness_mask(fit, "MNAR")
mar_mask <- missingness_mask(fit, "MAR")
absent_mask <- missingness_mask(fit, "absent")
stopifnot(identical(mnar_mask | mar_mask | absent_mask, is.na(fit$data)))
```

Each mask has the same dimensions, names, and order as the classified data.
Only remaining missing cells are `TRUE`. Observed and rescue values are
`FALSE`; dropped features have no rows. Keep absent blocks separate until
the analysis specifies how to handle them.

Create masks before normalisation or imputation. Align downstream replacement
matrices to the same feature and sample names before assigning selected cells.
Masks select targets, not the complete training data an imputation method needs.

For experiment results, the helper uses the recorded classified assay without
repeating its name. Keep the result's assay order unchanged. Older multi-assay
results without `assay_index` need reclassification before using this helper.

Normalisation and actual kNN, MinProb, MinDet, QRILC, BPCA, or other imputation
remain downstream choices, not ImputeFinder defaults or runtime dependencies.

</details>

## Further reading

The [vignette](vignettes/imputefinder.Rmd) covers automatic cutoffs, sensitivity
checks, downstream branching, and a bundled PXD020717 Tau-phosphosite extract.
The selected real-data example demonstrates classification, not disease effects.

For a rendered local vignette, install with `build_vignettes = TRUE` and its
suggested dependencies, then run `vignette("imputefinder")`. The first-run script
does not need a vignette build. Run `citation("imputefinder")` for citation details.

Regenerate the README figure from the example with
`Rscript tools/render-quickstart.R`. The package tests run the same example and
check its classifications and exported files.

## Verification

Deterministic regression tests compare package output with known classifications
under both manual and automatic cutoffs. They cover intensity-independent MAR,
intensity-dependent MNAR, condition-specific on/off features, unequal group
sizes, rescue reproducibility, and input-order invariance. The automatic-cutoff
tests also include profiles for which no defensible boundary can be estimated.

These simulations verify implementation behavior under their stated model; they
are not estimates of accuracy for arbitrary proteomics experiments.

Run the committed build and check from the repository root:

```sh
Rscript tools/check-package.R
Rscript tools/check-package.R --bioc
```

Install the package dependencies first. The second command also runs BiocCheck.
Logs remain under `.checks/`. The checks execute examples, tests, and vignettes.
CI is configured for Linux, macOS, Windows, and Linux with R-devel.
Hosted results are established only when those jobs actually run.

## Limitations

- ImputeFinder cannot determine whether values are log2 transformed or
  unnormalised. It checks their structure and finiteness, then relies on the
  user to supply suitable data.
- The intensity cutoff is a model assumption. Biology, acquisition, and
  preprocessing can produce missingness patterns that do not follow one
  descending transition.
- Arithmetic mean intensity is intentional: one strong observed signal can move
  a partial block toward MAR. It can also be sensitive to unusual observations.
- Automatic selection needs both incomplete and complete features. When the
  data do not support a clear cutoff, supply a manual value.
- A rescue value is artificial and remains visible in the returned data and
  log. Record its role and assess its downstream effect.
- With one sample in a condition, rescue fills the entire block, which becomes
  `complete`. Inspect original counts, `seeded`, and `seed_log`, or select
  `rescue = "none"` to keep the block absent. Replicated conditions are needed to assess
  within-condition missingness.
- The package classifies and coordinates a workflow. It neither proves causal
  missingness mechanisms nor chooses, runs, or validates an imputation method.

## License

ImputeFinder is distributed under GPL-3.0-or-later. See [COPYING](COPYING).
