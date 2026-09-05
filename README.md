# ImputeFinder

ImputeFinder classifies missing protein intensities separately within each
experimental condition. This preserves condition-specific on/off patterns that
can disappear when all samples are pooled. The result contains the filtered
data, the classification of each feature-condition block, and the reason each
feature was retained or dropped.

The package expects **unnormalised log2 protein intensities**. It classifies and
prepares data; it does not normalise or impute them.

## Why condition-specific classification?

Missing protein intensities do not all follow the same pattern. Scattered
values may be absent for approximately intensity-independent, MAR-like reasons,
while low-abundance proteins may fall below detection in an MNAR-like way. A
protein can also be present in one condition and absent in another. Pooling all
samples before classification can hide that on/off pattern.

ImputeFinder therefore assigns one state to each feature-condition block:

| State | Rule after the selected rescue policy | Downstream group |
|---|---|---|
| `complete` | No values are missing | `MAR_or_complete` |
| `MNAR` | Values are missing and the observed arithmetic mean is below the condition cutoff | `MNAR` |
| `MAR` | Values are missing, the mean is at or above the cutoff, and strictly more than half of the samples are observed | `MAR_or_complete` |
| `insufficient` | Values are missing, the mean is at or above the cutoff, and no strict observed majority exists | Feature is dropped |
| `absent` | No original values are observed and `rescue = "none"` | `absent`; no mechanism is assigned |

A mean equal to the cutoff is on the MAR side. The majority rule applies only
to MAR candidates: one low-intensity observation can still support a sparse
MNAR block. A feature is retained only when no condition is `insufficient` and
at least one condition is `complete` or `MAR`. Features absent globally and
features classified as `MNAR` or `absent` in every condition are dropped.

The labels mean **likely MNAR** and **likely MAR** given the observed data. They
are practical classifications, not proof of the underlying missingness
mechanism.

## Installation

ImputeFinder is under development and is not yet available from a Bioconductor
package repository. The current development package requires R 4.6 or later
and a compatible Bioconductor library. Install the GitHub version with
[`BiocManager`](https://bioconductor.org/install/):

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("eturkes/imputefinder")
```

After the package enters a Bioconductor release, the repository installation
will be:

```r
BiocManager::install("imputefinder")
```

After installation, run `vignette("imputefinder")` for an executable walkthrough
of cutoff selection, on/off rescue, diagnostics, and downstream branching. The
[vignette source](vignettes/imputefinder.Rmd) is also available in this
repository. Run `citation("imputefinder")` to cite the package and the thesis
that describes the method.

## Quick start

This example includes an on/off protein, MAR-like and MNAR-like blocks, an
unsupported sparse block, and a globally absent feature. Conditions A and B
each contain four samples, and both manual cutoffs are 12.

```r
library(imputefinder)

x <- rbind(
    on_off = c(NA, NA, NA, NA, 15, 16, 15, 17),
    mar_both = c(14, 15, NA, 16, 14, NA, 15, 16),
    sparse_mar = c(20, NA, NA, NA, 20, 21, 22, 23),
    sparse_mnar = c(8, NA, NA, NA, 14, 15, 16, 17),
    all_mnar = c(8, NA, NA, NA, 9, NA, NA, NA),
    globally_absent = rep(NA, 8),
    complete_low = c(8, 8, 8, 8, 9, 9, 9, 9)
)
colnames(x) <- paste0("s", seq_len(ncol(x)))
condition <- rep(c("A", "B"), each = 4)

fit <- classify_missingness(
    x,
    group = condition,
    cutoffs = c(A = 12, B = 12),
    seed = 1L
)
fit
```

```text
<imputefinder_result>
Retained 4 of 7 features; dropped 3.
Analysed 8 samples in 2 conditions; inserted 1 rescue value.
A: retained features include 2 MNAR, 1 MAR, 1 complete, 0 absent; cutoff 12 (manual)
B: retained features include 0 MNAR, 1 MAR, 3 complete, 0 absent; cutoff 12 (manual)
```

The feature-level summary explains both retention and exclusion:

| Feature | A | B | Retained | Drop reason |
|---|---|---|---:|---|
| `on_off` | MNAR after rescue | complete | yes | |
| `mar_both` | MAR | MAR | yes | |
| `sparse_mar` | insufficient | complete | no | `insufficient:A` |
| `sparse_mnar` | MNAR | complete | yes | |
| `all_mnar` | MNAR | MNAR | no | `MNAR_all_conditions` |
| `globally_absent` | not classified | not classified | no | `all_missing` |
| `complete_low` | complete | complete | yes | |

### What pooling conditions would lose

Apply the same cutoff and strict-majority rule to the original `on_off` row as
one pooled eight-sample block:

| Scope | Observed | Mean intensity | State | Consequence |
|---|---:|---:|---|---|
| Pooled original samples | 4/8 | 15.75 | `insufficient` | Drop |
| Condition A after rescue | 1/4 | 8.00 | `MNAR` | Retain because B is complete |
| Condition B | 4/4 | 15.75 | `complete` | Retain |

The pooled block is above the cutoff but has no strict observed majority, so it
would be dropped. ImputeFinder instead records the feature as MNAR in A and
complete in B, then retains it. The condition-specific states preserve an
on/off pattern that pooling would lose.

Inspect the classifications and condition-specific downstream groups directly:

```r
fit$classifications[
    , c("feature", "condition", "observed_count", "mean_intensity",
        "cutoff", "state", "seeded", "retained")
]

fit$groups$A$MNAR
#> [1] "on_off"     "sparse_mnar"

fit$groups$A$MAR_or_complete
#> [1] "mar_both"    "complete_low"
```

The input `x` remains unchanged. `fit$data` contains the four retained rows and
the recorded rescue value.

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
