# ImputeFinder

ImputeFinder classifies quantitative-proteomics missingness independently in
each experimental condition. It preserves plausible condition-specific on/off
proteins, filters blocks without enough evidence for downstream imputation, and
returns the modified data together with the evidence behind every decision.

The package expects **unnormalised log2 protein intensities**. It classifies and
prepares data; it does not normalise or impute them.

## Why condition-specific classification?

Quantitative proteomics can contain at least two useful missingness patterns.
Scattered observations may be absent for intensity-independent, MAR-like
reasons, while low-abundance proteins may drop below detection in an MNAR-like
way. The same protein can also be detected in one biological condition and
absent in another. Collapsing conditions into one global label loses that
on/off signal.

ImputeFinder therefore assigns one state to each feature-condition block:

| State | Evidence after rescue | Downstream group |
|---|---|---|
| `complete` | No values are missing | `MAR_or_complete` |
| `MNAR` | Values are missing and the observed arithmetic mean is below the condition cutoff | `MNAR` |
| `MAR` | Values are missing, the mean is at or above the cutoff, and strictly more than half of the samples are observed | `MAR_or_complete` |
| `insufficient` | Values are missing, the mean is at or above the cutoff, and no strict observed majority exists | Feature is dropped |

A mean equal to the cutoff is on the MAR side. The majority rule applies only
to MAR candidates: one low-intensity observation can still support a sparse
MNAR block. A feature is retained only when no condition is `insufficient` and
at least one condition is `complete` or `MAR`. Features absent globally and
features classified as `MNAR` in every condition are dropped.

These labels mean **likely MNAR** and **likely MAR** given the observed pattern.
They are auditable heuristics, not proof of a missingness mechanism.

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
repository.

## Minimal worked example

The fixture below exposes all central rules. Conditions A and B each contain
four samples, and both manual cutoffs are 12.

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
Features: 4/7 retained; 3 dropped
Samples: 8 across 2 conditions; seed insertions: 1
A: retained MNAR=2, MAR=1, complete=1; cutoff=12 (manual)
B: retained MNAR=0, MAR=1, complete=3; cutoff=12 (manual)
```

The complete feature audit explains both retention and exclusion:

| Feature | A | B | Retained | Drop reason |
|---|---|---|---:|---|
| `on_off` | MNAR after rescue | complete | yes | |
| `mar_both` | MAR | MAR | yes | |
| `sparse_mar` | insufficient | complete | no | `insufficient:A` |
| `sparse_mnar` | MNAR | complete | yes | |
| `all_mnar` | MNAR | MNAR | no | `MNAR_all_conditions` |
| `globally_absent` | not classified | not classified | no | `all_missing` |
| `complete_low` | complete | complete | yes | |

Inspect the long evidence table and condition-specific downstream groups
directly:

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
the one logged rescue value.

## Rescue and seed provenance

An on/off feature has no observed value from which to calculate its mean in the
absent condition. ImputeFinder handles that block before classification:

1. Features missing in every sample are marked `all_missing` and removed.
2. The minimum finite intensity in each condition is calculated once, before
   any artificial values are inserted.
3. For each surviving feature fully missing in a condition, one sample from
   that condition is selected uniformly and receives the condition minimum.
   The remaining cells stay `NA`.
4. Every insertion is recorded in `seed_log`.

For the worked fixture:

```r
fit$seed_log
#>   feature condition sample old_value inserted_value seed
#> 1  on_off         A     s1        NA              8    1
```

The seed is a low-value anchor for classification, not an imputation result.
Sampling is reproducible by feature, condition, and sample name; the default is
`seed = 1L`. The function uses a local RNG scope, so it leaves the caller's RNG
state unchanged. All originally observed values are preserved exactly.

## Manual and automatic cutoffs

Manual cutoffs are the most direct reproducibility control. Supply a named,
finite value for each condition with missing observations, as above. Partial
manual vectors are allowed: automatic selection fills the omitted conditions.
Values outside the observed feature-mean range are retained but produce a
diagnostic warning.

With `cutoffs = NULL`, the default, each condition containing missing values is
assessed independently for an automatic cutoff. A condition without missing
values needs no cutoff:

```r
fit_auto <- classify_missingness(
    intensity_matrix,
    group = sample_condition
)
```

Automatic selection operates on an explicit count-weighted density profile.
At intensity `x`, the plotted curve is the estimated proportion of features
that contain any missing observation:

```text
n_missing f_missing(x)
-----------------------------------------------
n_missing f_missing(x) + n_complete f_complete(x)
```

The detector locates the lower, right-hand boundary of the dominant descending
missingness transition. It records its method version, evidence and warnings in
`cutoff_diagnostics`. If a condition lacks enough evidence or has a flat or
ambiguous profile, the function raises a condition-specific error and asks for
a manual cutoff instead of fabricating one. Tiny fixtures such as the example
above are intentionally run with manual cutoffs.

Plotting reads only the stored profile and cutoff evidence:

```r
plot_missingness(fit, "A")
fit$cutoff_diagnostics$A
```

Bottom rug ticks identify rescued feature-condition blocks. The vertical line
is the recorded cutoff.

## Result contract

`classify_missingness()` returns an `imputefinder_result` with these components:

| Component | Contents |
|---|---|
| `data` | Filtered, seed-modified data; matrix in gives matrix out, while `SummarizedExperiment` input preserves that representation and metadata |
| `classifications` | One evidence row per surviving feature and condition, plus retention and drop status |
| `groups` | Retained `MNAR`, `MAR`, `complete`, and `MAR_or_complete` feature names for every condition |
| `feature_status` | One row per original feature with deterministic retention and drop reason |
| `cutoffs` | Named condition-specific numeric cutoffs |
| `cutoff_diagnostics` | Manual, automatic, or not-needed source; method metadata; quality evidence; warnings |
| `profiles` | Raw feature statistics and the explicit density grid for each condition |
| `seed_log` | Feature, condition, selected sample, old value, inserted value, and seed for every rescue |
| `groups_by_sample` | Condition vector named and ordered by output samples |
| `call` | Matched call for provenance |

`summary(fit)` reports all state and drop counts. `plot_missingness()` reconstructs
a condition profile without re-reading the input or using rendered plot
coordinates as scientific data.

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
2. Use `fit$data`, which excludes unsupported features and includes every logged
   rescue seed.
3. Apply the chosen normalisation to that complete returned object.
4. Split columns by `fit$groups_by_sample`.
5. Within each condition, apply an MNAR method only to rows named in
   `fit$groups[[condition]]$MNAR` and a MAR method only to rows named in
   `fit$groups[[condition]]$MAR_or_complete`.
6. Recombine conditions in the original sample order.

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

Normalisation and actual kNN, MinProb, MinDet, QRILC, BPCA, or other imputation
remain downstream choices, not ImputeFinder defaults or runtime dependencies.

## Limitations and evidence

- ImputeFinder cannot infer whether values are log2 transformed or suitably
  unnormalised. It validates structure and finite values, then trusts the input
  contract.
- The intensity cutoff is a model assumption. Biology, acquisition, and
  preprocessing can produce missingness patterns that do not follow one
  descending transition.
- Arithmetic mean intensity is intentional: one strong observed signal can move
  a partial block toward MAR. It can also be sensitive to unusual observations.
- Automatic selection requires both incomplete and complete feature evidence.
  Structured failure and a manual cutoff are expected for unidentifiable
  conditions.
- A rescue seed is artificial and remains visible in the returned data and log.
  Its role and downstream effect should be included in the analysis record.
- The package classifies and coordinates a workflow. It neither proves causal
  missingness mechanisms nor chooses, runs, or validates an imputation method.

The repository's [scientific validation
report](https://github.com/eturkes/imputefinder/blob/main/dev/scientific-validation.md)
uses deterministic simulations with 5% and 25% intensity-independent MAR,
intensity-dependent MNAR, on/off features, group sizes 4, 8, and 20,
condition-specific cutoffs, rescue seeds, and named permutations. The [cutoff
validation
report](https://github.com/eturkes/imputefinder/blob/main/dev/cutoff-validation.md)
records the automatic-method comparison and failure criteria. These synthetic
results support the tested behavior; they are not a claim of performance on
every proteomics dataset.

The normative behavior and release gates are maintained in
[PLAN.md](https://github.com/eturkes/imputefinder/blob/main/PLAN.md).

## License

ImputeFinder is distributed under GPL-3.0-or-later. See [COPYING](COPYING).
