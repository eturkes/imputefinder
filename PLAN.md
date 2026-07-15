# ImputeFinder completion plan

Status: canonical implementation plan for Codex

Baseline reviewed: `eturkes/imputefinder` `main` at commit `50c8154` (2026-07-14). Re-check the live repository before editing because this snapshot may become stale.

Source constraint: Codex will not receive the thesis. This document is intentionally self-contained. Treat the normative behavior below as the source of truth; thesis page references are provenance for the maintainer, not an implementation dependency.

## 1. Operating rules for multi-session work

`PLAN.md` defines stable product behavior and completion gates. Use `.agent/roadmap.md` as the live session ledger required by `AGENTS.md`; create it if absent. Do not weaken or reinterpret a normative requirement here merely to fit the current code.

At the start of every session:

1. Read `AGENTS.md`, `PLAN.md`, `.agent/roadmap.md`, and `.agent/memory.md` if present.
2. Run `git status`; inspect tracked source, package metadata, tests, and the latest commits.
3. Select one cohesive milestone or sub-milestone. Split work rather than combining unrelated changes.
4. Establish the failing test or explicit verification that defines the session outcome.

At the end of every session:

1. Run the narrow tests for the changed behavior, then the broadest affordable package checks.
2. Update `.agent/roadmap.md` with completed work, remaining blockers, commands run, and the exact next task.
3. Update this file only for status checkboxes, approved decisions, or genuinely improved acceptance criteria. Preserve its role as the method specification.
4. Make one scoped commit for the cohesive change. Leave the working tree clean unless a blocker is explicitly recorded.

If a requirement remains ambiguous after reading this plan and the tests, stop and ask the maintainer. Record the resulting decision here or in a short tracked decision record. Do not use the current prototype's behavior as an implicit answer where it conflicts with this plan.

## 2. Target outcome

Deliver a usable, tested R/Bioconductor package that implements ImputeFinder as a condition-aware missingness classifier and workflow coordinator for quantitative proteomics.

The package must:

- accept unnormalised log2 protein-intensity data plus sample-to-condition assignments;
- preserve biologically plausible condition-specific "on/off" proteins;
- classify each feature separately in each experimental condition;
- distinguish likely MNAR from MAR-or-complete cases using an intensity/missingness cutoff;
- apply the majority-observed rule only to MAR candidates;
- reconcile classifications across any number of conditions without losing condition-specific state;
- return filtered, seed-modified data plus auditable classifications, cutoffs, diagnostics, and seed provenance;
- leave normalisation and actual imputation to downstream tools;
- run cleanly from a fresh installation, with documentation and checks suitable for Bioconductor preparation.

Priority order:

1. Scientific and logical correctness.
2. Reproducibility and auditability.
3. Stable, clear public API.
4. Robust automatic cutoff selection.
5. Packaging, documentation, and release polish.

## 3. Explicit non-goals

The core package does not perform kNN, MinProb, MinDet, QRILC, BPCA, normalisation, differential-abundance testing, or method recommendation. It classifies and prepares data so users can apply those methods per condition.

The first release does not need:

- built-in normalisation pipelines;
- built-in imputation implementations;
- a web application;
- instrument-specific learned models;
- a universal/global cutoff shared across conditions;
- exact reproduction of unpublished or unavailable thesis datasets;
- broad support for every matrix backend before ordinary numeric matrices and `SummarizedExperiment` objects are correct.

Keep these as future work, not reasons to delay a correct classifier.

## 4. Current baseline and defects to replace

The current repository is a partial prototype, not a safe base for incremental patching without tests.

Known baseline defects:

- `classify_missingness()` calls an absent `plot_detect_custom()`.
- `ggplot_build()` is unimported and is incorrectly used as an algorithmic data source.
- the final MAR expression references `final_MAR` while defining it and uses an undefined `%||%` operator;
- the fixed `min_non_na = 5` only matches an eight-sample condition and is not a general majority rule;
- the function returns global MAR/MNAR unions, discarding the condition-specific classification that downstream mixed imputation requires;
- the fully missing-condition rescue/seed step is absent;
- proteins classified MNAR in every condition are retained rather than excluded;
- the `condition` metadata column is hard-coded;
- row order is repeatedly mutated inside the condition loop;
- automatic cutoff detection can produce empty candidates, `Inf`, division by zero, and unstable derivative results;
- input validation is minimal;
- package metadata is placeholder text and the declared licence is invalid;
- several imported packages are unused;
- there is no test suite, CI, vignette, worked example, or scientific validation fixture.

Disposition of prototype pieces:

| Prototype element | Required disposition |
|---|---|
| `classify_missingness` function name | Retain unless a documented API review finds a materially better name; preserving it minimizes needless churn. |
| Single 97-line implementation | Replace with small, pure helpers and one public orchestration function. |
| `plot_detect_custom` concept | Reimplement internally from explicit profile data; do not copy plot-built coordinates back into the algorithm. |
| Public `threshold` derivative parameter | Remove. Automatic selection must not expose an arbitrary second-derivative tuning knob. |
| Public fixed `min_non_na` default | Remove or deprecate. Derive strict majority from each condition's sample count. |
| Global `MAR`/`MNAR` output vectors | Replace with per-condition output plus an all-feature audit table. |
| `Reduce`/union reconciliation | Delete and replace with an explicit feature-by-condition state matrix. |

The package is unreleased and the current function is not cleanly executable, so correctness takes precedence over preserving broken behavior. Before deliberately breaking arguments, search public dependent code; add a concise deprecation shim only when real usage justifies it.

## 5. Normative method specification

### 5.1 Terminology

A feature is a protein row. A sample is a matrix column. A condition is an experimental group assigned to one or more samples.

ImputeFinder's MAR/MNAR labels are evidence-based heuristics, not proof of the true missingness mechanism. Human-facing documentation must use wording such as "likely MNAR" and "likely MAR" and expose the diagnostics behind the classification.

### 5.2 Input contract

Canonical data:

- numeric matrix-like object with features in rows and samples in columns;
- values are unnormalised log2 protein intensities;
- missing observations are `NA`;
- rows and columns have non-empty, unique names;
- condition labels align one-to-one with columns;
- condition labels contain no missing values;
- at least two conditions are expected for the complete ImputeFinder workflow;
- each condition must contain at least one finite intensity somewhere in the matrix, otherwise its minimum replacement value is undefined and the call must fail clearly.

Do not transform, log, normalise, centre, scale, or reinterpret zeros. The package cannot reliably infer whether input is log2 transformed; document the requirement rather than silently guessing.

Reject or clearly normalize unsupported special values:

- reject `Inf` and `-Inf`;
- reject `NaN` rather than silently treating it as a conventional `NA`;
- reject duplicated or absent feature/sample names;
- reject group/data mismatches;
- reject a missing requested assay or ambiguous multiple-assay selection.

Required interfaces:

1. A primary matrix interface accepting a group vector or an aligned design data frame and group column.
2. A `SummarizedExperiment` adapter with configurable assay selection and configurable group column; no hard-coded `condition` requirement.
3. Both interfaces must call the same pure matrix core and produce equivalent classifications.

The input object must not be modified by reference or through side effects.

### 5.3 Stable state vocabulary

For every feature `i` and condition `g`, calculate state after the rescue step:

- `complete`: zero missing values in the condition. Complete cases need no missingness mechanism, regardless of mean intensity.
- `MNAR`: at least one value is missing and the observed mean intensity is strictly below the condition cutoff.
- `MAR`: at least one value is missing, the observed mean is greater than or equal to the cutoff, and a strict majority of samples are observed.
- `insufficient`: at least one value is missing, the observed mean is greater than or equal to the cutoff, and no strict majority is observed.

This explicit `complete` state resolves an ambiguity in the original shorthand. MNAR is a missingness mechanism; a fully observed low-intensity feature must not be discarded merely because its mean lies below the cutoff.

For a condition with `n_g` samples, strict majority means:

```text
minimum observed = floor(n_g / 2) + 1
```

Examples: 3/4, 5/8, and 11/20 pass; 2/4, 4/8, and 10/20 fail.

### 5.4 Step A - remove globally absent features

Before any replacement, identify features missing in every sample across every condition.

- Mark them `all_missing` in the feature audit.
- Exclude them from all later calculations.
- Never seed them: the method requires evidence that the feature was detected in at least one condition.

### 5.5 Step B - rescue a feature fully missing within one condition

This is the defining "on/off protein" behavior and is mandatory.

For each condition `g`:

1. Compute `condition_min[g]` once as the minimum finite intensity in that condition before any seed values are inserted.
2. Find every non-globally-absent feature whose values are all `NA` in `g`.
3. For each such feature, choose one sample column from condition `g` uniformly at random. Sampling is with replacement across features, so the same sample may be selected for several features.
4. Replace exactly that one cell with `condition_min[g]`; leave the other cells in the feature-condition block as `NA`.

Required reproducibility behavior:

- default reproducibility seed is `1L`, matching the original analysis intent;
- use a local RNG scope: the caller's global RNG state must be identical before and after the function call;
- process condition labels, feature identifiers, and candidate sample identifiers in stable order so row/column permutation does not change classification or seed assignment by sample name under the same seed;
- return a seed log containing feature, condition, selected sample, old value, inserted value, and seed;
- preserve all original observed values exactly.

The old ad hoc script reset `set.seed(1)` while processing conditions. Exact old random cell choices are not a compatibility target; uniform selection, reproducibility, order stability, and provenance are.

### 5.6 Step C - calculate condition-specific feature statistics

After rescue, calculate for every non-globally-absent feature and condition:

- sample count;
- observed count;
- missing count;
- missing fraction;
- `has_missing` boolean;
- arithmetic mean of observed intensities;
- whether a seed was inserted.

Use the arithmetic mean, not the median. A high observed intensity in one sample is intended to pull the summary upward and make a MAR classification more likely.

No surviving feature-condition block may have a non-finite mean: globally absent features were removed and condition-specific all-missing blocks were seeded.

### 5.7 Step D - construct the missingness/intensity profile

The cutoff diagnostic is calculated independently for every condition.

The profile corresponding to the thesis figure is a stacked density of feature mean intensity split by whether the feature has any missing values in that condition. It is not computed from rendered plot coordinates.

For condition `g`:

1. Split feature means into `has_missing = TRUE` and `FALSE` groups.
2. Estimate count-weighted density curves on a common intensity grid.
3. Let `n_1`, `f_1(x)` describe features with missing values and `n_0`, `f_0(x)` complete features.
4. Calculate the plotted missing proportion explicitly:

```text
p_missing(x) = n_1 * f_1(x) / (n_1 * f_1(x) + n_0 * f_0(x))
```

This is the explicit equivalent of mapping density count and using a filled stacked density plot. The low-intensity region should generally be dominated by features with missing values; the profile then falls toward a noisier low-missing region at higher intensity.

Store the raw per-feature statistics and the derived profile grid. The plot is a presentation of those data, not a source of truth for them.

Diagnostics should also expose raw missing fractions, seeded features, feature counts in each class, density/smoothing choices, and any cutoff-quality warnings.

### 5.8 Step E - determine a cutoff per condition

Manual and automatic paths are both required.

Manual path:

- accept a named numeric cutoff for any condition with incomplete features;
- names must match condition labels exactly;
- validate finite values and report suspicious values outside the observed feature-mean range;
- record source as `manual`.

Automatic path:

- is the default only after the validation milestone in this plan passes;
- operates on the explicit `p_missing(x)` profile, never on `ggplot_build()` output;
- finds the dominant descending "MNAR cliff" and selects the right-hand/bottom boundary of that transition, not the first noisy curvature point;
- derives each condition independently;
- exposes no arbitrary public derivative threshold;
- records method version, quality statistics, and warnings;
- never silently returns `Inf`, an endpoint caused by no candidate, or a cutoff from a flat/unidentifiable profile.

When automatic selection is not credible, raise a structured, condition-specific error that tells the user to inspect the diagnostic profile and supply a manual cutoff. A deterministic failure is preferable to a plausible-looking fabricated cutoff.

A condition with no missing values requires no cutoff; mark it `NA` with source `not_needed`, and classify every feature as complete in that condition.

### 5.9 Step F - assign per-condition states

For every feature-condition block:

```text
if missing_count == 0:
    state = complete
else if mean_intensity < cutoff:
    state = MNAR
else if observed_count > sample_count / 2:
    state = MAR
else:
    state = insufficient
```

Boundary rule: a mean exactly equal to the cutoff is on the MAR side.

Sparse MNAR is valid. A single observed value can be sufficient for an MNAR feature because the intended downstream methods can operate from a low-value anchor. Do not apply the majority rule to MNAR.

### 5.10 Step G - reconcile across conditions

Generalize to any number of conditions using the full feature-by-condition state matrix. Do not write pairwise `Neg`/`Pos` special cases.

Retain feature `i` if and only if both conditions hold:

```text
1. Every condition state is one of complete, MAR, or MNAR.
2. At least one condition state is complete or MAR.
```

Equivalent interpretation:

- drop a feature if any condition is `insufficient`;
- drop a feature if it is MNAR in every condition;
- retain complete/MAR features when all other conditions are complete, MAR, or MNAR;
- retain condition-specific on/off features, such as MNAR in disease and complete/MAR in control.

Record an auditable drop reason:

- `all_missing`;
- `insufficient:<condition labels>`;
- `MNAR_all_conditions`.

If multiple reasons could apply, use deterministic precedence: `all_missing`, then `insufficient`, then `MNAR_all_conditions`.

### 5.11 Step H - return condition-specific classifications

For every condition, return disjoint retained-feature sets:

- `MNAR`;
- `MAR`;
- `complete`;
- `MAR_or_complete = union(MAR, complete)`.

For every retained feature in every condition, exactly one of `MNAR` and `MAR_or_complete` must be true. These are the lists a downstream mixed-imputation workflow needs.

Never collapse these lists into global unions as the primary result. A feature may be MNAR in one condition and MAR/complete in another.

### 5.12 Pipeline order after ImputeFinder

Document the intended external workflow:

1. Run ImputeFinder on unnormalised log2 intensities.
2. Use the returned filtered, seed-modified data.
3. Apply the user's chosen normalisation.
4. Split columns by condition.
5. Within each condition, apply an MNAR method to rows in that condition's `MNAR` set and a MAR method to rows in its `MAR_or_complete` set.
6. Recombine conditions without changing sample order.

The package may show examples using kNN for MAR and MinProb for MNAR, but these are examples, not imposed defaults or runtime dependencies.

## 6. Required public result contract

Return a documented list or lightweight S3 object, provisionally classed `imputefinder_result`, with at least:

```text
data
    Filtered, seed-modified object in the same broad representation as input:
    matrix in -> matrix out; SummarizedExperiment in -> SummarizedExperiment out.

classifications
    Long data frame with one row per original feature x condition, where applicable:
    feature, condition, sample_count, observed_count, missing_count,
    missing_fraction, mean_intensity, cutoff, state, seeded, retained, drop_reason.

groups
    Named list by condition. Each entry contains retained feature names:
    MNAR, MAR, complete, MAR_or_complete.

feature_status
    One row per original feature with retained flag and deterministic drop reason.

cutoffs
    Named numeric vector by condition.

cutoff_diagnostics
    Named list by condition containing source (manual/automatic/not_needed),
    method identifier/version, quality measures, warnings, and profile metadata.

profiles
    Named list by condition containing raw feature statistics and explicit profile grid.

seed_log
    Data frame of all artificial seed insertions.

groups_by_sample
    Named condition vector aligned to output columns.

call
    Matched call and relevant reproducibility metadata.
```

Recommended methods:

- `print.imputefinder_result()` - concise feature/group/cutoff summary;
- `summary.imputefinder_result()` - state and drop counts by condition;
- `plot_missingness(result, condition)` or a `plot()` method - reconstruct plot from stored profile data and add cutoff line.

Do not store a second full copy of the original matrix merely for convenience. The caller's input remains unchanged, while `seed_log` and `feature_status` provide provenance.

Public signature decided in M0:

```r
classify_missingness <- function(
    x,
    group = NULL,
    group_col = NULL,
    assay = NULL,
    cutoffs = NULL,
    seed = 1L
)
```

Argument routing is explicit:

- matrix `x` + atomic `group`: unnamed groups align positionally; named groups must have exactly the sample names and are aligned by name; `group_col` and `assay` must be `NULL`;
- matrix `x` + data-frame `group`: unique row names must equal the sample names, rows are aligned by name, `group_col` selects exactly one column, and `assay` must be `NULL`;
- `SummarizedExperiment` `x`: `group` is omitted, `group_col` explicitly selects `colData(x)`, and `assay` is optional only when exactly one assay exists;
- `cutoffs` is `NULL` for automatic selection or a named, optionally partial numeric vector of manual condition cutoffs; automatic selection fills conditions without a manual value;
- `seed` is one non-missing integer and defaults to `1L`.

`threshold`, `min_non_na`, and `return_plot` are removed. Stored profiles plus the plot method replace `return_plot`. M0 found no public dependent usage, releases, tags, or forks, so no compatibility shim is warranted.

Illustrative use:

```r
fit <- classify_missingness(
    x = intensity_matrix,
    group = design$condition,
    cutoffs = c(Neg = 12.4, Pos = 12.4)
)

fit$groups$Neg$MNAR
fit$groups$Neg$MAR_or_complete
plot_missingness(fit, "Neg")
```

For `SummarizedExperiment`, users must be able to select the assay and metadata column explicitly.

## 7. Normative worked fixture

Use this fixture in tests and documentation so the central semantics cannot drift. Conditions A and B each have four samples; manual cutoffs are `A = 12`, `B = 12`.

```r
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
group <- rep(c("A", "B"), each = 4)
```

Expected outcome:

| Feature | A state | B state | Retained | Reason |
|---|---|---|---:|---|
| `on_off` | MNAR after one A seed | complete | yes | condition-specific on/off protein |
| `mar_both` | MAR (3/4 observed) | MAR (3/4 observed) | yes | adequate MAR in both |
| `sparse_mar` | insufficient (1/4 observed, high mean) | complete | no | insufficient:A |
| `sparse_mnar` | MNAR (1/4 observed, low mean) | complete | yes | sparse MNAR allowed |
| `all_mnar` | MNAR | MNAR | no | MNAR_all_conditions |
| `globally_absent` | not seeded | not seeded | no | all_missing |
| `complete_low` | complete | complete | yes | complete overrides low-intensity mechanism shorthand |

Additional assertions:

- A strict majority for four samples is three.
- `on_off` receives exactly one seed in A, equal to A's pre-seed condition minimum.
- the original matrix remains byte-for-byte unchanged;
- output row/column order follows the original order, minus dropped rows;
- condition-specific lists agree with the state table.

## 8. Implementation architecture

Prefer small, pure functions and base R matrix operations. Avoid a tidyverse dependency for simple reshaping or aggregation.

Suggested source split:

```text
R/
  classify_missingness.R       public orchestration
  input.R                      validation + matrix/SE adapters
  seed_missing_conditions.R    global-absence and rescue logic
  profile_missingness.R        feature statistics + density profile
  detect_cutoff.R              manual validation + automatic detector
  reconcile.R                  state matrix + retention rules
  result.R                     constructor, print, summary
  plot.R                       plots generated from stored profile data
  package.R                    package-level documentation
```

Suggested tests:

```text
tests/testthat/
  helper-fixtures.R
  test-input.R
  test-seeding.R
  test-profile.R
  test-cutoff-manual.R
  test-cutoff-auto.R
  test-classification.R
  test-reconciliation.R
  test-output.R
  test-summarizedexperiment.R
  test-plot.R
  test-invariants.R
```

Architecture requirements:

- one internal matrix core;
- no feature reordering during condition loops;
- stable named indexing rather than implicit positional joins;
- no extraction of scientific data from ggplot internals;
- no undefined convenience operators;
- no mutable global state;
- no silent fallbacks on invalid cutoffs;
- no condition-name special cases;
- all generated `NAMESPACE` and Rd files maintained through roxygen2;
- runtime dependencies kept minimal and justified.

Target runtime dependencies should normally be limited to base/recommended R, `ggplot2`, and Bioconductor infrastructure needed for `SummarizedExperiment`. Remove `dplyr`, `tidyr`, and `assertthat` unless later work proves a real, documented need.

## 9. Milestones and suggested Codex sessions

The session boundaries below are suggestions. Split any milestone further when that produces a cleaner test-driven commit.

### M0 - Baseline, tests, and API decision record

- [x] Re-audit current `main`; run the existing package build/check and record actual failures.
- [x] Add `testthat` edition 3 infrastructure and the normative fixture.
- [x] Add red tests for the central current failures: absent plot helper, fully missing condition, dynamic majority, per-condition output, and all-MNAR exclusion.
- [x] Decide and document the exact public function signature while preserving all capabilities required in Sections 5-6.
- [x] Decide whether to add a temporary deprecation shim after searching public dependent code.
- [x] Create `.agent/roadmap.md` and map the remaining milestones.

Gate: tests express the intended method even though implementation tests are still red; package baseline failures are documented rather than guessed.

### M1 - Input core and result skeleton

- [x] Implement matrix validation and aligned group handling.
- [x] Implement `SummarizedExperiment` extraction/reconstruction with configurable assay and group column.
- [x] Reject unsupported values and ambiguous inputs with specific messages.
- [x] Create the internal state vocabulary and result constructor.
- [x] Guarantee input immutability and original order preservation.

Gate: matrix and `SummarizedExperiment` adapters feed identical matrices/groups into the core; validation tests pass.

### M2 - Global absence and condition-specific rescue

- [x] Implement `all_missing` detection and audit status.
- [x] Compute pre-seed condition minima.
- [x] Insert exactly one condition-minimum value for each fully missing feature-condition block that is observed elsewhere.
- [x] Implement local seeded randomness, stable processing order, and a complete seed log.
- [x] Verify the caller's RNG state is unchanged.

Gate: normative `on_off` and `globally_absent` tests pass; row/column permutations preserve named results under the same seed.

### M3 - Manual-cutoff classification and reconciliation

- [x] Calculate per-feature/per-condition statistics after seeding.
- [x] Validate named manual cutoffs.
- [x] Implement `complete`, `MNAR`, `MAR`, and `insufficient` states exactly as specified.
- [x] Derive strict majority separately for every condition.
- [x] Implement n-condition reconciliation and deterministic drop reasons.
- [x] Populate per-condition `MNAR`, `MAR`, `complete`, and `MAR_or_complete` lists.
- [x] Remove the broken global union logic.

Gate: all normative fixture expectations pass using manual cutoffs, including unequal group sizes and more than two conditions.

### M4 - Explicit profile data and plotting

- [x] Recreate the stacked missing-vs-complete density profile from explicit calculations.
- [x] Store raw statistics and profile grids in the result.
- [x] Implement plotting from those stored data, including condition title, percentage axis, and cutoff line.
- [x] Mark seeded features and warnings in diagnostics where useful without cluttering the default plot.
- [x] Remove all algorithmic use of `ggplot_build()`.

Gate: profile calculations are unit tested numerically; plots build without missing functions or namespace leaks.

### M5 - Automatic cutoff research and implementation

- [x] Build deterministic synthetic profiles covering clear, noisy, weak, absent, and multi-transition cliffs.
- [x] Benchmark at least segmented/change-point and derivative-based candidates against the criteria in Section 11.
- [x] Select the simplest method that reliably locates the right-hand/bottom cliff boundary.
- [x] Implement quality scoring and structured failure.
- [x] Record algorithm identity/version in results.
- [x] Remove the public `threshold` argument.
- [x] Keep manual cutoffs as a permanent reproducibility override.

Gate: automatic selection passes the scientific/stability acceptance tests and never fabricates a cutoff for an unidentifiable profile.

### M6 - Public API integration and compatibility cleanup

- [x] Integrate manual/automatic paths into the public function.
- [x] Implement `print`, `summary`, and missingness plot methods/helpers.
- [x] Ensure matrix and `SummarizedExperiment` outputs preserve metadata and order.
- [x] Remove obsolete prototype code and unused imports.
- [x] Add any justified deprecation behavior with tests and clear messages (none
  justified by the M0 public-usage audit).

Gate: the installed package's public API works from a clean R session and exposes no unqualified/missing symbols.

### M7 - Scientific regression suite

- [x] Add compact simulations of uniform MAR plus intensity-dependent MNAR.
- [x] Add condition-specific on/off simulations.
- [x] Test classification accuracy, retention, and cutoff error across seeds and permutations.
- [x] Test sensitivity to group sizes 4, 8, and 20.
- [x] Test a known cliff near intensity 12 and cutoff sweeps analogous to 8-14.
- [x] Keep long benchmarks outside routine unit tests; store concise deterministic summaries.

Gate: the package has evidence for its central scientific behavior, not only code-path coverage.

### M8 - Documentation and package metadata

- [x] Replace all `DESCRIPTION` placeholders with accurate title, author, description, GPL metadata, URLs, bug tracker, and current Bioconductor fields.
- [x] Audit `Imports`, `Suggests`, `biocViews`, and minimum R/Bioconductor versions against current guidance.
- [x] Rewrite README with method rationale, installation, minimal example, result interpretation, and limitations.
- [x] Add a vignette covering manual and automatic cutoffs, the on/off fixture, plotting, filtering, and per-condition downstream imputation.
- [x] Document that input is unnormalised log2 intensity and that classification is heuristic.
- [x] Document seed provenance and pipeline order before normalisation/imputation.
- [x] Add package-level help and `NEWS.md`.
- [x] Ensure every example runs offline or is explicitly non-evaluated for optional downstream packages.

Gate: a user with no thesis can understand and correctly use the package from installed documentation alone.

### M9 - CI and Bioconductor hardening

- [x] Add current R package CI using maintained actions and an appropriate Bioconductor environment.
- [x] Run unit tests, `R CMD build`, and `R CMD check --as-cran` from a clean source tarball.
- [x] Run `BiocCheck` using current official guidance.
- [ ] Fix all errors and warnings; resolve actionable notes rather than suppressing them.
- [x] Verify licence files, line endings, generated documentation, file sizes, examples, and vignette build.
- [x] Review code coverage; omit it because existing focused gates diagnose failures without fragile infrastructure.

Gate: clean package checks in CI and locally, with any unavoidable note documented precisely.

### M10 - Release-candidate adversarial review

- [x] Review every guarantee in README/Rd/vignette against tests and implementation.
- [x] Re-test no-missing data, one-class profiles, flat profiles, duplicated means, tiny groups, and all-special-value failures.
- [ ] Verify no global RNG, options, graphics, or working-directory side effects.
- [ ] Verify output invariance under row, column, and condition ordering by names.
- [ ] Verify observed cells never change except logged seed insertions.
- [ ] Benchmark a representative matrix (for example 10,000 features x 50 samples) for avoidable copying or pathological runtime.
- [ ] Remove dead development artefacts and stale comments.
- [ ] Set the first usable development/release version consistently, provisionally `0.1.0` or the Bioconductor-appropriate equivalent.

Gate: all Definition of Done items pass from a fresh checkout and clean library.

## 10. Required test matrix

At minimum, tests must cover the following behavior.

### Input and adapter tests

- [x] matrix + vector groups;
- [x] matrix + aligned design data frame;
- [x] `SummarizedExperiment` + explicit group column;
- [x] one assay vs multiple assays requiring explicit choice;
- [x] missing, duplicate, and mismatched feature/sample names;
- [x] missing group labels;
- [x] nonnumeric assay;
- [x] `Inf`, `-Inf`, and `NaN` rejection;
- [x] zeros and negative finite log2 intensities preserved;
- [x] input object unchanged.

### Rescue tests

- [x] feature all missing globally is dropped and never seeded;
- [x] feature all missing in one condition but observed elsewhere gets exactly one seed;
- [x] seed equals the condition-wide minimum computed before any seeds;
- [x] different rescued features may select the same sample;
- [x] seed assignments are reproducible and logged;
- [x] caller RNG state is unchanged;
- [x] a condition with no finite value anywhere errors clearly;
- [x] sample/feature permutation retains assignments by name under the same seed.

### Classification tests

- [x] complete overrides low mean;
- [x] incomplete mean below cutoff is MNAR;
- [x] incomplete mean equal to cutoff is on the MAR side;
- [x] incomplete mean above cutoff with strict majority is MAR;
- [x] incomplete mean above cutoff at exactly half observed is insufficient;
- [x] sparse MNAR with one observation remains eligible;
- [x] arithmetic mean, not median, drives boundary behavior;
- [x] no-missing condition produces complete states and needs no cutoff.

### Reconciliation tests

- [x] complete/MAR in all groups retained;
- [x] MNAR in one group and complete/MAR in another retained;
- [x] MNAR in every group dropped;
- [x] insufficient in any group dropped;
- [x] globally absent drop reason has precedence;
- [x] rules generalize to 3+ conditions;
- [x] convenience lists contain retained features only and agree with the long table;
- [x] no retained feature is simultaneously in `MNAR` and `MAR_or_complete` within a condition;
- [x] every retained feature appears in exactly one of those two sets per condition.

### Profile and cutoff tests

- [x] explicit profile matches the count-weighted density formula;
- [x] profile computation handles repeated means without division-by-zero artifacts;
- [x] manual cutoff validation and source metadata;
- [x] automatic clear-cliff accuracy;
- [x] automatic condition-specific cutoffs;
- [x] automatic stability under permutations and modest bootstrap noise;
- [x] flat/unidentifiable profile returns structured failure;
- [x] only-complete and only-missing profile behavior is explicit;
- [x] no candidate ever becomes `Inf`, `-Inf`, or silent endpoint fallback;
- [x] plot is generated from stored profile data and includes the recorded cutoff.

### Package tests

- [x] fresh-session namespace contains every called symbol;
- [x] examples run;
- [x] vignette builds;
- [x] matrix and `SummarizedExperiment` core results are equivalent;
- [x] result print/summary methods are concise and accurate;
- [x] source tarball check is clean.

## 11. Automatic cutoff validation protocol

Automatic cutoff selection is the least settled part of the original method. Treat it as a scientific component requiring comparative evidence, not a cosmetic refactor.

### 11.1 Candidate methods

Evaluate at least:

1. Segmented or piecewise regression/change-point detection on `p_missing(x)`.
2. A robust derivative/curvature method on the explicit profile, with smoothing and boundary rules fixed internally.
3. Optionally, monotone/isotonic smoothing followed by dominant transition-boundary detection when it materially improves robustness.

Reject an approach that only works because of an exposed hand-tuned threshold.

### 11.2 Synthetic scenarios

Generate deterministic scenarios with known ground truth:

- sharp single cliff;
- broad cliff;
- low-amplitude cliff;
- noisy right tail;
- heavy MAR background (approximately 25% random missingness);
- sparse MAR background (approximately 5%);
- no cliff/flat random missingness;
- two apparent transitions;
- unequal class counts;
- duplicate feature means;
- different cutoffs by condition;
- small feature count near the method's minimum viable input.

MNAR amputation should become more probable as intensity falls below a known boundary. MAR amputation should be independent of intensity.

### 11.3 Selection metrics

Measure:

- absolute cutoff error from the known right-hand cliff boundary;
- bias toward the left/top vs right/bottom of the transition;
- MNAR/MAR classification precision and recall on missing feature-condition blocks;
- retained-feature precision/recall under the cross-condition rule;
- stability across seeds, row order, column order, and bootstrap resampling;
- false confidence rate on no-cliff data;
- runtime and complexity.

### 11.4 Acceptance criteria

The chosen method must:

- select the right/bottom region on clear simulated cliffs;
- remain within a predeclared tolerance across realistic noise, with tolerance defined in the validation report before final comparison;
- fail explicitly on unidentifiable profiles rather than guessing;
- produce condition-specific results;
- use fixed internal behavior with no public derivative threshold;
- be simpler than a competing method when performance is materially equivalent;
- preserve a manual override for exact reproducibility.

Store the experiment and conclusion in a concise tracked development report, for example `dev/cutoff-validation.md` plus reproducible script. Do not put long stochastic benchmarks in ordinary package checks.

## 12. Scientific validation beyond cutoff placement

Create a compact, independent simulation inspired by the intended use case; no thesis file or private data may be required.

Recommended generator:

1. Generate complete log2 intensities with realistic feature-specific means/variances and multiple conditions.
2. Introduce condition-specific abundance differences, including suppressed "off" features.
3. Add intensity-independent MAR at 5% and 25% scenarios.
4. Add intensity-dependent MNAR below known condition cutoffs, with probability increasing as intensity decreases.
5. Run ImputeFinder with known manual cutoffs and automatic cutoffs.
6. Measure missingness classification and retained-feature behavior.
7. Optionally, in a non-core validation script with suggested packages, compare downstream differential-abundance recovery after per-condition kNN/MinProb against unimputed and complete-case analyses.

Do not hard-code or advertise the thesis's exact benchmark percentages unless the exact source data and full analysis are reproducibly included. Claims in package documentation must match evidence available in the repository.

## 13. Documentation requirements

README and vignette must explain, without relying on the thesis:

- why proteomics can contain both MAR-like and MNAR-like missingness;
- why classification is condition-specific;
- why a protein fully absent in one condition can be biologically informative;
- exactly how the one-cell condition-minimum seed works and why it is logged;
- what the stacked density profile means;
- how manual and automatic cutoffs differ;
- why arithmetic mean is used;
- why MAR requires strict majority but MNAR does not;
- why all-MNAR and insufficient-MAR features are excluded;
- the result structure and per-condition lists;
- how to normalize and impute separately by condition after classification;
- that the package does not prove missingness mechanisms or perform imputation itself;
- limitations of intensity-cutoff assumptions and automatic breakpoint selection.

At least one vignette example must show a feature that is MNAR in A and complete/MAR in B and then use different downstream imputation branches by condition.

## 14. Package and Bioconductor requirements

Complete and verify:

- valid `DESCRIPTION` fields and GPL declaration consistent with `COPYING`;
- accurate author/maintainer data without invented identifiers;
- concise title and description focused on condition-aware proteomics missingness classification;
- `URL` and `BugReports` pointing to the repository;
- current, valid `biocViews` selected from official vocabulary;
- minimal imports and complete namespace qualification;
- package-level documentation;
- unit tests and vignette in standard locations;
- installation from source without development-only files leaking into the tarball;
- `R CMD check --as-cran` and `BiocCheck` under current supported R/Bioconductor versions;
- maintained CI configuration rather than copied obsolete action versions.

Codex must consult current official R/Bioconductor documentation when choosing version fields, CI actions, and submission checks because these requirements change over time.

## 15. Definition of Done

The work is complete only when all of the following are true:

- [ ] The fully missing-condition rescue behavior is implemented exactly and audited.
- [ ] Globally absent features are never rescued.
- [ ] State is calculated per feature per condition and retained in output.
- [ ] Complete, MNAR, MAR, and insufficient semantics match Section 5.3.
- [ ] Majority is derived from each condition size.
- [ ] Cross-condition reconciliation works for arbitrary condition counts.
- [ ] All-MNAR features are excluded.
- [ ] Condition-specific output supports separate downstream MAR/MNAR imputation.
- [ ] Manual cutoffs work and automatic cutoffs satisfy Section 11.
- [ ] No algorithm reads data back from a ggplot build object.
- [ ] No missing/undefined namespace symbols remain.
- [ ] Input data and global RNG state are unchanged.
- [ ] Observed values are unchanged except for explicitly logged seeds in returned data.
- [ ] Output order and named results are stable under input permutation.
- [ ] The normative fixture and all required tests pass.
- [ ] README, Rd help, and vignette are sufficient without thesis access.
- [ ] Runtime dependencies are minimal and justified.
- [ ] `R CMD check --as-cran`, package examples, vignette build, and `BiocCheck` pass from a clean checkout.
- [ ] Claims in documentation do not exceed repository evidence.
- [ ] `.agent/roadmap.md` contains no unresolved P0/P1 item for this release.

## 16. Deferred roadmap after first correct release

Do not mix these into the core completion work unless the maintainer explicitly promotes them:

- optional helpers that execute common MAR/MNAR imputation methods;
- optional normalisation orchestration;
- data-driven imputation-method recommendation;
- global vs per-condition cutoff comparison;
- instrument-specific cutoff calibration;
- probabilistic feature-level missingness models beyond a hard cutoff;
- richer HTML reports;
- web/server workflow;
- broader matrix backends and very large delayed arrays;
- independent validation on additional public proteomics datasets.
