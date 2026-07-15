# ImputeFinder roadmap

Canonical method + gates → `PLAN.md`. This ledger owns session state, evidence, blockers, and exact next work.

## Current release path

- [x] M0 baseline + red contract + API decision
- [x] M1 input core + result skeleton
  - [x] M1a matrix validation + vector/design group alignment
  - [x] M1b `SummarizedExperiment` adapter + representation-preserving result skeleton
- [x] M2 global absence + condition rescue
  - [x] M2a deterministic rescue core + seed audit
  - [x] M2b local RNG + row/column permutation invariants
- [x] M3 manual classification + reconciliation
  - [x] M3a post-seed statistics + cutoff validation + four states
  - [x] M3b n-condition retention + groups + complete result contract
- [x] M4 explicit profiles + plotting
  - [x] M4a count-weighted profile calculation + result storage
  - [x] M4b profile plotting from stored data
- [x] M5 automatic cutoff research + implementation
  - [x] M5a deterministic simulations + predeclared benchmark protocol
  - [x] M5b candidate benchmark + method decision
  - [x] M5c pure detector + quality/failure contract
- [x] M6 public integration + obsolete dependency cleanup
  - [x] M6a result presentation + installed S3 dispatch
  - [x] M6b representation/order integration audit
  - [x] M6c obsolete code/import + compatibility cleanup
- [x] M7 scientific regression suite
  - [x] M7a independent simulation protocol + frozen gates
  - [x] M7b full scientific benchmark + assessment
  - [x] M7c compact routine regressions + gate closure
- [x] M8 documentation + package metadata
  - [x] M8a metadata + package overview
  - [x] M8b README + runnable public examples
  - [x] M8c vignette + NEWS + documentation gate
- [ ] M9 CI + Bioconductor hardening
  - [x] M9a maintained Bioconductor-aware CI
  - [x] M9b submission checks + clean `--as-cran`
  - [ ] M9c `BiocCheck` findings + hardening gate closure
- [ ] M10 release-candidate adversarial review
  - [x] M10a public-claim + boundary-case audit
  - [x] M10b side-effect + invariance audit
  - [x] M10c performance + release cleanup

## Session ledger

### 2026-07-14 - M0 baseline, tests, API

Scope → establish executable baseline; lock the required API; encode central semantics as expected-red tests.

Baseline at `1ae9106`:

- `R CMD build . --no-manual --no-build-vignettes` → pass.
- Initial `R CMD check` → 1 ERROR because every declared external dependency was absent from the container.
- Repository-local library provisioned against R 4.6 / Bioconductor 3.23; dependency-complete check → 1 WARNING + 2 NOTEs:
  - WARNING: invalid placeholder licence.
  - NOTE: unused `assertthat`, `dplyr`, `methods`, `tidyr` imports.
  - NOTE: undefined `plot_detect_custom` + unimported `ggplot_build`.
- Direct installed-package call → `could not find function "plot_detect_custom"`.
- Prototype also contains an undefined `%||%`, self-referential `final_MAR`, fixed `min_non_na = 5`, global MAR/MNAR unions, row mutation, and no rescue path; see `PLAN.md` Section 4.

API decision:

```r
classify_missingness(
    x,
    group = NULL,
    group_col = NULL,
    assay = NULL,
    cutoffs = NULL,
    seed = 1L
)
```

Routing + behavior are normative in `PLAN.md` Section 6. Compatibility shim → none: GitHub code searches for package/function/install reference returned zero; repository has no releases, tags, forks, or stars; CRAN registry lookup returned 404. Search performed 2026-07-14. Broken, unreleased prototype arguments provide no observed compatibility value.

Expected-red gate:

- normative fixture shared by tests;
- red assertions cover absent helper dependency, condition rescue/audit, dynamic strict majority, condition-specific state, and all-MNAR exclusion;
- failures are expected until M1-M3 implement the decided contract.

Verification + completion evidence:

- `testthat::test_local(".", reporter = "summary", load_package = "source")` → expected red: 1 helper-dependency failure + 4 new-API errors.
- `R CMD build . --no-manual --no-build-vignettes` after test addition → pass.
- source-tarball `R CMD check --no-manual --no-build-vignettes` → expected 1 ERROR (5 red tests), plus the baseline licence WARNING + 2 dependency/code NOTEs above.
- `git diff --check` → pass.

Exact next task after M0 → M1a: add focused red validation/alignment tests, then implement pure matrix input validation and vector/design group extraction without classification logic.

Blockers → none.

### 2026-07-14 - M1a matrix input core

Scope → validate ordinary numeric matrices and resolve vector/design conditions into one sample-aligned representation; classification and `SummarizedExperiment` routing remain outside this commit.

Implementation:

- pure `.prepare_matrix_input()` boundary returns unchanged matrix data, sample-named character conditions, and representation provenance;
- unnamed atomic groups align positionally; fully named groups align by exact sample-name set;
- design data frames require exact, unique sample row names and exactly one explicit `group_col`;
- rejects non-matrices/non-numeric matrices, absent/empty/duplicate feature or sample names, `NaN`, infinities, missing/empty conditions, mismatched identifiers, and matrix-inapplicable routing arguments;
- finite zero and negative values remain untouched.

Verification + completion evidence:

- focused red test → 7 errors, all due to absent `.prepare_matrix_input`.
- focused green test → 35 expectations pass.
- package-wide source tests → new input suite passes; expected M0 classifier contract remains red (1 failure + 4 errors).
- `R CMD build . --no-manual --no-build-vignettes` → pass.
- source-tarball `R CMD check --no-manual --no-build-vignettes` → expected 1 ERROR (5 M2-M3 red contract tests), baseline licence WARNING, and baseline unused-import + undefined-prototype-symbol NOTEs; installed input suite reports 35 passes.
- `git diff --check` → pass.

Exact next task after M1a → M1b: add focused red adapter/result tests, then implement `SummarizedExperiment` assay/group extraction, unified input routing, representation-preserving reconstruction, the stable state vocabulary, and the result skeleton without rescue/classification logic.

Blockers → none.

### 2026-07-14 - M1b SummarizedExperiment adapter + result skeleton

Scope → complete the shared input boundary and output scaffolding without implementing rescue or scientific classification.

Implementation:

- unified matrix/`SummarizedExperiment` dispatch; one selected ordinary numeric assay and one explicit `colData` condition column feed the matrix-core contract;
- implicit selection only for a single assay; named selection required and exact for multiple assays;
- representation reconstruction replaces only the selected assay after ordered row filtering, preserving other assays, `rowData`, `colData`, metadata, class, and sample order;
- stable `complete`/`MNAR`/`MAR`/`insufficient` vocabulary plus typed empty classification, feature-status, condition-group, cutoff, profile, and seed-log components;
- public signature now matches the M0 decision and fails explicitly at the still-unimplemented classification boundary rather than executing prototype logic;
- obsolete `plot_detect_custom`, `ggplot_build`, `%||%`, fixed-majority, and global-union code removed.

Verification + completion evidence:

- focused adapter/result red baseline → expected missing dispatch/reconstruction/constructor errors plus old-signature failure;
- focused matrix + adapter + result green suite → 95 expectations pass;
- package-wide source tests → 96 expectations pass; expected M2-M3 classifier contract remains red (4 errors);
- `R CMD build . --no-manual --no-build-vignettes` → pass;
- source-tarball `R CMD check --no-manual --no-build-vignettes` → expected 1 ERROR from the four M2-M3 contract tests, baseline placeholder-licence WARNING, and baseline unused-import NOTE; no code/documentation mismatch or undefined-symbol finding;
- `git diff --check` → pass.

Exact next task after M1b → M2a: add focused red global-absence/condition-minimum/rescue-audit tests, then implement the deterministic rescue core and seed log without classification or reconciliation.

Blockers → none.

### 2026-07-14 - M2a deterministic rescue core + seed audit

Scope → remove global absences; compute immutable pre-seed minima; seed fully missing condition blocks once; retain provenance. Classification/reconciliation remain outside this commit.

Implementation:

- `.seed_missing_conditions()` returns surviving data in original order, full initial feature audit, condition minima, and typed seed log;
- globally absent features → `retained = FALSE`, `drop_reason = all_missing`, excluded before rescue;
- rescue plan sorts conditions/features/candidate sample names, then samples independently with replacement and inserts each condition's pre-seed minimum;
- seed log records feature, condition, selected sample, old/inserted values, and normalized integer seed;
- caller matrix + all originally observed cells remain unchanged; invalid seeds and conditions without finite values fail explicitly;
- local RNG scope and name-stable plan implemented; formal side-effect/permutation regressions remain M2b.

Verification + completion evidence:

- focused red test → 6 errors, all due to absent `.seed_missing_conditions()`;
- focused green test → 35 expectations pass;
- package-wide source tests → 131 expectations pass; expected M3 classifier contract remains red (4 errors);
- `R CMD build . --no-manual --no-build-vignettes` → pass;
- source-tarball `R CMD check --no-manual --no-build-vignettes` → expected 1 ERROR from four M3 contract tests, baseline placeholder-licence WARNING, and baseline unused-import NOTE; installed suite reports 131 passes;
- direct RNG/permutation smoke check → existing and absent `.Random.seed` states preserved; row/column permutations retain seed log + named assignments;
- `git diff --check` → pass.

Exact next task after M2a → M2b: formalize RNG-state restoration (seed initially present/absent), row/column/condition-order named invariants, and original-order output assertions; adjust the rescue core if those adversarial tests expose a gap, then close the M2 gate.

Blockers → none.

### 2026-07-14 - M2b RNG + permutation invariants

Scope → formalize rescue side-effect and order guarantees, adversarially verify the M2a core, and close the M2 gate. Production code required no correction.

Regression coverage:

- caller `.Random.seed` present → exact state restored after rescue;
- caller `.Random.seed` absent → rescue leaves it absent;
- row order, column order, condition block order, and factor level order changes → identical canonical condition minima, seed log, and named seed cells under the same seed;
- shuffled named group vectors → identical rescue result after sample-name alignment;
- returned data + feature audit retain each input's feature/sample order while named reindexing reproduces the baseline result.

Verification + completion evidence:

- first focused adversarial run → 1 test-comparison failure from incidental data-frame row names after named reindexing; comparison normalized those non-contractual row names;
- focused M2 suite → 48 expectations pass;
- package-wide source tests → 144 expectations pass; expected M3 classifier contract remains red (4 errors);
- `R CMD build . --no-manual --no-build-vignettes` → pass;
- source-tarball `R CMD check --no-manual --no-build-vignettes` → expected 1 ERROR from four M3 contract tests, baseline placeholder-licence WARNING, and baseline unused-import NOTE; installed suite reports 144 passes;
- `git diff --check` → pass.

Exact next task after M2b → M3a: add focused red tests for post-seed feature-condition statistics, named manual-cutoff validation, no-missing cutoff behavior, and exact four-state boundaries; implement those pure helpers without cross-condition reconciliation.

Blockers → none.

### 2026-07-14 - M3a statistics + manual states

Scope → calculate post-rescue feature-condition evidence, resolve manual cutoffs, and assign the stable state vocabulary; cross-condition retention + public orchestration remain M3b.

Implementation:

- typed feature-major long statistics preserve input feature order + canonical condition order; counts, missing fraction, arithmetic mean, and seed provenance derive directly from rescued data;
- surviving all-missing blocks fail as a rescue-order invariant violation;
- manual cutoffs require unique exact condition names + finite values, allow partial overrides for later automatic filling, and preserve suspicious out-of-range values with condition-specific diagnostics;
- conditions without missing values force cutoff `NA`, diagnostic source `not_needed`, and `complete` states;
- state assignment applies `complete` first, strict `< cutoff` for MNAR, equality on the MAR side, and per-condition strict majority only to MAR candidates; unresolved incomplete conditions fail explicitly;
- classification output already matches the typed long-result schema, with retention fields intentionally unresolved for M3b.

Verification + completion evidence:

- focused red baseline → 12 errors from absent statistics/cutoff/state helpers;
- focused M3a suite → 60 expectations pass;
- package-wide source tests → 204 expectations pass; expected M3b public contract remains red (4 errors);
- `R CMD build . --no-manual --no-build-vignettes` → pass;
- source-tarball `R CMD check --no-manual --no-build-vignettes` → expected 1 ERROR from four M3b public integration tests, baseline placeholder-licence WARNING, and baseline unused-import NOTE; installed suite reports 204 passes;
- `git diff --check` → pass.

Exact next task after M3a → M3b: add focused red n-condition reconciliation/result tests; implement deterministic `all_missing` → `insufficient:<sorted conditions>` → `MNAR_all_conditions` precedence, retained per-condition groups, matrix/`SummarizedExperiment` reconstruction, and manual-cutoff public orchestration so the normative fixture + full suite pass.

Blockers → none.

### 2026-07-14 - M3b reconciliation + manual public pipeline

Scope → reconcile arbitrary condition counts, populate retained-only condition
groups, and connect the validated matrix core to the public matrix/SE result.
Automatic cutoffs and explicit profiles remain M4-M5 work.

Implementation:

- explicit feature-by-condition state matrix applies deterministic precedence:
  `all_missing` → `insufficient:<sorted labels>` → `MNAR_all_conditions`;
- feature audit + every surviving long-table row receive final retention/reason;
- per-condition `MNAR`, `MAR`, `complete`, and `MAR_or_complete` groups contain
  retained features only and preserve original feature order;
- public classifier now runs preparation → rescue → statistics → manual cutoff
  resolution → state assignment → reconciliation → representation restoration;
- filtered matrix output preserves sample order; `SummarizedExperiment` output
  preserves non-selected assays, row/column metadata, metadata, and class;
- result carries classifications, groups, feature audit, cutoffs/diagnostics,
  seed provenance, aligned groups, and matched call.

Verification + completion evidence:

- focused red baseline → 6 expected errors from absent reconciliation helpers +
  public skeleton stop;
- focused M3b suite → 45 expectations pass, including all drop rules, empty
  retained output, three conditions, all-complete/no-cutoff input, and matrix/SE
  equivalence;
- package-wide source suite → 266 expectations pass;
- `R CMD build . --no-manual --no-build-vignettes` → pass;
- source-tarball `R CMD check --no-manual --no-build-vignettes` → tests pass,
  1 known placeholder-licence WARNING + 1 known unused-import NOTE; no errors,
  code problems, documentation mismatch, or new finding;
- `git diff --check` → pass.

Exact next task after M3b → M4a: add focused red numerical tests for the common
intensity grid, count-weighted missing proportion, repeated/one-class means, and
stored per-condition raw/profile data; implement the pure profile builder and
wire profiles into results without cutoff detection or plotting.

Blockers → none.

### 2026-07-14 - M4a explicit profile calculation + storage

Scope → construct the missing-vs-complete intensity profile directly from
post-rescue evidence and store it for later plotting/cutoff detection. Plotting
and automatic cutoff selection remain M4b-M5 work.

Implementation:

- each condition stores full raw feature statistics plus `has_missing`, a
  512-point common intensity grid, and count-weighted density components;
- Gaussian KDE uses one pooled `nrd0` bandwidth per condition, shared across
  classes, with explicit singleton handling and a three-bandwidth grid extension;
- `missing_proportion = n_missing * f_missing / (n_missing * f_missing +
  n_complete * f_complete)` is calculated from stored numeric components;
- repeated means remain finite; complete-only and missing-only profiles expose
  exact 0/1 proportions and an explicit profile type;
- metadata records class/seed counts, smoothing choices, observed/grid ranges,
  unsupported points, and warnings; cutoff diagnostics reference the same
  metadata;
- matrix and `SummarizedExperiment` public paths produce identical profiles.

Verification + completion evidence:

- focused red baseline → 3 missing-helper errors + 5 empty-profile failures;
- focused profile suite → 49 expectations pass;
- package-wide source suite → 300 expectations pass;
- row/statistics, feature, sample, and condition-order profile invariants → pass;
- `R CMD build . --no-manual --no-build-vignettes` → pass;
- source-tarball `R CMD check --no-manual --no-build-vignettes` → tests pass,
  1 known placeholder-licence WARNING + 1 known unused-import NOTE; no errors,
  code problems, documentation mismatch, or new finding;
- `git diff --check` → pass.

Exact next task after M4a → M4b: add focused red plot tests, then implement
`plot_missingness(result, condition)` exclusively from stored profile data with
condition title, percentage y axis, recorded cutoff line when available, and
uncluttered seed/diagnostic cues; verify a clean installed namespace.

Blockers → none.

### 2026-07-14 - M4b stored-profile plotting

Scope → expose the M4a evidence as a condition-specific plot without deriving
scientific values from rendered objects or re-reading input data.

Implementation:

- exported `plot_missingness(result, condition)` validates the result,
  condition, stored grid/raw evidence, cutoff, and diagnostic schemas;
- filled proportion curve uses the stored 512-point grid directly; fixed
  percentage scale, condition title, and mean-log2-intensity axis make the
  diagnostic explicit;
- finite recorded cutoffs receive a labelled subtitle + vertical line;
  no-missing conditions omit the irrelevant cutoff layer;
- bottom rug ticks mark seeded feature-condition blocks; captions appear only
  for seed provenance or stored profile/cutoff warnings;
- plotting uses namespace-qualified `ggplot2` calls; algorithmic code remains
  independent of `ggplot_build()`.

Verification + completion evidence:

- focused red baseline → 5 expected errors from absent `plot_missingness()`;
- focused plot suite → 29 expectations pass;
- package-wide source suite → 329 expectations pass;
- installed-package smoke → export present; `ggplot_build()` renders;
  1200x750 PNG visual QA confirms title, percentage axis, cutoff, seed tick, and
  caption;
- `R CMD build . --no-manual --no-build-vignettes` → pass;
- source-tarball `R CMD check --no-manual --no-build-vignettes` → tests pass,
  1 known placeholder-licence WARNING + 1 known unused-import NOTE; no error,
  code problem, documentation mismatch, namespace leak, or new finding;
- `git diff --check` → pass.

Exact next task after M4b → M5a: create deterministic clear/noisy/weak/flat/
multi-transition synthetic profile scenarios plus a reproducible validation
harness/report; predeclare cutoff-error, boundary-bias, stability, false-
confidence, and runtime criteria before comparing detector candidates.

Blockers → none.

### 2026-07-14 - M5a frozen automatic-cutoff benchmark

Scope → freeze simulation truth, comparison metrics, and pass/failure gates before
any candidate detector is evaluated; provide the reusable benchmark machinery.

Implementation:

- `dev/cutoff-validation.R` generates sample-level log2 intensities + independent
  cell-level MAR and ramped MNAR masks for 12 deterministic scenarios / 13 target
  profiles, then uses the production rescue, statistics, reconciliation, and
  profile functions rather than a parallel approximation;
- a complete anchor condition preserves fully absent target-condition blocks for
  real one-cell rescue while eliminating accidental global absence;
- frozen cases cover sharp/broad/weak cliffs, noisy tail, 5%/25% MAR, flat/null,
  dominant dual transitions, class imbalance, exact duplicate means, two
  condition-specific boundaries, and 96-feature near-floor data;
- protocol fixes eight simulation seeds, right-boundary targets, scenario-specific
  error/bias/bootstrap tolerances, ≥80 total / ≥12-per-class eligibility floor,
  oracle-relative mechanism/retention metrics, runtime gates, and zero-tolerance
  false confidence on flat data;
- detector contract + benchmark, summary, machine gate, stratified bootstrap, and
  malformed-output handling are reusable for M5b; raw stochastic output remains
  regenerable rather than tracked;
- KDE inputs are radix-sorted before bandwidth/density calculation, removing
  machine-precision row-order drift exposed by the large scenarios; large-profile
  exact-order regression added.

Verification + completion evidence:

- focused red baseline → harness path absent, so the declared `--verify` command
  failed before implementation;
- `Rscript --vanilla dev/cutoff-validation.R --verify` → pass: 12 scenarios, 13
  target profiles, deterministic truth/RNG restoration, ≥12 members in both
  profile classes, exact raw-row/sample/condition-order profiles + rescue log,
  detector-contract/scoring self-checks;
- focused profile suite → 51 expectations pass;
- package-wide source suite → 331 expectations pass;
- `R CMD build . --no-manual --no-build-vignettes` → pass;
- source-tarball `R CMD check --no-manual --no-build-vignettes` → tests pass,
  1 known placeholder-licence WARNING + 1 known unused-import NOTE; no error,
  code problem, documentation mismatch, or new finding;
- `git diff --check` → pass.

Exact next task after M5a → M5b: keep the committed protocol frozen; review
primary method sources, implement at least one segmented/change-point and one
fixed derivative/boundary candidate in the development harness, run the full
eight-seed + bootstrap comparison, append every failed gate and concise result to
`dev/cutoff-validation.md`, then select the simplest passing method or record that
no candidate qualifies.

Blockers → none.

### 2026-07-14 - M5b cutoff candidates + method decision

Scope → compare fixed segmented and derivative families under the frozen M5a
protocol; select a production candidate only if every scientific, failure,
stability, bootstrap, and runtime gate passes. Production integration remains
M5c.

Implementation + research:

- reviewed primary piecewise-regression (Muggeo 2003) and local-polynomial
  smoothing/differentiation (Savitzky-Golay 1964) method sources;
- added clean-environment, base-R `segmented_plateau_v1` and
  `derivative_boundary_v1` candidates with fixed internal evidence/failure rules;
- extended the harness with the declared 200-resample stratified bootstrap,
  per-detector selection gate, deterministic elapsed-free hashes, and a
  `--benchmark` machine entry point;
- kept frozen scenarios, seeds, targets, tolerances, metrics, and gates unchanged;
  raw stochastic rows remain regenerable rather than tracked.

Decision:

- `derivative_boundary_v1` selected → 13/13 eight-seed rows + 13/13 bootstrap
  rows pass; flat false confidence = 0/8 + 0/200; exact repeat/order decisions;
- `segmented_plateau_v1` rejected → eight-seed failures on heavy MAR, sharp, and
  two-transition profiles; bootstrap failures on small + two-transition cases;
  every failed gate + value recorded in `dev/cutoff-validation.md`;
- selected method = cubic Savitzky-Golay derivative, first credible descending
  half-depth lobe, KDE-boundary correction + asymmetric-noise cap, with negative
  likelihood-ratio trend/evidence floors as structured-failure prerequisites.

Verification + completion evidence:

- preimplementation `Rscript --vanilla dev/cutoff-validation.R --benchmark` →
  expected usage failure (M5b path absent);
- `Rscript --vanilla dev/cutoff-validation.R --verify` → pass: frozen M5a
  manifest unchanged;
- final full `Rscript --vanilla dev/cutoff-validation.R --benchmark` → pass in
  111.4 s; selected derivative warmed medians 1-3 ms/profile, max p95 8 ms;
- deterministic hashes → decisions `6ada371e3669e8c4892d187c37fb6ee7`, benchmark
  assessment `6af8fc7e616cc162a3f3f183163eec85`, bootstrap assessment
  `7dbbbd7573f30ffbac2e77e5ebeab56c`;
- package-wide source suite → 331 expectations pass;
- `R CMD build . --no-manual --no-build-vignettes` → pass;
- source-tarball `R CMD check --no-manual --no-build-vignettes` → tests pass,
  1 known placeholder-licence WARNING + 1 known unused-import NOTE; no error,
  code problem, documentation mismatch, or new finding;
- `git diff --check` → pass.

Exact next task after M5b → M5c: promote only `derivative_boundary_v1` into a pure
runtime detector; add focused red tests for evidence floors, one-class/flat
profiles, condition-specific structured errors, finite/interior cutoffs, quality
statistics, algorithm/version metadata, exact order stability, and frozen clear/
stress fixtures; then wire automatic resolution into the internal cutoff path
without yet broadening public documentation work.

Blockers → none.

### 2026-07-14 - M5c production automatic cutoff

Scope → promote the selected M5b derivative boundary into the runtime,
preserve every frozen decision, expose its evidence gates, and connect automatic
plus partial-manual resolution to the public classifier.

Implementation:

- pure `derivative_boundary` version `1` detector preserves the frozen evidence,
  density-support, logistic-trend, cubic local-polynomial, half-depth boundary,
  KDE correction, and peak-offset rules without a public tuning parameter;
- diagnostics expose typed evidence/support/trend/derivative measures and the
  selected transition; quality remains the interpretable vector of predeclared
  gates rather than a synthetic composite score;
- automatic failure returns no finite/endpoint fallback and raises
  `imputefinder_cutoff_unidentifiable` / `imputefinder_cutoff_error` carrying the
  exact condition, reason, diagnostic, and stored profile for manual recovery;
- public `cutoffs = NULL` resolves each incomplete condition automatically;
  partial named cutoffs remain exact manual overrides; complete conditions remain
  `not_needed`;
- automatic decisions are stable under feature/sample/condition order, while
  output feature/sample order continues to follow input order.

Verification + completion evidence:

- focused red baseline → 5 absent-detector errors + unresolved automatic-path
  error + generic public failure;
- focused M5c suite → 74 expectations pass: clear/heavy-MAR/near-floor frozen
  fixtures, both one-class directions, evidence floors, flat/unsupported failure,
  finite interior boundaries, quality/version metadata, condition specificity,
  partial manual override, typed error payload, and exact order invariance;
- production vs selected M5b implementation → byte-identical status/cutoff
  for all 104 frozen eight-seed scenario-condition decisions;
- bootstrap parity → byte-identical status/cutoff/reason for all 2,600 frozen
  resamples, including 200/200 structured flat-profile failures;
- package-wide source suite → 405 expectations pass;
- `R CMD build . --no-manual --no-build-vignettes` → pass;
- source-tarball `R CMD check --no-manual --no-build-vignettes` → tests,
  installed namespace, code, and Rd checks pass; only the known placeholder-
  licence WARNING + obsolete-import NOTE remain for M6/M8;
- `git diff --check` → pass.

Exact next task after M5c → M6a: add focused red result-presentation tests,
then implement concise `print.imputefinder_result()` and
`summary.imputefinder_result()` methods with registered installed-namespace S3
dispatch; preserve the existing result contract and stored-profile plot helper.

Blockers → none.

### 2026-07-14 - M6a result presentation

Scope → expose concise result/summary presentation with machine-readable counts
and installed S3 dispatch. Representation/order audit + obsolete dependency
cleanup remain M6b-M6c; the stored-profile plot helper was completed in M4b.

Implementation:

- registered `print.imputefinder_result()`,
  `summary.imputefinder_result()`, and the summary print method;
- structured summary retains the call and exact feature/sample/seed counts,
  all-feature + retained state counts by condition, precedence-ordered drop
  counts, and condition cutoffs/sources;
- result print shows retained MNAR/MAR/complete group counts + cutoff provenance;
  summary print shows all classified states + drop reasons without dumping the
  underlying matrices, profiles, or diagnostics;
- no-cutoff complete conditions render as `not needed`; every print method
  returns its input invisibly.

Verification + completion evidence:

- focused red result suite → expected 4 failures + 2 errors from base list
  printing/default summary and absent structured fields;
- focused green result suite → 40 expectations pass;
- package-wide source suite → 424/424 expectations pass;
- roxygen namespace/help generation → three S3 registrations + one result-method
  help topic; `R CMD check` reports consistent generics/methods and Rd usage;
- clean installed-package smoke → all three `getS3method()` lookups, result
  summary dispatch, and concise print paths pass in a fresh R session;
- `R CMD build . --no-manual --no-build-vignettes` → pass;
- source-tarball `R CMD check --no-manual --no-build-vignettes` → tests,
  install/load, S3, code, and Rd checks pass; only the known placeholder-licence
  WARNING + obsolete `assertthat`/`dplyr`/`tidyr` import NOTE remain for M8/M6c;
- `git diff --check` → pass.

Exact next task after M6a → M6b: adversarially audit public matrix and
`SummarizedExperiment` outputs across manual + automatic paths; add focused
metadata, assay, feature/sample/condition-order, and core-equivalence regressions,
fix any exposed gap, then close the representation/order gate.

Blockers → none.

### 2026-07-14 - M6b representation/order integration audit

Scope → adversarially verify public matrix and `SummarizedExperiment` parity,
representation preservation, and ordering across manual, automatic, and empty
outputs. Obsolete dependency/prototype cleanup remains M6c.

Audit coverage:

- manual normative path permutes features, interleaves conditions/samples, uses
  named groups + reversed cutoff names, and exercises deterministic rescue;
- automatic path permutes 1,200-feature evidence, samples, condition blocks, and
  factor levels while preserving the frozen condition-specific cutoffs;
- both paths select the second of two assays and prove identical classifications,
  groups, feature audit, cutoffs/diagnostics, profiles, seed log, and aligned groups;
- output preserves exact class, assay names/order, non-selected assay values,
  selected seed-modified data, `rowData`, `colData`, object metadata, retained
  feature order, and original sample order while leaving both inputs unchanged;
- all-ineligible `SummarizedExperiment` output retains both zero-row assays,
  samples, column metadata, object metadata, and full audit evidence;
- shared automatic fixtures moved into the helper layer so cutoff and adapter
  tests use one deterministic generator.

Findings:

- production adapter/core required no correction;
- initial empty-output assertion exposed only valid container convention:
  empty base matrix row names = `NULL`; empty `SummarizedExperiment` row names =
  `character(0)`. Retained-feature audit defines the cross-representation
  invariant without conflating those representations.

Verification + completion evidence:

- pre-audit package-wide source suite → 424/424 expectations pass;
- focused M6b suite → 65/65 expectations pass;
- focused automatic cutoff + M6b + adapter suite → pass;
- package-wide source suite → 489/489 expectations pass;
- `R CMD build . --no-manual --no-build-vignettes` → pass;
- source-tarball `R CMD check --no-manual --no-build-vignettes` → tests,
  installation, namespace, code, and Rd checks pass; only the known placeholder-
  licence WARNING + obsolete `assertthat`/`dplyr`/`tidyr` import NOTE remain;
- `git diff --check` → pass.

Exact next task after M6b → M6c: prove no tracked code, tests, docs, or public
usage needs the prototype-only `assertthat`, `dplyr`, or `tidyr` dependencies;
remove obsolete imports and any remaining dead prototype artefacts, run the
installed API/namespace and source-tarball checks, then close M6.

Blockers → none.

### 2026-07-14 - M6c dependency + prototype cleanup

Scope → remove unused runtime declarations and dangerous scaffold generators,
confirm the compatibility decision, and close the installed-public-API gate.

Audit + implementation:

- pre-change source-tarball check reproduced the expected unused-import NOTE for
  `assertthat`, `dplyr`, and `tidyr`; exact tracked search found no runtime, test,
  documentation, or public use beyond `DESCRIPTION` and historical ledgers;
- removed those three imports; retained `methods`, `SummarizedExperiment`, and
  `ggplot2`, which are all namespace-qualified runtime dependencies;
- deleted four unreferenced `biocthis` bootstrap scripts that could recreate
  placeholder metadata, obsolete workflows, branches, tests, and documentation;
  retained the reproducible M5 cutoff-validation script/report;
- no deprecation shim added: M0 found no releases or public dependents, and the
  obsolete arguments never formed a working released API.

Verification + completion evidence:

- package-wide source suite → 489/489 expectations pass;
- post-change tracked search → removed dependencies appear only in normative or
  historical plan/ledger text;
- `R CMD build . --no-manual --no-build-vignettes` → pass;
- source-tarball `R CMD check --no-manual --no-build-vignettes` → tests,
  installation, dependency analysis, namespace, S3, code, and Rd checks pass;
  the unused-import NOTE is gone and only the known M8 placeholder-licence
  WARNING remains;
- isolated-library vanilla-session smoke → installed package path, dependency
  declarations, matrix classification, plotting, result/summary output, and all
  three S3 registrations pass;
- `git diff --check` → pass.

Exact next task after M6c → M7a: freeze a compact independent scientific
simulation protocol for manual + automatic cutoffs, including uniform 5%/25%
MAR, intensity-dependent MNAR, condition-specific on/off truth, group sizes
4/8/20, a cliff near 12, cutoff sweeps 8-14, seeds/permutations, and predeclared
classification/retention/cutoff metrics; keep long runs in `dev/` and establish
the smallest deterministic routine-test subset before judging results.

Blockers → none.

### 2026-07-14 - M7a frozen scientific-validation protocol

Scope → freeze an independent causal simulation, oracle, metrics, acceptance
gates, and routine/long-run split before inspecting any classifier result.

Protocol:

- sample-level two-condition generator retains complete intensities plus separate
  cell-level uniform MAR and ramped intensity-dependent MNAR masks; disjoint
  structural A-off/B-off sets force full condition dropout and exercise rescue;
- causal oracle treats any realized MNAR mask as MNAR, applies strict majority
  only to MAR-only blocks, then independently reconciles global absence,
  insufficiency, and all-MNAR state;
- 13 long scenarios pair 5%/25% MAR with group sizes 4/8/20 and a translated
  cutoff sweep 8-14; eight simulation seeds, three rescue seeds, and four exact
  name-aligned permutation audits are fixed;
- the 20-sample/25%-MAR case is explicitly an evidence-sensitivity audit because
  `.75^20` leaves too few complete blocks for the production detector; it must
  fail structurally below evidence floors rather than fabricate a cutoff;
- mechanism/state/retention/on-off/cutoff metrics + quantitative manual/automatic
  gates are predeclared in `dev/scientific-validation.md`; two compact cases are
  fixed for routine tests only after the long assessment;
- generator/scorer use local base-R RNG and no production helper; M7a verification
  examines generator/oracle invariants only, so this commit makes no performance
  claim.

Verification + completion evidence:

- initial red command → expected missing `dev/scientific-validation.R` failure;
- `Rscript --vanilla dev/scientific-validation.R --verify` → pass: 13 long + 2
  routine scenarios, exact RNG restoration/repetition, causal-mask union,
  structural-off truth, perfect-oracle scorer, translated 8/14 masks/states, and
  valid B-first feature/sample permutation; manifest MD5
  `4011e381bba2d0d747e91d277a45de5e`;
- package-wide source suite → 489/489 expectations pass;
- `R CMD build . --no-manual --no-build-vignettes` → pass;
- source-tarball `R CMD check --no-manual --no-build-vignettes` → tests,
  installation, namespace, code, and Rd checks pass; only the known M8
  placeholder-licence WARNING remains;
- `git diff --check` → pass.

Exact next task after M7a → M7b: preserve protocol version 1 unchanged; add the
public-path manual/automatic benchmark runner, invariant + alternative-rescue-
seed + exact permutation assessors, and deterministic summaries; run all 13 × 8
simulations, append every outcome/failed gate to the report, and correct the
runtime method only when evidence identifies a real implementation defect.

Blockers → none.

### 2026-07-15 - M7b full scientific benchmark + assessment

Scope → execute protocol v1 unchanged through the installed public classifier,
audit scientific performance + all invariants, and publish deterministic
scenario/condition evidence without adding long stochastic work to package
tests.

Implementation + assessment:

- `--benchmark` now runs all 13 scenarios × 8 simulation seeds through paired
  true-manual and automatic paths; a partial manual override isolates the other
  condition when a joint automatic call raises its structured first failure;
- every public invocation audits input/RNG immutability, original observed
  values, inserted cells, seed-log provenance, pre-seed minima, output axes, and
  aligned groups; structured failures audit target condition + profile evidence;
- alternative rescue seeds 7/29 compare statistics/states/retention/profiles/
  cutoffs with seed 1; translated sweeps compare masks, manual/automatic states,
  retention, status, and cutoff error; four B-first/factor-level permutations
  compare complete canonical named results including chosen seed samples;
- deterministic summaries + gate evaluator preserve raw per-seed regeneration,
  report every scenario/condition outcome, list failures, hash the result and
  harness, and fail the command on any gate;
- first pass found only a harness type mismatch: permutation truth held group
  labels as a factor while the public seed log correctly returned characters;
  normalizing that audit comparison changed provenance from 656/664 to 664/664
  without changing classifier output or a scientific metric;
- no runtime-method defect found; production R files remain unchanged.

Scientific outcome:

- all 270 gates pass; manual worst q10 MNAR/MAR/macro/retention F1 =
  `.9659/.9473/.9140/.9863` at 5% MAR and
  `.9316/.9606/.9620/.9564` at 25% MAR;
- all required/stress automatic cases succeed 8/8 per condition; maximum median
  absolute cutoff error = `.14755`, maximum q90 = `.33788`, signed medians span
  `[-.13456,.14755]`, worst q10 MNAR/retention delta =
  `-.04793/-.01076`, and all 192 estimates are interior;
- all 16 `group_n20_mar25` condition profiles contain only 1-7 complete blocks,
  remain below the 12-block evidence floor, and fail with the correct structured
  condition-specific error rather than a fabricated cutoff;
- exact audits → 664/664 public-call provenance, 208/208 manual + 416/416
  automatic rescue-seed invariants, 48/48 sweep masks/manual states, 96/96 each
  automatic sweep status/state/shift comparison, and 4/4 manual + 8/8 automatic
  permutations;
- protocol MD5 remains `4011e381bba2d0d747e91d277a45de5e`; deterministic
  result MD5 = `72142402979a25eccead24b5af4bba11`; harness MD5 =
  `e2e755fa36f69055e8acb23ce638816e`.

Verification + completion evidence:

- red baseline → `--benchmark` rejected because M7b execution was absent;
- fresh project-local `R CMD INSTALL --preclean .` → pass;
- fresh-installed full benchmark → 270/270 gates pass in 24.390 s; repeated
  deterministic result + harness hashes exact;
- package-wide source suite → 489/489 expectations pass;
- `R CMD build . --no-manual --no-build-vignettes` → pass;
- source-tarball `R CMD check --no-manual --no-build-vignettes` → tests,
  installation, namespace, code, Rd, and examples checks pass; only the known M8
  placeholder-licence WARNING remains;
- `git diff --check` → pass.

Exact next task after M7b → M7c: encode the frozen `manual_on_off` and
`automatic_cliff` routine cases in ordinary package tests with their predeclared
scientific thresholds + exact named permutation assertion; keep the long runner
in `dev/`, rerun broad checks, mark PLAN M7 coverage complete, and close M7.

Blockers → none.

### 2026-07-15 - M7c routine scientific regressions

Scope → move the two predeclared protocol-v1 routine cases into ordinary package
tests, preserve the long benchmark in `dev/`, and close the scientific-regression
gate without changing the runtime classifier.

Implementation + evidence:

- test-only protocol mirror regenerates `manual_on_off` (320 features, 4 samples
  per condition, 25% uniform MAR) and `automatic_cliff` (800 features, 8 samples
  per condition, 5% MAR) at simulation seed 104677; causal truth remains derived
  from independent MAR/MNAR masks + strict-majority/on-off reconciliation;
- mirror matrices, aligned groups, oracle states, and feature retention are exact
  against the frozen M7a generator; compact scores are exact against the M7b
  scorer while remaining self-contained when `dev/` is excluded from a source
  package;
- manual routine MNAR/MAR/macro/retention F1 =
  `.94382/.96736/.96805/.95094`; eligible on/off retention + off-block rescue =
  `1/1`, clearing all frozen 25%-MAR thresholds;
- automatic routine cutoffs A/B = `11.97274/11.88546`; automatic
  MNAR/MAR/macro/retention F1 = `.97770/.96703/.98158/.98983`, paired
  manual deltas = `-.00385/-.00082` for MNAR/retention F1, and both on/off gates
  = `1`; all frozen routine gates clear;
- even/odd feature order + B-first/reversed samples + reversed factor levels
  preserve the complete canonical named automatic result exactly, including
  filtered/seeded data, cutoffs, diagnostics, profiles, classifications,
  retention, groups, and seed assignments;
- production R files remain unchanged; the long 13-scenario × 8-seed benchmark,
  alternative rescue seeds, sweeps, and four-case permutation audit remain only
  in `dev/`.

Verification + completion evidence:

- focused red → 2 expected missing-helper errors;
- focused green → 16/16 expectations pass;
- package-wide source suite → 505/505 expectations pass;
- `Rscript --vanilla dev/scientific-validation.R --verify` → protocol v1 + MD5
  `4011e381bba2d0d747e91d277a45de5e` unchanged;
- `R CMD build . --no-manual --no-build-vignettes` → pass;
- source-tarball `R CMD check --no-manual --no-build-vignettes` → package tests,
  installation, namespace, code, Rd, and dependency checks pass; only the known
  M8 placeholder-licence WARNING remains;
- `git diff --check` → pass.

Exact next task after M7c → M8a: verify current official R/Bioconductor metadata
guidance, then replace `DESCRIPTION` identity/licence placeholders, audit runtime
and suggested dependencies + minimum platform fields + `biocViews`, and add the
package-level documentation needed to make that metadata checkable. Keep the
README/vignette workflow rewrite as the next cohesive documentation slice.

Blockers → none.

### 2026-07-15 - M8a Bioconductor metadata + package overview

Scope → replace scaffold identity/licence metadata, align the declared platform
with the current Bioconductor release, and add an installed package overview.
README workflow prose, runnable examples, vignette, and NEWS remain M8b-M8c.

Guidance + decisions:

- current official DESCRIPTION guidance inspected 2026-07-15 requires a
  three-sentence functional description, one active `Authors@R` maintainer,
  Bioconductor submission version `0.99.0`, and at least two valid same-branch
  `biocViews`; Bioconductor 3.23 is paired with R 4.6;
- `DESCRIPTION` now identifies Emir Turkes via the repository's durable contact,
  uses a behavior-bounded title/description, declares `GPL (>= 3)` to match the
  tracked version-3-or-later source headers + `COPYING`, and links the source +
  issue tracker;
- Software views = `MassSpectrometry`, `Proteomics`, `Preprocessing`; `Software`
  records package type. `Depends: R (>= 4.6.0)` targets the current platform;
- runtime audit retains namespace-used `SummarizedExperiment`, `ggplot2`, and
  `methods`; `Suggests` currently contains only test-used `testthat`. Vignette
  dependencies enter with the vignette rather than becoming speculative
  metadata;
- new roxygen package help states the unnormalised-log2 input contract,
  heuristic likely-MNAR/MAR interpretation, audit evidence, and boundary before
  external normalisation/imputation.

Primary guidance:

- `https://contributions.bioconductor.org/description.html`
- `https://contributions.bioconductor.org/license.html`
- `https://bioconductor.org/about/release-announcements/`

Verification + completion evidence:

- pre-change source-tarball check → expected single WARNING for the scaffold
  non-standard licence string;
- package-wide source suite → 505/505 expectations pass;
- roxygen generation + internal DESCRIPTION validation → pass;
- source build + project-local install → pass; installed version, GPL metadata,
  maintainer, and `imputefinder-package` help smoke → pass;
- final source-tarball `R CMD check --no-manual --no-build-vignettes` → status
  OK, including installed tests, dependency analysis, namespace, code, and Rd;
  examples/vignette remain intentionally absent until M8b-M8c;
- `git diff --check` → pass.

Exact next task after M8a → M8b: increment the Bioconductor development patch
version; replace the README scaffold with rationale, Bioconductor/development
installation, an offline normative-fixture walkthrough, result interpretation,
automatic/manual cutoff guidance, seed provenance, downstream pipeline order,
and evidence-bounded limitations; add compact runnable examples to the public
Rd pages and verify installed examples. Keep the full vignette + NEWS for M8c.

Blockers → none.

### 2026-07-15 - M8b README + runnable public examples

Scope → make the complete public workflow understandable without thesis access,
add executable help examples, and advance the Bioconductor development version.
The full vignette + NEWS remain M8c.

Guidance + implementation:

- current official Bioconductor README/install guidance confirms
  `BiocManager::install()` for repository packages and `owner/repository`
  GitHub development installs; README identifies the latter as current and the
  former as available only after acceptance;
- version advances `0.99.0` → `0.99.1` as the next Bioconductor development
  patch;
- README now explains condition-specific MAR/MNAR rationale, exact four-state +
  reconciliation rules, the normative on/off fixture, input assumptions,
  one-cell condition-minimum rescue, local seed provenance, result components,
  manual/automatic cutoff evidence + failure, matrix/SE routing, downstream
  normalisation/imputation order, limitations, and evidence boundaries;
- source-tarball-safe canonical links expose the tracked M5/M7 reports + plan;
- compact manual-cutoff examples on classifier, result, and plot help pages all
  exercise rescue + condition groups offline; roxygen-generated Rd stays synced;
- direct plot-example verification produces `Rplots.pdf`; the exact generated
  filename is now ignored and the working tree remains artefact-free;
- no runtime or suggested dependency was added.

Primary guidance:

- `https://contributions.bioconductor.org/readme.html`
- `https://bioconductor.org/install/`
- `https://bioconductor.github.io/BiocManager/reference/install.html`

Verification + completion evidence:

- pre-change acceptance audit → version `0.99.0`, six required README concepts
  absent, and no public Rd example blocks;
- project-local install + direct installed `example()` calls for
  `classify_missingness`, `imputefinder_result`, and `plot_missingness` → pass;
- package-wide source suite → 505/505 expectations pass;
- Pandoc → standalone HTML → Chromium PDF → six rendered pages inspected;
  headings, tables, code, links, and page flow remain legible with no overflow;
- `R CMD build . --no-manual --no-build-vignettes` → pass;
- source-tarball `R CMD check --no-manual --no-build-vignettes` → status OK,
  including examples, installed tests, dependencies, namespace, code, and Rd;
- `git diff --check` → pass.

Exact next task after M8b → M8c: consult current official Bioconductor vignette
and NEWS guidance; add one installed vignette covering manual + automatic cutoffs,
the on/off fixture, stored plotting evidence, filtering, seed provenance, and
condition-specific downstream branches; add `NEWS.md`, finish public/package
help cross-links, run the vignette/examples from a clean source tarball, and
close the M8 documentation gate.

Blockers → none.

### 2026-07-15 - M8c vignette + NEWS + documentation gate

Scope → supply the installed task-oriented vignette + release history, make the
long-form workflow discoverable from public help, and close every M8/package-
vignette gate. CI and Bioconductor submission checks remain M9.

Guidance + decisions:

- current official Bioconductor documentation guidance requires at least one
  task-oriented Rmd/Rnw vignette with executed non-trivial package code,
  introduction + Bioconductor motivation, installation instructions,
  references, and terminal `sessionInfo()`; only installation is justified as
  non-evaluated here;
- current guidance recommends `BiocStyle`; `DESCRIPTION` therefore declares
  `BiocStyle`, `knitr`, and `rmarkdown` in `Suggests` plus
  `VignetteBuilder: knitr`; no runtime dependency enters;
- current NEWS guidance accepts one top-level Markdown file; `NEWS.md` records
  the initial development behavior under version `0.99.1`;
- primary guidance:
  `https://contributions.bioconductor.org/docs.html`,
  `https://contributions.bioconductor.org/news.html`, and
  `https://bioconductor.org/help/package-vignettes/`.

Implementation:

- installed `BiocStyle` vignette executes 26/26 chunks offline except the
  explicitly unevaluated package-installation chunk;
- normative manual fixture demonstrates global absence, strict majority,
  all-MNAR/insufficient exclusion, a feature MNAR in A + complete in B,
  one-cell condition-minimum rescue, seed provenance, result components, and
  condition-specific downstream branches without choosing an imputer;
- deterministic 200-feature public-API example demonstrates two independent
  automatic cutoffs + quality metadata; manual + automatic figures read only
  stored profiles and recorded cutoffs;
- text covers input/pipeline boundaries, arithmetic-mean + majority rationale,
  count-weighted stacked-density meaning, automatic structured failure,
  reconciliation, limitations, related Bioconductor tools, and matrix/SE
  parity; the profile formula remains readable without browser-side MathJax;
- README + all public/package Rd pages point to `vignette("imputefinder")`;
  roxygen source + generated Rd remain synchronized;
- optional BiocStyle image cropping is disabled rather than adding `magick`.

Verification + completion evidence:

- pre-change documentation assertion → expected failure at absent
  `vignettes/`; `NEWS.md`, vignette metadata, and builder dependencies were
  absent too;
- public automatic-example probe across 200/300/400/600 features → both
  conditions succeed; selected 200-feature fixture cutoffs =
  `A 10.706691`, `B 13.978929`;
- direct `rmarkdown` build → 26/26 chunks pass without warnings; Chromium PDF
  review covered title/TOC, rule table, audit output, complete result table,
  both profile figures, downstream split, references, and session information;
  no horizontal overflow; review findings (optional crop dependency, reserved
  References heading, sparse-fixture caveat, online math rendering) corrected;
- package-wide source suite → pass (505 existing expectations; runtime
  unchanged);
- vignette-bearing source build → pass; tarball contains source + installed
  R/Rmd/HTML vignette and `NEWS.md`, while excluded development/agent files stay
  absent;
- clean tarball install + installed vignette R script + installed examples for
  `classify_missingness`, `imputefinder_result`, and `plot_missingness` → pass;
- final `R CMD check --no-manual imputefinder_0.99.1.tar.gz` → status OK,
  including examples, tests, installed vignette files/dependencies, and
  vignette re-build;
- `git diff --check` → pass.

Exact next task after M8c → M9a: consult live official Bioconductor CI guidance
and maintained action releases; add a minimal matched R/Bioconductor workflow
that installs dependencies and runs the package test/build/check path; validate
workflow syntax, permissions, caching, and command parity locally. Preserve
clean `--as-cran` + `BiocCheck` execution/fixes as M9b-M9c.

Blockers → none.

### 2026-07-15 - M9a Bioconductor-devel CI

Scope → add one maintained, least-privilege package-check workflow matched to
the live Bioconductor development environment. Submission-grade `--as-cran` and
`BiocCheck` execution remain M9b-M9c.

Guidance + decisions:

- current official guidance says contributed packages develop against
  Bioconductor devel; in the current release window that is Bioconductor 3.24
  on R 4.6;
- the official `bioconductor/bioconductor_docker:devel` image tracks that pair,
  supplies Bioconductor build-system dependencies, and approximates the Linux
  build machine; a floating devel tag is intentional for continuously current
  CI rather than archival reproducibility;
- maintained `r-lib/actions` major tag remains `v2`; its check action builds a
  source tarball before checking it and runs package tests through `R CMD check`;
- primary guidance:
  `https://contributions.bioconductor.org/general.html`,
  `https://contributions.bioconductor.org/use-devel.html`,
  `https://bioconductor.org/help/docker/`, and
  `https://github.com/r-lib/actions`.

Implementation:

- `.github/workflows/R-CMD-check.yaml` runs on main pushes, pull requests, and
  manual dispatch with read-only repository permission + duplicate-run
  cancellation;
- one Ubuntu job runs inside the official devel container, exports its live
  `BIOCONDUCTOR_VERSION` as `R_BIOC_VERSION`, then verifies both
  `BiocManager::version()` and installed `BiocVersion`; this matters because
  `pak` otherwise selects the released Bioconductor branch when one R version
  supports both release and devel;
- `setup-r-dependencies@v2` installs all declared/check dependencies plus
  `rcmdcheck` with its maintained package cache; forced Suggests keep vignette
  and example dependencies mandatory;
- `check-r-package@v2` builds, tests, and checks the source package with warnings
  fatal; `--as-cran` remains deliberately deferred until its clean local M9b
  run;
- `.github` is explicitly excluded from source packages via `.Rbuildignore`.

Verification + completion evidence:

- red baseline → expected absence of `.github/workflows/R-CMD-check.yaml`;
- `actionlint` 1.7.12 + independent YAML parse → pass;
- live refs resolve for `actions/checkout@v6` and all `r-lib/actions@v2`
  components; official devel container tag resolved to its current multi-arch
  manifest;
- package-wide source suite → pass (505 existing expectations);
- `R CMD build . --no-manual` → pass with vignette build; 486,030-byte source
  tarball excludes `.github`, `.agent`, `AGENTS.md`, and `PLAN.md`;
- source-tarball `R CMD check --no-manual` → status OK, including examples,
  installed tests, dependency analysis, and vignette rebuild;
- no local container engine is installed; hosted execution necessarily begins
  after the maintainer pushes, while workflow syntax, references, environment
  contract, and exact package commands are locally validated.

Exact next task after M9a → M9b: provision the project-local library against
Bioconductor 3.24, install current `BiocCheck`, then run the unit suite, a clean
source build, `R CMD check --as-cran`, `BiocCheckGitClone()`, and
`BiocCheck(..., new-package = TRUE)` under current official check guidance.
Promote `--as-cran` into CI only after the local gate is clean; carry every
hardening finding into M9c rather than suppressing it.

Blockers → none.

### 2026-07-15 - M9b submission checks + clean `--as-cran`

Scope → execute the current official new-package check sequence on
Bioconductor devel, make the full source-tarball CRAN check locally executable,
and promote its clean mode into CI. Repository-controlled `BiocCheck` findings
remain the cohesive M9c hardening scope.

Guidance + environment:

- current official guidance requires Bioconductor devel plus `R CMD build`,
  `R CMD check`, `BiocCheckGitClone()`, and
  `BiocCheck(..., \`new-package\` = TRUE)`; errors and warnings must be cleared
  and notes addressed or justified;
- upgraded the ignored project library from Bioconductor 3.23 to 3.24 on R
  4.6 and installed matching `BiocCheck` 1.49.29; `BiocManager::valid()` reports
  a valid installation;
- installed self-contained TinyTeX + required Courier/MakeIndex components and
  HTML Tidy under `.agent/`; neither host state nor user PATH was changed;
- primary guidance:
  `https://contributions.bioconductor.org/general.html`,
  `https://bioconductor.org/packages/devel/bioc/html/BiocCheck.html`, and
  `https://yihui.org/tinytex/`.

Verification + findings:

- package-wide source suite on Bioconductor 3.24 → pass (505 expectations);
- vignette-bearing `R CMD build .` → pass;
- source-tarball `R CMD check --as-cran` with PDF + HTML manuals enabled →
  0 errors, 0 warnings, 1 unavoidable NOTE (`New submission`); examples, tests,
  vignettes, manual PDF, and HTML validation pass;
- `BiocCheckGitClone()` on a clean local clone → 0 errors, 0 warnings,
  0 notes;
- `BiocCheck(..., new-package = TRUE)` on the source tarball → 1 error,
  1 warning, 9 notes:
  - external error: maintainer Support Site profile does not watch the
    `imputefinder` tag;
  - code warning: one locally scoped `set.seed()` call; implementation restores
    the caller RNG exactly, but the static check still requires resolution;
  - actionable review notes: consider `Classification` biocView; simplify six
    condition-message constructions + one redundant signal; review two
    warning suppressions, 12 functions over 50 lines, and formatting findings;
  - metadata notes requiring maintainer facts/actions rather than invention:
    ORCID, optional funder role, optional publication-backed `CITATION`, and
    unverifiable bioc-devel subscription;
- CI now runs `R CMD check --as-cran` with warnings fatal;
- generated check/tool directories are ignored uniformly via `.agent/*/` and
  `*.BiocCheck/`.

Exact next task after M9b → M9c: test-drive removal of the `set.seed()`
warning while preserving uniform rescue selection, exact local-RNG behavior,
and all name-order invariants; resolve every repository-controlled BiocCheck
finding, audit licence/line endings/generated docs/file sizes, then rerun the
clean clone + tarball checks. Add BiocCheck to CI only after that gate is clean.

External blocker → maintainer must add `imputefinder` to Watched Tags at
`https://support.bioconductor.org/accounts/edit/profile`; the full unskipped
new-package check cannot report zero errors until the account state changes.

### 2026-07-15 - M9c repository-controlled hardening

Scope → eliminate every repository-controlled `BiocCheck` diagnostic, audit
submission artefacts, and strengthen RNG + automatic-cutoff failure guarantees.
Full M9c closure remains gated only by Support Site propagation and the ensuing
CI promotion.

Implementation:

- replaced direct `set.seed()` handling with `withr` local scopes pinned to
  Mersenne-Twister/Inversion/Rejection; nested seed + RNG-kind preservation now
  restores even a non-default kind paired with an intentionally absent
  `.Random.seed`;
- added adversarial regressions for caller RNG-kind independence and exact
  absent/present RNG state restoration;
- captures `glm.fit()` numerical warnings in cutoff diagnostics and rejects
  warned trend fits instead of silently accepting them; removed the unnecessary
  bandwidth warning suppression;
- added the suggested `Classification` biocView, simplified condition signals,
  removed the redundant signal false positive, and refactored every function to
  at most 50 source lines;
- reformatted tracked R/vignette sources to the Bioconductor checks, migrated
  generated docs to roxygen2 8, and retained its accurate legacy marker because
  current `BiocCheck` detects generated Rd only through `RoxygenNote`;
- reviewed coverage infrastructure and omitted it: 83 focused tests already
  localise contract failures, while a new service would add no diagnostic value.

Verification + audit evidence:

- package-wide source suite → 83 tests / 513 expectations, 0 failures,
  0 warnings;
- vignette-bearing `R CMD build .` → pass;
- source-tarball `R CMD check --as-cran` on R 4.6.1 / Bioconductor 3.24 →
  0 errors, 0 warnings, 1 expected NOTE (`New submission`); examples, tests,
  vignette rebuild, PDF manual, and HTML validation pass;
- patched clean-clone `BiocCheckGitClone()` → pass with no findings;
- unskipped tarball `BiocCheck(..., new-package = TRUE)` → 1 external error,
  0 warnings, 2 metadata/admin notes; coding practice, warning handling,
  function length, formatting, dependencies, docs, tests, and file checks are
  clean;
- remaining error says the public checker cannot yet see `imputefinder` in the
  maintainer's Watched Tags; maintainer reports having added it during this
  session, so this is now a propagation wait rather than an unperformed action;
- GPL-3 text + `License: GPL (>= 3)` verified; no tracked CR bytes; largest
  tracked file is 86,505 bytes; generated docs report no pending roxygen work;
  source tarball excludes local tools/check artefacts and has no oversized files.

Exact next task → rerun the full unskipped tarball `BiocCheck` after Support Site
state propagates. Once it reports zero errors/warnings, add fatal
`BiocCheckGitClone()` + new-package `BiocCheck()` gates to CI, validate the
workflow, mark M9/M9c complete, and begin M10. Preserve the optional ORCID,
funder, and publication-backed CITATION omissions unless maintainer facts are
provided.

External blocker → Support Site Watched Tags propagation; checker still reports
the tag absent as of the final local run despite the maintainer's completed
profile update.

### 2026-07-15 - M10a public guarantees + boundary cases

Scope → audit every README, Rd, and vignette guarantee against implementation
and tests; re-run all named release-boundary cases. M9's Support Site lookup
remains external and does not block independent M10 work.

Findings + corrections:

- pre-reconciliation rescue provenance intentionally includes anchors on
  features later dropped, while filtered `data` contains anchors only for
  retained features; README/vignette now state both halves and a public
  regression locks the audit behavior;
- complete-only conditions need no cutoff, so vignette automatic-failure text
  now distinguishes that case from missing-only/flat/sparse/ambiguous failures;
- input-order reproducibility is scoped to a fixed named dataset + seed, and
  globally absent blocks are excluded from the state-table guarantee;
- one- and two-sample condition tests lock the normative post-rescue result:
  a seeded one-sample block becomes complete, while one of two observations is
  MNAR below the cutoff and insufficient on its MAR side;
- public-path regressions now reject `NaN`, `Inf`, `-Inf`, and conditions with
  no finite pre-seed minimum;
- first rebuilt tarball exposed a root `imputefinder.BiocCheck/` leak despite a
  stale M9 audit claim; `.Rbuildignore` now excludes all `.BiocCheck` outputs,
  and the rebuilt tarball contains no local check/tool artefact.

Claim-family evidence → input/routing (`test-input`, `test-summarizedexperiment`);
rescue/provenance (`test-seeding`, `test-reconciliation`); four states +
reconciliation (`test-classification*`, `test-reconciliation`); stored profile,
automatic/manual cutoffs, structured failure, and plotting (`test-profile`,
`test-cutoff-*`, `test-plot`); result/representation/order (`test-result`,
`test-representation-order`); evidence bounds (`dev/*-validation.md` + routine
scientific regressions). No other implementation/claim mismatch found.

Verification + completion evidence:

- package-wide source suite → 86 tests / 528 expectations, 0 failures,
  warnings, errors, or skips;
- vignette-bearing source build → pass;
- source-tarball `R CMD check --as-cran` → status OK, including examples,
  installed tests, vignette rebuild, PDF manual, and HTML manual;
- tarball content audit → no `.agent`, Git, `Rcheck`, or `BiocCheck` artefacts;
- `BiocCheck(..., new-package = TRUE)` → 0 repository-controlled errors or
  warnings; known Support Site propagation error + optional/admin notes only;
- `git diff --check` → pass.

Exact next task → M10b: adversarially snapshot caller RNG state/kind, options,
graphics devices/parameters, and working directory around the complete public
manual + automatic paths; prove named row/column/condition-order invariance and
observed-cell integrity at the final result boundary, correcting any gap before
closing those three M10 checks.

External blocker → Support Site Watched Tags propagation remains the sole M9
error; it does not block M10b-M10c.

### 2026-07-15 - M10b side effects + invariants

Scope → verify the complete public manual/automatic paths at the process-state
and final-cell boundaries, then close release-order invariance with canonical
full-result comparisons.

Adversarial coverage:

- exact snapshots around both classifier paths + stored plot construction cover
  caller RNG state and kind, every R option, working directory, active graphics
  device/list, and writable graphics parameters under a non-default RNG kind;
- one shared RNG-restoration helper now loads with the test fixtures rather than
  depending on test-file execution order;
- the frozen automatic scientific fixture is permuted simultaneously by feature
  order, sample order, condition-block order, named group order, factor-level
  order, and manual-cutoff name order; canonical manual and automatic results
  match their baselines across data, states, groups, status, cutoffs,
  diagnostics, profiles, seed log, and sample-condition mapping;
- final manual + automatic matrices are compared cell-by-cell with retained
  input rows: every originally observed value is byte-identical, every changed
  cell was originally `NA`, and changed coordinates/values equal the retained
  subset of `seed_log` exactly;
- no production defect surfaced.

Verification + completion evidence:

- focused RNG/invariance/side-effect suite → pass;
- package-wide source suite → 88 tests / 542 expectations, 0 failures,
  warnings, errors, or skips;
- vignette-bearing source build + source-tarball `R CMD check --as-cran` →
  status OK, including installed tests and vignette/manual rebuilds;
- new-package `BiocCheck` → 0 repository-controlled errors/warnings; only the
  known Support Site propagation error + optional/admin notes;
- `git diff --check` → pass.

Exact next task → M10c: benchmark a representative 10,000 x 50 matrix with
stage-level time/allocation evidence, remove dead development artefacts + stale
comments, resolve the Bioconductor submission-version requirement from current
official guidance, and close the release-candidate gate if clean-check evidence
supports it.

External blocker → Support Site Watched Tags propagation remains the sole M9
error; it does not block M10c.

### 2026-07-15 - M10c performance + release cleanup

Scope → benchmark the full public paths at representative release scale, audit
tracked artefact/comment provenance, and resolve the first release version from
current Bioconductor policy.

Version + cleanup decisions:

- current official policy requires new submissions to use `0.99.0`, increment
  `z` during development, and lets Bioconductor promote `0.99.z` to `1.0.0` at
  release; the release candidate is therefore `0.99.2`, not `0.1.0`;
- primary guidance:
  `https://contributions.bioconductor.org/versionnum.html`,
  `https://contributions.bioconductor.org/description.html`, and
  `https://contributions.bioconductor.org/git-version-control.html`;
- NEWS now separates dated 0.99.2/0.99.1 entries and renders through installed
  `utils::news()`; its first multi-version render exposed and corrected the
  previously latent undated-heading failure;
- tracked provenance audit retained all three development harness families:
  M5 candidate comparator + cutoff protocol, M7 scientific protocol, and M10
  performance gate all remain executable evidence and stay outside the source
  tarball; no dead tracked artefact was found;
- stale present-tense prototype wording in `PLAN.md` and the historical M7
  package-warning sentence now identify their exact baseline/milestone scope;
  copyright provenance, RStudio package config, canonical plan/agent ledger,
  frozen rejected comparator, and historical benchmark versions remain factual.

Performance evidence:

- `dev/performance-validation.R` SHA-256
  `74c53ccaf0f0c400693ac76607c575f686012c5b3980ae4d758a5dc4893fc2c2`;
  fixture + six gates were frozen before first execution;
- deterministic 10,000 x 50 matrix = five ten-sample conditions, known cliffs,
  5% background missingness, and 200 structural off blocks;
- manual/automatic median elapsed = `.139/.178` s; cumulative allocation =
  `204.4/299.7` MiB (`45.37x/66.53x` input); largest allocation = `3.81` MiB
  (`.85x`); result = `21.19/21.25` MiB (`4.70x/4.72x`); process high-water
  RSS <=`172.32` MiB; automatic maximum cutoff error = `.282`;
- all predeclared runtime, allocation, result-size, RSS, cutoff, unchanged-input,
  and seed-only-cell-change gates pass; no performance refactor is justified.

Verification + completion evidence:

- M5 cutoff + M7 scientific frozen-protocol verification → pass with unchanged
  manifests;
- package-wide source suite → 88 tests / 542 expectations, 0 failures,
  warnings, errors, or skips;
- installed version + dated NEWS rendering → pass at `0.99.2`;
- vignette-bearing source build → pass; 489,720-byte tarball excludes every
  development/agent/check/CI artefact;
- source-tarball `R CMD check --as-cran` → status OK, including examples,
  installed tests, vignette rebuild, and PDF/HTML manuals;
- new-package `BiocCheck` → 0 repository-controlled errors/warnings; known
  Support Site propagation error + optional/admin notes only;
- no tracked CR bytes; `git diff --check` → pass.

M10 work items are complete. Its release-candidate gate remains open only
because Definition of Done includes the externally gated M9 `BiocCheck`/CI
closure.

Exact next task → when the Support Site API exposes the watched `imputefinder`
tag, rerun unskipped new-package `BiocCheck`; then add fatal Git-clone + tarball
BiocCheck gates to CI, validate them, close M9/M10 + the two remaining Definition
of Done items, and declare the first release candidate complete.

External blocker → Support Site Watched Tags propagation remains the sole M9 /
release-candidate gate error.
