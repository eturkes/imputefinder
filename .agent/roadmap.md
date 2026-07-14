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
- [ ] M7 scientific regression suite
- [ ] M8 documentation + package metadata
- [ ] M9 CI + Bioconductor hardening
- [ ] M10 release-candidate adversarial review

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
