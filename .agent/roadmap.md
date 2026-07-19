# ImputeFinder roadmap

Canonical method + gates → `PLAN.md`. This ledger owns session state, evidence, blockers, and exact next work.

## Current innovation path

- [x] M11 freeze v1 + extension seams
- [x] M12 frozen A-C validation + gate registry
- [ ] M13 design core + confounding-sentinel study
- [ ] M14 robustness certificates
- [ ] M15 detectability contrast engine
- [ ] M15P release confirmation + next-card promotion review

Parked after review: D censor-aware evidence; E workflow stress-test lab; F
distributional handoff; G peptide-protein graph; H1 calibration; H2 feedback; I
delayed backend; J1 standards QC; J2 atlas. Parked work has no active checkbox
or release obligation.

Exact method, dependency graph, falsification gates, and Definition of Done →
`PLAN.md` Sections 2-14. M13c v3 contract hash =
`85d121ebd5ddd589be0f389553a3b7ffa6076a12a8fcbcb476a8969f7d987ae4`.
Exact next implementation task → build the semantic candidate-evidence resolver
and study harness, then green the guarded manifest authorizer and exact public
panel materializer/handoff. Open synthetic candidate replicates `1-32` only
after every resolver + manifest rail passes. Keep development replicates
`33-64`, HarmonizR development bytes, and all confirmation evidence sealed until
the candidate choice locks.

## Completed v1 release path

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
- [x] M9 CI + Bioconductor hardening
  - [x] M9a maintained Bioconductor-aware CI
  - [x] M9b submission checks + clean `--as-cran`
  - [x] M9c `BiocCheck` findings + hardening gate closure
- [x] M10 release-candidate adversarial review
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
- Prototype also contains an undefined `%||%`, self-referential `final_MAR`, fixed `min_non_na = 5`, global MAR/MNAR unions, row mutation, and no rescue path; see the historical completion spec at `30618e1:PLAN.md`, Section 4.

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

Routing + behavior are normative in the historical completion spec at
`30618e1:PLAN.md`, Section 6. Compatibility shim → none: GitHub code searches
for package/function/install reference returned zero; repository has no
releases, tags, forks, or stars; CRAN registry lookup returned 404. Search
performed 2026-07-14. Broken, unreleased prototype arguments provide no
observed compatibility value.

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
- committed-state fresh clone at `69fb2384c1b774d99bcc1089a74ed662a43c7472`
  → clean worktree; `BiocCheckGitClone()` + source build pass;
- fresh-clone source-tarball `R CMD check --as-cran` → status OK, including
  examples, installed tests, vignette rebuild, and PDF/HTML manuals;
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

External blocker → Support profile 78003 resolves to Emir Turkes, but the live
watched-tags API still returned `{"watched_tags": ""}` at 2026-07-15 09:25 UTC.
Propagation remains the sole M9 / release-candidate gate error.

### 2026-07-15 - M9c watched-tag gate recheck

Scope → re-run the sole external release gate against the current committed
0.99.2 source and distinguish delayed propagation from an unsaved Support Site
field. No package implementation change is warranted.

Evidence:

- fresh vignette-bearing `R CMD build .` → pass;
- unskipped `BiocCheck(..., new-package = TRUE)` → 1 error, 0 warnings,
  2 optional/admin notes; the Watched Tags requirement remains the only error;
- the live Support Site endpoint returned `{"watched_tags": ""}` at
  2026-07-15 09:30 UTC, including a cache-busting query;
- current deployed-site source returns `user.profile.watched_tags` directly
  from the account record in this endpoint, and the response carries no cache
  metadata; continued passive waiting is therefore unsupported by the
  available evidence;
- the Support Site had a resolved 2025 defect that erased My Tags/Watched Tags
  during profile edits, but the present failure cannot be corrected from the
  repository or without the maintainer's authenticated account.

Exact next task → maintainer reopens the Support Site profile, enters
`imputefinder` specifically in **Watched Tags**, commits the tag token with
comma/Enter/Space, saves the complete form, and confirms it persists after a
reload. Then re-query the public endpoint; once it returns `imputefinder`, run
the full clean-clone + tarball gates, add fatal `BiocCheckGitClone()` and
new-package `BiocCheck()` CI steps, validate the workflow, and close M9/M10 +
the two remaining Definition of Done items.

External blocker → authenticated Support Site profile correction. M9/M10
remain deliberately open; adding a CI check known to fail would not close the
release gate.

### 2026-07-16 - M9/M10 release-candidate closure

Scope → verify the corrected Support Site state through the same public API
used by `BiocCheck`, execute the complete clean-check sequence, and make both
Bioconductor checks fatal CI gates.

Evidence:

- live watched-tags API returned `genefunnel,imputefinder` at 00:19 UTC;
- clean committed-source clone → clean worktree; `BiocCheckGitClone()` reports
  0 errors, 0 warnings, 0 notes;
- vignette-bearing source build → pass; 489,554-byte tarball;
- clean-clone source-tarball `R CMD check --as-cran` → status OK, including
  examples, installed tests, vignette rebuild, and PDF/HTML manuals; sole note
  is the expected `New submission` incoming-feasibility note;
- new-package `BiocCheck` → 0 errors, 0 warnings, 2 already-reviewed notes:
  optional ORCID and locally unverifiable mailing-list subscription;
- Support Site check now explicitly reports that `imputefinder` is registered
  in Watched Tags;
- CI installs current devel `BiocCheck`, rebuilds a clean detached clone, and
  runs both checks through `dev/bioc-gates.R`; either error or warning is fatal;
- exact CI shell/R gate executed locally against the release candidate → pass;
  workflow YAML + gate-script parse checks and `git diff --check` → pass.

M9/M10 and both remaining Definition of Done items are complete. ImputeFinder
0.99.2 is the first release candidate.

Exact next external step → push the release-candidate commit, observe the
hosted devel workflow, then begin Bioconductor new-package submission.

Blockers → none.

### 2026-07-16 - PLAN status reconciliation

Scope → reconcile the milestone checklist with the completed M9/M10 release
gate; identify whether any in-scope implementation work remains.

Evidence:

- `PLAN.md`'s sole open box was M9's aggregate errors/warnings/actionable-notes
  item;
- the immediately preceding clean-clone gate recorded `R CMD check --as-cran`
  status OK, `BiocCheckGitClone()` at 0 errors/warnings/notes, and new-package
  `BiocCheck()` at 0 errors/warnings with only reviewed optional/admin notes;
- M9, M10, the Definition of Done, and the full required test matrix were
  already closed from that same evidence;
- no other open checklist item exists in `PLAN.md` or this ledger.

The stale M9 aggregate box is now closed. Core PLAN work is complete; Section
16 remains intentionally deferred until the maintainer promotes a specific
item.

Exact next external step → push `main`, observe the hosted devel workflow, then
begin Bioconductor new-package submission.

Blockers → none.

### 2026-07-17 - successor innovation plan

Scope → promote the completed v1 frontier into a research-backed, falsifiable
successor program while preserving the release-candidate classifier as an exact
stable rail.

Evidence:

- clean baseline at `30618e1`; M0-M10 and the v1 Definition of Done remain
  complete;
- repository audit found the principal limits: one biological grouping axis,
  hard point states, synthetic-first validation, no design-confounding audit,
  no downstream workflow comparison, and no uncertainty propagation;
- current primary research supports continuous detection models, condition-
  associated non-detection endpoints, downstream-centric workflow evaluation,
  multiple-imputation/post-imputation uncertainty, peptide hierarchy, batch-
  associated missingness, and explicit workflow-step interaction studies;
- `PLAN.md` now contains a compact self-contained v1 contract, one input-first
  sidecar architecture, active A-C method/gate cards, an initial LFQ-only
  validation charter, M11-M15P dependency path, parked
  D-G/H1-H2/I/J1-J2 option cards,
  risk/kill controls, research anchors, and active-slice Definition of Done;
- the historical completion specification remains recoverable at
  `30618e1:PLAN.md`; current sessions route directly to the active method;
- `.agent/roadmap.md` exposes only M11-M15P as active. Every option card remains
  parked and cannot block the first successor release.

Verification → primary-source DOI/title checks; link probes resolved to source
targets or expected publisher access controls (one legacy journal host verified
separately after its DOI redirect timed out); `git diff --check` pass; no runtime
code changed, so package tests are unchanged and intentionally not rerun.

Exact next implementation task → M11: add stable method/schema/profile/cutoff/
rescue/policy identifiers as internal/companion metadata plus differential
fixtures, then introduce the companion analysis-object seam with identical
`imputefinder_result` behavior.

Blockers → none.

### 2026-07-17 - M11a stable identifiers + sidecar spec

Scope → freeze machine-readable v1 scientific/component identities and the
minimum experimental companion `spec` while proving the stable result remains
structurally and scientifically unchanged.

Implementation:

- internal method/profile/cutoff/rescue/policy identifiers name the exact v1
  rail without adding a field, attribute, or class to `imputefinder_result`;
- the minimum deterministic sidecar `spec` records analysis-schema, lifecycle,
  sentinel/stability module, scientific-component, package, and live namespace-
  version identifiers;
- spec construction is dependency-free, has no global state, and is absent from
  the stable classifier execution path;
- an exact serialized differential fixture locks the ten-field/one-class v1
  boundary around sidecar-spec construction.

Verification + completion evidence:

- focused red baseline → 3 expected errors from absent identifier/spec helpers;
- focused green suite → 3 tests / 13 expectations;
- package-wide source suite → 91 tests / 555 expectations, 0 failures,
  warnings, errors, or skips;
- frozen M7 scientific protocol verification → pass; 13 long scenarios, routine
  subset, and protocol MD5 unchanged;
- 10,000 x 50 manual/automatic performance gate → all runtime, allocation,
  result-size, RSS, cutoff-error, unchanged-input, and seed-cell checks pass;
- vignette-bearing source build + source-tarball `R CMD check --as-cran` →
  status OK;
- new-package `BiocCheck` → 0 errors, 0 warnings, 2 previously reviewed
  optional/admin notes;
- `git diff --check` → pass.

Exact next task → M11b: run the public/internal API naming experiment, add red
alignment/role/lifecycle tests, then implement the minimum typed
`missingness_design` without changing `classify_missingness()`.

Blockers → none.

### 2026-07-18 - M11b typed design + API boundary

Scope → run the sidecar naming experiment; freeze and export only the minimum
typed sample-design constructor; prove alignment/lifecycle failures are portable
and `rules_v1` remains byte-exact.

API decision:

- current Bioconductor naming/namespace guidance + base design-matrix vocabulary
  informed the rubric; exact official CRAN/Bioconductor documentation searches
  and 145 installed namespaces found no candidate collision;
- selected golden-path names = `missingness_design()`, future
  `analyze_missingness()`, and future `test_detectability()`; only implemented
  functions are exported; schema constructors, validators, alignment,
  fingerprints, compatibility comparators, and module runners stay internal;
- protocol, alternatives, caveat, and public/internal boundary are retained in
  `dev/api-naming-experiment.md`.

Implementation:

- experimental `missingness_design_v1` object stores schema/lifecycle, declared
  role-only sample metadata, canonical role map, and canonical explicit
  interaction sets; unselected metadata is omitted;
- required condition + optional nuisance/block/acquisition roles require exact
  columns and complete atomic values; condition/block/acquisition identifiers
  canonicalize to character while nuisance types remain declared;
- explicit unique sample row names + exact name-set alignment preserve caller
  data/order and deterministically align a copied design to analysis order;
- typed schema/role/alignment/lifecycle errors carry machine-readable context;
  unknown schemas/lifecycles fail with reconstruction guidance rather than
  silently upgrading;
- compact print method, generated help, NEWS, export, and development-version
  increment complete the first experimental public seam;
- constructor/metadata/printing leave the stable classifier signature,
  serialized normative result, fields, class, and attributes exact.

Verification + completion evidence:

- focused red baseline → 9 expected errors from the absent constructor;
- focused green design suite → 9 tests / 66 expectations;
- package-wide source suite → 100 tests / 621 expectations, 0 failures,
  warnings, errors, or skips;
- frozen M5 protocol verification → 12 scenarios / 13 target profiles pass;
- frozen M7 protocol verification → 13 long scenarios + routine subset pass;
  protocol MD5 remains `4011e381bba2d0d747e91d277a45de5e`;
- 10,000 x 50 manual/automatic performance gate → median `.415/.638` s,
  allocation `45.37x/66.53x`, result `4.70x/4.72x`, peak RSS `172.72` MiB,
  maximum automatic cutoff error `.282`; every frozen gate passes;
- vignette-bearing source build + source-tarball `R CMD check --as-cran` →
  status OK;
- committed-state clean-clone `BiocCheckGitClone` → 0 errors, 0 warnings, 0
  notes; new-package tarball `BiocCheck` → 0 errors, 0 warnings, 2 previously
  reviewed optional/admin notes; initial function-length/line-width notes
  exposed by the check were removed before the final run;
- `git diff --check` → pass.

Exact next task → M11c: compare current primary fingerprint primitives and
mask encodings under determinism, collision/security, cross-session/platform,
size, and dependency constraints; record the choice before implementation, then
add red sidecar round-trip/mismatched-input/unknown-schema tests.

Blockers → none.

### 2026-07-18 - M11c analysis schema + input identity

Scope → choose the deterministic matrix fingerprint and inline original-mask
representation before use; freeze the minimum companion schema around those
identities without constructing or exporting the input-first analyzer.

Research decision:

- current R/NIST primary documentation + an executable candidate harness favor
  dependency-free `tools::sha256sum(bytes = ...)`; SHA-256 supplies 128-bit
  collision and 256-bit preimage strength, while the faster OpenSSL path does
  not justify another runtime/system dependency for one construction/recheck;
- arbitrary R serialization was rejected as the canonical payload because R
  warns that its format may change; `matrix_core_be_v1` instead domain-tags and
  encodes storage, 32-bit big-endian dimensions, length-prefixed UTF-8-or-byte
  names, packed missingness, and finite observed values in column-major order;
- `base_packbits_lsb0_v1` stores the original mask at one bit/cell with zero
  padding: 62,531 serialized bytes at 500,000 cells versus 2,000,031 for a
  logical mask; unlike sparse indices, its size is density-invariant;
- frozen harness SHA-256 =
  `0857e92b191457779a71f67545964533d106d93b27e9b455a85f64f0e885dd60`;
  full protocol, alternatives, measured results, limits, and primary links →
  `dev/fingerprint-mask-experiment.md`.

Implementation:

- internal input records now retain dimensions, canonical ordered names,
  storage/representation/assay declarations, SHA-256 fingerprint, packed mask,
  trusted log2-scale declaration, and optional aligned acquisition values - no
  original numeric matrix;
- exact decoder validates mask encoding, cell/byte counts, and zero padding;
  canonical hashing additionally proves that supplied mask bytes describe the
  matrix before hashing;
- input verification recomputes identity and emits a typed mismatch with only
  dimension/name/mask/storage/fingerprint families - never intensity values;
- experimental `imputefinder_analysis` freezes seven fields: exact classic
  success or portable stage/class/message/call failure; sidecar spec; aligned
  typed design + empty estimability; input identity; empty sentinel/stability;
  and call/seed/hash/warning/failure/assumption/training provenance;
- schema, lifecycle, fingerprint-canonicalization, and mask-encoding changes
  require explicit reconstruction; malformed nested records fail typed checks;
- all constructors/validators remain internal, and `classify_missingness()` is
  absent from their execution path.

Verification + completion evidence:

- focused red baseline → 12 expected errors from absent identity/schema
  helpers; focused green → 12 tests / 80 expectations;
- exact canonical fixture SHA-256 agrees in source and a separately installed
  package process; text identifiers agree across declared UTF-8/Latin-1 while
  byte-marked identifiers remain exact;
- 160 randomized integer/double matrices spanning non-byte-aligned dimensions,
  missingness, signed zero, and negative values round-trip masks exactly,
  agree with independent `digest` SHA-256, detect observed-cell changes, and
  leave caller RNG state exact;
- package-wide source suite → 112 tests / 701 expectations, 0 failures,
  warnings, errors, or skips;
- frozen M5 protocol → 12 scenarios / 13 profiles pass; frozen M7 protocol →
  13 long scenarios + routine subset pass, MD5 unchanged;
- 10,000 x 50 v1 performance gate → manual/automatic medians `.195/.238` s,
  allocation `45.37x/66.53x`, result `4.70x/4.72x`, peak RSS `177.77` MiB,
  maximum automatic cutoff error `.282`; every frozen gate passes;
- vignette-bearing source build + source-tarball `R CMD check --as-cran` →
  status OK; committed-state clean-clone `BiocCheckGitClone` → 0 errors,
  0 warnings, 0 notes; new-package `BiocCheck` → 0 errors, 0 warnings,
  2 previously reviewed optional/admin notes;
- exact v1 serialization, fields/class/attributes, and public signature remain
  unchanged; `git diff --check` passes.

Exact next task → M11d: freeze pre-rescue evidence fields and red input-first
matrix/SE/failure/alignment/unchanged-input behavior, then implement original-x
construction behind the internal API. Defer compatible `base_fit` reuse to the
following bounded slice and export only when the selected public signature is
complete.

Blockers → none.

### 2026-07-18 - M11d input-first analysis + pre-rescue evidence

Scope → freeze the smallest immutable evidence record available before rescue;
construct the sidecar from original matrix/`SummarizedExperiment` input and a
typed design; preserve that record when the stable classic attempt fails.
Compatible `base_fit` reuse and the public analyzer remain deferred.

Schema + implementation:

- `pre_rescue_evidence_v1` lives under the existing sentinel seam and stores
  exact feature-condition/sample counts, missing fractions, and observed means;
  zero-observation blocks retain explicit counts + `NA` means, and no original
  numeric matrix is serialized;
- matrix designs align by exact sample-name set before their condition drives
  the shared core; `SummarizedExperiment` designs additionally require their
  named condition role to exist in and equal selected-assay `colData` sample by
  sample; typed failures identify absent/mismatched conditions without values;
- original mask/fingerprint/acquisition declarations are built before the
  stable attempt; classic success remains a complete `imputefinder_result`,
  while rescue/cutoff failures become portable stage/class/message/call records
  and remain mirrored in provenance;
- one extracted prepared-input core now drives both stable and sidecar paths;
  stage observation is local, and stable calls/results/RNG behavior remain
  isolated;
- the internal constructor records the normalized classic rescue seed and stays
  unexported until `base_fit` + module routing complete the selected public API.

Verification + completion evidence:

- focused red baseline → 6 expected absent-constructor errors; focused green
  matrix/SE/failure/alignment/schema/side-effect suite → 7 tests / 71
  expectations;
- package-wide source suite → 119 tests / 772 expectations, 0 failures,
  warnings, errors, or skips;
- exact matrix + SE classic parity, serialized sidecar round-trip, automatic
  cutoff + no-finite-condition failure survival, original-mask equality,
  unchanged source objects, caller RNG restoration, and byte-exact v1 fixture
  regressions pass;
- frozen M5 protocol → 12 scenarios / 13 profiles pass; frozen M7 protocol →
  13 long scenarios + routine subset pass, MD5 unchanged;
- frozen 10,000 x 50 v1 performance gate → every runtime, allocation,
  result-size, RSS, cutoff-error, unchanged-input, and seed-cell check passes;
  two input-first manual sidecar probes elapsed `.424/.957` s, with a `2.717`
  MiB pre-rescue record versus the `4.505` MiB original matrix;
- vignette-bearing build + source-tarball `R CMD check --as-cran` → status OK
  apart from the expected new-submission note; new-package `BiocCheck` → 0
  errors, 0 warnings, 2 previously reviewed optional/admin notes;
- no future public symbol/export or stable-v1 result/schema change;
  `git diff --check` passes.

Exact next task → M11e: define compatibility categories + red tests for exact
compatible reuse, every identity/scientific mismatch family, conflicting
`base_fit`/`cutoffs`, source-object immutability, and fit-only limitations; then
implement recomputation-backed optional `base_fit` acceptance + portable
`unavailable` without exporting `analyze_missingness()` yet.

Blockers → none.

### 2026-07-18 - M11e compatible base fit + portable unavailable

Scope → accept an optional stable fit only when the original input-first path
can reproduce it under its recorded cutoff-source policy + requested seed;
freeze fit-only abstention without exporting the incomplete analyzer.

Contract + implementation:

- `base_fit_compatibility_v1` partitions comparisons into exact returned data,
  representation, discrete decisions, rescue/group/audit/call components;
  canonical cutoff/profile structure + roles/labels/warnings; and continuous
  classification/cutoff/diagnostic/profile doubles under absolute + relative
  `sqrt(.Machine$double.eps)` bounds;
- ordinary numeric result matrices and ordinary numeric SE assays add canonical
  SHA-256 identities to exact comparison, preserving signed-zero distinctions
  that R's `identical()` omits; special numeric values remain position-exact;
- manual cutoff sources become explicit recomputation inputs; automatic sources
  are re-estimated; complete-only sources remain `not_needed`; any mismatch or
  recomputation failure returns a redacted typed error + complete category
  report, while a compatible supplied result is retained byte-for-byte;
- simultaneous `base_fit` + `cutoffs` fails as a typed specification conflict;
  malformed fit/policy inputs have dedicated portable error classes;
- `imputefinder_unavailable` freezes `status`/`quantity`/`code`/`message`/
  `requires`; fit-only dropped-row, original-mask, and pre-rescue requests use
  `original_input_required` and point to original `x` + typed design;
- helper extraction keeps the internal analyzer below BiocCheck's 50-line
  recommendation; the future public symbol/export remains absent.

Verification + completion evidence:

- focused red baseline → 8 expected missing-API/unused-argument errors; focused
  green suite → 10 tests / 93 expectations covering manual/automatic/
  complete-only success, identity/design/representation/seed/signed-zero drift,
  schema/policy/conflict failures, exact/canonical/tolerance partitions, RNG +
  source-object immutability, round-trip, and fit-only abstention;
- package-wide source suite → 129 tests / 865 expectations, 0 failures, errors,
  warnings, or skips;
- frozen M5 protocol → 12 scenarios / 13 target profiles pass; frozen M7
  protocol → 13 long scenarios + routine subset pass, MD5 unchanged at
  `4011e381bba2d0d747e91d277a45de5e`;
- frozen 10,000 x 50 v1 performance gate → manual/automatic medians
  `.500/.615` s, allocation `45.37x/66.53x`, result `4.70x/4.72x`, peak RSS
  `170.00` MiB, maximum automatic cutoff error `.282`; every gate passes;
- vignette-bearing source build + final source-tarball `R CMD check --as-cran`
  → status OK; an earlier check process was externally interrupted without a
  terminal status and was discarded before clean terminal reruns;
- committed-state clean-clone `BiocCheckGitClone` → 0 errors, 0 warnings, 0
  notes; a working-tree probe correctly rejected ignored local check artifacts;
  final new-package `BiocCheck` → 0 errors, 0 warnings, 2 previously reviewed
  optional/admin notes; its intermediate function-length note was removed;
- exact stable-v1 isolation/public signature regressions + `git diff --check`
  pass; no stable result/schema or public namespace change.

Exact next task → M11f: freeze compact differential fixtures + overhead gates
for the four comparison categories; add the complete `modules` argument routing
and public-wrapper red contract, then export/document `analyze_missingness()`
only when selected and unavailable module behavior is explicit.

Blockers → none.

### 2026-07-18 - M11f differential gate + public sidecar

Scope → freeze the four-category compatibility/performance oracle; complete the
golden-path analyzer signature and honest module router; export the experimental
sidecar only after selected/skipped/pending behavior is executable and
documented.

Contract + implementation:

- exported `analyze_missingness(x, design, assay = NULL, base_fit = NULL,
  cutoffs = NULL, seed = 1L, modules = c("sentinel", "stability"))` is a thin
  matched-call wrapper around the input-first core; stable classifier signature,
  result schema, attributes, and execution path remain unchanged;
- `modules` requires a unique exact subset of `sentinel`/`stability`,
  canonicalizes routing order, and accepts `character()` as an explicit skip;
  unselected slots are `NULL`, selected sentinel stores schema-validated
  pre-rescue evidence, and selected stability returns portable
  `module_pending_validation` requiring `robustness_certificate_v1`;
- stability validation accepts only the exact pending record until an actual
  validated module schema ships; corrupted status and invented output fail at
  the typed analysis-schema boundary;
- compact `print.imputefinder_analysis()` reports schema, dimensions, classic
  branch, and module availability without dumping retained evidence;
- public help, README/vignette path, package metadata/NEWS, version 0.99.4, and
  the accurate `ExperimentalDesign` Bioconductor view complete the exported
  seam;
- compact routine fixtures freeze exact discrete, canonical profile-warning,
  within/beyond-tolerance, and byte-exact stable-isolation cases;
  `dev/sidecar-differential-validation.*` freezes the representative overhead
  protocol before first measurement.

Differential/performance evidence:

- harness SHA-256 =
  `63c4a943efccfe81e1c6445f363c6db0941ae99b087a6ab7ec91f7684332b471`;
  unchanged M10 harness SHA-256 remains
  `74c53ccaf0f0c400693ac76607c575f686012c5b3980ae4d758a5dc4893fc2c2`;
- unchanged 10,000 x 50 stable gate → manual/automatic medians `.146/.168` s,
  allocation `45.37x/66.53x`, result `4.70x/4.72x`, peak RSS `176.57` MiB,
  maximum automatic cutoff error `.282`; every gate passes;
- compatible-base-fit default sidecar median `.426` s (`2.92x` stable),
  allocation `465.55` MiB = `103.34x` input / `2.28x` stable, largest allocation
  `.90x` input, result `5.48x` input / `1.16x` embedded classic; every frozen
  exact/input/mask/comparator/runtime/allocation/object-size gate passes.

Verification + completion evidence:

- focused red baseline → 6 expected absent-public-symbol errors; focused schema,
  analyzer, compatibility, input-first, spec, and differential suites green;
- package-wide source suite → 137 tests / 918 expectations, 0 failures, errors,
  warnings, or skips;
- frozen M5 protocol → 12 scenarios / 13 target profiles pass; frozen M7
  protocol → 13 long scenarios + routine subset pass, MD5 unchanged at
  `4011e381bba2d0d747e91d277a45de5e`;
- vignette-bearing source build succeeds; source-tarball `R CMD check --as-cran`
  completes with only the expected new-submission note; an initial unpublished
  report-link note was removed before the clean terminal rerun;
- new-package `BiocCheck` → 0 errors, 0 warnings, 2 previously reviewed
  optional/admin notes; its new actionable `ExperimentalDesign` suggestion was
  adopted before the terminal rerun;
- exact classic reuse, stable-v1 byte isolation, round-trip/schema corruption,
  caller-object immutability, generated documentation, harness hashes, and
  `git diff --check` pass.

M11 is complete. Exact next task → M12a: freeze generator-extension and
gate/data-card schemas, then inventory the minimum controlled + public LFQ
DDA/DIA candidates and their legal/scientific split feasibility before any data
role or numeric candidate result is opened.

Blockers → none.

### 2026-07-18 - M12a validation schemas + metadata-only data inventory

Scope → freeze the Section 9.1 generator-extension envelope, every A-C claim
family, executable gate/data/artifact schemas, and the smallest plausible
controlled/public DDA/DIA inventory without assigning roles or opening result
artifacts.

Frozen contract:

- dependency-free `dev/m12-validation-contract.R` owns exact ordered schemas,
  13 specified-unimplemented Tier 1-2 scenario identities, 25 A-C claims, an
  intentionally empty numeric gate registry, four data cards, six checksummed
  artifacts, four documented exclusions, and primary-source provenance;
- generator scope covers 2-4 conditions; balanced/unequal + independent/paired
  designs; run/batch effects with crossed/partial/perfect confounding; every
  Section 9.1 missingness pattern; outlier/heteroscedastic/support/low-evidence
  stress; DDA/DIA strata; null, flat-profile, causal-wording, and grouped-sibling
  leakage controls;
- gate rows are accepted only as complete Section 9.4 tuples with a finite
  numeric threshold, positive sample size, explicit uncertainty/failure rules,
  known data IDs, and SHA-256 data/protocol hashes - provisional rows cannot
  masquerade as frozen gates;
- role/opening state machine requires card/artifact role agreement + split
  SHA-256 before assignment and a local SHA-256 before artifact opening;
  M12a additionally fixes all roles to `unassigned`, all artifacts to
  `metadata_only`, local hashes absent, and gate-row count zero;
- canonical aggregate contract SHA-256 =
  `3f9a53d73495b0b7e3ce9ccf6a36ddcbabe4928f39cd5297e9a7cf79848b9e77`;
  component hashes + rationale live in `dev/m12-validation-contract.md`.

Candidate feasibility:

- MultiPro HCC1806/HS578T PXD041391 → CC0, balanced deliberate-batch DDA/DIA,
  3 biological replicates/class with technical/instrument siblings; 36 runs per
  acquisition; biological replicate = minimum split unit; natural missing cells
  retain no causal labels;
- MS-DAP PXD036134 → CC0 matched DDA/DIA, 3 known yeast levels x technical
  triplicates; whole matched family = one role; no biological-unit claim;
- PXD002099 → CC0 DDA UPS1/yeast series + compact non-normalized protein CSV;
  whole artifact = one role and an exact <=4-level subset/contrast is required;
- LFQbench PXD002952 → checksummed 2-condition DIA analysis under reviewed
  EMBL-EBI terms; all master-mixture/software/instrument derivatives stay one
  family; owner-rights caveat retained;
- exact PRIDE SHA-1 + byte identities recorded for every candidate artifact;
  SHA-1 is upstream transfer identity only and a post-role local SHA-256 is
  mandatory before opening;
- PXD041421, PXD028735, the rheumatoid-arthritis SWATH study, and PXD003881 are
  excluded with explicit artifact/independence/licence/scope reconsideration
  conditions rather than silently dropped.

Research boundary:

- only PRIDE/API/project metadata, repository terms, and primary papers were
  inspected; no result archive/CSV bytes or confirmation artifact were
  downloaded or parsed, no archive contents were listed, and no role was
  assigned;
- technical replicates never become independent split units; linked DDA/DIA,
  instruments, software outputs, and master mixtures cannot manufacture an
  independent confirmation;
- the manifest freezes scenario coverage, not acquisition physics, numeric
  generator parameters, estimator behavior, or evidence that A-C passes.

Verification + completion evidence:

- frozen contract verifier → schemas/catalog/hash exact; 4 adversarial rejection
  rails pass (`generator_scope`, `role_requires_split`,
  `opening_requires_role_hash`, `checksum_shape`);
- package-wide source suite → 137 tests / 918 expectations, 0 failures,
  errors, warnings, or skips;
- frozen M5 protocol → 12 scenarios / 13 target profiles pass; frozen M7
  protocol → 13 long scenarios + routine subset pass, MD5 unchanged at
  `4011e381bba2d0d747e91d277a45de5e`;
- vignette-bearing source build succeeds; source-tarball
  `R CMD check --as-cran --no-manual` completes with only the expected
  new-submission note;
- `git diff --check` passes.

Exact next task → M12b: freeze numeric generator parameters/seeds, implement
only the manifest scenarios/negative controls, and assign whole/grouped data
families to development vs sealed confirmation with exact split/protocol hashes
after rechecking terms. Validate generator truth + unchanged v1 behavior only;
external results remain unopened until the A-C numeric registry/protocol freeze.

Blockers → none.

### 2026-07-18 - M12b executable generators + sealed evidence roles

Scope → freeze numeric Section 9.1 generators/seed streams, validate generator
truth + unchanged stable-v1 behavior, and assign indivisible public families to
development or sealed confirmation before any external result bytes open.

Generator protocol:

- dependency-free `dev/m12-generator-validation.R` implements exactly the 13
  manifest scenarios with separate synthetic DDA/DIA parameter strata; latent
  complete intensities, component/union masks, generating probabilities,
  condition/design truth, and descriptive rule oracle remain explicit;
- model axes → 2-4 balanced/unequal conditions, independent/paired/technical-
  sibling units, correlated feature modules, run/batch/subject effects,
  crossed/partial/perfect confounding, intensity-independent/monotone/mixed/
  block/batch/structural-off-compatible/no-cliff masks, heteroscedasticity,
  outliers, differing support, and all predeclared negative controls;
- 64 explicit SHA-256-derived development seeds/scenario = 832 unique streams;
  local RNG kind = `Mersenne-Twister/Inversion/Rejection`, rescue seeds cycle
  `1/7/29`; endpoints 1/64 receive exact full-output fingerprints;
- generator protocol SHA-256 =
  `cdea1bf874152e63fba08c49e390da1dd18690105e43a2dc7a03fc6866d0d080`;
  endpoint audit =
  `4d216966c8cc54b1b23685d8801c7c65b2afcc5879dd15430d57a88f765e8a00`;
  all-832 compact audit =
  `a9f05c4b2beea35fbd2ac7729f10c95bb93a38c4494f24a4e8b7bb221ad70671`;
- all 832 seeds satisfy design/mask/probability/truth/RNG invariants; realized
  missingness stays within each frozen family (`.0300-.3619` overall); no
  generator result chooses or evaluates an A-C method.

Data-role freeze:

- live PRIDE v3 recheck confirms all project licences and every inventoried
  byte/SHA-1 identity; EMBL-EBI terms remain revised 2024-02-05 with
  attribution + owner-rights caveats; metadata-only inspection added the
  independently sourced HarmonizR mouse four-batch family PXD027467;
- A association: HarmonizR PXD027467 whole-family `development`; MultiPro HCC
  PXD041391 whole DDA/DIA family `confirmation`;
- B/C: MS-DAP PXD036134 matched DDA/DIA family `development`; UPS1/yeast
  PXD002099 DDA and LFQbench PXD002952 DIA `confirmation`;
- UPS1 exact initial-scope subset = 2/4/10/25 fmol triplicates; 50 fmol excluded
  before role assignment; each card's canonical role row supplies its exact
  split SHA-256;
- all 11 artifacts remain `metadata_only`, local SHA-256 absent; no archive,
  CSV, workbook, or confirmation result bytes were downloaded, listed, or
  parsed; shared-family opening remains synchronized per linked claims;
- v2 aggregate validation-contract SHA-256 =
  `bb5624fecd9cd52980bd901c1b75eac0d58dcb1869d95f92014e0d9e086a53ed`;
  six adversarial role/split/open/checksum/subset/scope rejection rails pass.

Stable compatibility + verification:

- all 13 first-seed simulations pass exact repeat, input/RNG preservation,
  reversed named-order, observed-value/rescue-log, and matrix/SE semantic gates;
  stable-v1 audit SHA-256 =
  `8463f4306fe9b6192336487e4ca365c8c7e3525dc06aeca5b66b1a99b1a10250`;
- package-wide source suite → 137 tests / 918 expectations, 0 failures, errors,
  warnings, or skips;
- frozen M5 protocol → 12 scenarios / 13 target profiles pass; frozen M7
  protocol → 13 long scenarios + routine subset pass, MD5 unchanged at
  `4011e381bba2d0d747e91d277a45de5e`;
- vignette-bearing source build succeeds; source-tarball
  `R CMD check --as-cran --no-manual` → only expected new-submission note;
  `git diff --check` passes.

Exact next task → M12c: populate complete numeric Section 9.4 rows for every
retained A-C claim, including explicit failure treatment, then freeze A/C
candidate comparisons + B perturbation protocol. Development bytes remain
closed until their linked hashes pass; confirmation remains sealed through
M15P.

Blockers → none.

### 2026-07-18 - M12c numeric gates + frozen candidate protocols

Scope → populate every retained A-C claim with executable numeric Section 9.4
gates/failure treatment, then freeze A/C candidate comparisons and B
perturbation/calibration mechanics before any candidate or external-result run.

Frozen evaluation contract:

- simulation replicates `1-32` = candidate selection/calibration; `33-64` =
  untouched development; all protocol states remain `frozen_unrun`;
- 27 claims map to 66 numeric rows (`A=20`, `B=18`, `C=28`), each binding the
  exact metric/unit/comparator/strata/n/uncertainty/operator/threshold/failure,
  canonical evidence hash, and track-protocol hash;
- result evaluator requires every requested gate exactly once + exact frozen-row
  hash; measured endpoints are finite, while failed/unavailable rows carry
  `NA` endpoints and fail - fabricated unavailable measurements are rejected;
- 20 contract/result rails pass, including nonempty staged subsets +
  adversarial drift;
  aggregate validation-contract SHA-256 =
  `2cd5ecd8bf5763da1c5d9e1d8994207e701bcea9a7efcf4ae539dfd9b52d7431`.

Candidate/perturbation freeze:

- A canonicalizes labels/terms, uses SVD rank + an ordered null-projector basis
  for rotation-invariant aliases + contrast row-space estimability, preserves
  biological/block units, and compares
  OLS-HC3/CR2, studentized restricted Freedman-Lane, and independent-unit
  quasibinomial association; support/permutation-resolution floors, exact vs
  `9999`-draw tails, Holm families, ranking, and no-winner outcome are explicit;
- B separates 999 biological/block bootstraps + exhaustive leave-one-unit-out,
  999 module/feature bootstraps with cutoff refit, and six deterministic policy
  scenarios; multiplicity weights preserve names; cutoff success floor =
  `900/999`, range = type-8 `.025/.975`; display labels use draws `1-499` for
  scoring and `500-999` for validation, remain optional, and cannot kill
  continuous panels;
- C separates empirical-row standardized detection-risk difference from
  conditional observed-abundance contrast; component truth/nonzero alternatives
  are distinct (five detection vs seven abundance scenarios), preventing
  abundance-only generators from inflating detection power; detection candidates
  = GLM-HC3, mean bias reduction, GEE+Mancl-DeRouen, and odds-only exact
  conditional auxiliary; abundance = OLS-HC3/CR2 + ordinary/robust `limma`;
  BY families, support/separation rules, and 999 whole-block bootstrap mechanics
  are frozen separately;
- track-specific SHA-256 seed keys exactly match executable stream generation;
  eight candidate-protocol mutation/input rails pass; A/B/C/bundle hashes =
  `e20bdcc8e040eed65457480be7ae2ce2542761165c3d27161b953f86f39e5edc` /
  `57ffa1220c740dcb69508d3ad1b9beb05e7013c3a92af528293211e41a13f256` /
  `58b6e3c6d7f35f5561549a3a52420e49b0f3eb296752353d1bb3c9d433cf9d96` /
  `8eb5592131407ae48d4b4caa06ef13eea84371471215eae41997faf0e5ec0821`.

Evidence boundary + review:

- public B evidence was corrected to finite-dataset external behavior -
  controlled abundance ratios do not identify a true automatic cutoff; C alone
  retains ratio/direction recovery truth;
- all 11 external artifacts remain `metadata_only`, local SHA-256 absent; no
  result archive/CSV/workbook bytes were downloaded, listed, or parsed;
- primary/current review covers Holm/BY, type-8 quantiles, mean bias reduction,
  exact conditional likelihood, GEE + Mancl-DeRouen, `limma`, restricted
  permutation, cluster bootstrap, CR2/Satterthwaite, and risk-coverage;
  package availability selects no winner.

Verification + completion evidence:

- M12 contract/candidate verifiers → 66 frozen rows; all 20 + 8 protocol rails
  pass; generator protocol/audit unchanged at
  `cdea1bf874152e63fba08c49e390da1dd18690105e43a2dc7a03fc6866d0d080` /
  `4d216966c8cc54b1b23685d8801c7c65b2afcc5879dd15430d57a88f765e8a00`;
- expanded stable-v1 audit passes unchanged at
  `8463f4306fe9b6192336487e4ca365c8c7e3525dc06aeca5b66b1a99b1a10250`;
- package-wide source suite → 137 tests / 918 expectations, 0 failures,
  errors, warnings, or skips;
- frozen M5 → 12 scenarios / 13 target profiles pass; frozen M7 → 13 long +
  routine scenarios pass, MD5 unchanged at
  `4011e381bba2d0d747e91d277a45de5e`;
- vignette-bearing source build + source-tarball
  `R CMD check --as-cran --no-manual` → status OK; `git diff --check` passes.

M12 is complete. Exact next task → M13a: write red mandatory-core tests, then
implement the canonical declared model matrix, SVD rank/null-space alias report,
contrast row-space estimability, and block/technical-unit accounting. Keep A
association candidates + development bytes closed until the mandatory core is
green; confirmation remains sealed through M15P.

Blockers → none.

### 2026-07-18 - M13a mandatory algebraic design core

Scope → implement and validate only the mandatory shared algebraic core needed
by A-C: canonical declared model, SVD rank/null-space aliases, row-space
contrast estimability, and declared block/technical-unit accounting. A
association candidates and every external result artifact remain closed.

Implementation:

- every new sidecar populates its existing `design$estimability` seam with
  schema `design_estimability_v1`; older experimental objects with the prior
  `NULL` placeholder remain structurally readable;
- canonical model rows sort sample identities; observed categorical labels sort
  independently of factor levels/character encoding; integer/double nuisance
  re-encoding agrees; intercept + treatment contrasts cover condition,
  declared nuisance main effects, block fixed effects, and only explicit
  interactions; acquisition is retained as declared metadata and enters only
  an explicit interaction;
- exact variable/term/coefficient maps accompany the numeric matrix; no formula
  environment or ambient contrast option enters the stored result;
- SVD rank threshold = `max(nrow, ncol) * .Machine$double.eps * d[1]`;
  null projector = `I - V_r V_r^T`; canonical null basis scans coefficient
  axes with modified Gram-Schmidt + fixed first-nonzero sign; affected
  coefficients/terms retain alias magnitudes;
- declared-level and raw coefficient contrasts use the frozen row-space
  residual rule `L2 <= sqrt(.Machine$double.eps) * max(1, contrast_L2)` and
  return named affected terms for non-estimable directions;
- declared block = independent unit; otherwise sample = unit; sample weights,
  unit/condition multiplicities, and within-condition technical siblings stay
  explicit and name-aligned;
- sidecar validation recomputes the core with continuous tolerance, rejecting
  schema drift or inconsistent model/rank/alias/unit evidence; inputs, stable
  classic branches, RNG, options, and working directory remain untouched;
- package version advances to `0.99.5`; README, vignette, NEWS, package help,
  and `analyze_missingness()` help document the always-on algebraic core and
  its non-association/non-correction boundary.

Frozen development evidence:

- dependency-free `dev/m13-design-core-validation.R` binds only the 11 M12
  mandatory gate rows and regenerates all 13 DDA/DIA scenarios at untouched
  replicates `33-64` = 416 instances plus the normative v1 parity case;
- all required sensitivity/specificity/retention/accounting/invariance
  estimates = `1`; exact-alias false-positive fraction = `0`; global-state
  mutation and classic exact-content drift counts = `0`;
- complete audit SHA-256 =
  `37528216590dcfff3338501ef61b766697ea63c14b7324912fcb00f84d58790e`;
  registry-bound gate-result SHA-256 =
  `cb72a69301100287cee569ae7a0d709f42edc10e3b10a62dbff3a86ab88aec6f`;
- verifier asserts the upstream generator hash plus candidate/external sealed
  state before running; A association remains `frozen_unrun`, all 11 external
  artifacts remain `metadata_only`, and no external byte was downloaded,
  listed, parsed, or assigned a new role;
- `dev/m13-design-core-validation.md` records exact denominators, comparison
  normalization, evidence identities, boundaries, and reproduction command.

Adversarial review + verification:

- red M13a suite initially failed on the absent core; final package-wide source
  suite → 144 tests / 996 expectations, 0 failures, errors, warnings, or skips;
- constructed crossed vs perfect condition-batch designs, `p > n` null spaces,
  one-sample/one-condition support, multi-level contrasts, paired technical
  siblings, unequal replication, schema mutation, and factor/sample/role-order
  variants pass;
- frozen M13 verifier reproduced both hashes after install and after the final
  helper refactor; all 11 numeric registry gates pass exactly;
- frozen M5 → 12 scenarios / 13 target profiles pass; frozen M7 → 13 long +
  routine scenarios pass, MD5 unchanged at
  `4011e381bba2d0d747e91d277a45de5e`;
- vignette-bearing source build succeeds; source-tarball
  `R CMD check --as-cran --no-manual` → status OK;
- staged-tree synthetic commit `BiocCheckGitClone` → 0 errors, 0 warnings, 0
  notes; new-package tarball `BiocCheck` → 0 errors, 0 warnings, only the two
  reviewed external/admin notes; initial function-length/line-width notes were
  removed by helper extraction and source wrapping;
- `git diff --check` passes.

Exact next task → M13b: write red coverage-schema/order/zero-cell tests, then
implement static pre-rescue sample/condition/nuisance/block summaries from the
typed design + stored original evidence. Keep association fitting and all
external result bytes closed until the static sentinel schema is green.

Blockers → none.

### 2026-07-18 - M13b static sentinel coverage

Scope → implement deterministic descriptive pre-rescue coverage from the typed
design + original evidence. Fit no association candidate; open no candidate,
development, external, or confirmation result evidence.

Implementation:

- selected sentinel output retains unchanged `pre_rescue_evidence_v1` and adds
  `sentinel_static_coverage_v1` with exact sample, condition, role-level,
  condition-role, and pairwise feature-overlap tables;
- sample detection fractions use globally observable proteins - finite in at
  least one original cell - matching the frozen A response denominator; full
  input counts and observed intensity mean/minimum/median/maximum remain beside
  them;
- condition summaries retain sample/full/eligible/detected/complete support,
  observed + eligible cells, and detection fraction;
- role-level evidence includes every declared condition/nuisance/block/
  acquisition value, numeric/categorical encoding, numeric value where
  applicable, sample/condition counts, and singleton flag;
- condition-role evidence materializes the full product of canonical conditions
  with every nuisance/block/acquisition value; zero-sample cells remain explicit
  with zero support, `NA` detection fraction, and `empty = TRUE`;
- feature overlap reports shared/union/side-only/neither counts + Jaccard for
  every unordered condition pair; global absences never inflate denominators;
- sample/features/conditions/roles/labels canonicalize independently of input,
  factor-level, nuisance-selector, and integer/double representation order;
- deserialization validation reconstructs structural tables exactly from the
  packed original mask + typed design, cross-checks pre-rescue means and
  intensity-range consistency, and rejects lifecycle/schema/count corruption;
  prior pre-rescue-only experimental objects remain structurally readable;
- one redundant constructor validation pass was removed; terminal representative
  sidecar gates remain inside every frozen runtime/allocation/object-size bound;
- package version advances to `0.99.6`; README, vignette, NEWS, package help,
  generated Rd, and `dev/m13-static-coverage-validation.md` document exact
  denominators, zero cells, descriptive scope, and integrity limits.

Adversarial review + verification:

- focused red baseline → 57 failures + one lifecycle error at the absent seam;
  final focused suite → 8 tests / 81 expectations; package-wide source suite →
  152 tests / 1,077 expectations, no failures/errors/warnings/skips;
- singleton levels, categorical + numeric nuisance grids, empty condition-role
  cells, zero-support samples, wholly absent input, feature-overlap arithmetic,
  named order/re-encoding invariance, corruption, legacy readability, and
  unchanged caller inputs/global state pass;
- installed frozen M13 verifier → 416 instances, all 11 gates exact; audit/
  result hashes unchanged at
  `37528216590dcfff3338501ef61b766697ea63c14b7324912fcb00f84d58790e` /
  `cb72a69301100287cee569ae7a0d709f42edc10e3b10a62dbff3a86ab88aec6f`;
  association remains `frozen_unrun`, all external artifacts `metadata_only`;
- frozen M5 → 12 scenarios / 13 target profiles; frozen M7 → 13 long + routine,
  MD5 `4011e381bba2d0d747e91d277a45de5e`; all pass;
- uncontended 10,000 × 50 differential gate → sidecar median `2.257` s =
  `4.57x` stable; allocation `622.85` MiB = `3.05x` stable; largest allocation
  `.90x` input; object `1.17x` embedded classic; every frozen gate passes;
- vignette-bearing build succeeds; source-tarball
  `R CMD check --as-cran --no-manual` → only expected new-submission note;
  staged-tree `BiocCheckGitClone` → 0 errors/warnings/notes; new-package
  `BiocCheck` → 0 errors/warnings + three reviewed notes; automatic `Coverage`
  view suggestion rejected as genomic/read-coverage metadata, not design/data
  support; one intermediate clone check was discarded when temporary-worktree
  cleanup raced its yielded process, then the explicit terminal rerun passed;
- `git diff --check` passes. No A candidate fit or external-result access
  occurred.

Exact next task → M13c: write red association output/support/abstention tests,
implement the frozen OLS-HC3/CR2, restricted Freedman-Lane, and independent-unit
quasibinomial candidates, then run only simulation candidate replicates `1-32`
for selection/no-winner. Keep development replicates `33-64`, HarmonizR bytes,
and confirmation sealed until selection locks.

Blockers → none.

### 2026-07-19 - M12 pre-result literature audit + v2 correction

Scope → use newly available authenticated/OA full text to adversarially audit
the frozen M12 sources and method claims before M13c opens any A candidate
result. Preserve the M13a/b implementation/evidence and every external seal.

Source audit:

- 24 pre-existing contract/PLAN DOI anchors were title-resolved: 23 matched;
  `cluster_bootstrap_2013` incorrectly used
  `10.1016/j.jmva.2012.10.006`, an unrelated almost-periodic time-series paper;
  corrected source = Cheng/Yu/Huang cluster-bootstrap consistency,
  `10.1016/j.jmva.2012.09.003`;
- full CR2 author manuscript + official 2023 corrigendum constrain any
  absorbed-fixed-effect shortcut to ordinary unweighted LS, identity working
  model, and rank conditions; leverage-dependent Satterthwaite df can be much
  lower than cluster count;
- Winkler exchangeability/variance-group conditions + Helwig robust-W results
  leave small-n Freedman-Lane calibration empirical; independent robust-W
  evidence does not directly validate the paired CR2 extension;
- MSqRob's author version models peptide detections within each protein and
  moderates protein-specific dispersions - not one all-protein count/sample;
  the global quasibinomial candidate is now explicitly analogy-only;
- full RA SWATH methods show disease-category/batch overlap + ad hoc count
  adjustments, strengthening its exclusion from confounder-adjusted validation;
- Taylor & Francis original CR2 and ACS final MSqRob copies remained gated;
  final DOI metadata + titles were checked, while claims use author versions,
  open full text, and the open official correction. No human access action is
  outstanding.

Pre-result freeze:

- contract/candidate/gate versions advance to v4/v2/v2; A/C CR2 uses the full
  design + identity working model and returns named unavailable below computed
  Satterthwaite df `5`;
- restricted Freedman-Lane recomputes scalar HC3/CR2 robust Wald statistics and
  permits only null-invariant transforms within compatible nuisance,
  exchangeability, and variance groups; incompatible layouts abstain;
- candidate protocol hashes A/B/C/bundle =
  `ca0cf8dbbf082446d9116ce280f0acf2c6517d35ec6137edb7b415585ce92683` /
  `20912c4dfadeebb0c1da7100cb54dfa27f202ea39c3e0d991fb0ee38bab383c2` /
  `f28b01a4cf5373807b4c71241a3563bc64cf55aa8d446cb62e242c9ecae76b04` /
  `3bc4ef531d6c84475114befc5c153e50680459f1d8125fdbabac5fd840a1de65`;
  source/aggregate contract hashes =
  `d81aae9a43a48b692bec5d248d4d78c711a2918be9632c2c06f3a4ff27fea691` /
  `45d1cda936b21a42d6735d900917734398eecb3ff1c136cab486cf90ee2b21e5`;
- M13a's 416-row audit remains byte-identical at
  `37528216590dcfff3338501ef61b766697ea63c14b7324912fcb00f84d58790e`;
  registry-only rebinding changes its 11-row result hash to
  `816303e18fb165d90891833e46e3d056e3d97289e4e2d267460d9e61b9f85144`
  with all estimates/thresholds unchanged.

Verification:

- M12 contract → 20/20 schema/gate rails; candidate protocol → 9/9 mutation
  rails; generator hashes/audit unchanged; all pass;
- installed M13 verifier → all 11 gates pass across 416 instances; association
  stays `frozen_unrun`, external artifacts stay `metadata_only`;
- package-wide source suite → 152 tests / 1,077 expectations, no failures,
  errors, warnings, or skips; vignette-bearing source build succeeds;
- `git diff --check` passes. No A/B/C candidate result, external result, or
  sealed artifact byte was opened.

Exact next task → M13c: write red association output/support/abstention tests,
implement the v2-frozen OLS-HC3/CR2, restricted Freedman-Lane, and
independent-unit quasibinomial candidates, then run only simulation candidate
replicates `1-32` for selection/no-winner. Keep A candidate computations on
development replicates `33-64`, HarmonizR bytes, and confirmation sealed until
selection locks.

Blockers → none.

### 2026-07-19 - M13c pre-result association addendum v3

Scope → close Track A ambiguities before any candidate computation while
preserving every M12 threshold, candidate, generator, evidence role, and seal.

Freeze:

- standalone executable `m13_a_association_protocol_v3` binds the immutable
  M12 contract/gate/A-protocol/bundle and generator protocol/audit hashes;
- globally observable denominator remains fixed over the complete input;
  acquisition levels define separate rebuilt model/Holm strata;
- hypotheses are canonical 1-df condition/nuisance coefficients and explicit
  condition/nuisance-only product columns; intercept, block, and acquisition
  remain adjustment/descriptive-only;
- support now has exact categorical sides, numeric unit-position median sides,
  Cartesian interaction cells, independent `4`/side and blocked `6`-complete-
  block floors; red-test derivation exposed multi-numeric interaction loss, so
  typed component encodings + exact per-component numeric support were added
  before candidate computation and adversarially rereviewed;
- production red tests exposed incidental `vapply()` names that rejected every
  declared-acquisition result; the comparison now strips only those names and
  a declared-stratum fixture locks the valid path before candidate computation;
- HC3, full-design identity-working-model CR2/Satterthwaite, constrained-null
  Freedman-Lane, empirical independent quasibinomial, abstention, covariance,
  transformation-count, and canonical seed/hash rules are executable;
- Monte Carlo transforms now have an exact seeded stream-to-map algorithm,
  identity/duplicate rejection, retry audit, map hash, and RNG restoration;
- exact nested association/result/unavailable schemas + candidate-agnostic
  opportunity denominator prevent silent successful-fit denominator changes;
- top-level abstention represents zero observable features/no hypotheses;
  winner panels bind v3 gate, implementation, and candidate-evidence hashes;
- fixed alternative targets = DDA crossed `batch[batch_2]`, DDA unequal
  `condition[D]`, DIA partial `batch[batch_3]`, DIA paired `condition[B]`;
- corrected effective Track A scope for `dda_monotone_unequal`; candidate gates
  bind replicates `1-32`; every screening-qualified candidate faces the full
  gate before deterministic ranking; development opens only after winner lock;
  confirmation remains M15P-only;
- adversarial architecture/contract/mathematical rereview closed all pre-seal
  findings; historical M12 objects remain byte-identical.

Verification:

- executable contract self-tests reject upstream/component/hash/allocation/
  stage drift, incomplete gates, malformed result bindings, and prove canonical
  collision-safe seed/map ordering plus RNG restoration;
- protocol hash =
  `254e9cb03d3f981bbbcf621b8e13985c8abb43916ab9b79901b32bac13b55c59`;
- aggregate v3 hash =
  `04e824321ac68c73633277064397f1a3fa96ba73ff342fbfe575df7fc6b1f4a6`;
- `Rscript --vanilla dev/m13-association-contract.R --verify` passes all 27
  named self-tests with exact component hashes;
- `git diff --check` passes; no candidate simulation ran and no external result
  artifact was read.

Exact next task → red tests for response, strata, coefficient eligibility,
support boundaries, available/unavailable records, Holm, canonical mask/design
hash, implementation/evidence manifests, candidate opportunity denominators,
HC3/CR2, restricted transformations, and quasibinomial effects. Keep all result
allocations sealed until those rails are green and the implementation manifest
locks.

Blockers → none.

### 2026-07-19 - M13d association preparation + identity rails

Scope → implement the candidate-agnostic response/hypothesis/support boundary
under the resealed v3 contract; retain all candidate/development/confirmation
result allocations sealed.

Implementation:

- `association_mask_design_be_v1` now has an independent production encoder;
  its gold digest matches the executable contract and is invariant to row/
  column order, factor levels, integer/double representation, role order,
  interaction order, and finite intensity values while remaining sensitive to
  mask, design, roles, and interactions;
- `association_preparation_v1` stores only packed original-mask/global typed-
  design identity, global sample detection fractions, acquisition-specific
  rebuilt cores, eligible coefficient axes, and common-support records - no
  intensity or candidate result enters it;
- acquisition strata remove acquisition-containing terms before rebuilding;
  condition/nuisance main and eligible explicit interaction axes retain exact
  raw-coefficient estimability, including nonestimable rows;
- support counts unique biological-unit positions, numeric medians/ties/ranges,
  coefficient-specific categorical targets, Cartesian interaction cells, and
  six-complete-block requirements; technical sibling multiplicity cannot
  inflate support;
- top-level zero-observable/no-hypothesis states use the frozen exact
  `imputefinder_unavailable` shape;
- validation reconstructs hash, response, stratum IDs, exact sample bijections,
  designs, cores, hypotheses, and support from the packed global identity;
  detached hashes/masks/designs, sample omission, forged/reordered strata,
  empty available state, missing hypotheses, and support drift reject.

Adversarial corrections:

- math/architecture review found all multilevel main coefficients initially
  targeted the first nonreference level; coefficient-level lookup + A/B/C
  regression now distinguish A/B from A/C support;
- contract/architecture review found the initial preparation validator trusted
  a shape-only SHA and subset-directional strata; packed identity + exact
  reconstruction closes both trust gaps;
- implementation tests exposed the frozen result validator's incidental
  `vapply()` names bug; contract protocol/aggregate identities were resealed to
  `254e9cb03d3f981bbbcf621b8e13985c8abb43916ab9b79901b32bac13b55c59` /
  `04e824321ac68c73633277064397f1a3fa96ba73ff342fbfe575df7fc6b1f4a6`
  in commit `59e0ccd` before any candidate computation.

Verification:

- focused preparation suite → 9 tests / 85 expectations; full source suite →
  161 tests / 1,162 expectations; zero failures/errors/warnings/skips;
- three-agent rereview → no remaining identity, lifecycle, hypothesis,
  support-math, or stratum-rebuild blocker;
- frozen executable contract → 27/27 named self-tests with exact hashes;
  vignette-bearing source build/check → `Status: OK`; `git diff --check` passes;
- no candidate simulation ran; development replicates `33-64`, HarmonizR
  bytes, and confirmation remain sealed.

Exact next task → red/green the shared rank-coordinate + HC3/CR2 engine and
final unavailable handoff, then restricted Freedman-Lane and independent
quasibinomial. Freeze exact source/environment/test implementation manifest
only after all engine/handoff rails pass; then open candidate replicates `1-32`.

Blockers → none.

### 2026-07-19 - M13e robust association candidate

Scope → implement and adversarially validate only the frozen
`a_fraction_ols_hc3_cr2` candidate; keep every candidate replicate,
development/confirmation allocation, external byte, final panel, and
implementation/evidence manifest sealed.

Implementation:

- `association_candidate_artifact_v1` retains the exact response,
  hypothesis/support joins, robust available/unavailable outcomes,
  stratum-local Holm families, typed-empty seed record, and diagnostics;
- trusted validation requires the genuine `association_preparation_v1`, checks
  provenance before dereference, restricts this slice to robust-only failure
  stages, enforces estimability/support/residual-df/degenerate precedence, and
  rejects detached hashes, impossible codes, Holm/family drift, typed-seed
  drift, and diagnostic joins;
- retained-SVD coordinates use `Z=U D`, `a=D^-1 U'y`, and cancellation-stable
  HC3/CR2 sandwiches without a second `solve()` rank cutoff;
- independent fits apply HC3 + the registered leverage denominator floor;
  blocked fits keep the complete fixed-effect design, identity-working-model
  symmetric MP CR2 block correction, and scalar Satterthwaite reference df;
- covariance cleanup symmetrizes, rejects materially negative eigenvalues, and
  removes only tolerance-scale eigencomponents, preserving every retained
  direction unchanged when no clamp is required; CR2 df projections normalize
  their homogeneous scale before Gram arithmetic;
- the exact registered `nu>=5` floor remains a literal computed-double rule:
  a mathematical boundary that rounds below `5` abstains. Any tolerance policy
  change requires a pre-result protocol reseal.

Adversarial corrections:

- the initial global `1-h` gate wrongly rejected valid blocked leverage-one
  rows; it now applies only to HC3 while CR2 carries null directions through
  its MP correction;
- default `solve(crossprod(Z))` introduced an undeclared squared-condition
  cutoff and discarded SVD-retained scaled axes; direct retained-singular-value
  algebra removes it;
- unconditional eigen reconstruction damaged legitimate small covariance
  directions beside a large one; low-rank-only correction restores the frozen
  clamp semantics;
- artifact validation originally allowed provenance-free calls, cross-candidate
  failure codes, one-way precedence, and logical typed-empty seed columns; all
  four trust gaps now reject under mutation tests;
- validation proves preparation provenance, schema, joins, and internal
  arithmetic; it intentionally does not independently recompute a coherently
  forged numerical fit. The later genuine-run implementation/evidence hash
  locks output drift before any candidate result is accepted.

Verification:

- focused robust suite → 7 tests / 90 expectations; independent hard HC3,
  blocked hard CR2, blocked leverage-one, rank/scaling, covariance cleanup,
  low-reference, acquisition-Holm, abstention precedence, and validator
  mutations all pass;
- full source suite → 168 tests / 1,252 expectations, zero failures/errors/
  warnings/skips; three-agent final rereview → no remaining math, contract, or
  architecture blocker;
- frozen executable contract → 27/27 named self-tests with exact protocol/
  contract identities; vignette-bearing source build/check → `Status: OK`;
  `git diff --check` passes;
- no candidate simulation ran; replicates `1-32`, development `33-64`,
  HarmonizR bytes, and confirmation remain sealed.

Exact next task → red/green the constrained-null transformation identity,
exact/Monte Carlo map/seed manifest, and restricted Freedman-Lane candidate on
local unit fixtures. Then implement the independent quasibinomial candidate +
final handoff before locking the implementation manifest and opening candidate
replicates `1-32`.

Blockers → none.

### 2026-07-19 - M13f restricted Freedman-Lane candidate

Scope → implement and independently validate only the frozen
`a_fraction_freedman_lane` candidate on local fixtures. Preserve all synthetic
candidate/development allocations, external result bytes, final panel, and the
implementation/evidence manifest as sealed.

Pre-result amendment:

- implementation exposed three impossible/ambiguous frozen requirements before
  any result allocation: delimiter-joined seeds excluded valid acquisition
  labels, the blocked helper contradicted radix block ordering, and finite
  doubles could not hold every exact transformation count;
- the authorized v3 reseal uses a domain-separated typed length-prefixed seed
  key, caller-preserved radix block order, `log(T)` branch decisions, and
  finite-double saturation only for Monte Carlo diagnostics;
- global seed collision resolution still covers every planned Monte Carlo
  instance, while result manifests retain only available outcomes; every
  byte-identical residual stabilizer refits the literal original response;
- the final seal also binds the complete-QR CR2 residual projector and rejects
  malformed seed types/keys, duplicate per-hypothesis map hashes, and retry
  counts beyond the declared limit in both executable-contract and production
  validation;
- protocol/contract/gate identities now equal
  `5f484f614d72d9560d96e2b8f75b71cc1027950958ea1459be93eece595069f5` /
  `85d121ebd5ddd589be0f389553a3b7ffa6076a12a8fcbcb476a8969f7d987ae4` /
  `f2fd8d0171003e0ae9b31757acf2d43acf7477965d517f1757f4d895c6abe731`.

Implementation:

- each coefficient restriction is constructed in the retained full-model rank
  coordinates; the constrained-null projection produces `y0_hat`, residuals,
  and the exact frozen `P Z0` tolerance;
- independent transforms use pairwise-complete, transitive equivalence groups;
  exact maps are globally lexicographic. Blocked transforms use radix block
  swaps of complete condition/nuisance-matched sibling bundles;
- `20 <= T <= 100000` enumerates every distinct map; larger spaces generate
  9,999 distinct seeded nonidentity maps with exact RNG restoration, retry +
  big-endian map-hash audit, and frozen `+1` Monte Carlo arithmetic;
- every transform recomputes the full scalar HC3/CR2 Wald statistic through a
  cached-response refit. A complete-QR orthogonal-complement projector prevents
  roundoff-only CR2 PSD failures at the first blocked Monte Carlo size;
- hypotheses plan/evaluate/discard maps sequentially. Artifact validation
  regenerates maps, refits all transformed responses, recomputes exceedance,
  p-value, disposition, Holm families, and the available-only seed manifest;
- low resolution, incompatible layouts, retry exhaustion, and transformed-fit
  failure enter the exact frozen per-hypothesis unavailable path.

Adversarial corrections:

- independent connected-component traversal initially reduced the wrong matrix
  margin; a multi-group exact-map fixture now locks the correct partition;
- identity-only response preservation undercounted nonidentity residual
  stabilizer ties; literal residual equality now selects original response
  bytes for every stabilizer map;
- a shape-valid artifact could originally forge coherent exceedance/p-value
  arithmetic; provenance validation now reruns every transformation;
- all-hypothesis map retention and exact planning duplicated avoidable memory;
  enumeration/generation is now delayed to one evaluated hypothesis at a time;
- retry exhaustion originally aborted the full candidate and unavailable Monte
  Carlo rows leaked seed records; both now obey per-hypothesis disposition and
  result-manifest rules.

Verification:

- focused restricted-permutation suite → 13 tests / 149 expectations;
  focused robust suite → 7 tests / 90 expectations; both pass without
  failure, warning, error, or skip;
- package-wide source suite → 181 tests / 1,401 expectations, zero failures,
  warnings, errors, or skips;
- executable frozen contract → 30/30 named self-tests; exact new protocol,
  component, contract, gate, allocation, and stage identities reproduce;
- vignette-bearing source build succeeds; source-tarball
  `R CMD check --as-cran --no-manual` passes every substantive rail with only
  the expected new-submission note; `git diff --check` passes;
- no candidate simulation ran; replicates `1-32`, development `33-64`,
  HarmonizR bytes, and confirmation remain sealed.

Exact next task → red/green the empirical independent-unit quasibinomial
candidate, then final panel/handoff + implementation manifest. Open candidate
replicates `1-32` only after all three candidates and manifest rails lock.

Blockers → none.

### 2026-07-19 - M13g empirical quasibinomial candidate

Scope → implement and adversarially validate only the frozen empirical
independent-unit `a_fraction_quasibinomial` candidate. Preserve synthetic
candidate/development allocations, external bytes, public panel, and the
implementation/evidence manifests as sealed.

Implementation:

- grouped `(k,G-k)` counts fit through the exact frozen `stats::glm.fit`
  controls in retained SVD coordinates; every count/probability boundary,
  convergence result, final weighted-rank check, and literal finite `phi>0`
  rule has an exact structured disposition;
- Pearson dispersion and the weighted information inverse pass the common
  covariance PSD/scalar tolerance rails; the primary 0→1 encoded-column delta
  effect uses its analytic gradient and unclamped residual-df t interval;
  log-odds effect, SE, interval, statistic, and raw p remain secondary;
- only independent units are licensed. Candidate records retain common
  estimability/support/residual-df precedence, paired-scope abstention,
  per-contrast continuation, stratum-local Holm families, typed-empty seed
  provenance, and exact preparation-bound fit/effect replay;
- common covariance/failure/multiplicity infrastructure now lives in the
  candidate boundary rather than the robust engine. Unavailable validation
  replays the complete canonical object, including human message/requirements;
- the declared `testthat` floor is `3.1.7`, the first release containing the
  mocked-binding API used by the deterministic failure rails.

Adversarial corrections:

- the first detection-mask fixture accidentally made `G=max(k)` and therefore
  forced a legitimate count boundary; rotating observed feature positions
  preserves the fixed global denominator while keeping intended rows interior;
- malformed GLM output was initially conflated with a fitted-probability
  boundary; nonfinite/shape failures now use numerical failure while only
  finite probabilities outside the inclusive frozen interval use the boundary
  code;
- an apparent zero-dispersion fixture has positive roundoff-scale `phi` and
  correctly reaches scalar singular covariance. A separate exact-zero fixture
  now locks the literal `phi<=0` numerical path;
- initial tests reconstructed delta effects through the implementation's own
  coordinates. Raw-design grouped GLM/covariance plus finite-difference
  gradient oracles now independently cover numeric and interaction inference,
  alongside mixed contrast failure, both count limits, inclusive probability
  floors, unclamped CI, acquisition strata, and rank deficiency.

Verification:

- focused quasibinomial suite → 9 tests / 77 expectations; focused robust and
  restricted-permutation suites remain 7 / 90 and 13 / 149; all pass;
- package-wide source suite → 190 tests / 1,478 expectations, zero failures,
  warnings, errors, or skips; a separate 100-fit synthetic raw-coordinate
  equivalence sweep had maximum relative error `3e-15`;
- three-agent final rereview → no remaining mathematical, contract, schema,
  provenance, dependency, or architecture finding; executable frozen contract
  remains 30/30 with unchanged v3 identities and `frozen_unrun` lifecycle;
- vignette-bearing source build succeeds; source-tarball
  `R CMD check --as-cran --no-manual` passes every substantive rail with only
  the expected new-submission note; usage scan and `git diff --check` pass;
- no candidate simulation ran; replicates `1-32`, development `33-64`,
  HarmonizR bytes, final panel, manifests, and confirmation remain sealed.

Exact next task → M13h: implement the exact public-panel handoff and
implementation-manifest constructor/validator, including executable-contract
quasibinomial panel coverage. Lock the manifest before opening candidate
replicates `1-32`.

Blockers → none.

### 2026-07-19 - M13h fail-closed panel + implementation-seal rails

Scope → implement the pre-allocation panel/manifest rails; adversarially prove
they cannot mint a winner or seal detached code/evidence.

Implemented:

- exact frozen public-panel shape reader reconstructs the internal candidate
  artifact and replays candidate-specific numerics; robust replay now rebuilds
  the deterministic HC3/CR2 artifact exactly;
- executable contract gains a quasibinomial panel fixture outside the normative
  source boundary; all frozen v3 identities remain unchanged;
- manifest core ports the frozen escaped-frame/self-hash protocol and requires a
  canonical complete inventory of package, M12/M13 contract, test, and harness
  sources with live-byte hashes;
- fixed subprocess registry derives pass status + positive counts from command
  exit/output rather than caller claims; focused receipt binds 443 expectations,
  stable pre/post source-inventory SHA-256, and the child package-namespace
  fingerprint;
- environment seal binds every live non-base namespace plus a child/parent
  fingerprint over all package-owned bindings; stored, tested, live source, and
  executable namespace drift all reject;
- candidate-evidence tables/materializers from the initial draft were deleted
  after review proved opaque aggregates could fabricate a winner from entirely
  unavailable outcomes. The sentinel association slot remains absent and
  `SEAL_READY` remains `FALSE` until semantic resolution exists.

Verification:

- focused panel/manifest/robust rail → 15 tests / 132 expectations; package-wide
  source suite → 198 tests / 1,520 expectations; zero failures, warnings, errors,
  or skips;
- registered child harness → 443 focused expectations with child/parent
  namespace and source-snapshot matches; frozen contract → 31/31 with unchanged
  protocol/contract/gate/effective-manifest hashes;
- vignette-bearing source build succeeds; source-tarball `R CMD check --as-cran
  --no-manual` passes every substantive rail with only the expected new-
  submission note; `git diff --check` passes;
- independent math + contract rereviews → no current finding after robust replay,
  source snapshot, namespace, receipt-schema, and fail-closed lifecycle fixes;
- no candidate simulation, candidate/development allocation, external data,
  HarmonizR bytes, public panel, implementation manifest, or confirmation result
  opened or created.

Exact next task → M13i: implement a companion resolver that hashes and validates
every exact input/result artifact, derives outcome bindings + six screening
metrics from frozen truth/opportunity allocations, resolves raw full-gate rows,
runs the frozen evaluator, and alone may authorize ranking/materialization. Then
flip the implementation-seal guard, lock the exact manifest, and only afterward
open candidate replicates `1-32`.

Blockers → none; authorization is intentionally red pending M13i.
