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
- [ ] M5 automatic cutoff research + implementation
  - [ ] M5a deterministic simulations + predeclared benchmark protocol
  - [ ] M5b candidate benchmark + method decision
  - [ ] M5c pure detector + quality/failure contract
- [ ] M6 public integration + obsolete dependency cleanup
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
