# ImputeFinder roadmap

Canonical method + gates → `PLAN.md`. This ledger owns session state, evidence, blockers, and exact next work.

## Current release path

- [x] M0 baseline + red contract + API decision
- [x] M1 input core + result skeleton
  - [x] M1a matrix validation + vector/design group alignment
  - [x] M1b `SummarizedExperiment` adapter + representation-preserving result skeleton
- [ ] M2 global absence + condition rescue
  - M2a deterministic rescue core + seed audit
  - M2b local RNG + row/column permutation invariants
- [ ] M3 manual classification + reconciliation
  - M3a post-seed statistics + cutoff validation + four states
  - M3b n-condition retention + groups + complete result contract
- [ ] M4 explicit profiles + plotting
- [ ] M5 automatic cutoff research + implementation
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
