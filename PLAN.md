# ImputeFinder plan - stable rules + evidence horizon

Status: canonical scientific contract and active Codex implementation plan

Reviewed baseline: clean `main` at `30618e1` (2026-07-17).

The first-release completion plan is finished. Its full historical specification
remains in git at `30618e1:PLAN.md`; completed session evidence remains in
`.agent/roadmap.md`; validation evidence remains under `dev/`. This successor
plan keeps the stable scientific contract self-contained while concentrating
active work on the next shippable horizon.

## 1. Operating model

`PLAN.md` owns stable scientific behavior, promoted successor scope, and gates.
`.agent/roadmap.md` owns live status, session evidence, blockers, and exact next
work. `.agent/memory.md` holds only durable facts that earn repetition beyond
code, tests, this plan, and git history.

At session start:

1. Read `AGENTS.md`, Sections 1-5 here, the active direction/milestone, the
   roadmap head + latest ledger entry, and memory.
2. Run `git status`; inspect relevant tracked source, tests, metadata, and recent
   commits.
3. Select one cohesive milestone/sub-milestone and define its failing test or
   explicit verification.
4. Re-check current primary sources before a new method, dependency, standard,
   or dataset decision.

At session close:

1. Run narrow tests, then the broadest proportional package/scientific gates.
2. Update the roadmap with evidence, failures, blockers, and exact next task.
3. Update this plan only for approved scope, improved gates, promotion state, or
   durable scientific decisions.
4. Make one scoped commit; record an intentional dirty tree as a blocker.

Ambiguous scientific behavior, estimands, or scope → stop and ask. A failed
research direction is a valid result. Preserve failed candidates and frozen
holdout misses rather than weakening gates.

## 2. Stable `rules_v1` contract

`classify_missingness()` is a released-shape, deterministic heuristic for
condition-aware quantitative-proteomics missingness. The stable rail classifies
and prepares data; it neither proves a causal mechanism nor normalises or
imputes values.

### 2.1 Input

Canonical input:

- ordinary numeric matrix data: features × samples;
- unnormalised log2 protein intensities; missing observations = `NA`;
- finite observed values, unique non-empty row/column names;
- condition labels aligned one-to-one with samples; the full comparative
  workflow expects at least two conditions;
- each condition contains at least one finite intensity somewhere;
- matrix + atomic group vector or aligned design data frame; or
  `SummarizedExperiment` + explicit `group_col`, with `assay` optional only when
  exactly one assay exists.

Reject `NaN`, infinities, duplicate/absent identifiers, missing group labels,
ambiguous assays, and alignment mismatches. Preserve finite zeros and negative
log2 values. Trust the documented scale; the package cannot infer whether data
were transformed correctly.

### 2.2 Algorithm

For every call:

1. **Global absence** - mark features missing in every sample `all_missing`;
   exclude them before later calculations and never seed them.
2. **Condition rescue** - compute each condition's minimum finite intensity
   before insertion. For each surviving feature fully missing in a condition,
   choose one condition sample uniformly and insert that minimum in exactly one
   cell. Default seed = `1L`; sampling uses stable feature/condition/sample-name
   order and local RNG scope. Log feature, condition, sample, old value,
   inserted value, and seed.
3. **Block statistics** - after rescue, calculate sample/observed/missing counts,
   missing fraction, arithmetic observed mean, and seeded flag for each
   feature-condition block.
4. **Missingness profile** - independently by condition, estimate complete and
   incomplete feature-mean densities on one grid. The explicit plotted profile
   is:

   ```text
   p_missing(x) = n_missing f_missing(x)
                  -----------------------------------------------
                  n_missing f_missing(x) + n_complete f_complete(x)
   ```

   Store raw statistics, grid, density metadata, support, and warnings. Plots
   consume stored data; rendered graphics never become algorithmic input.
5. **Cutoff** - use a finite named manual cutoff when supplied. Otherwise locate
   the right/bottom boundary of the dominant descending profile transition.
   Record method/version/quality; raise a structured condition-specific failure
   for weak, flat, or ambiguous evidence. A complete-only condition needs no
   cutoff.
6. **State** - for feature `i`, condition `g`:

   ```text
   missing_count == 0                         → complete
   mean_intensity < cutoff                    → MNAR
   mean_intensity >= cutoff and
     observed_count > sample_count / 2        → MAR
   otherwise                                  → insufficient
   ```

   Thus equality lies on the MAR side; strict majority =
   `floor(sample_count / 2) + 1`; sparse MNAR remains eligible.
7. **Reconciliation** - retain a feature exactly when every condition is
   `complete`, `MAR`, or `MNAR`, and at least one condition is `complete` or
   `MAR`. Drop precedence: `all_missing` → `insufficient:<conditions>` →
   `MNAR_all_conditions`.

Labels mean evidence-compatible **likely MAR/MNAR**, not identified truth.

### 2.3 Public API and result

Stable signature:

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

Stable `imputefinder_result` fields, in order:

1. `data` - retained, seed-modified data in the input's broad representation;
2. `classifications` - long block evidence/state/retention table;
3. `groups` - retained `MNAR`, `MAR`, `complete`, `MAR_or_complete` per condition;
4. `feature_status` - one original-feature retention/drop record;
5. `cutoffs`;
6. `cutoff_diagnostics`;
7. `profiles`;
8. `seed_log`, including rescues later dropped;
9. `groups_by_sample` aligned to output columns;
10. `call`.

Public helpers: `print()`, `summary()`, `plot_missingness()`.

### 2.4 Invariants and maintenance

- Input object remains unchanged. On returned retained rows, originally observed
  values remain exact; only logged rescue cells differ.
- Caller RNG kind/state, options, graphics state, and working directory remain
  unchanged.
- Named scientific results survive feature/sample/condition permutation.
- Matrix and `SummarizedExperiment` paths share one matrix core and agree.
- Manual/automatic failures are explicit; automatic selection never fabricates
  an endpoint or non-finite cutoff.
- The normative fixture, full unit suite, scientific regressions, 10,000 × 50
  performance gate, source check, and Bioconductor gates remain the v1 oracle.

`rules_v1` freezes scientific semantics, result fields/order/class, named
discrete outputs, masks, and rescue assignments. Versioned `rules_v1.x`
maintenance may improve messages, structured condition metadata, security, or
performance while preserving fixture behavior and documented meaning.

Comparison vocabulary for future tests:

- **exact** - names, masks, rescue assignments, discrete states, retained rows,
  result components/class, and calls;
- **canonical** - stored plot data, layer roles, labels, and diagnostics rather
  than environment-bearing plot object identity;
- **tolerance** - continuous fits, probabilities, intervals, and floating-point
  backend reductions;
- **threshold/noninferiority** - runtime, allocation, RSS, and scientific rates;
- **semantic parity** - equivalent scientific contents when representations
  intentionally differ.

### 2.5 Stable downstream handoff

The v1 workflow remains:

```text
classify unnormalised log2 data
  → use retained seed-modified data
  → user-chosen normalisation
  → split by condition
  → route condition MNAR vs MAR_or_complete sets
  → external imputation/analysis
```

Actual normalisation, imputation, differential testing, and method selection
remain outside the stable rail.

## 3. Successor brief

### 3.1 First persona and regime

Primary user: an R/Bioconductor analyst with bulk label-free protein-level
intensities, a biological condition, and optional batch and/or subject metadata.

Initial successor claims cover only validation-backed bulk LFQ regimes. DDA and
DIA are recorded as separate acquisition strata. TMT, precursor/peptide
hierarchies, single-cell proteomics, multi-omics, cross-study priors, and delayed
backends remain parked until their option cards are promoted. Stable
`classify_missingness()` retains its current documented input contract.

### 3.2 Program target and release floor

Build one opt-in sidecar workflow that attempts to answer:

1. Which biological and technical design effects are estimable, aliased, or
   inseparable from missingness?
2. Which classic states/retention decisions are stable to sampling, cutoff
   estimation, and declared assumption changes?
3. Does detection frequency differ by condition after declared design terms,
   and how is that distinct from observed-abundance change?

The guaranteed release floor is the sidecar seam plus validated typed-design/
estimability core. The A association panel, B, and C ship only when their frozen
development and confirmation gates pass; killed/parked tracks remain documented
negative results. The active program studies three directions:

- **A - design/confounding sentinel**;
- **B - robustness certificate**;
- **C - detectability contrast**.

### 3.3 Scientific boundary

Observed intensities plus missing indicators generally do not identify the true
causal missingness mechanism. Biological absence, abundance-dependent detection,
stochastic identification, interference, processing failure, and batch coverage
can overlap.

Therefore the successor:

- reports association, estimability, model-conditional evidence, stability, and
  unresolved cases;
- reserves causal/structural-absence claims for designs or auxiliary-variable
  assumptions sufficient to identify the stated estimand;
- treats abstention/`unavailable` as valid outcomes;
- keeps evidence separate from versioned decision policy;
- freezes protocols, metrics, thresholds, and holdouts before final candidate
  results;
- distinguishes synthetic truth, experimental abundance/detection truth, and
  naturally missing cells whose causal labels remain unknown.

## 4. Portfolio governance

Track states:

| State | Meaning |
|---|---|
| `stable` | Supported contract; compatibility rules apply |
| `active` | Approved implementation path with roadmap milestone + release gate |
| `candidate` | Bounded research proposal awaiting promotion evidence |
| `parked` | Valuable horizon; absent from current release/DoD |
| `killed` | Frozen evidence rejected the direction; negative result retained |
| `shipped` | Active gate passed and release evidence closed |

Current portfolio:

| ID | Direction | State | Dependency |
|---|---|---|---|
| v1 | Condition-aware hard-state classifier | `stable` | complete |
| A | Design/confounding sentinel core | `active` | M11-M12 |
| A-panel | Statistical missingness association | `killed` | M13: no full-gate passer |
| B | Robustness certificate | `active` | mandatory A design core + M12 |
| C | Detectability contrast | `active` | mandatory A design core + M12 |
| D | Censor-aware continuous evidence | `parked` | A design core + B + identification study |
| E | Leakage-safe workflow stress-test lab | `parked` | release 1 + bounded adapter contract |
| F | Distributional uncertainty handoff | `parked` | promoted D or E distributional path |
| G | Precursor/peptide/protein evidence graph | `parked` | modality-specific validation |
| H1 | Acquisition calibration | `parked` | B-D + transfer evidence |
| H2 | Value-of-information feedback | `parked` | H1 + action-specific evidence |
| I | Delayed/on-disk backend | `parked` | explicit stable/research adapter decision |
| J1 | Standards-aware QC exchange | `parked` | A + standards lifecycle decision |
| J2 | Cross-study atlas | `parked` | J1 + dataset governance + transfer evidence |

WIP cap: one active implementation milestone plus its validation work. Parallel
agents may research/review independent pieces inside that milestone. After the
M13-M15 studies close, M15P may activate at most one of D or E; every other card
stays parked until a later explicit decision. Parked cards never block the active
release Definition of Done.

## 5. Sidecar product and architecture

### 5.1 Golden path

M11 naming experiment selected this experimental API; only implemented
functions are exported:

```r
design <- missingness_design(
    sample_data,
    condition = "condition",
    nuisance = "batch",
    block = "subject",
    acquisition = "acquisition_mode"
)

analysis <- analyze_missingness(
    x,
    design,
    assay = NULL,
    base_fit = NULL,
    cutoffs = NULL,
    seed = 1L,
    modules = c("sentinel", "stability")
)

detectability <- test_detectability(
    analysis,
    x = x,
    contrast = list(condition = c(treated = 1, control = -1))
)
```

`analyze_missingness()` owns the original `x` during computation and either
attempts the classic fit or checks an optional `base_fit` for scientific and
structural compatibility with `x`, design, cutoff-source policy, and seed. This
check recomputes the classic path and compares the defined exact/canonical/
tolerance components; it establishes compatibility, not historical provenance.
Conflicting `base_fit`/`cutoffs` specifications fail. `test_detectability()`
rechecks the supplied `x` fingerprint because the sidecar does not retain numeric
intensities. A fit-only call cannot reconstruct dropped-row cells, the original
mask, or pre-rescue sample evidence; any limited fit-only accessor must return
structured `unavailable` for analyses needing them.

The frozen M11 comparator report is `base_fit_compatibility_v1`. Exact
components = result schema/class, returned data/representation, discrete block
and reconciliation evidence, groups, feature audit, rescue log, aligned sample
groups, and call. Canonical components = cutoff-source/diagnostic and stored
profile schemas, roles, labels, warnings, and other non-continuous values.
Tolerance components = finite doubles in classifications, cutoffs, diagnostics,
and profiles under absolute + relative `sqrt(.Machine$double.eps)` bounds;
special-value positions remain exact. Any non-empty difference category rejects
reuse. Manual sources are rerun as explicit cutoffs; automatic sources are
re-estimated. The supplied call is mirrored only to compare structure, and the
compatible supplied fit is retained unchanged; neither fact proves its history.

Portable abstention records use class `imputefinder_unavailable` and fields
`status`, `quantity`, `code`, `message`, `requires`. Fit-only requests for
dropped-row cells, the original mask, or pre-rescue evidence return code
`original_input_required` and require original `x` + typed design.

The M11 module router accepts a unique subset of `sentinel` and `stability`;
unselected slots are `NULL`. Selected `sentinel` stores the input-first
`pre_rescue` record. Until M12/M14 validation and implementation close, selected
`stability` returns `module_pending_validation` requiring
`robustness_certificate_v1`. The default requests both, making incomplete
research capability explicit without fabricating output.

For matrix input, the typed design supplies the aligned condition. For
`SummarizedExperiment`, `assay` follows the stable routing rule and the design's
named condition column must exist in, align with, and equal `colData(x)`; users
with external-only condition metadata pass the selected assay as a matrix.

### 5.2 Execution

```text
x + typed design
  ├─→ immutable pre-rescue evidence + original mask
  ├─→ stable `rules_v1` attempt: result | structured failure
  ├─→ A: estimability/association diagnostics
  └─→ B: named perturbation runs of the stable evidence/decision path
          ↓
     companion analysis + provenance
          + fingerprint-matched x + contrast
          ↓
     C: standalone detection + observed-abundance result
```

Stable output stays untouched. Shared A/B evidence lives in an experimental
companion class, provisionally `imputefinder_analysis`:

```text
classic
    Exact `imputefinder_result` on success; otherwise a portable structured
    failure record with stage/class/message/call. Never a partial v1 result.
spec
    Sidecar schema/module/policy/software identifiers and lifecycle status.
design
    Aligned supplied roles, declared terms, estimability, and unavailable roles.
input
    Dimensions, ordered names, deterministic fingerprint, original mask,
    scale/acquisition declarations; no duplicate full numeric matrix.
sentinel
    A outputs or NULL. Input-first objects begin with a schema-tagged
    `pre_rescue` record: feature-condition + sample counts, missing fractions,
    and observed means computed before any rescue/classic attempt.
stability
    B outputs or NULL, separated by uncertainty family.
provenance
    Call, seeds/streams, hashes, warnings, failures, assumptions, training scope.
```

The caller owns numeric `x`; serialized sidecars retain computed evidence and
the original mask, not a second full intensity matrix. Re-running a new module
requires `x` again and verifies its fingerprint. The M11 fingerprint is SHA-256
over the versioned `matrix_core_be_v1` protocol: storage tag, big-endian
dimensions/observed values, length-prefixed UTF-8-or-byte identifiers, and the
column-major missing mask. The inline mask uses zero-padded LSB-first base-R
bit packing (`base_packbits_lsb0_v1`). Protocol identifiers travel with the
record; unknown identifiers and missing/mismatched artifacts fail explicitly.

Pre-rescue evidence is constructed before the classic attempt. Supported A/C
work may therefore continue after an automatic cutoff failure; B records the
failure frequency rather than requiring a fabricated fit. Modules that truly
depend on a successful classic result return structured `unavailable`.
`pre_rescue_evidence_v1` preserves zero-observation blocks as explicit counts
with `NA` observed means; it stores summaries rather than cell intensities.

`test_detectability()` returns an immutable standalone
`imputefinder_detectability` object referencing the sidecar schema/input
fingerprint, design hash, contrast, and its own provenance. Multiple contrasts
therefore neither mutate nor ambiguously overwrite the sidecar.

### 5.3 Design object

One row per sample, names exactly aligned to `x`. Roles:

- required `condition`;
- optional declared `nuisance` terms such as batch/plate/run order;
- optional `block` for subject/pair/repeated measures;
- optional acquisition mode;
- optional explicitly declared interactions;
- resampling unit derived from block/biological replicate semantics.

The sidecar audits every supplied and estimand-required role. Absent optional
metadata is recorded as unavailable and constrains claims; it does not invalidate
an otherwise supported analysis. Automatic interaction search, sample exclusion,
batch correction, or data mutation is outside A-C.

The design/estimability core in A is mandatory shared infrastructure for B/C.
A's statistical missingness-association panel remains an independently
falsifiable release track; its failure cannot remove the validated design core.

### 5.4 Lifecycle and exactness

- Research class/schema is experimental from M11; schema ID and upgrade/error
  behavior ship with the first object.
- M11 identifiers live only in internal constants + sidecar `spec`; v1 receives
  no new field, attribute, or class.
- Inline artifacts only for the first sidecar schema. External content-addressed
  references await a separate lifecycle design.
- Exact/canonical/tolerance/threshold/semantic comparisons follow Section 2.4.
- Sidecar modules read stored scientific data, never rendered plot coordinates.
- Local RNG scopes and stable instance IDs preserve caller state and auditability.

## 6. Active direction A - design/confounding sentinel

### 6.1 Hypothesis and output

Retaining declared design context will reveal missingness associations, aliasing,
and non-estimable biological/technical separation before users over-interpret a
condition label.

Report, without automatic correction:

- design rank, contrast estimability, exact aliasing, singleton levels, and
  condition-nuisance overlap;
- sample detection rate, intensity support, feature overlap, and pre-rescue
  coverage;
- condition × declared nuisance/block coverage tables;
- missingness association with declared main effects/interactions;
- practical overlap/uncertainty diagnostics for near-confounding;
- unsupported acquisition/metadata claims and actionable wording.

Language: “associated with,” “aliased with,” “cannot be separated under this
design.” The sentinel never claims which factor caused missingness.

### 6.2 Method constraints

- Algebraic design/contrast checks precede statistical association.
- Statistical terms are user-declared; multiplicity and low support are visible.
- Pre-rescue evidence drives QC; rescue anchors remain separately marked.
- A owns static summaries. Perturbation-based influence belongs to B.
- DDA and DIA evidence stays stratified until validation supports pooling.

### 6.3 Gate

M12 freezes numeric entries in the gate registry before implementation results:

Mandatory design-core gates:

- exact detection of constructed rank deficiency, perfect aliasing, and
  non-estimable contrasts;
- correct block/paired and unequal-replication accounting;
- label-preserving order/re-encoding invariance;
- no mutation of input, classic result/failure, or global state.

Independently disposable association-panel gates:

- predeclared false-flag/error behavior on null simulations;
- one frozen multi-batch public case meeting its predeclared aliasing/coverage
  finding and wording rubric, with zero causal attribution or unsupported
  correction claims;
- calibrated uncertainty/multiplicity behavior under the frozen alternatives.

## 7. Active direction B - robustness certificate

### 7.1 Three distinct quantities

Every output is assigned to exactly one family:

1. **Sampling stability** - biological/block resampling and leave-one-unit-out.
2. **Estimator uncertainty** - automatic profile/cutoff feature resampling and
   structured cutoff failure frequency.
3. **Assumption/policy sensitivity** - fixed versus re-estimated cutoff, cutoff
   sweep, and predeclared summary/rule scenarios.

The sidecar never pools these into a single “confidence” score. Resampling
frequencies are not posterior probabilities. Classic rescue sample placement is
classification-invariant because the inserted value/count are unchanged;
downstream placement sensitivity belongs to parked workflow Direction E.

### 7.2 Perturbation contract

A manifest records:

- perturbation family, ID, seed stream, source sample/feature IDs;
- draw-instance IDs, multiplicities or weights, and resampling unit;
- fixed/re-estimated parameters and policy version;
- success/failure class and canonical alignment back to original names.

Internal weighted/instance evidence handles bootstrap multiplicity; duplicate-
named matrices never pass through the v1 public validator.
Feature resampling uses cluster/block structure when declared dependence makes
ordinary exchangeability implausible. Its output is finite-dataset stability
unless a frozen sampling model justifies inferential coverage.

### 7.3 Outputs

By feature-condition and retained feature:

- classic intensity and strict-majority margins;
- sampling state/retention frequencies + transition counts;
- cutoff empirical resampling distribution/quantile range + failure frequency;
- assumption scenario matrix and decision boundaries;
- influential biological/block units from leave-one-unit-out;
- structured `unavailable` when support/design cannot estimate the quantity.

Any `stable`/`fragile` display policy is derived only after M12 calibration and
retains the continuous components beside the label.

### 7.4 Gate

- Where the frozen sampling model supports it, cutoff quantile-range coverage of
  the generating boundary + false-confidence thresholds from M12 pass;
  otherwise the output stays explicitly descriptive.
- Flat/no-cliff controls trigger weak-identification/failure in the cutoff panel;
  sparse designs return `unavailable` for unsupported quantities. Other panels
  retain their independently measured stability.
- Perturbation instance order, input order, and seed stream reproduce canonical
  named outputs.
- Caller RNG/global state is unchanged.
- Stable default calls incur zero sidecar work; v1 performance gates still pass.

## 8. Active direction C - detectability contrast

### 8.1 Estimands

Condition-associated non-detection can carry predictive or inferential signal.
Preserve two distinct effects per feature:

1. **Model-adjusted detection contrast** - condition difference in
   probability/odds of detection after declared nuisance/block structure.
2. **Conditional observed-abundance contrast** - condition intensity difference
   conditional on detection, with its selection limitation explicit.

The first release reports both components separately. A combined decision is
allowed only if M12 predeclares the estimand and M15 validates its error control.
Neither component alone proves structural biological absence.

### 8.2 Candidate study

Before implementation, compare the simplest viable exact, penalized,
beta-binomial/quasi-binomial, block-aware, or Bayesian candidates needed for:

- unequal groups and paired/repeated designs;
- overdispersion/correlation across replicates;
- complete/quasi separation;
- effect size + interval, rather than only a p-value;
- thousands of correlated feature tests.

A current primary-source/tooling review precedes the decision. Package
availability alone does not choose the method.

### 8.3 Support, multiplicity, and abstention

- Contrast must be algebraically estimable under A.
- M12 freezes per-group/block evidence floors by estimand.
- Sparse or separated cases use the selected valid estimator or return
  `unavailable`; infinite silent coefficients are invalid.
- Family definition and FDR/false-sign control are explicit; adjusted and raw
  evidence remain distinguishable.
- Naturally missing cells receive no invented abundance.
- Batch-exclusive coverage yields non-estimability, not a biological hit.

### 8.4 Outputs and gate

Standalone output: analysis/input/design references, design/contrast,
model-adjusted detection contrast + interval/test evidence, conditional
observed-abundance contrast + interval/test evidence, support, separation
status, multiplicity adjustment, assumptions, abstention reason, and provenance.

Gate registry must include:

- null rejection/FDR and false-sign behavior under correlated features;
- effect bias and interval coverage;
- power/precision-recall on frozen alternatives;
- paired, unequal, separated, low-support, and confounded cases;
- experimentally identified reference detection/recovery quantities from
  dilution, spike-in, mixed-species, or replicate data;
- untouched development test data and named-order invariance;
- wording audit separating differential detection, observed abundance, and
  causal absence.

## 9. Initial validation charter

### 9.1 Scope for M12-A-C

Initial generator/data scope:

- bulk LFQ protein-level matrices and `SummarizedExperiment` equivalents;
- 2-4 conditions, balanced/unequal replication, small-to-moderate cohorts;
- independent and paired/subject-blocked designs;
- declared batch/run effects, partial and perfect confounding;
- intensity-independent, monotone abundance-associated, mixed, blockwise,
  batch-associated, structural-off-compatible, and no-cliff patterns;
- outliers, heteroscedasticity, differing intensity support, and low evidence;
- DDA and DIA as distinct strata.

TMT, single-cell, peptide hierarchy, multi-omics, transfer priors, and delayed
storage add their own validation versions only after promotion. M12 does not
simulate contracts that have not been designed.

### 9.2 Evidence tiers

| Tier | Evidence | Claim supported |
|---:|---|---|
| 0 | Existing normative fixture + v1 tests/reports | Compatibility |
| 1 | Deterministic latent-truth simulation | Method behavior/failure |
| 2 | Several semi-synthetic mask generators | Misspecification robustness |
| 3 | Controlled dilution/spike-in/mixed-species/replicates | Abundance/detection recovery |
| 4 | Public bulk LFQ studies with design/batch metadata | External behavior/shift |

Tier 1-2 causal labels are generator-specific. Tier 3 truth is used only for
quantities the experiment identifies. Tier 4 naturally missing cells remain
causally unlabeled.

### 9.3 Split and holdout ownership

Each dataset card declares deployment target, unit, related samples, licence,
checksum, acquisition, truth, exclusions, and split hash.

Separate data roles:

- training/resampling data for fit/tuning;
- frozen development test data opened for a milestone's candidate decision;
- independently sourced, sealed release-confirmation data reserved for M15P.

All learned preprocessing fits inside training data. Related subject/technical
units stay together when deployment requires unseen units. A sealed-holdout miss
is a failed claim. A later holdout requires a materially revised hypothesis or
method plus independently sourced data, never another split chosen from the
exhausted dataset.

### 9.4 Gate registry

Before candidate results, freeze for every gate:

```text
claim + metric + unit + comparator + strata + sample size
+ uncertainty method + numeric threshold + failure treatment + data/protocol hash
```

Report distributions and strata. Reconstruction error alone cannot support a
mechanism, contrast, or workflow recommendation.

### 9.5 Evidence package

Each method milestone produces:

- executable `dev/` script;
- frozen protocol + gate registry;
- concise tracked report with all candidates/negative controls/failures;
- dependency/session provenance and data manifest;
- lightweight deterministic regression distilled into routine tests.

Long stochastic/external runs stay outside routine package checks; summaries and
hashes remain tracked.

## 10. Active milestones

Dependency:

```text
M11 sidecar seam
  → M12 frozen initial validation
      → M13 design core + sentinel study ─┬→ M14 robustness study
                                          └→ M15 detectability study
M13 + M14 + M15 → M15P release confirmation + promotion review
```

### M11 - Stable seam + minimum sidecar schema

- [x] Add internal stable method/profile/cutoff/rescue/policy identifiers and
      sidecar `spec`; preserve every v1 field, attribute, class, and output.
- [x] Define minimum `missingness_design` roles, validation, and lifecycle.
- [x] Define `imputefinder_analysis` schema with inline original mask,
      fingerprint, provenance, empty module slots, `classic` result/failure sum
      type, and explicit upgrade/error behavior.
- [x] Implement input-first construction + optional compatible `base_fit`;
      pre-rescue evidence survives classic cutoff failure and fit-only
      limitations return structured `unavailable`.
- [x] Add differential fixtures for exact/canonical/tolerance/performance
      categories.
- [x] Run the API naming experiment and freeze the experimental public/internal
      boundary before export.
- [x] Research and record fingerprint/mask implementation choice before use.

Gate: all v1 unit/scientific/package gates pass; existing performance thresholds
remain satisfied; no component/attribute/class drift; sidecar round-trip and
mismatched-input failures pass.

### M12 - Frozen validation + gate registry

- [x] Extend generators only across Section 9.1 initial scope.
- [x] Select the minimum controlled + public datasets needed for A-C; record
      cards, licences, checksums, and data roles.
- [x] Freeze development and sealed release-confirmation identities.
- [x] Freeze numeric gate registry for A-C before candidate implementation.
- [x] Freeze candidate protocols for A association and C contrast estimators,
      plus B resampling units/counts/replacement, weighting, seed streams,
      cutoff-range definition, failure handling, and display calibration.
- [x] Assign sealed confirmation artifacts per track/claim. Declare shared
      artifacts up front; their synchronized opening retires them for every
      linked claim.
- [x] Add null/confounded/leakage/low-support negative controls.

Frozen M12 allocation = simulation replicates `1-32` candidate
selection/calibration, `33-64` untouched development. Registry = 27 retained
claims + 66 numeric gates bound to exact data/protocol hashes. Track protocols
`m12_a_candidate_protocol_v2`, `m12_b_perturbation_protocol_v2`, and
`m12_c_candidate_protocol_v2` remain `frozen_unrun`; external artifacts remain
`metadata_only`. The pre-result v2 literature correction requires full-design
identity-working-model CR2 with computed Satterthwaite df `>=5`, compatible
exchangeability/variance groups for residual permutation, and treats the global
quasibinomial model as empirical rather than MSqRob-equivalent.
`dev/m12-validation-contract.md` owns thresholds, hashes, audit evidence, and
the exact candidate/perturbation summary.

M13's pre-result Track A v4 addendum binds the immutable M12 hashes while
freezing one coefficient per hypothesis, acquisition-specific rebuilt models,
support cells, constrained-null permutation, canonical mask/design seeds,
typed per-component numeric support, result schema, candidate-agnostic
opportunity denominators, four fixed alternative targets, and exact
candidate/development/confirmation stage seals. V4 additionally freezes the
128-target coverage fraction and 9,999-draw scenario-stratified whole-replicate
bootstrap mechanics, preserving within-replicate feature-module dependence.
`dev/m13-association-contract.md` owns these superseding details; its aggregate
hash is `422954c88ed7683f58ac296cb4d9e785d61c3f80a276914235d4b285bbc16547`
(protocol `e9d3b3799bb160d65feb0ae63b2e49a92bb4e8e005b89232ada96b83cfabed27`;
gate registry `65f01b20e74b216399e17433df393f32e9bda189f73992ab419069cf7d5c343b`).

The production `association_preparation_v1` rail now binds the packed global
mask/design identity to the exact response, acquisition-specific rebuilt
cores, coefficient hypotheses, estimability, and biological-unit support.
The three internal candidate engines now implement the frozen OLS-HC3/
full-design identity-working-model CR2/Satterthwaite path, studentized
restricted Freedman-Lane with exact/Monte Carlo constrained-null maps, and the
empirical independent-unit grouped quasibinomial path. They bind exact
unavailable handoff, per-stratum Holm families, preparation-bound artifact
validation, and candidate-specific numerical provenance replay. The exact
final-panel shape reader and fail-closed implementation-manifest core now bind
fixed passing-test receipts, stable source snapshots, the full live dependency
set, and child/parent namespace fingerprints. The semantic evidence resolver
now owns exact input/result reads, generator and engine replay, outcome audits,
screening/full gates, ranking, and self-hash; production accepts no caller
provenance callback and fixture evidence is non-authoritative. The fixed study
harness has no development/confirmation command. Implementation manifest
`d821cae0...c1cba4` locked before allocation; all 416 canonical inputs and 1,248
candidate results were then generated. Independently replayed evidence
`a5a242c3...60cbc6` resolved `no_winner`: all candidates failed screening, so
the full-gate/ranking sets are empty. The association panel is `killed`; public
materialization, development replicates `33-64`, HarmonizR bytes, and
confirmation remain unopened. The mandatory design core survives independently.

Gate: protocol hashes exist; every A-C claim has an executable numeric gate and
failure rule; v1 behavior/limitations reproduce on the expanded harness.

### M13 - A: design core + confounding sentinel study

- [x] Implement typed design alignment and algebraic estimability/aliasing.
- [x] Implement pre-rescue sample/condition/nuisance coverage summaries.
- [x] Freeze the executable M13c association addendum before candidate results.
- [x] Implement only predeclared association terms and uncertainty.
- [x] Retain exact candidate diagnostics/evidence; omit the public panel,
      print/summary, and wording surface after the killed disposition.
- [x] Run candidate A gates and retain every failure. No candidate passed
      screening, so development was correctly not opened and the association
      panel is `killed`; confirmation remains sealed for synchronized M15P.

Gate: the algebraic design/estimability core must pass because B/C depend on it.
The association panel closes through a recorded disposition; only a passing
panel becomes a confirmation candidate. Neither component mutates data or makes
causal attribution/automatic correction.

M13 candidate evidence → `dev/m13-association-validation.md`. The killed panel
is not a failed mandatory-core claim and does not block B/C.

### M14 - B: robustness certificate

- [x] Implement weighted/instance perturbation manifests + deterministic streams.
- [ ] Implement sampling, cutoff-estimator, and assumption/policy panels
      separately.
- [ ] Add margins, state transitions, retention/cutoff summaries, failures, and
      influence outputs.
- [ ] Calibrate optional display labels through risk-coverage analysis.
- [ ] Run the frozen B study; retain failures, distill passing regressions, and
      assign `confirmation_candidate`, `killed`, or `parked`.

Gate: the study and disposition are complete. Only a Section 7.4 pass becomes a
confirmation candidate; no disposition permits a pooled confidence score or
global side effect, and stable calls retain existing performance.

### M15 - C: detectability contrast

- [ ] Complete frozen estimator comparison; select simplest passing method.
- [ ] Implement estimability/support/separation/abstention behavior.
- [ ] Implement detection + conditional-observed-abundance contrasts and
      multiplicity.
- [ ] Add effect/interval diagnostics and wording tests.
- [ ] Run the frozen C study; retain disagreement/failures and assign
      `confirmation_candidate`, `killed`, or `parked`.

Gate: the study and disposition are complete. Only a Section 8.4 development
pass becomes a confirmation candidate. No result is called structural absence
or unconditional differential abundance.

### M15P - Release confirmation + next promotion

- [ ] Open only confirmation-candidate artifacts under their per-track registry;
      open shared artifacts synchronously and retire them for every linked claim.
- [ ] Adversarially map each public claim to code, test, dataset, and interval.
- [ ] Complete docs/vignette/NEWS, privacy/redaction, accessibility, and failure
      UX reviews.
- [ ] Run unit tests, build, examples, vignette, `R CMD check --as-cran`,
      `BiocCheck`, minimal-library/optional-dependency-absence, and performance
      gates.
- [ ] Mark the A association panel, B, and C `shipped`, `killed`, or still
      `candidate` independently; preserve the mandatory design core separately.
- [ ] Use evidence/user need to promote at most one of D or E, or pause.

Gate: Section 14 passes for shipped tracks. A failed optional track does not
weaken or block a validated one.

## 11. Parked option cards

These are innovative horizons, not active commitments or release gates.
Before promotion, every card receives a Section 9.4 registry whose numeric gate
and comparator directly test its stated hypothesis.

### D - Censor-aware continuous evidence

Hypothesis: a monotone detection/selection model with explicit sensitivity and
abstention can outperform one hard cutoff in low-information compatible strata.

- Candidate outputs: abundance-associated detection evidence, residual pattern,
  calibration, support, drift, sensitivity; “structural” only under designs or
  auxiliary-variable assumptions sufficient to identify the stated estimand.
- Compare logistic detection curves, selection likelihoods, hierarchical factor
  models, and simple calibrated baselines.
- First gate: measurable calibration/selective-risk gain over the v1 hard cutoff
  in a frozen compatible low-information stratum, external-shift refusal, and v1
  noninferiority on its home scenarios.
- Promotion requires the A design core, B shipped, and a focused identification
  study.

### E - Leakage-safe workflow stress-test lab

Hypothesis: dataset/task-specific comparison of complete workflows can beat a
universal imputation recipe.

- Candidates include no-imputation/direct models and external imputation
  adapters; evaluate effect bias, coverage, FDR/power, rank stability, runtime,
  and cell error by evidence tier.
- Training-fold-only fit applies to cutoff, normalisation, masking, feature
  selection, batch handling, and imputation.
- Registered adapters can be audited. Arbitrary callbacks are `unverified` and
  cannot support package recommendations without restricted/replayable proof.
- First gate: nested held-out regret (excess frozen loss versus the best member
  of the frozen candidate set) + intentional leakage detection + abstention when
  candidates tie. Masked observed cells never become truth for naturally
  missing cells; one benchmark reports FDP/TPR, while FDR/power require repeated
  studies or simulations.

### F - Distributional uncertainty handoff

Hypothesis: draws, precision weights, or direct missingness-aware contracts yield
better calibrated downstream inference than one completed matrix.

- Preserve original mask, draws/distribution, typed predictive/credible/
  confidence intervals or estimator SEs, method/seed/scope, within/between-
  imputation variance, and selective unfilled cells.
- Support Rubin-style pooling only for proper imputations + a congenial analysis;
  mark one-matrix views explicitly lossy.
- First gate: latent-value coverage in simulation and reference-value/effect
  coverage in experiments, accounting for reference measurement error, by
  intensity and mask stratum, versus frozen single-matrix and no-imputation
  baselines.
- Depends on a promoted distributional D/E path.

### G - Precursor/peptide/protein evidence graph

Hypothesis: lower-level support, sibling covariance, and mapping ambiguity improve
protein missingness decisions.

- Prototype `QFeatures`/long-table input in `dev/`; preserve shared/ambiguous
  peptide edges and singleton warnings.
- Optional transcript evidence requires paired-data and negative-transfer tests.
- First gate: held-out or genuinely external lower-level reference evidence
  improves detection decisions without degraded calibration/retention.
- Modality-specific validation precedes package-facing API/dependency decisions.

### H1 - Acquisition calibration

Hypothesis: compatible reference runs can stabilize small studies.

- Transfer artifact records acquisition, intensity support, training population,
  curve/uncertainty, version, and drift.
- First transfer gate: leave-lab/instrument-out benefit + refusal outside support.

### H2 - Value-of-information feedback

Hypothesis: an explicit utility/cost and acquisition-outcome model can rank a
useful next acquisition action.

- First gate: hidden-replicate replay beats random/coverage-only action
  for the same hidden action type. Biological-replicate recommendations require
  their own utility/cost and acquisition-outcome model.
- Depends on validated H1 calibration in the target acquisition regime.

### I - Delayed/on-disk backend

Hypothesis: an internal block contract can preserve scientific semantics with
bounded memory on large matrices.

- Promotion first decides: additive stable adapter in `rules_v1.x` or research-
  rail-only input.
- Candidate: `DelayedArray`; `HDF5Array` only as optional benchmark backend.
- Serial named rescue plan stays independent of chunks/workers.
- First gate: semantic matrix parity, no full realization, bounded RSS, backing
  file immutability, chunk/worker invariance within numeric tolerance.

### J1 - Standards-aware QC exchange

Hypothesis: versioned SDRF/mzQC exchange can preserve package evidence losslessly
across conforming tools and improve auditability.

- Static offline reports redact cell values/local paths by default.
- First gate: round-trip/version/controlled-vocabulary conformance without loss
  of package-native evidence; experimental metrics remain namespaced until
  accepted by the relevant standard.

### J2 - Cross-study atlas

Hypothesis: curated experimental truth can support guarded transfer across
compatible studies.

- Dataset cards: source, licence, checksum, acquisition, truth, split, claim,
  exclusions.
- Atlas prior uploads no user data and requires leave-study/lab/instrument-out
  validation + shift refusal.
- Depends on J1 provenance exchange plus explicit dataset governance.

## 12. Risk register

| Risk | Consequence | Control / kill rule |
|---|---|---|
| Mechanism non-identifiability | False causal certainty | Association/model-conditional wording, sensitivity, abstention |
| V1 drift | Existing analyses break | Dual rail + differential oracle + internal/sidecar-only IDs |
| Fit-only reconstruction assumption | Invalid QC/resampling | Input-first sidecar; structured unavailable |
| Confounding attribution | Batch called biology | Algebraic estimability + “inseparable” output |
| Bootstrap false precision | Fragile decisions look stable | Three separate panels + null calibration |
| Multiplicity/low support | Detectability false discoveries | Frozen family/FDR/support/abstention rules |
| Generator favoritism | Method wins its simulated world | Several generators + controlled + external truth |
| Holdout shopping | Inflated evidence | Sealed independent confirmation; a miss remains failure |
| Scope explosion | No shippable release | A-C only + WIP cap + one-card promotion checkpoint |
| Dependency expansion | Fragile Bioconductor release | `dev/` research → boundary decision → `Suggests`/companion/import |
| Artifact/privacy leakage | Nonportable or sensitive result | Inline mask only, caller-owned values, redacted reports |
| Acquisition overgeneralisation | Invalid TMT/single-cell claims | Initial LFQ persona + modality-specific promotion |

## 13. Research anchors

Sources motivate questions and gates; they do not preselect implementations.
Re-check current versions, licences, follow-up evidence, and official guidance
at the implementing milestone.

Core scientific boundary and detection models:

- [Nonignorable missingness in label-free proteomics](https://doi.org/10.1214/18-AOAS1144)
  motivates mechanism-aware interval evaluation and cautions against single
  imputation as complete truth.
- [protDP](https://doi.org/10.1093/bioinformatics/btad200) estimates
  intensity-dependent detection probability beyond a pure random/censored split.
- [Sequential identification of nonignorable
  mechanisms](https://doi.org/10.5705/ss.202016.0328) shows that observed data
  alone require identifying restrictions for MNAR models.
- [limpa](https://smythlab.github.io/limpa/) uses a detection-probability curve,
  recovers complete protein estimates with uncertainty, and carries precision
  information rather than treating filled values as observed.
- [Pirat](https://doi.org/10.1093/biostatistics/kxaf006) fits an assumed global
  instrument-censoring model using peptide covariance and optional transcript
  evidence.
- [msBayesImpute](https://doi.org/10.1038/s42004-026-02106-3) combines
  protein-specific dropout with Bayesian factorisation and current controlled
  benchmark designs.

Active A-C rationale:

- [Batch-associated missing values](https://doi.org/10.1093/bib/bbaf168)
  motivate design/coverage diagnostics before interpretation or imputation.
- [Missingness as disease-discriminative signal](https://doi.org/10.1093/bioinformatics/btz898)
  shows condition-associated missingness was disease-discriminative in one
  small longitudinal rheumatoid-arthritis SWATH cohort with disease/batch
  overlap and ad hoc count adjustments; it is neither proof of causation nor
  validation of a confounder-adjusted association estimator.
- [MSqRob's missing hurdle](https://doi.org/10.1021/acs.analchem.9b04375)
  directly motivates paired detection/intensity components for label-free DDA,
  with associational interpretation and empirical FDP/power evaluation on a
  controlled spike-in benchmark. Its grouped response is peptides within each
  protein, not one all-protein detected fraction per sample.
- [Restricted permutation](https://doi.org/10.1016/j.neuroimage.2014.01.060)
  requires null-invariant exchangeability blocks compatible with distinct
  variance groups; the later [robust-W comparison](https://doi.org/10.1016/j.neuroimage.2019.116030)
  supplies asymptotic heteroscedastic support but warns that small-sample
  residual-permutation inflation can persist.
- [CR2/Satterthwaite inference](https://doi.org/10.1080/07350015.2016.1247004)
  makes reference df leverage-dependent; its
  [corrigendum](https://doi.org/10.1080/07350015.2023.2174123) limits the
  absorbed-fixed-effect shortcut to ordinary unweighted LS with identity
  working model + rank conditions.
- [Cluster-bootstrap consistency](https://doi.org/10.1016/j.jmva.2012.09.003)
  supports whole-cluster resampling while documenting sparse-binary failure
  modes at small cluster counts.
- [Downstream-centric imputation criteria](https://doi.org/10.1021/acs.jproteome.3c00205)
  show that reconstruction error could mislead downstream method choice in the
  evaluated benchmarks.

Parked workflow/uncertainty/hierarchy rationale:

- [mi4p](https://doi.org/10.1371/journal.pcbi.1010420) propagates multiple-
  imputation variability into moderated differential analysis.
- [PEPerMINT](https://doi.org/10.1093/bioinformatics/btae389) uses a peptide-
  protein graph; predicted uncertainty ranked errors and supported filtering,
  without establishing nominal interval coverage.
- [Doubly robust post-imputation inference](https://doi.org/10.1214/25-AOAS2012)
  motivates bias-aware inference around flexible imputers under its MAR,
  positivity, independence, rank, and nuisance-rate assumptions.
- [Proteomics workflow interactions](https://doi.org/10.1038/s41467-024-47899-w)
  motivate whole-workflow rather than isolated-step evaluation.
- [DataSAIL](https://doi.org/10.1038/s41467-025-58606-8) shows in molecular-
  property and drug-target tasks that random splits can overestimate OOD
  performance when similar entities cross folds, motivating deployment-matched
  grouping/similarity audits where that risk applies.
- [`QFeatures`](https://bioconductor.org/packages/release/bioc/html/QFeatures.html),
  [SDRF-Proteomics v1.1.0](https://sdrf.quantms.org/specification.html), and
  [HUPO-PSI mzQC](https://hupo-psi.github.io/mzQC/) are current candidate
  infrastructure/standards for parked G/J1/J2.

## 14. Active-release Definition of Done

Evaluate only `stable`, `active`, or explicitly promoted tracks. Parked/killed
cards close through recorded disposition and never count as unresolved P0/P1.

- [ ] Every v1 exact/canonical/tolerance/performance/package gate passes.
- [ ] V1 result fields, attributes, class, discrete decisions, rescue assignments,
      and scientific semantics remain compatible.
- [ ] The sidecar starts from original `x` + typed design, preserves the original
      mask without duplicating the numeric matrix, and rejects mismatched input.
- [ ] Supplied/estimand-required roles align; unavailable optional metadata is
      recorded and constrains claims.
- [ ] The mandatory design core reports estimability/aliasing and supports B/C;
      any shipped A panel reports declared association without causal
      attribution, mutation, or automatic correction.
- [ ] Any shipped B reports sampling stability, cutoff-estimator uncertainty,
      and assumption/policy sensitivity separately, with structured unavailable.
- [ ] Any shipped C reports detection and conditional-observed-abundance
      estimands separately, with support, separation, multiplicity,
      FDR/false-sign, and abstention.
- [ ] All learned operations fit only on their declared training data.
- [ ] Gate registries, protocol/data hashes, negative controls, candidate
      failures, and sealed-holdout misses remain tracked.
- [ ] Claims name evidence tier, acquisition/design scope, estimand, uncertainty,
      and limitations.
- [ ] Static output is offline, auditable, accessible, and redaction-aware.
- [ ] Routine tests hold compact scientific regressions; long external/stochastic
      reports remain reproducible under `dev/`.
- [ ] Unit tests, source build, examples, vignette, `R CMD check --as-cran`,
      `BiocCheck`, minimal-library optional-dependency, and performance gates pass.
- [ ] Documentation explains non-identifiability, fit-only limits, abstention,
      design confounding, and exact downstream handoff without overclaim.
- [ ] The A association panel, B, and C each has an independent shipped/killed/
      candidate disposition; the validated design core is tracked separately.
- [ ] Promotion review activates at most one next card or records a pause.
- [ ] `.agent/roadmap.md` has no unresolved P0/P1 item for the shipped active slice.
