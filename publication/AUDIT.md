# Fixed publication audit

Scope = the ten rows fixed in `PROTOCOL.md`; article = bounded software/evaluation claim, not universal imputation superiority. Findings outside those claims = deferred below.

| ID | Ruling | Evidence + adjudication |
|---|---|---|
| A1 | Pass | `manuscript.Rmd` defines cutoff labels as operational, distinguishes incomplete-feature density from cell missingness, describes method v3 safeguards, cites DEP/DAPAR/msImpute/sector routing + prior Fowler use. No causal/novelty/superiority claim. |
| A2 | Pass | `benchmark.R` records independently generated cell causes before package invocation; eight scenarios include covariate-MAR confounding + inverse MNAR. Forced structural absence has its own excluded mechanism denominator. |
| A3 | Pass | Every method-dataset has a row; failures have zero unconditional availability/coverage/power, undefined unexecuted FDP/RMSE; dataset-level intervals + successful counts + paired contrasts retained. |
| A4 | Pass | Pooled/condition intensity, count, all-MAR, DEP absence routes; uniform kNN/MinProb; observed-only limma + actual proDA inference. Shared-imputer DEP comparator explicitly is not a reproduction of the full DEP workflow. |
| A5 | Pass | Three replication levels, three fixed manual cutoffs, rescue/no-rescue, real-data ±1 cutoffs, 100 feature bootstraps per rescue policy, 20 rescue seeds, donor/run sensitivity. Unsupported + ambiguous profiles remain visible. |
| A6 | Pass | Both full public source tables have SHA256-pinned snapshots + extraction recipes. PXD020717 flags + donor aggregation explicit. iPRG known spikes survive contaminant filtering; six technical-mixture contrasts retain all six positives. Complete-row masking is labelled semi-synthetic and selected. |
| A7 | Pass | 720 paired downstream outcomes + iPRG84 pipeline/contrast/normalization outcomes +504 spike estimates; `check-results.R` passes denominators/failures. proDA converges in90/90 simulations and12/12 iPRG fits with cap100;89 simulations stop within default20. |
| A8 | Pass | Original-observation assertions pass in imputation and rescue checks. Fresh43-artifact replay from empty destinations is byte-identical; independent repeated HTML render is byte-identical. |
| A9 | Pass | Standalone Rmd, self-contained HTML, editable Word,20 references, three main research figures + four diagnostic profiles. Visual QA covers all three main figures, article opening, and external-calibration table; declarations identify outstanding author decisions. |
| A10 | Pass | `Rscript tools/check-package.R --bioc` passes package build/check and BiocCheck. `reproduce.py` + `check-results.R` are the tracked numerical/regeneration gates; `AGENTS.md` points to both. |

## Adjudicated fixes

- Shared A/B imputation RNG → distinct deterministic condition/imputer streams; paired method seeds preserved.
- Existence-only result reuse → every analysis executes afresh; only checksum-validated immutable source snapshots can be reused.
- Orphan artifacts → empty generated destinations + exact declared inventory; real-file orphan/missing negative controls + valid positive control pass.
- Manifest IDs → POSIX separators on every platform; PureWindowsPath control passes without claiming a Windows execution.
- Diagnostic states of dropped features → separate actionable coverage/correctness in every simulation row and main figure.
- Paired availability difference → Student-t difference interval; Wilson interval reserved for binary availability.
- Real rescue-seed interpretation → report null median effect because every seeded block was subsequently dropped; infer no observed normalization effect.
- proDA silently returns nonconverged fits → explicit numerical validity fields; primary availability requires convergence. Fixed random initialization isolates solver-cap probes and prevents dependence on preceding methods' RNG use.
- Generated diagnostic whitespace → idempotent `tidy-generated.py`; immutable upstream metadata remains untouched.

Adjudicated =10/10; pass=10; unresolved in-scope findings=0. `VERIFICATION.md` records gate outcomes. These rulings apply to the bounded claims, not universal scientific guarantees.

## Deferred claims + editorial dependencies

- General superiority over complete DEP, DAPAR, msImpute, or newer model-based imputation workflows = unsupported and unclaimed; additional comparisons required if that claim is pursued.
- Mechanism identification, calibrated cutoff probability, universal FDR control, and disease discovery = unsupported and unclaimed.
- Article author list, affiliations, funding, competing interests, target journal format, repository release DOI, and actual submission/publication require final author/editorial decisions.
- iPRG = one relatively complete controlled experiment with technical replication; no generalisation to every acquisition regime.
- BiocCheck advisory notes = function length, line length, generated-document indentation, and admin-only mailing-list subscription verification;0 errors/0 warnings. Existing style/account checks do not establish journal readiness or acceptance.
