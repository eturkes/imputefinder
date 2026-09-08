# ImputeFinder publication protocol

## Scope + estimands

- Software article = auditable condition-aware operational routing, not identified missingness mechanisms or validated disease effects.
- Primary simulation estimands = automatic availability; original missing-cell routing coverage + correctness against recorded simulated causes; mechanism-specific recall; retention + on/off retention. No threshold tuning against outcomes.
- Independent simulation units = datasets. Interval = two-sided 95% Student-t Monte Carlo interval across dataset replicates; bounded proportions clipped to [0,1]. Availability = Wilson binomial interval. Not patient/population uncertainty.
- Automatic failures = denominator-preserving rows; score availability separately; failed classifications receive zero unconditional routing coverage/correctness. Selective accuracy excludes unassigned cells, reports its denominator.
- Correlated features + cells = never treated as independent replicates for intervals. Numerical seeds = fixed; estimands may be undefined when a cause is absent.
- Real data = full protein-group intensity table from BA46 soluble PXD020717 archive; exclude reverse, contaminant, site-only rows; retain every remaining row/sample. No outcome-driven site selection. Use observed complete-row masking only for reconstruction; explicitly acknowledge selection bias + unrealistic masking.
- Downstream experiment = paired identical latent datasets + fixed imputation algorithms. Evaluate reconstruction, effect error, observed-only and imputed differential analysis, realized FDP + power. BH 0.05 labels do not guarantee FDR after single imputation.

## Fixed experiments

1. Eight missingness scenarios × replication (3, 6, 10) × 30 dataset seeds; 600 features × two conditions; initial20% nonzero effects; overwrite5% as forced on/off blocks, with final effect truth recorded. Scenarios = MCAR, observed-covariate MAR, smooth low-tail MNAR, hard censoring, dual transition, high-tail MNAR, heteroscedastic/batch, mixed MAR/MNAR.
2. Routes = automatic minimum/no rescue; fixed log2 cutoff 12 with ±1 perturbation; pooled 25th-percentile intensity; condition 25th percentile; observed-majority; all MAR; DEP-style fully absent condition rule. Comparator labels are operational heuristics, not reproduced published full pipelines.
3. Downstream = 30 seeds × MCAR/mixed/heteroscedastic × six samples/condition; shared imputers = MinProb + kNN; observed-only limma + proDA; uniform MinProb/kNN; DEP-style and ImputeFinder mixed routing. Comparisons count pipeline failures + retained/tested fractions. Selection-aware power denominator = every true changed feature; FDP denominator = discoveries. Methods receive distinct deterministic condition/imputer random streams. proDA = inference comparator, reconstruction not applicable.
4. Real data = automatic success/ambiguity; cutoff bootstrap (100 feature resamples, failure-inclusive); rescue/no-rescue and ±1 cutoff sensitivity; complete-row held-out MCAR and low-tail masking (30 seeds). Held-out truth only assesses measured cells.
5. External controlled calibration = iPRG2015 complete raw-Intensity protein table, four mixtures from one yeast background × three technical runs. Six declared spikes with known nonzero effects in every pairwise comparison; six contrasts, seven pipelines, observed-median normalization + unnormalized sensitivity. Count all six spikes even if dropped/untested; common original-feature BH universe. Report background discoveries + realized FDP; six correlated contrasts are not independent experiments. Spike amounts establish effects, not cell missingness mechanisms. Preserve known spikes before generic contaminant exclusion. PXD028038 reprocessing enables match-between-runs; upstream identification is not rerun.
6. proDA numerical validity = explicit fixed initialization seed; uniform cap100 + unchanged tolerance0.001; record convergence/iterations/error + whether stopping occurs within default20. Nonconvergence means primary unavailability. Solver-only seed1 probes compare caps20/100 with identical initialization, without effect/discovery outcomes. Cap100 accommodates a default-cap failure that converges at iteration22. No solver-unconverged superiority claims.

## Fixed acceptance/review set

| ID | Required evidence |
|---|---|
| A1 | API rules + limits match source; no causal/novelty/superiority overclaim |
| A2 | Independent generator + fixed parameters; no implementation-derived truth |
| A3 | Every run/failure represented; dataset-level uncertainty + denominator |
| A4 | Pooled, count, uniform, and existing condition-aware routing comparators |
| A5 | Rescue, sample count, cutoff, unsuitable-profile sensitivity |
| A6 | Real-data source integrity + full extraction + held-out selection caveat |
| A7 | Downstream estimates + power/FDP + retention; identical paired inputs |
| A8 | Original observations remain unchanged; deterministic rerun output |
| A9 | Standalone rendered article + references + figures/tables + supplement |
| A10 | Package tests/build/check + committed regeneration/check entry point |

Adjudicate these ten rows once; fixes receive one finding-specific recheck. Outside findings = deferred, not silent scope expansion. Unmet rows remain explicit.
