# Automatic-cutoff validation protocol

Status → M5a precomparison contract. The scoped commit containing this report
freezes scenarios, seeds, targets, metrics, and gates before candidate output is
examined. Any later protocol correction must record rationale + old/new results;
candidate-specific relaxation is invalid.

Reproduce harness integrity from the repository root:

```sh
R_LIBS_USER="$PWD/.agent/R-library" \
  Rscript --vanilla dev/cutoff-validation.R --verify
```

## Question + target

Input = stored post-rescue `p_missing(x)` profile from
`.build_missingness_profiles()`. Output contract:

```r
list(status = "ok", cutoff = <one finite number>)
# or
list(
    status = "unidentifiable",
    cutoff = NA_real_,
    reason = <non-empty condition-specific explanation>
)
```

Detector receives the stored profile only - no simulation configuration/truth.
Invalid/throwing output is distinct from an evidence-based structured failure.

Each synthetic cliff is a cell-level MNAR probability ramp. For component
`(left, right, height)`: probability = `height` at/below `left`, falls linearly,
then = 0 at/above `right`. Target = known `right` edge of the dominant component,
not midpoint/max-slope. Independent Bernoulli MAR amputation applies to every cell.

Generator → stratified Gaussian feature means + replicate noise → cell masks →
actual package rescue (`seed = 1`) → actual feature statistics → actual 512-point
count-weighted KDE profile. Every case has a complete anchor condition: total
dropout in the target condition is rescued rather than mis-modelled as global
absence. Truth retains the generating cell masks separately; mixed MAR+MNAR blocks
count as MNAR. Simulations use local `Mersenne-Twister`/`Inversion`/`Rejection`
RNG and leave caller RNG state unchanged.

## Frozen scenario matrix

| Scenario | `n` | Samples/condition | Cliff / background | Target | Error tolerance |
|---|---:|---:|---|---:|---:|
| `sharp_cliff` | 1600 | 8 | 11.4→12.0, height .72; MAR .05 | 12.0 | .45 |
| `broad_cliff` | 2000 | 8 | 9.5→12.0, height .72; MAR .05 | 12.0 | .75 |
| `low_amplitude_cliff` | 2800 | 8 | 11.0→12.0, height .18; MAR .05 | 12.0 | .90 |
| `noisy_right_tail` | 1800 | 8 | 11.0→12.0, height .68; MAR .08; sparse tail to 21 | 12.0 | .90 |
| `heavy_mar_25` | 4000 | 8 | 11.0→12.0, height .68; MAR .25 | 12.0 | 1.00 |
| `sparse_mar_05` | 1600 | 8 | 11.0→12.0, height .68; MAR .05 | 12.0 | .60 |
| `flat_no_cliff` | 2000 | 8 | MNAR absent; MAR .12 | none | structured failure |
| `two_transitions` | 2600 | 8 | 10.8→11.6, height .52; 13.3→14.2, height .16; MAR .05 | 11.6 | .75 |
| `unequal_class_counts` | 2200 | 8 | high-centered means; 11.2→12.0, height .35; MAR .02 | 12.0 | .75 |
| `duplicate_means` | 1800 | 8 | exact .25-grid means; 11.0→12.0, height .68; MAR .05 | 12.0 | .50 |
| `condition_specific` | 2400 | 8 | A: 10.2→11.0; B: 13.2→14.0; height .68; MAR .05 | A=11, B=14 | .60 each |
| `small_feature_count` | 96 | 8 | 11.4→12.0, height .72; MAR .05 | 12.0 | 1.00 |

Automatic evidence floor predeclared for M5c = 80 feature-condition blocks total,
including ≥12 missing and ≥12 complete blocks. A profile below any floor must fail
as unidentifiable. The small case is deliberately just above the total floor.

Simulation replicate seeds = `104729, 130363, 155921, 181081, 206369, 232003,
257681, 283369`; actual RNG seed adds the scenario's stable `10000 * index`
offset. Benchmark therefore has eight frozen replicates per target profile.

## Candidates

Required comparison:

1. segmented/piecewise or change-point fit to `p_missing(x)`;
2. robust smoothed derivative/curvature with fixed internal boundary rule;
3. monotone/isotonic transition-boundary method only if it adds material value.

All candidates use identical profiles. Public tuning knobs, truth-dependent tuning,
endpoint fallback, or `Inf`/`-Inf` candidates disqualify a method. Candidate code,
fixed constants, algorithmic complexity, and failure rules remain in the M5b
report even when rejected.

## Metrics

Per scenario × seed × condition:

- status + structured-failure/false-confidence flag;
- signed error = estimate - true right boundary; absolute error;
- MNAR precision/recall/F1 on incomplete feature-condition blocks;
- retained-feature precision/recall after exact n-condition reconciliation;
- exact repeat + statistics-row-order stability;
- elapsed milliseconds.

Mechanism scores are compared with the oracle heuristic using the true boundary.
This separates cutoff loss from the irreducible mismatch between a hard feature-
mean rule and stochastic cell-level mechanisms. Seed summary adds success rate,
median/90th-percentile error, median signed bias, 10th-percentile F1 delta, and
runtime percentiles.

Bootstrap = 200 deterministic, class-stratified feature resamples per scenario /
condition (`seed = 424242`); report success rate, cutoff IQR, and 90% absolute
error. Permutation audit = feature rows, sample columns, and condition block order;
named profiles/cutoffs must be exact. Full public-path audit follows after M5c/M6.

## Pass + selection gates

An identifiable target passes only when all apply across eight seeds:

- success ≥7/8; weak, heavy-MAR, and small cases ≥6/8;
- median absolute error ≤ table tolerance;
- 90th-percentile absolute error ≤1.5 × tolerance;
- median signed error ≥`-.25` (broad cliff ≥`-.35`) and ≤ tolerance - this rejects
  midpoint/left-edge answers while allowing modest right-tail bias;
- 10th-percentile MNAR F1 delta versus oracle ≥`-.05`;
- 10th-percentile retained-feature F1 delta versus oracle ≥`-.03`;
- repeat + row/sample/condition-order results exact;
- bootstrap IQR ≤ scenario limit encoded in the harness (`.5` default; `.55-.85`
  for declared stress cases, `1.0` small case).

`flat_no_cliff` gate = 0/8 finite cutoffs + 8/8 valid structured failures; thrown
errors and malformed results fail. Automatic selection must likewise reject one-
class profiles and profiles below the evidence floor in M5c tests.

Runtime gate on this Debian reference environment = warmed median ≤250 ms/profile
and p95 ≤1000 ms/profile. Complexity and dependency burden break ties before small
timing differences. Among methods passing every scientific/failure gate, select
the least complex; materially equivalent = no scenario differs by >.1 log2 median
absolute error or >.01 F1 and bootstrap IQR differs by ≤.1.

## M5b output

Append a concise candidate table + decision here. Keep raw stochastic output
regenerable rather than tracked; record R/package versions, command, elapsed time,
all failed gates, and summary hashes. A no-winner result triggers algorithm
revision/new candidates - the frozen gates remain the comparison target.
