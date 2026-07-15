# Scientific-regression validation protocol

Status → protocol v1 frozen in M7a; full M7b assessment appended below. This
report and `dev/scientific-validation.R` freeze simulation truth, metrics, gates,
seeds, and the routine/long-run split before any ImputeFinder result is inspected.
Protocol version = `1`; later scientific corrections must record rationale and
old/new outcomes. Result-specific threshold relaxation is invalid.

Verify the generator/oracle contract from the repository root:

```sh
R_LIBS_USER="$PWD/.agent/R-library" \
  Rscript --vanilla dev/scientific-validation.R --verify
```

M7a verification hash for the serialized protocol manifest =
`4011e381bba2d0d747e91d277a45de5e`. The enclosing scoped git commit freezes the
generator and scoring implementation as well. M7a verification invokes no
classifier and establishes no performance claim.

## Question + independence

Question → does condition-specific ImputeFinder classification recover causal
missingness and retention behavior when data contain both uniform MAR and
intensity-dependent MNAR, including biological on/off features? Separately, how
closely does automatic selection reproduce the known cutoff and the manual
true-cutoff result?

This simulation is independent of the M5 detector benchmark:

- M5 validates cutoff placement from purpose-built profile shapes and uses a
  complete anchor condition.
- M7 generates complete sample-level data for two experimental conditions,
  amputates individual cells, includes structural on/off blocks, and scores the
  public pipeline against causal masks and cross-condition retention truth.
- production helpers generate neither data nor truth; the public classifier will
  see only the amputated matrix, aligned groups, cutoff arguments, and rescue seed.

## Frozen data-generating model

Every scenario has conditions `A` and `B`, equal replicate counts within that
scenario, unique feature/sample names, and unnormalised log2-scale intensities.
For feature `i`:

1. Stratified normal quantiles generate a relative feature mean with center
   `cutoff + 1` and SD `2.2`; seed-specific shuffling removes intensity ordering.
2. A condition contrast drawn from `N(0, 0.35^2)` is split equally/oppositely
   across A/B.
3. Replicate error is independent `N(0, 0.18^2)`.
4. Four percent of features in each long scenario are forced off in A and four
   percent off in B. Off/on means are `cutoff - 2.5` / `cutoff + 2.5`; every cell
   in the off block is causally MNAR. Sets are disjoint.
5. MAR is an independent Bernoulli draw for every cell with scenario rate `.05`
   or `.25`.
6. For non-structural cells, causal MNAR probability is

```text
0.8 * clamp((cutoff - complete intensity) / 1.2, 0, 1)
```

7. A cell is missing when either causal mask is true. MAR/MNAR overlap remains
   explicit; observed complete intensities and both masks are retained only in
   the development truth object.

Simulation RNG = local `Mersenne-Twister` / `Inversion` / `Rejection`; caller RNG
kind/state must be byte-identical before/after generation. Fixed simulation seeds:

```text
104677, 130387, 155939, 181123, 206407, 232019, 257707, 283397
```

Classifier rescue seeds = `1, 7, 29`. Default-seed results define performance;
the alternatives audit that seed-cell choice changes neither statistics, states,
retention, nor cutoffs.

## Frozen scenario matrix

All group-size scenarios use cutoff `12` and 2,000 features. Each has 80 features
off per condition (160 total).

| Scenario | Samples/condition | Uniform MAR | Automatic gate |
|---|---:|---:|---|
| `group_n04_mar05` | 4 | 5% | required |
| `group_n08_mar05` | 8 | 5% | required |
| `group_n20_mar05` | 20 | 5% | required |
| `group_n04_mar25` | 4 | 25% | stress |
| `group_n08_mar25` | 8 | 25% | stress |
| `group_n20_mar25` | 20 | 25% | evidence sensitivity |

The cutoff sweep uses 1,600 features, 8 samples/condition, 5% MAR, and cutoffs
`8, 9, ..., 14`. For a fixed seed, all latent values and masks are identical
after subtracting the cutoff shift. Thus the sweep tests location equivariance,
not seven accidentally different random datasets.

The 20-sample/25%-MAR case is predeclared as an evidence-sensitivity scenario:
even without MNAR, a feature's probability of being complete is
`.75^20 = .00317`, or roughly six of 2,000 features. The production detector
requires at least 12 missing and 12 complete feature blocks. A cutoff is therefore
often not identifiable from its missing-vs-complete profile. This case measures
safe coverage/failure rather than imposing a scientifically impossible success
rate.

## Causal oracle

Truth is assigned from realized masks before rescue:

```text
globally absent across A+B       -> all_missing
no missing cell in condition    -> complete
at least one causal MNAR cell    -> MNAR
MAR-only + strict majority seen -> MAR
MAR-only + no strict majority   -> insufficient
```

MNAR dominates an overlapping MAR mask. This asks whether a block needing an
MNAR branch is found; it does not relabel stochastic causal truth to make the
intensity heuristic look exact. A condition block missing entirely by chance
from MAR is therefore `insufficient`, even though production rescue will insert
one low anchor and may make an error. Such cases remain in the score.

Oracle reconciliation applies the production specification independently:
`all_missing` → any `insufficient` condition → `MNAR_all_conditions` → retain.
Structural on/off retention is scored only when the oracle says the feature is
usable (the on condition is complete/MAR and the feature is not globally absent).

## Metrics

Each public result is aligned to truth by `feature + condition`; positional
agreement is insufficient. Scores:

- accuracy and macro-F1 over `complete`, `MNAR`, `MAR`, `insufficient`;
- MNAR and MAR precision/recall/F1 from causal states;
- insufficient F1;
- retention precision/recall/F1 over every original feature;
- eligible structural on/off retention recall;
- structural off-block rescue recall;
- automatic signed/absolute cutoff error from the known boundary;
- automatic-minus-manual MNAR-F1 and retention-F1 deltas.

For a binary class, precision = `TP/(TP+FP)`, recall = `TP/(TP+FN)`, and F1 is
their harmonic mean. Undefined class scores remain `NA`; they are not converted
to perfect values. Macro-F1 averages finite state F1 values. Long-run quantiles
are calculated separately per scenario and condition across eight simulation
seeds, so one condition cannot hide failure in the other.

## Predeclared gates

Manual runs use the true named cutoff in both conditions. Every scenario must
meet its MAR-rate row:

| Uniform MAR | q10 MNAR F1 | q10 MAR F1 | q10 macro-F1 | q10 retention F1 |
|---:|---:|---:|---:|---:|
| 5% | ≥.72 | ≥.72 | ≥.76 | ≥.94 |
| 25% | ≥.60 | ≥.60 | ≥.66 | ≥.90 |

For every individual run, eligible on/off retention recall and structural
off-block rescue recall must both equal `1`. The input object and caller RNG must
remain exact. Every returned non-seed observed value must equal its complete
input value; every difference must have one matching seed-log row and the
condition's pre-seed minimum.

Automatic runs are paired with the manual result from the same generated matrix:

- required cases → success ≥7/8 for each condition;
- 25%-MAR stress cases at 4/8 samples → success ≥6/8 for each condition;
- median absolute cutoff error ≤`.75`, q90 ≤`1.25`;
- median signed error ∈`[-.35, .75]` to reject a left/mid-transition answer;
- q10 MNAR-F1 delta from manual ≥`-.05`;
- q10 retention-F1 delta from manual ≥`-.03`;
- every finite cutoff lies strictly inside observed feature-mean support.

For `group_n20_mar25`, a profile below either production evidence floor must
raise the structured, condition-specific unidentifiable-cutoff error. An eligible
profile is scored by the ordinary automatic gates. Any finite estimate, even from
an otherwise ineligible run, must be within `1.25` of truth. Coverage and exact
class counts are reported rather than hidden.

For each seed, the 8-14 sweep requires exact missing masks, manual states, and
retention after translating intensities/cutoffs. Automatic decisions must have
the same status, automatic state/retention, and cutoff error within `1e-10`.

Permutation audit scenarios = `group_n04_mar25`, `group_n08_mar05`, `sweep_c08`,
and `sweep_c14`. Even/odd feature permutation plus B-first/reversed sample blocks
must preserve named cutoffs, diagnostics, profiles, classifications, retention,
groups, and seed assignments exactly after canonical name alignment. Exactness,
not a numerical tolerance, is the gate.

## Routine checks vs long benchmark

Routine package tests will encode only these version-1 cases after M7b judges
them:

| Case | Features | Samples/condition | MAR | Public path |
|---|---:|---:|---:|---|
| `manual_on_off` | 320 | 4 | 25% | true cutoff `12` |
| `automatic_cliff` | 800 | 8 | 5% | automatic + one permutation |

Both use simulation seed `104677` and rescue seed `1`. The manual case must meet
the 25%-MAR per-run thresholds and exact on/off gates above. The automatic case
must meet the 5%-MAR thresholds, cutoff error ≤`1`, MNAR-F1 delta ≥`-.05`,
retention-F1 delta ≥`-.03`, and exact named permutation output. These cases are
the smallest predeclared subset that exercises causal scoring, dynamic majority,
condition rescue, reconciliation, automatic placement, and order stability.

The full 13-scenario × 8-seed manual/automatic benchmark, alternative rescue-seed
audit, and four-scenario permutation audit remain in `dev/`; they do not burden
routine package checks. M7b adds the benchmark runner and the concise summaries,
failed-gate ledger, R/package context, hashes, and elapsed time below. Raw
stochastic rows remain regenerable and untracked.

## M7b full assessment

Run from a current project-local installation:

```sh
R_LIBS_USER="$PWD/.agent/R-library" R CMD INSTALL --preclean .
R_LIBS_USER="$PWD/.agent/R-library" \
  Rscript --vanilla dev/scientific-validation.R --benchmark
```

Outcome → all `270/270` frozen gates pass. No threshold, generator, oracle, or
production-method change followed result inspection. The only first-pass
failure was a harness-only character/factor comparison in the permutation
provenance audit; normalizing the simulated factor label to character changed
that audit from `656/664` to `664/664` and changed no classifier output or
scientific metric.

### Manual true-cutoff results

Each cell is condition `A / B`; values are q10 across eight simulation seeds.
The translated sweep has exactly the same outcome at every cutoff 8-14.

| Scenario | MNAR F1 | MAR F1 | Macro-F1 | Retention F1 |
|---|---:|---:|---:|---:|
| `group_n04_mar05` | .9727 / .9765 | .9473 / .9486 | .9474 / .9587 | .9921 / .9921 |
| `group_n08_mar05` | .9771 / .9749 | .9683 / .9651 | .9584 / .9140 | .9904 / .9904 |
| `group_n20_mar05` | .9659 / .9677 | .9722 / .9737 | .9794 / .9805 | .9863 / .9863 |
| `group_n04_mar25` | .9333 / .9316 | .9639 / .9606 | .9629 / .9620 | .9564 / .9564 |
| `group_n08_mar25` | .9553 / .9540 | .9776 / .9778 | .9754 / .9706 | .9756 / .9756 |
| `group_n20_mar25` | .9638 / .9620 | .9812 / .9795 | .9740 / .9725 | .9854 / .9854 |
| `sweep_c08`-`sweep_c14` | .9695 / .9744 | .9586 / .9651 | .9760 / .9650 | .9890 / .9890 |

Worst q10 margins remain well inside the frozen gates: at 5% MAR, MNAR/MAR/
macro/retention F1 minima = `.9659/.9473/.9140/.9863` against
`.72/.72/.76/.94`; at 25% MAR, minima = `.9316/.9606/.9620/.9564` against
`.60/.60/.66/.90`. Every one of 208 condition-level manual summaries has
eligible on/off retention and structural off-block rescue recall = `1`.

### Automatic results

Each cell is condition `A / B`. Success is out of eight seeds; errors and deltas
summarize successful runs. `MedAE`/`q90AE` = absolute cutoff error; `MedSE` =
signed error; `ΔMNAR`/`Δretain` = q10 F1 delta from the paired manual result.

| Scenario | Success | MedAE | q90AE | MedSE | ΔMNAR | Δretain |
|---|---:|---:|---:|---:|---:|---:|
| `group_n04_mar05` | 8/8 / 8/8 | .06792 / .04734 | .11705 / .09660 | -.06792 / -.04734 | -.00399 / -.00718 | -.00189 / -.00189 |
| `group_n08_mar05` | 8/8 / 8/8 | .06491 / .05629 | .12277 / .15006 | .06491 / .01847 | -.00379 / -.00320 | -.00147 / -.00147 |
| `group_n20_mar05` | 8/8 / 8/8 | .11741 / .14755 | .22364 / .25216 | .11741 / .14755 | -.00016 / -.00633 | -.00148 / -.00148 |
| `group_n04_mar25` | 8/8 / 8/8 | .13476 / .04945 | .29373 / .22908 | -.13456 / -.04062 | -.01245 / -.00349 | -.00248 / -.00248 |
| `group_n08_mar25` | 8/8 / 8/8 | .11676 / .08150 | .17972 / .33788 | -.04819 / .04182 | -.00954 / -.04793 | -.01076 / -.01076 |
| `group_n20_mar25` | 0/8 / 0/8 | NA / NA | NA / NA | NA / NA | NA / NA | NA / NA |
| `sweep_c08`-`sweep_c14` | 8/8 / 8/8 | .04437 / .04889 | .07606 / .08457 | .01830 / .04889 | -.00548 / -.00104 | -.00087 / -.00087 |

Every one of 192 finite cutoffs is strictly inside observed feature-mean
support. Across eligible profiles, maximum median absolute error = `.14755`,
maximum q90 absolute error = `.33788`, signed medians span
`[-.13456, .14755]`, worst q10 MNAR-F1 delta = `-.04793`, and worst q10
retention-F1 delta = `-.01076`. All clear/stress coverage is `8/8` per
condition, exceeding required `7/8` and stress `6/8` gates.

The predeclared evidence-sensitivity case behaves structurally: all 16
condition profiles have fewer than 12 complete blocks and return the exact
condition-specific unidentifiable-cutoff error. Complete-block counts by seed,
shown as `A/B`, are `4/1, 3/2, 5/3, 3/6, 7/2, 5/2, 7/6, 2/5`; no cutoff is
fabricated.

### Exact audits + provenance

- 664/664 public calls preserve input, caller RNG kind/state, observed values,
  and seed provenance exactly.
- Rescue-seed invariance: 208/208 manual and 416/416 condition-specific
  automatic comparisons exact for seeds 7/29 versus seed 1.
- Translation sweep: 48/48 masks, 48/48 manual state/retention comparisons,
  and 96/96 automatic statuses, state/retention comparisons, and cutoff-shift
  errors within the frozen `1e-10` exactness tolerance.
- Name-aligned permutation: 4/4 manual and 8/8 condition-specific automatic
  results exact, including cutoffs, diagnostics, profiles, classifications,
  retention, groups, data, and seed assignments.
- Structured automatic failures report the target condition and carry matching
  diagnostic/profile evidence in 16/16 cases.

Execution context → R `4.6.1`, `x86_64-pc-linux-gnu`; installed imputefinder
`0.0.0.9000` from `.agent/R-library`. Protocol MD5 =
`4011e381bba2d0d747e91d277a45de5e`; deterministic result MD5 =
`72142402979a25eccead24b5af4bba11`; benchmark-harness MD5 =
`e2e755fa36f69055e8acb23ce638816e`. Final measured elapsed time = `24.390 s`.
Raw rows remain intentionally untracked and regenerate to the same result hash.

Conclusion → protocol-v1 evidence supports the manual classifier, automatic
boundary detector, safe evidence-floor failure, condition-specific on/off
retention, order invariance, and seed provenance across the frozen scope. No
runtime defect was identified.

## M7c routine regression closure

Ordinary package tests now regenerate only the two predeclared seed-104677 cases
from a compact test-only protocol mirror. The mirror's matrices, aligned groups,
oracle states, retention truth, and scores were checked exact against this
harness; `dev/` remains excluded from the source package and owns every long run.

| Case/path | MNAR F1 | MAR F1 | Macro-F1 | Retention F1 | Cutoff A/B |
|---|---:|---:|---:|---:|---:|
| `manual_on_off`, manual | .94382 | .96736 | .96805 | .95094 | 12 / 12 |
| `automatic_cliff`, manual | .98155 | .97222 | .98459 | .99065 | 12 / 12 |
| `automatic_cliff`, automatic | .97770 | .96703 | .98158 | .98983 | 11.97274 / 11.88546 |

Both cases retain every eligible on/off feature and rescue every structural off
block. Automatic-minus-manual MNAR/retention F1 = `-.00385/-.00082`; the full
canonical automatic result is exact after the frozen even/odd-feature,
B-first/reversed-sample permutation. All 16 routine expectations, 505 package
expectations, protocol verification, build, and source-tarball tests passed. At
M7c closure, the only package-check warning was the licence placeholder assigned
to M8. M7 is closed.
