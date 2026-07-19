# M13c association addendum v3

Status: pre-result contract frozen; all three candidates remain unrun. Synthetic
replicates 1-32 are the only association evidence authorized next. Association
computations on replicates 33-64, HarmonizR result bytes, and all confirmation
artifacts remain sealed.

## Scope and bindings

M12 fixed the candidates and thresholds before results, but Track A v2 left the
scalar hypothesis, acquisition handling, support cells, permutation groups,
runtime schema, and stage denominators under-specified. This standalone
`m13_a_association_protocol_v3` supersedes those details without rewriting the
historical M12 objects. It binds:

- M12 contract/gate registry:
  `45d1cda936b21a42d6735d900917734398eecb3ff1c136cab486cf90ee2b21e5` /
  `70d101b890de687be8973c9ff263caa454b61505f99470e5ec6061cb344db92f`;
- M12 Track A/protocol bundle:
  `ca0cf8dbbf082446d9116ce280f0acf2c6517d35ec6137edb7b415585ce92683` /
  `3bc4ef531d6c84475114befc5c153e50680459f1d8125fdbabac5fd840a1de65`;
- generator protocol/audit:
  `cdea1bf874152e63fba08c49e390da1dd18690105e43a2dc7a03fc6866d0d080` /
  `4d216966c8cc54b1b23685d8801c7c65b2afcc5879dd15430d57a88f765e8a00`;
- `design_estimability_v1`, `canonical_treatment_contrasts_v1`,
  `svd_relative_v1`, and the other frozen design-core identifiers.

The executable tables, encoder, seed manifest, schema validator, and gate
evaluator live in `dev/m13-association-contract.R`. Its marked normative source
region and this full Markdown file are themselves SHA-256 inputs to the protocol
hash. The verifier prints component/protocol/aggregate identities; this file
does not embed a self-referential v3 hash.

## Response and acquisition strata

Let `D_fi = 1` iff the original pre-rescue value for feature `f`, sample `i` is
finite. Freeze the globally observable feature set over the complete supplied
input before acquisition splitting:

```
F+  = {f: sum_i D_fi > 0}
G   = |F+|
k_i = sum_{f in F+} D_fi
y_i = k_i / G
```

`G = 0` yields `association_no_observable_features`. Intensities never enter a
candidate. Store probability units; percentage points are presentation only.

When acquisition is declared, each canonical acquisition label defines a
separate analysis. Subset samples first, remove acquisition and every
acquisition-containing interaction, then rebuild encoding, SVD, aliases, units,
support, and hypotheses. A globally built model must not be subset. Declared
stratum IDs are `s_` plus the full SHA-256 of the length-prefixed acquisition
label; without acquisition, use stratum `all` and display label `undeclared`.
Holm families never cross strata, and acquisition itself receives no
hypothesis.

## Coefficient hypotheses and identities

Each eligible canonical encoded column is one 1-df raw coefficient-axis
hypothesis:

- condition and nuisance main effects - categorical nonreference level versus
  radix-first reference, or plus one numeric encoded unit;
- explicit interaction product columns - only when every component role is
  condition or nuisance; each product column is separate, never omnibus.

Intercept, block coefficients, block-containing interactions, acquisition, and
acquisition-containing interactions remain adjustment/descriptive-only and are
absent from the hypothesis table. Every eligible role coefficient remains in
that table even if nonestimable or unsupported. Its hypothesis ID is `a_` plus
the full SHA-256 of the version tag and length-prefixed stratum, coefficient,
label, and term ID. The outcomes list uses exactly those IDs, in hypothesis
order. Available and existing `imputefinder_unavailable` records set `quantity`
to the same list key, so every abstention joins deterministically.
Each hypothesis also records a named `numeric`/`treatment` encoding vector whose
names and order exactly equal its canonical components.

Raw coefficient axes must pass the existing row-space estimability rule. Only
finite available raw p-values enter Holm. Nonestimable/unsupported outcomes
remain explicit; they never enter a successful-fit denominator.

These are conditional model-encoded associations - not causal effects,
missingness mechanisms, structural-absence claims, or automatic corrections.

## Shared rank-coordinate algebra

Within a rebuilt stratum, use the existing rank rule and coordinates:

```
X = U_r D_r V_r'
Z = X V_r = U_r D_r
K = (Z'Z)^-1
a_hat = K Z' y
beta_hat = V_r a_hat
H = U_r U_r'
M = I - H                         [independent]
M = Q_perp Q_perp'                [blocked]
e = M y
d0 = n - r
```

For blocked CR2, obtain `[Q_r,Q_perp]` from the complete LAPACK Householder QR
of `U_r`, form `M = Q_perp Q_perp'`, and symmetrize it. This is mathematically
`I-H`; the orthogonal-complement construction is frozen because direct
subtraction can create tolerance-scale negative cluster-block eigenvalues.

For raw coefficient axis `l_j`, set `c_j = V_r' l_j` and
`theta_j = c_j'a_hat`. This retains every estimable direction without inverting
a singular raw `X'X`. Require `d0 >= 3`.

For any symmetric numeric object `A`, the general backward-error tolerance is

```
tau(A) = max(dim(A)) * .Machine$double.eps *
         max(1, max(abs(A)))
```

Covariances are symmetrized. An eigenvalue below `-tau(V)` fails; values in
`[-tau(V), tau(V)]` become zero; scalar contrast variance must exceed `tau(V)`.
No ridge, substitute variance estimator, or silent singular-hypothesis
pseudoinverse is allowed.

### Independent HC3

Clamp leverage into `[0,1]` only when its excursion is within `tau(H)` and
require `1-h_i > tau(H)`:

```
q_i = e_i / (1-h_i)
V_HC3 = K Z' diag(q_i^2) Z K
v_j = c_j' V_HC3 c_j
t_j = theta_j / sqrt(v_j)
```

Use a two-sided `t_d0` p-value and 95% interval.

### Full-design CR2

With a declared block, retain the complete design including block fixed
effects; clusters are blocks. For cluster selector `C_g`:

```
B_g = C_g M C_g'
A_g = B_g^(+1/2)
s_g = Z_g' A_g e_g
V_CR2 = K [sum_g s_g s_g'] K
```

Compute the symmetric Moore-Penrose inverse square root with

```
tau_g = tau(B_g) = max(dim(B_g)) * .Machine$double.eps *
        max(1, max(abs(B_g)))
```

Materially negative eigenvalues fail; values within tolerance define the
mathematical null space. Scalar Satterthwaite reference df is

```
p_g = M C_g' A_g Z_g K c_j
nu_j = [sum_g p_g'p_g]^2 / sum_g sum_h (p_g'p_h)^2
```

Use `t_nu_j`; `nu_j < 5` yields `association_low_reference_df`.

## Exact common support

Biological units are declared blocks when present and samples otherwise.
Technical duplicates never count independently. A unit-position is:

- independent - canonical sample ID;
- blocked - block ID plus canonical values of every condition and nuisance
  column.

Rows with the same complete blocked key are technical siblings and collapse
once. Declared numeric values must agree after integer-to-double normalization;
otherwise use `association_incompatible_unit_positions`. Each unique position
gets equal median weight, regardless of sibling count.

- Categorical main coefficient - reference versus encoded target; at least 4
  independent units per side, or 6 blocks containing both sides and 6 blocks
  per side.
- Numeric main coefficient - median over unique unit-positions; exact ties are
  on neither side; at least 4 independent units per side, or 6 blocks with at
  least one position on each side. Retain range and whether encoded 0-to-1
  evaluation extrapolates.
- Interaction coefficient - form every categorical reference/target and every
  numeric below/above-median side, then materialize all `2^k` Cartesian cells.
  Every cell requires 4 independent units or 6 blocks. A block counts once per
  occupied cell, regardless of repeated positions, and must occupy every cell
  to be complete. Stored low/high counts are the minima across component
  marginal sides. The support row carries a canonical nested table with one
  median, range, and 0-to-1 extrapolation flag per numeric component; a scalar
  summary is insufficient when an interaction contains multiple numeric
  variables.
  Numeric support rows must equal that recorded numeric-component subset
  exactly; their extrapolation flags are derived from the stored ranges.

The role-eligible manifest retains all eligible columns. The candidate-agnostic
opportunity manifest is its subset with observable response, estimability,
common support, and residual-df floor. Nonestimable/support failures remain in
the disposition audit but outside scope coverage; method failures and
candidate-out-of-scope outcomes remain opportunity misses.

## Restricted Freedman-Lane

For estimable direction `c_j`, build the deterministic ordered-axis
modified-Gram-Schmidt basis `N_j` for `{a: c_j'a=0}`:

```
Z0 = Z N_j
H0 = Z0 Z0+
y0_hat = H0 y
e0 = (I-H0) y
y_b* = y0_hat + P_b e0
```

Every transform refits the full model and recomputes the complete HC3/CR2 Wald
statistic `W_b* = theta_b*^2/v_b*`. When `P_b e0` is byte-identical to `e0`,
refit the original response bytes rather than numerically recomposing the same
`y`; this preserves the literal tie for every distinct residual-stabilizer map.
A failed transform makes the hypothesis unavailable; no draw is discarded.
Require

```
max(abs(P_b Z0-Z0)) <= sqrt(.Machine$double.eps) * max(1,max(abs(Z0)))
```

`P_b` acts on reduced-model residual sources; the design and its labels remain
fixed at target rows. Independent sources stay inside their constrained-null
row group. Blocked sources stay inside the same block and matched nuisance
signature. The tested condition pivot may exchange reference and target -
requiring its label vector to be invariant would forbid the intended test.

Independent maps permute canonical unit IDs inside constrained-null row groups.
Two rows share a group only when every pair meets the tolerance above; a
nontransitive tolerance relation is incompatible. Radix-sort groups and units,
then enumerate distinct canonical row-index vectors in global lexicographic
order with identity first.
`T = product_g factorial(n_g)`.

Blocked maps require a categorical condition pivot: the radix-first condition
component in the hypothesis. Within each block, reference and target bundles
are matched by the canonical values of every nuisance column; sibling
multiplicity must match and siblings map positionally. Other condition levels
and incomplete blocks stay fixed. Radix-sort blocks, put identity first, then
use ascending binary swap masks; `T = 2^B`. No condition pivot, an untouched
tested direction, or a failed invariant yields
`association_incompatible_permutation`.

Distinct means distinct canonical row-index permutation vectors - even if two
maps produce equal responses/statistics. Deduplicate maps before counting.

- `T < 20` - `association_low_permutation_resolution`;
- `20 <= T <= 100000` - enumerate all maps and use literal
  `sum(W_b* >= W_obs)/T`;
- `T > 100000` - draw 9,999 distinct nonidentity maps and use
  `(1 + sum(W_b* >= W_obs))/10000`.

The branch uses `log(T)`, so every finite map space above 100,000 reaches Monte
Carlo even when `T` exceeds double range. The diagnostic
`allowable_transformations` is exact for enumerated spaces; for Monte Carlo it
stores the nearest finite-double `T`, saturated at `.Machine$double.xmax`.

Finite doubles use literal `>=`; there is no tie jitter. The Monte Carlo seed
function accepts the complete candidate/acquisition/hypothesis manifest in one
call, sorts it canonically, and resolves collisions globally before returning
any seed. Its SHA-256 input is a `seed-v2` domain tag followed by the typed,
length-prefixed text vector of protocol, input hash, candidate, acquisition,
hypothesis, decimal draw ID, and decimal collision nonce; every valid
acquisition label is therefore unambiguous. For each sorted draw ID, its seed
starts a fresh
Mersenne-Twister/Inversion/Rejection stream. An independent proposal calls
`sample.int(n_g)` once per radix group and maps the group to that sampled order;
a blocked proposal calls `sample.int(2, B, replace=TRUE)-1` and composes the
selected disjoint canonical block swaps. Reject identity or a previously
accepted row-index vector, then continue the same draw stream. Record the retry
count and SHA-256 of the big-endian integer map. More than 1,000,000 retries for
one draw is `association_numerical_failure`. This is sequential rejection over
the finite product map space; unique seeds alone are never treated as proof of
unique transforms. The result seed manifest retains rows only for available
Monte Carlo outcomes, although collision resolution always used the complete
pre-map manifest. RNG kind, seed existence, and seed bytes are restored exactly.

## Empirical quasibinomial candidate

This candidate is licensed only for independent units. Fit grouped counts
`(k_i,G-k_i)` in rank coordinates with `stats::glm.fit`,
`quasibinomial()`, `start = etastart = NULL`, the family default `mustart`, and
`glm.control(epsilon=1e-10,maxit=100,trace=FALSE)`. Final weighted SVD rank must
remain `r`.

```
W_i = G mu_i(1-mu_i)
K_q = (Z'WZ)^-1
phi = sum_i (k_i-G mu_i)^2/[G mu_i(1-mu_i)]/(n-r)
V_Q = phi K_q
```

Require every `0 < k_i < G`, fitted probabilities in
`[sqrt(.Machine$double.eps), 1-sqrt(.Machine$double.eps)]`, finite convergence,
and finite `phi > 0`. Symmetrized `V_Q` and its scalar delta variance then pass
the same `tau(A)` PSD/variance checks as every other covariance. The tested raw
p-value uses the log-odds direction `c_j'a_hat` with residual-df t reference.

Primary effect toggles only encoded column `j`:

```
z0_i = z_i - x_ij c_j
z1_i = z0_i + c_j
mu0_i = plogis(z0_i'a_hat)
mu1_i = plogis(z1_i'a_hat)
Delta_j = mean_i(mu1_i-mu0_i)
g_j = mean_i[mu1_i(1-mu1_i)z1_i - mu0_i(1-mu0_i)z0_i]
SE(Delta_j) = sqrt(g_j' V_Q g_j)
```

The 95% delta interval is not clamped. Retain log-odds effect/interval
secondarily.

## Generator oracle and candidate screening

For each replicate, use its realized `F+` but average the generator's marginal
`1-missing_probability` over that set. This is a predeclared plug-in oracle -
not `E[y|F+]`, because selection into the realized set is not conditioned in
the marginal probabilities.

Fixed targets are:

- `dda_batch_crossed` - `batch[batch_2]` versus `batch_1`;
- `dda_monotone_unequal` - `condition[D]` versus `A`;
- `dia_batch_partial` - `batch[batch_3]` versus `batch_1`;
- `dia_monotone_paired` - `condition[B]` versus `A`.

OLS and Freedman-Lane truth is the target coefficient from the rank-coordinate
OLS projection of plug-in sample fractions. Quasibinomial truth is its
standardized 0-to-1 effect from a quasibinomial projection of fractional
plug-in counts. These are candidate-specific projection calibrations on one
response; the contract does not pretend their effect estimands are identical.

Screening uses the M12 thresholds inside each candidate's licensed scope, with
in-scope failures retained as misses. OLS and Freedman-Lane screen on 64 null
families + 128 targets. Quasibinomial screens on the same 64 independent null
families + 96 independent targets; paired DIA is out of its performance
denominator but remains a miss in common scope coverage. Every screened null
family must contain at least one available p-value - zero-p abstention rejects
screening instead of improving false-flag error.

Every screening-qualified candidate then faces the full fixed v3 synthetic
registry: 64 null families and all 128 targets. A zero-p null family,
unavailable target, or numeric miss eliminates that candidate. Only full-gate
passers enter winner ranking, so failure of a narrow candidate cannot suppress
a passing broader candidate. The panel is killed/parked only when the
full-passer set is empty. No successful-fit denominator is allowed. Thus a
partial-scope method may screen well but cannot produce a release candidate by
hiding unsupported rows.

`m13a_rank_candidates()` freezes a total order. Repeatedly anchor the remaining
candidate with greatest common opportunity coverage, breaking an exact tie by
candidate ID. Its comparison cluster contains candidates within 0.02 of that
anchor in each named metric: common opportunity coverage, Holm target power,
and median absolute candidate-specific projection bias. Choose the cluster's
lowest frozen simplicity rank, then median runtime, then candidate ID; remove
it and repeat. Anchoring each iteration prevents nontransitive pairwise
near-tie decisions.

## Versioned gates and bookkeeping correction

The 11 completed design-core gates retain their exact M12 IDs and binding
hashes. The nine association/public rows receive new `_v3` IDs in
`m13a_gate_registry_v3`, full numeric rows, v3 protocol/data hashes, lineage to
the original M12 gate binding, stage, replicate range, and a new binding hash.
`m13a_evaluate_gate_results()` accepts exactly one complete evidence stage and
rejects omissions, detached bindings, malformed or inverted interval endpoints,
and unavailable endpoints carrying numeric values. The three calibration
comparators name the candidate-specific plug-in projection rather than the
historical M12 expectation wording.

The M12 confirmation comparator included acquisition as a tested term, which
contradicts acquisition-specific models. Its v3 successor keeps the same
threshold, sample size, data, claim, and failure treatment but tests
per-acquisition condition/instrument/nuisance coefficients and labels
acquisition descriptive. The historical M12 row remains unchanged.

M12 also omitted Track A from `dda_monotone_unequal.required_tracks` despite
using that scenario in all three A alternative gates. The v3 effective manifest
sets `A;B;C`; v3 evidence hashes derive from that manifest. Historical M12
manifest/data hashes remain unchanged and are retained as upstream lineage.

Candidate association gates bind exact replicates 1-32. After a winner and its
evidence hash lock, association computations on replicates 33-64 and HarmonizR
development may open. MultiPro confirmation remains M15P-only. M13 closes with
a development disposition, not a claim that confirmation ran.

## Runtime and execution seals

The association slot is a sum type: either exact
`sentinel_association_v1` or an exact top-level `imputefinder_unavailable` with
quantity `association` and code `association_no_observable_features` or
`association_no_testable_hypotheses`. Thus `G=0` and pre-hypothesis abstention
are representable without fabricating a positive denominator.

`sentinel_association_v1` is an exact
`imputefinder_association_panel` with protocol/method records, response,
hypothesis/support tables, one keyed available/unavailable outcome per
hypothesis, Holm summary, and diagnostics. The executable catalog freezes field
order, type, cardinality, key, class, nesting, NA rules, canonical order, and
schema values. `m13a_validate_result_shape()` checks joins, available and legacy
unavailable identities, finite evidence, support keys, exact Holm arithmetic,
and diagnostic/seed schemas on synthetic fixtures. Its protocol record binds
the v3 contract and gate-registry hashes, locked implementation-manifest hash,
candidate-evidence hash, winner state, and historical M12 hashes. A merely
candidate-shaped object without those winner bindings is invalid.

Legacy sentinel objects with `pre_rescue` alone or `pre_rescue + coverage`
remain readable. New production objects add `association` only after method
lock; unknown association schema fails, and candidate comparison tables never
enter the public object.

Before any candidate computation, freeze a separate implementation manifest:
exact engine/evaluator/schema/test/harness source hashes,
R/platform/package/loaded-nonbase-dependency environment, focused-test
evidence, and v3 protocol/contract/gate/effective-manifest hashes. The exact
nested schemas and `m13a_evidence_bundle_hash_v1` component/self-hash algorithm
are implementation rails whose validators must turn green before any result
allocation. Candidate evidence contains the complete canonical 3-candidate x
13-scenario x 32-replicate run grid, outcome dispositions, six screening
metrics per candidate, the complete selected-winner-stage gate set for every
screening-qualified candidate, full-passer ranking metrics/order, and the
selected winner or no-winner state. Run hashes cover internal candidate
artifacts without the later public winner wrapper; only after evidence locks is
the public panel materialized with that evidence hash. This one-way binding
avoids a self-referential result/evidence hash. Any drift rejects execution and
requires a new pre-result implementation seal.

## Verification

```sh
Rscript --vanilla dev/m13-association-contract.R --verify
```

The verifier checks upstream identities, normative source/detail digests,
component and aggregate hashes, effective-manifest/gate lineage, gate evaluator,
allocation/stage seals, global seed collision handling, exact runtime joins and
Holm arithmetic, mask/design byte invariance/sensitivity, and mutation
rejection. It computes no candidate and reads no external result artifact.
