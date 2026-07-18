# M12 validation contract - frozen generators, evidence, gates, and candidates

Status: M12 complete. All generators, evidence roles, numeric gates, candidate
comparisons, perturbations, support rules, and calibration procedures are
hash-frozen before results. Every external artifact remains unopened.

Run from the repository root:

```sh
Rscript --vanilla dev/m12-validation-contract.R --verify
Rscript --vanilla dev/m12-generator-validation.R --verify
Rscript --vanilla dev/m12-candidate-protocol.R --verify
R_LIBS_USER="$PWD/.agent/R-library" \
  Rscript --vanilla dev/m12-generator-validation.R --v1
```

The contract verifier requires exact schemas, generator linkage, assigned
whole-family roles/split hashes, licence/upstream identities, protocol linkage,
complete claim coverage, executable gate bindings, and this state:

```text
claims = 27 (A=9; B=9; C=9)
numeric gate rows = 66 (A=20; B=18; C=28)
candidate protocols = frozen_unrun
data/artifact roles = development | confirmation
artifact open state = metadata_only
local artifact SHA-256 = absent
```

Twenty contract/gate rails prove nonempty staged subsets and reject scope,
split, role,
opening, checksum, subset, claim-coverage, numeric-threshold, data/protocol
binding, result-seal,
omission, detached-result, unavailable, fabricated unavailable endpoints, and
wrong-side or altered-registry failures. Eight candidate rails reject allocation,
replacement, cutoff-floor, multiplicity, result-opening, seed-key, and
C-resampling drift plus empty seed manifests. This proves protocol integrity -
not that a candidate passes.

## Frozen hashes

Canonical frame encoding = sorted first key, ordered fields, UTF-8, explicit
types/`NA`, tab/newline escaping; digest = SHA-256. Aggregate hashes cover
ordered component names + hashes.

| Contract component | SHA-256 |
|---|---|
| Schema catalog | `3c117b427cda433a58548855d25ebe7b6bed0f235cbae7a12aab69a3b2e1689c` |
| Generator manifest | `9d1381833cd59e8775bd99e730d5c783b98c5a82afff6d52276e83a46a8be76a` |
| Generator protocol linkage | `9f2166bc7a31112d5529d81adcd68cae70b99d12374f3d77bc45844d3d16bd7a` |
| Claim inventory | `c31032238677acc5e57478057852481f29dd24bc3f6fb88d629653be5ba8d296` |
| Candidate protocol registry | `61f7c1d9c42e80e1dafb13693e7a1818c7835ee9fd219e8a953167fe0c84ab6d` |
| Internal evidence | `7c2deb16729dc06720bea8ce14b1ed4024dcc51f836c03e00b3ccc7696464a1a` |
| Numeric gate registry | `12e1fb5f9886b4ca6a4accb18909b3aa40cf33ac220800d970d79d456d9fcd87` |
| Data cards | `0ba7fa86286afd4cc72dfea6eda4b0bc1324553229c2a3b7c412d41c777d684c` |
| Data roles | `500d08b89b44c3adde9782176ae1065bf20c2eb135b89e666c8623ec2c83665c` |
| Artifact manifest | `7c3267020d484ad9bc1bf4b8ac442cf0916d554e9e467e5c696ab2a20fb1c652` |
| Candidate exclusions | `7ef14dd357171ffaad2e3793dfd1e038d1cc4108d56d052a2f751e4061e71af8` |
| Source manifest | `e413eb42bb5f173dde1f994664b2dd9b081eee304323cea3cf8ac390ac2651b4` |
| **Aggregate contract** | **`2cd5ecd8bf5763da1c5d9e1d8994207e701bcea9a7efcf4ae539dfd9b52d7431`** |

Generator subprotocol:

| Component | SHA-256 |
|---|---|
| Descriptor | `8076ad5a990841a1fcd37d5145b3f01eaae9cb9f9dd2c922712c150bc221d132` |
| Acquisition profiles | `d868a3bdf9d2daaf2c9ed482833416cbbd07b83956a205056bc2d371af731fbe` |
| Scenario parameters | `019a81a0438958e7080243955e1e6673999a29e82fbe2ed70b5d8fc0d100b86e` |
| 832-row seed manifest | `c3f550c378b3d679e44b472ac01178776ef7e3f386ee4d5a43dde928d752de79` |
| **Aggregate generator protocol** | **`cdea1bf874152e63fba08c49e390da1dd18690105e43a2dc7a03fc6866d0d080`** |
| Two-seed generator audit | `4d216966c8cc54b1b23685d8801c7c65b2afcc5879dd15430d57a88f765e8a00` |
| Expanded stable-v1 audit | `8463f4306fe9b6192336487e4ca365c8c7e3525dc06aeca5b66b1a99b1a10250` |

Candidate/perturbation protocols:

| Track | Protocol | SHA-256 |
|---|---|---|
| A | Design core + association comparison | `e20bdcc8e040eed65457480be7ae2ce2542761165c3d27161b953f86f39e5edc` |
| B | Perturbations + cutoff/display calibration | `57ffa1220c740dcb69508d3ad1b9beb05e7013c3a92af528293211e41a13f256` |
| C | Detection + observed-abundance comparison | `58b6e3c6d7f35f5561549a3a52420e49b0f3eb296752353d1bb3c9d433cf9d96` |
| **Bundle** | ordered A/B/C hashes | **`8eb5592131407ae48d4b4caa06ef13eea84371471215eae41997faf0e5ec0821`** |

Intentional changes require version/hash updates before linked results are
inspected.

## Frozen evaluation and numeric gates

Each generator's 64 streams now have a second, sealed allocation: replicates
`1-32` = candidate selection/display calibration; `33-64` = untouched
development test. A winner and every calibrated threshold lock before the
second half opens. Confirmation data never participate in selection.

Every Section 9.4 row includes claim, metric, unit, comparator, strata, sample
size, uncertainty, numeric operator/value/scale, failure disposition, canonical
data IDs, evidence hash, and track-protocol hash. Later result rows must match
the frozen gate-row hash; omission, `unavailable`, candidate failure, or the
wrong side of a threshold fails. Failed/unavailable rows carry `NA` endpoint
fields, so no synthetic finite measurement can conceal abstention. High-level
bounds:

| Track | Frozen gate families |
|---|---|
| A mandatory | rank/alias/estimability/block/order/state exact `1` or mutation/drift count `0` |
| A panel | null upper bound `.15`; per-acquisition false flags `<=.125`; 95% coverage `>=.90`; power `>=.70`; median absolute probability bias `<=.03`; public completion `>=.95`; wording violations `0` |
| B | cutoff 95% range coverage `>=.90`, median error `<=.35` log2; no-cliff false-confidence upper `<=.10` and abstention `>=.95`; sparse/reproducibility/side-effect exact rails; stable runtime ratio `<=1.05`; optional display held-out risk upper `<=.15`, coverage `>=.10` |
| C | separate detection/abundance null FDR upper `<=.10`, false-sign upper `<=.10`, power `>=.50`, precision `>=.90`; median detection bias `<=.075` probability points; abundance bias `<=.25` log2; 95% coverage `>=.90`; edge-case exact rails; controlled-reference errors/directions; wording violations `0` |

These are frozen falsification criteria, not observed performance. Stratum
distributions remain required; pooled values cannot hide a DDA/DIA miss.

## Frozen candidate protocols

A uses the original mask and one common denominator: proteins observed at least
once globally. It canonicalizes the declared design, computes SVD rank/null
space and row-space estimability, and derives a rotation-invariant null basis
from the null projector + ordered projected axes, then compares three
sample-detection-fraction
association candidates: OLS with independent HC3 or paired/repeated CR2,
studentized restricted Freedman-Lane, and independent-unit quasibinomial.
Per-side/block/residual-df and permutation-resolution floors precede fitting.
Freedman-Lane freezes exchangeability, two-sided tails, exact enumeration, and
the `9999`-draw Monte Carlo correction. Holm `0.05` controls the
acquisition-specific declared-term family. Candidate selection is gate-first,
then supported scope, then explicit simplicity/noninferiority; there is no
silent fallback.

B keeps three uncertainty families separate:

- sampling = 999 condition-stratified biological/block bootstraps + exhaustive
  leave-one-unit-out, baseline cutoff fixed;
- estimator = 999 feature/module bootstraps, automatic cutoff re-estimated;
- assumption = six named fixed/re-estimated/cutoff/summary/majority policies.

Multiplicities become internal integer weights - duplicate public matrix names
never reach v1. Successful-cutoff floor = `900/999`; range = type-8
`.025/.975` quantiles. Public ranges are descriptive. Sampling/estimator labels
use draws `1-499` for a score and `500-999` for held-out risk; a label ships only
if calibration and untouched-development risk/coverage pass. Assumption output
stays continuous because six deterministic policies supply no independent draw
law. Track-specific SHA-256 keys over the frozen protocol/input/stream/instance/
draw identity seed local restored base-R RNG scopes; open addressing operates
over the complete manifest for each stream. A display-gate miss disables labels;
it cannot kill otherwise passing continuous B panels.

C freezes two separate estimands: standardized model-adjusted detection-risk
difference and design-adjusted observed log2 contrast conditional on detection.
Detection truth is the empirical-row standardized generating-probability
contrast; conditional-abundance truth is the corresponding
`sum(p * intensity) / sum(p)` contrast. Gates derive component alternatives
from nonzero component truth: five detection scenarios and seven abundance
scenarios, so an abundance-only generator cannot masquerade as detection power.
No combined decision exists. Detection candidates = GLM+HC3, mean bias-reduced
logit, small-sample-corrected GEE, plus exact conditional logistic as an
odds-only paired comparator. Paired GLM/bias-reduced intervals use `999`
whole-block bootstrap draws from the frozen C stream. Observed-abundance
candidates = OLS+HC3 and ordinary/robust `limma` empirical-Bayes fits on finite
original cells; paired/repeated OLS uses CR2/Satterthwaite. Support floors,
separation, nonestimability, and component-specific abstention precede fitting.
BY `0.05` families are separate by contrast/component and count
unsupported original features conservatively as `p=1`.

Cell-level feature-fixed association, beta-binomial/quasibinomial feature
dispersion, and hierarchical Bayesian candidates were excluded before results
for estimand loss, absent replication, or uncalibrated-prior complexity. Their
reconsideration conditions remain executable protocol data.

## Executable generator

`dev/m12-generator-validation.R` implements exactly the 13 Section 9.1
manifest scenarios. It owns 64 explicit, unique, SHA-256-derived development
seeds per scenario and rescue seeds cycling `1/7/29`. RNG = local
`Mersenne-Twister` / `Inversion` / `Rejection`; caller kind/state is exact.
Audit replicates `1` and `64` freeze both ends of every stream.

Acquisition names select separate synthetic parameter strata only:

| Stratum | Mean | Feature SD | Detection midpoint | Logit slope | Measurement SD |
|---|---:|---:|---:|---:|---:|
| DDA | 13.0 | 2.2 | 12.0 | 1.25 | .22 |
| DIA | 13.2 | 2.0 | 11.8 | 1.50 | .16 |

These are stress profiles, not acquisition-physics estimates. Each simulation
stores complete latent intensities, amputated data, component/union masks,
cell-level generating probabilities, feature-condition effects, design truth,
and a rule-aligned descriptive oracle. Generator-specific truth never labels
naturally missing public cells.

Model components:

- stratified-normal feature baselines + 25-feature correlated modules;
- predeclared condition alternatives, subject effects, run-order drift, and
  crossed/partial/perfect batch layouts;
- monotone cell missingness
  `m_max * logistic(slope * (midpoint - latent intensity))`;
- independent cell loss, abundance-associated whole-batch/subject loss,
  structural-off-compatible forced blocks, heteroscedastic noise, and outliers;
- missing probability = complement of component survival; realized mask =
  component union; structural blocks force probability `1`.

Null, flat/no-cliff, perfect-confounding, causal-wording, low-support, and
grouped-technical-sibling controls are executable. Validation checks exact
design counts/ranks, component eligibility, mask/probability identity, finite
support, paired coverage, structural blocks, grouped resampling units, and RNG
restoration. Across audit endpoints, realized missing fractions span `.0300` to
`.3460`; generator output hashes freeze every latent value, mask, design, and
truth record.

The stable-v1 expansion runs each first-seed matrix with its named generating
midpoint as a manual cutoff. Repeat calls, reversed feature/sample/condition
orders, and matrix/`SummarizedExperiment` paths have identical canonical exact
outputs; inputs, caller RNG, and originally observed values remain exact. The
13 cases retain `1,000-3,500` features, rescue `0-755` blocks, and reproduce
`0-176` global absences. These are compatibility observations, not accuracy or
mechanism claims.

## Frozen evidence roles

External training/tuning data = none; deterministic simulations own that role.
Every public family is indivisible across evidence roles.

| Dataset | Role | Linked track | Frozen inclusion | Split SHA-256 |
|---|---|---|---|---|
| HarmonizR mouse `PXD027467` | development | A public association | tumour/control; all four DDA preservation/timepoint batches; 25 animals | `e980a7de4f1f144e4ae086295008481beeaec6dd9659f630eda48da07d2725c8` |
| MultiPro HCC `PXD041391` | confirmation | A public association | HCC1806/HS578T; complete DDA+DIA instrument/replicate family | `da552730ff792b5e470c9031c433d112a557e515cb060c42c6266b60afa4bfba` |
| MS-DAP `PXD036134` | development | B external behavior + C controlled recovery | all 12.5/15.625/18.75 ng yeast mixtures; linked DDA+DIA | `ffb15df98f3bc56f8149068f956495eb8f782ac4559cb31ccfdc36d10ec30b54` |
| UPS1/yeast `PXD002099` | confirmation | B external behavior + C DDA recovery | exact near-log-spaced 2/4/10/25 fmol subset; 50 fmol excluded | `7ce5e2da20ba241610699f104caf52f0eedee0a8a6d69a58c174ec403b7dce1f` |
| LFQbench `PXD002952` | confirmation | B external behavior + C DIA recovery | HYE124 A/B; inventoried 6600/64-variable-window family | `62cbd5e44c28491ec60126dfb5e7a7501fa5705f6633be8c90270a201701692d` |

HarmonizR supplies a whole-family multibatch development case independently of
MultiPro, allowing MultiPro to remain sealed confirmation. Technical repeats,
acquisition modes, instruments, software outputs, and master mixtures never
manufacture independent splits. Shared confirmation opening retires the family
for every linked claim simultaneously.

Development artifacts may open only after their linked numeric gate and
candidate-protocol hashes freeze. Confirmation artifacts remain sealed through
M15P. A later download must match upstream SHA-1/size, receive local SHA-256,
and transition the manifest before parsing.

## Artifact identities

PRIDE `checksum` = upstream SHA-1 transfer identity. No bytes below were
downloaded in M12.

| Artifact | Bytes | Upstream SHA-1 |
|---|---:|---|
| `Mouse_Metadata.xlsx` | 12,361 | `8392a710eeeff2b801336697d2134ec74e85ea01` |
| `Mouse_FFPE_2018_SEARCH.zip` | 97,768,338 | `d6788a617342ffc676dcc7a81ce371e40e67e3b5` |
| `Mouse_FF_2018_SEARCH.zip` | 168,611,708 | `634314ff9fa1bf5eb4c53aa22fac91ecf863d4c4` |
| `Mouse_FFPE_2020_SEARCH.zip` | 184,092,139 | `8379e4c92ce2bcdabcc6f8d264efe38681d1a604` |
| `Mouse_FF_2020_SEARCH.zip` | 304,630,271 | `3fc9ee403270c9e6b147ee8ac279c3469f83d71d` |
| `HCC1806_HS578T_DDA_FragPipe_Output.zip` | 2,942,612,390 | `4e372285a6ee8c8301ee3f1e7ff9bc548e0ee4bd` |
| `HCC1806_HS578T_DIA_DIANN_Output.zip` | 3,083,088,287 | `d027d5f89b2eea04f6c90b8189760856f70f039b` |
| `timstofpro2_30SPD_two-proteome_spike-in_series_DDA.zip` | 101,116,244 | `9113637ddfdcd9c97d10ab0358bc918e948eee19` |
| `timstofpro2_30SPD_two-proteome_spike-in_series_DIA.zip` | 127,676,636 | `f03cb375f06ba243f73d02580881ee44a9ccd734` |
| `YEAST-Data-NonNormalized.csv` | 273,919 | `90e76ee3fa6238ba5325aac7d82bd2b8b6bcafb0` |
| `lfqbench_analysis_HYE124_TTOF6600_64var.zip` | 263,968,636 | `a4420dd9f5c5833371b5623cebd6742d36888241` |

## Sources and exclusions

Live PRIDE recheck on 2026-07-18 confirmed CC0 for PXD027467, PXD041391,
PXD036134, and PXD002099; PXD002952 still records EMBL-EBI terms. The terms
page is revised 2024-02-05: EMBL-EBI adds no owner-independent data restriction,
expects attribution, and retains the original-owner rights caveat.

Primary records: [PRIDE API](https://www.ebi.ac.uk/pride/markdownpage/prideapi),
[checksum utility](https://github.com/PRIDE-Archive/pride-checksum),
[EMBL-EBI terms](https://www.ebi.ac.uk/about/terms-of-use/), HarmonizR
[paper](https://doi.org/10.1038/s41467-022-31007-x) +
[record](https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD027467), MultiPro
[paper](https://doi.org/10.1038/s41597-023-02779-8) +
[record](https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD041391), MS-DAP
[paper](https://doi.org/10.1021/acs.jproteome.2c00513) +
[record](https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD036134), UPS1
[paper](https://doi.org/10.1021/acs.jproteome.5b00183) +
[record](https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD002099), and
LFQbench [paper](https://doi.org/10.1038/nbt.3685) +
[record](https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD002952).

The candidate review also froze primary/current references for R
[Holm/BY adjustment](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html)
and [quantile types](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/quantile.html),
[bias-reduced logistic regression](https://cran.r-project.org/package=brglm2),
[conditional logistic likelihood](https://stat.ethz.ch/R-manual/R-devel/library/survival/html/clogit.html),
[GEE](https://cran.r-project.org/package=geepack), and
[limma](https://bioconductor.org/packages/release/bioc/html/limma.html), plus
the primary [CR2/Satterthwaite paper](https://doi.org/10.1080/07350015.2016.1247004),
small-sample covariance, permutation, cluster-bootstrap, and risk-coverage
papers listed in the source manifest. Package availability did not select a
winner.

PXD041421 remains related to MultiPro; PXD028735 exposes only >700 raw files in
the checked metadata; the rheumatoid-arthritis SWATH case still lacks a stable
deidentified/licensed protein artifact; PXD003881 still lacks a compact frozen
protein derivation. Their explicit reconsideration conditions remain in the
contract.

## Exact next slice

M13a implements only the mandatory algebraic design core: canonical declared
model matrix, SVD rank/null-space alias report, contrast row-space
estimability, block/technical-unit accounting, and exact named-order/side-effect
rails. A association candidates and development artifact bytes remain closed
until that core is green. Confirmation stays sealed through M15P.
