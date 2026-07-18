# M12b validation contract - executable generators + sealed data roles

Status: generator/data-role freeze complete. Every external result artifact
remains unopened; A-C numeric gates and candidate protocols remain the next
freeze.

Run from the repository root:

```sh
Rscript --vanilla dev/m12-validation-contract.R --verify
Rscript --vanilla dev/m12-generator-validation.R --verify
R_LIBS_USER="$PWD/.agent/R-library" \
  Rscript --vanilla dev/m12-generator-validation.R --v1
```

The contract verifier requires exact schemas, generator linkage, assigned
whole-family roles/split hashes, licence and upstream identities, and this
state:

```text
gate rows = 0
data/artifact roles = development | confirmation
artifact open state = metadata_only
local artifact SHA-256 = absent
```

Six adversarial rails reject an out-of-scope generator, missing split hash,
card/artifact role disagreement, opening without local SHA-256, malformed
upstream checksum, and drift from the exact UPS1 four-level subset. This proves
protocol integrity - no A-C candidate has been run.

## Frozen hashes

Canonical frame encoding = sorted first key, ordered fields, UTF-8, explicit
types/`NA`, tab/newline escaping; digest = SHA-256. Aggregate hashes cover
ordered component names + hashes.

| Contract component | SHA-256 |
|---|---|
| Schema catalog | `24e4aa3c93f5b0d87ff5ed4f92e70b08d0e258359a57b7c88fd1bd931c8b18f3` |
| Generator manifest | `9d1381833cd59e8775bd99e730d5c783b98c5a82afff6d52276e83a46a8be76a` |
| Generator protocol linkage | `9f2166bc7a31112d5529d81adcd68cae70b99d12374f3d77bc45844d3d16bd7a` |
| Claim inventory | `c6265115e6e2abb9ed30ef52a93d1605427de6869beed70f642a358bb9e3ce40` |
| Empty gate registry | `0f10cd353879c95a851be3925a3c2f9570179cdf4592134215cccfdc58d80305` |
| Data cards | `ea8605fe7f0b5c0765cb13eecb9851fee53ea527038984a44b8839d20d576188` |
| Data roles | `650b6492d2f2a8d5542bc7097602c67f1f68db71d706a18e0bfd4dcfb4637aa7` |
| Artifact manifest | `7c3267020d484ad9bc1bf4b8ac442cf0916d554e9e467e5c696ab2a20fb1c652` |
| Candidate exclusions | `7ef14dd357171ffaad2e3793dfd1e038d1cc4108d56d052a2f751e4061e71af8` |
| Source manifest | `b56898ee3ce5cd8ad2f2656a0b1cfb94792d8f1fb7f51b54472555279d8033ed` |
| **Aggregate contract** | **`bb5624fecd9cd52980bd901c1b75eac0d58dcb1869d95f92014e0d9e086a53ed`** |

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

Intentional changes require version/hash updates before linked results are
inspected.

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
| MS-DAP `PXD036134` | development | B/C controlled recovery | all 12.5/15.625/18.75 ng yeast mixtures; linked DDA+DIA | `189611f39df498f3f99470b24a4d21d09703f58a141594c9d21cff2fb8dc7e62` |
| UPS1/yeast `PXD002099` | confirmation | B/C DDA | exact near-log-spaced 2/4/10/25 fmol subset; 50 fmol excluded | `a90bdc8c98e5ed33593f8b949521e2c8e78bbc130caf3f1857afc593edc498eb` |
| LFQbench `PXD002952` | confirmation | B/C DIA | HYE124 A/B; inventoried 6600/64-variable-window family | `02860b2f126937b89742f67b6fac2fe971b6d9b6267c8b9508e9a636f725a8cd` |

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
downloaded in M12b.

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

PXD041421 remains related to MultiPro; PXD028735 exposes only >700 raw files in
the checked metadata; the rheumatoid-arthritis SWATH case still lacks a stable
deidentified/licensed protein artifact; PXD003881 still lacks a compact frozen
protein derivation. Their explicit reconsideration conditions remain in the
contract.

## Exact next slice

M12c populates every retained A-C claim with complete numeric gate rows and
failure rules, then freezes A/C estimator candidates and B perturbation units,
counts, weighting, seed streams, cutoff-range/failure policy, and display
calibration. External result bytes remain unopened until their linked protocol
hashes pass; confirmation stays sealed through M15P.
