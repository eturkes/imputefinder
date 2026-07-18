# M12a validation contract - scope, schemas, candidate inventory

Status: frozen metadata-only contract. Candidate results remain unopened; data
roles, split identities, numeric gates, and estimator protocols remain
unassigned.

Run from the repository root:

```sh
Rscript --vanilla dev/m12-validation-contract.R --verify
```

The verifier requires exact schema/order/type vocabularies, Section 9.1 scope
coverage, valid upstream SHA-1 identities, legal-status records, grouped split
rules, and the M12a state:

```text
gate rows = 0
data/artifact roles = unassigned
artifact open state = metadata_only
local artifact SHA-256 = absent
```

It also proves four rejection rails: a five-condition generator, role assignment
without a split hash, artifact opening without role/local hash, and malformed
upstream checksum all fail. This is a protocol guard, not evidence that any A-C
candidate works.

## Frozen hashes

Canonical frame encoding = sorted first key, ordered fields, UTF-8, explicit
types/`NA`, tab/newline escaping; digest = SHA-256. The aggregate hash covers
the ordered component names + hashes.

| Component | SHA-256 |
|---|---|
| Schema catalog | `70e8953371cd81f2d189e08c975a201dee8420236daa9df18cc43b074ad7e8a6` |
| Generator manifest | `aba19431f8c4f4a51f29df8471596b274d157f56f5dd6fb9218d8c09907f1e49` |
| Claim inventory | `c6265115e6e2abb9ed30ef52a93d1605427de6869beed70f642a358bb9e3ce40` |
| Empty gate registry | `0f10cd353879c95a851be3925a3c2f9570179cdf4592134215cccfdc58d80305` |
| Data cards | `b27e2c64e6e2190ead46e6534616444bb83875ac1ddd3de7d778916363c46ecd` |
| Artifact manifest | `58cd567e81e36e2c9af611bf65596b49e43a45524172e2d6116fc80e14153ddd` |
| Candidate exclusions | `7ef14dd357171ffaad2e3793dfd1e038d1cc4108d56d052a2f751e4061e71af8` |
| Source manifest | `65f29787525ff4420c1820ecab86cd8879b4402ebe7d5f3a35492a7ef470543b` |
| **Aggregate contract** | **`3f9a53d73495b0b7e3ce9ccf6a36ddcbabe4928f39cd5297e9a7cf79848b9e77`** |

Any intentional change requires a version/hash update before results are
inspected.

## Generator extension

Thirteen specified-unimplemented Tier 1-2 scenarios cover only `PLAN.md`
Section 9.1:

- DDA and DIA remain separate synthetic parameter strata. Their labels do not
  claim to simulate acquisition physics.
- Conditions = 2-4; balanced and unequal independent designs plus paired
  subject blocks are present.
- Declared run/batch effects include crossed, partial-confounding, and
  perfect-aliasing cases.
- Missingness patterns cover intensity-independent, monotone-abundance, mixed,
  blockwise, batch-associated, structural-off-compatible, and no-cliff cases.
- Stress axes include outliers, heteroscedasticity, differing intensity support,
  adequate/low evidence, a flat-profile control, a causal-wording control, and
  a grouped-sibling leakage trap.

The manifest freezes coverage and scenario identities, not numeric generator
parameters or outcomes. Those belong to the next protocol-freeze slice.

## Executable schemas

The script owns six relational schemas plus a source manifest:

```text
generator_manifest -> bounded simulation/semi-synthetic scope
claim_inventory     -> every mandatory/independent A-C claim family
gate_registry       -> complete numeric rows only; partial rows invalid
data_cards          -> design, unit, relation, licence, truth, exclusions, split
artifact_manifest   -> exact file/size/upstream hash + role/open state
candidate_exclusions-> negative inventory with reconsideration condition
```

A gate row is admissible only when it has the full Section 9.4 tuple: claim,
metric, unit, comparator, strata, positive sample size, uncertainty method,
finite numeric threshold/operator/scale, failure treatment, known data IDs,
protocol ID, and SHA-256 data/protocol hashes. The registry is intentionally
empty in M12a: inserting provisional metrics or thresholds would falsely imply
they were frozen.

A candidate artifact can move from `metadata_only` to `opened` only after its
card and artifact share a non-`unassigned` role, the card has an exact split
SHA-256, and the downloaded bytes have a local SHA-256. Role assignment must
therefore precede result inspection.

## Minimum candidate inventory

The inventory uses only repository/project metadata and primary papers. No
result archive/CSV bytes were downloaded or parsed, and no archive contents
were listed.

| ID | Acquisition/design | Identified truth | Split constraint | Licence |
|---|---|---|---|---|
| `multipro_hcc_pxd041391` | DDA + DIA; 2 cell lines; 3 biological replicates; 3 technical replicates; 2 instruments; deliberate batch | Cell class, instrument batch, replicate hierarchy | Biological replicate is the smallest split unit; all injections/modes/instruments stay together; only 3 units/class means the family cannot supply multiple evidence roles | PRIDE record: CC0 |
| `msdap_hy_pxd036134` | Matched DDA + DIA; 3 yeast spike levels; triplicates | Fixed human background; 12.5/15.625/18.75 ng yeast | Technical repeats only; DDA/DIA are linked; whole family gets one role | PRIDE record: CC0 |
| `ups1_yeast_pxd002099` | DDA; 5 UPS1 levels; triplicates | Fixed yeast background; known UPS1 inputs | Whole artifact gets one role; exact <=4-level subset or pairwise contrast must be frozen because the full design exceeds initial scope | PRIDE record: CC0 |
| `lfqbench_pxd002952` | DIA; inventoried HYE124 TripleTOF 6600 64-variable-window analysis; A/B triplicates | Known species mixture/ratios | Master mixtures, instruments, software outputs, and triplicates are related; whole family gets one role | EMBL-EBI terms; owner-rights caveat retained |

MultiPro supplies the needed multibatch public design without pretending that
its naturally missing cells have causal labels. MS-DAP supplies matched modern
DDA/DIA concentration evidence but no biological replication. The compact
UPS1 and LFQbench candidates offer independently sourced acquisition strata;
neither may be split internally to manufacture confirmation.

### Artifact identities

The PRIDE API file `checksum` values are SHA-1; PRIDE's checksum utility defines
that algorithm. SHA-1 is retained only as the upstream integrity identity.
After role freeze, local downloads receive SHA-256 before opening.

| Artifact | Bytes | Upstream SHA-1 |
|---|---:|---|
| `HCC1806_HS578T_DDA_FragPipe_Output.zip` | 2,942,612,390 | `4e372285a6ee8c8301ee3f1e7ff9bc548e0ee4bd` |
| `HCC1806_HS578T_DIA_DIANN_Output.zip` | 3,083,088,287 | `d027d5f89b2eea04f6c90b8189760856f70f039b` |
| `timstofpro2_30SPD_two-proteome_spike-in_series_DDA.zip` | 101,116,244 | `9113637ddfdcd9c97d10ab0358bc918e948eee19` |
| `timstofpro2_30SPD_two-proteome_spike-in_series_DIA.zip` | 127,676,636 | `f03cb375f06ba243f73d02580881ee44a9ccd734` |
| `YEAST-Data-NonNormalized.csv` | 273,919 | `90e76ee3fa6238ba5325aac7d82bd2b8b6bcafb0` |
| `lfqbench_analysis_HYE124_TTOF6600_64var.zip` | 263,968,636 | `a4420dd9f5c5833371b5623cebd6742d36888241` |

## Excluded candidates

| Candidate | Current exclusion | Reconsideration condition |
|---|---|---|
| MultiPro A549/K562, PXD041421 | Same laboratory/protocol family as inventoried MultiPro; processed-archive checksums absent in checked metadata | Non-independent exploratory role + immutable artifact hashes |
| Multiplatform LFQ, PXD028735 | More than 700 runs; no minimum checksummed processed protein artifact selected | Platform-shift claim promoted + exact <=4-condition artifact/protocol frozen |
| Rheumatoid-arthritis SWATH study | Strong blocked longitudinal design, but no stable deidentified protein artifact/accession + reuse licence identified | Artifact, subject/batch metadata, privacy, and licence verified |
| IonStar, PXD003881 | Strong DDA mixed-species truth, but no compact deposited protein matrix | Versioned protein derivation + preprocessing/hash frozen |

An exclusion is a present feasibility result, not a judgment on scientific
quality.

## Primary metadata and legal sources

- [PRIDE API](https://www.ebi.ac.uk/pride/markdownpage/prideapi) and
  [PRIDE checksum utility](https://github.com/PRIDE-Archive/pride-checksum)
- [EMBL-EBI terms of use](https://www.ebi.ac.uk/about/terms-of-use/)
- MultiPro [paper](https://doi.org/10.1038/s41597-023-02779-8) and
  [PXD041391 record](https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD041391)
- MS-DAP [paper](https://doi.org/10.1021/acs.jproteome.2c00513) and
  [PXD036134 record](https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD036134)
- UPS1 series [paper](https://doi.org/10.1021/acs.jproteome.5b00183) and
  [PXD002099 record](https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD002099)
- LFQbench [paper](https://doi.org/10.1038/nbt.3685) and
  [PXD002952 record](https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD002952)

Checked 2026-07-18. Repository terms can change; role assignment rechecks the
record and freezes the applicable text/version.

## Exact next slice

M12b must freeze generator parameters/seeds and the minimum whole-family data
roles, including any <=4-level subset and split/protocol hashes, while keeping
all result artifacts unopened. It then populates every retained A-C claim with
complete numeric gate rows and explicit failure treatment. Only after those
hashes pass may a development artifact be downloaded, locally SHA-256 hashed,
and opened. Confirmation artifacts remain sealed through M15P.
