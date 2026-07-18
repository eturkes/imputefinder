# M13b static sentinel coverage

Status: deterministic static coverage schema green; A association candidates
remain `frozen_unrun`; every external artifact remains `metadata_only`.

## Contract

Selected sentinel output now contains the unchanged
`pre_rescue_evidence_v1` record plus `sentinel_static_coverage_v1`:

- `sample` - canonical sample/condition order, full + globally observable
  feature denominators, detected-feature fraction, observed mean, and observed
  intensity minimum/median/maximum;
- `condition` - sample/full/eligible support, detected + complete features,
  observed/eligible cells, and detection fraction;
- `role_level` - every declared condition/nuisance/block/acquisition level,
  encoding, numeric value where applicable, sample/condition support, and
  singleton status;
- `condition_role` - full condition product with every declared nuisance,
  block, and acquisition level, retaining zero-sample cells with zero counts +
  `NA` detection fraction;
- `feature_overlap` - every unordered condition pair with shared/union/side-only
  and neither-detected counts plus Jaccard overlap.

Globally observable = finite in at least one original cell. It is the frozen A
association-response denominator; the original feature count remains explicit.
All summaries precede rescue/classic fitting. Numeric nuisance values retain
doubles plus locale-independent labels; categorical labels, roles, samples,
features, and conditions use canonical ordering.

The validator regenerates every structural count/table from the packed original
mask + typed design. It cross-checks stored sample means against pre-rescue
evidence and enforces finite ordered intensity support; exact extrema cannot be
regenerated without the deliberately omitted original numeric matrix. Existing
experimental objects containing only `pre_rescue` remain readable; newly
constructed objects always contain coverage.

## Verification

- Red baseline → 57 failed expectations + one lifecycle error at the absent
  coverage seam; final focused suite → 8 tests / 81 expectations, all green.
- Package-wide source suite → 152 tests / 1,077 expectations, no failures,
  errors, warnings, or skips.
- Named feature/sample/role order, factor/character and integer/double
  re-encoding, singleton/empty cells, zero-support samples, wholly absent input,
  schema corruption, legacy readability, original-object immutability, and
  matrix/`SummarizedExperiment` integration pass.
- Frozen M13 mandatory verifier → all 11 gates pass across 416 development
  instances; audit/result hashes remain
  `37528216590dcfff3338501ef61b766697ea63c14b7324912fcb00f84d58790e` /
  `cb72a69301100287cee569ae7a0d709f42edc10e3b10a62dbff3a86ab88aec6f`.
- Frozen M5 → 12 scenarios / 13 target profiles pass; frozen M7 → 13 long +
  routine scenarios pass, MD5
  `4011e381bba2d0d747e91d277a45de5e`.
- Uncontended 10,000 x 50 sidecar gate → median 2.257 s = 4.57x stable;
  allocation 622.85 MiB = 3.05x stable; largest allocation .90x input; object
  size 1.17x embedded classic; every exact/input/mask/comparator/runtime/
  allocation/object-size gate passes.
- Vignette-bearing source build passes; source-tarball
  `R CMD check --as-cran --no-manual` has only the expected new-submission note;
  staged-tree `BiocCheckGitClone` = 0 errors/warnings/notes; new-package
  `BiocCheck` = 0 errors/warnings + three reviewed notes. Automatic `Coverage`
  biocView suggestion was rejected because this package reports design/data
  support, not genomic/read coverage.

No A candidate was fitted, no candidate/development result allocation changed,
and no external result byte was downloaded, listed, or parsed.

## Next work

M13c implements the three frozen A association candidates + support/abstention
schema, then opens only synthetic candidate replicates 1-32 for gate-first
selection or a no-winner result. Untouched development replicates 33-64 and the
HarmonizR development artifact remain closed until that choice locks;
confirmation remains sealed through M15P.
