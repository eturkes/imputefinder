# M13a mandatory design-core validation

Status: all mandatory frozen development gates pass; A association candidates
remain `frozen_unrun`; every external artifact remains `metadata_only`.

## Scope

Protocol `m13_design_core_validation_v1` executed the M12-frozen mandatory core
against all 13 synthetic DDA/DIA scenarios, using only development replicates
33-64 = 416 scenario-replicates. It additionally checked classic-branch parity
on the normative v1 fixture. This run consumed algebraic/design and compatibility
endpoints only: it fitted no A association candidate and computed no association
candidate comparison or external-data result.

The package implementation:

- canonicalizes sample and observed categorical-label order;
- encodes an intercept, condition and declared nuisance main effects, a declared
  block fixed effect, and only explicit interactions with treatment contrasts;
- retains exact term/coefficient maps and one named model row/weight per sample;
- applies the frozen SVD rank threshold
  `max(nrow, ncol) * .Machine$double.eps * d[1]`;
- forms `I - V_r V_r^T`, then derives a rotation-invariant null basis by scanning
  coefficient axes with modified Gram-Schmidt and a fixed sign convention;
- assesses named declared-level or coefficient contrasts by their row-space
  residual under the frozen `sqrt(.Machine$double.eps)` rule;
- maps a declared block to the independent resampling unit; samples sharing one
  block, including within-condition technical siblings, remain inseparable.

No automatic batch correction, interaction search, sample exclusion, causal
attribution, or association test occurs in this core.

## Frozen gate results

| Gate | n | Estimate | Required | Result |
|---|---:|---:|---:|---|
| rank-deficiency sensitivity | 32 | 1 | `== 1` | pass |
| full-rank specificity | 96 | 1 | `== 1` | pass |
| exact-alias sensitivity | 32 | 1 | `== 1` | pass |
| false exact-alias fraction | 64 | 0 | `== 0` | pass |
| non-estimable contrast rejection | 32 | 1 | `== 1` | pass |
| estimable contrast retention | 96 | 1 | `== 1` | pass |
| paired/block unit accounting | 96 | 1 | `== 1` | pass |
| unequal-replication accounting | 96 | 1 | `== 1` | pass |
| named order/re-encoding invariance | 416 | 1 | `== 1` | pass |
| caller-state mutation count | 416 | 0 | `== 0` | pass |
| classic exact-content drift count | 417 | 0 | `== 0` | pass |

The classic comparison nulls only each branch's call before exact comparison,
because direct `classify_missingness()` and sidecar construction necessarily
record different matched calls. Every other result field, class, attribute, and
value remains in the comparison.

## Evidence identities

- complete 416-row audit SHA-256:
  `37528216590dcfff3338501ef61b766697ea63c14b7324912fcb00f84d58790e`;
- registry-bound 11-row gate-result SHA-256:
  `cb72a69301100287cee569ae7a0d709f42edc10e3b10a62dbff3a86ab88aec6f`;
- upstream M12 generator protocol:
  `cdea1bf874152e63fba08c49e390da1dd18690105e43a2dc7a03fc6866d0d080`;
- upstream M12 validation contract:
  `2cd5ecd8bf5763da1c5d9e1d8994207e701bcea9a7efcf4ae539dfd9b52d7431`.

Reproduce after installing the current source package locally:

```sh
R_LIBS_USER="$PWD/.agent/R-library" \
  Rscript --vanilla dev/m13-design-core-validation.R --verify
```

The verifier regenerates every instance, binds each endpoint to its frozen M12
gate row, requires all 11 thresholds to pass, and rejects evidence-hash drift.

## Boundary and next work

This evidence validates the mandatory algebraic seam used by B/C; it does not
validate the independently disposable A association panel. M13b next adds
pre-rescue sample/condition/nuisance/block coverage and overlap summaries from
stored evidence. Only after those static summaries are green may the frozen A
association candidates run. Confirmation remains sealed through M15P.
