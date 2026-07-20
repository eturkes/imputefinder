# M13 association candidate validation

Status: verified `no_winner`; statistical association panel `killed`. The
mandatory algebraic design/coverage core remains active and unchanged.

## Frozen identities + allocation

- protocol = `e9d3b3799bb160d65feb0ae63b2e49a92bb4e8e005b89232ada96b83cfabed27`
- contract = `422954c88ed7683f58ac296cb4d9e785d61c3f80a276914235d4b285bbc16547`
- gate registry = `65f01b20e74b216399e17433df393f32e9bda189f73992ab419069cf7d5c343b`
- implementation manifest = `d821cae04b689c6ddb3b5843602e146e3bdecada485faf7afcd83544cac1cba4`
- artifact inventory = `9aefb8ec8d2906bbb6530d48a6e8778cf95fe207760aa31dc0c7ab324c1cd551`
- candidate evidence = `a5a242c3adeba3843cc9a91546543794f54b924785db831781fe8c4ae360cbc6`
- allocation = 13 scenarios x replicates `1-32` = 416 inputs; three candidates
  = 1,248 result artifacts. Development `33-64`, HarmonizR, and confirmation
  remained sealed.

## Candidate result

All 1,248 runs completed as measured panels. All candidates passed both null
negative-control screens, but each failed interval coverage, target power, and
the candidate-specific bias metric. No candidate screened into the fixed five-
gate stage; full-gate and ranking tables therefore contain zero rows by frozen
protocol.

| candidate | common coverage | null CP upper | max null stratum | interval coverage | target power | bias | disposition |
|---|---:|---:|---:|---:|---:|---|---|
| `a_fraction_ols_hc3_cr2` | 0.6667 | 0.1167 pass | 0.0625 pass | 0.7266 fail | 0.2188 fail | unavailable/fail | screened out |
| `a_fraction_freedman_lane` | 0.4583 | 0.1167 pass | 0.0625 pass | 0.4766 fail | 0.0625 fail | unavailable/fail | screened out |
| `a_fraction_quasibinomial` | 0.6250 | 0.0951 pass | 0.0312 pass | 0.6354 fail | 0.3125 fail | unavailable/fail | screened out |

Thresholds remained fixed: null CP upper `<=0.15`; maximum acquisition false-
flag fraction `<=0.125`; coverage `>=0.90`; power `>=0.70`; median absolute
bias `<=0.03`. No gate was weakened or reinterpreted.

## Failure + negative-control audit

The 2,592 hypothesis dispositions comprise 1,344 available and 1,248 structured
unavailable outcomes. Unavailable codes: degenerate response 288; singular
covariance 252; nonestimable 192; low permutation resolution 160;
quasibinomial paired-scope exclusion 128; low reference df 128; low independent
support 96; numerical failure 4. The four numerical failures were two
hypotheses in quasibinomial `dda_batch_crossed` replicates 15 and 17; their runs
remained measured and fail-closed. Structural-off, perfect-confounding, low-
support, paired-scope, grouped-leakage, and permutation-resolution controls
therefore remained explicit rather than improving successful-fit denominators.

Runtime medians per run: OLS/HC3/CR2 0.472 s; Freedman-Lane 4.060 s;
quasibinomial 0.323 s. Frozen elapsed values are descriptive only and do not
enter selection.

## Reproduction + disposition

From the repository root with the project R library:

```sh
Rscript --vanilla dev/m13-association-study.R --status
Rscript --vanilla dev/m13-association-study.R --verify-evidence
```

Independent complete replay reproduced the evidence hash and returned
`verified_no_winner`; final status reported manifest present, `416/416` inputs,
`1248/1248` results, and evidence present. Under `PLAN.md` portfolio semantics,
frozen evidence rejected the current statistical association direction, so the
panel is `killed`, not materialized or advanced to development. This negative
result does not weaken or remove the mandatory design/estimability core.
