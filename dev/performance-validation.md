# Release-candidate performance validation

Status → M10c premeasurement contract. Fixture, metrics, and gates below were
fixed before the first benchmark execution. The harness is deterministic and
uses only base R plus the installed package.

Run after installing the current source tree:

```sh
R_LIBS_USER="$PWD/.agent/R-library" \
  Rscript --vanilla dev/performance-validation.R
```

## Representative workload

- matrix = 10,000 protein features x 50 samples = five conditions x ten
  samples;
- deterministic condition-specific cliffs at log2 intensities 11-13, 5%
  intensity-independent background missingness, and 40 structural off-features
  per condition;
- manual and automatic public paths each run three timed repetitions plus one
  allocation-profiled repetition;
- every run verifies unchanged input and proves output cell differences equal
  the retained seed-log coordinates and values.

## Predeclared gates

On the canonical Debian/R/Bioconductor development environment, each path must:

- complete in median <=10 seconds;
- allocate <=100 x input size cumulatively in one profiled call;
- make no single allocation larger than 2 x input size;
- return an object <=12 x input size;
- keep process peak resident memory <=768 MiB when `/proc/self/status` exposes
  the high-water mark;
- place every automatic condition cutoff within one log2-intensity unit of its
  known boundary.

These broad gates detect accidental quadratic behavior, repeated whole-matrix
copying, runaway retained evidence, and scientifically invalid automatic output
without treating machine-scale timing noise as a regression.

## Assessment

Harness SHA-256 =
`74c53ccaf0f0c400693ac76607c575f686012c5b3980ae4d758a5dc4893fc2c2`.
The script remained unchanged after the gates were frozen.

Reference run → Debian 13 x86_64; Intel Core Ultra 7 268V; R 4.6.1;
Bioconductor 3.24; imputefinder 0.99.2.

| Path | Median s (range) | Allocated MiB (x input) | Largest MiB (x input) | Result MiB (x input) | Retained | Seeds | Max cutoff error |
|---|---:|---:|---:|---:|---:|---:|---:|
| Manual | .139 (.135-.151) | 204.4 (45.37x) | 3.81 (.85x) | 21.19 (4.70x) | 9,635 | 200 | 0 |
| Automatic | .178 (.174-.182) | 299.7 (66.53x) | 3.81 (.85x) | 21.25 (4.72x) | 9,635 | 200 | .282 |

Process peak resident memory = 171.63 MiB at the gate snapshot and 172.32 MiB at
process exit under `/usr/bin/time -v`. All six gates pass. Runtime is sub-second,
no allocation exceeds the input matrix, cumulative allocation stays below 67x
input, and the retained audit object stays below 5x input. No avoidable
pathological copying or runtime behavior was identified at the
release-candidate scale.
