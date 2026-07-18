# Sidecar differential validation

Status → M11f premeasurement contract. The fixture, comparison categories, and
thresholds below were fixed before the first execution. The harness first runs
the unchanged M10 stable performance gate, then measures the compatible
`base_fit` sidecar path on the same deterministic workload.

Run after installing the current source tree:

```sh
R_LIBS_USER="$PWD/.agent/R-library" \
  Rscript --vanilla dev/sidecar-differential-validation.R
```

## Differential fixture

The compact routine fixture and representative-scale harness require:

- exact → compatible supplied result retained byte-for-byte; discrete-state
  mutation reported only in the exact category;
- canonical → stored profile-warning mutation reported only in the canonical
  category;
- tolerance → a one-quarter absolute-tolerance perturbation accepted and a
  100-fold perturbation reported only in the tolerance category;
- performance → original input unchanged, packed mask exactly one bit/cell,
  and bounded sidecar runtime, allocation, and retained-object overhead.

## Representative workload and frozen gates

The workload is the M10 10,000-feature x 50-sample deterministic LFQ matrix.
The sidecar uses a compatible manual `base_fit`, default modules, and seed 101.
It therefore exercises input fingerprinting, mask packing, pre-rescue evidence,
full stable recomputation, exact/canonical/tolerance comparison, and structured
stability abstention.

On the canonical Debian/R/Bioconductor development environment, the sidecar
must:

- complete in median <=30 seconds and <=8x the same-process stable median;
- allocate <=250x input size cumulatively and <=4x the stable profiled call;
- make no single allocation larger than 2x input size;
- return an object <=8x input size and <=1.5x its embedded classic result;
- preserve the supplied result exactly, preserve input exactly, pass all four
  compact comparator checks, and store exactly `ceiling(cells / 8)` mask bytes.

The wide timing and allocation ratios are regression tripwires, not performance
claims: they tolerate instrumentation and machine noise while rejecting an
accidental second classic recomputation, duplicated full evidence trees, or
quadratic copying. The object-size gates directly bound the serialized sidecar
surface that runtime ratios cannot protect.

## Assessment

Harness SHA-256 =
`63c4a943efccfe81e1c6445f363c6db0941ae99b087a6ab7ec91f7684332b471`.
The script remained unchanged after the thresholds were frozen.

Reference run → canonical Debian x86_64 development container; R 4.6.1;
imputefinder 0.99.4. The unchanged M10 gate passed first: stable manual and
automatic medians were .146/.168 seconds, allocation ratios 45.37x/66.53x,
result ratios 4.70x/4.72x, peak RSS 176.57 MiB, and maximum automatic-cutoff
error .282.

| Sidecar measure | Result | Gate |
|---|---:|---:|
| Median elapsed (range) | .426 s (.414-.461) | <=30 s |
| Elapsed/stable manual | 2.92x | <=8x |
| Cumulative allocation | 465.55 MiB (103.34x input) | <=250x input |
| Allocation/stable manual | 2.28x | <=4x |
| Largest allocation | 4.05 MiB (.90x input) | <=2x input |
| Sidecar size | 24.68 MiB (5.48x input) | <=8x input |
| Sidecar/embedded classic size | 1.16x | <=1.5x |

All exact-input/result/mask, exact/canonical/tolerance comparator, runtime,
allocation, and object-size gates passed. The sidecar adds one fingerprint/mask
pass, pre-rescue summaries, and comparison work without a second classic
recomputation or a duplicated original numeric matrix.
