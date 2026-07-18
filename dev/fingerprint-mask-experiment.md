# M11c fingerprint + original-mask decision

Decision date: 2026-07-18. Scope: experimental
`imputefinder_analysis_v1`; stable `rules_v1` remains untouched.

## Decision

- Fingerprint = SHA-256 via `tools::sha256sum(bytes = ...)` over the
  versioned `matrix_core_be_v1` byte protocol.
- Protocol = domain tag + storage tag + 32-bit big-endian dimensions +
  length-prefixed names + packed missing mask + big-endian observed values,
  all in matrix column-major order. Text names canonicalize to UTF-8; names
  explicitly marked `bytes` retain byte identity. Missing numeric payloads
  canonicalize to zero because their locations are already encoded by the
  mask. Integer and double storage remain distinct.
- Original mask = inline raw bitset from `packBits(..., type = "raw")`, LSB
  first, column-major, zero-padded to a complete byte. The record stores its
  encoding and cell count; dimensions and ordered names live beside it in the
  input record.
- Fingerprints identify the validated numeric matrix core, including order,
  names, storage, values, and missingness. The matrix versus
  `SummarizedExperiment` wrapper, assay routing, typed design, and declarations
  remain separately validated provenance.

The fingerprint is a collision-resistant mismatch sentinel, not a signature,
proof of historical provenance, or authorization primitive. Unknown algorithm,
canonicalization, mask encoding, analysis schema, or lifecycle identifiers fail
with typed reconstruction guidance.

## Why this pair

`tools::sha256sum()` accepts raw bytes at the package's R floor, returns a
64-hex-character digest, adds no dependency or system library, and uses the
NIST-approved SHA-256 primitive. SHA-256 provides 128-bit collision strength
and 256-bit preimage strength. OpenSSL was materially faster in the local
microbenchmark, but hashing runs once per sidecar construction/recheck and its
speedup does not justify another runtime/system dependency. `digest` likewise
adds a dependency without improving the protocol. MD5 is absent from the
candidate set because collision resistance is part of the requirement.

The byte protocol avoids hashing arbitrary R serialization. R documents XDR as
portable, but also warns that its serialization format may change. Explicit
big-endian numeric bytes, tagged storage, length-delimited identifiers, and a
domain/version prefix make the compatibility boundary ours. The protocol is
injective over accepted matrix-core values before SHA-256's collision bound.

A logical matrix costs four bytes per cell in current R; one raw byte per cell
still wastes seven bits. Sparse integer indices win near 1% missingness but
grow with missingness and lose at realistic moderate/high proteomics rates.
The packed bitset is density-invariant, 32x smaller than logical storage, and
uses documented base primitives with single-digit-millisecond pack/unpack
medians at the 500,000-cell release scale.

## Frozen experiment

Harness: `dev/fingerprint-mask-experiment.R`, SHA-256
`0857e92b191457779a71f67545964533d106d93b27e9b455a85f64f0e885dd60`.
Debian x86-64, R 4.6.1; 10,000 x 50 double matrix, 20% missing for hashing;
500,000-cell masks. SHA-256 candidates agreed exactly. Representative first-run
medians (machine load makes timing descriptive, not a gate):

| SHA-256 implementation | Median ms / 4.05 MiB payload | Runtime surface |
|---|---:|---|
| `tools` | 43.10 | base R |
| `digest` | 81.85 | added R package |
| `openssl` | 4.85 | added R package + OpenSSL |

Serialized representation sizes and base pack/unpack timing:

| Missing | Logical | Raw/cell | Packed | Sparse indices | Pack ms | Unpack ms |
|---:|---:|---:|---:|---:|---:|---:|
| 1% | 2,000,031 | 500,031 | 62,531 | 20,043 | 1.40 | 5.12 |
| 20% | 2,000,031 | 500,031 | 62,531 | 400,031 | 1.16 | 4.98 |
| 50% | 2,000,031 | 500,031 | 62,531 | 1,000,047 | 1.40 | 5.62 |

Bytes include each representation's R serialization envelope solely to compare
sidecar storage. The production fingerprint hashes the explicit byte protocol,
not those envelopes. Exact canonical fixture:
`a09c1f5501caee9b7d9a6dcce7322e56c2028a19c103c1f3bba7764d6338a24a`.

## Primary sources

- R `tools::sha256sum()` raw-byte contract:
  <https://stat.ethz.ch/R-manual/R-patched/library/tools/html/sha256sum.html>
- NIST approved hashes + security strengths:
  <https://csrc.nist.gov/Projects/hash-functions>
- R `packBits()` / `rawToBits()` LSB-first contract:
  <https://stat.ethz.ch/R-manual/R-devel/library/base/help/rawConversion.html>
- R binary transfer + explicit endianness:
  <https://stat.ethz.ch/R-manual/R-devel/library/base/html/readBin.html>
- R serialization XDR behavior + long-term warning:
  <https://stat.ethz.ch/R-manual/R-devel/library/base/help/serialize.html>
- R character-encoding normalization:
  <https://stat.ethz.ch/R-manual/R-devel/library/base/html/Encoding.html>
