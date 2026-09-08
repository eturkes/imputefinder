#!/usr/bin/env python3
"""Extract the complete iPRG2015 raw-intensity matrix and declared spike truth."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.util
import io
import json
import math
import sys
import urllib.request
import zipfile
from pathlib import Path


URL = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/04/PXD028038/iPRG2015_MQ_txt.zip.zip"
ARCHIVE_BYTES = 580636238
MEMBER = "txt/proteinGroups.txt"
MEMBER_METADATA = (7658671, 2449219, 0xDDA01138)
SOURCE_SHA256 = "4ab473a7a1bd942ac8d95e29e5993718f581a0bec06f250a7c7f604c79515792"
AUTHOR_COMMIT = "962d154ffe20b4eea27a341f56f3d67f8a6cd580"
SPIKES = {
    "P44015": ("Ovalbumin", (65, 55, 15, 2)),
    "P55752": ("Myoglobin", (55, 15, 2, 65)),
    "P44374": ("Phosphorylase_b", (15, 2, 65, 55)),
    "P44983": ("BetaGalactosidase", (2, 65, 55, 15)),
    "P44683": ("BSA", (11, 0.6, 10, 500)),
    "P55249": ("CarbonicAnhydrase", (10, 500, 11, 0.6)),
}
SAMPLES = tuple(f"{group}_{rep}" for group in range(1, 5) for rep in "ABC")


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def fetch(package_root: Path) -> bytes:
    # Reuse the package's seekable range reader; keep the source tree read-only.
    sys.dont_write_bytecode = True
    source = package_root / "data-raw/build-pxd020717-example.py"
    spec = importlib.util.spec_from_file_location("range_reader", source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    with urllib.request.urlopen(urllib.request.Request(URL, method="HEAD"), timeout=60) as response:
        size = int(response.headers["Content-Length"])
    if size != ARCHIVE_BYTES:
        raise ValueError(f"Archive size changed: {size}")
    with module.HTTPRangeReader(URL, size) as reader, zipfile.ZipFile(reader) as archive:
        info = archive.getinfo(MEMBER)
        if (info.file_size, info.compress_size, info.CRC) != MEMBER_METADATA:
            raise ValueError("Archive member metadata changed")
        return archive.read(info)


def write_csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def accessions(text: str) -> set[str]:
    return {part.split("|")[1] if "|" in part else part for part in text.split(";")}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("publication/data"))
    parser.add_argument("--source", type=Path, help="Pinned proteinGroups TSV or gzip for offline regeneration.")
    parser.add_argument("--package-root", type=Path, default=Path.cwd())
    args = parser.parse_args()
    if args.source:
        data = args.source.read_bytes()
        if args.source.suffix == ".gz":
            data = gzip.decompress(data)
    else:
        data = fetch(args.package_root)
    if digest(data) != SOURCE_SHA256:
        raise ValueError("Protein-table SHA256 mismatch")
    reader = csv.DictReader(io.StringIO(data.decode("utf-8-sig")), delimiter="\t")
    rows = list(reader)
    columns = [f"Intensity {sample}" for sample in SAMPLES]
    if len(rows) != 3314 or not set(columns).issubset(reader.fieldnames):
        raise ValueError("Protein-table shape changed")
    matrices, features, exclusions, truth = [], [], [], []
    found = set()
    for index, row in enumerate(rows, 1):
        feature = f"iprg{index:04d}"
        spike_ids = accessions(row["Protein IDs"]) & SPIKES.keys()
        if len(spike_ids) > 1:
            raise ValueError("Mixed spike groups require separate interpretation")
        spike = next(iter(spike_ids), "")
        reasons = [field for field in ("Reverse", "Only identified by site", "Potential contaminant")
                   if row[field].strip() == "+" and not (spike and field == "Potential contaminant")]
        if reasons:
            exclusions.append({"feature": feature, "protein_ids": row["Protein IDs"], "reason": ";".join(reasons)})
            continue
        values = [float(row[col] or 0) for col in columns]
        if any(not math.isfinite(x) or x < 0 for x in values):
            raise ValueError("Invalid raw intensity")
        matrices.append({"feature": feature, **{
            sample: "NA" if value == 0 else format(math.log2(value), ".17g")
            for sample, value in zip(SAMPLES, values)}})
        features.append({"feature": feature, "protein_ids": row["Protein IDs"],
                         "majority_protein_ids": row["Majority protein IDs"],
                         "spike_id": spike, "spike_name": SPIKES[spike][0] if spike else "",
                         "truth_class": "spike" if spike else "yeast_background",
                         "observed": sum(x > 0 for x in values)})
        if spike:
            if spike in found:
                raise ValueError("Duplicate spike protein group")
            found.add(spike)
            for condition, amount in enumerate(SPIKES[spike][1], 1):
                truth.append({"feature": feature, "spike_id": spike, "protein": SPIKES[spike][0],
                              "condition": str(condition), "amount_fmol": amount})
    if found != SPIKES.keys():
        raise ValueError(f"Missing spike identities: {SPIKES.keys() - found}")
    args.output.mkdir(parents=True, exist_ok=True)
    write_csv(args.output / "iprg2015-intensities.csv", ["feature", *SAMPLES], matrices)
    write_csv(args.output / "iprg2015-features.csv", list(features[0]), features)
    write_csv(args.output / "iprg2015-exclusions.csv", ["feature", "protein_ids", "reason"], exclusions)
    write_csv(args.output / "iprg2015-truth.csv", list(truth[0]), truth)
    write_csv(args.output / "iprg2015-samples.csv", ["sample", "condition", "technical_replicate"],
              [{"sample": s, "condition": s[0], "technical_replicate": s[-1]} for s in SAMPLES])
    # An immutable local source makes offline extraction independent of the remote ZIP.
    source_path = args.output / "iprg2015-proteinGroups.tsv.gz"
    source_path.write_bytes(gzip.compress(data, mtime=0))
    manifest = {
        "original_study": "iPRG2015; four mixtures from one yeast lysate, three technical runs per mixture",
        "original_accession": "MSV000079843",
        "reprocessing_accession": "PXD028038",
        "reprocessing_article": "https://doi.org/10.1021/acs.jproteome.1c00683",
        "source_url": URL, "archive_bytes": ARCHIVE_BYTES, "member": MEMBER,
        "member_bytes": MEMBER_METADATA[0], "member_compressed_bytes": MEMBER_METADATA[1],
        "member_crc32": f"{MEMBER_METADATA[2]:08x}", "source_sha256": SOURCE_SHA256,
        "spike_amounts_source": f"https://raw.githubusercontent.com/MassDynamics/lfq_benchmark/{AUTHOR_COMMIT}/resources/IPRGSpikes.csv",
        "spike_id_source": f"https://raw.githubusercontent.com/MassDynamics/lfq_benchmark/{AUTHOR_COMMIT}/app/benchmarkers/iprg_benchmarker.py",
        "identifiers": "Blinded surrogate spike accessions supplied by the study; not current organism annotations",
        "intensity": "MaxQuant 1.6.17.0 Intensity columns, not LFQ intensity; zero becomes NA; positive values become log2",
        "upstream_matching": "Match between runs enabled; this analysis does not repeat peptide identification or alignment",
        "filtering": "Exclude reverse, site-only, and contaminant groups; retain declared spikes even if contaminant-flagged",
        "source_rows": len(rows), "retained_rows": len(matrices), "excluded_rows": len(exclusions),
        "samples": len(SAMPLES), "spike_groups": len(found),
        "truth_limit": "Known mixture effects and constant yeast background; no known causes for individual missing cells",
        "outputs_sha256": {p.name: digest(p.read_bytes()) for p in sorted(args.output.glob("iprg2015-*.csv"))},
    }
    (args.output / "iprg2015-manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"iPRG2015: {len(matrices)} retained groups, {len(exclusions)} excluded groups, 12 runs, six spikes")


if __name__ == "__main__":
    main()
