#!/usr/bin/env python3
"""Integrity-pinned full protein-group extraction; only ZIP member is transferred."""
import contextlib
import csv
import hashlib
import gzip
import importlib.util
import json
import math
from pathlib import Path
import sys
import zipfile

ROOT = Path(__file__).resolve().parent.parent
sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location(
    "fixture", ROOT / "data-raw/build-pxd020717-example.py"
)
fixture = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fixture)
MEMBER = "txt/proteinGroups.txt"
SHA = "d7defecb4f3a0acf2e031844e2f17521df7aa57a8d89d966c4162e2c550a1a7f"
cache = ROOT / "publication/cache/proteinGroups.txt"
cache.parent.mkdir(exist_ok=True)
if not cache.exists():
    snapshot = ROOT / "publication/data/pxd020717-proteinGroups.tsv.gz"
    if snapshot.exists():
        cache.write_bytes(gzip.decompress(snapshot.read_bytes()))
    else:
        with contextlib.closing(fixture.archive_reader(None)) as reader:
            with zipfile.ZipFile(reader) as archive:
                info = archive.getinfo(MEMBER)
                if (info.file_size, info.compress_size, info.CRC) != (
                    9642977, 2956620, 0x116D7C34
                ):
                    raise ValueError("PXD020717 protein-table member metadata changed")
                cache.write_bytes(archive.read(info))
if hashlib.sha256(cache.read_bytes()).hexdigest() != SHA:
    raise ValueError("PXD020717 protein-table SHA256 mismatch")
with cache.open(encoding="utf-8-sig", newline="") as stream:
    reader = csv.DictReader(stream, delimiter="\t")
    columns = [x for x in reader.fieldnames if x.startswith("Intensity ")]
    rows = list(reader)
flags = ("Only identified by site", "Reverse", "Potential contaminant")
clean = [r for r in rows if all(r[f] != "+" for f in flags)]
samples = [c.removeprefix("Intensity ") for c in columns]
donors = sorted({s.split("_")[0] for s in samples})
out = ROOT / "publication/data"
out.mkdir(exist_ok=True)
(out / "pxd020717-proteinGroups.tsv.gz").write_bytes(gzip.compress(cache.read_bytes(), mtime=0))

def emit(name, fields, records):
    with (out / name).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(records)

def log2(x):
    return "" if not x else format(math.log2(x), ".10f")

def positive_mean(values):
    positive = [x for x in values if x > 0]
    return sum(positive) / len(positive) if positive else 0

emit("pxd020717-runs.csv", ["feature", *samples], (
    ["PG" + r["id"], *(log2(float(r[c] or 0)) for c in columns)] for r in clean
))
emit("pxd020717-donors.csv", ["feature", *donors], (
    ["PG" + r["id"], *(log2(positive_mean([
        float(r[c] or 0) for c, s in zip(columns, samples)
        if s.split("_")[0] == donor
    ])) for donor in donors)] for r in clean
))
emit("pxd020717-features.csv", ["feature", "protein_ids", "gene_names"], (
    ["PG" + r["id"], r["Protein IDs"], r["Gene names"]] for r in clean
))
emit("pxd020717-samples.csv", ["sample", "donor", "condition"], (
    [s, s.split("_")[0], "AD" if s.startswith("AD") else "control"] for s in samples
))
manifest = dict(url=fixture.ARCHIVE_URL, archive_bytes=fixture.ARCHIVE_SIZE,
    member=MEMBER, member_sha256=SHA, input_rows=len(rows),
    retained_rows=len(clean), excluded_flagged_rows=len(rows)-len(clean),
    runs=len(samples), donors=len(donors),
    transform="Mean positive linear intensities within donor; zero-only = missing; log2.",
    files={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
           for p in sorted(out.glob("pxd020717-*.csv"))})
(out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
print(json.dumps({k: manifest[k] for k in ("input_rows", "retained_rows", "runs", "donors")}))
