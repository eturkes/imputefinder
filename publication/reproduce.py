#!/usr/bin/env python3
"""Regenerate evidence from empty output directories; --verify checks exact inventory and hashes."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent.parent
DATA = [f"pxd020717-{name}.csv" for name in ("donors", "runs", "features", "samples")] + [
    f"iprg2015-{name}.csv" for name in ("intensities", "features", "exclusions", "truth", "samples")]
RESULTS = [f"{name}.csv" for name in (
    "simulation", "simulation-summary", "downstream", "downstream-summary", "downstream-paired-vs-DEP",
    "real-summary", "real-bootstrap", "real-quality", "real-run-sensitivity", "real-rescue-seeds",
    "holdout", "holdout-summary", "spikein", "spikein-effects", "spikein-normalization", "proda-solver-probe")]
RESULTS += [f"real-states-{rescue}-{delta}.csv" for rescue in ("none", "minimum") for delta in ("auto", "-1", "0", "1")]
FIGURES = [f"{name}.{ext}" for name in ("simulation", "downstream", "real-cutoffs") for ext in ("svg", "png")]
FIGURES += [f"profile-{rescue}-{condition}.svg" for rescue in ("none", "minimum") for condition in ("AD", "control")]
EXPECTED = {f"publication/{folder}/{name}" for folder, names in (("data", DATA), ("results", RESULTS), ("figures", FIGURES)) for name in names}


def relative_id(path, root):
    return path.relative_to(root).as_posix()


def evidence_inventory(root):
    pub = root / "publication"
    paths = [*list((pub / "data").glob("*.csv")), *list((pub / "results").glob("*.csv")),
             *list((pub / "figures").glob("*.svg")), *list((pub / "figures").glob("*.png"))]
    return {relative_id(p, root) for p in paths}


def check_inventory(observed, expected):
    if observed != expected:
        raise ValueError(f"Evidence inventory mismatch; missing={sorted(expected-observed)}, orphaned={sorted(observed-expected)}")


def evidence_hashes(root):
    check_inventory(evidence_inventory(root), EXPECTED)
    return {name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in sorted(EXPECTED)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    os.chdir(ROOT)
    pub = ROOT / "publication"
    manifest_path = pub / "artifact-manifest.json"
    expected = json.loads(manifest_path.read_text()) if args.verify else None
    if expected is not None:
        check_inventory(set(expected), EXPECTED)
    environment = dict(os.environ, TZ="UTC", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    if environment.get("R_LIBS_USER"):
        environment["R_LIBS_USER"] = os.pathsep.join(
            str(Path(p).expanduser().resolve()) for p in environment["R_LIBS_USER"].split(os.pathsep)
        )
    cache = pub / "cache"
    cache.mkdir(exist_ok=True)
    # Both destinations contain script-owned generated output only.
    # Empty destinations prevent orphans from surviving removed generation steps.
    for folder in ("results", "figures"):
        destination = pub / folder
        if destination.exists():
            shutil.rmtree(destination)
        destination.mkdir()
    commands = [
        [sys.executable, "publication/test-reproduction.py"],
        [sys.executable, "publication/fetch-data.py"],
        [sys.executable, "publication/fetch-spikein.py", "--source", "publication/data/iprg2015-proteinGroups.tsv.gz"],
        *[["Rscript", f"publication/{name}.R"] for name in (
            "benchmark", "proda-diagnostics", "downstream", "real-data", "holdout", "spikein")],
        [sys.executable, "publication/tidy-generated.py"],
        *[["Rscript", f"publication/{name}.R"] for name in ("check-results", "figures", "render")],
    ]
    for command in commands:
        name = Path(command[1]).stem
        print(f"Running {name}", flush=True)
        with (cache / f"reproduce-{name}.log").open("w") as log:
            subprocess.run(command, env=environment, check=True, stdout=log, stderr=subprocess.STDOUT)
    observed = evidence_hashes(ROOT)
    if expected is not None:
        changed = sorted(k for k in EXPECTED if expected[k] != observed[k])
        if changed:
            raise SystemExit("Reproduction mismatch: " + ", ".join(changed))
        print(f"PASS: {len(observed)} evidence artifacts byte-identical after fresh analyses.")
    else:
        manifest_path.write_text(json.dumps(observed, indent=2, sort_keys=True) + "\n")
        print(f"Wrote {len(observed)} evidence hashes; run --verify for independent replay.")


if __name__ == "__main__":
    main()
