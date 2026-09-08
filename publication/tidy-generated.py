#!/usr/bin/env python3
"""Canonicalise generated R text only; preserve every upstream byte and numerical artifact."""
from pathlib import Path

root = Path(__file__).resolve().parent
files = ("session-info.txt", "downstream-session-info.txt", "real-diagnostics.R")
changed = 0
for name in files:
    path = root / "results" / name
    original = path.read_text()
    canonical = "\n".join(line.rstrip(" \t") for line in original.splitlines()) + "\n"
    if original != canonical:
        path.write_text(canonical)
        changed += 1
print(f"Generated diagnostic text canonicalised: {changed}/{len(files)} files changed.")
