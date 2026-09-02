#!/usr/bin/env python3
"""Build the compact PXD020717 Tau-phosphosite vignette fixture."""

from __future__ import annotations

import argparse
import contextlib
import csv
import hashlib
import io
import math
import urllib.request
import zipfile
from pathlib import Path
from typing import BinaryIO


ARCHIVE_URL = (
    "https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/12/PXD020717/"
    "MaxQuant_BA46_SOLUBLE.zip"
)
ARCHIVE_SIZE = 933_534_392
MEMBER = "txt/Phospho (STY)Sites.txt"
MEMBER_SIZE = 870_975
MEMBER_COMPRESSED_SIZE = 164_619
MEMBER_CRC = 0xB6967597
SAMPLES = (
    "CTR03_SOLUBLE_01",
    "CTR06_SOLUBLE_01",
    "CTR45_SOLUBLE_01",
    "AD03_SOLUBLE_01",
    "AD50_SOLUBLE_01",
    "AD51_SOLUBLE_01",
)
SITES = (
    ("MAPT_pS404", "S", "404"),
    ("MAPT_pS199", "S", "199"),
    ("MAPT_pS202", "S", "202"),
    ("MAPT_pT403", "T", "403"),
    ("MAPT_pT217", "T", "217"),
    ("MAPT_pT181", "T", "181"),
    ("MAPT_pT231", "T", "231"),
)
DEFAULT_OUTPUT = Path("inst/extdata/PXD020717_tau_phosphosites.csv")
EXPECTED_OUTPUT_SHA256 = "eb65d2e543da1ea97010f76f7f2d3c0e234001a312adb429d8c71113b6c55e1d"


class HTTPRangeReader(io.RawIOBase):
    """Seekable reader that transfers only byte ranges requested by zipfile."""

    def __init__(self, url: str, size: int, timeout: float = 120.0) -> None:
        self.url = url
        self.size = size
        self.timeout = timeout
        self.position = 0

    def readable(self) -> bool:
        return True

    def seekable(self) -> bool:
        return True

    def tell(self) -> int:
        return self.position

    def seek(self, offset: int, whence: int = io.SEEK_SET) -> int:
        if whence == io.SEEK_SET:
            position = offset
        elif whence == io.SEEK_CUR:
            position = self.position + offset
        elif whence == io.SEEK_END:
            position = self.size + offset
        else:
            raise ValueError(f"Unsupported seek mode: {whence}")
        if position < 0:
            raise ValueError("Cannot seek before the archive start")
        self.position = position
        return position

    def read(self, size: int = -1) -> bytes:
        if self.position >= self.size:
            return b""
        if size < 0:
            size = self.size - self.position
        end = min(self.size - 1, self.position + size - 1)
        request = urllib.request.Request(
            self.url,
            headers={"Range": f"bytes={self.position}-{end}"},
        )
        with urllib.request.urlopen(request, timeout=self.timeout) as response:
            if response.status != 206:
                raise RuntimeError(
                    f"Expected HTTP 206 for a range request; got {response.status}"
                )
            data = response.read()
        expected = end - self.position + 1
        if len(data) != expected:
            raise RuntimeError(
                f"Range request returned {len(data)} bytes; expected {expected}"
            )
        self.position += len(data)
        return data

    def readinto(self, buffer: bytearray) -> int:
        data = self.read(len(buffer))
        buffer[: len(data)] = data
        return len(data)


def verify_remote_size(url: str) -> None:
    request = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(request, timeout=120.0) as response:
        observed = int(response.headers["Content-Length"])
    if observed != ARCHIVE_SIZE:
        raise RuntimeError(
            f"Remote archive size changed: observed {observed}, expected {ARCHIVE_SIZE}"
        )


def archive_reader(local_archive: Path | None) -> BinaryIO:
    if local_archive is None:
        verify_remote_size(ARCHIVE_URL)
        return HTTPRangeReader(ARCHIVE_URL, ARCHIVE_SIZE)
    if local_archive.stat().st_size != ARCHIVE_SIZE:
        raise RuntimeError(
            f"Local archive size changed: observed {local_archive.stat().st_size}, "
            f"expected {ARCHIVE_SIZE}"
        )
    return local_archive.open("rb")


def read_selected_sites(reader: BinaryIO) -> dict[str, dict[str, str]]:
    targets = {(amino_acid, position): label for label, amino_acid, position in SITES}
    selected: dict[str, dict[str, str]] = {}
    with zipfile.ZipFile(reader) as archive:
        info = archive.getinfo(MEMBER)
        observed = (info.file_size, info.compress_size, info.CRC)
        expected = (MEMBER_SIZE, MEMBER_COMPRESSED_SIZE, MEMBER_CRC)
        if observed != expected:
            raise RuntimeError(
                f"Source member metadata changed: observed {observed}, expected {expected}"
            )
        with archive.open(info) as binary:
            text = io.TextIOWrapper(binary, encoding="utf-8-sig", newline="")
            rows = csv.DictReader(text, delimiter="\t")
            for row in rows:
                if row["Gene names"] != "MAPT":
                    continue
                proteins = row["Proteins"].split(";")
                positions = row["Positions within proteins"].split(";")
                if "P10636-8" not in proteins:
                    continue
                position = positions[proteins.index("P10636-8")]
                label = targets.get((row["Amino acid"], position))
                if label is not None:
                    selected[label] = row
    expected_labels = {label for label, _, _ in SITES}
    if set(selected) != expected_labels:
        raise RuntimeError(
            f"Selected sites changed: observed {sorted(selected)}, "
            f"expected {sorted(expected_labels)}"
        )
    return selected


def transformed_intensity(row: dict[str, str], sample: str) -> str:
    value = float(row[f"Intensity {sample}"] or 0.0)
    return "" if value <= 0 else f"{math.log2(value):.3f}"


def write_fixture(
    output: Path,
    selected: dict[str, dict[str, str]],
) -> str:
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["sample", "condition", "subject", *(label for label, _, _ in SITES)]
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for sample in SAMPLES:
            writer.writerow(
                {
                    "sample": sample,
                    "condition": "control" if sample.startswith("CTR") else "AD",
                    "subject": sample.split("_", maxsplit=1)[0],
                    **{
                        label: transformed_intensity(selected[label], sample)
                        for label, _, _ in SITES
                    },
                }
            )
    return hashlib.sha256(output.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    with contextlib.closing(archive_reader(args.archive)) as reader:
        selected = read_selected_sites(reader)
    digest = write_fixture(args.output, selected)
    if digest != EXPECTED_OUTPUT_SHA256:
        raise RuntimeError(
            f"Generated fixture changed: observed {digest}, "
            f"expected {EXPECTED_OUTPUT_SHA256}"
        )
    print(f"Wrote {args.output} (sha256={digest})")


if __name__ == "__main__":
    main()
