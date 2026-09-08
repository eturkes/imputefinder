#!/usr/bin/env python3
"""Independent negative controls for the inventory guarantee, without running analysis."""
import importlib.util
from pathlib import Path, PureWindowsPath
import tempfile
import unittest

spec = importlib.util.spec_from_file_location("reproduction", Path(__file__).with_name("reproduce.py"))
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class InventoryContract(unittest.TestCase):
    def setUp(self):
        cache = Path(__file__).parent / "cache"
        cache.mkdir(exist_ok=True)
        self.temporary = tempfile.TemporaryDirectory(dir=cache)
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        for name in ("data", "results", "figures"):
            (self.root / "publication" / name).mkdir(parents=True)
        self.current = self.root / "publication/results/current.csv"
        self.current.write_text("value\n1\n")
        self.expected = {"publication/results/current.csv"}

    def test_valid_inventory_passes(self):
        module.check_inventory(module.evidence_inventory(self.root), self.expected)

    def test_orphan_is_rejected(self):
        self.current.with_name("stale.csv").write_text("obsolete\n1\n")
        with self.assertRaisesRegex(ValueError, "orphaned=.*stale"):
            module.check_inventory(module.evidence_inventory(self.root), self.expected)

    def test_missing_output_is_rejected(self):
        self.current.unlink()
        with self.assertRaisesRegex(ValueError, "missing=.*current"):
            module.check_inventory(module.evidence_inventory(self.root), self.expected)

    def test_windows_ids_use_manifest_separators(self):
        self.assertEqual(module.relative_id(
            PureWindowsPath("C:/work/publication/results/current.csv"),
            PureWindowsPath("C:/work")), "publication/results/current.csv")


if __name__ == "__main__":
    unittest.main(verbosity=1)
