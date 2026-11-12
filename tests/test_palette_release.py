"""End-to-end consistency of the downloadable scientific release."""
import hashlib
import json
from pathlib import Path
import unittest

from openpyxl import load_workbook

from analysis.palette_audit import audit
from analysis.palette_evaluation import records


DATA = Path(__file__).resolve().parents[1] / "docs" / "data"


def read(name):
    return json.loads((DATA / name).read_text(encoding="utf8"))


class PublishedPaletteTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.catalog = read("palette.json")
        cls.raw = read("palette-raw.json")
        cls.manifest = read("palette-manifest.json")

    def test_per_repeat_values_and_addresses_match_the_two_workbooks(self):
        for system, source in self.catalog["sources"].items():
            path = DATA / Path(source["local"]).name
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), source["sha256"])
            book = load_workbook(path, read_only=False, data_only=True)
            try:
                sheet = book[source["sheet"]]
                for dye in self.catalog["dyes"]:
                    panel = dye["systems"][system]
                    self.assertEqual(sheet[panel["column"] + "1"].value, panel["worksheetDye"])
                    for key, series in panel["metrics"].items():
                        scale = 100 if key in ("dropPct","destroyedPct") else 1
                        for point in series["points"]:
                            self.assertAlmostEqual(sheet[point["cell"]].value * scale, point["value"])
            finally:
                book.close()

    def test_portable_audit_and_evaluation_match_the_released_data(self):
        report = audit(self.catalog, self.raw, self.manifest)
        self.assertEqual(report, read("palette-audit.json"))
        self.assertEqual(records(self.catalog, self.raw), read("palette-evaluation.json"))
        self.assertEqual(report["primaryFieldDigestDuplicates"], [])
        self.assertTrue(report["cy3bFourthMainFieldAlsoOccursAs"])
        self.assertEqual(report["principalCases"], 108)
        self.assertEqual(report["mainFigureCases"], 109)
        self.assertTrue(self.raw["coverage"]["complete"])
        self.assertEqual(self.raw["coverage"]["totalFiles"], report["mainFigureHdf5"])

    def test_raw_role_identifiers_match_the_verified_manifest(self):
        listing = {x["path"]:x for x in self.manifest["files"]}
        for record in self.raw["replicates"]:
            for field in record["files"]:
                member=listing[field["path"]]
                self.assertEqual(field["archiveSha256"], member["sha256"])
                self.assertEqual(field["archiveCrc32"], member["crc32"])
                sidecar = Path(field["path"]).with_suffix(".yaml").as_posix()
                if sidecar in listing:
                    self.assertEqual(field["metadataPath"], sidecar)
                    self.assertEqual(field["metadataSha256"], listing[sidecar]["sha256"])
                else:
                    self.assertNotIn("metadataPath", field)
                self.assertEqual(field["rows"], member["datasets"]["locs"]["shape"][0])
