import json
import math
import unittest
from hashlib import sha256
from statistics import mean, median

from openpyxl import load_workbook

from analysis.probe_design import DATA, DYE_NAMES, HASHES, build, file_name


class ProbeEvidence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = build()
        cls.published = json.loads((DATA / "probe-design.json").read_text(encoding="utf-8"))
        cls.by_name = {r["id"]: r for r in cls.data["dyes"]}

    def test_bundle_and_structures(self):
        self.assertEqual(list(self.by_name), DYE_NAMES)
        for n, digest in HASHES.items():
            self.assertEqual(sha256((DATA / file_name(n)).read_bytes()).hexdigest(), digest)
        for d in self.data["dyes"]:
            image = DATA.parent / d["structure"]["image"]
            self.assertTrue(image.read_bytes().startswith(b"\x89PNG\r\n\x1a\n"))
        self.assertEqual(self.data, self.published)

    def test_titration_blanks_zero_and_coordinates(self):
        self.assertEqual([len(d["titration"]) for d in self.data["dyes"]], [25]*4 + [31]*5)
        self.assertNotIn(1.14, [p["pH"] for p in self.by_name["JF639b"]["titration"]])
        self.assertEqual(self.by_name["JFX650b"]["titration"][-1]["readings"], [0, 0])
        book = load_workbook(DATA / file_name(9), data_only=True)
        try:
            sheet = book[self.data["sources"]["9"]["sheet"]]
            for d in self.data["dyes"]:
                for p in d["titration"]:
                    self.assertEqual([c.value for c in sheet[p["cells"]][0]], p["readings"])
                    self.assertEqual(sheet[p["pHCell"]].value, p["pH"])
                    self.assertAlmostEqual(mean(p["readings"]), p["mean"])
        finally:
            book.close()

    def test_context_identity_and_experimental_units(self):
        er, lys = (self.data["contexts"][k] for k in ("er", "lysosome"))
        self.assertEqual((er["nCells"], lys["nCells"]), (4, 5))
        self.assertEqual(set(lys["dyes"]), {"JF639b", "JF614b"})
        self.assertNotIn("useful", lys["metrics"])
        for state, expected in (("live", [41.1591, 29.0161]), ("fixed", [15.306, 52.1853])):
            self.assertEqual([er["dyes"][n][state]["useful"]["mean"] for n in ("JF639b", "JF630b")], expected)
        for context in [er, lys]:
            for d in context["dyes"].values():
                for state in d.values():
                    self.assertEqual(set(state), set(context["metrics"]))
                    self.assertTrue(all(r["sem"] >= 0 for r in state.values()))

    def test_same_ligand_holdout_for_all_three_methods(self):
        c = self.data["calibration"]
        pairs = {p["id"]: p for p in c["pairs"]}
        self.assertEqual(set(self.by_name) - set(pairs), {"HM-SiR", "JF646b", "JF614b"})
        self.assertEqual(c["fit"]["pKaRange"], [1.22, 4.48])
        self.assertEqual(len(c["folds"]), 6)
        for fold in c["folds"]:
            self.assertEqual(set(fold["training"]), set(pairs) - {fold["held"]})
            train = [pairs[d] for d in fold["training"]]
            held = pairs[fold["held"]]
            low, high = min(p["pKa"] for p in train), max(p["pKa"] for p in train)
            if not low <= held["pKa"] <= high:
                self.assertEqual(fold["status"], "outside training range")
                self.assertNotIn("error", fold)
                continue
            x = [p["pKa"] for p in train]
            y = [math.log10(p["ratio"]) for p in train]
            beta = sum((a - mean(x)) * (b - mean(y)) for a,b in zip(x,y)) / sum((a - mean(x))**2 for a in x)
            predicted = mean(y) + beta * (held["pKa"] - mean(x))
            self.assertAlmostEqual(math.log10(fold["predicted"]), predicted)
            self.assertAlmostEqual(math.log10(fold["baseline"]), median(y))
            ref = min(train, key=lambda p: abs(p["pKa"] - held["pKa"]))
            self.assertEqual((fold["nearestDye"], fold["nearestRatio"]), (ref["id"], ref["ratio"]))
        self.assertEqual((c["validation"]["scored"], c["validation"]["winsVsNearest"]), (4, 1))
        for name, expected in (("mae", .31771044776), ("nearestMae", .21620068859), ("baselineMae", .49878348336)):
            self.assertAlmostEqual(c["validation"][name], expected, places=10)


if __name__ == "__main__":
    unittest.main()
