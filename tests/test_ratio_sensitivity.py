import csv
import gzip
from itertools import product
import json
from math import sqrt
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from analysis.analogue_benchmark import benchmark, predict
from analysis.dyes import DATA
from analysis.ratio_audit import PAIRINGS, read_audit
from analysis.ratio_sensitivity import build


class RatioSensitivity(unittest.TestCase):
    def test_source_cells_and_estimands(self):
        rows = read_audit()
        self.assertEqual(len(rows), 9)
        outside = [r for r in rows if abs(r["ratio_of_means_diff_pct"]) > .1]
        self.assertEqual([r["dye"] for r in outside], ["JF626b"])
        anomaly = outside[0]
        self.assertAlmostEqual(anomaly["ratio_of_means_diff_pct"], 3.654770739, places=8)
        self.assertGreater(anomaly["paired_min"], anomaly["ratio_of_means"])
        self.assertEqual(anomaly["on_cells_ms"], "Extended Data Figure 4b-c!B9:D9")
        self.assertEqual(anomaly["j_cells_pka_ratio_sd_n"], "Extended Data Figure 4j!B8:E8")
        for row in rows:
            self.assertEqual(row["reported_n"], 3)
            if row["j_sd"] is not None:
                self.assertEqual(row["j_sd"], row["reported_sd"])
                self.assertAlmostEqual(row["j_sd_div_sqrt_n"], row["reported_sd"]/sqrt(3))
            for pairing in PAIRINGS:
                exact = sum(row["on_"+on+"_ms"] / row["off_"+off+"_ms"] for on, off in zip("BCD", pairing))/3
                self.assertAlmostEqual(row["paired_"+pairing], exact, places=15)

    def test_every_scenario_and_independent_scalar_oracle(self):
        raw = [r for r in read_audit() if r["free_pka"] is not None]
        with TemporaryDirectory() as tmp:
            audit = build(DATA, Path(tmp))
            checked = audit["pairings"]
            self.assertEqual(checked["overall"]["rankCountsLowestFirst"], {"local < nearest < regression": 46656})
            self.assertEqual(checked["scoredDyes"], 4)
            self.assertEqual(checked["calOrder"][0]["scenarios"], 46656)
            self.assertEqual(audit["registered"]["boundaries"], audit["estimands"]["j_reported"]["calBoundaries"])
            published = json.loads((DATA/"ratio-sensitivity.json").read_text(encoding="utf-8"))
            self.assertEqual(audit, published)
            with gzip.open(Path(tmp)/"ratio-scenarios.csv.gz", "rt") as handle:
                cases = csv.DictReader(handle)
                for index, (mix, case) in enumerate(zip(product(PAIRINGS, repeat=6), cases, strict=True)):
                    self.assertEqual(tuple(case[r["dye"]+"_off"] for r in raw), mix)
                    if index not in (0, 1, 15, 575, 3622, 12345, 23004, 46655):
                        continue
                    ref = [dict(id=r["dye"], dye=r["dye"], pKa=r["free_pka"], ratio=r["paired_"+p])
                           for r, p in zip(raw, mix)]
                    result = benchmark("test", "", "holland", "development", ref)
                    self.assertEqual(result["summary"]["scored"], 4)
                    cal = predict(ref, 3.2)["predictions"]
                    for m, value in result["summary"]["mae"].items():
                        self.assertAlmostEqual(float(case["mae_"+m]), value, places=9)
                        self.assertAlmostEqual(float(case["cal_"+m]), cal[m], places=12)
                    for a, b in (("nearest", "local"), ("regression", "nearest")):
                        self.assertAlmostEqual(float(case[f"tie_{a}_{b}"]), sqrt(cal[a]*cal[b]), places=12)

    def test_external_denominator_and_missing_scores(self):
        evidence = json.loads((DATA/"ratio-evidence.json").read_text(encoding="utf-8"))
        self.assertEqual(evidence["countsByStudy"], {"Uno 2014": 14, "Chen 2024": 4})
        self.assertEqual(len({r["dye"] for r in evidence["records"]}), evidence["distinctScreened"])
        self.assertEqual((evidence["matchedScored"], evidence["frozenTransferScored"]), (0, 0))
        inside = [r for r in evidence["records"] if r["predictions"]]
        self.assertEqual([r["dye"] for r in inside], ["AMAcRG"])
        self.assertTrue(all(v is None for row in evidence["records"] for v in row["absLogErrors"].values()))
        self.assertEqual({r["dye"] for r in evidence["chenByCondition"]}, {"SMI-NIR1"})
        self.assertEqual(len({(r["pH"], r["powerKwCm2"]) for r in evidence["chenByCondition"]}), 10)
        self.assertTrue(all(r["opticalMidpointApproxPH"] > 4.48 for r in evidence["records"] if r["study"] == "Chen 2024"))
        for bar in evidence["uncertaintyAudit"]["visibleWhiskers"]:
            self.assertLess(abs(bar["approxHalfBar"]-bar["sdOverSqrtN"]), .0015/118)
            self.assertGreater(abs(bar["approxHalfBar"]-bar["workbookSD"]), .0015/118)
