import copy
import csv
import json
import math
from pathlib import Path
import tempfile
import unittest

from analysis.analogue_benchmark import PROTOCOL, SCOPE, benchmark, build, independent_panel, observations, predict
from analysis.dyes import DATA
from analysis.probe_design import build as holland


class AnalogueBenchmarks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = build()
        cls.sets = {s["id"]: s for s in cls.result["sets"]}
        cls.public = observations(DATA / "analogue-observations.csv")

    def test_release_is_reproducible_and_discloses_untested_primary(self):
        self.assertEqual(self.result, json.loads((DATA / "analogue-benchmark.json").read_text(encoding="utf-8")))
        self.assertEqual((self.result["uniqueExternalDyes"], self.result["externalConditions"]), (11, 25))
        s = self.sets["matched"]["summary"]
        self.assertEqual((s["total"], s["scored"], s["mae"], s["pairedBootstrap95"]), (0, 0, None, None))
        self.assertEqual(self.result["primaryVerdict"], "advantage unconfirmed")
        for group in self.sets.values():
            for row in group["rows"]:
                self.assertNotIn(row["id"], row["trainingIds"])
                if row["status"] == "scored":
                    for method, y in row["predictions"].items():
                        self.assertAlmostEqual(row["errors"][method], abs(math.log10(y / row["observed"])), places=10)
                else:
                    self.assertNotIn("predictions", row)

    def test_external_results_and_change_of_training_set(self):
        frozen = self.sets["remmel-igg-frozen"]
        within = self.sets["remmel-igg-local"]
        for group, expected_ids, maes, wins in (
            (frozen, [1, 2, 3, 4], [.292964, .307467, .207599], 2),
            (within, [2, 3, 4, 5], [.162958, .360401, .291799], 4),
        ):
            scored = [r for r in group["rows"] if r["status"] == "scored"]
            self.assertEqual([r["id"] for r in scored], [f"remmel-{i}" for i in expected_ids])
            self.assertEqual(group["summary"]["wins"], wins)
            self.assertEqual(group["summary"]["total"], 6)
            for method, value in zip(("local", "nearest", "regression"), maes):
                self.assertAlmostEqual(group["summary"]["mae"][method], value, places=6)
        self.assertTrue(all(len(r["trainingIds"]) == 6 and "JF630b" in r["trainingIds"] for r in frozen["rows"]))
        self.assertTrue(all(len(r["trainingIds"]) == 5 and all(x.startswith("remmel-") for x in r["trainingIds"])
                            for r in within["rows"]))
        lo, hi = frozen["summary"]["pairedBootstrap95"]
        self.assertAlmostEqual(lo, -.0080472534)
        self.assertAlmostEqual(hi, .0424571136)
        self.assertIsNone(within["summary"]["pairedBootstrap95"])

    def test_source_units_and_same_parent_in_distinct_formulations(self):
        datasets = [self.sets[f"remmel-{k}-frozen"] for k in ("igg", "halotag7", "nanobody")]
        same_parent = [next(r for r in s["rows"] if r["id"] == "remmel-4") for s in datasets]
        for row, on_ms, off_s in zip(same_parent, [3, 1.6, 2.5], [1.02, 8.95, 1.16]):
            self.assertAlmostEqual(row["observed"], on_ms / (1000 * off_s), places=12)
            self.assertEqual(row["pKa"], 4.38)
            self.assertEqual(row["predictions"], same_parent[0]["predictions"])
        self.assertEqual(self.sets["tyson-glass-frozen"]["rows"][0]["status"], "outside training range")
        self.assertEqual(self.sets["remmel-halotag7-local"]["rows"][0]["status"], "insufficient training dyes")
        for power in (500, 1000, 1500, 2000):
            local = self.sets[f"zheng-{power}-local"]
            self.assertEqual({r["id"] for r in local["rows"] if r["status"] == "scored"}, {"stmr", "srhb"})
            self.assertEqual(local["summary"]["wins"], 1)
            for r in self.sets[f"zheng-{power}-frozen"]["rows"]:
                self.assertEqual((r["status"], r["rangeStatus"]), ("different response window", "outside training range"))
        stmr = next(r for r in self.public["zheng", "500"] if r["id"] == "stmr")
        self.assertEqual(stmr["pKa"], 7.18)
        self.assertAlmostEqual(stmr["ratio"], .026 / (1 - .026), places=12)

    def test_held_out_response_cannot_change_its_own_prediction(self):
        rows = copy.deepcopy(self.public["remmel", "IgG"])
        rows[2]["ratio"] *= 100
        alt = benchmark("test", "", "", "transfer-local", rows)["rows"][2]
        original = self.sets["remmel-igg-local"]["rows"][2]
        self.assertEqual(original["predictions"], alt["predictions"])
        self.assertNotEqual(original["errors"], alt["errors"])
        with self.assertRaises(ValueError):
            benchmark("overlap", "", "", "frozen", rows, reference=[rows[0]])

    def test_optional_bound_titration_uses_paired_dyes_and_identical_training_folds(self):
        covariate = self.result["exploratoryTitrations"]
        self.assertIn("post hoc", covariate["status"])
        self.assertEqual([r["power"] for r in covariate["byPower"]], [500, 1000, 1500, 2000])
        for group in covariate["byPower"]:
            self.assertEqual({r["id"] for r in group["rows"] if r["paired"]}, {"srhb", "stmr"})
            for r in group["rows"]:
                self.assertEqual(r["free"]["observed"], r["bound"]["observed"])
                self.assertEqual(r["free"]["trainingIds"], r["bound"]["trainingIds"])
                self.assertNotIn(r["id"], r["bound"]["trainingIds"])
                self.assertNotEqual(r["free"]["pKa"], r["bound"]["pKa"])
        wins = [g["power"] for g in covariate["byPower"] if g["mae"]["bound"]["local"] < g["mae"]["free"]["local"]]
        self.assertEqual(wins, [1000])

    def test_ties_exact_match_and_abstentions(self):
        rows = [dict(id=n, pKa=p, ratio=r) for n, p, r in [("a", 2, .01), ("b", 2, .0001), ("c", 4, .1)]]
        mid = predict(rows, 3)
        self.assertEqual(mid["nearestDyes"], ["a", "b"])
        self.assertAlmostEqual(mid["predictions"]["nearest"], .001)
        self.assertAlmostEqual(mid["predictions"]["local"], .01)
        self.assertEqual(predict(rows, 2)["predictions"]["local"], .001)
        self.assertEqual(predict(rows[:2], 2)["status"], "insufficient distinct pKa values")
        self.assertEqual(predict(rows, 1)["status"], "outside training range")

    def test_invalid_or_ambiguous_source_values(self):
        raw = DATA / "analogue-observations.csv"
        with raw.open(encoding="utf-8") as f:
            example = next(csv.DictReader(f))
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.csv"
            for edit in (dict(on_ms="0"), dict(pka="nan"), dict(duty="0.9"), dict(on_ms="", off_s="", duty="1")):
                with path.open("w", newline="", encoding="utf-8") as f:
                    writer = csv.DictWriter(f, fieldnames=list(example))
                    writer.writeheader()
                    writer.writerow({**example, **edit})
                with self.assertRaises(ValueError):
                    observations(path)


class MatchedPanelIntake(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.catalog = holland()

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / "unit-test-panel.json"
        self.dye = {"parentInchiKey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C", "name": "test-only-01", "aliases": [],
                    "monotonicFreeTitration": True, "pKa": 3.2, "preparationRatios": [.001, .002, .003],
                    "sourceLocator": "Synthetic unit test, row 1"}
        self.panel = {"schema": "blink-atlas/matched-h2b/1", "scope": SCOPE, "protocol": dict(PROTOCOL),
                      "source": "https://example.org/unit-test-only", "independentAcquisition": True,
                      "identitiesReviewedAgainstHolland": True, "dyes": [self.dye]}

    def save(self):
        self.path.write_text(json.dumps(self.panel), encoding="utf-8")

    def test_import_preserves_units_replicates_provenance_and_original_training(self):
        self.save()
        result = build(matched=self.path)
        self.assertEqual(result["matchedProtocol"], PROTOCOL)
        row = result["sets"][0]["rows"][0]
        self.assertEqual(row["trainingIds"], [p["id"] for p in self.catalog["calibration"]["pairs"]])
        self.assertEqual(row["measurement"]["preparationRatios"], self.dye["preparationRatios"])
        self.assertEqual(row["observed"], .002)
        self.assertIsNone(result["sets"][0]["summary"]["pairedBootstrap95"])
        del self.dye["preparationRatios"]
        self.dye["reportedMeanRatio"] = .004
        self.save()
        rows, source, _ = independent_panel(self.path, self.catalog)
        self.assertEqual((rows[0]["ratio"], source), (.004, self.panel["source"]))
        self.assertEqual(rows[0]["measurement"]["sourceLocator"], self.dye["sourceLocator"])

    def test_reject_known_dyes_duplicate_parents_and_incompatible_protocol(self):
        original = copy.deepcopy(self.panel)
        edits = [lambda p: p["dyes"][0].update(aliases=["jf-639B"]),
                 lambda p: p["dyes"].append({**p["dyes"][0], "name": "alias of same parent"}),
                 lambda p: p["protocol"].update(offCycleMs=100), lambda p: p.update(independentAcquisition=False),
                 lambda p: p.update(source="https://doi.org/10.1038/s41592-026-03062-5"),
                 lambda p: p["dyes"][0].update(reportedMeanRatio=.2),
                 lambda p: p["dyes"][0].update(preparationRatios=[0]),
                 lambda p: p["dyes"][0].update(pKa=None),
                 lambda p: p["dyes"][0].update(monotonicFreeTitration=False),
                 lambda p: p["dyes"][0].update(sourceLocator=" ")]
        for edit in edits:
            self.panel = copy.deepcopy(original)
            edit(self.panel)
            self.save()
            with self.assertRaises(ValueError):
                independent_panel(self.path, self.catalog)


if __name__ == "__main__":
    unittest.main()
