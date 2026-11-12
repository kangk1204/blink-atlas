"""Check source reproduction, finite exposure physics and Python/browser parity."""

import json
import math
from pathlib import Path
import shutil
import subprocess
import unittest

import numpy as np

from analysis.dyes import DATA, FILENAME, extract
from analysis.model import (
    DEFAULTS, blink_probabilities, configuration, evaluate_prior, event_rate,
    scenario_probabilities,
)

ROOT = Path(__file__).resolve().parents[1]
DYE_DATA = json.loads((DATA / "dyes.json").read_text(encoding="utf-8"))
DYES = DYE_DATA["dyes"]
SETTINGS = [{}, {"densityMin": 10000, "densityMax": 10500}, {"minPhotons": 2200},
            {"sampleState": "live"}, {"densityMin": 1, "densityMax": 100000, "minPhotons": 1500},
            {"candidatePanel": "validation", "densityMin": 10000, "densityMax": 10500}]


class PhotophysicsTest(unittest.TestCase):
    def test_publisher_records_rebuild(self):
        self.assertEqual(DYE_DATA, extract(DATA / FILENAME))
        self.assertEqual(len(DYES), 9)
        for i, dye in enumerate(DYES, 1):
            self.assertEqual(dye["ligand"], f"{i}HTL")
            self.assertEqual(dye["metrics"]["onMs"]["cells"], f"B{i+1}:D{i+1}")
            self.assertEqual(dye["metrics"]["onOffRatio"]["n"], 3)
        jf630 = DYES[5]["metrics"]
        self.assertAlmostEqual(jf630["onMs"]["mean"], 10.14, delta=0.005)
        self.assertAlmostEqual(jf630["photonsPerMs"]["mean"], 67.2, delta=0.05)
        jf635 = DYES[4]["metrics"]
        self.assertAlmostEqual(jf635["onMs"]["mean"], 14.53, delta=0.005)
        self.assertAlmostEqual(jf635["photonsPerMs"]["mean"], 51.5, delta=0.05)
        rebuilt = configuration(DYES)
        published = json.loads((DATA / "model.json").read_text(encoding="utf-8"))
        # exp() and expm1() can differ by a few ulps across system math libraries.
        # Keep nonnumeric settings and candidate ordering exact.
        lhs = rebuilt["referenceResult"].pop("ranking")
        rhs = published["referenceResult"].pop("ranking")
        self.assertEqual(rebuilt, published)
        self.assertEqual(len(lhs), len(rhs))
        for row, reference in zip(lhs, rhs):
            self.assertEqual(row.keys(), reference.keys())
            for key, value in reference.items():
                if isinstance(value, float):
                    self.assertTrue(math.isclose(row[key], value, rel_tol=1e-13, abs_tol=1e-15), key)
                else:
                    self.assertEqual(row[key], value)

    def test_simulated_stationary_first_burst(self):
        rng = np.random.default_rng(442)
        n, on, off, flux, exposure, threshold = 350_000, 11., 90., 50., 25., 220.
        in_on = rng.random(n) < on / (on + off)
        start = np.where(in_on, 0., rng.exponential(off, n))
        dwell = rng.exponential(on, n)
        lit_time = np.maximum(0, np.minimum(exposure - start, dwell))
        observations = [(lit_time * flux >= threshold).mean(), (lit_time > 0).mean()]
        predictions = blink_probabilities(on, off, flux, exposure, threshold)
        for empirical, expected in zip(observations, predictions):
            self.assertLess(abs(empirical - expected), 5 * math.sqrt(expected * (1 - expected) / n))
        self.assertEqual(blink_probabilities(on, off, flux, exposure, 1251)[0], 0)
        self.assertEqual(*blink_probabilities(on, off, flux, exposure, 0))

    def test_previous_snapshot_overlap_example(self):
        # Independent Bernoulli emitters in one footprint. This checks the
        # earlier frame simulation; it is a different denominator from the
        # candidate-conditioned overlap probability in the recommendation.
        rng = np.random.default_rng(1112)
        for labels, reference in [(10, 0.0042662), (100, 0.264238)]:
            exact = 1 - 0.99**labels - labels * 0.01 * 0.99**(labels - 1)
            self.assertAlmostEqual(exact, reference, delta=1e-6)
            observed = (rng.binomial(labels, 0.01, 200000) >= 2).mean()
            self.assertLess(abs(observed - exact), 5 * math.sqrt(exact * (1 - exact) / 200000))

    def test_scenario_dominance_and_bands(self):
        result = evaluate_prior(DYES, SETTINGS[1])
        pair = result["recommendation"]
        self.assertEqual(set(pair), {"JF630b", "JF629b"})
        for density in np.linspace(10000, 10500, 7):
            for scale in (0.5, 1.0):
                intervals = {}
                for dye in DYES:
                    ps = scenario_probabilities(dye, 25, 100)[scale]
                    rates = [event_rate(density, p, 27) for p in ps]
                    intervals[dye["id"]] = (min(rates), max(rates))
                    row = next(row for row in result["ranking"] if row["id"] == dye["id"])
                    self.assertGreaterEqual(min(rates), row["low"] - 1e-10)
                    self.assertLessEqual(max(rates), row["high"] + 1e-10)
                self.assertGreater(min(intervals[i][0] for i in pair),
                                   max(bounds[1] for i, bounds in intervals.items() if i not in pair))
        self.assertEqual(evaluate_prior(DYES)["status"], "uncertain")
        self.assertEqual(evaluate_prior(DYES, SETTINGS[2])["status"], "low_signal")
        self.assertEqual(evaluate_prior(DYES, SETTINGS[3])["status"], "out_of_scope")

    @unittest.skipUnless(shutil.which("node"), "Node.js is required for the browser parity check")
    def test_browser_against_python(self):
        js = """
          import fs from 'node:fs';
          import {evaluate, curve} from './docs/model-core.js';
          const read = path => JSON.parse(fs.readFileSync(path, 'utf8'));
          const d=read('docs/data/dyes.json'), m=read('docs/data/model.json');
          const settings=JSON.parse(fs.readFileSync(0,'utf8'));
          console.log(JSON.stringify(settings.map(s=>({r:evaluate(d,m,s),c:curve(d,m,s,[100,10000])}))));
        """
        items = json.loads(subprocess.run(["node", "--input-type=module", "-e", js],
                           cwd=ROOT, input=json.dumps(SETTINGS), encoding="utf-8", check=True, capture_output=True).stdout)
        for settings, item in zip(SETTINGS, items):
            expected, actual = evaluate_prior(DYES, settings), item["r"]
            for key in ("status", "recommendation", "contenders"):
                self.assertEqual(expected[key], actual[key])
            self.assertEqual(len(actual["ranking"]), len(expected["ranking"]))
            for a, b in zip(expected["ranking"], actual["ranking"]):
                self.assertEqual(a["id"], b["id"])
                for key in ("rate", "low", "high", "detectProbability", "overlapProbability"):
                    self.assertAlmostEqual(a[key], b[key], delta=max(1e-10, abs(a[key]) * 1e-12))
            if expected["status"] != "out_of_scope":
                s = DEFAULTS | settings
                for point in item["c"]:
                    self.assertEqual(set(point["rates"]), {row["id"] for row in expected["ranking"]})
                    dye = next(d for d in DYES if d["id"] == expected["ranking"][0]["id"])
                    p = blink_probabilities(*(dye["metrics"][k]["mean"] for k in ("onMs","offMs","photonsPerMs")),
                                            s["exposureMs"], s["minPhotons"])
                    self.assertAlmostEqual(point["rates"][dye["id"]], event_rate(point["density"], p, s["frameIntervalMs"]))


if __name__ == "__main__":
    unittest.main()
