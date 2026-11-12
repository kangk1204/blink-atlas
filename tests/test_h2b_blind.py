import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest

from analysis.analogue_benchmark import PROTOCOL, primary_verdict, summarize
from analysis.dyes import DATA
from analysis.h2b_blind import (
    ROSTER, csv_rows, make_lock, public_lock, save, score, score_files,
)


class BlindedFixedH2B(unittest.TestCase):
    def setUp(self):
        self.tmp = TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.folder = Path(self.tmp.name)
        self.assay = {"study": "test-only", "acquisitionRef": "synthetic test fixture",
                      "independentAcquisition": True, "identitiesReviewedAgainstHolland": True,
                      "outcomesWithheldUntilLock": True, "freeParentTitrationsMatchHolland": True,
                      "illuminationAndDetectionMatched": True,
                      "actualProtocol": dict(PROTOCOL)}

    def cohort(self, *rows):
        return make_lock(self.assay["study"], [
            {"blind_id": dye, "free_parent_pka": value} for dye, value in rows])

    @staticmethod
    def preps(blind_id):
        return [dict(blind_id=blind_id, preparation_id=f"test-{i+1}", mean_on_ms=on, mean_off_s=off,
                     source_locator=f"unit-test-row-{i+1}") for i, (on, off) in enumerate([(1, 2), (2, 2), (3, 1)])]

    def test_frozen_cal664sb_matches_registered_document_and_checksum(self):
        path = DATA / "h2b-cal664sb-lock.json"
        public = public_lock()
        self.assertEqual(json.loads(path.read_bytes()), public)
        digest = save(self.folder / "reproduced.json", public)
        self.assertEqual(digest, "c4e9ae62123abff64ca703a72b8a3c63e29e2ccec312a85f02f969d19eaaca8f")
        self.assertEqual(digest, hashlib.sha256(path.read_bytes()).hexdigest())
        dye = public["dyes"][0]
        self.assertNotIn("observed", dye)
        self.assertEqual(dye["upper"], ["JF639b"])
        self.assertEqual(dye["nearestDyes"], ["JF635b"])
        self.assertEqual(dye["predictions"], {"local": .000579845004467, "nearest": .000483915,
                                               "regression": .000429728198699})
        self.assertAlmostEqual(public["equalErrorBoundaries"][-1]["ratio"], .000529712842337, places=15)
        self.assay.update(study="Cal664sb", freeParentTitrationsMatchHolland=False)
        result = score(public, self.preps("Cal664sb"), self.assay)
        self.assertEqual((result["arm"], result["verdict"]), ("literature buffer transfer", "conditional test only"))

    def test_outcome_free_roster_and_whole_dye_identity(self):
        p = self.folder / "roster.csv"
        p.write_text("blind_id,free_parent_pka,mean_on_ms\ntest-X,3.2,13\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "exactly"):
            csv_rows(p, ROSTER)
        for rows in [[("JF-630b", "3.2")], [("test-X", "3.1"), ("Test_X", "2.8")], [("T9", "nan")]]:
            with self.assertRaises(ValueError):
                self.cohort(*rows)
        lock = self.cohort(("T1", 6.5), ("T2", ""))
        self.assertEqual([d["status"] for d in lock["dyes"]], ["outside training range", "no exact monotonic pKa"])

    def test_receipt_and_acquisition_match_before_unblinding(self):
        lock = self.folder / "locked.json"
        src = self.cohort(("T1", 3.1))
        sha = save(lock, src)
        self.assertEqual(save(lock, src), sha)
        with self.assertRaisesRegex(ValueError, "preserve"):
            save(lock, {**src, "study": "changed"})
        for bad in ("f"*64, ""):
            with self.assertRaisesRegex(ValueError, "SHA-256"):
                score_files(lock, bad, self.folder / "unopened.csv", self.folder / "unopened.json")
        assayfile = self.folder / "assay.json"
        self.assay["actualProtocol"]["offCycleMs"] = 100
        assayfile.write_text(json.dumps(self.assay), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "offCycleMs"):
            score_files(lock, sha, self.folder / "unopened.csv", assayfile)
        self.assay["actualProtocol"] = dict(PROTOCOL)
        self.assay["freeParentTitrationsMatchHolland"] = False
        with self.assertRaisesRegex(ValueError, "freeParentTitrations"):
            score(src, [], self.assay)
        self.assay.update(freeParentTitrationsMatchHolland=True, study="a different study")
        with self.assertRaisesRegex(ValueError, "study ID"):
            score(src, [], self.assay)

    def test_three_paired_preparations_and_all_enrolled_dyes(self):
        lock = self.cohort(("T1", 3.2), ("T2", 2.25), ("T3", 8))
        test = score(lock, self.preps("T1"), self.assay)
        self.assertEqual(test["summary"]["scored"], 1)
        self.assertEqual(test["summary"]["total"], 3)
        self.assertEqual(test["verdict"], "incomplete panel")
        self.assertEqual(test["pendingDyes"], 1)
        self.assertAlmostEqual(test["rows"][0]["observed"], .0015)
        self.assertIsNone(test["summary"]["localVsRegression"]["pairedBootstrap95"])
        for changed in (self.preps("different"), [self.preps("T1")[0]]*2,
                        self.preps("T1") + [dict(self.preps("T1")[0], preparation_id="fourth")],
                        [dict(self.preps("T1")[0], mean_off_s=0)]):
            with self.assertRaises(ValueError):
                score(lock, changed, self.assay)

    def test_command_line_private_panel_round_trip(self):
        roster, outcome, assay, lock, output = (self.folder / n for n in (
            "roster.csv", "out.csv", "assay.json", "locked.json", "scored.json"))
        roster.write_text("blind_id,free_parent_pka\nT1,3.2\n", encoding="utf-8")
        outcome.write_text("blind_id,preparation_id,mean_on_ms,mean_off_s,source_locator\n"
                           "T1,day-A,1,2,test-A\nT1,day-B,2,2,test-B\nT1,day-C,3,1,test-C\n", encoding="utf-8")
        assay.write_text(json.dumps(self.assay), encoding="utf-8")
        cli = [sys.executable, "-m", "analysis.h2b_blind"]
        run = subprocess.run(cli + ["lock", "--study", "test-only", "--roster", str(roster), "--output", str(lock)],
                             check=True, text=True, capture_output=True)
        sha = run.stdout.split()[0]
        self.assertEqual(sha, hashlib.sha256(lock.read_bytes()).hexdigest())
        subprocess.run(cli + ["score", "--lock", str(lock), "--lock-sha256", sha, "--assay", str(assay),
                              "--measurements", str(outcome), "--output", str(output)], check=True, capture_output=True)
        result = json.loads(output.read_bytes())
        self.assertEqual(result["lockSha256"], sha)
        self.assertEqual((result["pendingDyes"], result["summary"]["scored"]), (0, 1))
        self.assertAlmostEqual(result["rows"][0]["observed"], .0015)
        self.assertEqual(result["verdict"], "advantage unconfirmed")

    def test_beating_only_nearest_cannot_pass_the_joint_decision(self):
        rows = [{"status": "scored", "delta": .1, "errors": {"nearest": .2, "regression": .05, "local": .1}}]*4
        summary = summarize(rows, independent=True)
        self.assertGreater(summary["pairedBootstrap95"][0], 0)
        self.assertLess(summary["localVsRegression"]["pairedBootstrap95"][1], 0)
        self.assertEqual(primary_verdict(summary), "advantage unconfirmed")
        both = copy.deepcopy(rows)
        for item in both:
            item["errors"]["regression"] = .3
        self.assertEqual(primary_verdict(summarize(both, independent=True)), "independent review warranted")


if __name__ == "__main__":
    unittest.main()
