"""Scientific release checks specific to the structure-only endpoint (optional RDKit environment)."""

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest

try:
    import rdkit  # noqa: F401
except ImportError:
    rdkit = None

if rdkit is not None:
    import numpy as np
    from rdkit import Chem, DataStructs
    from analysis import structure_pk as pk


@unittest.skipIf(rdkit is None, "Install analysis/requirements-structure.txt for the structure benchmark")
class StructureEndpoint(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows, cls.train, cls.folds, cls.selected = pk.context()
        cls.by_name = {(r["study"], r["dye"]): r for r in cls.rows if r["method"] == "Absorbance"}
        cls.published = json.loads((pk.OUTPUT / "benchmark.json").read_text(encoding="utf-8"))

    def test_original_aliases_isotopes_and_exclusions(self):
        def row(study, dye):
            return self.by_name[study, dye]
        self.assertEqual(row("Uno", "HMSiR")["inchiKey"], row("Tyson", "HMSiR")["inchiKey"])
        self.assertNotEqual(row("Uno", "HMSiR")["inchiKey"], row("Holland", "HM-SiR")["inchiKey"])
        jfx = pk.smiles_mol(row("Holland", "JFX650b")["smiles"])
        self.assertEqual(sum(a.GetIsotope() == 2 for a in jfx.GetAtoms()), 16)
        self.assertEqual(row("Tyson", "Cal664sb")["formula"], "C33H34F4N2O3Si")
        self.assertNotEqual(row("Tachibana", "3MHMRG")["inchiKey"], row("Tachibana", "5MHMRG")["inchiKey"])
        prior = json.loads((pk.OUTPUT.parent / "ratio-evidence.json").read_text(encoding="utf-8"))
        fluorescence = {r["dye"]: str(r["siFluorescencePKcycl"]) for r in prior["records"] if r["study"] == "Uno 2014"}
        self.assertEqual(fluorescence, {r["dye"]: r["reported"] for r in self.rows
                                      if r["study"] == "Uno" and r["method"] == "Fluorescence"})
        measured = next(p for p in self.train if p["dye"] == "HMRG")
        self.assertEqual(measured["studies"], ["Sakabe", "Uno"])
        self.assertFalse(any(r["used"] for r in self.rows if r["status"] != "exact" or r["method"] == "Fluorescence"))
        self.assertEqual(len(pk.training(self.rows)), 26)
        self.assertEqual(len([r for r in self.rows if r["study"] == "Tachibana" and r["evidence"] == "quotation"]), 17)

    def test_split_includes_all_target_identities(self):
        for f in self.folds:
            test_keys = {r["inchiKey"] for r in self.rows if r["study"] == f["study"]}
            self.assertEqual(set(f["blockedKeys"]), test_keys)
            self.assertFalse(set(p["inchiKey"] for p in f["train"]) & test_keys)
            self.assertTrue(all(f["study"] not in p["studies"] for p in f["train"]))
            self.assertEqual(len(f["rows"]), len({r["inchiKey"] for r in f["rows"]}))
            self.assertEqual(f["mae"]["nearest"], pk.rounded(np.mean([abs(r["nearest"] - r["measured"]) for r in f["rows"]])))
        self.assertEqual([f["eligible"] for f in self.folds], [6, 5])
        self.assertTrue(all(f["mae"]["ridge"] < f["mae"]["nearest"] for f in self.folds))
        self.assertEqual([r["omit"] for r in self.folds[1]["leaveOneOut"] if r["nearest"] < r["ridge"]], ["Yale676sb"])
        self.assertEqual(self.selected, "nearest")

    def test_tied_structures_and_independent_ridge_solve(self):
        parents = [dict(dye=str(i), inchiKey=str(i), studies=["test"], smiles=s, measured=y)
                   for i, (s, y) in enumerate([("CO", 2), ("CO", 4), ("CCO", 7)])]
        model = pk.Predictor(parents)
        mol = pk.smiles_mol("CO")
        got = model.predict(mol)
        self.assertEqual(got["nearest"], 3)
        self.assertEqual(len(got["neighbors"]), 2)
        m = np.array([DataStructs.BulkTanimotoSimilarity(fp, model.fps) for fp in model.fps])
        # Solve the intercept-augmented KRR directly, independently of kernel centering.
        augmented = np.block([[m + np.eye(3), np.ones((3, 1))], [np.ones((1, 3)), np.zeros((1, 1))]])
        result = np.linalg.solve(augmented, [2, 4, 7, 0])
        self.assertAlmostEqual(got["ridge"], m[0] @ result[:3] + result[3], places=9)

    def test_support_and_known_complex_traces(self):
        def score(smiles):
            return pk.structure_result(pk.smiles_mol(smiles), "sample", pk.Predictor(self.train),
                                       self.selected, pk.model_id(self.train), pk.assay_observations(self.rows))
        for name in ("HM-SiR", "JF646b", "JF614b"):
            result = score(self.by_name["Holland", name]["smiles"])
            self.assertIsNone(result["ridge"])
            self.assertIsNone(result["nearest"])
            self.assertTrue(result["observations"])
        for smi in ("c1ccccc1", self.by_name["Holland", "JF626b"]["smiles"] + ".[Na+]",
                    self.by_name["Uno", "2MeSiR"]["smiles"]):
            with self.assertRaises(ValueError):
                score(smi)
        control = self.by_name["Holland", "JF626b"]
        hydrated = pk.parent_input(pk.smiles_mol(control["smiles"] + ".O.[Na+].[Cl-]"))
        self.assertEqual(pk.key(hydrated), control["inchiKey"])
        with self.assertRaises(ValueError):
            pk.parent_input(pk.smiles_mol(control["smiles"] + ".c1ccccc1"))

    def test_rebuild_and_cli_round_trip(self):
        with TemporaryDirectory() as work:
            out = Path(work) / "release"
            pk.build(out)
            self.assertEqual(sorted(f.name for f in out.iterdir()), sorted(f.name for f in pk.OUTPUT.iterdir()))
            for f in out.iterdir():
                self.assertEqual(f.read_bytes(), (pk.OUTPUT / f.name).read_bytes(), f.name)
            with (out / "measurements.csv").open(encoding="utf-8", newline="") as handle:
                by_key = {r["inchiKey"] for r in csv.DictReader(handle)}
            supplier = Chem.SDMolSupplier(str(out / "parents.sdf"))
            self.assertEqual({pk.key(m) for m in supplier}, by_key)
            self.assertEqual(len(supplier), len(by_key))
            dest = Path(work) / "input.json"
            for row in self.published["candidates"]:
                n = row["dye"]
                base = [sys.executable, "-m", "analysis.structure_pk", "predict", "--output", str(dest)]
                for args in (["--sdf", str(out / (n + ".sdf"))], ["--smiles", row["smiles"], "--name", n]):
                    subprocess.run(base + args, cwd=pk.ROOT, check=True)
                    self.assertEqual(dest.read_bytes(), (out / (n + ".json")).read_bytes())
            with self.assertRaises(subprocess.CalledProcessError):
                subprocess.run(base + ["--sdf", str(out / "parents.sdf")], cwd=pk.ROOT, capture_output=True, check=True)

    def test_quantum_protocol_and_registration(self):
        q = self.published["quantum"]
        self.assertEqual(len(q), 17)
        self.assertEqual(sum(r["matchedPrimary"] for r in q), 7)
        self.assertEqual([(r["originalThreeWater"], r["quantum"]) for r in q[:2]], [(11.3, 7.9), (10.1, 6.2)])
        expected = {"h2b-cal664sb-lock.json": "c4e9ae62123abff64ca703a72b8a3c63e29e2ccec312a85f02f969d19eaaca8f",
                    "registration.json": "5a0fb0e189d4ac9a36dbda179301b5b52583bad504beea1eae0ec38e699f3af9"}
        for file in (pk.OUTPUT.parent / "h2b-cal664sb-lock.json", pk.REG / "registration.json"):
            self.assertEqual(hashlib.sha256(file.read_bytes()).hexdigest(), expected[file.name])
