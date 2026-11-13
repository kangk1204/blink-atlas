import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

try:
    import rdkit
except ImportError:
    rdkit = None

if rdkit:
    import numpy as np
    from rdkit import Chem
    from analysis import structure_transfer as transfer
    from analysis.structure_catalog import ACID, NME, SI, arm, graph


@unittest.skipUnless(rdkit, "Install analysis/requirements-structure.txt")
class TransferModel(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.base, cls.rows = transfer.registry()
        cls.bundle = json.loads((transfer.OUT / "benchmark.json").read_text(encoding="utf-8"))

    def test_graph_features_and_donor_boundaries(self):
        def fs(s):
            return dict(zip(transfer.NAMES, transfer.features(s), strict=True))
        plain = fs(graph(NME))
        silicon = fs(graph(NME, bridge=SI, pendant=ACID))
        self.assertEqual((plain["siliconBridge"], silicon["siliconBridge"], silicon["pendantCO2H"]), (0, 1, 1))
        self.assertEqual(silicon["donorWithdrawal"], 0)  # pendant acid is outside the N branches
        rows = {r["dye"]: r for r in self.rows if r["used"]}
        for name, ring4, fused in [("JF635b", 1, 0), ("HMSiRTHQ", 0, 1), ("HMTMR", 0, 0)]:
            x = fs(rows[name]["smiles"])
            self.assertEqual((x["fourMemberedDonors"], x["fusedDonors"]), (ring4, fused))
        self.assertEqual(fs(rows["3MHMRG"]["smiles"])["crowdedAryl3"], 1)
        self.assertEqual(fs(rows["5MHMRG"]["smiles"])["crowdedAryl3"], 0)
        self.assertEqual(fs(graph(arm("N3CC(C#N)C3"), bridge=SI))["donorWithdrawal"], (3.04-2.55)/8)
        mol = Chem.MolFromSmiles(rows["JF639b"]["smiles"])
        reverse = Chem.RenumberAtoms(mol, list(reversed(range(mol.GetNumAtoms()))))
        self.assertEqual(fs(Chem.MolToSmiles(reverse, canonical=False)), fs(rows["JF639b"]["smiles"]))

    def test_novel_identities_and_study_block(self):
        b = self.bundle
        self.assertEqual([b["baseParents"], b["trainingParents"]], [26, 29])
        self.assertEqual([len(f["rows"]) for f in b["folds"]], [6, 5, 0, 1])
        self.assertEqual({r["dye"] for r in b["newParents"]["rows"]}, {"HM-DS655","Ac-HMDiEtR","Ac-HMDiMeR"})
        new_keys = {r["inchiKey"] for r in b["newParents"]["rows"]}
        self.assertFalse(new_keys.intersection(r["inchiKey"] for r in transfer.v1.training(self.base)))
        for fold in b["folds"] + b["oldDataFolds"] + b["opticalSensitivity"] + [b["newParents"], b["hollandWithoutChi"]]:
            train = fold["train"]
            self.assertFalse(set(fold["blocked"]).intersection(p["inchiKey"] for p in train))
            self.assertFalse(set(fold["excludedStudies"]).intersection(s for p in train for s in p["studies"]))
        for row in b["withinHolland"]["rows"]:
            self.assertNotIn(row["inchiKey"], row["trainingKeys"])
        statuses = {r["dye"]: r["used"] for r in b["additions"]}
        self.assertFalse(statuses["HMbisTFER"] or statuses["HM-DS531"] or statuses["HMtetTFER"])
        self.assertEqual(b["opticalSensitivity"][-1]["summary"]["n"], 4)

    def test_training_only_standardization_and_ablation(self):
        orig = self.bundle["folds"][0]
        varied = [dict(r, reported="80") if r["study"] == "Holland" and r["used"] else r for r in self.rows]
        other = transfer.parent_fold(varied, "Holland")
        self.assertEqual(orig["model"], other["model"])
        self.assertEqual([r["predictions"] for r in orig["rows"]], [r["predictions"] for r in other["rows"]])
        # A feature present only in the test molecule cannot acquire a coefficient.
        m = transfer.ChemicalRidge([{"smiles":graph(NME), "measured":6}, {"smiles":graph(arm("N")), "measured":8}])
        self.assertEqual(m.weights[0], 0)
        self.assertEqual(m.predict(transfer.features(graph(NME))), m.predict(transfer.features(graph(NME,bridge=SI))))

    def test_paired_scoring_and_removal(self):
        b = self.bundle
        audit = b["audit"]
        self.assertEqual((audit["siTrain"],audit["azetidineTrain"],audit["belowTrain"]), (5,0,4))
        self.assertTrue(all(r["nearestBridge"] == ["O"] and r["kernelSigned"] > 0 for r in audit["rows"]))
        self.assertTrue(all([n["dye"] for n in r["neighbors"]] == ["HM-DS655"] for r in b["folds"][0]["rows"]))
        self.assertEqual(b["hollandWithoutChi"]["summary"]["metrics"]["nearest"]["mae"], audit["priorMAE"]["nearest"])
        self.assertEqual(b["preferred"], "nearest")
        for f in [b["newParents"], *b["folds"], b["withinHolland"]]:
            if not f["rows"]:
                continue
            for key in transfer.METHODS:
                err = [r["predictions"][key] - r["measured"] for r in f["rows"]]
                self.assertAlmostEqual(float(np.mean(np.abs(err))), f["summary"]["metrics"][key]["mae"])
                self.assertAlmostEqual(float(np.mean(err)), f["summary"]["metrics"][key]["bias"])
        self.assertEqual(b["newParents"]["summary"]["wins"], 3)
        self.assertLess(b["newParents"]["summary"]["metrics"]["trainMean"]["mae"], b["newParents"]["summary"]["metrics"]["coreDonor"]["mae"])

    def test_regeneration_queries_and_chemical_scope(self):
        b = self.bundle
        fit = transfer.Fit(transfer.v1.training(self.rows))
        self.assertEqual(fit.model_id(), b["modelId"])
        with TemporaryDirectory() as d:
            out = Path(d)
            transfer.build(out)
            for f in transfer.OUT.iterdir():
                self.assertEqual(f.read_bytes(), (out / f.name).read_bytes(), f.name)
        for released in b["candidates"]:
            mol = Chem.SDMolSupplier(str(transfer.v1.OUTPUT/(released["dye"] + ".sdf")))[0]
            self.assertEqual(transfer.query_result(Chem.MolToSmiles(mol), released["dye"], b, self.rows, fit), released)
        for name in ("HM-SiR", "JF614b", "HMbisTFER"):
            r = next(r for r in self.rows if r["dye"] == name)
            self.assertIsNone(transfer.query_result(r["smiles"], name, b, self.rows, fit)["result"])
        charged = next(r for r in self.rows if r["dye"] == "HM-DS531")
        self.assertEqual(Chem.GetFormalCharge(Chem.MolFromSmiles(charged["smiles"])), 2)
        with self.assertRaises(ValueError):
            transfer.query_result(charged["smiles"], "charged", b, self.rows, fit)
