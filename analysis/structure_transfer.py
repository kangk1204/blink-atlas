"""Audit and improve free-dye pKcycl transfer across original experiments using CPU and RDKit."""

import argparse
from functools import lru_cache
import json
from pathlib import Path
from statistics import mean

import numpy as np
from rdkit import Chem

from analysis import structure_pk as v1
from analysis import transfer_sources as src


OUT = v1.OUTPUT.parent / "structure-transfer"
PLAN = v1.REG / "transfer-registration.json"
NAMES = ("siliconBridge", "aminoClosure", "pendantCO2H", "crowdedAryl3", "crowdedAryl6",
         "donorSubstitution", "acylatedDonors", "cyclicDonors", "fusedDonors", "fourMemberedDonors",
         "donorWithdrawal")
METHODS = ("nearest", "kernel", "coreDonor", "withoutWithdrawal", "trainMean")
ELECTRONEGATIVITY = {7: 3.04, 8: 3.44, 9: 3.98, 17: 3.16, 35: 2.96, 53: 2.66}
CARBON_EN = 2.55
# Aromatic donor ipso atoms are labeled explicitly so assignment is independent of atom order.
CORE = Chem.MolFromSmarts("[O,Si:1]1c2c[c:2]ccc2[C;X4:3]2([c:4]3[c:5]ccc[c:6]3[CH2:7][O,N,S:8]2)c2cc[c:9]cc21")
CARBOXYL = Chem.MolFromSmarts("c-[C](=O)[O;H1,-]")


def branch_atoms(donor, core_indices):
    """Traverse the donor side of the aryl–N bond. Stop on either fused aryl boundary."""
    frontier = [(donor, 0)]
    seen = set(core_indices) | {donor.GetIdx()}
    for atom, distance in frontier:
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in seen:
                frontier.append((neighbor, distance + 1))
                seen.add(neighbor.GetIdx())
    return frontier


def is_nitrile_n(atom):
    return atom.GetAtomicNum() == 7 and any(b.GetBondType() == Chem.BondType.TRIPLE for b in atom.GetBonds())


def is_carbonyl(atom):
    return atom.GetAtomicNum() == 6 and any(b.GetOtherAtom(atom).GetAtomicNum() == 8
        and b.GetBondType() == Chem.BondType.DOUBLE for b in atom.GetBonds())


@lru_cache(maxsize=512)
def features(smiles):
    mol = v1.smiles_mol(smiles)
    matches = mol.GetSubstructMatch(CORE)
    if not matches or Chem.GetFormalCharge(mol) != 0 or len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("Expected one neutral, closed (Si-)rhodamine free parent.")
    core = set(matches)
    labeled = {atom.GetAtomMapNum(): mol.GetAtomWithIdx(matches[atom.GetIdx()])
               for atom in CORE.GetAtoms() if atom.GetAtomMapNum()}
    donors = []
    for pos in (2, 9):
        found = [n for n in labeled[pos].GetNeighbors() if n.GetIdx() not in core and n.GetAtomicNum() == 7]
        if len(found) != 1:
            raise ValueError("Expected two aryl amino donors at the rhodamine 3,6 positions.")
        donors.extend(found)
    # Atoms adjacent to each end of the pendant aryl/spiro fusion, excluding the two ring junctions.
    pos3 = next(n for n in labeled[6].GetNeighbors() if n.GetIsAromatic() and n.GetIdx() != labeled[4].GetIdx())
    crowded = lambda a: float(any(n.GetIdx() not in core and n.GetAtomicNum() > 1 for n in a.GetNeighbors()))
    pendant_ring = next(set(ring) for ring in mol.GetRingInfo().AtomRings()
                        if labeled[4].GetIdx() in ring and all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring))
    acid = any(m[0] in pendant_ring for m in mol.GetSubstructMatches(CARBOXYL))
    donor_rows = []
    for n in donors:
        arms = [a for a in n.GetNeighbors() if not a.GetIsAromatic() and a.GetAtomicNum() > 1]
        cycles = [ring for ring in mol.GetRingInfo().AtomRings() if n.GetIdx() in ring]
        withdrawal = sum((ELECTRONEGATIVITY[a.GetAtomicNum()] - CARBON_EN) * (0.5 ** (d - 1))
                         for a, d in branch_atoms(n, core) if d and (a.GetAtomicNum() in (8, 9, 17, 35, 53) or is_nitrile_n(a)))
        donor_rows.append([len(arms), any(is_carbonyl(a) for a in arms), bool(cycles),
                           any(core.intersection(ring) for ring in cycles), any(len(ring) == 4 for ring in cycles),
                           withdrawal])
    return (float(labeled[1].GetAtomicNum() == 14), float(labeled[8].GetAtomicNum() == 7), float(acid),
            crowded(pos3), crowded(labeled[5]), *np.mean(donor_rows, axis=0).tolist())


class ChemicalRidge:
    def __init__(self, parents, columns=NAMES):
        self.columns = [NAMES.index(n) for n in columns]
        matrix = np.array([features(p["smiles"]) for p in parents])[:, self.columns]
        self.center = matrix.mean(axis=0)
        scale = matrix.std(axis=0)
        self.scale = np.where(scale > 1e-12, scale, 1.)
        x = (matrix - self.center) / self.scale
        x[:, scale < 1e-12] = 0
        labels = np.array([p["measured"] for p in parents])
        self.mean = float(labels.mean())
        self.weights = np.linalg.solve(x.T @ x + np.eye(len(self.columns)), x.T @ (labels - self.mean))

    def predict(self, values):
        query = np.asarray(values)[self.columns]
        return v1.rounded(self.mean + ((query - self.center) / self.scale) @ self.weights)

    def snapshot(self):
        return dict(names=[NAMES[j] for j in self.columns], center=list(map(v1.rounded, self.center)),
                    scale=list(map(v1.rounded, self.scale)), weights=list(map(v1.rounded, self.weights)),
                    intercept=v1.rounded(self.mean))


class Fit:
    def __init__(self, parents):
        self.parents = parents
        self.morgan = v1.Predictor(parents)
        self.chemical = ChemicalRidge(parents)
        self.ablated = ChemicalRidge(parents, NAMES[:-1])
        self.by_key = {p["inchiKey"]: p for p in parents}

    def model_id(self):
        compiled = dict(train=self.parents, descriptors=self.chemical.snapshot(), rdkit=v1.rdBase.rdkitVersion,
                        inputs=[features(p["smiles"]) for p in self.parents])
        return v1.digest(PLAN.read_bytes() + json.dumps(compiled, sort_keys=True).encode())

    def verification(self):
        m = self.morgan
        return dict(bits=[list(fp.GetOnBits()) for fp in m.fps], intercept=v1.rounded(m.mean),
                    means=list(map(v1.rounded, m.column_means)), grandMean=v1.rounded(m.grand_mean),
                    weights=list(map(v1.rounded, m.weights)), ablated=self.ablated.snapshot())

    def score(self, smiles):
        mol = v1.smiles_mol(smiles)
        old = self.morgan.predict(mol)
        x = features(smiles)
        pred = dict(nearest=old["nearest"], kernel=old["ridge"], coreDonor=self.chemical.predict(x),
                    withoutWithdrawal=self.ablated.predict(x), trainMean=v1.rounded(self.chemical.mean))
        neighbors = [dict(**p, silicon=bool(features(self.by_key[p["inchiKey"]]["smiles"])[0])) for p in old["neighbors"]]
        return dict(predictions=pred, similarity=old["similarity"], neighbors=neighbors,
                    descriptors=dict(zip(NAMES, x, strict=True)))


def comparison(test):
    def stats(method):
        signed = [r["predictions"][method] - r["measured"] for r in test]
        return dict(mae=v1.rounded(mean(abs(e) for e in signed)), bias=v1.rounded(mean(signed)))
    errors = {name: stats(name) for name in METHODS}
    differences = [r["errors"]["nearest"] - r["errors"]["coreDonor"] for r in test]
    reduced = [{"omit": r["dye"], "gain": v1.rounded(mean(differences[:i] + differences[i+1:]))}
               for i, r in enumerate(test) if len(test) > 1]
    gain = v1.rounded(mean(differences))
    return dict(n=len(test), metrics=errors, gain=gain,
                wins=sum(d > 1e-9 for d in differences), ties=sum(abs(d) <= 1e-9 for d in differences),
                losses=sum(d < -1e-9 for d in differences), deletion=reduced,
                stableGain=gain > 1e-9 and all(r["gain"] > 1e-9 for r in reduced))


def parent_fold(rows, study, external=False, within=False, *, excluded=None, targets=None):
    excluded_studies = set(excluded or (src.EXTERNAL if external else (study,)))
    all_blocked = {r["inchiKey"] for r in rows if r["study"] in excluded_studies and r["inchiKey"]}
    if targets is None:
        targets = [r for r in rows if r["study"] == study and r["used"]]
    if len(set(r["inchiKey"] for r in targets)) != len(targets):
        raise ValueError("Repeated test identity")
    out = {"study": study, "excludedStudies": [] if within else sorted(excluded_studies), "within": within, "rows": []}
    def fit_excluding(blocked):
        train = v1.training(rows, blocked=blocked, omit=() if within else excluded_studies)
        if blocked.intersection(p["inchiKey"] for p in train):
            raise ValueError("Train/test parent leakage.")
        if not within and set(s for p in train for s in p["studies"]) & excluded_studies:
            raise ValueError("Test paper in training.")
        return Fit(train)
    fold_model = None if within else fit_excluding(all_blocked)
    if fold_model is not None:
        out.update(train=fold_model.parents, blocked=sorted(all_blocked), model=fold_model.chemical.snapshot(), overlap=0)
    for r in targets:
        fit = fold_model if fold_model is not None else fit_excluding({r["inchiKey"]})
        measured = float(r["reported"])
        s = fit.score(r["smiles"])
        info = dict(study=r["study"], dye=r["dye"], inchiKey=r["inchiKey"], measured=measured, **s,
                    errors={m: v1.rounded(abs(pred - measured)) for m, pred in s["predictions"].items()})
        if within:
            info["trainingKeys"] = sorted(fit.by_key)
        out["rows"].append(info)
    out["summary"] = comparison(out["rows"]) if targets else None
    return out


def audit_original(rows):
    old = v1.study_fold(rows, "Holland")
    vals = [r["measured"] for r in old["train"]]
    tm = [r["measured"] for r in old["rows"]]
    bykey = {r["inchiKey"]: r for r in rows}
    for r in old["rows"]:
        r["nearestSigned"] = v1.rounded(r["nearest"] - r["measured"])
        r["kernelSigned"] = v1.rounded(r["ridge"] - r["measured"])
        r["nearestBridge"] = sorted({"Si" if features(bykey[p["inchiKey"]]["smiles"])[0] else "O" for p in r["neighbors"]})
    return dict(rows=old["rows"], nTrain=len(vals), trainMean=mean(vals), trainRange=[min(vals), max(vals)],
                testRange=[min(tm), max(tm)], similarityRange=[min(r["similarity"] for r in old["rows"]),
                max(r["similarity"] for r in old["rows"])], belowTrain=sum(x < min(vals) for x in tm),
                siTrain=sum(bool(features(p["smiles"])[0]) for p in old["train"]),
                azetidineTrain=sum(features(p["smiles"])[NAMES.index("fourMemberedDonors")] > 0 for p in old["train"]),
                priorMAE=old["mae"])


def registry():
    base = v1.catalogue()
    return base, base + src.additions()


def new_parents(rows, base_keys):
    studies = ("Sakabe", *src.EXTERNAL)
    new = [r for r in rows if r["used"] and r["study"] in studies and r["inchiKey"] not in base_keys]
    return parent_fold(rows, "New parents", excluded=studies, targets=new)


def build(out=OUT):
    out.mkdir(parents=True, exist_ok=True)
    base, rows = registry()
    folds = [parent_fold(rows, study, external=study in src.EXTERNAL) for study in ("Holland", "Tyson", *src.EXTERNAL)]
    general = (any(f["summary"] for f in folds if f["study"] in src.EXTERNAL)
               and all(f["summary"]["stableGain"] for f in folds if f["summary"]))
    preferred = "coreDonor" if general else "nearest"
    ablation = [parent_fold(base, study) for study in ("Holland", "Tyson")]
    within = parent_fold(rows, "Holland", within=True)
    train = v1.training(rows)
    fit = Fit(train)
    model_id = fit.model_id()
    assays = src.additions()
    bundle = dict(schema="blink-atlas/structure-transfer/1", modelId=model_id, preferred=preferred,
                  baseParents=len(v1.training(base)), training=train, trainingParents=len(train),
                  audit=audit_original(base), folds=folds, oldDataFolds=ablation, withinHolland=within,
                  fullModel=fit.chemical.snapshot(), additions=assays, sources=src.DOIS, candidates=[],
                  verification=fit.verification(), observed=observations(rows))
    # Controlled removal: how much of the Holland improvement survives without the closest new reference?
    bundle["hollandWithoutChi"] = parent_fold([r for r in rows if r["study"] != src.CHI], "Holland")
    base_keys = {p["inchiKey"] for p in v1.training(base)}
    bundle["newParents"] = new_parents(rows, base_keys)
    extra = [dict(r, used=True) if r["status"] == "readout unresolved" else r for r in rows]
    bundle["opticalSensitivity"] = [parent_fold(extra, s, external=s in src.EXTERNAL)
                                    for s in ("Holland", "Tyson", *src.EXTERNAL)] + [new_parents(extra, base_keys)]
    spec = json.loads((v1.REG / "registration.json").read_text())
    for c in spec["candidates"]:
        query = query_result(c["closedSmiles"], c["id"], bundle, rows, fit)
        bundle["candidates"].append(query)
        v1.dump(out / (c["id"] + ".json"), query)
    v1.dump(out / "benchmark.json", bundle)
    v1.csv_file(out / "assays.csv", assays)
    per_dye = []
    cases = [("expanded", f) for f in folds] + [("original", f) for f in ablation]
    cases.extend([("within series", within), ("new parents combined", bundle["newParents"]),
                  ("expanded, Chi removed", bundle["hollandWithoutChi"])])
    cases.extend(("optical sensitivity", f) for f in bundle["opticalSensitivity"])
    for data, case in cases:
        for r in case["rows"]:
            per_dye.append(dict(panel=data, split=case["study"], study=r["study"], dye=r["dye"], inchiKey=r["inchiKey"],
                    measured=r["measured"], similarity=r["similarity"], neighbors=[n["dye"] for n in r["neighbors"]],
                    **r["predictions"], **{k+"AbsError": e for k,e in r["errors"].items()}))
    v1.csv_file(out / "per-dye.csv", per_dye)
    return bundle


def observations(rows):
    found = v1.assay_observations(rows)
    for row in rows:
        if row["status"] == "readout unresolved":
            found.setdefault(row["inchiKey"], []).append({k: row[k] for k in
                ("study", "dye", "reported", "status", "locator")})
    return found


def query_result(smiles, name, bundle, rows, fit):
    mol = v1.smiles_mol(smiles)
    # Use the same free-parent screening as the existing CLI.
    legacy = v1.structure_result(mol, name, fit.morgan, "nearest", bundle["modelId"], observations(rows))
    predictions = None
    if legacy["status"] != "observed complex or limited":
        predictions = fit.score(smiles)
    return dict(schema="blink-atlas/transfer-query/1", modelId=bundle["modelId"], dye=name, **v1.identity(mol),
                status=legacy["status"], knownMeasured=legacy["knownMeasured"], observations=legacy["observations"],
                preferred=bundle["preferred"], fingerprint=legacy["fingerprint"], result=predictions)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    make = sub.add_parser("build")
    make.add_argument("--out", type=Path, default=OUT)
    query = sub.add_parser("predict")
    g = query.add_mutually_exclusive_group(required=True)
    g.add_argument("--smiles")
    g.add_argument("--sdf", type=Path)
    query.add_argument("--name")
    query.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.command == "build":
            b = build(args.out)
            for f in b["folds"]:
                s = f["summary"]
                print(f["study"], {m: (s["n"], round(v["mae"], 4)) for m, v in s["metrics"].items()} if s else "0 eligible")
            print("Selected", b["preferred"])
        else:
            if args.sdf:
                mols = list(Chem.SDMolSupplier(str(args.sdf)))
                if len(mols) != 1 or mols[0] is None:
                    raise ValueError("Expected one valid SDF record.")
                mol = mols[0]
            else:
                mol = v1.smiles_mol(args.smiles)
            name = args.name or (mol.GetProp("_Name") if mol.HasProp("_Name") else "User structure")
            smi = Chem.MolToSmiles(v1.parent_input(mol))
            _, rows = registry()
            bundle = json.loads((OUT / "benchmark.json").read_text(encoding="utf-8"))
            fit = Fit(v1.training(rows))
            if fit.model_id() != bundle["modelId"]:
                raise ValueError("Release and source model differ. Rebuild the benchmark before predicting.")
            v1.dump(args.output, query_result(smi, name, bundle, rows, fit))
    except (OSError, ValueError) as exc:
        parser.exit(2, str(exc) + "\n")


if __name__ == "__main__":
    main()
