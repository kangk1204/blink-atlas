"""Structure-only free-dye pKcycl benchmark, release builder and local prediction CLI."""

import argparse
import csv
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
from statistics import mean

import numpy as np
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import AllChem, Draw, rdFingerprintGenerator, rdMolDescriptors

from analysis.structure_catalog import DOI, catalog, quantum_rows


ROOT = Path(__file__).resolve().parents[1]
REG = ROOT / "validation/structure_pk"
OUTPUT = ROOT / "docs/data/structure-pk"
SCHEMA = "blink-atlas/free-parent-pk/1"
GENERATOR = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True)
# Closed ortho-(hydroxymethyl/aminomethyl/mercaptomethyl) xanthene or Si-xanthene.
SCAFFOLD = Chem.MolFromSmarts("[O,Si]1c2ccccc2[C;X4]2(c3ccccc3C[O,N,S]2)c2ccccc21")


def dump(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                    encoding="utf-8", newline="\n")


def digest(value):
    return hashlib.sha256(value).hexdigest()


def rounded(n):
    return round(float(n), 10)


def key(mol):
    return Chem.MolToInchiKey(mol)


def smiles_mol(s):
    mol = Chem.MolFromSmiles(s)
    if mol is None or not mol.GetNumAtoms():
        raise ValueError("Invalid or empty molecular graph.")
    return mol


def parent_input(mol):
    fragments = Chem.GetMolFrags(mol, asMols=True)
    # Only remove explicit spectators. An organic mixture needs a user-selected parent.
    spectators = {"[Cl-]", "[Br-]", "[I-]", "[Na+]", "[K+]", "[Li+]", "O"}
    parents = [f for f in fragments if Chem.MolToSmiles(f) not in spectators]
    if len(parents) != 1:
        raise ValueError("Expected one free-dye parent; select a single component from organic mixtures.")
    return parents[0]


def identity(mol):
    return dict(inchiKey=key(mol), smiles=Chem.MolToSmiles(mol),
                formula=rdMolDescriptors.CalcMolFormula(mol, separateIsotopes=True))


def catalogue():
    rows = catalog()
    for dye, smiles, measured, _ in quantum_rows():
        # The published comparison combines upper branches, earlier measurements and protocol-development rows.
        origin = "Uno or earlier; row-level provenance unresolved"
        reason = "Parent or original numeric assay not reconciled with the Table 1 quotation."
        if dye in ("HMRG", "AMRG", "HMTMR", "AMTMR", "HMRB"):
            origin, reason = "Uno", "Matching primary Fig. 2 measurement catalogued separately."
        elif dye in ("HMDiMeR", "HMDiEtR"):
            origin, reason = "Sakabe", "Original spectra in SI Fig. S1; value quoted again."
        elif dye in ("HMSiR620h", "HMJSiR"):
            origin, reason = "Iwatate 2018", "SI Table S1 has two transitions. Tachibana quotes the upper one."
        elif dye.startswith(("HMJ", "HMDiMeF", "HMDiMeC")):
            origin, reason = "Iwatate 2016", "Secondary Table 1 verified. Primary numeric panel unavailable for independent checking."
        elif dye in ("HMSiR", "AMSiR"):
            reason = "Table 1 draws an unsubstituted pendant aryl. Uno Fig. 2 draws a carboxylated parent with a complex curve."
        rows.append(dict(study="Tachibana", dye=dye, smiles=smiles, reported=str(measured), status="quoted",
                         reason=reason, method="Literature pKcycl", buffer="Original conditions (see origin)",
                         locator="Table 1, " + dye, source="https://doi.org/" + DOI["Tachibana"],
                         origin=origin, evidence="quotation", form="closed"))
    for r in rows:
        mol = smiles_mol(r["smiles"])
        r.update(identity(mol))
        # Explicit isotope labels are also represented in the canonical graph and SDF.
        isotopes = Counter(f"{a.GetIsotope()}{a.GetSymbol()}" for a in mol.GetAtoms() if a.GetIsotope())
        r["isotopes"] = ", ".join(f"{label} ×{count}" for label, count in sorted(isotopes.items())) or "No explicit isotope label"
        r.setdefault("substitution", "Atom-resolved canonical SMILES and cited structure drawing")
        r["used"] = r["status"] == "exact" and r["method"] == "Absorbance" and r["evidence"] == "original"
        if r["method"] == "Fluorescence" and r["status"] == "exact":
            r["reason"] = "Same paper and parent as its preferred absorbance measurement."
    return rows


def training(rows, omit=(), blocked=()):
    grouped = defaultdict(dict)
    for r in rows:
        if not r["used"] or r["study"] in omit or r["inchiKey"] in blocked:
            continue
        by_study = grouped[r["inchiKey"]]
        if r["study"] in by_study:
            raise ValueError("Repeated training label within an original parent/study.")
        by_study[r["study"]] = r
    out = []
    for k, sources in sorted(grouped.items()):
        first = next(iter(sources.values()))
        out.append(dict(inchiKey=k, dye=first["dye"], smiles=first["smiles"],
                        studies=sorted(sources), measured=mean(float(s["reported"]) for s in sources.values())))
    return out


class Predictor:
    def __init__(self, parents):
        self.parents = parents
        self.fps = [GENERATOR.GetFingerprint(smiles_mol(p["smiles"])) for p in parents]
        targets = np.array([p["measured"] for p in parents])
        self.mean = targets.mean()
        kernel = np.array([DataStructs.BulkTanimotoSimilarity(fp, self.fps) for fp in self.fps])
        self.column_means = kernel.mean(axis=0)
        self.grand_mean = kernel.mean()
        centered = kernel - self.column_means - self.column_means[:, None] + self.grand_mean
        self.weights = np.linalg.solve(centered + np.eye(len(parents)), targets - self.mean)

    def predict(self, mol):
        fp = GENERATOR.GetFingerprint(mol)
        similarities = np.array(DataStructs.BulkTanimotoSimilarity(fp, self.fps))
        best = similarities.max()
        ties = [p for p, sim in zip(self.parents, similarities, strict=True) if abs(sim - best) <= 1e-12]
        centered = similarities - self.column_means - similarities.mean() + self.grand_mean
        return dict(nearest=rounded(mean(p["measured"] for p in ties)),
                    ridge=rounded(self.mean + centered @ self.weights), similarity=rounded(best),
                    neighbors=[{k: p[k] for k in ("dye", "inchiKey", "studies", "measured")} for p in ties])


def study_fold(rows, study):
    assay = [r for r in rows if r["study"] == study]
    blocked = {r["inchiKey"] for r in assay}
    train = training(rows, omit=[study], blocked=blocked)
    overlap = blocked & {p["inchiKey"] for p in train}
    if overlap:
        raise ValueError("Target structure in training.")
    model = Predictor(train)
    test = []
    for row in assay:
        if not row["used"]:
            continue
        label = float(row["reported"])
        p = model.predict(smiles_mol(row["smiles"]))
        test.append(dict(dye=row["dye"], inchiKey=row["inchiKey"], measured=label, **p,
                         nearestError=rounded(abs(p["nearest"] - label)), ridgeError=rounded(abs(p["ridge"] - label))))

    def maes(sample):
        return {m: rounded(mean(x[m + "Error"] for x in sample)) for m in ("nearest", "ridge")}

    overall = maes(test)
    deletions = [{"omit": r["dye"], **maes([t for t in test if t is not r])} for r in test]
    return dict(study=study, blockedKeys=sorted(blocked), train=train, overlap=len(overlap), rows=test,
                eligible=len(test), total=len(blocked), mae=overall, leaveOneOut=deletions,
                stableRidge=all(m["ridge"] < m["nearest"] for m in [overall, *deletions]))


def quantum(rows):
    out = []
    for i, (name, smi, observed, published) in enumerate(quantum_rows()):
        mol = smiles_mol(smi)
        p = Predictor(training(rows, blocked=[key(mol)])).predict(mol)
        exact = any(r["used"] and r["inchiKey"] == key(mol) and float(r["reported"]) == observed for r in rows)
        stage = ("Hydration changed to five waters after comparison" if i < 2 else
                 "Hydration protocol development" if i < 8 else "Subsequent application in Table 1")
        out.append(dict(dye=name, inchiKey=key(mol), measured=observed, quantum=published,
                        originalThreeWater=([11.3, 10.1][i] if i < 2 else published), stage=stage,
                        tablePrintedError=[0.2, 0.0, 0.0, 0.3, 0.1, 0.1, 0.5, 0.6, 0.3, 0.5,
                                           0.4, 0.1, 0.4, 0.0, 0.4, 0.2, 0.4][i],
                        matchedPrimary=exact, **p, quantumError=rounded(abs(published - observed)),
                        nearestError=rounded(abs(p["nearest"] - observed)), ridgeError=rounded(abs(p["ridge"] - observed))))
    return out


def csv_file(path, records):
    def scalar(v):
        return json.dumps(v, ensure_ascii=False, separators=(",", ":")) if isinstance(v, (list, dict)) else v
    with path.open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(records[0]), lineterminator="\n")
        w.writeheader()
        for r in records:
            w.writerow({k: scalar(v) for k, v in r.items()})


def sdf(path, records):
    writer = Chem.SDWriter(str(path))
    for r in records:
        mol = smiles_mol(r["smiles"])
        AllChem.Compute2DCoords(mol)
        mol.SetProp("_Name", r["dye"])
        for k in ("inchiKey", "formula", "form", "aliases", "source", "locator"):
            if k in r:
                mol.SetProp(k, str(r[k]))
        writer.write(mol)
    writer.close()
    # RDKit pads SD property headers. End-of-line padding carries no molfile columns.
    lines = path.read_text(encoding="utf-8").splitlines()
    path.write_text("".join(line.rstrip() + "\n" for line in lines), encoding="utf-8", newline="\n")


def model_id(parents):
    return digest((REG / "registration.json").read_bytes() + json.dumps(parents, sort_keys=True).encode())


def assay_observations(rows):
    observed = defaultdict(list)
    for row in rows:
        if row["evidence"] == "original" and row["method"] == "Absorbance":
            observed[row["inchiKey"]].append({field: row[field] for field in
                ("study", "dye", "reported", "status", "locator")})
    return dict(observed)


def structure_result(mol, name, model, preferred, version, observations):
    if len(Chem.GetMolFrags(mol)) != 1 or Chem.GetFormalCharge(mol) != 0 or not mol.HasSubstructMatch(SCAFFOLD):
        raise ValueError("Supply one neutral, ring-closed free (Si-)rhodamine parent with an ortho CH2–O/N/S spiro ring. Remove salts and ligands.")
    if any(a.GetAtomicNum() not in {1, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53} for a in mol.GetAtoms()):
        raise ValueError("Unsupported element in free dye.")
    if mol.HasSubstructMatch(Chem.MolFromSmarts("ClCCCCCC")):
        raise ValueError("A chloroalkyl HaloTag ligand is outside the free-parent assay.")
    p = model.predict(mol)
    known = next((r for r in model.parents if r["inchiKey"] == key(mol)), None)
    measurements = observations.get(key(mol), [])
    status = "known exact" if known else "observed complex or limited" if measurements else "unmeasured in corpus"
    abstain = status == "observed complex or limited"
    return dict(schema=SCHEMA, modelId=version, dye=name, **identity(mol), preferred=preferred,
                status=status, observations=measurements,
                fingerprint=list(GENERATOR.GetFingerprint(mol).GetOnBits()),
                nearest=None if abstain else p["nearest"], ridge=None if abstain else p["ridge"],
                similarity=p["similarity"], neighbors=p["neighbors"],
                modelGap=None if abstain else rounded(abs(p["ridge"] - p["nearest"])),
                knownMeasured=known["measured"] if known else None,
                endpoint="Free-parent pKcycl", trainingParents=len(model.parents))


def context():
    rows = catalogue()
    train = training(rows)
    folds = [study_fold(rows, study) for study in ("Holland", "Tyson")]
    selected = "ridge" if all(f["stableRidge"] for f in folds) else "nearest"
    return rows, train, folds, selected


def build(out):
    out.mkdir(parents=True, exist_ok=True)
    rows, train, folds, selected = context()
    model = Predictor(train)
    observed = assay_observations(rows)
    version = model_id(train)
    spec = json.loads((REG / "registration.json").read_text(encoding="utf-8"))
    rounds = []
    control = next(r for r in rows if r["dye"] == "JF626b")
    for name, smi, desc in [(c["id"], c["closedSmiles"], c["substituent"]) for c in spec["candidates"]] + [
            ("JF626b", control["smiles"], "3-(trifluoromethyl)azetidine on both arms")]:
        mol = smiles_mol(smi)
        result = structure_result(mol, name, model, selected, version, observed)
        dump(out / (name + ".json"), result)
        is_control = name == "JF626b"
        sdf(out / (name + ".sdf"), [dict(dye=name, **identity(mol), form="closed",
            source="https://doi.org/" + DOI["Holland"] if is_control else "Blink Atlas frozen virtual design")])
        reloaded = list(Chem.SDMolSupplier(str(out / (name + ".sdf"))))
        if len(reloaded) != 1 or structure_result(reloaded[0], name, model, selected, version, observed) != result:
            raise ValueError("SDF round-trip changed identity or prediction.")
        AllChem.Compute2DCoords(mol)
        drawer = Draw.MolDraw2DSVG(520, 265)
        drawer.drawOptions().clearBackground = False
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        (out / (name + ".svg")).write_text(drawer.GetDrawingText(), encoding="utf-8", newline="\n")
        rounds.append(dict(**result, substitution=desc, role="measured control" if is_control else "virtual candidate"))
    unique = {}
    for r in rows:
        unique.setdefault(r["inchiKey"], []).append(r)
    parents = []
    for rs in unique.values():
        parents.append(dict(**rs[0], aliases=" | ".join(dict.fromkeys(r["study"] + ": " + r["dye"] for r in rs))))
    sdf(out / "parents.sdf", parents)
    csv_file(out / "measurements.csv", rows)
    csv_file(out / "holdouts.csv", [dict(study=f["study"], **r) for f in folds for r in f["rows"]])
    qr = quantum(rows)
    csv_file(out / "quantum.csv", qr)
    qr_summary = {}
    for label, sample in (("publishedTable", qr), ("primaryMatched", [r for r in qr if r["matchedPrimary"]])):
        qr_summary[label] = {"n": len(sample), **{m: rounded(mean(r[m + "Error"] for r in sample))
                                                  for m in ("nearest", "ridge", "quantum")}}
    snapshot = dict(parents=[dict(**p, bits=list(fp.GetOnBits())) for p, fp in zip(train, model.fps, strict=True)],
                    weights=[rounded(n) for n in model.weights], mean=rounded(model.mean),
                    columnMeans=[rounded(n) for n in model.column_means], grandMean=rounded(model.grand_mean),
                    observed=observed)
    counts = {}
    for s in dict.fromkeys(r["study"] for r in rows):
        sub = [r for r in rows if r["study"] == s]
        counts[s] = dict(represented=len({r["inchiKey"] for r in sub}),
                         exactOriginal=len({r["inchiKey"] for r in sub if r["used"]}),
                         exclusions=[dict(dye=r["dye"], status=r["status"], reason=r["reason"]) for r in sub if not r["used"]])
    payload = dict(schema="blink-atlas/structure-release/1", modelId=version, preferred=selected, rdkit=rdBase.rdkitVersion,
                   uniqueParents=len(parents), trainingParents=len(train), counts=counts, folds=folds,
                   quantum=qr, quantumSummary=qr_summary, candidates=rounds, verification=snapshot,
                   sources={s: "https://doi.org/" + doi for s, doi in DOI.items()})
    dump(out / "benchmark.json", payload)
    return payload


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    make = sub.add_parser("build")
    make.add_argument("--out", type=Path, default=OUTPUT)
    predict = sub.add_parser("predict")
    group = predict.add_mutually_exclusive_group(required=True)
    group.add_argument("--smiles")
    group.add_argument("--sdf", type=Path)
    predict.add_argument("--name")
    predict.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.command == "build":
            bundle = build(args.out)
            for f in bundle["folds"]:
                print(f["study"], f["eligible"], f["mae"], "overlap", f["overlap"])
            print("Selected:", bundle["preferred"])
            return
        if args.sdf:
            molecules = list(Chem.SDMolSupplier(str(args.sdf)))
            if len(molecules) != 1 or molecules[0] is None:
                raise ValueError("Expected exactly one valid SDF record.")
            mol = molecules[0]
        else:
            mol = smiles_mol(args.smiles)
        name = args.name or (mol.GetProp("_Name") if mol.HasProp("_Name") else "User structure")
        mol = parent_input(mol)
        rows, train, _, selected = context()
        dump(args.output, structure_result(mol, name, Predictor(train), selected, model_id(train), assay_observations(rows)))
    except (OSError, ValueError) as exc:
        parser.exit(2, str(exc) + "\n")


if __name__ == "__main__":
    main()
