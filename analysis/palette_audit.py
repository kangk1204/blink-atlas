"""Join every Fig. 2/3 replicate to its deposited field, its source cells and its raw-table estimates."""
import argparse
from collections import Counter, defaultdict
import json
import math
from pathlib import Path
from statistics import mean, median


DATA = Path(__file__).resolve().parents[1] / "docs" / "data"
# Explicit reporting grid for the descriptive agreement count, in displayed units.
# These are not error bounds on the experiment or thresholds for selecting a dye.
GRID = {"nenaNm": .01, "centerPhotons": 1, "brightTimeS": .01, "centerSbr": .01,
        "dropPct": .1, "destroyedPct": .1, "specificity": .01}


def audit(catalog, raw, manifest):
    datafiles = {f["path"]: f for f in manifest["files"]}
    if len(datafiles) != len(manifest["files"]):
        raise ValueError("Duplicate path in extraction manifest")
    groups = {s: {d["id"] for d in catalog["dyes"]} for s in ("origami", "npc")}
    expected = {(s, d, n) for s, dyes in groups.items() for d in dyes for n in (1, 2, 3)}
    records = {(r["system"],r["dye"],r["repeat"]): r for r in raw["replicates"]}
    if len(records) != len(raw["replicates"]) or set(records) != expected | {("npc","cy3b",4)}:
        raise ValueError("Principal repeat panel or additional Cy3B field is incomplete")
    referenced = set()
    main_hashes = defaultdict(list)
    fields = {}
    for (s,d,n), record in records.items():
        for f in record["files"]:
            member = datafiles[f["path"]]
            if f["path"] in referenced or member["bytes"] != f["bytes"]:
                raise ValueError("Incorrect file identity or repeat role: " + f["path"])
            if member["datasets"]["locs"]["shape"] != [f["rows"]]:
                raise ValueError("Localization-row count mismatch: " + f["path"])
            referenced.add(f["path"])
            if f["role"] == "target_field":
                fields[s,d,n] = member["sha256"]
                if n <= 3:
                    main_hashes[member["sha256"]].append([s,d,n])
    candidates = {p for p,f in datafiles.items() if f["type"]=="hdf5" and f["section"] in ("fig_2","fig_3")}
    if candidates != referenced or set(fields) != set(records):
        raise ValueError("The analyzed main-panel roles do not cover the source files")
    duplicate_primary = [g for g in main_hashes.values() if len(g)>1]
    all_hdf5 = defaultdict(list)
    for f in datafiles.values():
        if f["type"] == "hdf5":
            all_hdf5[f["sha256"]].append(f["path"])
    principal_in_reference = [
        path for path in all_hdf5[fields["npc","cy3b",4]]
        if "/reference_cy3b/" in path
    ]
    if any(len(x) != 64 for x in all_hdf5):
        raise ValueError("Missing file digest")

    differences = defaultdict(list)
    summary_cells = 0
    for dye in catalog["dyes"]:
        for s, source in dye["systems"].items():
            for metric, worksheet in source["metrics"].items():
                for key, value in (("mean","computedMean"),("sd","computedSd")):
                    a,b = worksheet[key], worksheet[value]
                    if a is None and b is None:
                        continue
                    if a is None or b is None or not math.isclose(a,b,abs_tol=1e-9,rel_tol=1e-10):
                        raise ValueError(f"Published {key} does not reproduce: {dye['id']}/{s}/{metric}")
                    summary_cells += 1
                for p in worksheet["points"]:
                    row = records[s, dye["id"], p["repeat"]]
                    observed = row["metrics"].get(metric)
                    published = p["value"]
                    delta = None if observed is None else observed - published
                    differences[s,metric].append({"dye":dye["id"],"repeat":p["repeat"],"cell":p["cell"],
                                                  "published":published,"recomputed":observed,"difference":delta})
    summaries=[]
    for (s,m), values in sorted(differences.items()):
        usable=[v for v in values if v["difference"] is not None]
        deltas=[v["difference"] for v in usable]
        summaries.append({"system":s,"metric":m,"sourceCells":len(values),"comparable":len(usable),
                          "withinHalfGrid":sum(abs(v)<=GRID[m]/2+1e-9 for v in deltas),
                          "gridUnit":GRID[m],"medianAbsoluteDifference":median(map(abs,deltas)) if deltas else None,
                          "meanSignedDifference":mean(deltas) if deltas else None,
                          "largestDifferences":sorted(usable,key=lambda v:-abs(v["difference"]))[:3],
                          "missing":[{"dye":v["dye"],"repeat":v["repeat"],"cell":v["cell"]} for v in values if v["difference"] is None]})
    return {"schema":"dna-paint-audit/1", "masterArchiveMd5":manifest["archive_md5"],
            "members":len(datafiles),"extractedBytes":sum(f["bytes"] for f in datafiles.values()),
            "hdf5Members":manifest["counts"]["hdf5"],"mainFigureHdf5":len(referenced),
            "mainFigureFileCounts":dict(Counter(datafiles[p]["section"] for p in referenced)),
            "mainFigureCases":len(records), "principalCases":len(expected),
            "primaryFieldDigestDuplicates":duplicate_primary,"cy3bFourthMainFieldAlsoOccursAs":sorted(principal_in_reference),
            "summaryCellsReproduced":summary_cells,
            "agreementRule":"Absolute deviation ≤ half of the stated fixed reporting grid. Descriptive agreement with worksheet values.",
            "metrics":summaries}


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data",type=Path,default=DATA)
    args=parser.parse_args()
    def load(name):
        return json.loads((args.data/name).read_text(encoding="utf8"))
    result=audit(load("palette.json"),load("palette-raw.json"),load("palette-manifest.json"))
    (args.data/"palette-audit.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf8",newline="\n")
