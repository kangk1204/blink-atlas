"""Evaluate frozen two-dye policies on Steen et al. NPC repeat three."""
import argparse
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean


DATA = Path(__file__).resolve().parents[1] / "docs" / "data"
PROTOCOL = Path(__file__).with_name("palette-protocol.md")
WEIGHTS = (0, 0.5, 1)
EPSILON_NM = 1e-10


def ranks(values):
    ordered = sorted(values, key=values.get)
    output = {}
    start = 0
    while start < len(ordered):
        stop = start + 1
        while stop < len(ordered) and values[ordered[stop]] == values[ordered[start]]:
            stop += 1
        for dye in ordered[start:stop]:
            output[dye] = (start + 1 + stop) / 2
        start = stop
    return output


def training_ranks(groups, observations, repeats):
    result = {}
    for wavelength, dyes in groups.items():
        result[wavelength] = {
            system: ranks({dye: mean(observations[system, dye, r] for r in repeats) for dye in dyes})
            for system in ("npc", "origami")
        }
    return result


def select(training, weight):
    ranked = {}
    for wavelength, specimens in training.items():
        scores = {dye: (1 - weight) * npc_rank + weight * specimens["origami"][dye]
                  for dye, npc_rank in specimens["npc"].items()}
        ranked[wavelength] = sorted(scores, key=lambda dye: (scores[dye], dye))[:2]
    return ranked


def heldout(groups, observations, repeat, predictions, baselines):
    channels = []
    for wavelength, dyes in sorted(groups.items()):
        values = {dye: observations["npc", dye, repeat] for dye in dyes}
        pair, baseline = predictions[wavelength], baselines[wavelength]
        oracle = sorted(values, key=lambda dye: (values[dye], dye))[:2]
        selected, fixed, best = (mean(values[dye] for dye in two) for two in (pair, baseline, oracle))
        channels.append({"wavelength": wavelength, "pair": pair, "baselinePair": baseline, "oraclePair": oracle,
                         "meanNm": selected, "baselineMeanNm": fixed, "oracleMeanNm": best,
                         "regretNm": selected - best, "baselineRegretNm": fixed - best})
    return channels


def evaluate(groups, observations):
    missing = []
    for dyes in groups.values():
        for dye in dyes:
            for system, repeats in (("origami", (1, 2)), ("npc", (1, 2, 3))):
                for repeat in repeats:
                    v = observations.get((system, dye, repeat))
                    if isinstance(v, bool) or not isinstance(v, (float, int)) or not math.isfinite(v) or v <= 0:
                        missing.append(f"{system}/{dye}/n{repeat}")
    if missing:
        return {"status": "unavailable", "gatePassed": False, "channels": [],
                "reasons": ["Incomplete or invalid NeNA values: " + ", ".join(missing)]}
    first = training_ranks(groups, observations, (1,))
    tuning = []
    for w in WEIGHTS:
        rows = heldout(groups, observations, 2, select(first, w), select(first, 0))
        tuning.append({"weight": w, "macroRegretNm": mean(c["regretNm"] for c in rows), "channels": rows})
    weight = min(tuning, key=lambda t: (round(t["macroRegretNm"], 12), t["weight"]))["weight"]
    pooled = training_ranks(groups, observations, (1, 2))
    pairs, baseline = select(pooled, weight), select(pooled, 0)
    channels = heldout(groups, observations, 3, pairs, baseline)
    regret = mean(c["regretNm"] for c in channels)
    base = mean(c["baselineRegretNm"] for c in channels)
    improved = sum(c["regretNm"] + EPSILON_NM < c["baselineRegretNm"] for c in channels)
    passed = regret + EPSILON_NM < base and improved >= 2
    return {"status": "evaluated", "weight": weight, "tuning": tuning, "gatePassed": passed,
            "channels": channels, "macroRegretNm": regret, "baselineMacroRegretNm": base, "improvedChannels": improved,
            "reasons": [] if passed else ["The held-out results did not meet both prespecified improvement criteria."]}


def records(catalog, raw):
    groups, published = defaultdict(list), {}
    for dye in catalog["dyes"]:
        groups[dye["wavelength"]].append(dye["id"])
        for system, source in dye["systems"].items():
            for point in source["metrics"]["nenaNm"]["points"]:
                published[system, dye["id"], point["repeat"]] = point["value"]
    if {k: len(v) for k, v in groups.items()} != {488: 4, 560: 4, 640: 10}:
        raise ValueError("Wrong dye panel or excitation channels")
    observations = {}
    for row in raw["replicates"]:
        key = (row["system"], row["dye"], row["repeat"])
        if key in observations:
            raise ValueError(f"Duplicate replicate {key}")
        observations[key] = row["metrics"].get("nenaNm")
    return {"schema": "dna-paint-evaluation/1",
            "protocol": {"file": "analysis/palette-protocol.md", "primary": "NPC NeNA (nm)",
                         "sha256": hashlib.sha256(PROTOCOL.read_bytes().replace(b"\r\n", b"\n")).hexdigest(),
                         "tuning": "n1 predicts n2", "test": "n1+n2 predict n3", "weights": list(WEIGHTS),
                         "baseline": "NPC training rank, weight 0", "fourthRepeatExcluded": True,
                         "gate": "Lower macro regret and at least two of three channels strictly improved"},
            "raw": evaluate(groups, observations), "published": evaluate(groups, published)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DATA)
    args = parser.parse_args()
    def load(name):
        return json.loads((args.data / name).read_text(encoding="utf8"))
    out = records(load("palette.json"), load("palette-raw.json"))
    (args.data / "palette-evaluation.json").write_text(json.dumps(out, indent=2, allow_nan=False) + "\n", encoding="utf8")
