"""Finite-exposure forward model and publication-prior recommendations.

Run `python -m analysis.model` to rebuild the public model configuration and
its reference calculation. The browser evaluates these same equations for new
settings and can condition the density range on a matched pilot localization
table (docs/model-core.js).
"""

import argparse
from itertools import product
import json
import math
from pathlib import Path


DATA = Path(__file__).resolve().parents[1] / "docs" / "data"
DEFAULTS = {
    "densityMin": 100, "densityMax": 1000, "minPhotons": 100,
    "exposureMs": 25, "frameIntervalMs": 27, "excitationNm": 647,
    "sampleState": "fixed", "cellLine": "COS-7", "target": "H2B", "tag": "HaloTag",
}
LIMITS = {"density": [1, 100000], "minPhotons": [1, 5000]}
RADIUS_NM = 250
AREA_UM2 = math.pi * (RADIUS_NM / 1000) ** 2
OFF_SCALES = (0.5, 1.0)


def blink_probabilities(on_ms, off_ms, photons_per_ms, exposure_ms, threshold_photons):
    """Return qualifying-burst and any-ON probabilities in one frame.

    Stationary alternating exponential ON/OFF times. For a qualifying
    burst the first ON segment intersecting the frame must supply at least
    threshold/rate ms of ON time. A second burst within one exposure and
    Poisson variation in observed photons are omitted from this proxy.
    """
    if min(on_ms, off_ms, photons_per_ms) <= 0 or min(exposure_ms, threshold_photons) < 0:
        raise ValueError("Dwell times and photon rates must be positive; exposure and threshold nonnegative")
    if not all(math.isfinite(v) for v in (on_ms, off_ms, photons_per_ms, exposure_ms, threshold_photons)):
        raise ValueError("Nonfinite photophysics")
    on_fraction = on_ms / (on_ms + off_ms)

    def first_burst_by(time_ms):
        return on_fraction - (1 - on_fraction) * math.expm1(-time_ms / off_ms)

    any_on = first_burst_by(exposure_ms)
    needed_ms = threshold_photons / photons_per_ms
    if needed_ms > exposure_ms:
        return 0.0, any_on
    qualified = math.exp(-needed_ms / on_ms) * first_burst_by(exposure_ms - needed_ms)
    return qualified, any_on


def event_rate(density, probabilities, frame_interval_ms):
    """Expected isolated, qualifying detections per square micron per second."""
    qualified, any_on = probabilities
    return density * qualified * math.exp(-density * AREA_UM2 * any_on) * 1000 / frame_interval_ms


def kinetics(dye, off_scale, measured_scenarios):
    m = dye["metrics"]
    keys = ("onMs", "offMs", "photonsPerMs")
    if measured_scenarios:
        combinations = product(*(m[k]["replicates"] for k in keys))
    else:
        combinations = [tuple(m[k]["mean"] for k in keys)]
    return [(on, off * off_scale, rate) for on, off, rate in combinations]


def scenario_probabilities(dye, exposure_ms, threshold_photons):
    by_scale = {}
    for scale in OFF_SCALES:
        by_scale[scale] = [blink_probabilities(*parameters, exposure_ms, threshold_photons)
                           for parameters in kinetics(dye, scale, measured_scenarios=True)]
    return by_scale


def domain_supported(settings):
    norm = lambda s: "".join(c for c in str(s).lower() if c.isalnum())
    bio = all(norm(settings[k]) == norm(DEFAULTS[k]) for k in ("sampleState", "cellLine", "target", "tag"))
    return (bio and abs(settings["excitationNm"] - 647) <= 1
            and abs(settings["exposureMs"] - 25) <= 0.2 and abs(settings["frameIntervalMs"] - 27) <= 0.2)


def evaluate_prior(dyes, overrides=None):
    """Rank at the geometric midpoint. Certify an unordered pair over its full range."""
    s = DEFAULTS | (overrides or {})
    lo, hi = s["densityMin"], s["densityMax"]
    h, exposure, cycle = s["minPhotons"], s["exposureMs"], s["frameIntervalMs"]
    if not all(isinstance(v, (int, float)) and math.isfinite(v)
               for v in (lo, hi, h, exposure, cycle, s["excitationNm"])):
        raise ValueError("Nonfinite or nonnumeric model setting")
    if not LIMITS["density"][0] <= lo <= hi <= LIMITS["density"][1]:
        raise ValueError("Density outside model grid")
    if not LIMITS["minPhotons"][0] <= h <= LIMITS["minPhotons"][1] or not 0 < exposure <= cycle:
        raise ValueError("Invalid exposure, camera cycle or required photon count")
    if not domain_supported(s):
        return {"status": "out_of_scope", "recommendation": [], "ranking": [], "contenders": []}
    center = math.sqrt(lo * hi)
    ranking, bounds = [], {}
    for dye in dyes:
        ps = scenario_probabilities(dye, exposure, h)
        endpoints, minima, maxima = {}, [], []
        for scale, scenarios in ps.items():
            at_lo = [event_rate(lo, p, cycle) for p in scenarios]
            at_hi = [event_rate(hi, p, cycle) for p in scenarios]
            endpoints[scale] = [(min(at_lo), max(at_lo)), (min(at_hi), max(at_hi))]
            minima.extend((min(at_lo), min(at_hi)))
            for probabilities in scenarios:
                peak_density = 1 / (AREA_UM2 * probabilities[1])
                maxima.append(event_rate(min(hi, max(lo, peak_density)), probabilities, cycle))
        bounds[dye["id"]] = endpoints
        nominal = blink_probabilities(*kinetics(dye, 1.0, False)[0], exposure, h)
        ranking.append({
            "id": dye["id"], "name": dye["name"], "rate": event_rate(center, nominal, cycle),
            "low": min(minima), "high": max(maxima),
            "detectProbability": nominal[0], "anyOnProbability": nominal[1],
            "overlapProbability": -math.expm1(-center * AREA_UM2 * nominal[1]),
        })
    ranking.sort(key=lambda item: -item["rate"])

    def dominates(i, j):
        # The log rate ratio for two kinetic settings is affine in density.
        # Checking both endpoints proves dominance over the entire interval.
        return all(bounds[i][scale][k][0] > bounds[j][scale][k][1] + 1e-12
                   for scale in OFF_SCALES for k in range(2))

    contenders = [a["id"] for a in ranking if a["high"] > 0 and
                  sum(dominates(b["id"], a["id"]) for b in ranking if b["id"] != a["id"]) < 2]
    pair = [d["id"] for d in ranking[:2]]
    robust_pair = (ranking[1]["low"] > 0 and
                   all(dominates(i, x["id"]) for i in pair for x in ranking if x["id"] not in pair))
    status = "stable" if robust_pair else "uncertain"
    if sum(row["high"] > 0 for row in ranking) < 2:
        status = "low_signal"
    return {
        "status": status, "recommendation": pair if robust_pair else [], "contenders": contenders,
        "ranking": ranking, "density": [lo, hi], "centerDensity": center,
    }


def configuration(dyes):
    return {
        "schema": "smlm-isolation/1", "defaults": DEFAULTS, "limits": LIMITS,
        "collisionRadiusNm": RADIUS_NM, "offTimeScales": list(OFF_SCALES),
        "pilotConfig": {"gridPoints": 181, "credibleMass": 0.9, "maxEffectiveTiles": 512,
                        "maxDeviance": 13.82, "minInformationNats": 0.05},
        "rateUnit": "isolated detections per square micrometre per second",
        "methods": [
            "Two-state blinking at stationary equilibrium; exponential ON and OFF waits represent each reported mean. Effective projected label density is an input assumption.",
            "A dye qualifies when its first ON segment overlaps the exposure long enough to supply the chosen expected photon count. Photon shot noise, camera background, and bleaching during a long run are outside this short-exposure proxy.",
            "Another label with any ON segment during the same exposure and within a 250 nm centre-to-centre disc interferes. Locally uniform independent positions and blinking give a Poisson empty-disc probability. Defining the neighbourhood this way is a modeling choice.",
            "The nominal estimate uses mean dwell times and emission rate. Sensitivity crosses the three measured means for each of the three kinetic inputs within each dye (27 combinations). Replicates from differing acquisition protocols are not paired. A common dark-time scale of 1 or 0.5 tests protocol sensitivity across dyes; dye-specific shortening or unmeasured extremes remain open.",
            "Bands span the specified sensitivities and density interval. The rate reaches its maximum at density 1/(disc area × probability of any ON time). Top two are reported only when every other candidate is outscored throughout the interval at both dark-time settings. These are sensitivity bounds.",
            "A registered fixed-condition pilot is accepted only when its acquisition and photon filter match. For each trial density the likelihood jointly assigns observed tiles to empty, single and multiple localization categories. The total count adds no second likelihood. The working Poisson tile model approximates residual spatial dependence after overlap rejection. One field worth of tiles, capped at 512 effective observations, limits pseudo-replication.",
            "The pilot starts from a log-uniform density grid and marginalizes equally over the stated reference-dye sensitivity scenarios. The 90% density interval is conditional on the observation model. Multimodal support is retained in the interval. Pilot data inconsistent with all sensitivity scenarios trigger abstention.",
        ],
        "referenceResult": evaluate_prior(dyes),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--dyes", type=Path, default=DATA / "dyes.json")
    parser.add_argument("--output", type=Path, default=DATA / "model.json")
    args = parser.parse_args()
    dyes = json.loads(args.dyes.read_text(encoding="utf-8"))["dyes"]
    args.output.write_text(json.dumps(configuration(dyes), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Saved model assumptions and default result to {args.output}")


if __name__ == "__main__":
    main()
