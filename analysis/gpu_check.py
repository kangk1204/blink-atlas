"""Optional CUDA Monte Carlo verification of the finite-exposure formula.

Run with CUDA-enabled PyTorch: python -m analysis.gpu_check
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import subprocess
import time

from analysis.model import DATA, DEFAULTS, RADIUS_NM, blink_probabilities, event_rate


SEED = 20261112
SINGLE_EXPOSURES = 8_388_608
NEIGHBORHOODS = 65_536
DENSITIES = (100.0, 1000.0)
Z_LIMIT = 5.0


def gpu_snapshot(device_index):
    fields = "name,utilization.gpu,memory.free,driver_version"
    result = subprocess.run(
        ["nvidia-smi", "-i", str(device_index), f"--query-gpu={fields}",
         "--format=csv,noheader,nounits"],
        check=True, capture_output=True, text=True,
    )
    name, util, free, driver = next(csv.reader([result.stdout.strip()], skipinitialspace=True))
    return {"name": name, "utilization_percent": int(util),
            "free_memory_mib": int(free), "driver": driver}


def first_segment(torch, n, on_ms, off_ms, exposure_ms, generator, device):
    """Sample the stationary initial state and first intersecting bright segment."""
    starts_on = torch.rand(n, generator=generator, device=device) < on_ms / (on_ms + off_ms)
    wait_to_on = torch.empty(n, device=device).exponential_(1 / off_ms, generator=generator)
    wait_to_on.masked_fill_(starts_on, 0)
    bright_lifetime = torch.empty(n, device=device).exponential_(1 / on_ms, generator=generator)
    return torch.minimum((exposure_ms - wait_to_on).clamp_min_(0), bright_lifetime)


def trajectory_counts(torch, inputs, generator, device):
    hits = qualifying = 0
    on, off, exposure, rate, threshold = inputs
    for start in range(0, SINGLE_EXPOSURES, 262_144):
        n = min(262_144, SINGLE_EXPOSURES - start)
        length = first_segment(torch, n, on, off, exposure, generator, device)
        hits += (length > 0).sum().item()
        qualifying += (length * rate >= threshold).sum().item()
    return hits, qualifying


def neighbor_counts(torch, density, inputs, generator, device):
    on, off, exposure, _, _ = inputs
    radius_um = RADIUS_NM / 1000
    square_mean = density * (2 * radius_um) ** 2
    empty_disks = 0
    for start in range(0, NEIGHBORHOODS, 2048):
        n = min(2048, NEIGHBORHOODS - start)
        # Simulate all possible neighbors in an enclosing square. The radius and
        # XY positions determine interference; no analytical thinning is used.
        counts = torch.poisson(torch.full((n,), square_mean, device=device),
                               generator=generator).long()
        ends = counts.cumsum(0)
        num_labels = int(ends[-1].item())
        xy = (torch.rand((2, num_labels), device=device, generator=generator) * 2 - 1) * radius_um
        inside_disc = xy.square().sum(0) <= radius_um ** 2
        lengths = first_segment(torch, num_labels, on, off, exposure, generator, device)
        interference = inside_disc & (lengths > 0)
        prefix = torch.cat((torch.zeros(1, dtype=torch.long, device=device),
                            interference.cumsum(0)))
        per_disc = prefix[ends] - prefix[ends - counts]
        empty_disks += (per_disc == 0).sum().item()
    return empty_disks


def check(estimate, expected, standard_error):
    error = abs(estimate - expected)
    return {"monte_carlo": estimate, "analytic": expected,
            "standard_error": standard_error, "absolute_error": error,
            "z_score_abs": error / standard_error if standard_error else None,
            "passed": error <= Z_LIMIT * standard_error}


def proportion(hits, n, expected):
    return {**check(hits / n, expected, math.sqrt(expected * (1 - expected) / n)),
            "successes": hits, "trials": n}


def verify(dyes, device_index=0):
    # Availability is sampled immediately before the compute context is opened.
    before = gpu_snapshot(device_index)
    if before["utilization_percent"] > 10 or before["free_memory_mib"] < 1024:
        raise RuntimeError(f"Selected GPU is busy or lacks room for this check: {before}")
    import torch

    if not torch.cuda.is_available():
        raise RuntimeError("This verification requires CUDA-enabled PyTorch.")
    device = torch.device(f"cuda:{device_index}")
    torch.cuda.set_device(device_index)
    torch.empty(1, device=device)  # Initialize the device allocator before resetting its counters.
    generator = torch.Generator(device=device).manual_seed(SEED)
    torch.cuda.reset_peak_memory_stats(device)
    start = time.perf_counter()
    checked = datetime.now(timezone.utc).isoformat(timespec="seconds")
    exposure, interval, threshold = (DEFAULTS[k] for k in ("exposureMs", "frameIntervalMs", "minPhotons"))
    conditions = []
    checks = []
    by_id = {d["id"]: d for d in dyes}
    for name in ("JF630b", "JF635b"):
        metrics = by_id[name]["metrics"]
        on, off, rate = (metrics[k]["mean"] for k in ("onMs", "offMs", "photonsPerMs"))
        inputs = (on, off, exposure, rate, threshold)
        expected_qualified, expected_any = blink_probabilities(on, off, rate, exposure, threshold)
        any_count, qualified_count = trajectory_counts(torch, inputs, generator, device)
        any_result = proportion(any_count, SINGLE_EXPOSURES, expected_any)
        qualified_result = proportion(qualified_count, SINGLE_EXPOSURES, expected_qualified)
        checks.extend([any_result, qualified_result])
        scenarios = []
        for density in DENSITIES:
            empty_count = neighbor_counts(torch, density, inputs, generator, device)
            expected_empty = math.exp(-density * math.pi * (RADIUS_NM / 1000) ** 2 * expected_any)
            empty_result = proportion(empty_count, NEIGHBORHOODS, expected_empty)
            # Two independent simulation streams form a factorized estimator.
            # This gives useful precision for rare qualifying dyes. Include
            # both sampling variances and their product in the standard error.
            v1, v2 = qualified_result["standard_error"] ** 2, empty_result["standard_error"] ** 2
            product_variance = expected_empty ** 2 * v1 + expected_qualified ** 2 * v2 + v1 * v2
            scale = density * 1000 / interval
            yield_result = check(
                scale * qualified_result["monte_carlo"] * empty_result["monte_carlo"],
                event_rate(density, (expected_qualified, expected_any), interval),
                scale * math.sqrt(product_variance),
            )
            checks.extend([empty_result, yield_result])
            scenarios.append({"density_per_um2": density, "no_active_neighbor": empty_result,
                              "isolated_rate_per_um2_s": yield_result})
        conditions.append({"dye": name, "mean_on_ms": on, "mean_off_ms": off,
                           "photons_per_ms": rate, "any_on": any_result,
                           "qualified_first_burst": qualified_result, "densities": scenarios})
    torch.cuda.synchronize(device)
    elapsed = time.perf_counter() - start
    peak = torch.cuda.max_memory_allocated(device) / 1024 ** 2
    torch.cuda.empty_cache()
    return {
        "schema": "finite-exposure-cuda-check/1", "checked_utc": checked,
        "device": torch.cuda.get_device_name(device), "cuda_runtime": torch.version.cuda,
        "torch_version": str(torch.__version__), "before": before,
        "wall_seconds": round(elapsed, 3), "peak_memory_mib": round(peak, 3),
        "peak_memory_scope": "PyTorch allocated tensors; CUDA context excluded",
        "seed": SEED, "single_exposures_per_dye": SINGLE_EXPOSURES,
        "disks_per_density_per_dye": NEIGHBORHOODS, "exposure_ms": exposure,
        "frame_interval_ms": interval, "threshold_photons": threshold,
        "exclusion_radius_nm": RADIUS_NM, "off_time_scale": 1, "z_limit": Z_LIMIT,
        "check_count": len(checks), "passed": all(c["passed"] for c in checks),
        "conditions": conditions, "source_doi": "10.1038/s41592-026-03062-5",
        "method": (
            "Independent stationary ON/OFF exponential trajectories sample the time to the "
            "first bright segment and its exposure-truncated duration. First-burst detection "
            "uses deterministic photons accrued in that segment. Independent unthinned Poisson "
            "labels are sampled uniformly in a square. Any active neighbor whose coordinates "
            "fall in the 250 nm disk counts as interference. Yield combines independent estimates "
            "of qualifying and clear-disk probabilities. Each check requires agreement within "
            "five null standard errors. Checks share component estimates."
        ),
        "interpretation": (
            "This tests the analytic implementation under its own stationary kinetics, "
            "uniform spatial density and constant photon rate assumptions. Testing dye selection "
            "in cells requires matched experimental data."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dyes", type=Path, default=DATA / "dyes.json")
    parser.add_argument("--out", type=Path, default=DATA / "gpu-check.json")
    parser.add_argument("--device", type=int, default=0)
    args = parser.parse_args()
    dyes = json.loads(args.dyes.read_text(encoding="utf-8"))["dyes"]
    result = verify(dyes, args.device)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f'{result["device"]}: {result["check_count"]} checks; passed={result["passed"]}; '
          f'{result["wall_seconds"]} seconds, {result["peak_memory_mib"]} MiB peak tensors')
    if not result["passed"]:
        raise SystemExit("One or more Monte Carlo checks failed.")


if __name__ == "__main__":
    main()
