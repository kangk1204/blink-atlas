"""Numerical tests for reconstruction from deposited localizations.

The optional raw reanalysis environment includes SciPy, h5py, Numba and lmfit.
"""

import itertools
import hashlib
import json
import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

try:
    import lmfit
    from analysis import palette_raw as raw
    from analysis import palette_raw_run as runner
except ImportError:
    raw = runner = None


@unittest.skipIf(
    raw is None, "Install the optional DNA-PAINT raw analysis dependencies"
)
class PaletteRawTests(unittest.TestCase):
    def test_summary_is_identical_in_worker_and_resume_orders(self):
        rows = [
            {
                "system": "npc",
                "dye": "cy3b",
                "repeat": i,
                "principal": True,
                "files": [],
                "metrics": {"dropPct": v},
            }
            for i, v in enumerate([1e16, -1e16, 1.0], 1)
        ]
        self.assertEqual(
            runner.make_summary(rows, {}), runner.make_summary(rows[::-1], {})
        )

    def test_public_manifest_and_staging_index_agree(self):
        entries = [
            {"path": f"39.dye_analysis/fig_2/cy3b/n2/{name}.hdf5", "bytes": 50 + i}
            for i, name in enumerate(("full_field", "binding_sbs"))
        ]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            file = root / "index.json"
            file.write_text(
                json.dumps(
                    {
                        "files": [
                            {**e, "crc32": "000000ff", "sha256": "a" * 64}
                            for e in entries
                        ]
                    }
                )
            )
            original = runner.read_cases(root, file)
            key = ("origami", "cy3b", 2)
            self.assertEqual(original[key]["target_sites"]["crc"], 255)
            full_version = runner.identity(key, original[key], root)
            file.write_text(json.dumps([{**e, "crc": 255} for e in entries]))
            quick = runner.read_cases(root, file)
            self.assertEqual(
                quick[key]["target_field"]["path"],
                original[key]["target_field"]["path"],
            )
            self.assertNotEqual(full_version, runner.identity(key, quick[key], root))

    def test_localization_and_metadata_digests_are_verified(self):
        import h5py

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "field.hdf5"
            with h5py.File(path, "w") as h:
                h["locs"] = np.zeros(3, dtype=[("frame", "u4")])
            sidecar = path.with_suffix(".yaml")
            sidecar.write_text("---\nFrames: 20000\n---\nPixelsize: 130\n")
            entry = {
                "path": path.name,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "metadataSha256": hashlib.sha256(sidecar.read_bytes()).hexdigest(),
            }
            info = runner.file_info(root, "target_field", entry)
            self.assertEqual(info["archiveSha256"], entry["sha256"])
            self.assertEqual(info["pixelSizeNm"], 130)
            with self.assertRaisesRegex(ValueError, "SHA-256"):
                runner.file_info(root, "target_field", {**entry, "sha256": "0" * 64})
            sidecar.write_text("---\nFrames: 19999\n")
            with self.assertRaisesRegex(ValueError, "Sidecar SHA-256"):
                runner.file_info(root, "target_field", entry)

    def test_consecutive_events_and_pixel_integral(self):
        cols = [
            ("frame", "i4"),
            ("group", "i4"),
            ("photons", "f8"),
            ("bg", "f8"),
            ("sx", "f8"),
            ("sy", "f8"),
        ]
        data = np.array(
            [
                (f, g, p, 10, 1, 1)
                for g, f, p in [
                    (0, 1, 1),
                    (0, 2, 100),
                    (0, 3, 2),
                    (0, 7, 3),
                    (0, 8, 4),
                    (1, 1, 5),
                    (1, 2, 200),
                    (1, 3, 300),
                    (1, 4, 6),
                ]
            ],
            dtype=cols,
        )
        metrics, counts = raw.picked_metrics(data[::-1], 10, False)
        self.assertEqual(counts["groups"], 2)
        self.assertEqual(counts["events"], 3)
        self.assertEqual(counts["centerFrames"], 3)
        self.assertAlmostEqual(metrics["centerPhotons"], 200.0)
        self.assertAlmostEqual(
            metrics["centerSbr"], 20 * math.erf(0.5 / math.sqrt(2)) ** 2
        )
        self.assertAlmostEqual(metrics["eventPhotons"], sum(data["photons"]) / 3)
        self.assertIsNone(metrics["destroyedPct"])
        self.assertAlmostEqual(metrics["dropPct"], 75.0)
        self.assertIsNone(metrics["darkTimeS"])

    def test_cdf_matches_direct_single_parameter_optimization(self):
        from scipy.optimize import minimize_scalar

        for vals in [[3, 3, 3, 3], [1, 1, 1, 2, 3, 4, 5], [1, 2, 5, 8, 35, 450]]:
            a = np.array(vals)
            target = np.linspace(0.0, 1.0, len(a))
            fit = minimize_scalar(
                lambda lt: np.sum((1 - np.exp(-a / math.exp(lt)) - target) ** 2),
                bounds=(-5, 10),
                method="bounded",
            )
            self.assertAlmostEqual(raw.fitted_tau(a), math.exp(fit.x), delta=2e-4)
        self.assertTrue(np.isnan(raw.fitted_tau(np.array([2]))))

    def test_sparse_dwell_fit_matches_released_lmfit_initialization(self):
        def cdf(t, tau):
            return 1.0 - np.exp(-(t / tau))

        for values in ((2, 2500), (0, 1, 10, 800, 2500), (8, 80, 100, 150, 200, 250)):
            times = np.asarray(values)
            ranks = np.arange(len(times)) / (len(times) - 1)
            expected = lmfit.Model(cdf).fit(ranks, t=times, tau=200).best_values["tau"]
            self.assertAlmostEqual(raw.fitted_tau(times), expected, places=8)

    def test_neighbor_histogram_counts_all_adjacent_pairs_and_respects_groups(self):
        rng = np.random.default_rng(481)
        n = 120
        dt = [("frame", "i4"), ("x", "f8"), ("y", "f8"), ("group", "i4")]
        locs = np.zeros(n, dtype=dt)
        locs["frame"] = rng.integers(0, 5, n)
        locs["x"], locs["y"] = rng.uniform(-2, 2, (2, n))
        locs["group"] = rng.integers(1, 4, n)
        expected = np.zeros(1000, dtype=int)
        for a, b in itertools.product(locs, repeat=2):
            if a["frame"] + 1 == b["frame"] and a["group"] == b["group"]:
                dist = np.hypot(a["x"] - b["x"], a["y"] - b["y"])
                if dist < 1:
                    expected[int(dist * 1000)] += 1
        np.testing.assert_array_equal(raw.next_frame_histogram(locs), expected)
        self.assertGreater(expected.sum(), 10)

    def test_picasso_historical_source_batch_remainder(self):
        # Original Picasso's 100 progress chunks omit 1..99 source candidates,
        # while searches can still reach remainder rows as next-frame targets.
        dt = [("frame", "i4"), ("x", "f4"), ("y", "f4")]
        locs = np.zeros(105, dtype=dt)
        locs["frame"] = [0] * 99 + [1] * 2 + [2] * 4
        a = raw.next_frame_histogram(locs)
        b = raw.next_frame_histogram(locs, historical_batches=True)
        self.assertEqual(a[0], 99 * 2 + 2 * 4)
        self.assertEqual(b[0], 99 * 2 + 1 * 4)

    def test_historical_nena_fit_recovers_source_model(self):
        from scipy.special import i0

        d = (np.arange(1000) + 0.5) * 0.001
        a, s, ac, dc, sc = 380.0, 0.047, 230.0, 0.075, 0.11
        y = a * (d / s**2) * np.exp(-0.5 * d**2 / s**2) + ac * (d / sc**2) * np.exp(
            -0.5 * (d**2 + dc**2) / sc**2
        ) * i0(d * dc / sc)
        estimate, diag = raw.nena_fit(np.rint(y).astype(int), 0.043)
        self.assertAlmostEqual(estimate, 130 * s, delta=0.02)
        self.assertLess(diag["relativeResidualL2"], 0.002)

    def test_specificity_equal_areas_cancel_within_corresponding_rounds(self):
        roles = (
            "target_npc",
            "target_cytoplasm",
            "reference_npc",
            "reference_cytoplasm",
        )
        infos = [
            {"role": r, "rows": n, "nonemptyPicks": k}
            for r, n, k in zip(roles, [1000, 400, 600, 400], [4, 2, 3, 4])
        ]
        value, _ = runner.specificity(infos)
        self.assertAlmostEqual(value, (250 / 200) / (200 / 100))
        for item, diameter in zip(infos, [1.0, 4.0, 1.0, 4.0]):
            item["pickDiameterPixels"] = diameter
        self.assertAlmostEqual(value, runner.specificity(infos)[0])
        self.assertIsNone(runner.specificity(infos[:2])[0])


if __name__ == "__main__":
    unittest.main()
