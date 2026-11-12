"""Regression checks for coordinate matching, uncertainty and score denominators."""
import io
import itertools
import json
import math
from pathlib import Path
import unittest
from zipfile import ZipFile

import numpy as np

try:
    from analysis import epfl_precision as epfl
except ModuleNotFoundError as exc:
    if exc.name != "scipy":
        raise
    epfl = None


@unittest.skipUnless(epfl, "Install the optional EPFL SciPy environment to run these numerical checks")
class GeometryTests(unittest.TestCase):
    def test_lexicographic_matching_survives_greedy_trap(self):
        # Closest-first lets prediction 0 take the only in-range truth of prediction 1.
        pred = np.array([[0.5, 0.0], [-1.0, 0.0]])
        truth = np.array([[0.0, 0.0], [2.0, 0.0]])
        i, j, ambiguous = epfl.one_to_one(pred, truth, 1.5)
        self.assertEqual(list(zip(i, j)), [(0, 1), (1, 0)])
        self.assertEqual(ambiguous, 1)
        self.assertEqual(len(epfl.one_to_one([], truth, 1.5)[0]), 0)

    def test_hungarian_result_against_brute_force(self):
        gen = np.random.default_rng(11)
        for _ in range(30):
            a, b = gen.uniform(0, 6, size=(2, 4, 2))
            i, j, _ = epfl.one_to_one(a, b, 2)
            best = (0, 0.0)
            for n in range(1, 5):
                for left in itertools.combinations(range(4), n):
                    for right in itertools.permutations(range(4), n):
                        ds = np.linalg.norm(a[list(left)] - b[list(right)], axis=1)
                        if np.all(ds <= 2):
                            best = max(best, (n, -float(ds.sum())))
            self.assertEqual(len(i), best[0])
            self.assertAlmostEqual(np.linalg.norm(a[i] - b[j], axis=1).sum(), -best[1])

    def test_camera_variance_independent_finite_differences(self):
        # Independently differentiate the pixel means numerically. EPFL's
        # photoelectron-equivalent noise includes read and multiplication terms.
        theta = np.array([3.11, 3.38, 1200., 40., 1.1])  # x,y,N,bg,sigma
        grid = np.arange(7)

        def image(t):
            x, y, photons, bg, width = t
            cdf = np.vectorize(lambda v: math.erf(v/math.sqrt(2)))
            ax = (cdf((grid+.5-x)/width)-cdf((grid-.5-x)/width))/2
            ay = (cdf((grid+.5-y)/width)-cdf((grid-.5-y)/width))/2
            return (photons * np.outer(ay, ax) + bg).ravel()

        jac = []
        for k, step in enumerate((1e-4, 1e-4, .01, 1e-4, 1e-4)):
            diff = np.zeros(5)
            diff[k] = step
            jac.append((image(theta+diff)-image(theta-diff))/(2*step))
        jac = np.stack(jac, axis=-1)
        weights = 1/(25 + 1.41**2 * image(theta))
        exact = np.linalg.inv(jac.T @ (weights[:, None]*jac))[:2, :2]*100**2
        # MLE theta has row first, then column.
        model = epfl.camera_covariance(np.array([theta[[1, 0, 2, 3, 4]]]))[0]
        np.testing.assert_allclose(model, exact, rtol=2e-6, atol=1e-6)

    def test_predicted_uncertainty_filter_does_not_imply_accuracy(self):
        truth = [np.array([[0., 0.], [1000., 0.]])]
        fits = np.array([[25., 0.], [1000., 0.], [5000., 0.]])
        cov = np.tile(np.eye(2)*25., (3, 1, 1))
        report = epfl.evaluate(truth, np.array([1, 1, 1]), fits, cov, np.full((3, 2), 25.))
        self.assertAlmostEqual(report["fdr"], 1/3)
        self.assertEqual(report["recall"], 1)
        counts = report["precision20"]
        self.assertEqual(counts["passing_detections"], 3)
        self.assertEqual(counts["passing_matched"], 2)
        self.assertEqual(counts["matched_actual_radial_le20"], 1)
        self.assertEqual(counts["fraction_actual_radial_le20_given_matched_pass"], .5)

    def test_frame_ids_are_not_lexically_reindexed(self):
        buffer = io.BytesIO()
        with ZipFile(buffer, "w") as z:
            z.writestr("sequence/00011.tif", b"fake")
            z.writestr("sequence/00001.tif", b"fake")
            z.writestr("oracle/00001.tif", b"other")
        with ZipFile(buffer) as z:
            self.assertEqual(epfl.frame_index(z, r"sequence/(\d+)\.tif"),
                             {1:"sequence/00001.tif", 11:"sequence/00011.tif"})


class PublishedOutputTests(unittest.TestCase):
    def test_all_archived_frames_and_counts(self):
        path = Path(__file__).resolve().parents[1] / "docs/data/epfl-precision.json"
        result = json.loads(path.read_text(encoding="utf-8"))
        self.assertEqual([r["frames"] for r in result["datasets"]], [12000, 361])
        for r in result["datasets"]:
            self.assertEqual(r["truth_count"], 81049)
            self.assertEqual(r["matches"] + r["false_positives"], r["detection_count"])
            self.assertEqual(r["matches"] + r["false_negatives"], r["truth_count"])
            self.assertEqual(r["candidates"] - r["rejected_fits"], r["detection_count"])
            self.assertEqual(r["matches"], r["radius_sensitivity"][2]["matches"])
            self.assertLess(r["median_radial_error_nm"], 150.)
            self.assertAlmostEqual(r["psf"]["measured_x_fwhm_nm"], 258.21, delta=1.0)
            self.assertEqual(sorted(x["matches"] for x in r["radius_sensitivity"]),
                             [x["matches"] for x in r["radius_sensitivity"]])


if __name__ == "__main__":
    unittest.main()
