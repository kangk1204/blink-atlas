import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from analysis.pilot import evaluate, summarize
from analysis.s21 import SELECTED_SHA256


class PilotTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.path = Path(self.directory.name) / "counts.csv"
        self.field = dict(width_nm=800, height_nm=500, exposure_ms=25,
                          frame_interval_ms=27, total_frames=1000, sample_state="fixed",
                          cell_line="COS-7", target="H2B", tag="HaloTag",
                          excitation_nm=647, reference_dye="JF630b")

    def write_rows(self, rows):
        self.path.write_text("frame,x_nm,y_nm,photons\n" + "\n".join(rows), encoding="utf-8")

    def test_zero_frames_and_incomplete_edge(self):
        self.write_rows(["0,20,20,125", "0,90,30,110", "0,300,30,101",
                         "199,760,10,100", "199,200,200,99", "200,10,10,120",
                         "999,10,10,-2"])
        result = summarize(self.path, self.field)
        stats = result["pilot"]
        self.assertEqual(stats["localizations"], 4)
        self.assertEqual(stats["frames"], 200)
        tiles = stats["tiles"]
        self.assertEqual(tiles["observations"], 200 * 3 * 2)
        self.assertEqual(tiles["empty"], 1198)
        self.assertEqual(tiles["single"], 1)
        self.assertEqual(tiles["multiple"], 1)
        self.assertEqual(tiles["localizations"], 3)
        self.assertEqual(result["evaluation"]["blocks"][0]["observed"], 1)
        self.assertEqual(result["quality_control"]["negative_photon_fits_in_evaluation_window"], 1)

    def test_explicit_roi_with_half_open_edges(self):
        self.write_rows(["0,249,200,120", "0,250,200,120", "0,749.9,0,125", "0,750,0,125"])
        result = summarize(self.path, dict(self.field, roi_nm=[250, 0, 750, 500]))
        self.assertEqual(result["pilot"]["tiles"]["tiles_per_frame"], 4)
        self.assertEqual(result["pilot"]["localizations"], 2)
        self.assertEqual(result["quality_control"]["outside_selected_roi_in_source"], 2)

    def test_incompatible_metadata_and_short_exposures(self):
        self.write_rows(["0,100,20,100"])
        with self.assertRaisesRegex(ValueError, "shorter"):
            summarize(self.path, dict(self.field, frame_interval_ms=1))
        self.write_rows(["3,0,0,120"])
        result = summarize(self.path, dict(self.field, total_frames=3), frame_base=1)
        self.assertEqual(result["pilot"]["localizations"], 1)
        self.assertIsNone(result["evaluation"])
        with self.assertRaisesRegex(ValueError, "outside"):
            summarize(self.path, dict(self.field, total_frames=3))

    def test_empty_table_and_invalid_entries(self):
        self.write_rows([])
        result = summarize(self.path, self.field)
        self.assertEqual(result["pilot"]["tiles"]["empty"], result["pilot"]["tiles"]["observations"])
        self.assertEqual(result["evaluation"]["model_mae"], 0)
        self.assertFalse(result["evaluation"]["model_improved"])
        for row in ["0.5,0,0,10", "1,-1,0,10", "0,0,0,nan", "x,0,0,10"]:
            with self.subTest(row=row), self.assertRaises(ValueError):
                self.write_rows([row])
                summarize(self.path, self.field)

    def test_predictions_ignore_held_out_observations(self):
        pilot = np.array([10] * 50 + [12] * 50 + [8] * 50 + [7] * 50)
        a = evaluate(np.r_[pilot, np.zeros(800, dtype=int)], len(pilot))
        b = evaluate(np.r_[pilot, np.full(800, 100)], len(pilot))
        self.assertEqual([g["prediction"] for g in a["blocks"]], [g["prediction"] for g in b["blocks"]])
        self.assertEqual(a["blocks"][0]["baseline"], b["blocks"][0]["baseline"])
        self.assertNotEqual(a["model_mae"], b["model_mae"])

    def test_picasso_centers_match_normalized_csv(self):
        try:
            import h5py
        except ModuleNotFoundError:
            self.skipTest("h5py required for the HDF5 comparison")
        path = self.path.with_suffix(".hdf5")
        dtype = [("frame", "i4"), ("x", "f4"), ("y", "f4"), ("photons", "f4")]
        data = np.array([(0, 1, 1, 120), (5, 2, 2, 130)], dtype=dtype)
        with h5py.File(path, "w") as hf:
            hf.create_dataset("locs", data=data)
            hf.create_dataset("metadata", data=json.dumps([{"Width": 8, "Height": 5, "Frames": 1000}, {"Pixelsize": 100}]))
        other = summarize(path, self.field)
        self.write_rows(["0,150,150,120", "5,250,250,130"])
        csv = summarize(self.path, self.field)
        self.assertEqual(other["pilot"], csv["pilot"])
        path.with_suffix(".yaml").write_text("---\nWidth: 8\nHeight: 5\nFrames: 1000\n---\nPixelsize: 100\nROI: [[0, 0], [3, 3]]\n", encoding="utf-8")
        cropped = summarize(path, self.field)
        self.assertEqual(cropped["acquisition"]["camera_width_nm"], 800)
        self.assertEqual(cropped["acquisition"]["width_nm"], 300)
        with self.assertRaisesRegex(ValueError, "disagrees"):
            summarize(path, dict(self.field, width_nm=500))

    def test_published_case(self):
        root = Path(__file__).resolve().parents[1]
        result = json.loads((root / "docs/data/pilot.json").read_text(encoding="utf-8"))
        self.assertEqual(result["acquisition"]["sample_state"], "live")
        self.assertEqual(result["pilot"]["localizations"], 12513)
        self.assertEqual(result["pilot"]["four_block_counts"], [3287, 3072, 3093, 3061])
        self.assertEqual([b["observed"] for b in result["evaluation"]["blocks"]],
                         [6103, 5986, 5873, 5592, 5550, 4968, 5116, 4742])
        self.assertAlmostEqual(result["evaluation"]["baseline_mae"], 765.25)
        self.assertAlmostEqual(result["evaluation"]["model_mae"], 456.3202121513, places=5)
        for key, expected in SELECTED_SHA256.items():
            self.assertEqual(result["source"]["members"][key]["sha256"], expected)


if __name__ == "__main__":
    unittest.main()
