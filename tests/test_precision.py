import copy
import json
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np

from analysis.precision import compare, field, validate_field


ROOT = Path(__file__).resolve().parents[1]


class PrecisionTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.path = Path(self.directory.name) / "locs.csv"
        self.meta = dict(width_nm=2000, height_nm=1000, roi_nm="full_frame", total_frames=10, exposure_ms=25,
                         frame_interval_ms=27, sample_state="fixed", cell_line="COS-7", target="H2B",tag="HaloTag",
                         excitation_nm=647, reference_dye="JF639b", prep_id="P1", field_id="A", illumination_state="readout",
                         density_labels_per_um2=100, density_origin="measured",density_calibration_id="independent count",
                         acquisition_id="M1",analysis_id="fitter 1.0",roi_protocol_id="preset rectangle")

    def make(self,rows,meta=None):
        self.path.write_text("frame,x_nm,y_nm,sigma_x_nm,sigma_y_nm\n"+"\n".join(rows),encoding="utf-8")
        return field(self.path, {**self.meta,**(meta or {})})

    def test_threshold_invalid_missing_and_geometry(self):
        f = self.make(["0,0,0,20,20","1,13,12,20.0001,1","2,12,14,,6","3,12,14,0,2",
                       "4,13,12,4,inf","5,800,12,2,2","6,-8,20,2,2","9,1,1,2,2"],
                      {"roi_nm":[0,0,500,1000],"window_frames":9})
        c=f["counts"]
        self.assertEqual((c["qualified"],c["valid_precision"],c["missing_precision"],c["invalid_precision"]),(1,2,1,2))
        self.assertEqual((c["outside_roi"],c["invalid_geometry"],c["outside_window"]),(1,1,1))
        self.assertEqual(f["denominators"]["empty_frames"],4)
        self.assertAlmostEqual(f["denominators"]["elapsed_s"],.241)
        self.assertAlmostEqual(f["result"]["per_um2"],2)
        self.assertEqual(f["result"]["status"],"partial")
        validate_field(f)

    def test_valid_zero_empty_and_no_precision(self):
        for rows in [[],["0,0,0,21,21"]]:
            f=self.make(rows)
            self.assertEqual(f["result"]["status"],"measured")
            self.assertEqual(f["result"]["per_um2_s"],0)
            validate_field(f)
        bad=self.make(["0,0,0,1,"])
        self.assertEqual(bad["result"]["status"],"unavailable")
        self.assertIsNone(bad["counts"]["qualified"])
        self.path.write_text("frame,x_nm,y_nm\n")
        bad=field(self.path,self.meta)
        self.assertIsNone(bad["counts"]["qualified"])
        validate_field(bad)

    def test_scalar_requires_definition_kept_separate(self):
        self.path.write_text("frame,x_nm,y_nm,precision_nm\n1,0,0,20\n")
        with self.assertRaisesRegex(ValueError,"scalar-name"):
            field(self.path,self.meta,frame_base=1)
        s=field(self.path,self.meta,frame_base=1,scalar_name="radial RMS")
        self.assertEqual(s["counts"]["qualified"],1)
        a=self.make(["0,0,0,10,10"],{"reference_dye":"JF630b"})
        self.assertIn("Precision definition, units or threshold differ.",compare([s,a],["JF639b","JF630b"])["reasons"])

    def test_hdf5_calibration_and_reversed_norm(self):
        try:
            import h5py
        except ModuleNotFoundError:
            self.skipTest("h5py required for the Picasso calibration test")
        path=self.path.with_suffix(".hdf5")
        a=np.array([(0,1,2,.2,.1),(1,1,2,.3,.1)],dtype=[("frame","i4"),("x","f8"),("y","f8"),("lpx","f8"),("lpy","f8")])
        with h5py.File(path,"w") as h:
            h.create_dataset("locs",data=a)
            h.create_dataset("metadata",data=json.dumps([{"Width":20,"Height":10,"Frames":10},{"Pixelsize":100}]))
        f=field(path,self.meta)
        self.assertEqual(f["counts"]["qualified"],1)
        self.assertEqual(f["denominators"]["empty_frames"],8)
        self.assertEqual(f["acquisition"]["pixel_size_nm"],100)
        validate_field(f)
        with self.assertRaisesRegex(ValueError,"position unit"):
            field(path,dict(self.meta,position_unit="nm"))

    def make_comparison(self):
        a=self.make(["0,0,0,10,10","1,0,0,5,5"],{"width_nm":1000})
        b=self.make(["0,0,0,10,10"],{"reference_dye":"JF630b","field_id":"B","width_nm":2000})
        c=self.make(["0,0,0,10,10"],{"field_id":"C","width_nm":1000})
        return a,b,c

    def test_paired_preparations_and_abstention(self):
        a,b,c=self.make_comparison()
        dyes=["JF639b","JF630b"]
        matched=compare([a,b,c],dyes)
        row=matched["preparations"][0]
        self.assertEqual(matched["status"],"matched")
        self.assertAlmostEqual(row["numerator"]["per_um2"],1.5)
        self.assertAlmostEqual(row["denominator"]["per_um2"],0.5)
        self.assertAlmostEqual(row["ratio_per_um2_s"],3)
        self.assertEqual(matched["model_domain"],"reference_acquisition")
        for key,val in [("density_labels_per_um2",200),("prep_id","P2"),("frame_interval_ms",30),
                        ("density_origin","nominal"),("acquisition_id","another setup"),("sample_state","live")]:
            other=self.make(["0,0,0,10,10"],{"reference_dye":"JF630b",key:val})
            self.assertEqual(compare([a,other],dyes)["status"],"unavailable",key)
        self.assertEqual(compare([a,a,b],dyes)["status"],"unavailable")
        zero=self.make([],{"reference_dye":"JF630b"})
        row=compare([a,zero],dyes)["preparations"][0]
        self.assertIsNone(row["ratio_per_um2_s"])
        self.assertGreater(row["difference_per_um2_s"],0)
        a2=self.make(["0,0,0,2,2"],{"excitation_nm":640,"exposure_ms":20})
        b2=self.make(["0,0,0,2,2"],{"reference_dye":"JF630b","excitation_nm":640,"exposure_ms":20})
        self.assertEqual(compare([a2,b2],dyes)["model_domain"],"outside_reference_acquisition")
        corrupt=copy.deepcopy(a)
        corrupt["counts"]["qualified"]=900
        with self.assertRaisesRegex(ValueError,"qualified counts"):
            validate_field(corrupt)

    def test_browser_python_parity(self):
        node=shutil.which("node")
        if not node:
            self.skipTest("Node required for the browser parity harness")
        a,b,c=self.make_comparison()
        q=self.make([],{"reference_dye":"JF630b"})
        bad=self.make(["0,0,0,1,2"],{"prep_id":"P2","reference_dye":"JF630b"})
        d=self.make(["0,250,20,11,20"],{"prep_id":"P2","field_id":"A"})
        ambiguous=self.make(["0,0,0,1,"])
        cases=[[a,b,c],[a,q],[b],[a,bad],[a,b,bad,d],[a,a,b],[a,b,ambiguous],[]]
        dyes=["JF639b","JF630b"]
        payload=[{"fields":x,"dyes":dyes,"expected":compare(x,dyes)} for x in cases]
        script=ROOT/"tests/precision-parity.mjs"
        p=subprocess.run([node,str(script)],input=json.dumps(payload),text=True,capture_output=True,cwd=ROOT)
        self.assertEqual(p.returncode,0,p.stderr)

    def test_public_data_and_scope(self):
        j=json.loads((ROOT/"docs/data/precision.json").read_text(encoding="utf-8"))
        self.assertEqual(j["holland"]["qualifying_count_bounds"],[928067,938285])
        self.assertEqual(j["s21"]["counts"]["qualified"],110694)
        self.assertEqual(j["independent_check"]["qualified"],110694)
        self.assertEqual(len(j["manifest"]),21)
        self.assertEqual(len(j["auxiliary"]),19)
        self.assertEqual(sum(m["status"]=="quarantined" for m in j["manifest"]),1)
        for record in [j["s21"],*j["auxiliary"]]:
            validate_field(record)
        self.assertEqual(j["matched_fixed_comparison"]["status"],"unavailable")
