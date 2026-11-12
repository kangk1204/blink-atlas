# Reproduce or import a pilot

From the root of the source checkout, rebuild the S21 example with:

```bash
python -m analysis.s21 --cache .cache/s21 --docs docs
```

The fetch reads indexed byte ranges from [Zenodo record 16941576](https://zenodo.org/records/16941576). It retrieves the complete 9.46 MB Picasso table, the 25 kB YAML and a 32 MiB uncompressed prefix of the first OME-TIFF. CRCs verify complete ZIP members. SHA-256 pins all three local inputs, including the TIFF prefix. The full ZIP is 7.36 GB. A rebuild may reuse a verified cache with `--reuse`. The reference build used Python 3.11, NumPy 1.26.4, h5py 3.12.1, PyYAML 6.0.2, tifffile 2024.9.20 and Pillow 12.3.0.

The S21 preview covers the first 40 frames, with matching author fits. Pixel scale is 106 nm from the Picasso YAML. Exposure is 100.014 ms and the average frame interval in the clip is 100.042 ms from the TIFF timestamps. The full-field and cropped displays use separate fixed contrast ranges recorded in the JSON. The 224 × 224 pixel crop uses the region with the most localizations in the pilot. It is shown at 2×, and the scale bar is 2 µm. Data are from Walterspiel et al., CC BY 4.0. [Original experiment](https://doi.org/10.1002/anie.202424955)

## Process your experiment

CSV accepts `frame,x_nm,y_nm,photons` plus any additional columns. Use the upper-left outer edge of the detector as the coordinate origin. A frame ID of zero is the first acquired frame. For one-based exports add `--frame-base 1`. Picasso input accepts its structured `/locs` HDF5 dataset and matching multi-document YAML, or embedded JSON metadata. The importer converts pixel-center coordinates to edge-relative nanometers. YAML is detected by basename or can be supplied with `--yaml`.

Create a JSON file containing the actual acquisition values. The numbers shown here illustrate a fixed, 100 × 100 pixel acquisition at 100 nm/px; replace them with the values for your own sample and camera. The total number of acquired frames includes frames where no spots were localized. Use an ROI around a representative labeled region for an estimate of local density.

```json
{
  "sample_state": "fixed", "cell_line": "COS-7", "target": "H2B",
  "tag": "HaloTag", "reference_dye": "JF630b", "excitation_nm": 647,
  "exposure_ms": 25, "frame_interval_ms": 27,
  "width_nm": 10000, "height_nm": 10000, "total_frames": 1000
}
```

```bash
python -m analysis.pilot locs.csv --metadata acquisition.json --output my-pilot.json --min-photons 100
python -m analysis.pilot recording_locs.hdf5 --metadata acquisition.json --output my-pilot.json
```

For Picasso, the YAML or embedded metadata can supply image dimensions, `Pixelsize` in nm, frame count and exposure. For a raw HDF5 table without calibration, provide all the CSV acquisition values plus `pixel_size_nm`. The known excitation, specimen, target, reference dye and camera timing define a compatible fixed-cell pilot. Leave genuinely unknown conditions unspecified for quality inspection. The browser checks the conditions before using a summary to update fixed-cell predictions.

Optional `roi_nm` in the settings is a rectangle `[left, top, right, bottom]` relative to the full detector. A single Picasso ROI is honored automatically. For several localization ROIs, choose one rectangle lying inside a measured region. The spatial bins are complete 250 nm square tiles by default. Cells at incomplete right and bottom edges are omitted from the tile observations; all retained localizations in the selected ROI count toward the temporal forecast. Use the same photon threshold as the ranking model. Local processing writes only acquisition parameters and aggregate counts to JSON.

## Forecast check

The first 200 frames fit both candidate forecasts. The constant forecast uses their mean. The trend forecast divides those frames into four consecutive 50-frame blocks, computes the least-squares slope of the log count rates with a 0.5-count continuity correction, shrinks the slope by half, normalizes to the observed pilot mean and exponentially damps future log-rate change over 200 frames. Its parameters are fixed before evaluating frames 201 through 1,000. The score is the mean absolute error of counts in eight disjoint 100-frame blocks. A shorter input yields the available complete blocks; data that end during calibration can still yield a pilot-only summary.

The public S21 case keeps the authors' nonnegative photon fits. Fifteen of 214,794 source entries had negative fitted amplitudes, including zero in calibration and six in evaluation. The 12,513 pilot localizations yielded a baseline of 6,256.5 per 100 frames. The held-out MAEs were 765.25 for the constant baseline and 456.32 for the damped trend. The last blocks decayed faster than the forecast. This is a within-recording diagnostic in live U2OS cells with 100 ms exposure. It supplies no measured comparison of candidate dyes in fixed COS-7.
