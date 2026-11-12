# Measure localization precision

The **field** command counts localizations with finite, positive uncertainties of at most 20 nm on **both** lateral axes. Each accepted localization contributes one event; a fluorophore can contribute again in a subsequent frame. The endpoint measures fitted position uncertainty. Resolution of reconstructed structures is a separate measurement.

## Prepare one field

Accepted inputs are a Picasso HDF5 `/locs` table with `frame,x,y,lpx,lpy`, or a CSV with `frame,x_nm,y_nm,sigma_x_nm,sigma_y_nm`. Picasso reports lateral location and uncertainty in camera pixels. Its companion YAML or HDF5 metadata supplies the physical pixel size, or you can specify `pixel_size_nm`. Coordinates in a CSV are in nm from the outer left/top detector edges. Set `--frame-base 1` for one-based CSV exports. Read the [Picasso data format](https://picassosr.readthedocs.io/en/latest/files.html).

Create one acquisition metadata file per field. The numbers below are an illustrative input template. Replace the density with an **independently measured** label density and record its calibration protocol. A reagent concentration in nM does not measure labels per area.

```json
{
  "prep_id": "prep-01", "field_id": "well-A-field-03",
  "sample_state": "fixed", "cell_line": "COS-7", "target": "H2B",
  "tag": "HaloTag", "reference_dye": "JF639b",
  "density_labels_per_um2": 100, "density_origin": "measured",
  "density_calibration_id": "calibration-protocol-and-reference", "density_unit": "labels/um2",
  "illumination_state": "readout", "excitation_nm": 647,
  "exposure_ms": 25, "frame_interval_ms": 27,
  "acquisition_id": "scope-readout-calibration-v1", "analysis_id": "fit-filter-settings-v1",
  "roi_protocol_id": "preregistered-nuclear-rectangles-v1",
  "width_nm": 20000, "height_nm": 20000,
  "roi_nm": [1000, 2000, 8000, 9000],
  "total_frames": 5000, "start_frame": 0, "window_frames": 1000
}
```

```bash
python -m analysis.precision field positions.csv --metadata acquisition.json --output field.json
python -m analysis.precision field sample_locs.hdf5 --yaml sample_locs.yaml --metadata acquisition.json --output field.json
```

Use a rectangular ROI `[left,top,right,bottom]` or declare `"roi_nm":"full_frame"` for a full-field localization pass. The right and bottom edges are exclusive. Picasso localization ROIs constrain the usable region. A single recorded Picasso ROI can be used directly; choose a single rectangle inside the measured region when several exist. Include every acquired frame in the declared length, including frames without detections. `start_frame` is always zero-based. The elapsed acquisition window spans the first shutter opening to the last shutter closing, `((N−1) × frame_interval_ms + exposure_ms) / 1000`. Missing cadence permits an area-only summary for a multi-frame recording.

For a table with a named scalar uncertainty in nm, specify its definition explicitly:

```bash
python -m analysis.precision field scalar.csv --metadata acquisition.json --scalar-column precision_nm --scalar-name "Reported lateral RMS" --output scalar-field.json
```

Different scalar definitions and per-axis endpoints stay separate. Zero, negative and infinite uncertainties are invalid. Blank or NaN values are missing. Per-field output includes every excluded category, precision-qualified count, fraction of valid fits, events/µm², and events/(µm²·s). Valid uncertainties with no qualified fits yield a measured zero. A verified empty acquisition with the expected columns also yields zero events and an undefined fraction. A table lacking usable precision yields an unavailable result. Partial fields disclose the measured subset and coverage. No photon cutoff is added to the precision endpoint; record acquisition-side and processing filters in `analysis_id`.

## Compare fields

```bash
python -m analysis.precision compare day1-639-a.json day1-639-b.json day1-630-a.json --dyes JF639b JF630b --output matched.json
```

Field and comparison JSON files can be opened in the precision view of the static website. The browser can also compare multiple field summaries locally. No uploaded summary contains individual coordinates or source paths. A field summary includes the source-table hash and the calibration, ROI and frame window used to produce it.

- The first dye argument is the numerator. Pair within each `prep_id`. The comparison requires complete, positive, finite precision coverage, except for verified empty fields. A duplicate preparation/dye/field identity, an unpaired preparation, or a missing comparison condition suppresses the contrast.
- Match fixed COS-7 H2B–HaloTag, measured label density and calibration ID, illumination, excitation, exposure, cycle, frame window, acquisition/analysis/ROI protocol IDs, precision definition and units. The comparison uses matched density values with only numerical roundoff tolerance. Externally matched acquisitions can inform an observed comparison even if their acquisition parameters differ from the dye prediction model.
- Within each preparation and dye, pool counts over the sum of ROI areas and over the sum of ROI area × elapsed seconds. The experiment can therefore use different ROI areas. Compute the difference and, where the denominator is nonzero, the ratio between dye arms. Across preparations summarize the **mean of preparation-level differences**. Fields from the same preparation contribute to the same comparison.

## Reproduce public evidence

```bash
python -m analysis.precision public --cache .cache/precision --output docs/data/precision.json
python -m analysis.precision public --offline --cache .cache/precision --output precision-rebuilt.json
```

The acquisition code verifies published ZIP identities, member CRCs and an audited membership fingerprint. It selectively fetches fit tables, YAML acquisition metadata and compressed TIFF prefixes to determine camera cadence. Each eligible Walterspiel table uses the full field and all acquired frames; elapsed times are estimates from the mean of the first 39 timestamp intervals. The public manifest gives the precise member, hashes, measured frame interval and timing method.

- Holland et al. publish a **scalar, binned** fixed-cell H2B–HaloTag distribution for JF630b. The 19.75–20.25 nm bin straddles the threshold. Of 1,048,575 displayed histogram counts, 928,067 certainly qualify and 10,218 are ambiguous, yielding 88.51–89.48%. The figure does not supply a matched low-density JF639b field. [Source workbook](https://pmc-oa-opendata.s3.amazonaws.com/PMC13167466.1/41592_2026_3062_MOESM11_ESM.xlsx)
- Walterspiel et al. supply live-cell QC. An independent whole-table query and the chunked importer both find 110,694 qualifying fits of 214,794 in the S21 pre-illumination control. S21 used 106 nm/px. Its projected area is 4,044.96 µm² and the estimated time is 500.210 s. [Data and license](https://zenodo.org/records/16941576)
- Across the four archives there are 21 actual fit tables. One S20 before-illumination member encodes ordinary HaloTag in its filename but its parent archive specifies psHaloTag1a; that recording is quarantined. The other 20 records are identified by tag and illumination phase. For the separate within-movie count forecast, the damped model had a lower MAE in 12 records and a higher MAE in eight. Multiple phases of a single biological specimen can be represented, so the table describes recordings. [Experimental protocol](https://doi.org/10.1002/anie.202424955)
