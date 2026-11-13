# Blink Atlas

**Design a free dye for titration.** Probe design now accepts local structure-result JSON and includes two fluorinated Si-rhodamine candidates, downloadable SDFs, held-out study errors and the matched JF626b control. [Open the workbench](https://kangk1204.github.io/blink-atlas/#probe-design) · [Reproduce the structure benchmark](analysis/structure-pk.md) · [Paired titration protocol](validation/structure_pk/titration.md)

An interactive research demo for selecting the first two spontaneously blinking HaloTag dyes to test in a fixed-cell SMLM experiment. Begin with fixed COS-7 cells expressing H2B–HaloTag, 647 nm excitation, a 25 ms exposure and a 27 ms camera cycle. The nine dyes, measured properties, site assets and scientific checks load without an account. [Open the comparison](https://kangk1204.github.io/blink-atlas/) · [Watch the two-minute walkthrough](https://kangk1204.github.io/blink-atlas/media/walkthrough.mp4) · [One-page collaboration proposal](https://kangk1204.github.io/blink-atlas/media/validation-proposal.pdf)

## Use the comparison

The time-ratio audit reconstructs nine source rows and compares three models under all 46,656 whole-dye pairing scenarios. It also records the screened Uno and Chen dyes, including each reason a score could not be calculated. [Methods, results and reproduction](analysis/ratio-audit.md) · [Open the validation screen](https://kangk1204.github.io/blink-atlas/#unseen-dyes)

The probe workbench links the nine free-dye structures and pH titrations to the corresponding HaloTag-ligand measurements. Compare live and fixed cellular SOFI, then lock and score a new analogue’s predicted H2B ON:OFF time ratio. [Open Probe design](https://kangk1204.github.io/blink-atlas/#probe-design) · [Data and prospective protocol](analysis/probe-design.md)

The unseen-dye benchmark adds 25 measurements across 11 additional parent dyes. On four eligible external IgG-conjugated dyes, the frozen local model had mean absolute log₁₀ error 0.293, nearest pKa 0.307 and the original regression 0.208. Local interpolation beat nearest for two dyes. Independent validation on new dyes in the fixed H2B assay remains open. [Explore all predictions and assay changes](https://kangk1204.github.io/blink-atlas/#unseen-dyes) · [Frozen benchmark, results and next experiments](analysis/analogue-benchmark.md)

The first prospective candidate is Cal664sb. Its published free-parent pKcycl of 3.2 yields a conditional local forecast of 0.000580, nearest-pKa forecast of 0.000484 and regression forecast of 0.000430. Matched cellular observations are pending. A laboratory can enroll additional blinded parents and keep all unpublished files local. [Dated forecasts](https://github.com/kangk1204/blink-atlas/blob/4c68d5a39e2559d05e4b8b8144f4fba5bf2341aa/validation/h2b_blind/registration.md) · [Laboratory files and scoring instructions](validation/h2b_blind/README.md)

1. Set a plausible *local* label density interval in labels/µm² and the minimum expected photon count. The slider shifts the interval along a logarithmic scale. Broad range and Dense range are illustrative assumptions.
2. Inspect the two reference curves or display all nine. The table shows nominal isolated detections per µm²/s at the geometric midpoint, bounds over the selected interval and sensitivity assumptions, and the conditional probability of an active neighbor within 250 nm. Estimates refer to detections, which can include repeated emission from the same molecule.
3. A pair is recommended only if it outranks every alternative for the full density interval at all tested kinetic settings. At 100–1,000 labels/µm² and 100 expected photons, the pair is unresolved. At 10,000–10,500 with the same threshold the model selects JF630b and JF629b. These are input-dependent predictions. Cellular performance requires an independent matched experiment.
4. To condition the estimate on a known reference dye, process a representative ROI locally and import its summary JSON. The browser checks the specimen, protocol, photon filter and tile counts. A matched informative pilot updates the prior density interval, and its result can be compared with the publication-only estimate. An incompatible pilot withholds the recommendation. Export analysis saves the settings, pilot summary and computed results. Importing a saved analysis recomputes it. The browser processes the local file in memory.

Under **Candidate yields**, choose **Validation panel · 3 ligands** to score the three proposed study reagents only. This recomputes the uncertainty check with JF630b, JF635b and JF639b as the eligible candidates. Set `"candidatePanel":"validation"` in the same settings JSON for CLI use.

The reference-pair images are actual first-frame crops from the published supplementary movie. They retain the authors' display scaling. The source workbook, recorded clip, acquisition metadata and source cells are accessible under **Methods & data**. [Holland et al.](https://doi.org/10.1038/s41592-026-03062-5) · [Crop provenance](docs/data/reference-frames.json)

## Public results and interpretation

| Check | Result | Supported interpretation |
|---|---|---|
| Fixed-cell source reproduction | All nine ligands. Mean ON times for JF630b and JF635b are 10.14 and 14.53 ms, with 67.2 and 51.5 photons/ms, respectively. | Three cellular imaging experiments per dye. Off-times are upper-bound estimates under the measurement protocol. Detailed units, SDs, replicate values and workbook cells are included. |
| S21 before-illumination forecast | 12,513 retained author localizations in the first 200 frames. In eight consecutive 100-frame blocks, the damped trend scored 456.32 mean absolute count error and the constant mean scored 765.25. | One live U2OS JF630b recording with approximately 100 ms exposures. It tests within-recording forecasting. The holdout contains 800 frames. |
| EPFL overlap geometry | At 250 nm, 9.482% of the 100,000 published emitter-frame coordinates have another emitter within the radius. Whole-field homogeneous density predicts 0.553%. Exhaustive and spatial-bin pair checks agree. | Tubulin clustering materially changes the overlap risk. A global density estimate can greatly understate the risk in dense structures. |
| Independent CUDA check | All 12 checks of finite-exposure and neighbor probabilities passed a prespecified five-standard-error tolerance. | The stochastic check verifies the analytical implementation under the simulated assumptions. |

The camera model uses a first ON segment, exponential dwell times, a hard threshold on expected photons and a spatial Poisson approximation. Individual optically separable events at low SNR, heterogeneous labeling, detector noise, photobleaching and the biological number of distinct molecules can differ from the proxy. The pilot forecast and fixed-cell recommendations have separate inputs and evaluation domains. [Model derivation](analysis/photophysics.md) · [External pilot record](https://zenodo.org/records/16941576) · [EPFL benchmark](docs/data/epfl.json) · [CUDA results](docs/data/gpu-check.json)

## Reproduce the analysis

### Precision evidence

The precision view reports actual fitted uncertainties. Holland's fixed-cell JF630b histogram puts the passing fraction at 88.51–89.48%. The available bins cannot resolve the 20 nm boundary further. Public matched low-density JF639b fields with the same endpoint are still needed. For the live JF630b S21 control, 110,694 of 214,794 fits pass the two-axis cutoff. Across the 20 eligible live recordings, the pilot trend forecast improves count MAE in 12 and worsens it in eight. Illumination segments can share the same specimen. [Precision evidence](https://kangk1204.github.io/blink-atlas/#precision) · [Source search](validation/h2b_precision/README.md)

One local command imports calibrated CSV or Picasso per-event uncertainties. A second compares two dyes within matched preparations and preserves the area and time denominators. The browser can open either output or a group of field summaries. An absent precision column, a partial field, an unmatched acquisition and a measured count of zero have separate meanings. [Field metadata and commands](analysis/precision-input.md)

The independent EPFL test processes complete sparse and dense sensor sequences. The single-emitter fitter recovered 91.31% and 7.37% of true emitter positions, respectively, at a 150 nm matching radius. The nominal 95% uncertainty ellipses contained 80.85% and 3.28% of the matched truth positions. [EPFL fitting and matching](analysis/epfl-precision.md) · [Machine-readable outcomes](docs/data/epfl-precision.json)

### DNA-PAINT palette

The palette compares 18 imager dyes within three excitation channels on DNA origami and nuclear pore complexes. Each measurement links to a workbook cell and each recomputed value to a set of deposited Picasso localization tables. In the browser, choose a channel, specimen, metric and dye. Open **Source agreement** to see where the reanalysis differs. The NPC evaluation selects two dyes using repeats 1 and 2 and measures their mean NeNA on repeat 3. Its fixed baseline uses the historical NPC ranks. The test unit is a dye acquisition. [Open the palette](https://kangk1204.github.io/blink-atlas/#dna-paint) · [Frozen decision protocol](analysis/palette-protocol.md) · [Steen et al.](https://doi.org/10.1038/s41592-024-02374-8)

The complete archive has 1,167 files, including 575 HDF5 tables. The two 18-dye panels use 435 of those tables and their 435 YAML sidecars. The analysis covers 108 principal dye/specimen/repeat groups plus an additional Cy3B group that duplicates a reference acquisition. All 653 reported repeat cells have corresponding raw-table estimates, and all 432 workbook means and SDs reproduce from the workbook repeats. Some HDF5 estimates differ from the published cells. For JF585 NPC repeat 3, the deposited target and historical Picasso implementation give 4.3674 nm, and the workbook reports 19.9 nm. The cause is unresolved. [Archive inventory](docs/data/palette-manifest.json) · [Per-metric agreement](docs/data/palette-audit.json) · [Source identity and fit checks](analysis/palette-raw.md#two-nena-discrepancies)

On the tuning repeat, weights 0 and 0.5 tied. The prespecified rule chose 0, so the selected pairs matched the NPC-only baseline in all three channels. Relative to the best two dyes observed in each held-out channel, both policies averaged 0.3182 nm additional NeNA in the HDF5 reanalysis and 0.0633 nm in the separate workbook replay. The display withholds a new pair recommendation. [Held-out results and individual pairs](docs/data/palette-evaluation.json)

For a small replay, the two source workbooks and all computed JSON are included. Install `openpyxl==3.1.5` with Python 3.11 or later and run these commands at the repository root.

```sh
python -m analysis.palette_catalog
python -m analysis.palette_audit
python -m analysis.palette_evaluation
```

For a complete reanalysis, fetch the 38,156,495,377-byte archive to a POSIX data volume with room for both the RAR and its contents. Use Python 3.11 or later with `libarchive-c==5.3`, `rarfile==4.2`, a system libarchive with RAR5 support and h5py for inventory. The fetch command verifies Zenodo's MD5 and can resume from its range journal. The manifest records each extracted file's path, byte size, CRC32 and SHA-256, plus dataset shape and field names for HDF5. Then run the panel fitter with the environment in the raw analysis guide, followed by the three replay commands above. [Raw fitter, definitions and dependency versions](analysis/palette-raw.md) · [Source archive and license](https://zenodo.org/records/10960858)

```sh
python -m analysis.palette_archive fetch /data/blink-palette
python -m analysis.palette_archive extract /data/blink-palette
python -m analysis.palette_archive inventory /data/blink-palette --output docs/data/palette-manifest.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m analysis.palette_raw_run \
  --root /data/blink-palette/extracted --index docs/data/palette-manifest.json \
  --cache /data/blink-palette/case-cache --workers 4 --output docs/data/palette-raw.json
```

### Original comparison and recording

The static site is served from `docs/`. Start a local preview with `python -m http.server 8000 --directory docs`, then open `http://localhost:8000`. The published site uses the same files.

For analytical regeneration, use Python 3.11 and Node.js 22.7 or later. Install the CPU dependencies in a virtual environment.

```sh
python -m pip install -r requirements.txt
python -m analysis.dyes
python -m analysis.model
node analysis/predict.mjs
python -m analysis.s21 --cache .cache/s21 --docs docs
python -m analysis.epfl
python -m analysis.make_proposal
python -m unittest discover -s tests
node --test tests/model-pilot.mjs
```

The small publisher workbook is included with its SHA-256. S21 reproduction downloads the complete author HDF5 and YAML and a verified 32 MiB TIFF prefix using HTTP range requests. The public example includes 40 raw camera frames with their author-fitted localization overlays. The EPFL script fetches the 5.5 MB ground-truth archive and checks its hash. No training data or GPU is required to run the site. On a CUDA-enabled PyTorch installation, `python -m analysis.gpu_check` repeats the optional Monte Carlo audit. [File hashes and input formats](analysis/pilot-input.md)

For a CSV or Picasso table, prepare an acquisition metadata JSON, run the importer, and choose **Import pilot or analysis JSON** on the site. For CLI inference supply a JSON containing the acquisition settings, such as `{"densityMin":100,"densityMax":1000,"minPhotons":100}`.

```sh
python -m analysis.pilot locs.csv --metadata acquisition.json --min-photons 100 --output pilot.json
node analysis/predict.mjs settings.json pilot.json
```

The importer supports zero- or one-based frame numbering, explicitly calibrated pixel or nanometre coordinates and a representative ROI. Use actual acquired frame counts, including frames with no detections. [CSV, Picasso and SMAP export instructions](analysis/pilot-input.md)

### Media

- To rebuild the two published fixed-cell micrographs, run `python tools/reference_frames.py` with FFmpeg on PATH. The script retrieves the original supplementary movie and records the source digest and crop coordinates.
- To rebuild the English film, serve the site locally and use Playwright with an installed Chromium browser and Poppler's `pdftoppm` on PATH. Run `node tools/capture-tour.cjs http://localhost:8000 .qa/tour.webm chrome`, then `pwsh -File tools/narrate.ps1 -OutputDirectory .qa/voice` on Windows with an English SAPI voice. Finally run `python tools/build_video.py .qa/tour.webm .qa/voice`. The narrated MP4 has a selectable English subtitle track and a separate WebVTT track for the web player. The script and on-screen actions are in `tools/`.

## Prior work and proposed experiment

Fluorophore screening has already established the relevance of photon yield, duty cycle and photostability. SMIS simulates fluorophore photophysics and imaging conditions, and autonomous illumination control estimates active-emitter density from images. The question for this demo is whether a short pilot changes the choice of two dyes enough to improve the outcome on matched samples. [Dempsey et al.](https://doi.org/10.1038/nmeth.1768) · [Bittel et al.](https://doi.org/10.1038/srep29687) · [SMIS](https://doi.org/10.1038/s42003-023-04432-x) · [Autonomous illumination control](https://doi.org/10.1364/OE.26.030882)

Start with existing matched fixed-cell recordings or localization tables from at least two candidate dyes. The proposed follow-up uses three dyes, two labeling densities, three independent preparation days and three fields per condition, for 54 movies. Compare the model-selected pair against a fixed JF630b/JF635b pair using the same pilot-time budget. The primary outcome is localizations per µm² passing the prespecified ≤20 nm localization-uncertainty criterion at equal acquisition times. Aggregate fields within each preparation day. [Collaboration proposal](docs/media/validation-proposal.pdf)

## Credits

The dye measurements, workbook and reference movie crops are from Holland et al., *Nature Methods* (2026), under CC BY 4.0. The S21 files are from Walterspiel et al., Zenodo (2025), under CC BY 4.0. The site displays a cropped, contrast-scaled and annotated excerpt. The external experiment is described in *Angewandte Chemie International Edition*. EPFL's Tubulins II coordinates are fetched from the authors' site on demand. The bundled DM Sans and IBM Plex Mono fonts include their SIL Open Font License texts in `docs/assets/`. [Holland et al.](https://doi.org/10.1038/s41592-026-03062-5) · [S21 dataset](https://doi.org/10.5281/zenodo.16941576) · [S21 paper](https://doi.org/10.1002/anie.202424955) · [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
