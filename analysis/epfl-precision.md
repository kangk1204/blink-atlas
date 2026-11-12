# Image-to-precision benchmark

Blink Atlas processes all sensor TIFFs and all corresponding projected ground-truth coordinates from two EPFL 2013 training sets. Both simulations use bundled tubules, 100 nm pixels and 64 × 64 pixel frames. Each archive contains 81,049 ground-truth emitter-frame coordinates. The sparse sequence spans 12,000 frames. The dense sequence contains 361 matching TIFF/CSV frame IDs, despite the 168-frame headline on its web page. The six archive SHA-256 digests and original download locations are in the [results JSON](../docs/data/epfl-precision.json). [Sparse dataset](https://bigwww.epfl.ch/srm/dataset/challenge-2D-simulation/Bundled_Tubes_Long_Sequence/index.html) · [Dense dataset](https://bigwww.epfl.ch/srm/dataset/challenge-2D-simulation/Bundled_Tubes_High_Density/index.html) · [Challenge paper](https://doi.org/10.1038/nmeth.3442)

## Reproduce on a CPU

The full reference run used Linux, Python 3.11.16, NumPy 1.26.4, SciPy 1.13.1, Numba 0.61.0, llvmlite 0.44.0, tifffile 2025.6.11, PyYAML 6.0.2 and PicassoSR 0.7.5. Use an isolated environment for this optional analysis. Only the published numerical MLE module is imported from Picasso. Its GUI and visualization dependencies are unnecessary for the benchmark.

```bash
python3.11 -m venv .venv-epfl
. .venv-epfl/bin/activate
python -m pip install numpy==1.26.4 scipy==1.13.1 numba==0.61.0 llvmlite==0.44.0 tifffile==2025.6.11 PyYAML==6.0.2
python -m pip install --no-deps picassosr==0.7.5
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python -m analysis.epfl_precision --workers 8
python -m unittest discover -s tests -p 'test_epfl_precision.py' -v
```

The script downloads the archives once and checks their digests on each run. Use `--cache DIRECTORY` for another location or existing downloads and `--out PATH` for a separate result. The full versioned reference tables are under `datasets[].ground_truth_crowding`, `uncertainty_bins`, `radius_sensitivity` and `precision20`.

## Detection, matching and uncertainty

Both conditions use the same finder and fitter. The difference of two Gaussian filters with widths of 1 and 3 pixels emphasizes spots. The finder keeps positive 3 × 3 local maxima more than six robust noise units above the frame median. On the unfiltered image, it subtracts the 100 ADU camera baseline and passes each valid 7 × 7 pixel ROI to Picasso's circular, pixel-integrated Gaussian Poisson MLE with a constant local background. The result file lists the acceptance bounds, number of rejected fits and fits reaching the iteration limit. [Picasso implementation](https://github.com/jungmannlab/picasso/blob/v0.7.5/picasso/gaussmle.py) · [Picasso methods](https://doi.org/10.1038/nprot.2017.024)

Frame IDs come from the one-based filenames. Picasso returns pixel-center indices, with row first in its parameter array. EPFL coordinates use the field edge as origin. A shift of half a pixel in each axis and multiplication by 100 convert the recovered XY centers to nm. A 10 nm-sampled calibration PSF gave a focal-plane x FWHM of 258.30 nm for both datasets (documented FWHM 258.21 nm). The PSF serves as an independent spatial calibration check. The fitted width comes from each ROI. [EPFL optical and camera settings](https://bigwww.epfl.ch/srm/dataset/challenge-2D-simulation/Bundled_Tubes_Long_Sequence/index.html)

Frame-wise assignments maximize the number of one-to-one matches inside 150 nm, then minimize their total Euclidean XY distance. The result includes alternative matching radii of 50, 100 and 250 nm. Recall divides matched positions by all truth positions. FDR divides unmatched detections by all accepted fits. An additional control matches each set of detections to a truth frame half an acquisition away, preserving the common tubulin structure while breaking the temporal correspondence.

For each fitted spot the primary uncertainty estimate uses the inverse of the five-parameter information matrix, including fitted position, brightness, background and width. A working noise variance of 1.41² times the fitted mean plus 5² accounts for the benchmark camera's nominal multiplication and read noise. The benchmark retains the raw Picasso Poisson CRLB for comparison. The bivariate 95% coverage uses the full XY marginal covariance and a squared Mahalanobis radius of 5.991. This is a first-order, model-based precision estimate. Bias, unresolved emitters and model mismatch affect its calibration.

## Results

| Full-sequence metric | Sparse | Dense |
|---|---:|---:|
| Accepted fits / matches within 150 nm | 75,220 / 74,002 | 6,028 / 5,972 |
| Recall | 91.31% | 7.37% |
| FDR at 150 nm | 1.62% | 0.93% |
| Median / RMS radial error among matches | 7.76 / 28.74 nm | 42.89 / 60.08 nm |
| Empirical coverage of a nominal 95% ellipse | 80.85% | 3.28% |
| True radial error ≤20 nm among matched fits passing the 20 nm per-axis uncertainty filter | 89.21% | 18.55% |

The dense case has 224.5 truth positions per frame on average. A wrong frame still yields 4,648 assignments within 150 nm, compared with 5,972 for the correct frames. These control counts show how readily a fitted peak can fall near an unrelated emitter in a crowded field. Error and coverage statistics condition on the matched set, and the matching radius censors its error tail. This benchmark characterizes the single-emitter analysis of a simulated camera sequence. Independent JF dye rankings require matched cell experiments.
