# Photophysics and density inference

The catalog and score use the fixed COS-7 H2B–HaloTag experiments in [Holland et al., *Nature Methods* (2026)](https://doi.org/10.1038/s41592-026-03062-5). The source workbook under `docs/data` is the unmodified Extended Data Fig. 4 statistical source data, licensed CC BY 4.0. Its SHA-256, per-value sheet coordinates, optical settings and ligand mapping are in `dyes.json`. Panel h–i supplies the reported mean-based and median-based time ratios. The arithmetic mean of three measured means and the arithmetic mean of three measured medians are stored separately.

## Event calculation

Density ρ is the effective projected density of labeled molecules within a locally homogeneous region (labels/µm²). The optical assumption is a circular exclusion neighborhood of area A = π(0.25 µm)², measured by distance between emitter centers. Consider a stationary two-state molecule with mean bright dwell τ₁, mean dark dwell τ₀, and mean photon rate r while ON. T and Δ are exposure and frame interval in ms. Default T = 25 ms and Δ = 27 ms.

1. The equilibrium ON probability is p = τ₁/(τ₀ + τ₁). The probability that at least one ON interval intersects the exposure by time u is B(u) = p + (1 − p)(1 − exp(−u/τ₀)).
2. A burst needs at least a = h/r milliseconds of overlap to reach a threshold of h expected photons. The probability that its first intersecting bright segment qualifies is q = exp(−a/τ₁) B(T − a) for a ≤ T, and zero otherwise. The model neglects aggregation of multiple bursts from the same molecule in one exposure. Dark waits in the source data are much longer than the 25 ms exposure.
3. Any neighbor within the collision disc that is ON during the exposure can interfere. Under independent positions and blinking, the conditional isolation probability for the candidate is exp(−ρ A B(T)). The expected isolated yield per µm² per second is **ρ q exp(−ρ A B(T)) × 1000/Δ**. Per-frame observations can include the same emitter in consecutive frames. The model does not estimate the number of distinct molecules.

The photon cutoff uses expected counts. Photon shot noise, camera background, spatially varying labeling efficiency, PSF fitting, transport, bleaching over a long acquisition and structural clustering can change observed yields. Spot overlap is a geometric proxy. The EPFL check documents the effect of spatial clustering. The earlier 1%-on Bernoulli calculation asks for the fraction of all footprints containing at least two simultaneously active molecules. The website reports the chance of at least one interfering neighbor conditional on a candidate dye event at finite exposure.

## Sensitivity and the decision

For each dye we cross all three independent experimental means for ON dwell, OFF dwell and photon rate (27 sensitivity combinations). These quantities came from differing acquisition protocols. We also use two **shared** scales on OFF dwell, ×1 and ×0.5, to probe the dark-time upper bound reported for this protocol. They are stress scenarios. The measurements do not establish a lower dark-time bound or a universal shortening factor for all dyes.

Each row gives a nominal yield at the geometric center of the requested density range, with min/max yield over all scenarios and all densities in that interval. The curve can peak inside an interval, at ρ = 1/(A B(T)), so its upper bound checks that interior maximum. To certify a top-two pair, each of those dyes must exceed every outside candidate under all within-dye combinations at each shared dark-time scale, across the entire density interval. Log yield differences for a fixed pair of kinetics are affine in ρ; checking both endpoints proves this condition. Candidate overlap bands do not by themselves determine the recommendation. A dye is an unresolved contender if fewer than two others uniformly outrank it.

## Pilot summary

The importer in `analysis/pilot.py` records equal-area square tile counts for empty, single and multiple localizations, including empty frames. The square tile area belongs to the observation model and differs from the circular exclusion area. A matched fixed-cell pilot must identify the reference dye and use the same experimental domain and recorded photon cutoff as the requested comparison. S21 is a separate live-cell verification case.

For each density and reference-dye sensitivity scenario, the expected detections per tile in a frame are λ = tile area × ρ q exp(−ρ A B(T)). The three predicted category probabilities are exp(−λ), λ exp(−λ) and 1 − (1 + λ) exp(−λ). The joint multinomial likelihood uses all three categories once. Counts beyond two are coarsened. The raw total does not contribute a second likelihood. This Poisson tile observation model approximates spatially dependent detections following overlap rejection. Use a homogeneous field or crop, and retain zero-count frames.

The grid starts log-uniform inside the supplied density range with 181 points and trapezoidal weights. The likelihood is marginalized over equally weighted kinetic scenarios. It uses at most one field of tile observations, capped at 512, to limit pseudo-replication from repeat frames. A best-scenario multinomial deviance above 13.82 rejects a mismatched summary. This is a working quality-control cutoff. An informative update returns the equal-tail 90% density interval. Information gain below 0.05 nats retains the full input interval. These are conditional model estimates with no empirical coverage guarantee. Similar observed spot densities can arise on either side of the yield peak. The reported interval spans such alternative solutions.

## Reproduce

```sh
python -m analysis.dyes
python -m analysis.model
node analysis/predict.mjs
node analysis/predict.mjs settings.json pilot-summary.json
python -m unittest discover -s tests -p test_photophysics.py
node --test tests/model-pilot.mjs
```

For a high-density sensitivity example, use `{"densityMin":10000,"densityMax":10500,"minPhotons":100}` in `settings.json`. The default density interval is 100–1,000 labels/µm². Both intervals are illustrative assumptions. The Python extraction requires openpyxl, the stochastic tests use NumPy, and the JavaScript CLI requires Node.js 22.7 or later. The published web application runs in a modern browser.
