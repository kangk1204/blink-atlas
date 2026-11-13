# How the definition of the ON/OFF ratio affects dye prediction

Local interpolation has the lowest held-out average error under every tested aggregation of the Holland measurements. Nearest pKa wins three of the four individual comparisons. The independent fixed-cell H2B claim is still open. Results and sources were reviewed on 13 November 2026. [Interactive results](https://kangk1204.github.io/blink-atlas/#unseen-dyes)

## Recovering the published response

Columns B–D of Extended Data Fig. 4b–c and 4e–f contain three experiment means in milliseconds. Columns E–G describe medians. We used the mean columns and matched ligand identities across sheets. The reference is the full-precision, dimensionless ratio in 4h–i column B. The six calibration means in 4j column C repeat those values. Fig. 4n reports rounded percentages. The CSV preserves those displayed values and converts them to ratios. [Source workbook](https://pmc-oa-opendata.s3.amazonaws.com/PMC13167466.1/41592_2026_3062_MOESM10_ESM.xlsx) · [Original figure including 4n](https://pmc-oa-opendata.s3.amazonaws.com/PMC13167466.1/41592_2026_3062_Fig6_ESM.jpg)

For mean ON divided by mean OFF, divide the sum of the three ON means by the sum of the three OFF means. For the mean of ratios, fix the ON order B, C, D and permute the OFF columns. BDC pairs ON B with OFF B, ON C with OFF D and ON D with OFF C. BCD is the same-column assumption. Signed percentage difference is 100 × (reconstructed / reported − 1). The CSV includes the absolute difference in ratio units for every calculation. [Nine dyes, source cells, units and per-pair differences](../docs/data/ratio-audit.csv)

Eight ratio-of-means values reproduce the reference within 0.1%. JF626b gives 0.000300271 versus 0.000289684, a difference of +3.655%. Its six mean-of-ratios values lie 4.056–4.222% above the reference. The printed 4n value is 0.029%. The source of the discrepancy is unresolved. For JFX650b the same-column mean is 7.31% above the reference, although the ratio of means differs by only −0.038%.

The paper reports three independent cellular experiments. ON acquisition used 2 ms exposures and OFF acquisition used 25 ms exposures at approximately 100-fold lower labeling concentrations. The authors interpret the OFF times as upper bounds under these conditions. The released sheets and methods supply no identifier linking an ON column to an OFF preparation. The original prospective scorer expects identified preparation pairs and averages their individual ratios. A subsequent matched experiment can report both estimands and specify its primary rule before viewing outcomes. [Holland methods and caption](https://doi.org/10.1038/s41592-026-03062-5) · [Registered scoring rule](https://github.com/kangk1204/blink-atlas/blob/4c68d5a39e2559d05e4b8b8144f4fba5bf2341aa/validation/h2b_blind/registration.md)

### SD and SEM

The worksheet and caption agree on SD for 4h–i. For 4j the worksheet header D2 reads SD, and the caption specifies SEM. The six numeric SDs exactly match their 4h–i counterparts. A pixel check gives half-bar heights of about 0.000203 for JFX650b and 0.000140 for JF639b. Their SDs divided by √3 are 0.000205 and 0.000136. One ordinate pixel represents about 0.0000127. These two bars support the caption’s SEM interpretation. The smaller whiskers cannot reliably distinguish the two conventions. The 4n graphic also prints JF639b as 0.10 ± 0.2% while the underlying 4h–i SD is 0.0236 percentage points. All three fits use the means. [Original figure](https://pmc-oa-opendata.s3.amazonaws.com/PMC13167466.1/41592_2026_3062_Fig6_ESM.jpg) · [Pixel coordinates and source hashes](../docs/data/ratio-evidence.json)

## Whole-dye predictions

All methods use the same training dyes and log₁₀ response. Local interpolation uses the two bracketing pKas. Nearest pKa copies the closest training response, with a lower-pKa tie break. The original regression fits a line through all training pKas and log responses. An omitted dye is scored only within the remaining pKa range. JFX650b and JF626b are the endpoints, leaving four eligible folds out of six. Each target is fully removed before fitting. Mean absolute log₁₀ error weights each eligible dye equally.

| Response used for fitting and scoring | Local | Nearest pKa | Original regression |
|---|---:|---:|---:|
| Published ratios | 0.20035 | 0.21620 | 0.31771 |
| Mean ON divided by mean OFF | 0.20145 | 0.21622 | 0.32047 |
| Same-column mean of three ratios | 0.20182 | 0.22339 | 0.32197 |
| All 46,656 paired cases, min–max | 0.19934–0.20327 | 0.22000–0.22458 | 0.31936–0.32378 |

The enumeration allows six OFF permutations separately for each of six dyes. It includes the six cases where all dyes share a permutation. Each scenario uses the assigned permutations to reconstruct training labels and the held-out observation. The average error advantage over nearest pKa ranges from 0.01836 to 0.02361. JF635b supplies the only local win over nearest in every scenario. Nearest wins for JF639b, JF630b and JF629b. Regression has the smallest individual error for JF635b. These individual rankings also hold across the three named estimands. The ranges quantify sensitivity to possible column correspondences. The case counts carry no sampling probabilities. [Per-dye errors, shared-permutation subset and ranks](../docs/data/ratio-sensitivity.json) · [Every scenario and its forecasts](../docs/data/ratio-scenarios.csv.gz)

### What would change a decision for Cal664sb?

The input is the published pKcycl of 3.2 in sodium phosphate. The forecast remains conditional on transfer to the Holland assay. With absolute log error, adjacent predictions tie at their geometric mean. Below the lower boundary regression wins, between the boundaries nearest pKa wins, and above the upper boundary local interpolation wins.

| Quantity | Original registration | All 46,656 mean-of-ratios cases |
|---|---:|---:|
| Local interpolation | 0.000579845 | 0.000589278–0.000597359 |
| Nearest pKa | 0.000483915 | 0.000489600–0.000496886 |
| Regression | 0.000429728 | 0.000441501–0.000445282 |
| Regression / nearest boundary | 0.000456017 | 0.000464929–0.000470377 |
| Nearest / local boundary | 0.000529713 | 0.000537132–0.000544811 |

For a hypothetical observation of 0.000540, the registered forecasts favor local interpolation. In the sensitivity analysis either local or nearest can win depending on pairing. There is no matched Cal664sb outcome. The original registration JSON has SHA-256 c4e9ae62123abff64ca703a72b8a3c63e29e2ccec312a85f02f969d19eaaca8f. [Original forecasts](../docs/data/h2b-cal664sb-lock.json) · [Measured parent titration](https://doi.org/10.1021/acscentsci.1c00670)

## Additional public measurements

The two new panels cover 18 distinct parents, including HMSiR from the Holland nine. To score a frozen model we require an independently measured free-parent pKa within 1.22–4.48 and an independently observed compatible response. Per-model missing errors are JSON null. [Dye records and reasons](../docs/data/ratio-evidence.json)

- **Uno 2014.** Fourteen distinct parents occur in supplementary Fig. S2. Table 1 gives six parent titrations and open-state relaxation times or limits, alongside the flash photolysis in Fig. S4. AMAcRG is in range at 3.7. Local, nearest and regression forecast 0.00083253, 0.00101210 and 0.00058702 at that input. Its 3.09 µs flash relaxation lacks an independent mean OFF dwell for computing an error. The supplementary single-molecule experiments use adsorbed dye–IgG and the fixed-cell imaging targets microtubules. HMSiR is a known parent. Its absorption titration in Table 1 is 5.8 and its supplementary fluorescence readout is 5.9. [Uno et al., Table 1](https://www.nature.com/articles/nchem.2002/tables/1) · [Supporting information](https://static-content.springer.com/esm/art%3A10.1038%2Fnchem.2002/MediaObjects/41557_2014_BFnchem2002_MOESM177_ESM.pdf)
- **Chen 2024.** We inspected the publisher abstract and complete author-deposited supplement. The subscription article body was unavailable. Free-dye fluorescence plots S3A/C/E/G have descending half-heights near pH 6.6, 6.0, 6.2 and 5.8 for SMI-NIR1–4. We digitized the midpoint between displayed extrema and rounded to 0.1 pH. These optical estimates sit above the frozen range. Exact fitted parent pKas were unavailable in the inspected supplement. Fig. S7c supplies ten pre-bleach duty fractions for one biotin conjugate at five 721 nm powers in each of two buffers. At 3.93 kW/cm² the plotted fractions are about 0.0019 in PBS (pH 7.4) and 0.083 in MES (pH 4.7). S6 gives example traces for the other three parents. The JSON records dye, buffer and power for each digitized measurement. [Chen et al.](https://doi.org/10.1021/acs.analchem.4c02445) · [Author-deposited supplement](https://doi.org/10.1021/acs.analchem.4c02445.s001)

**A score on a different physical measurement is possible.** We refitted the three algorithms separately to Uno’s solution relaxation times in seconds after a 308 nm pulse at pH 7.4. Four Table 1 dyes provide exact inputs and lifetimes, two provide only upper limits. Holding out each exact-value parent gives two in-range tests. On HMAcRG the local, nearest and regression absolute log₁₀ errors are 0.931, 1.051 and 0.698. On AMSiR they are 1.866, 3.309 and 1.828. The averages are **1.398 / 2.180 / 1.263**. HMAcRG and AMSiR span different ring-closing chemistries, and the aminomethyl lifetime uses the slow decay component. The result concerns this solution assay and its own training set. [Table 1](https://www.nature.com/articles/nchem.2002/tables/1) · [Predictions, training IDs and per-parent errors](../docs/data/ratio-evidence.json)

Neither panel contributes a matched fixed COS-7 H2B dwell-ratio measurement or a scorable pair for the frozen transfer comparison. The earlier frozen Remmel IgG check scores four in-range dyes. Its mean errors are 0.293 local, 0.307 nearest and 0.208 regression. An independent H2B advantage and a change in the number of localizations passing a 20 nm threshold remain unverified. [Remmel measurements](https://pure.mpg.de/rest/items/item_3489427_3/component/file_3506946/content) · [Earlier benchmark](../docs/data/analogue-benchmark.json)

## Reproduce

Use Python 3.11 with the project requirements and Poppler on PATH. The model reads the bundled, checksum-verified workbook and original forecast lock. The external extractor downloads the original XML, figure and two supplementary PDFs and verifies their hashes. Further runs can use the populated cache offline. The Uno Table 1 constants and rounded 4n percentages are manual transcriptions identified in the source. Code extracts the Uno SI titrations and reads Chen’s plots from the pinned figures. All computations run on CPU.

~~~sh
python -m pip install -r requirements.txt
python -m analysis.ratio_sensitivity --output-dir docs/data
python -m tools.extract_ratio_evidence --cache .cache/ratio-originals
python -m tools.extract_analogue_sources --cache .cache/analogue-sources --verify-igg
python -m analysis.analogue_benchmark
python -m unittest discover -s tests
~~~

The gzipped CSV orders its first six columns by calibration dye. Each three-letter value specifies the OFF order relative to ON B, C, D. The remaining columns give the three fold-averaged errors, three Cal664sb forecasts and both equal-error thresholds. The JSON contains named-estimand folds, per-dye errors, extrema and rank counts. Its intervals retain more digits than the browser display.
