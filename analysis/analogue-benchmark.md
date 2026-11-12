# Unseen-dye benchmark and development plan

## Evaluation contract frozen before extracting the additional outcomes

The primary question is whether a new hydroxymethyl Si-rhodamine HaloTag ligand produces the predicted **mean ON:OFF time ratio in fixed COS-7 H2B–HaloTag** under the Holland acquisition and analysis protocol. The covariate is the independently measured apparent pKa of the free parent. The existing six exact Holland pairs and the previous leave-one-dye-out results are known development data.

**Candidate.** Use the two training dyes whose pKa values bracket the query. Linearly interpolate their measured **log10 time ratios** as a function of pKa. The local slope allows the relationship to vary across a chemical series. There are no fitted hyperparameters or outcome-selected interpolation weights. The primary comparator copies the measured time ratio of the training dye nearest in pKa. Also report the existing global log-linear regression. Ties at the nearest distance favor the lower pKa. If training dyes share a pKa, average their log ratios at that pKa for local predictions. An exact match uses that pKa's mean log ratio.

**Split and metric.** Hold out the entire chemical identity, including its aliases, replicates and repeated conditions. Never split spot, molecule or frame rows as independent dyes. Train every method on the same eligible dyes. Score only pKa values within the training interval. Require a positive empirical ratio or convert an explicitly defined single-molecule bright fraction f to f/(1−f). An absorbance titration used to obtain pKa cannot also supply its own prediction target. Preserve all exclusions. Report paired per-dye absolute log10 errors, mean errors, individual wins/ties/losses and the denominator before abstention. A positive improvement is baseline error minus candidate error. Keep null and negative results.

**Three evaluations.**

1. Report the Holland leave-one-dye-out result as development evidence. Those six dyes informed this design.
2. With the six-dye Holland model frozen, score every eligible additional dye from an independent matching H2B acquisition. Check chemical identity, sample and fixation, HaloTag, acquisition, analysis, pKa definition, and response definition against the primary assay. Report the count even when zero.
3. Where an independent publication reports paired experimental acidity and spontaneous blinking kinetics in a different assay, test transfer separately. Report predictions from the frozen Holland fit only for commensurate inputs and response definitions inside its range, explicitly labeling the assay shift. In addition, if the publication contains its own calibration series, run leave-one-dye-out **within one matched assay at a time**, using the same frozen algorithm. These fits use only the other dyes of that assay. Treat distinct acidity definitions and acquisition conditions as separate tasks. Never pool their errors with the primary endpoint.

**Decision.** Report whether the candidate's mean error is lower on paired eligible dyes. Report uncertainty and coverage before interpreting a small difference. Keep the existing reference visible. Promote a primary-assay replacement only on independent matched data with positive paired improvement and a 95% uncertainty interval above zero across independent dyes. A cross-assay win can motivate acquisition of matched data. A release claiming state-of-the-art performance additionally needs independent laboratories or batches and comparison with relevant contemporary methods on a common benchmark. Public retrospective scoring is not a prospective laboratory test.

## Realistic path to a stronger model

| Stage | Work and exit condition |
|---|---|
| Public-data pilot, this release | Trace new primary-source tables to their exact cells or printed values. Freeze the rule, run the three evaluations and publish every dye-level prediction, absence and abstention. The computation fits on a CPU. |
| Matched calibration, approximately 2–4 weeks after materials are available | Obtain existing fixed COS-7 H2B–HaloTag recordings or prepare additional monotonic hydroxymethyl SiR analogues spanning the calibration interval. Blind the structures' experimental outcomes. Fix reference dyes, exposure and sparse OFF acquisition before imaging. Record the per-preparation ratio, microscope settings and a control repeated on every day. |
| Modeling, approximately 1–2 weeks after a larger development panel is available | Compare the local model with nearest pKa and the original regression on grouped development folds. Test one additional independent covariate at a time, starting with protein-bound versus free-dye titration shift. Collect it before measuring cellular kinetics. Allow only training-fold selection. Estimate how many independent dyes a confirmation study needs from pilot paired-error variability and the minimum gain the lab considers actionable. |
| Independent confirmation, approximately 3–6 weeks conditional on instrument and dye access | Lock one model and protocol. Evaluate unseen dye identities in a separately acquired blinded panel, preferably at a second laboratory. Report per-dye error differences, group-level uncertainty, coverage and failure cases. Recheck performance against the strongest reproducible comparator and the same-budget nearest-pKa rule. |

## Public evaluation on 12 November 2026

The original evaluation contract is archived at [33951f3](https://github.com/kangk1204/blink-atlas/blob/33951f380301398a0f5e6bfff3d1a1dcca232b98/analysis/analogue-benchmark.md). The three additional publications supply 11 free-parent dyes and 25 dye–condition measurements. There were **zero independent unseen dyes with the matched fixed H2B endpoint** in the reviewed public sources. The default assay and existing reference models remain available. [Interactive results](https://kangk1204.github.io/blink-atlas/#unseen-dyes) · [Complete per-dye predictions, training IDs and errors](../docs/data/analogue-benchmark.json)

Mean absolute log₁₀ errors use exactly the same test dyes for all three methods within each row. Coverage counts the eligible and total parent dyes in that condition. Each comparison with study-specific training removes the test parent before fitting. Repeated powers and conjugates share parent identities.

| Test assay and training | Coverage | Local interpolation | Nearest pKa | Original regression | Local wins / ties / losses vs nearest |
|---|---:|---:|---:|---:|---:|
| Independent fixed H2B, frozen Holland fit | 0 available | — | — | — | — |
| Holland, leave one development dye out | 4/6 | 0.200 | 0.216 | 0.318 | 1 / 0 / 3 |
| Remmel IgG, frozen Holland fit | 4/6 | 0.293 | 0.307 | 0.208 | 2 / 0 / 2 |
| Remmel IgG, calibrate with the other five IgG dyes | 4/6 | 0.163 | 0.360 | 0.292 | 4 / 0 / 0 |
| Remmel HaloTag7, frozen Holland fit | 1/1 | 0.775 | 0.780 | 0.701 | 1 / 0 / 0 |
| Remmel nanobody, frozen Holland fit | 1/1 | 0.306 | 0.301 | 0.381 | 0 / 0 / 1 |
| Zheng purified HaloTag, fit at 500 W/cm² | 2/4 | 0.190 | 0.245 | 0.236 | 1 / 0 / 1 |
| Zheng purified HaloTag, fit at 1,000 W/cm² | 2/4 | 0.273 | 0.268 | 0.333 | 1 / 0 / 1 |
| Zheng purified HaloTag, fit at 1,500 W/cm² | 2/4 | 0.303 | 0.306 | 0.365 | 1 / 0 / 1 |
| Zheng purified HaloTag, fit at 2,000 W/cm² | 2/4 | 0.279 | 0.275 | 0.338 | 1 / 0 / 1 |

**What transfers.** The frozen IgG evaluation tests Remmel compounds 1–4. Compounds 5 and 6 have pKa₂ values of 4.60 and 5.00, above the Holland upper limit of 4.48. Mean paired improvement for the local rule was +0.0145 log₁₀ units. The descriptive paired-dye bootstrap interval was −0.0080 to +0.0425 (2.5th and 97.5th percentiles of all 4⁴ resamples with replacement). The original regression had the lowest mean error. Its own paired interval against nearest pKa was −0.0545 to +0.2542. The within-IgG result uses compounds 2–5 as tests with five other IgG-conjugated dyes as references in each fold. Its folds share training dyes. Neither panel estimates performance in fixed H2B. [Remmel et al., Table 1](https://doi.org/10.1002/smll.202206026)

**How the endpoint changes.** Remmel's upper pKa₂ represents the open to spirolactam equilibrium. A separate protonation at low pH produces a second dark state. We use pKa₂ as an explicitly approximate input for transfer. Divide mean ON time in ms by 1,000 times mean OFF time in seconds. The authors report both times rounded. Compound 4 has the same parent pKa₂ of 4.38 but observed ratios of 0.00294 on IgG, 0.000179 on purified HaloTag7 and 0.00216 on a nanobody. The first ratio is 16.5 times the second. The frozen local prediction is 0.00106 in all three cases. The linkers and protein environments both change. The IgG row repeated in Table S1 is counted once. [Article](https://doi.org/10.1002/smll.202206026) · [Authors’ supplement, Table S1 and Fig. S7](https://www.biorxiv.org/content/biorxiv/early/2022/08/29/2022.08.29.505670/DC1/embed/media-1.pdf?download=true)

Zheng's Table S2 records a mean bright fraction over all 6,000 frames at pH 7.4. We convert that fraction to occupancy odds `f/(1−f)`. A full observation window includes time after bleaching. All four free-parent titrations are outside the H2B training range. Only SRhB and STMR are interior when the local model is trained within one surface assay. Yale676sb has measured pKcycl 5.9, mean ON time 4.5 ms and OFF time 3.8 s on glass. It is outside the H2B range and has no same-assay training series in the curated pairs. [Zheng et al., Tables S1–S2](https://doi.org/10.1021/jacs.2c11395.s001) · [Tyson et al., Fig. 3 and single-molecule measurements](https://doi.org/10.1021/acscentsci.1c00670)

### Exploratory comparison of free and protein-bound titrations

After evaluating the frozen methods, we used the four independently measured bound titrations in Zheng Table S1 to change one input in the within-assay benchmark. The interpolator, nearest rule, held-out identities and training dyes stay fixed. Scoring requires support for both inputs. This question was added after the first outcomes were inspected. It has no part in the frozen H2B evaluation.

| Power (W/cm²) | Local, free pKa | Local, bound pKa | Nearest, free pKa | Nearest, bound pKa |
|---|---:|---:|---:|---:|
| 500 | 0.190 | 0.289 | 0.245 | 0.353 |
| 1,000 | 0.273 | 0.258 | 0.268 | 0.472 |
| 1,500 | 0.303 | 0.345 | 0.306 | 0.543 |
| 2,000 | 0.279 | 0.283 | 0.275 | 0.488 |

Bound-pKa interpolation improved on free-pKa interpolation at one of four powers. All four rows score the same two dyes. Binding-sensitive features merit a **paired measurement in the target assay**. Changing to the bound titration alone did not consistently improve this different endpoint. [Source titrations and illumination conditions](https://doi.org/10.1021/jacs.2c11395.s001)

## Discovery and source handling

On 12 November 2026 the Europe PMC H2B query returned 15 indexed records and the chemistry/kinetics query returned 31. The DataCite JF/fixed query returned zero records and the H2B/blinking/HaloTag query returned two Spot-On records. Publisher supplements, the authors’ preprints and Zenodo were checked for experimental labels. These are search yields for the stated queries. They place no bound on unpublished data or unindexed repositories. The machine-readable ledger records duplicates, other targets, censored titrations, computed quantities and incompatible acquisitions. [Screening ledger and primary sources](../docs/data/analogue-studies.json)

~~~text
Europe PMC
("H2B" OR "histone") AND ("spontaneously blinking" OR "blinkogenic") AND ("HaloTag" OR "Halo-Tag")
("hydroxymethyl" OR "Si-rhodamine" OR "spirolactam") AND ("duty cycle" OR "on:off ratio") AND ("blinking" OR "spirocyclization")
DataCite
("JF630b" OR "JF639b" OR "HM-DS655") AND ("H2B" OR "fixed")
("H2B" AND ("blinking" OR "HMSiR") AND "HaloTag")
~~~

The extractor verifies the SHA-256 of each original before reading it. It reads Remmel's Table 1 from the published PDF in the Max Planck repository and conjugate Table S1 from the authors’ supplement, links free and bound titrations to matching rows in Zheng Tables S1 and S2, and extracts Tyson's measured values from the full-text XML. The published Remmel Table 1 agrees with the preprint. The separate Holland workbook provides the frozen six-dye training data. The public release includes row-level measurements and references. [Remmel publisher PDF](https://pure.mpg.de/rest/items/item_3489427_3/component/file_3506946/content) · [Zheng supporting data record](https://api.figshare.com/v2/articles/22147942)

~~~sh
python -m pip install -r requirements.txt
# Poppler's pdftotext must be on PATH to re-extract the PDFs.
python -m tools.extract_analogue_sources
python -m analysis.analogue_benchmark
python -m unittest discover -s tests -p test_analogue_benchmark.py
node --test tests/probe-trial.mjs
~~~

The checked-in [observations CSV](../docs/data/analogue-observations.csv) suffices to rebuild predictions without downloading the PDFs. [Source extraction with checksums](../tools/extract_analogue_sources.py) · [Model, eligibility and paired scoring](analogue_benchmark.py)

## Next blinded measurements

The first prospectively registered parent is Cal664sb. The published pKcycl of 3.2 comes from a different buffer. Its cross-buffer forecasts are fixed before a matched H2B measurement and are available separately from the primary panel. [Dated registration](https://github.com/kangk1204/blink-atlas/blob/4c68d5a39e2559d05e4b8b8144f4fba5bf2341aa/validation/h2b_blind/registration.md) · [Blinded enrollment and scoring](../validation/h2b_blind/README.md)

Start with existing matched recordings for new parent identities. Include reference ligands on each preparation day. If new data are acquired, span the calibrated pKa interval with several chemically distinct series and measure matched free-parent titrations. Only measure a bound midpoint if the protein remains functional throughout the transition. Holland et al. encountered HaloTag instability at the required low pH. Retain linker identity and the method used to estimate dwell times. Use a reference ligand to match excitation intensity and photon detection across instruments. Keep the model developers blind to the independent panel until model and baselines are frozen. [Holland et al.](https://doi.org/10.1038/s41592-026-03062-5)

During development compare nearest pKa, global log-linear regression, local interpolation and a regression adding the bound-minus-free titration shift. A measured bound titration costs an extra assay. Report the extra information and measurement budget alongside each comparison. Split by chemical family when estimating transfer to a new series. In the final panel compute paired errors on the same eligible dyes and disclose abstentions, preparation variation and differences by laboratory. Choose sample size from a pre-agreed useful error reduction and between-dye paired variance from the calibration panel. Re-evaluate published structure/energy models if they can predict this exact experimental endpoint using the same available inputs. Existing spirocyclization and recruiting-rate methods address related quantities. [Chi et al.](https://doi.org/10.1002/anie.202010169) · [DFT recruiting-rate model](https://doi.org/10.1021/acs.jpca.4c06408)

For prospective individual trials, the workbench locks four forecasts before imaging and later scores three independent preparation ratios. Previously downloaded three-method cards retain their original forecasts. For a curated independent **public panel** run the command below. The panel JSON uses schema `blink-atlas/matched-h2b/1`, scope `hydroxymethyl-sir-monotonic-fixed-cos7-h2b`, the acquisition fields in `PROTOCOL`, a public HTTPS `source`, and the boolean declarations `independentAcquisition` and `identitiesReviewedAgainstHolland`. Each entry in `dyes` contains `name`, standardized free-parent `parentInchiKey`, optional `aliases`, `pKa`, `monotonicFreeTitration`, `sourceLocator`, and either `preparationRatios` or the directly `reportedMeanRatio`. An independent curator must check the molecular identities and experimental conditions against the original source. The code validates the inputs, retains their provenance and checksums, and predicts with the original training set. The unit tests use synthetic fixtures only to exercise this import.

~~~sh
python -m analysis.analogue_benchmark --matched-panel /path/to/curated-panel.json --output /path/to/scored-panel.json
~~~
