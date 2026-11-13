# Explaining the Holland error and testing an expanded structure model

The deployed nearest-structure predictor now uses 29 unique, measured free parents. The new additive ridge model uses an explicit description of the bridge, ring closure and amino substituents. It improves on nearest structure on three newly curated, completely held-out parents. Whole-study tests in Holland and Tyson favor the expanded nearest-structure predictor. Both forecasts are available for new structures. [Web workbench](https://kangk1204.github.io/blink-atlas/#probe-design) · [Per-dye forecasts and errors](../docs/data/structure-transfer/per-dye.csv)

## Reproducing the original error

Holding out Holland left 20 eligible parents, of which five contained a silicon bridge and zero had a four-membered azetidine donor. The six Holland measurements span 1.22–4.48. The training labels span 3.20–10.80 and average 7.42. Four test values fall below the lowest training label. All six nearest neighbors have oxygen bridges. Their Tanimoto similarity to the test parent is only 0.383–0.424. The original fingerprint ridge compresses the predicted range to 7.58–7.75. [Original folds](../docs/data/structure-pk/benchmark.json)

| Holland free dye | Measured | Original nearest | Original error, predicted minus measured | Expanded nearest |
|---|---:|---|---:|---:|
| JFX650b | 4.48 | HMTMR 9.50 | +5.02 | 5.30 |
| JF639b | 3.97 | HMTMR 9.50 | +5.53 | 5.30 |
| JF635b | 2.95 | HMRB 9.20 | +6.25 | 5.30 |
| JF630b | 2.58 | HMTMR 9.50 | +6.92 | 5.30 |
| JF629b | 2.25 | HMTMR 9.50 | +7.25 | 5.30 |
| JF626b | 1.22 | HMRB 9.20 | +7.98 | 5.30 |

Chi's HM-DS655 has both a silicon bridge and azetidine donors. Its carboxylated pendant aryl distinguishes it from the Holland parents. It becomes the nearest reference for all six. Tanimoto increases to 0.492–0.525 and nearest MAE falls from 6.492 to 2.392. **Removing this single reference while retaining the other additions restores 6.492.** The expanded baseline still overpredicts all six and assigns them one common rank. The change identifies a cause in the model and dataset. Structural differences and buffer differences are confounded across studies, so separate experimental contributions cannot be estimated from this result. [Chi structures and aqueous absorbance curves](https://doi.org/10.1002/anie.202010169) · [Holland measurements](https://doi.org/10.1038/s41592-026-03062-5)

## Plan, sources and model

1. Reconstruct signed errors against the archived 26-parent release. Freeze the challenger, comparisons and rule for selecting a default before fitting. Holland and Tyson results were already known from the earlier release. [Primary protocol](../validation/structure_pk/transfer-registration.json)
2. Trace new free-parent structures and their optical measurements. Add Sakabe's Ac-HMDiMeR (6.2) and Ac-HMDiEtR (6.7), and Chi's neutral HM-DS655 (5.3). Sakabe's Ac-HMRG (5.3) maps to the known HMAcRG graph and receives one chemical-identity weight. Chi Fig. 2e distinguishes the experimental 5.3 from the calculated 5.65. Fig. 3c shows the directly observed decreasing absorbance. The other measured Chi dye HM-DS531 is a permanent dication outside the neutral free-parent domain. RDM2–4 have calculated values only. [Sakabe original study](https://doi.org/10.1021/ja309688m) · [Sakabe original spectra](https://doi.org/10.1021/ja309688m.s001) · [Chi Fig. 2](https://file1.lookchem.com/html/freearticle/2022081713/80859049/img_04.png) · [Chi Fig. 3](https://file1.lookchem.com/html/freearticle/2022081713/80859049/img_05.png)
3. Preserve unresolved labels. Uno reports 7.2 for HMbisTFER. Supplementary Fig. S3 includes both absorption and fluorescence, and the methods permit either fit. We verified the value, parent graph and reported single transition. We could not independently establish the optical channel for that exact number. Count it only in the more inclusive optical sensitivity analysis. Uno's HMtetTFER and HEtetTFER both have two resolved transitions. Zheng's HM-JF526 synthesis and protein imaging and the lactone equilibrium constants in the Grimm/Zheng series supply relevant structural precedent. A measured single-transition free-parent pKcycl could not be established there. The independent Chi quantum equation also uses prior measurements for calibration. Its output never supplies a training label. [Uno Table 1](https://doi.org/10.1039/C7CC07783A) · [Uno supporting methods and spectra](https://www.rsc.org/suppdata/c7/cc/c7cc07783a/c7cc07783a1.pdf) · [Zheng](https://doi.org/10.1021/acscentsci.9b00676) · [Grimm](https://doi.org/10.1038/nmeth.4403) · [Assay and exclusion ledger](../docs/data/structure-transfer/assays.csv)
4. Train a ridge regression on 11 features with an unpenalized intercept and alpha 1. The features encode the bridge, amino spiro closure, pendant acid, 3′ and 6′ crowding, N substitution, acylation, ring constraint, ring fusion, four-membered donors and a donor-branch electron-withdrawing score. Branches stop at the aromatic core. The last score sums the Pauling electronegativity excess over carbon of O, halogens and nitrile N, halving the contribution for each additional bond of graph distance. Average the two donor environments. This score is a chemical heuristic. The previous kernel and both nearest methods use the original radius-2, 2,048-bit Morgan fingerprints. All comparators share exactly the same training data per fold. Feature centering, scaling, target means and model weights are computed on training parents only. [Implementation](structure_transfer.py) · [Exact assay transcriptions](transfer_sources.py)

## Unseen dyes

The unit of analysis is one chemical identity. Splits block every structurally identified parent on the held-out titration roster, including reference dyes. Each row below is a separate retrospective experiment. The combined new-parent test excludes Sakabe, Uno 2018 and Chi 2020 simultaneously and scores only the three parents absent from the original 26. It shares dyes with the study analyses. The combined check and channel sensitivity were specified after the first model fits and before computing the combined test. They do not select the deployed model. [Evaluation extension](../validation/structure_pk/transfer-external-extension.md)

| Data and split | Train | Test | Nearest MAE | Fingerprint ridge | Core and donor ridge | Training mean |
|---|---:|---:|---:|---:|---:|---:|
| Original corpus, Holland held out | 20 | 6 | 6.492 | 4.751 | 4.367 | 4.512 |
| Expanded corpus, Holland held out | 23 | 6 | 2.392 | 4.400 | 2.892 | 4.335 |
| Expanded corpus, Tyson held out | 24 | 5 | 1.820 | 1.349 | 2.632 | 1.368 |
| Chi and Uno 2018 both held out, Chi scored | 26 | 1 | 1.600 | 0.046 | 1.330 | 1.060 |
| Same split, Uno 2018 scored | 26 | 0 | — | — | — | — |
| Three newly curated parents jointly held out | 22 | 3 | 2.467 | 0.477 | 0.617 | 0.499 |
| One Holland dye out, siblings in training | 28 | 6 folds | 1.180 | 1.517 | 1.051 | 3.561 |

All errors are in pH units. On the primary Holland fold the new regression beats expanded nearest on two of six dyes. Removing the electronic feature worsens its MAE from 2.892 to 4.515. This helps account for the measured shift but does not resolve the cross-study gap. In the sibling-assisted analysis it wins on three of six and its MAE advantage changes sign when either JF635b or JF626b is omitted from the average. Full-study and within-family generalization answer different experimental questions.

| Newly curated test parent | Measured | Nearest | Fingerprint ridge | Core and donor ridge |
|---|---:|---:|---:|---:|
| Ac-HMDiMeR | 6.20 | 9.50 | 7.20 | 6.18 |
| Ac-HMDiEtR | 6.70 | 9.20 | 7.12 | 6.18 |
| HM-DS655 | 5.30 | 6.90 | 5.30 | 3.98 |

The challenger beats nearest on 3/3 jointly unseen identities from two studies. The unseen labels cover a narrow interval. A constant estimated from the training parents has lower error than the challenger on this set, and so does the original fingerprint ridge. In the broader optical sensitivity, adding HMbisTFER as a fourth test parent changes the paired MAEs to 2.350 (nearest), 1.307 (new regression), 0.402 (fingerprint ridge) and 0.649 (training mean). The new regression wins 3/4. The source checks do not establish superiority across unseen studies or prospective accuracy. The default is the nearest measured dye in the **expanded** corpus. [Per-dye results](../docs/data/structure-transfer/per-dye.csv) · [Split membership, weights, individual deletions and ablations](../docs/data/structure-transfer/benchmark.json)

## Use and reproduce

~~~sh
python3.11 -m venv .venv
. .venv/bin/activate
python -m pip install -r analysis/requirements-structure.txt
python -m tools.structure_sources --cache /tmp/blink-originals
python -m tools.transfer_sources --cache /tmp/blink-originals
python -m analysis.structure_transfer build --out /tmp/blink-transfer
diff -rq docs/data/structure-transfer /tmp/blink-transfer
python -m unittest discover -s tests -p test_structure_transfer.py
python -m analysis.structure_transfer predict --sdf docs/data/structure-pk/HM-SiR-CH2F.sdf --output result.json
python -m analysis.structure_transfer predict --smiles 'C[Si]1(C)c2cc(N3CC(CF)C3)ccc2C2(c3ccccc3CO2)c2ccc(N3CC(CF)C3)cc21' --name HM-SiR-CH2F --output same.json
~~~

Open a result JSON in the web workbench. RDKit computes the fingerprint and descriptor vector locally. The browser verifies the released numerical predictions from these inputs. Previously registered forecasts and candidate structures remain available for comparison. For –CH₂F and –CHF₂, expanded nearest predicts 2.95 for both. The new regression predicts 3.22 and 2.30, respectively. These are exploratory forecasts for the same unsynthesized parents. The original matched titration proposal and the fixed-cell H2B registry use their own frozen protocols. [Matched titration](../validation/structure_pk/titration.md) · [H2B registration](../validation/h2b_blind/registration.md)
