# A prospective test for a blinking HaloTag analogue

The workbench asks whether a measured free-dye apparent pKa can predict the corresponding ligand’s cellular mean ON:OFF time ratio. Its other readouts let a chemist inspect the structural change, buffer response and cellular context behind the choice of comparison dyes. The data are from the nine hydroxymethyl Si-rhodamines reported by Holland, Plutkis, English and colleagues. [Open the workbench](https://kangk1204.github.io/blink-atlas/#probe-design) · [Article and experimental methods](https://doi.org/10.1038/s41592-026-03062-5)

## Evidence chain

| Measurement | Material and acquisition | Source and transformation |
|---|---|---|
| Structure and spectral peaks | Free parent compounds 1–9. Absorption and emission maxima in citrate buffer at pH 3, 150 mM NaCl and 0.1% SDS. | Fig. 1a gives the wavelengths and displayed pKa values. Closed structures are cropped from Supplementary Note 2, physical PDF pages 41–44 (printed pages 40–43). |
| Quantum yield and pH titration | Two separate solutions for quantum yield. Two technical absorbance readings at each measured pH for each parent. The titration uses HCl/KCl below pH 2.2 and citrate/phosphate above pH 2.2. | Fig. 1a workbook B:C and Extended Data Fig. 2j–k workbook B:S. The display averages available duplicate readings. Both original readings, pH and cell addresses are retained. Missing pairs are omitted and measured zeroes are preserved. |
| Cellular blinking | Ligands 1HTL–9HTL bound to H2B–HaloTag in fixed COS-7 cells. Three independent measurements. ON and OFF statistics come from separate acquisitions. | Extended Data Fig. 4h–i gives the authors’ mean-based ratios and SDs. Extended Data Fig. 4j supplies six exact parent pKa / ligand ratio pairs. The endpoint is dimensionless. |
| ER SOFI | HaloTag–Sec61β in U2OS, 100 nM of each ligand for 30 min, 637 nm HILO, 10,000 frames at 100 Hz. Separate live and fixed groups. | Fig. 2b–d workbook. Nine dyes, three metrics, mean ± SEM across four cells. Percentage-point differences describe group means. |
| Lysosomal SOFI | U2OS, 20 nM JF639b or 100 nM JF614b for 90 min, 100 Hz. Methods name CD63 and captions name LAMP. | Extended Data Fig. 8e–f workbook. SNR and decorrelation time, mean ± SEM across five cells. The interface labels the compartment. The unequal doses and aggregate values limit interpretation. |

The ER example gives 41.16% and 29.02% useful pixels in live cells for JF639b and JF630b, respectively. In fixed cells the values are 15.31% and 52.19%. These SOFI readouts describe image fluctuations. The fixed H2B question about localization counts at ≤20 nm requires matched localization tables with acquisition area, duration and density. [Fig. 2 source workbook](https://pmc-oa-opendata.s3.amazonaws.com/PMC13167466.1/41592_2026_3062_MOESM8_ESM.xlsx) · [H2B precision evidence](https://kangk1204.github.io/blink-atlas/#precision)

## What the six dyes predict

Fit an unweighted least-squares line to `log10(mean ON:OFF ratio)` versus measured parent pKa. Each ligand contributes one point. The calibration set is JFX650b, JF639b, JF635b, JF630b, JF629b and JF626b. HM-SiR and JF646b have complex titrations and JF614b has an upper limit (<1). Their paired numeric values are absent from Extended Data Fig. 4j. [Published calibration data](https://pmc-oa-opendata.s3.amazonaws.com/PMC13167466.1/41592_2026_3062_MOESM10_ESM.xlsx)

For each retrospective test, remove the entire ligand, refit with the other five and score only pKa values within the resulting training interval. Compare the regression with copying the ratio of the nearest measured pKa and with the training median on the log scale. An equidistant nearest-reference tie selects the lower pKa. Apply the same eligibility and absolute log₁₀ error to all three methods.

| Method | Mean absolute log₁₀ error on the same four dyes |
|---|---:|
| Nearest measured pKa | 0.216 |
| Regression | 0.318 |
| Training median | 0.499 |

The regression had smaller error than the median on four of four dyes and smaller error than the nearest reference on one of four. Withholding JFX650b or JF626b placed the test pKa outside the training range. This is a retrospective comparison within one chemical series and one publication. A small mean error here does not establish transfer to another scaffold or instrument.

## Run one blinded analogue

1. Assign an ID to an unmeasured hydroxymethyl Si-rhodamine analogue. Measure the free parent's apparent pKa from a monotonic titration in the source buffer series. A pKa outside 1.22–4.48, a different tag/assay, or a complex titration gives an abstention.
2. Before reviewing cellular data, download the prediction card and retain a dated copy in the lab record. It fixes the regression, both references, their predictions, the model identifier, source-workbook hash and scoring rule. Include the ligands bracketing the new pKa as controls. The local timestamp is supplied by the browser.
3. In three independent preparations, measure the matched HaloTag ligand under the fixed COS-7 H2B protocol. Randomize the order of analogue and controls within each preparation. Keep illumination, labeling regime, detection threshold and analysis settings matched across them. The source ON protocol uses 20,000 frames at 2 ms exposure and approximately 3 ms cycle time with five-frame gap closing. The sparse OFF protocol uses 100,000 frames at 25 ms exposure and 27 ms cycle time. The authors treat OFF times as upper bounds. [Labeling, imaging and processing methods](https://doi.org/10.1038/s41592-026-03062-5)
4. For each preparation divide mean ON time by mean OFF time, after converting both to the same unit. Enter the three positive ratios. The card takes the arithmetic mean, computes absolute log₁₀ error for each frozen prediction and reports the closest method, retaining ties. Check reference-ligand agreement with the source and document assay deviations with the experiment. The independent statistical unit for accumulating prospective evidence is the new compound.
5. Save the scored JSON alongside the original plan and raw observations. Importing a card in a fresh browser recomputes every stored prediction and score. The files stay on the user's computer. No prospective measurements ship with the demo.

## Reproduce

From the repository root, install the CPU dependencies in `requirements.txt` and run:

```sh
python -m analysis.probe_design
python -m unittest discover -s tests -p test_probe_design.py
node --test tests/probe-trial.mjs
```

Rebuild the structure crops with `python -m pip install pypdfium2==5.9.0` and `python -m tools.probe_structures`. To reuse a local copy, pass `--pdf /path/to/supplement.pdf`. The renderer checks the PDF SHA-256 and crops the identified page rectangles. The five statistical workbooks are bundled unchanged and verified before extraction. Their URLs and checksums, plot coordinates, sheet cells and full held-out predictions are in the [machine-readable bundle](../docs/data/probe-design.json). Structure crops and derived experimental data are attributed to Holland et al. under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
