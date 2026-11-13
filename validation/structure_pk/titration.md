# A matched test of two fluorinated Si-rhodamines

**Decision for a collaborating synthesis lab.** Measure the free-parent spirocyclization midpoint of the two proposed dyes alongside JF626b. The two candidates have the same nearest measured analogue, JF635b. The structure models have little power to resolve their order. A paired titration can test whether adding a fluorine on each arm lowers the midpoint.

| Free parent | Group on each azetidine C3 | Role |
|---|---|---|
| HM-SiR-CH2F | –CH₂F | Prospective free parent |
| HM-SiR-CHF2 | –CHF₂ | Prospective free parent |
| JF626b | –CF₃ | Matched scaffold control, published apparent pKa 1.22 |

Freeze the predicted experimental order CH₂F > CHF₂ > CF₃ before acquiring spectra. The public registration contains atom-level structures. Holland reports preparation of control 8 from ditriflate S4 and 3-(trifluoromethyl)azetidine hydrochloride by Pd-catalyzed C–N coupling (method S1, reported 92%). The two alternative azetidines are recorded building blocks. Test their suitability for the same coupling and characterize each isolated free parent by LC–MS and NMR before the assay. Candidate synthesis, yield, purity, titrations and biological performance are pending. [Frozen structures](registration.json) · [Holland SI, compound 8](https://pmc-oa-opendata.s3.amazonaws.com/PMC13167466.1/41592_2026_3062_MOESM1_ESM.pdf#page=44) · [Grimm et al.](https://doi.org/10.1038/nmeth.4403)

### Measurements

1. On each of three independent preparation days, make fresh dilutions of both free dyes and the control. Use matched dye concentrations, solvent fractions and plate geometry. Randomize coded sample positions on the plate. Take buffer blanks at each pH and duplicate instrumental readings for each dilution. The independent unit is the preparation day.
2. Use the Holland HCl/KCl series below pH 2.2 at ionic strength 0.2 M, and mixtures of 0.1 M citric acid and 0.2 M disodium phosphate above pH 2.2. Use the published Extended Data 2j–k pH grid covering both plateaus where the buffer permits. Record measured pH, temperature, cosolvent and absorbance spectrum for every sample. Identify the visible absorption band of the open dye and record the wavelength used. Holland used two measurements for its titration. This proposed assessment uses three independent days. [Holland methods](https://doi.org/10.1038/s41592-026-03062-5) · [Titration workbook](https://pmc-oa-opendata.s3.amazonaws.com/PMC13167466.1/41592_2026_3062_MOESM9_ESM.xlsx)
3. For each dye and day, plot blank-corrected absorbance against pH. If the trace is monotonic with two observable plateaus, fit a decreasing sigmoid with free plateaus and positive slope magnitude. Define pKcycl at half the absorbance interval. Save the raw spectra, fitted curve and residuals. For a bound, report the observable interval. For a nonmonotonic trace, report the curve and turning points and withhold the single-midpoint score.
4. Unblind after the fits and eligibility calls are saved. Report each day's CH₂F minus CHF₂ and CHF₂ minus CF₃ values and the mean paired differences across eligible days. Record how many days agree with the frozen order, any reversal, and each control's difference from 1.22. Compare the candidate measurements with both frozen numerical forecasts in pH units. Three days support a repeatability check. Keep this assay distinct from HaloTag conjugation and the fixed-cell H2B blind registry.
