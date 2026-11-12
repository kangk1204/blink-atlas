# Run the blinded fixed-H2B experiment

Cal664sb has one registered conditional forecast and zero matched cellular measurements. Tyson et al. provide its free-parent pKcycl in sodium phosphate. The imaging in that study used Yale676sb-CA. A collaborating laboratory would need a Cal664sb HaloTag ligand and a fixed H2B acquisition, or existing, independently acquired records for a new parent. [Original titration](https://pmc-oa-opendata.s3.amazonaws.com/PMC8393207.1/oc1c00670_0003.jpg) · [Dated protocol and numerical predictions](https://github.com/kangk1204/blink-atlas/blob/4c68d5a39e2559d05e4b8b8144f4fba5bf2341aa/validation/h2b_blind/registration.md)

## Files to exchange

| Who and when | File | Contents |
|---|---|---|
| Sample custodian, before unblinding | [Enrollment CSV](../../docs/data/h2b-blind-roster.csv) | One blinded ID per distinct new parent and its measured free-parent pKa in the reference buffers. Keep the ID in the roster with a blank pKa if the titration is missing, censored or nonmonotonic. Keep structures and the identity key with the custodian. |
| Analyst, before unblinding | Prediction lock and SHA-256 | The whole roster, all three predictions, reference dyes, calibration source, scope and abstentions. Record the checksum independently with a date, for example in the electronic lab notebook. |
| Acquisition team, after the checksum is recorded | [Preparation measurements CSV](../../docs/data/h2b-blind-measurements.csv) | One row per dye and independent preparation. Supply mean ON in milliseconds, mean OFF in seconds and a source file, sheet or trace locator. |
| Acquisition team, with the measurements | [Assay record](../../docs/data/h2b-blind-assay.json) | Enter the lock's study ID, a reference to the actual acquisition, confirmations and the observed settings. The blank template has no default claim of protocol agreement. A private lab reference is sufficient. |

The numerical registration fixes the three methods for all cohorts. The individual cohort lock binds the actual IDs and inputs to an independently dated receipt. A checksum alone cannot establish when anyone first saw an outcome. Document labeling, excitation dose, detector and filtering thresholds with the acquisition reference before confirming `illuminationAndDetectionMatched`. Keep unpublished inputs and results outside the public site.

## Reproduce and score

Run the following commands from the repository root using the CPU environment in `requirements.txt`. Replace `/private/panel` with the existing directory containing your filled forms.

~~~sh
python -m analysis.h2b_blind lock --study H2B-panel-A --roster /private/panel/roster.csv --output /private/panel/locked.json
python -m analysis.h2b_blind score --lock /private/panel/locked.json --lock-sha256 RECORDED_SHA256 --measurements /private/panel/measurements.csv --assay /private/panel/assay.json --output /private/panel/scored.json
~~~

The first command prints the digest to record. The second validates it before opening the outcomes and scores the locked numbers. It averages the three within-preparation ratios for each dye, reports paired absolute log₁₀ error against both comparators, and retains missing and out-of-range parents. The candidate and the two comparators use the same eligible dyes. Incomplete primary cohorts remain marked incomplete.

To test the **original cross-buffer Cal664sb forecast**, use the [downloadable lock](../../docs/data/h2b-cal664sb-lock.json), the study ID `Cal664sb` and SHA-256 `c4e9ae62123abff64ca703a72b8a3c63e29e2ccec312a85f02f969d19eaaca8f`. To reproduce the locked file independently, run `python -m analysis.h2b_blind public --output cal664sb.json`. A subsequent matched-buffer titration belongs in a new primary lock prepared before release of cellular outcomes. Preserve and report both results under their original arms.

## Public source follow-up

On 12 November 2026 Europe PMC listed 14 citing records for the final Holland paper and 13 for its preprint, with overlap between lists. Searching for Cal664sb with H2B or histone returned one review. DataCite returned zero records for Cal664sb with H2B or HaloTag. Recent nanobody and expansion studies use JF635b, a training dye. The accessible Walterspiel archive uses live cells. No independent matched Cal664sb measurement was identified in this search. [Final paper citations](https://www.ebi.ac.uk/europepmc/webservices/rest/MED/41986687/citations?format=json) · [Preprint citations](https://www.ebi.ac.uk/europepmc/webservices/rest/MED/38766149/citations?format=json) · [Screened studies](../../docs/data/analogue-studies.json) · [Walterspiel archive](https://zenodo.org/records/16941576)

## Request for a collaborating laboratory

We have frozen three forecasts of the mean ON:OFF dwell-time ratio for a Cal664sb HaloTag ligand on H2B in fixed COS-7 cells. Do you have this ligand, or existing matched measurements for another hydroxymethyl Si-rhodamine absent from the Holland series? A sample custodian can first supply a coded parent ID and a free-parent titration from the reference buffer series. Once the analyst returns dated predictions, the custodian can release preparation-level ON and OFF times and the acquisition settings. An existing matched dataset would let us test the forecasts without new imaging.
