# From a free dye structure to a titration experiment

The Probe design page estimates the apparent midpoint of spirocyclization of a free (Si-)rhodamine in solution. It uses a binary molecular fingerprint. Distances to measured dyes and a kernel regression give two testable estimates in pH units. The separate H2B predictor and its blind scoring protocol still require a **measured** free-dye pKa. [Web workbench](https://kangk1204.github.io/blink-atlas/#probe-design) · [Design frozen before scoring](https://github.com/kangk1204/blink-atlas/commit/5e43841) · [Titration proposal](../validation/structure_pk/titration.md)

### Reproduce

Use Python 3.11 on a CPU with an RDKit wheel (the reference environment is Linux x86-64). From the repository root:

~~~sh
python3.11 -m venv .venv
. .venv/bin/activate
python -m pip install -r analysis/requirements-structure.txt
python -m tools.structure_sources --cache /tmp/blink-atlas-primary
python -m analysis.structure_pk build --out /tmp/blink-atlas-release
diff -rq docs/data/structure-pk /tmp/blink-atlas-release
python -m unittest discover -s tests -p test_structure_pk.py
~~~

The source verifier fetches twelve SHA-256-pinned originals and checks Holland workbook cells and the two Tachibana tables. Drawings and raster-figure labels are human transcriptions in the catalog. Building and using the model from those committed transcriptions works offline. SDFs contain 2D coordinates and identifiers. The reference collection includes the open cation of the 2MeSiR noncyclizing control, tagged with its form. All other entries use closed graphs. [Structure transcriptions](structure_catalog.py) · [Primary files and checksums](../tools/structure_sources.py) · [Source rows and reasons](../docs/data/structure-pk/measurements.csv) · [Reference SDF](../docs/data/structure-pk/parents.sdf)

~~~sh
python -m analysis.structure_pk predict --sdf docs/data/structure-pk/HM-SiR-CH2F.sdf --output my-dye.json
python -m analysis.structure_pk predict --name HM-SiR-CH2F --smiles 'C[Si]1(C)c2cc(N3CC(CF)C3)ccc2C2(c3ccccc3CO2)c2ccc(N3CC(CF)C3)cc21' --output from-smiles.json
~~~

Supply one neutral, closed free parent and preserve explicit isotopes. The CLI accepts one record and checks the core and charge. Disconnected water and simple alkali-metal or halide counterions are stripped. For organic mixtures, salts requiring protonation changes, conjugates and open dyes, supply the appropriate neutral closed parent. The supported chemical series is narrow and predictions outside the training data may have large errors. In Probe design, expand **Open a result from your own structure** and select its JSON. The browser checks both numerical outputs against the released weights and fingerprint. RDKit assigns identity locally. For known complex, censored or undetermined measurements the CLI retains the observation and withholds a single numerical estimate. Known exact parents carry their original measurements and an explicit training-reference flag.

### Which measurements enter the benchmark?

| Study | Structures represented | Eligible original parents | Measurement decision |
|---|---:|---:|---|
| Holland | 9 | 6 | Absorbance, Fig. 1a and ED 2j–k. The six exact values are also in ED 4j B3:B8. Preserve the complex HM-SiR trace, undetermined JF646b and JF614b <1 separately. |
| Tyson | 7 | 5 | Figures 2 and 3, SI Figs. S2a–g and methods §D. HMSiR and HMSiR2-FlEt have two titratable groups. HMSiRTHQ shows one decreasing transition in both figures. |
| Uno | 14 | 9 | Fig. 2 absorbance is primary. SI Fig. S2 fluorescence is an alternate readout of the same assay. Two bell-shaped SiR curves, two bounds and one noncyclizing control are retained. |
| Tachibana | 21 | 4 | Table 2 / SI Figs. 16–19 give original measurements for four substituted HMRGs. Table 1 reuses 17 literature measurements. |
| Sakabe | 3 | 3 | HMRG, HMDiMeR and HMDiEtR have original absorption curves in SI Fig. S1. The HMRG identity also occurs in Uno. |
| Iwatate 2018 | 2 | 0 | HMSiR620h and HMJSiR both have two measured transitions in SI Table S1. Tachibana Table 1 uses their upper branches. |

Together these records describe 44 unique chemical identities including the open control, with 26 eligible training parents. InChIKey groups aliases and retained isotopes. Each original study contributes at most one absorbance value per parent. Independent study values are averaged within a parent before fitting. Quoted values carry their original or unresolved attribution and contribute no training label. The Iwatate 2016 paper identifies the asymmetric series. Its numerical supplement could not be independently inspected in this audit. Five entries sourced to that study remain secondary-only. AMRB and two discrepant unsubstituted SiR entries also retain unresolved row-level provenance. [Holland](https://doi.org/10.1038/s41592-026-03062-5) · [Tyson](https://doi.org/10.1021/acscentsci.1c00670) · [Uno](https://doi.org/10.1038/nchem.2002) · [Tachibana](https://doi.org/10.1038/s42004-020-0326-x) · [Sakabe](https://doi.org/10.1021/ja309688m) · [Iwatate 2016](https://doi.org/10.1002/chem.201503426) · [Iwatate 2018](https://doi.org/10.1021/acs.bioconjchem.7b00776)

Tyson Fig. 2 omits the para carboxyl drawn in Fig. 3. The SI synthesis schemes and HRMS formulas establish the additional CO₂ unit for the synthesized parents. Uno Fig. 2 draws the carboxylated HMSiR too. Holland compound 1 is the corresponding carboxyl-free parent. JFX650b retains sixteen deuteriums. Halogens and amines on asymmetric xanthene arms and 3′/5′ pendant-aryl isomers are distinguished explicitly in the graphs.

### Retrospective validation

Test on all eligible dyes from each held-out study. All identities in that study's titration roster, including ineligible controls and aliases, are removed from the other studies before fitting. The test midpoint and its magnitude play no part in fold assignment. Morgan radius 2, 2,048 bits and chirality are fixed. Nearest structure uses binary Tanimoto and averages all maximum-similarity ties. Kernel ridge uses the same Tanimoto kernel, centered training targets and kernel, an unpenalized intercept and alpha 1. The model receives only the parent graph. Buffer differences remain part of the transfer error.

| Held out | Scored | Nearest MAE | Ridge MAE | Shared train/test structures |
|---|---:|---:|---:|---:|
| Holland | 6 / 9 | 6.492 | 4.751 | 0 |
| Tyson | 5 / 7 | 1.868 | 1.360 | 0 |

Removing Yale676sb from the Tyson paired average reverses the ranking (nearest 1.598, ridge 1.687). The preregistered deployment rule therefore retains nearest structure. Holland’s held-out azetidine series is poorly predicted by both models. This retrospective benchmark supplies no prospective accuracy estimate for the two proposed compounds. [All per-dye errors](../docs/data/structure-pk/holdouts.csv) · [Training identities, removal checks and model](../docs/data/structure-pk/benchmark.json)

For Tachibana Table 1, each structure model excludes the query parent globally. Seven rows match an original single-transition measurement with the same graph and value. Their MAEs are 0.900 for nearest, 1.015 for ridge and 0.200 for the published quantum calculation. All 17 reported table values are also preserved. The first eight parents informed hydration/protocol development. Two had their hydration treatment changed after comparison with experiment. The five-water and initial three-water numbers are both provided. Errors use the displayed measured and computed numbers, so rounding can yield a different last decimal from the paper's error column. There is no published quantum calculation for these two virtual candidates in this comparison. [Full crosswalk](../docs/data/structure-pk/quantum.csv) · [Original table](https://www.nature.com/articles/s42004-020-0326-x/tables/1)

### Candidate identity and search coverage

The query pair was fixed before scoring. On 2026-11-13 PubChem PUG REST returned no CID for either full closed-parent InChIKey or for the corresponding open and closed same-connectivity searches. Each building block was positively identified and a connectivity-service control returned results. Exact keys and names combined with rhodamine/SiR/fluorophore in Europe PMC returned no indexed matching record. Targeted literature searches for the two azetidine names, Si-rhodamine and hydroxymethyl led to the published JF synthesis series and other scaffolds. Neither candidate was identified as a reported full dye in the materials inspected. Structures drawn only in schemes, patents using generic substituents and unindexed chemistry can escape these searches. [Timestamped queries and responses](../validation/structure_pk/identity-search.json) · [Search script](../tools/structure_identity_search.py) · [Candidate/control titration protocol](../validation/structure_pk/titration.md)
