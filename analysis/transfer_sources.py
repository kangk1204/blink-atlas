"""Additional figure-verified parent assays for the v2 structure transfer experiment."""

from analysis import structure_pk as v1
from analysis.structure_catalog import ACETAMIDE, ACID, NET, NH2, NME, SI, arm, graph


UNO = "Uno 2018"
CHI = "Chi 2020"
EXTERNAL = (UNO, CHI)
DOIS = {"Sakabe": "10.1021/ja309688m", UNO: "10.1039/C7CC07783A", CHI: "10.1002/anie.202010169"}
CHI_FIG2 = "https://file1.lookchem.com/html/freearticle/2022081713/80859049/img_04.png"
CHI_FIG3 = "https://file1.lookchem.com/html/freearticle/2022081713/80859049/img_05.png"
UNO_SI = "https://www.rsc.org/suppdata/c7/cc/c7cc07783a/c7cc07783a1.pdf"


def additions():
    rows = []

    def row(study, dye, smiles, reported, status, locator, buffer, reason="", *, method="Absorbance"):
        mol = v1.smiles_mol(smiles) if smiles else None
        ident = v1.identity(mol) if mol else {"inchiKey": "", "smiles": "", "formula": ""}
        if status == "exact" and mol is None:
            raise ValueError("Missing eligible parent")
        rows.append(dict(study=study, dye=dye, **ident, reported=str(reported), status=status, locator=locator,
                         source="https://doi.org/" + DOIS[study], method=method, buffer=buffer, reason=reason,
                         used=status == "exact", evidence="original" if status != "reference" else "quotation"))

    for dye, donor, y in [("Ac-HMRG", NH2, 5.3), ("Ac-HMDiMeR", NME, 6.2), ("Ac-HMDiEtR", NET, 6.7)]:
        cosolvent = "0.1% DMSO" if dye == "Ac-HMRG" else "<0.5% acetonitrile"
        row("Sakabe", dye, graph(ACETAMIDE, donor), y, "exact", "Fig. 1a,b and Table 1; SI Fig. S2",
            "0.2 M sodium phosphate; " + cosolvent)

    phosphate = "0.2 M sodium phosphate; 0.5 µM in Fig. 2c"
    tfe = "CC(F)(F)F"
    row(UNO, "HMbisTFER", graph(arm("N" + tfe)), 7.2, "readout unresolved", "Table 1; Figs. 1, 2c; SI compound 4, Fig. S3c,d",
        phosphate, "A single midpoint is reported. Absorbance and fluorescence spectra are both supplied; the fit channel could not be verified. Included in optical sensitivity only.",
        method="Optical pH fit, channel unverified")
    four = arm("N(" + tfe + ")" + tfe)
    row(UNO, "HMtetTFER", graph(four, pendant=ACID), "3.3 / 5.0", "complex",
        "Table 1; SI compound 8, spectroscopic methods and Fig. S3a,b", phosphate, "Two resolved transitions.")
    he = graph(four, pendant=ACID).replace(ACID + "CO2", ACID + "CCO2")
    row(UNO, "HEtetTFER", he, "3.1 / 5.5", "complex", "Table 1; SI compound 13 and spectroscopic methods",
        phosphate, "Hydroxyethyl spiro ring and two resolved transitions.")
    for name, smi, value in [("HMRG", graph(NH2), 8.1), ("HMAcRG", graph(ACETAMIDE, NH2), 5.3),
                             ("HMSiR", graph(NME, bridge=SI, pendant=ACID), 5.8)]:
        row(UNO, name, smi, value, "reference", "Table 1 prior-work reference", "See original assay",
            "Blocks the reference identity in the external split; no added label.", method="Literature")
    for name, value in [("HMRef", 10.2), ("HMRet", 10.9)]:
        row(UNO, name, "", value, "outside scope", "SI Fig. S2, rhodol side-series", "See SI",
            "Rhodols lack the two aryl amino donors of the registered model.", method="SI photophysical panel")

    water = "Aqueous pH titration in Fig. 3; buffer composition unspecified in figure"
    row(CHI, "HM-DS655", graph(arm("N3CCC3"), bridge=SI, pendant=ACID), 5.3, "exact",
        "RDM1, Figs. 2c,e and 3a,c", water, "Original monotonic peak absorbance; para carboxyl.")
    row(CHI, "HM-DS531", graph(arm("N3CC[N+](C)(C)CC3")), 5.3, "outside scope",
        "RDM5, Figs. 2c,e and 3d,f", water,
        "Two permanent quaternary ammonium charges on the closed free parent; outside neutral-parent model.")
    for name in ("RDM2", "RDM3", "RDM4"):
        row(CHI, name, "", "", "calculated only", "Fig. 2c,e", "", "No original measured midpoint.", method="Quantum calculation")
    return rows
