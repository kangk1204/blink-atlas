"""Source transcriptions for free dyes. SMILES are closed neutral graphs redrawn from the cited figures."""

from analysis.dyes import DYE_NAMES
from analysis.probe_design import CHEMISTRY


DOI = {
    "Holland": "10.1038/s41592-026-03062-5", "Tyson": "10.1021/acscentsci.1c00670",
    "Uno": "10.1038/nchem.2002", "Tachibana": "10.1038/s42004-020-0326-x",
    "Sakabe": "10.1021/ja309688m", "Iwatate 2016": "10.1002/chem.201503426",
    "Iwatate 2018": "10.1021/acs.bioconjchem.7b00776",
}
BUFFER = {
    "Holland": "HCl/KCl below pH 2.2 (0.2 M ionic strength); 0.1 M citric acid / 0.2 M disodium phosphate series above",
    "Tyson": "§D: 2 µM dye, 0.2 M sodium phosphate, pH 2–12 in 0.5 steps, room temperature, triplicate; SI Fig. S2: 2–5 µM, 0.1% DMSO",
    "Uno": "0.2 M sodium phosphate, 0.1% DMSO (Fig. 2 and SI Fig. S2)",
    "Tachibana": "0.1 M sodium phosphate, 0.1% DMSO (SI Figs. 16–19)",
    "Sakabe": "0.2 M sodium phosphate, 0.1% DMSO (SI Fig. S1)",
}
PARENT = "c3ccccc3"
ACID = "c3ccc(C(=O)O)cc3"
SI = "C[Si]1(C)"


def arm(donor, *, before="", after=""):
    def carbon(group):
        return "c(" + group + ")" if group else "c"
    return carbon(before), carbon(donor), carbon(after)


def julolidine(first, second):
    return arm(f"N{first}{second}", before=f"CCC{first}", after=f"CCC{second}")


def fused(n_sub, tether, closure):
    return arm(f"N{closure}{n_sub}", after=tether + closure)


def graph(left, right=None, bridge="O1", closure="O", pendant=PARENT):
    """Walk one aromatic arm from bridge to spiro center and the other in reverse."""
    right = left if right is None else right
    return (bridge + "c2" + "".join(left) + "cc2C2(" + pendant + "C" + closure + "2)c2c"
            + "".join(reversed(right)) + "c21")


NH2, NME, NET = arm("N"), arm("N(C)C"), arm("N(CC)CC")
ACETAMIDE = arm("NC(=O)C")


def catalog():
    rows = []

    def add(study, dye, smiles, value, locator, *, status="exact", reason="", method="Absorbance",
            buffer=None, form="closed", source=None, evidence="original", origin=None):
        rows.append(dict(study=study, dye=dye, smiles=smiles, reported=str(value), status=status,
                         reason=reason, method=method, buffer=buffer or BUFFER.get(study, ""),
                         locator=locator, source=source or "https://doi.org/" + DOI[study],
                         origin=origin or study, evidence=evidence, form=form))

    pyr_d8 = "N3" + "C([2H])([2H])" * 4 + "3"
    groups = ["N(C)C", pyr_d8, "N3CCC3", "N3CC(OC)C3", "N3CC(F)C3",
              "N3CC(S(=O)(=O)C)C3", "N3CC(C#N)C3", "N3CC(C(F)(F)F)C3", "N3CC(F)(F)C3"]
    for i, (name, spec, group) in enumerate(zip(DYE_NAMES, CHEMISTRY, groups, strict=True), 1):
        status = {"HM-SiR": "complex", "JF646b": "undetermined", "JF614b": "upper-bound"}.get(name, "exact")
        reason = {"complex": "Nonmonotonic titration; reported apparent midpoint retained.",
                  "undetermined": "Complex titration; no defined midpoint.",
                  "upper-bound": "Reported only as <1."}.get(status, "")
        add("Holland", name, graph(arm(group), bridge=SI), spec[3],
            f"Fig. 1a, compound {i}; ED Fig. 2j–k; Supplementary Note 2, PDF p. {spec[4]}",
            status=status, reason=reason)
        rows[-1]["substitution"] = f"Two {spec[0].lower()} donors; ortho hydroxymethyl, no pendant CO2H"

    acid_si = graph(NME, bridge=SI, pendant=ACID)
    uno = [
        ("HMAcRG", graph(ACETAMIDE, NH2), 5.3, 5.3),
        ("AMAcRG", graph(ACETAMIDE, NH2, closure="N"), 3.7, 3.7),
        ("HMSiR", acid_si, 5.8, 5.9),
        ("AMSiR", graph(NME, bridge=SI, pendant=ACID, closure="N"), 4.7, 4.7),
        ("MMSiR", graph(NME, bridge=SI, pendant=ACID, closure="S"), "<2", "<2"),
        ("2MeSiR", "C[Si]1(C)c2cc(N(C)C)ccc2[C+](c2ccc(C(=O)O)cc2C)c2ccc(N(C)C)cc21", ">11", ">11"),
        ("HMRG", graph(NH2), 8.1, 8.0), ("AMRG", graph(NH2, closure="N"), 6.2, 6.2),
        ("HMTMR", graph(NME), 9.5, 9.5), ("AMTMR", graph(NME, closure="N"), 7.8, 7.8),
        ("HySOx", graph(NME, closure="S"), "<2", "<2"), ("HMRB", graph(NET), 9.2, 9.3),
        ("HMR6G", graph(arm("NCC", after="C")), 10.1, 10.0),
        ("HMR101", graph(julolidine(4, 5), julolidine(6, 7)), 10.8, 10.8),
    ]
    for name, smi, absorbance, fluorescence in uno:
        status = ("complex" if name in ("HMSiR", "AMSiR") else
                  "open-control" if name == "2MeSiR" else "upper-bound" if absorbance == "<2" else "exact")
        reason = {"complex": "Bell-shaped signal in Fig. 2a.", "upper-bound": "Only an upper bound.",
                  "open-control": "Ortho methyl cannot cyclize; only open cation supplied."}.get(status, "")
        for value, method, loc in [(absorbance, "Absorbance", "Fig. 2a–b"),
                                   (fluorescence, "Fluorescence", "SI Fig. S2")]:
            add("Uno", name, smi, value, loc, method=method, status=status, reason=reason,
                form="open cation" if name == "2MeSiR" else "closed")
            rows[-1]["substitution"] = {
                "HMAcRG": "Hydroxymethyl; acetamido/amino", "AMAcRG": "Aminomethyl; acetamido/amino",
                "HMSiR": "Hydroxymethyl; para CO2H; dimethylamino", "AMSiR": "Aminomethyl; para CO2H; dimethylamino",
                "MMSiR": "Mercaptomethyl; para CO2H; dimethylamino", "2MeSiR": "Ortho methyl; para CO2H; dimethylamino",
                "HMRG": "Hydroxymethyl; two amino", "AMRG": "Aminomethyl; two amino",
                "HMTMR": "Hydroxymethyl; two dimethylamino", "AMTMR": "Aminomethyl; two dimethylamino",
                "HySOx": "Mercaptomethyl; two dimethylamino", "HMRB": "Hydroxymethyl; two diethylamino",
                "HMR6G": "Hydroxymethyl; two ethylamino; adjacent aryl methyls toward spiro center",
                "HMR101": "Hydroxymethyl; two julolidines"}[name]

    def tyson_fused(sub, chain):
        return graph(fused(sub, chain, "4"), fused(sub, chain, "5"), bridge=SI, pendant=ACID)

    thq = tyson_fused("C", "CCC")
    tyson = [("HMSiR", acid_si, 6.0, "3a,e; SI Scheme S1"),
             ("HMSiRindol", tyson_fused("C", "CC"), 6.4, "2; SI Scheme S2"),
             ("HMSiRjulol", graph(julolidine(4, 5), julolidine(6, 7), SI, pendant=ACID), 9.0, "2; SI Scheme S3"),
             ("HMSiRTHQ", thq, 6.9, "2, 3b,e; SI Scheme S4"),
             ("HMSiR2-FlEt", graph(arm("N(C)CCF"), bridge=SI, pendant=ACID), 5.0, "3a,e; SI Scheme S6"),
             ("Yale676sb", tyson_fused("CCF", "CCC"), 5.9, "3b,e; SI Scheme S5"),
             ("Cal664sb", tyson_fused("CC(F)F", "CCC"), 3.2, "3b,e; SI Scheme S7")]
    for name, smi, value, fig in tyson:
        complex_curve = name in ("HMSiR", "HMSiR2-FlEt")
        add("Tyson", name, smi, value, "Fig. " + fig + "; SI §D, Fig. S2",
            status="complex" if complex_curve else "exact",
            reason="Upper midpoint of two-group fit." if complex_curve else "")
        rows[-1]["substitution"] = "Para CO2H and ortho hydroxymethyl; " + {
            "HMSiR": "dimethylamino", "HMSiRindol": "N-methylindoline",
            "HMSiRjulol": "julolidine", "HMSiRTHQ": "N-methyltetrahydroquinoline",
            "HMSiR2-FlEt": "N-methyl-N-(2-fluoroethyl)amino", "Yale676sb": "N-(2-fluoroethyl)tetrahydroquinoline",
            "Cal664sb": "N-(2,2-difluoroethyl)tetrahydroquinoline"}[name] + " on both arms"
    add("Tyson", "HMSiR", acid_si, 5.8, "Fig. 2a comparator; compare Uno Fig. 2",
        status="comparator", evidence="alternate report", origin="Tyson Fig. 2; same number in Uno",
        reason="Independence of the Fig. 2 value is unclear; Fig. 3 reports 6.0 for a complex curve.")

    for name, substituent, position, value, fig in [("3FHMRG", "F", 3, 8.0, 17),
                                                   ("3MHMRG", "C", 3, 6.6, 16),
                                                   ("3CF3HMRG", "C(F)(F)F", 3, 5.3, 18),
                                                   ("5MHMRG", "C", 5, 8.2, 19)]:
        pendant = (f"c3cccc({substituent})c3" if position == 3 else f"c3cc({substituent})ccc3")
        add("Tachibana", name, graph(NH2, pendant=pendant), value, f"Table 2; SI Fig. {fig}")
        rows[-1]["substitution"] = f"Pendant aryl {position}′-{substituent}; amino donors on both arms"

    for name, smi, value in [("HMRG", graph(NH2), 8.1), ("HMDiMeR", graph(NME, NH2), 8.9),
                             ("HMDiEtR", graph(NET, NH2), 9.3)]:
        add("Sakabe", name, smi, value, "Main text, p. 410–411; SI Fig. S1")

    for name, smi, values, loc, cosolvent in [
        ("HMSiR620h", graph(NME, NH2, SI), "3.3, 5.0", "S3", "<0.12% acetonitrile and 0.68% DMSO"),
        ("HMJSiR", graph(julolidine(4, 5), NH2, SI), "1.3, 6.6", "S4", "0.5% DMSO")]:
        add("Iwatate 2018", name, smi, values, f"SI Table S1, Fig. {loc}", status="complex",
            buffer=f"0.2 M sodium phosphate, {cosolvent}", reason="Two reported transitions.")

    return rows


def quantum_rows():
    """All 17 rows, including table-only values. None become new training measurements."""
    plain = [(u, graph(a, closure=c), m, q) for u, a, c, m, q in [
        ("HMRG", NH2, "O", 8.1, 7.9), ("AMRG", NH2, "N", 6.2, 6.2),
        ("HMTMR", NME, "O", 9.5, 9.5), ("AMTMR", NME, "N", 7.8, 8.1),
        ("HMRB", NET, "O", 9.2, 9.3), ("AMRB", NET, "N", 8.2, 8.1)]]
    plain += [("HMSiR", graph(NME, bridge=SI), 5.7, 6.2),
              ("AMSiR", graph(NME, bridge=SI, closure="N"), 4.2, 4.8)]
    asym = [("HMDiMeR", NME, "", 8.9, 9.2), ("HMDiEtR", NET, "", 9.3, 8.8),
            ("HMJR", julolidine(4, 5), "", 10.3, 10.7), ("HMDiMeFR", NME, "F", 8.2, 8.3),
            ("HMDiMeCR", NME, "Cl", 7.7, 7.3), ("HMJFR", julolidine(4, 5), "F", 9.8, 9.8),
            ("HMJCR", julolidine(4, 5), "Cl", 9.1, 8.7)]
    plain += [(n, graph(a, arm("N", after=z)), m, q) for n, a, z, m, q in asym]
    plain += [("HMSiR620h", graph(NME, NH2, SI), 5.0, 4.8),
              ("HMJSiR", graph(julolidine(4, 5), NH2, SI), 6.6, 7.1)]
    return plain
