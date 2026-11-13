"""Capture full-parent public identity lookups. HTTP failures remain distinguishable from no matches."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import json
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from rdkit import Chem


BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
PMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?"
ROOT = Path(__file__).resolve().parents[1]


def queries():
    spec = json.loads((ROOT / "validation/structure_pk/registration.json").read_text())
    for row in spec["candidates"]:
        mol = Chem.MolFromSmiles(row["closedSmiles"])
        key = Chem.MolToInchiKey(mol)
        yield row["id"], "full closed InChIKey", BASE + f"inchikey/{key}/cids/JSON"
        spiro, oxygen = mol.GetSubstructMatch(Chem.MolFromSmarts("[C;X4;R2]-[O;R]"))
        opened = Chem.RWMol(mol)
        opened.RemoveBond(spiro, oxygen)
        opened.GetAtomWithIdx(spiro).SetFormalCharge(1)
        Chem.SanitizeMol(opened)
        for label, parent in (("closed", mol), ("open cation", opened)):
            query = urlencode({"smiles": Chem.MolToSmiles(parent), "identity_type": "same_connectivity"})
            yield row["id"], f"full {label} connectivity", BASE + "fastidentity/smiles/cids/JSON?" + query
        word = "difluoromethyl" if "CHF2" in row["id"] else "fluoromethyl"
        for label, terms in [("literature full key", '"' + key + '"'),
                             ("literature chemical text", f'("{word}azetidine" OR "({word})azetidine" OR "{word} azetidine") AND (rhodamine OR SiR OR fluorophore)')]:
            yield row["id"], label, PMC + urlencode({"query": terms, "format": "json", "pageSize": 50})
        yield row["id"], "building block control", BASE + f"cid/{row['buildingBlockCID']}/property/Title,InChIKey/JSON"
    yield "control", "identity service positive control", BASE + "fastidentity/smiles/cids/JSON?smiles=CCO&identity_type=same_connectivity"


def lookup(item):
    dye, kind, url = item
    result = dict(dye=dye, query=kind, url=url)
    try:
        with urlopen(Request(url, headers={"User-Agent": "BlinkAtlas/1.0"}), timeout=35) as resp:
            body = json.load(resp)
            result["http"] = resp.status
        if "hitCount" in body:
            result.update(count=body["hitCount"], articles=[
                {"title": r.get("title"), "doi": r.get("doi"), "id": r.get("id")}
                for r in body["resultList"]["result"]])
        elif "IdentifierList" in body:
            result["cids"] = body["IdentifierList"]["CID"]
        else:
            result["response"] = body
    except HTTPError as exc:
        result.update(http=exc.code, response=exc.read(1600).decode("utf-8", "replace"))
    except (OSError, ValueError) as exc:
        result["failure"] = str(exc)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with ThreadPoolExecutor(max_workers=2) as executor:
        results = list(executor.map(lookup, queries()))
    args.output.write_text(json.dumps({"retrievedUTC": datetime.now(timezone.utc).date().isoformat(),
        "queries": results}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
