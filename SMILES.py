#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, MACCSkeys
from rdkit.DataStructs import ConvertToNumpyArray

try:
    import selfies as sf
except ImportError as e:
    raise ImportError("SELFIES is required. Install with: pip install selfies") from e


OUTDIR = Path("/fs/scratch/PAS2942/Alejandro/RNA/drugs/SMILES")
OUTDIR.mkdir(parents=True, exist_ok=True)

FP_RADIUS = 2
FP_BITS = 512
MAX_SELFIES_LEN = 160
PUBCHEM_SLEEP = 0.25

PAD = "<PAD>"
BOS = "<BOS>"
EOS = "<EOS>"
UNK = "<UNK>"

RAW_TREATMENTS = [
    "AE37 Peptide/GM-CSF Vaccine",
    "Anastrozole",
    "Bevacizumab",
    "Brachytherapy, High Dose",
    "Brachytherapy, NOS",
    "Capecitabine",
    "Carboplatin",
    "Chemotherapy",
    "Cisplatin",
    "Clinical Trial",
    "Clinical Trial Agent",
    "Clodronate Disodium",
    "Clodronic Acid",
    "Cyclophosphamide",
    "Denosumab",
    "Docetaxel",
    "Doxorubicin",
    "Doxorubicin Hydrochloride",
    "Epirubicin",
    "Etoposide",
    "Everolimus",
    "Exemestane",
    "Fluorouracil",
    "Fulvestrant",
    "Gemcitabine",
    "Gemcitabine Hydrochloride",
    "Goserelin",
    "Goserelin Acetate",
    "Hormone Therapy",
    "Ibandronate Sodium",
    "Ifosfamide",
    "Ixabepilone",
    "Lapatinib",
    "Letrozole",
    "Leuprolide",
    "Leuprolide Acetate",
    "Megestrol Acetate",
    "Mesna",
    "Metformin",
    "Methotrexate",
    "Mitomycin",
    "Mitoxantrone",
    "Nab-paclitaxel",
    "Nelipepimut-S",
    "Not Reported",
    "Paclitaxel",
    "Palonosetron Hydrochloride",
    "Pamidronate Disodium",
    "Pamidronic Acid",
    "Pegfilgrastim",
    "Pegylated Liposomal Doxorubicin Hydrochloride",
    "Pemetrexed",
    "Pharmaceutical Therapy, NOS",
    "Prednisone",
    "Radiation Therapy, NOS",
    "Radiation, External Beam",
    "Radiation, Implants",
    "Radiation, Radioisotope",
    "Radiation, Stereotactic/Gamma Knife/SRS",
    "Rituximab",
    "Surgery, NOS",
    "Tamoxifen",
    "Tamoxifen Citrate",
    "Taxane Compound",
    "Tesetaxel",
    "Toremifene Citrate",
    "Trabectedin",
    "Trastuzumab",
    "Triptorelin",
    "Unknown",
    "Vinblastine",
    "Vincristine",
    "Vinorelbine",
    "Vinorelbine Tartrate",
    "Zoledronic Acid",
]

EXCLUDE = {
    "AE37 Peptide/GM-CSF Vaccine": "peptide vaccine / immunotherapy label",
    "Bevacizumab": "monoclonal antibody biologic",
    "Brachytherapy, High Dose": "procedure",
    "Brachytherapy, NOS": "procedure",
    "Chemotherapy": "generic treatment class",
    "Clinical Trial": "administrative label",
    "Clinical Trial Agent": "ambiguous label",
    "Denosumab": "monoclonal antibody biologic",
    "Hormone Therapy": "generic treatment class",
    "Nelipepimut-S": "peptide vaccine",
    "Not Reported": "missing label",
    "Pegfilgrastim": "protein biologic",
    "Pharmaceutical Therapy, NOS": "generic treatment class",
    "Radiation Therapy, NOS": "procedure",
    "Radiation, External Beam": "procedure",
    "Radiation, Implants": "procedure",
    "Radiation, Radioisotope": "procedure",
    "Radiation, Stereotactic/Gamma Knife/SRS": "procedure",
    "Rituximab": "monoclonal antibody biologic",
    "Surgery, NOS": "procedure",
    "Taxane Compound": "class label, not a compound",
    "Trastuzumab": "monoclonal antibody biologic",
    "Unknown": "missing label",
}

NORMALIZE = {
    "Clodronate Disodium": "Clodronic Acid",
    "Doxorubicin Hydrochloride": "Doxorubicin",
    "Gemcitabine Hydrochloride": "Gemcitabine",
    "Goserelin Acetate": "Goserelin",
    "Ibandronate Sodium": "Ibandronic Acid",
    "Leuprolide Acetate": "Leuprolide",
    "Nab-paclitaxel": "Paclitaxel",
    "Palonosetron Hydrochloride": "Palonosetron",
    "Pamidronate Disodium": "Pamidronic Acid",
    "Pegylated Liposomal Doxorubicin Hydrochloride": "Doxorubicin",
    "Tamoxifen Citrate": "Tamoxifen",
    "Toremifene Citrate": "Toremifene",
    "Vinorelbine Tartrate": "Vinorelbine",
}

QUERY_OVERRIDES = {
    "Fluorouracil": "5-Fluorouracil",
    "Ibandronic Acid": "Ibandronic acid",
    "Pamidronic Acid": "Pamidronic acid",
    "Zoledronic Acid": "Zoledronic acid",
    "Megestrol Acetate": "Megestrol acetate",
}


def safe_name(x: str) -> str:
    x = re.sub(r"[^\w\-.]+", "_", str(x).strip())
    x = re.sub(r"_+", "_", x).strip("_")
    return x


def write_json(path: Path, obj: Any) -> None:
    with open(path, "w") as f:
        json.dump(obj, f, indent=2)


def write_text(path: Path, text: str) -> None:
    with open(path, "w") as f:
        f.write(str(text).strip() + "\n")


def pubchem_lookup(name: str) -> Dict[str, Any]:
    query = QUERY_OVERRIDES.get(name, name)
    encoded = urllib.parse.quote(query)
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
        f"{encoded}/property/CanonicalSMILES,IsomericSMILES,InChIKey,IUPACName/JSON"
    )

    time.sleep(PUBCHEM_SLEEP)

    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            data = json.loads(response.read().decode("utf-8"))

        props = data["PropertyTable"]["Properties"][0]

        return {
            "status": "ok",
            "query": query,
            "canonical_smiles_pubchem": props.get("CanonicalSMILES", ""),
            "isomeric_smiles_pubchem": props.get("IsomericSMILES", ""),
            "inchikey": props.get("InChIKey", ""),
            "iupac_name": props.get("IUPACName", ""),
            "error": "",
        }

    except Exception as e:
        return {
            "status": "pubchem_failed",
            "query": query,
            "canonical_smiles_pubchem": "",
            "isomeric_smiles_pubchem": "",
            "inchikey": "",
            "iupac_name": "",
            "error": str(e),
        }


def canonicalize_smiles(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None

    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None, None

    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    return canonical, mol


def smiles_to_selfies(smiles: str):
    try:
        return sf.encoder(smiles)
    except Exception:
        return None


def selfies_tokens(selfies: str) -> List[str]:
    try:
        return list(sf.split_selfies(selfies))
    except Exception:
        return []


def bounded_tokens(tokens: List[str]) -> List[str]:
    toks = [BOS] + list(tokens) + [EOS]
    if len(toks) > MAX_SELFIES_LEN:
        toks = toks[:MAX_SELFIES_LEN]
        toks[-1] = EOS
    return toks


def build_vocab(token_lists: List[List[str]]):
    vocab = {PAD: 0, BOS: 1, EOS: 2, UNK: 3}

    for toks in token_lists:
        for tok in toks:
            if tok not in vocab:
                vocab[tok] = len(vocab)

    inv_vocab = {str(v): k for k, v in vocab.items()}
    return vocab, inv_vocab


def tokens_to_ids(tokens: List[str], vocab: Dict[str, int]) -> np.ndarray:
    ids = [vocab.get(t, vocab[UNK]) for t in tokens]

    if len(ids) < MAX_SELFIES_LEN:
        ids = ids + [vocab[PAD]] * (MAX_SELFIES_LEN - len(ids))

    return np.asarray(ids[:MAX_SELFIES_LEN], dtype=np.int64)


def token_counts(tokens: List[str], vocab: Dict[str, int]) -> np.ndarray:
    arr = np.zeros((len(vocab),), dtype=np.float32)

    for tok in tokens:
        arr[vocab.get(tok, vocab[UNK])] += 1.0

    s = arr.sum()
    if s > 0:
        arr /= s

    return arr


def morgan_fp(mol) -> np.ndarray:
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, FP_RADIUS, nBits=FP_BITS)
    arr = np.zeros((FP_BITS,), dtype=np.float32)
    ConvertToNumpyArray(fp, arr)
    return arr


def maccs_fp(mol) -> np.ndarray:
    fp = MACCSkeys.GenMACCSKeys(mol)
    arr = np.zeros((fp.GetNumBits(),), dtype=np.float32)
    ConvertToNumpyArray(fp, arr)
    return arr


def descriptor_vec(mol) -> np.ndarray:
    vals = [
        Descriptors.MolWt(mol),
        Descriptors.MolLogP(mol),
        Descriptors.TPSA(mol),
        Descriptors.NumHDonors(mol),
        Descriptors.NumHAcceptors(mol),
        Descriptors.NumRotatableBonds(mol),
        Descriptors.RingCount(mol),
        Descriptors.FractionCSP3(mol),
        rdMolDescriptors.CalcNumAromaticRings(mol),
        rdMolDescriptors.CalcNumAliphaticRings(mol),
        rdMolDescriptors.CalcNumHBA(mol),
        rdMolDescriptors.CalcNumHBD(mol),
        rdMolDescriptors.CalcExactMolWt(mol),
    ]
    return np.asarray(vals, dtype=np.float32)


def main() -> None:
    lookup_rows = []
    skipped_rows = []
    valid_rows = []

    for raw_name in RAW_TREATMENTS:
        if raw_name in EXCLUDE:
            skipped_rows.append({
                "raw_name": raw_name,
                "normalized_name": NORMALIZE.get(raw_name, raw_name),
                "reason": EXCLUDE[raw_name],
            })
            continue

        normalized_name = NORMALIZE.get(raw_name, raw_name)

        lookup = pubchem_lookup(normalized_name)
        lookup_rows.append({
            "raw_name": raw_name,
            "normalized_name": normalized_name,
            **lookup,
        })

        if lookup["status"] != "ok":
            skipped_rows.append({
                "raw_name": raw_name,
                "normalized_name": normalized_name,
                "reason": "pubchem_lookup_failed",
                "pubchem_query": lookup["query"],
                "error": lookup["error"],
            })
            continue

        smiles_candidate = lookup["isomeric_smiles_pubchem"] or lookup["canonical_smiles_pubchem"]

        if not smiles_candidate:
            skipped_rows.append({
                "raw_name": raw_name,
                "normalized_name": normalized_name,
                "reason": "pubchem_returned_no_smiles",
                "pubchem_query": lookup["query"],
            })
            continue

        smiles, mol = canonicalize_smiles(smiles_candidate)

        if smiles is None or mol is None:
            skipped_rows.append({
                "raw_name": raw_name,
                "normalized_name": normalized_name,
                "reason": "rdkit_parse_failed",
                "pubchem_query": lookup["query"],
                "smiles_candidate": smiles_candidate,
            })
            continue

        selfies = smiles_to_selfies(smiles)

        if selfies is None:
            skipped_rows.append({
                "raw_name": raw_name,
                "normalized_name": normalized_name,
                "reason": "selfies_encoding_failed",
                "pubchem_query": lookup["query"],
                "smiles": smiles,
            })
            continue

        toks = bounded_tokens(selfies_tokens(selfies))

        if len(toks) <= 2:
            skipped_rows.append({
                "raw_name": raw_name,
                "normalized_name": normalized_name,
                "reason": "empty_selfies_tokens",
                "pubchem_query": lookup["query"],
                "smiles": smiles,
            })
            continue

        valid_rows.append({
            "raw_name": raw_name,
            "normalized_name": normalized_name,
            "pubchem_query": lookup["query"],
            "inchikey": lookup["inchikey"],
            "iupac_name": lookup["iupac_name"],
            "canonical_smiles_pubchem": lookup["canonical_smiles_pubchem"],
            "isomeric_smiles_pubchem": lookup["isomeric_smiles_pubchem"],
            "smiles": smiles,
            "selfies": selfies,
            "selfies_tokens": toks,
            "n_selfies_tokens_unpadded": len(toks),
        })

    pd.DataFrame(lookup_rows).to_csv(OUTDIR / "pubchem_lookup_report.tsv", sep="\t", index=False)
    pd.DataFrame(skipped_rows).to_csv(OUTDIR / "smiles_selfies_skipped.tsv", sep="\t", index=False)

    if not valid_rows:
        print("[FAIL] No valid molecules were produced.")
        print(f"[FAIL] PubChem report: {OUTDIR / 'pubchem_lookup_report.tsv'}")
        print(f"[FAIL] Skipped report: {OUTDIR / 'smiles_selfies_skipped.tsv'}")
        raise RuntimeError("No valid molecules. Check PubChem/network failures in pubchem_lookup_report.tsv.")

    vocab, inv_vocab = build_vocab([r["selfies_tokens"] for r in valid_rows])

    index_rows = []

    for r in valid_rows:
        raw_name = r["raw_name"]
        stem = safe_name(raw_name)
        smiles = r["smiles"]
        selfies = r["selfies"]
        toks = r["selfies_tokens"]

        smiles, mol = canonicalize_smiles(smiles)
        if smiles is None or mol is None:
            continue

        morgan = morgan_fp(mol)
        maccs = maccs_fp(mol)
        desc = descriptor_vec(mol)
        ids = tokens_to_ids(toks, vocab)
        counts = token_counts(toks, vocab)

        smiles_path = OUTDIR / f"{stem}.smiles.txt"
        selfies_path = OUTDIR / f"{stem}.selfies.txt"
        token_json_path = OUTDIR / f"{stem}.selfies_tokens.json"
        token_ids_path = OUTDIR / f"{stem}.selfies_ids.npy"
        token_counts_path = OUTDIR / f"{stem}.selfies_token_counts.npy"
        morgan_path = OUTDIR / f"{stem}.morgan.npy"
        maccs_path = OUTDIR / f"{stem}.maccs.npy"
        desc_path = OUTDIR / f"{stem}.desc.npy"
        meta_path = OUTDIR / f"{stem}.json"

        np.save(morgan_path, morgan)
        np.save(maccs_path, maccs)
        np.save(desc_path, desc)
        np.save(token_ids_path, ids)
        np.save(token_counts_path, counts)

        write_text(smiles_path, smiles)
        write_text(selfies_path, selfies)

        write_json(token_json_path, {
            "raw_name": raw_name,
            "normalized_name": r["normalized_name"],
            "tokens": toks,
            "token_ids": ids.tolist(),
            "max_selfies_len": MAX_SELFIES_LEN,
        })

        write_json(meta_path, {
            "raw_name": raw_name,
            "normalized_name": r["normalized_name"],
            "pubchem_query": r["pubchem_query"],
            "inchikey": r["inchikey"],
            "iupac_name": r["iupac_name"],
            "canonical_smiles_pubchem": r["canonical_smiles_pubchem"],
            "isomeric_smiles_pubchem": r["isomeric_smiles_pubchem"],
            "smiles": smiles,
            "selfies": selfies,
            "morgan_radius": FP_RADIUS,
            "morgan_bits": int(morgan.shape[0]),
            "maccs_bits": int(maccs.shape[0]),
            "descriptor_dim": int(desc.shape[0]),
            "selfies_vocab_size": int(len(vocab)),
            "selfies_max_len": MAX_SELFIES_LEN,
            "smiles_file": str(smiles_path),
            "selfies_file": str(selfies_path),
            "selfies_ids_file": str(token_ids_path),
            "selfies_token_counts_file": str(token_counts_path),
            "morgan_file": str(morgan_path),
            "maccs_file": str(maccs_path),
            "descriptor_file": str(desc_path),
            "metadata_file": str(meta_path),
        })

        index_rows.append({
            "raw_name": raw_name,
            "normalized_name": r["normalized_name"],
            "pubchem_query": r["pubchem_query"],
            "inchikey": r["inchikey"],
            "iupac_name": r["iupac_name"],
            "smiles": smiles,
            "selfies": selfies,
            "selfies_unpadded_len": len(toks),
            "selfies_vocab_size": len(vocab),
            "selfies_max_len": MAX_SELFIES_LEN,
            "morgan_dim": int(morgan.shape[0]),
            "maccs_dim": int(maccs.shape[0]),
            "descriptor_dim": int(desc.shape[0]),
            "smiles_file": str(smiles_path),
            "selfies_file": str(selfies_path),
            "selfies_tokens_file": str(token_json_path),
            "selfies_ids_file": str(token_ids_path),
            "selfies_token_counts_file": str(token_counts_path),
            "morgan_file": str(morgan_path),
            "maccs_file": str(maccs_path),
            "descriptor_file": str(desc_path),
            "metadata_file": str(meta_path),
        })

    pd.DataFrame(index_rows).to_csv(OUTDIR / "smiles_selfies_index.tsv", sep="\t", index=False)

    write_json(OUTDIR / "selfies_vocab.json", vocab)
    write_json(OUTDIR / "selfies_inv_vocab.json", inv_vocab)

    write_json(OUTDIR / "smiles_selfies_config.json", {
        "source": "PubChem PUG-REST",
        "output_dir": str(OUTDIR),
        "fp_radius": FP_RADIUS,
        "fp_bits": FP_BITS,
        "max_selfies_len": MAX_SELFIES_LEN,
        "vocab_size": len(vocab),
        "n_raw_treatments": len(RAW_TREATMENTS),
        "n_valid_molecules": len(index_rows),
        "n_skipped": len(skipped_rows),
    })

    print(f"[OK] Wrote {len(index_rows)} valid molecules")
    print(f"[OK] Output directory: {OUTDIR}")
    print(f"[OK] Index: {OUTDIR / 'smiles_selfies_index.tsv'}")
    print(f"[OK] Lookup report: {OUTDIR / 'pubchem_lookup_report.tsv'}")
    print(f"[OK] Skipped report: {OUTDIR / 'smiles_selfies_skipped.tsv'}")


if __name__ == "__main__":
    main()
