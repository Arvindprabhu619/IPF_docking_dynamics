#!/usr/bin/env python3
"""
05_imppat_screen.py (v2.0.0)
Strict Two-Stage Disease-Driven Phytochemical Screening:
  Stage 1: Filter all IMPPAT 2.0 phytochemicals (with valid Canonical SMILES)
           for strict multi-rule drug-likeness (Lipinski 0 viol, Ghose 0 viol,
           Veber, Egan, GSK 4/400) using RDKit.
  Stage 2: Screen ONLY the strictly drug-like subset against the 12 upstream
           M2 macrophage regulators (5 primary + 7 secondary) using
           experimentally curated IMPPAT 2.0 compound-gene associations.
  Stage 3: Apply medicinal chemistry / safety triage (excluding endogenous
           steroid hormones, industrial xenoestrogens, and overt cytotoxic/
           mitochondrial poisons) to define the Top 10 docking candidates.
"""
from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors

ROOT = Path(__file__).resolve().parents[1]
COMP_DIR = ROOT / "data" / "compounds"
TAB_DIR = ROOT / "results" / "tables"
TAB_DIR.mkdir(parents=True, exist_ok=True)

PRIMARY_REGULATORS = {"NFKB1", "PPARG", "HIF1A", "CCR2", "STAT3"}
SECONDARY_REGULATORS = {"CSF1R", "MMP12", "TNF", "IL6", "MRC1", "CHIT1", "LPL"}
ALL_REGULATORS = PRIMARY_REGULATORS | SECONDARY_REGULATORS

NON_LEAD_EXCLUSIONS = {
    "progesterone": "Endogenous steroid hormone; not a plant-derived lead",
    "estrone": "Endogenous estrogenic steroid hormone; not a plant-derived lead",
    "estradiol": "Endogenous estrogenic steroid hormone; not a plant-derived lead",
    "testosterone": "Endogenous androgenic steroid hormone; not a plant-derived lead",
    "4-nonylphenol": "Endocrine-disrupting industrial xenoestrogen",
    "catechol": "Non-specific reactive PAINS fragment",
    "camptothecin": "Cytotoxic DNA topoisomerase I inhibitor (chemotherapeutic poison)",
    "rotenone": "Mitochondrial complex I neurotoxin/pesticide",
}


def evaluate_druglikeness(smiles: str):
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        return None

    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    rotb = Lipinski.NumRotatableBonds(mol)
    mr = Crippen.MolMR(mol)
    n_atoms = mol.GetNumAtoms()

    lip_viol = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    ghose_ok = (160 <= mw <= 480) and (-0.4 <= logp <= 5.6) and (40 <= mr <= 130) and (20 <= n_atoms <= 70)
    veber_ok = (rotb <= 10) and (tpsa <= 140)
    egan_ok = (logp <= 5.88) and (tpsa <= 131.6)
    gsk_ok = (mw <= 400) and (logp <= 4.0)

    strict_pass = (lip_viol == 0) and ghose_ok and veber_ok and egan_ok and gsk_ok
    return {
        "MW": round(mw, 2),
        "LogP": round(logp, 2),
        "HBD": int(hbd),
        "HBA": int(hba),
        "TPSA": round(tpsa, 2),
        "RotBonds": int(rotb),
        "Lipinski_Viol": int(lip_viol),
        "Strict_Druglike": bool(strict_pass),
    }


def main():
    chem_df = pd.read_csv(COMP_DIR / "Chemical_Information_IMPPAT_Phytochemicals.tsv", sep="\t", dtype=str)
    targ_df = pd.read_csv(COMP_DIR / "Target_IMPPAT_Phytochemicals.tsv", sep="\t", dtype=str)
    plant_df = pd.read_csv(COMP_DIR / "IMPPAT_Phytochemical_Plant_Association.tsv", sep="\t", dtype=str)

    chem_df = chem_df.drop_duplicates(subset=["IMPPAT_Phytochemical_identifier"]).copy()
    props = chem_df["Canonical_SMILES"].apply(evaluate_druglikeness)
    valid_mask = props.notnull()
    chem_full = pd.concat([chem_df[valid_mask].reset_index(drop=True),
                           pd.DataFrame(list(props[valid_mask]))], axis=1)

    druglike_df = chem_full[chem_full["Strict_Druglike"]].copy()
    druglike_df.to_csv(TAB_DIR / "imppat_druglike_library.tsv", sep="\t", index=False)

    targ_df["HGNC_Symbol"] = targ_df["HGNC_Symbol"].str.upper().str.strip()
    dl_ids = set(druglike_df["IMPPAT_Phytochemical_identifier"])

    sub_targ = targ_df[
        targ_df["IMPPAT_Phytochemical_identifier"].isin(dl_ids) &
        targ_df["HGNC_Symbol"].isin(ALL_REGULATORS)
    ].copy()

    grouped = (
        sub_targ.groupby("IMPPAT_Phytochemical_identifier")["HGNC_Symbol"]
        .apply(lambda s: sorted(set(s)))
        .reset_index(name="Regulators_List")
    )
    grouped["Total_Targets"] = grouped["Regulators_List"].apply(len)
    grouped["Primary_Targets"] = grouped["Regulators_List"].apply(lambda lst: len(set(lst) & PRIMARY_REGULATORS))
    grouped["Regulators_Hit"] = grouped["Regulators_List"].apply(lambda lst: ", ".join(lst))

    hits = grouped.merge(druglike_df, on="IMPPAT_Phytochemical_identifier", how="inner")
    plant_col = [c for c in plant_df.columns if "plant" in c.lower() or "species" in c.lower() or "name" in c.lower()][0]
    plants_grouped = (
        plant_df.groupby("IMPPAT_Phytochemical_identifier")[plant_col]
        .apply(lambda s: "; ".join(sorted(set(s.dropna()))[:5]))
        .reset_index(name="Source_Plants_Top5")
    )
    hits = hits.merge(plants_grouped, on="IMPPAT_Phytochemical_identifier", how="left")
    hits["Exclusion_Note"] = hits["Phytochemical name_standardised"].str.lower().str.strip().map(
        lambda n: NON_LEAD_EXCLUSIONS.get(n, "")
    )
    hits = hits.sort_values(
        by=["Total_Targets", "Primary_Targets", "MW"],
        ascending=[False, False, True]
    ).reset_index(drop=True)

    hits.to_csv(TAB_DIR / "druglike_m2_regulator_hits_all.tsv", sep="\t", index=False)

    eligible = hits[(hits["Total_Targets"] >= 2) & (hits["Exclusion_Note"] == "")].copy().reset_index(drop=True)
    eligible.insert(0, "Rank", range(1, len(eligible) + 1))
    top10 = eligible.head(10)
    top10.to_csv(TAB_DIR / "top10_docking_candidates.tsv", sep="\t", index=False)

    print("=== Excluded Multi-Target Compounds (Safety / Non-Plant Triage) ===")
    excl = hits[(hits["Total_Targets"] >= 2) & (hits["Exclusion_Note"] != "")]
    print(excl[["Phytochemical name_standardised", "Total_Targets", "Exclusion_Note"]].to_string(index=False))

    print("\n=== Final Top 10 Drug-Like Phytochemicals Selected for Docking ===")
    cols = ["Rank", "IMPPAT_Phytochemical_identifier", "Phytochemical name_standardised",
            "Pubchem_CID", "Total_Targets", "Primary_Targets", "MW", "LogP", "Regulators_Hit"]
    print(top10[cols].to_string(index=False))


if __name__ == "__main__":
    main()
