#!/usr/bin/env python3
"""
06_prepare_docking.py (v2.0.0)
Prepares clean receptor PDBQT structures, exact binding-site grid boxes,
reference crystal ligand PDBs (for redocking RMSD validation), and 3D
MMFF94-minimized ligand PDBQT files for:
  - Top 10 strictly drug-like IMPPAT phytochemicals (salt-stripped)
  - Clinical reference drug: Pirfenidone (CID_40632)
  - Co-crystallized native positive controls: BRL (2PRG), KQV (6NJS), F7N (6GPS)
"""
import json
import subprocess
import urllib.request
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize

ROOT = Path(__file__).resolve().parents[1]
TAB_DIR = ROOT / "results" / "tables"
DOCK_DIR = ROOT / "results" / "docking"
REC_DIR = DOCK_DIR / "receptors"
LIG_DIR = DOCK_DIR / "ligands"
REF_DIR = DOCK_DIR / "crystal_ligands"
for d in [REC_DIR, LIG_DIR, REF_DIR]:
    d.mkdir(parents=True, exist_ok=True)

TARGETS = {
    "PPARG": {
        "pdb": "2PRG",
        "chain": "A",
        "native_resname": "BRL",
        "description": "Human PPARG LBD (2.30 A) co-crystallized with rosiglitazone (BRL)",
        "box_size": [22.0, 22.0, 22.0],
    },
    "STAT3": {
        "pdb": "6NJS",
        "chain": "A",
        "native_resname": "KQV",
        "description": "Human STAT3 SH2 domain (2.70 A) co-crystallized with SD-36 (KQV)",
        "box_size": [22.0, 22.0, 22.0],
    },
    "CCR2": {
        "pdb": "6GPS",
        "chain": "A",
        "native_resname": "F7N",
        "description": "Human CCR2A (3.30 A) co-crystallized with orthosteric antagonist MK-0812 (F7N)",
        "box_size": [22.0, 22.0, 22.0],
    },
    "NFKB1": {
        "pdb": "1SVC",
        "chain": "P",
        "native_resname": None,
        "description": "Human NF-kB p50 Rel-homology domain (2.60 A); DNA-binding groove pocket",
        "pocket_residues": ["TYR60", "CYS62", "GLU63", "HIS67", "LYS147", "LYS148", "ASP209", "LYS244", "GLN277", "ARG308"],
        "box_size": [24.0, 24.0, 24.0],
    },
    "HIF1A": {
        "pdb": "4ZPR",
        "chain": "B",
        "native_resname": None,
        "description": "Murine HIF-1a:ARNT PAS-B heterodimer (3.90 A, Chain B HIF-1a); PAS-B pocket",
        "pocket_residues": ["TYR276", "GLN299", "HIS313", "MET315", "THR327", "HIS374"],
        "box_size": [24.0, 24.0, 24.0],
    },
}

EXTRA_LIGANDS = [
    {"id": "Pirfenidone", "name": "Pirfenidone (Clinical Control)", "smiles": "CC1=CC=C(C=C1)N2C=CC=CC2=O", "role": "Clinical_Control"},
    {"id": "Native_BRL",  "name": "Rosiglitazone (BRL, 2PRG Native)", "smiles": "CN(CCOc1ccc(CC2SC(=O)NC2=O)cc1)c1ccccn1", "role": "Native_Control"},
    {"id": "Native_KQV",  "name": "SD-36 (KQV, 6NJS Native)", "smiles": "CC(C)(C)c1ccc(NC(=O)C2CCCN2C(=O)C(CC(=O)O)NC(=O)c2ccc3ccccc3c2)cc1", "role": "Native_Control"},
    {"id": "Native_F7N",  "name": "MK-0812 (F7N, 6GPS Native)", "smiles": "CO[C@@H]1CCOC[C@H]1N[C@H]1CC[C@@](C(C)C)(C(=O)N2CCc3cc(C(F)(F)F)cnc3C2)C1", "role": "Native_Control"},
]


def download_pdb(pdb_id: str, out_path: Path):
    if not out_path.exists():
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        print(f"[06] Downloading {pdb_id} from RCSB...")
        urllib.request.urlretrieve(url, out_path)


def parse_receptor_and_grid(gene: str, cfg: dict):
    pdb_id = cfg["pdb"]
    raw_pdb = REC_DIR / f"{pdb_id}_raw.pdb"
    clean_pdb = REC_DIR / f"{gene}_{pdb_id}_clean.pdb"
    rec_pdbqt = REC_DIR / f"{gene}_{pdb_id}.pdbqt"
    download_pdb(pdb_id, raw_pdb)

    atom_lines = []
    native_lines = []
    native_coords = []
    pocket_coords = []
    pocket_set = set(cfg.get("pocket_residues", []))

    # Identify first residue number for co-crystallized ligand so only one copy is extracted
    first_nat_resseq = None

    with open(raw_pdb) as fh:
        for line in fh:
            rec = line[:6].strip()
            chain = line[21].strip() if len(line) > 21 else ""
            resname = line[17:20].strip() if len(line) > 20 else ""
            resseq = line[22:26].strip() if len(line) > 26 else ""

            if rec == "ATOM" and chain == cfg["chain"]:
                altloc = line[16]
                if altloc in (" ", "A"):
                    atom_lines.append(line[:16] + " " + line[17:])
                    if f"{resname}{resseq}" in pocket_set:
                        x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                        pocket_coords.append([x, y, z])

            elif rec == "HETATM" and cfg["native_resname"] and resname == cfg["native_resname"]:
                if first_nat_resseq is None:
                    first_nat_resseq = (chain, resseq)
                if (chain, resseq) == first_nat_resseq:
                    native_lines.append(line)
                    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    native_coords.append([x, y, z])

    with open(clean_pdb, "w") as out:
        out.writelines(atom_lines)
        out.write("END\n")

    ref_lig_pdb = None
    if native_lines:
        ref_lig_pdb = REF_DIR / f"{gene}_{pdb_id}_{cfg['native_resname']}_crystal.pdb"
        with open(ref_lig_pdb, "w") as out:
            out.writelines(native_lines)
            out.write("END\n")

    if native_coords:
        center = np.mean(native_coords, axis=0)
        method = f"Co-crystallized ligand ({cfg['native_resname']}) centroid"
    elif pocket_coords:
        center = np.mean(pocket_coords, axis=0)
        method = f"Curated pocket residue centroid ({len(pocket_set)} residues)"
    else:
        all_xyz = [[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in atom_lines]
        center = np.mean(all_xyz, axis=0)
        method = "Chain centroid"

    cmd = ["obabel", str(clean_pdb), "-xr", "-O", str(rec_pdbqt), "--partialcharge", "gasteiger"]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    return {
        "gene": gene,
        "pdb": pdb_id,
        "chain": cfg["chain"],
        "native_resname": cfg["native_resname"],
        "crystal_ligand_pdb": str(ref_lig_pdb) if ref_lig_pdb else "",
        "description": cfg["description"],
        "grid_method": method,
        "center_x": round(float(center[0]), 3),
        "center_y": round(float(center[1]), 3),
        "center_z": round(float(center[2]), 3),
        "size_x": cfg["box_size"][0],
        "size_y": cfg["box_size"][1],
        "size_z": cfg["box_size"][2],
        "receptor_pdbqt": str(rec_pdbqt),
    }


def smiles_to_pdbqt(lig_id: str, smiles: str) -> Path:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES for {lig_id}: {smiles}")
    lfc = rdMolStandardize.LargestFragmentChooser()
    mol = lfc.choose(mol)
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    AllChem.EmbedMolecule(mol, params)
    AllChem.MMFFOptimizeMolecule(mol, maxIters=500)

    sdf_path = LIG_DIR / f"{lig_id}.sdf"
    pdbqt_path = LIG_DIR / f"{lig_id}.pdbqt"
    Chem.MolToMolFile(mol, str(sdf_path))

    cmd = ["obabel", str(sdf_path), "-O", str(pdbqt_path), "--partialcharge", "gasteiger"]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return pdbqt_path


def main():
    grids = {}
    for gene, cfg in TARGETS.items():
        info = parse_receptor_and_grid(gene, cfg)
        grids[gene] = info
        print(f"[06] Receptor {gene} ({cfg['pdb']}:{cfg['chain']}) -> center=({info['center_x']}, {info['center_y']}, {info['center_z']}) [{info['grid_method']}]")

    with open(DOCK_DIR / "grid_configs.json", "w") as fh:
        json.dump(grids, fh, indent=2)

    top10 = pd.read_csv(TAB_DIR / "top10_docking_candidates.tsv", sep="\t")
    manifest = []

    for _, row in top10.iterrows():
        raw_name = row["Phytochemical name_standardised"].strip()
        safe_id = "".join(c if c.isalnum() else "_" for c in raw_name).strip("_")
        pdbqt = smiles_to_pdbqt(safe_id, row["Canonical_SMILES"])
        manifest.append({
            "lig_id": safe_id,
            "name": raw_name,
            "imppat_id": row["IMPPAT_Phytochemical_identifier"],
            "pubchem_cid": row["Pubchem_CID"],
            "role": "Lead_Candidate",
            "pdbqt": str(pdbqt),
        })

    for item in EXTRA_LIGANDS:
        pdbqt = smiles_to_pdbqt(item["id"], item["smiles"])
        manifest.append({
            "lig_id": item["id"],
            "name": item["name"],
            "imppat_id": "REF",
            "pubchem_cid": "REF",
            "role": item["role"],
            "pdbqt": str(pdbqt),
        })

    pd.DataFrame(manifest).to_csv(DOCK_DIR / "ligand_manifest.tsv", sep="\t", index=False)
    print(f"[06] Prepared {len(manifest)} 3D MMFF94 ligands in {LIG_DIR}")


if __name__ == "__main__":
    main()
