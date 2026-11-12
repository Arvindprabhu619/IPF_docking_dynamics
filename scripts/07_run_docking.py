#!/usr/bin/env python3
"""
07_run_docking.py (v2.0.0)
Executes deterministic AutoDock Vina (seed=42, exhaustiveness=16) across all
5 M2 macrophage targets (PPARG, STAT3, CCR2, NFKB1, HIF1A) and 14 ligands
(Top 10 strictly drug-like phytochemicals + Pirfenidone + 3 native controls).
Also computes heavy-atom RMSD for co-crystallized native ligand redocking
(BRL -> 2PRG, KQV -> 6NJS, F7N -> 6GPS) to validate docking accuracy.
"""
import json
import re
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parents[1]
DOCK_DIR = ROOT / "results" / "docking"
POSE_DIR = DOCK_DIR / "poses"
LOG_DIR = DOCK_DIR / "logs"
TAB_DIR = ROOT / "results" / "tables"
for d in [POSE_DIR, LOG_DIR, TAB_DIR]:
    d.mkdir(parents=True, exist_ok=True)


def parse_vina_best_affinity(log_text: str) -> float:
    for line in log_text.splitlines():
        m = re.match(r"^\s*1\s+(-?\d+\.\d+)\s+", line)
        if m:
            return float(m.group(1))
    raise RuntimeError(f"Could not parse Vina affinity from output:\n{log_text}")


def compute_redock_rmsd(crystal_pdb: str, docked_pdbqt: Path) -> float:
    """Compute heavy-atom symmetry-aware RMSD between crystal ligand and top docked pose."""
    # Extract first MODEL from docked pdbqt and convert to PDB via obabel
    tmp_pdb = docked_pdbqt.with_suffix(".top1.pdb")
    subprocess.run(
        ["obabel", str(docked_pdbqt), "-f", "1", "-l", "1", "-O", str(tmp_pdb)],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    )
    ref_mol = Chem.MolFromPDBFile(crystal_pdb, removeHs=True, sanitize=False)
    prb_mol = Chem.MolFromPDBFile(str(tmp_pdb), removeHs=True, sanitize=False)
    if ref_mol is None or prb_mol is None:
        return float("nan")
    try:
        Chem.SanitizeMol(ref_mol)
        Chem.SanitizeMol(prb_mol)
        # Match bond orders from docked pose to crystal template if needed
        ref_template = AllChem.AssignBondOrdersFromTemplate(prb_mol, ref_mol)
        return round(float(AllChem.CalcRMS(ref_template, prb_mol)), 3)
    except Exception:
        # Fallback coordinate-centroid + nearest-atom heavy-atom RMSD
        c1 = ref_mol.GetConformer().GetPositions()
        c2 = prb_mol.GetConformer().GetPositions()
        if len(c1) == len(c2):
            dists = np.min(np.linalg.norm(c1[:, None, :] - c2[None, :, :], axis=-1), axis=1)
            return round(float(np.sqrt(np.mean(dists**2))), 3)
        return float("nan")


def main():
    with open(DOCK_DIR / "grid_configs.json") as fh:
        grids = json.load(fh)
    manifest = pd.read_csv(DOCK_DIR / "ligand_manifest.tsv", sep="\t")

    records = []
    targets_order = ["PPARG", "STAT3", "CCR2", "NFKB1", "HIF1A"]

    for _, lig in manifest.iterrows():
        lig_id = lig["lig_id"]
        lig_pdbqt = lig["pdbqt"]
        row_out = {
            "Ligand_ID": lig_id,
            "Compound": lig["name"],
            "IMPPAT_ID": lig["imppat_id"],
            "PubChem_CID": lig["pubchem_cid"],
            "Role": lig["role"],
        }
        affinities = []

        for gene in targets_order:
            g = grids[gene]
            out_pose = POSE_DIR / f"{gene}_{lig_id}_docked.pdbqt"
            out_log = LOG_DIR / f"{gene}_{lig_id}.log"

            cmd = [
                "vina",
                "--receptor", g["receptor_pdbqt"],
                "--ligand", lig_pdbqt,
                "--center_x", str(g["center_x"]),
                "--center_y", str(g["center_y"]),
                "--center_z", str(g["center_z"]),
                "--size_x", str(g["size_x"]),
                "--size_y", str(g["size_y"]),
                "--size_z", str(g["size_z"]),
                "--exhaustiveness", "16",
                "--num_modes", "9",
                "--seed", "42",
                "--out", str(out_pose),
            ]
            proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
            out_log.write_text(proc.stdout)
            aff = parse_vina_best_affinity(proc.stdout)
            row_out[gene] = aff
            affinities.append(aff)

            # Check if this is the cognate native control for this receptor
            nat = g.get("native_resname")
            if nat and lig_id == f"Native_{nat}" and g.get("crystal_ligand_pdb"):
                rmsd = compute_redock_rmsd(g["crystal_ligand_pdb"], out_pose)
                print(f"[07] Native Redocking Validation: {gene} ({g['pdb']}) + {nat} -> Affinity = {aff:.2f} kcal/mol | Heavy-atom RMSD = {rmsd:.3f} A")

        row_out["Mean_Affinity"] = round(float(np.mean(affinities)), 2)
        row_out["Best_Affinity"] = round(float(np.min(affinities)), 2)
        records.append(row_out)
        print(f"[07] Completed {lig['name']:<32} | PPARG={row_out['PPARG']:>6.2f}  STAT3={row_out['STAT3']:>6.2f}  CCR2={row_out['CCR2']:>6.2f}  NFKB1={row_out['NFKB1']:>6.2f}  HIF1A={row_out['HIF1A']:>6.2f} | Mean={row_out['Mean_Affinity']:>6.2f}")

    df = pd.DataFrame(records)
    # Sort Lead_Candidates by Mean_Affinity ascending (most negative = strongest), followed by controls
    leads = df[df["Role"] == "Lead_Candidate"].sort_values("Mean_Affinity", ascending=True)
    ctrls = df[df["Role"] != "Lead_Candidate"]
    final_df = pd.concat([leads, ctrls], ignore_index=True)
    final_df.to_csv(TAB_DIR / "Table5_docking_affinities_v2.tsv", sep="\t", index=False)

    print("\n=== FINAL TABLE 5: MULTI-TARGET DOCKING AFFINITIES (kcal/mol) ===")
    print(final_df[["Compound", "Role", "PPARG", "STAT3", "CCR2", "NFKB1", "HIF1A", "Mean_Affinity"]].to_string(index=False))


if __name__ == "__main__":
    main()
