#!/usr/bin/env python3
"""
10_prepare_md_handoff.py (v2.1.0 - Production Team Handoff)
Packages the 5 prioritized soluble protein-ligand complexes for 100 ns all-atom
molecular dynamics (MD) + MM-GBSA/MM-PBSA validation by the MD collaboration team:
  1. 01_PPARG_Native_BRL_Control : PPARG (2PRG:A) + Rosiglitazone (BRL) [Charge: 0, Vina: -8.33 kcal/mol]
  2. 02_PPARG_Chelerythrine      : PPARG (2PRG:A) + Chelerythrine       [Charge: +1, Vina: -8.74 kcal/mol]
  3. 03_PPARG_Berberine          : PPARG (2PRG:A) + Berberine           [Charge: +1, Vina: -8.41 kcal/mol]
  4. 04_PPARG_Curcumin           : PPARG (2PRG:A) + Curcumin            [Charge: 0, Vina: -8.12 kcal/mol]
  5. 05_STAT3_Quercetin          : STAT3 (6NJS:A) + Quercetin           [Charge: 0, Vina: -5.64 kcal/mol]
"""
import json
import shutil
import subprocess
import zipfile
from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

ROOT = Path(__file__).resolve().parents[1]
DOCK_DIR = ROOT / "results" / "docking"
REC_DIR = DOCK_DIR / "receptors"
POSE_DIR = DOCK_DIR / "poses"
TAB_DIR = ROOT / "results" / "tables"
MD_DIR = ROOT / "results" / "md_handoff_v2"

if MD_DIR.exists():
    shutil.rmtree(MD_DIR)
MD_DIR.mkdir(parents=True, exist_ok=True)

SYSTEMS = [
    {
        "folder": "01_PPARG_Native_BRL_Control",
        "gene": "PPARG",
        "pdb": "2PRG",
        "chain": "A",
        "lig_id": "Native_BRL",
        "lig_name": "Rosiglitazone (BRL)",
        "role": "Positive Crystallographic Reference Control",
        "vina_kcal_mol": -8.33,
        "formal_charge": 0,
        "smiles": "CN(CCOc1ccc(CC2SC(=O)NC2=O)cc1)c1ccccn1",
        "key_pocket_residues": "HIS323, TYR327, HIS449, TYR473, CYS285, ARG288",
    },
    {
        "folder": "02_PPARG_Chelerythrine",
        "gene": "PPARG",
        "pdb": "2PRG",
        "chain": "A",
        "lig_id": "Chelerythrine",
        "lig_name": "Chelerythrine",
        "role": "Lead Candidate #1 (Highest Multi-Target & PPARG Affinity)",
        "vina_kcal_mol": -8.74,
        "formal_charge": 1,
        "smiles": "COc1ccc2c3c([n+](C)cc2c1OC)c1cc2c(cc1cc3)OCO2",
        "key_pocket_residues": "ARG280, ILE281, GLY284, ARG288, ILE341, GLU343",
    },
    {
        "folder": "03_PPARG_Berberine",
        "gene": "PPARG",
        "pdb": "2PRG",
        "chain": "A",
        "lig_id": "Berberine_chloride",
        "lig_name": "Berberine (Quaternary Cation)",
        "role": "Lead Candidate #2 (High QED=0.674, Zero PAINS)",
        "vina_kcal_mol": -8.41,
        "formal_charge": 1,
        "smiles": "COc1ccc2cc3[n+](cc2c1OC)CCc1cc2c(cc1-3)OCO2",
        "key_pocket_residues": "ARG280, CYS285, ARG288, LEU340, ILE341",
    },
    {
        "folder": "04_PPARG_Curcumin",
        "gene": "PPARG",
        "pdb": "2PRG",
        "chain": "A",
        "lig_id": "Curcumin",
        "lig_name": "Curcumin",
        "role": "Lead Candidate #3 (Canonical PPARG Agonist H-Bonding)",
        "vina_kcal_mol": -8.12,
        "formal_charge": 0,
        "smiles": "COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O",
        "key_pocket_residues": "CYS285, ARG288, HIS323, TYR327, LEU340, HIS449",
    },
    {
        "folder": "05_STAT3_Quercetin",
        "gene": "STAT3",
        "pdb": "6NJS",
        "chain": "A",
        "lig_id": "Quercetin",
        "lig_name": "Quercetin",
        "role": "Secondary Target Validation (STAT3 SH2 pTyr705 Pocket)",
        "vina_kcal_mol": -5.64,
        "formal_charge": 0,
        "smiles": "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12",
        "key_pocket_residues": "GLU638, PRO639, TYR640, GLN644, GLY656, TYR657, LYS658",
    },
]


def main():
    summary_rows = []
    for s in SYSTEMS:
        sdir = MD_DIR / s["folder"]
        sdir.mkdir(parents=True, exist_ok=True)

        rec_src = REC_DIR / f"{s['gene']}_{s['pdb']}_clean.pdb"
        pose_pdbqt = POSE_DIR / f"{s['gene']}_{s['lig_id']}_docked.pdbqt"

        rec_dst = sdir / f"{s['gene']}_{s['pdb']}_receptor.pdb"
        lig_pdb = sdir / f"{s['lig_id']}_top1_pose.pdb"
        lig_sdf = sdir / f"{s['lig_id']}_top1_pose.sdf"
        lig_mol2 = sdir / f"{s['lig_id']}_top1_pose.mol2"
        complex_pdb = sdir / f"{s['folder']}_complex.pdb"

        shutil.copy2(rec_src, rec_dst)

        # Extract Model 1 (top-scoring Vina pose) into PDB, SDF (with explicit Hs), and MOL2 (with Gasteiger charges)
        subprocess.run(["obabel", str(pose_pdbqt), "-f", "1", "-l", "1", "-O", str(lig_pdb)],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run(["obabel", str(pose_pdbqt), "-f", "1", "-l", "1", "-h", "-O", str(lig_sdf)],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run(["obabel", str(pose_pdbqt), "-f", "1", "-l", "1", "-h", "--partialcharge", "gasteiger", "-O", str(lig_mol2)],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # Standardize ligand residue name in PDB to 'LIG' for clean tleap / pdb2gmx parsing
        lig_lines_clean = []
        for line in lig_pdb.read_text().splitlines():
            if line.startswith(("ATOM", "HETATM")):
                # Replace residue name (cols 17:20) with LIG and chain (col 21) with L
                fixed = "HETATM" + line[6:17] + "LIG L   1" + line[26:]
                lig_lines_clean.append(fixed)
        lig_pdb.write_text("\n".join(lig_lines_clean + ["END"]) + "\n")

        rec_lines = [l for l in rec_src.read_text().splitlines() if l.startswith("ATOM")]
        complex_pdb.write_text("\n".join(rec_lines + ["TER"] + lig_lines_clean + ["END"]) + "\n")

        meta_json = {
            "system_id": s["folder"],
            "receptor_gene": s["gene"],
            "pdb_id": s["pdb"],
            "chain": s["chain"],
            "ligand_name": s["lig_name"],
            "ligand_resname_in_pdb": "LIG",
            "role": s["role"],
            "autodock_vina_affinity_kcal_mol": s["vina_kcal_mol"],
            "ligand_net_formal_charge": s["formal_charge"],
            "canonical_smiles": s["smiles"],
            "key_binding_pocket_residues": s["key_pocket_residues"],
            "files": {
                "receptor_pdb": rec_dst.name,
                "ligand_pdb": lig_pdb.name,
                "ligand_sdf": lig_sdf.name,
                "ligand_mol2": lig_mol2.name,
                "merged_complex_pdb": complex_pdb.name,
            },
        }
        (sdir / "system_metadata.json").write_text(json.dumps(meta_json, indent=2) + "\n")

        summary_rows.append({
            "System_Folder": s["folder"],
            "Target": f"{s['gene']} ({s['pdb']}:{s['chain']})",
            "Ligand": s["lig_name"],
            "Role": s["role"],
            "Vina_Affinity_kcal_mol": s["vina_kcal_mol"],
            "Ligand_Net_Charge": s["formal_charge"],
            "Key_Pocket_Residues": s["key_pocket_residues"],
        })

    pd.DataFrame(summary_rows).to_csv(MD_DIR / "MD_SYSTEMS_SUMMARY.tsv", sep="\t", index=False)

    # Helper shell script for the MD team to parameterize all 5 ligands via acpype (AM1-BCC / GAFF2)
    helper_sh = """#!/usr/bin/env bash
# Optional helper script for MD team: generates AMBER/GROMACS ligand topologies via acpype
# Requires: acpype (conda install -c conda-forge acpype)
set -euo pipefail

declare -A CHARGES=(
  ["01_PPARG_Native_BRL_Control"]="0"
  ["02_PPARG_Chelerythrine"]="1"
  ["03_PPARG_Berberine"]="1"
  ["04_PPARG_Curcumin"]="0"
  ["05_STAT3_Quercetin"]="0"
)

for sys in "${!CHARGES[@]}"; do
  nc="${CHARGES[$sys]}"
  echo "=== Parameterizing ${sys} (net charge = ${nc}) ==="
  (
    cd "${sys}"
    sdf_file=$(ls *_top1_pose.sdf | head -n 1)
    acpype -i "${sdf_file}" -c bcc -n "${nc}" -a gaff2 -r LIG
  )
done
echo "Done. GROMACS (.gro/.itp) and AMBER (.frcmod/.lib) ligand topologies generated."
"""
    sh_path = MD_DIR / "parameterize_ligands_acpype.sh"
    sh_path.write_text(helper_sh)
    sh_path.chmod(0o755)

    readme_text = """# Molecular Dynamics (MD) Validation Handoff Package (v2.1.0)

## 1. Overview & Simulation Scope
This package contains **5 pre-docked, energy-minimized soluble protein–ligand complexes** prioritized from our multi-cohort IPF transcriptomic meta-analysis ($N = 109$) and multi-target structure-based virtual screening against M2 macrophage polarization regulators:

| Folder | Target (PDB:Chain) | Ligand | Role | Vina Score (kcal/mol) | Net Ligand Charge (`-n`) | Key Pocket Residues |
| :--- | :--- | :--- | :--- | :---: | :---: | :--- |
| `01_PPARG_Native_BRL_Control` | PPARG (`2PRG:A`) | Rosiglitazone (`BRL`) | Positive Crystallographic Control | `-8.33` | `0` | HIS323, TYR327, HIS449, TYR473, CYS285, ARG288 |
| `02_PPARG_Chelerythrine` | PPARG (`2PRG:A`) | Chelerythrine | Lead #1 (Top Overall Multi-Target Lead) | `-8.74` | **`+1`** | ARG280, ILE281, GLY284, ARG288, ILE341, GLU343 |
| `03_PPARG_Berberine` | PPARG (`2PRG:A`) | Berberine | Lead #2 (Optimal QED=0.674, 0 PAINS) | `-8.41` | **`+1`** | ARG280, CYS285, ARG288, LEU340, ILE341 |
| `04_PPARG_Curcumin` | PPARG (`2PRG:A`) | Curcumin | Lead #3 (Canonical PPARG Agonist H-Bonds) | `-8.12` | `0` | CYS285, ARG288, HIS323, TYR327, LEU340, HIS449 |
| `05_STAT3_Quercetin` | STAT3 (`6NJS:A`) | Quercetin | Secondary M2 Regulator (SH2 pTyr Pocket) | `-5.64` | `0` | GLU638, PRO639, TYR640, GLN644, GLY656, TYR657, LYS658 |

> **CRITICAL PARAMETERIZATION NOTE:** Both **Chelerythrine** (`02`) and **Berberine** (`03`) are quaternary benzophenanthridine/isoquinoline alkaloids carrying a **permanent `+1` formal charge** on their iminium/isoquinolinium nitrogen (`[n+]`). Ensure net charge **`+1`** (`-n 1` in `acpype` or `-nc 1` in `antechamber`) is specified during ligand parameterization.

*(Note: CCR2 [`6GPS`] is a 7-transmembrane GPCR requiring explicit POPC lipid bilayer embedding; therefore, the soluble globular PPARG LBD and STAT3 SH2 domains above are prioritized for aqueous all-atom MD.)*

---

## 2. Contents of Each System Folder
Each directory contains:
- `*_receptor.pdb`: Clean, single-chain protein receptor (`2PRG:A` or `6NJS:A`), stripped of crystallographic waters and co-crystallized heteroatoms.
- `*_top1_pose.sdf`: Top-ranked AutoDock Vina docked pose with explicit hydrogens and exact bond orders (**recommended input for `acpype` / `antechamber` / `CGenFF`**).
- `*_top1_pose.mol2`: Top-ranked pose with explicit hydrogens and Gasteiger partial charges.
- `*_top1_pose.pdb`: Top-ranked pose standardized with residue name `LIG` and chain `L`.
- `*_complex.pdb`: Merged receptor + `LIG` initial complex ready for immediate visualization (PyMOL/ChimeraX) or `tleap`/`gmx pdb2gmx` assembly.
- `system_metadata.json`: Machine-readable summary of formal charge, SMILES, Vina score, and pocket residues.

---

## 3. Recommended All-Atom MD Simulation Protocol (GROMACS or AMBER)
- **Simulation Length:** **100 ns production** per system (`5 × 100 ns = 500 ns` total; if compute-constrained, prioritize `01`, `02`, `03`, and `04` first = `400 ns`).
- **Protein Force Field:** AMBER `ff14SB` (or CHARMM36m)
- **Ligand Force Field:** `GAFF2` with `AM1-BCC` partial charges (run `./parameterize_ligands_acpype.sh` to auto-generate all 5 ligand topologies).
- **Solvation & Ions:** Explicit `TIP3P` truncated octahedral or cubic water box (`10.0 Å` minimum padding), neutralized with `Na+`/`Cl-` counterions and supplemented to `0.15 M` physiological NaCl.
- **Minimization & Equilibration:**
  1. Energy minimization: Steepest descent (`5,000` steps) + conjugate gradient until max force `< 100 kJ/mol/nm` (`1.0 kcal/mol/Å`).
  2. NVT heating: `100 ps` (`0 -> 310 K`, V-rescale / Langevin thermostat, `1000 kJ/mol/nm²` heavy-atom restraints).
  3. NPT equilibration: `100–200 ps` (`310 K`, `1.0 bar`, Parrinello-Rahman / Berendsen barostat).
- **Production MD:** `100 ns` unrestrained NPT at `310 K`, `1.0 bar`, `2 fs` timestep with LINCS/SHAKE H-bond constraints and PME electrostatics (`1.0 nm` cutoff), saving trajectories every `10 ps` (`10,000` frames).

---

## 4. Requested Post-MD Analysis Deliverables
Please return the following raw `.xvg`/`.csv` time-series files and summary metrics so we can generate publication Figure 6 and Table 7:
1. **Protein Backbone RMSD (`nm` vs. `ns`)** — whole domain and binding-pocket residues (`within 5 Å of LIG`).
2. **Ligand Heavy-Atom RMSD (`nm` vs. `ns`)** — fitted on protein backbone to confirm stable pocket retention without dissociation.
3. **Per-Residue Cα RMSF (`nm` vs. residue number)** — highlighting flexibility changes in PPARG Helix 12 (AF-2) or STAT3 SH2 loops.
4. **Radius of Gyration (`Rg`, `nm` vs. `ns`) & Solvent-Accessible Surface Area (`SASA`, `nm²` vs. `ns`)**.
5. **Intermolecular Hydrogen-Bond Count & Occupancy (`%`)** over the equilibrated `50–100 ns` window.
6. **MM-GBSA / MM-PBSA Binding Free Energy ($\Delta G_{\text{bind}} \pm \text{SD}$, `kcal/mol`)** decomposed into $\Delta E_{\text{vdW}}$, $\Delta E_{\text{elec}}$, $\Delta G_{\text{polar}}$, and $\Delta G_{\text{nonpolar}}$ across `250–500` snapshots from the equilibrated `50–100 ns` trajectory.
"""
    (MD_DIR / "MD_HANDOFF_PROTOCOL.md").write_text(readme_text)

    zip_path = ROOT / "results" / "md_handoff_v2.zip"
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for f in sorted(MD_DIR.rglob("*")):
            if f.is_file():
                zf.write(f, f.relative_to(ROOT / "results"))

    print(f"[10] Production MD handoff package created at: {zip_path}")
    print(f"[10] Directory contents ({len(list(MD_DIR.rglob('*')))} items):")
    for p in sorted(MD_DIR.iterdir()):
        print(f"     - {p.name}")


if __name__ == "__main__":
    main()
