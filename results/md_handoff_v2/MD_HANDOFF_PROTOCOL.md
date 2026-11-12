# Molecular Dynamics (MD) Validation Handoff Package (v2.1.0)

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
- **Simulation Length:** **200 ns production** per system (`5 × 200 ns = 1,000 ns (1.0 us)` total; if compute-constrained, prioritize `01`, `02`, `03`, and `04` first = `800 ns`).
- **Protein Force Field:** AMBER `ff14SB` (or CHARMM36m)
- **Ligand Force Field:** `GAFF2` with `AM1-BCC` partial charges (run `./parameterize_ligands_acpype.sh` to auto-generate all 5 ligand topologies).
- **Solvation & Ions:** Explicit `TIP3P` truncated octahedral or cubic water box (`10.0 Å` minimum padding), neutralized with `Na+`/`Cl-` counterions and supplemented to `0.15 M` physiological NaCl.
- **Minimization & Equilibration:**
  1. Energy minimization: Steepest descent (`5,000` steps) + conjugate gradient until max force `< 100 kJ/mol/nm` (`1.0 kcal/mol/Å`).
  2. NVT heating: `100 ps` (`0 -> 310 K`, V-rescale / Langevin thermostat, `1000 kJ/mol/nm²` heavy-atom restraints).
  3. NPT equilibration: `100–200 ps` (`310 K`, `1.0 bar`, Parrinello-Rahman / Berendsen barostat).
- **Production MD:** `200 ns` unrestrained NPT at `310 K`, `1.0 bar`, `2 fs` timestep with LINCS/SHAKE H-bond constraints and PME electrostatics (`1.0 nm` cutoff), saving trajectories every `10 ps` (`20,000` frames).

---

## 4. Requested Post-MD Analysis Deliverables
Please return the following raw `.xvg`/`.csv` time-series files and summary metrics so we can generate publication Figure 6 and Table 7:
1. **Protein Backbone RMSD (`nm` vs. `ns`)** — whole domain and binding-pocket residues (`within 5 Å of LIG`).
2. **Ligand Heavy-Atom RMSD (`nm` vs. `ns`)** — fitted on protein backbone to confirm stable pocket retention without dissociation.
3. **Per-Residue Cα RMSF (`nm` vs. residue number)** — highlighting flexibility changes in PPARG Helix 12 (AF-2) or STAT3 SH2 loops.
4. **Radius of Gyration (`Rg`, `nm` vs. `ns`) & Solvent-Accessible Surface Area (`SASA`, `nm²` vs. `ns`)**.
5. **Intermolecular Hydrogen-Bond Count & Occupancy (`%`)** over the equilibrated `50–200 ns` window.
6. **MM-GBSA / MM-PBSA Binding Free Energy ($\Delta G_{	ext{bind}} \pm 	ext{SD}$, `kcal/mol`)** decomposed into $\Delta E_{	ext{vdW}}$, $\Delta E_{	ext{elec}}$, $\Delta G_{	ext{polar}}$, and $\Delta G_{	ext{nonpolar}}$ across `250–500` snapshots from the equilibrated `50–200 ns` trajectory.
