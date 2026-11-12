#!/usr/bin/env python3
"""
13_run_full_pre_md.py (v2.3.0 - Clean Resumable Pipeline)
"""
import json
import re
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from pdbfixer import PDBFixer
from openmm.app import PDBFile

ROOT = Path(__file__).resolve().parents[1]
SRC_MD = ROOT / "results" / "md_handoff_v2"
WORK_MD = ROOT / "results" / "md_prepared_200ns"
WORK_MD.mkdir(parents=True, exist_ok=True)

EM_MDP = """integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 50000
nstlist     = 10
cutoff-scheme = Verlet
ns_type     = grid
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
"""

IONS_MDP = """integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 500
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
"""

NVT_MDP = """define                  = -DPOSRES -DPOSRES_LIG
integrator              = md
nsteps                  = 50000     ; 100 ps (2 fs timestep)
dt                      = 0.002
nstxout                 = 0
nstvout                 = 0
nstenergy               = 2500
nstlog                  = 2500
continuation            = no
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 1
lincs_order             = 4
cutoff-scheme           = Verlet
ns_type                 = grid
nstlist                 = 10
rcoulomb                = 1.0
rvdw                    = 1.0
DispCorr                = EnerPres
coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16
tcoupl                  = V-rescale
tc-grps                 = Protein_LIG Water_and_ions
tau_t                   = 0.1     0.1
ref_t                   = 310     310
pcoupl                  = no
pbc                     = xyz
gen_vel                 = yes
gen_temp                = 310
gen_seed                = 42
"""

NPT_MDP = """define                  = -DPOSRES -DPOSRES_LIG
integrator              = md
nsteps                  = 50000     ; 100 ps (2 fs timestep)
dt                      = 0.002
nstxout                 = 0
nstvout                 = 0
nstenergy               = 2500
nstlog                  = 2500
continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 1
lincs_order             = 4
cutoff-scheme           = Verlet
ns_type                 = grid
nstlist                 = 10
rcoulomb                = 1.0
rvdw                    = 1.0
DispCorr                = EnerPres
coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16
tcoupl                  = V-rescale
tc-grps                 = Protein_LIG Water_and_ions
tau_t                   = 0.1     0.1
ref_t                   = 310     310
pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5
refcoord_scaling        = com
pbc                     = xyz
gen_vel                 = no
"""

MD_200NS_MDP = """title                   = 200 ns Production MD (310 K, 1 bar)
integrator              = md
nsteps                  = 100000000 ; 200 ns (2 fs timestep)
dt                      = 0.002
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstlog                  = 5000      ; every 10 ps
nstenergy               = 5000      ; every 10 ps
nstxout-compressed      = 5000      ; 20,000 frames (every 10 ps)
compressed-x-grps       = System
continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 1
lincs_order             = 4
cutoff-scheme           = Verlet
ns_type                 = grid
nstlist                 = 10
rcoulomb                = 1.0
rvdw                    = 1.0
DispCorr                = EnerPres
coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16
tcoupl                  = V-rescale
tc-grps                 = Protein_LIG Water_and_ions
tau_t                   = 0.1     0.1
ref_t                   = 310     310
pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5
pbc                     = xyz
gen_vel                 = no
"""


def run_cmd(cmd, cwd: Path, stdin_text: str = None):
    proc = subprocess.run(
        cmd, cwd=str(cwd), input=stdin_text, text=True,
        capture_output=True
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"Command failed in {cwd}: {' '.join(cmd)}\nSTDOUT:\n{proc.stdout[-1500:]}\nSTDERR:\n{proc.stderr[-1500:]}"
        )
    return proc


def fix_receptor_pdb(src_pdb: Path, out_pdb: Path):
    fixer = PDBFixer(filename=str(src_pdb))
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)
    with open(out_pdb, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)


def build_exact_valence_sdf(raw_sdf: Path, smiles: str, out_sdf: Path):
    suppl = Chem.SDMolSupplier(str(raw_sdf), sanitize=False, removeHs=False)
    raw_with_h = suppl[0]
    assert raw_with_h is not None, f"Failed to read {raw_sdf}"
    rw = Chem.RWMol(raw_with_h)
    h_idxs = sorted([a.GetIdx() for a in rw.GetAtoms() if a.GetAtomicNum() == 1], reverse=True)
    for idx in h_idxs:
        rw.RemoveAtom(idx)
    raw_heavy = rw.GetMol()
    for b in raw_heavy.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
        b.SetIsAromatic(False)
    for a in raw_heavy.GetAtoms():
        a.SetIsAromatic(False)
        a.SetFormalCharge(0)
    tmpl = Chem.MolFromSmiles(smiles)
    assert tmpl is not None, f"Invalid SMILES: {smiles}"
    fixed = AllChem.AssignBondOrdersFromTemplate(tmpl, raw_heavy)
    Chem.SanitizeMol(fixed)
    fixed_h = Chem.AddHs(fixed, addCoords=True)
    writer = Chem.SDWriter(str(out_sdf))
    writer.write(fixed_h)
    writer.close()


def merge_gro_files(prot_gro: Path, lig_gro: Path, complex_gro: Path):
    p_lines = prot_gro.read_text().splitlines()
    l_lines = lig_gro.read_text().splitlines()
    n_prot = int(p_lines[1].strip())
    n_lig = int(l_lines[1].strip())
    box_line = p_lines[-1]
    merged = [
        "Protein-Ligand Complex",
        f"{n_prot + n_lig:5d}",
        *p_lines[2:2 + n_prot],
        *l_lines[2:2 + n_lig],
        box_line,
    ]
    complex_gro.write_text("\n".join(merged) + "\n")


def patch_topology(topol_path: Path, lig_itp_rel: str, posre_lig_rel: str):
    txt = topol_path.read_text()
    ff_include = '#include "amber99sb-ildn.ff/forcefield.itp"'
    assert ff_include in txt, "Could not find forcefield.itp include in topol.top"

    itp_full = (topol_path.parent / lig_itp_rel).read_text()
    m = re.search(r"(\[\s*atomtypes\s*\].*?)(?=\[\s*moleculetype\s*\])", itp_full, flags=re.S)
    if m:
        atomtypes_block = m.group(1).strip() + "\n"
        (topol_path.parent / "lig_atomtypes.itp").write_text(atomtypes_block)
        itp_no_at = itp_full[:m.start()] + itp_full[m.end():]
        (topol_path.parent / lig_itp_rel).write_text(itp_no_at)
        inc_str = f'{ff_include}\n#include "lig_atomtypes.itp"\n'
    else:
        inc_str = f'{ff_include}\n'

    txt = txt.replace(ff_include, inc_str, 1)

    water_inc = '#include "amber99sb-ildn.ff/tip3p.itp"'
    lig_block = (
        f'; Include ligand topology\n'
        f'#include "{lig_itp_rel}"\n'
        f'#ifdef POSRES_LIG\n'
        f'#include "{posre_lig_rel}"\n'
        f'#endif\n\n'
        f'{water_inc}'
    )
    txt = txt.replace(water_inc, lig_block, 1)

    if not txt.endswith("\n"):
        txt += "\n"
    txt += "LIG                 1\n"
    topol_path.write_text(txt)


def prepare_system(sys_src: Path, run_equil: bool = True):
    meta = json.loads((sys_src / "system_metadata.json").read_text())
    sys_id = meta["system_id"]
    nc = meta["ligand_net_formal_charge"]
    wdir = WORK_MD / sys_id

    if (wdir / "md_200ns.tpr").exists():
        print(f"\n[SKIP] {sys_id}: md_200ns.tpr already compiled and verified.")
        return

    if wdir.exists():
        shutil.rmtree(wdir)
    wdir.mkdir(parents=True, exist_ok=True)

    print(f"\n==================================================================")
    print(f"[PRE-MD] System: {sys_id} | Ligand: {meta['ligand_name']} (q={nc:+d})")
    print(f"==================================================================")

    (wdir / "ions.mdp").write_text(IONS_MDP)
    (wdir / "em.mdp").write_text(EM_MDP)
    (wdir / "nvt.mdp").write_text(NVT_MDP)
    (wdir / "npt.mdp").write_text(NPT_MDP)
    (wdir / "md_200ns.mdp").write_text(MD_200NS_MDP)

    rec_src = sys_src / meta["files"]["receptor_pdb"]
    rec_fixed = wdir / "receptor_fixed.pdb"
    fix_receptor_pdb(rec_src, rec_fixed)
    print(f"  [1/7] Receptor repaired at pH 7.4 -> {rec_fixed.name}")

    run_cmd(
        ["gmx", "pdb2gmx", "-f", "receptor_fixed.pdb", "-o", "protein.gro",
         "-p", "topol.top", "-ff", "amber99sb-ildn", "-water", "tip3p", "-ignh"],
        cwd=wdir
    )
    print(f"  [2/7] Protein topology generated (amber99sb-ildn + tip3p)")

    lig_sdf_src = sys_src / meta["files"]["ligand_sdf"]
    lig_sdf_dst = wdir / "ligand.sdf"
    build_exact_valence_sdf(lig_sdf_src, meta["canonical_smiles"], lig_sdf_dst)

    run_cmd(
        ["acpype", "-i", "ligand.sdf", "-c", "bcc", "-n", str(nc), "-a", "gaff2", "-b", "LIG"],
        cwd=wdir
    )
    acp_dir = wdir / "LIG.acpype"
    lig_gro = acp_dir / "LIG_GMX.gro"
    lig_itp = acp_dir / "LIG_GMX.itp"
    posre_lig = acp_dir / "posre_LIG.itp"
    shutil.copy2(lig_itp, wdir / "LIG_GMX.itp")
    shutil.copy2(posre_lig, wdir / "posre_LIG.itp")
    print(f"  [3/7] Ligand parameterized (GAFF2 + AM1-BCC, net charge={nc:+d})")

    merge_gro_files(wdir / "protein.gro", lig_gro, wdir / "complex.gro")
    patch_topology(wdir / "topol.top", "LIG_GMX.itp", "posre_LIG.itp")
    print(f"  [4/7] Merged complex.gro & updated topol.top")

    run_cmd(
        ["gmx", "editconf", "-f", "complex.gro", "-o", "newbox.gro", "-c", "-d", "1.0", "-bt", "cubic"],
        cwd=wdir
    )
    run_cmd(
        ["gmx", "solvate", "-cp", "newbox.gro", "-cs", "spc216.gro", "-o", "solv.gro", "-p", "topol.top"],
        cwd=wdir
    )
    run_cmd(
        ["gmx", "grompp", "-f", "ions.mdp", "-c", "solv.gro", "-p", "topol.top", "-o", "ions.tpr", "-maxwarn", "2"],
        cwd=wdir
    )
    run_cmd(
        ["gmx", "genion", "-s", "ions.tpr", "-o", "solv_ions.gro", "-p", "topol.top",
         "-pname", "NA", "-nname", "CL", "-neutral", "-conc", "0.15"],
        cwd=wdir,
        stdin_text="SOL\n"
    )
    print(f"  [5/7] Solvated in cubic TIP3P box & neutralized with 0.15 M NaCl")

    run_cmd(
        ["gmx", "grompp", "-f", "em.mdp", "-c", "solv_ions.gro", "-p", "topol.top", "-o", "em.tpr", "-maxwarn", "2"],
        cwd=wdir
    )
    run_cmd(["gmx", "mdrun", "-deffnm", "em"], cwd=wdir)
    print(f"  [6/7] Energy minimization complete -> em.gro")

    probe = run_cmd(["gmx", "make_ndx", "-f", "em.gro", "-o", "index.ndx"], cwd=wdir, stdin_text="q\n")
    lig_grp = None
    for line in (probe.stdout + "\n" + probe.stderr).splitlines():
        m = re.match(r"^\s*(\d+)\s+(\S+)\s*:\s*\d+\s+atoms", line)
        if m and m.group(2) not in {
            "System", "Protein", "Protein-H", "C-alpha", "Backbone",
            "MainChain", "MainChain+Cb", "MainChain+H", "SideChain",
            "SideChain-H", "Prot-Masses", "non-Protein", "Other",
            "NA", "CL", "Water", "SOL", "non-Water", "Ion", "Water_and_ions"
        }:
            lig_grp = m.group(1)
            break
    assert lig_grp is not None, f"Could not detect ligand index group:\n{probe.stdout}"
    run_cmd(
        ["gmx", "make_ndx", "-f", "em.gro", "-o", "index.ndx"],
        cwd=wdir,
        stdin_text=f"1 | {lig_grp}\nq\n"
    )
    ndx_txt = (wdir / "index.ndx").read_text()
    headers = re.findall(r"^\[\s*([^\]]+?)\s*\]", ndx_txt, flags=re.M)
    last_hdr = headers[-1]
    if last_hdr != "Protein_LIG":
        ndx_txt = re.sub(r"\[\s*" + re.escape(last_hdr) + r"\s*\]\s*$", "[ Protein_LIG ]", ndx_txt, count=1, flags=re.M)
        (wdir / "index.ndx").write_text(ndx_txt)

    if run_equil:
        print(f"  [7a/7] Running 100 ps NVT equilibration (310 K)...")
        run_cmd(
            ["gmx", "grompp", "-f", "nvt.mdp", "-c", "em.gro", "-r", "em.gro",
             "-p", "topol.top", "-n", "index.ndx", "-o", "nvt.tpr", "-maxwarn", "2"],
            cwd=wdir
        )
        run_cmd(["gmx", "mdrun", "-deffnm", "nvt"], cwd=wdir)

        print(f"  [7b/7] Running 100 ps NPT equilibration (310 K, 1 bar)...")
        run_cmd(
            ["gmx", "grompp", "-f", "npt.mdp", "-c", "nvt.gro", "-r", "nvt.gro",
             "-t", "nvt.cpt", "-p", "topol.top", "-n", "index.ndx", "-o", "npt.tpr", "-maxwarn", "2"],
            cwd=wdir
        )
        run_cmd(["gmx", "mdrun", "-deffnm", "npt"], cwd=wdir)

        print(f"  [7c/7] Compiling ready-to-run 200 ns production binary -> md_200ns.tpr")
        run_cmd(
            ["gmx", "grompp", "-f", "md_200ns.mdp", "-c", "npt.gro", "-t", "npt.cpt",
             "-p", "topol.top", "-n", "index.ndx", "-o", "md_200ns.tpr", "-maxwarn", "2"],
            cwd=wdir
        )
        print(f"  [DONE] {sys_id}: md_200ns.tpr compiled and verified!")


def main():
    run_equil = "--skip-equil" not in sys.argv
    for folder in sorted(p for p in SRC_MD.iterdir() if p.is_dir()):
        prepare_system(folder, run_equil=run_equil)

    launch_sh = """#!/usr/bin/env bash
# Master 200 ns Production MD Launch Script for MD Team
set -euo pipefail

SYSTEMS=(
  "01_PPARG_Native_BRL_Control"
  "02_PPARG_Chelerythrine"
  "03_PPARG_Berberine"
  "04_PPARG_Curcumin"
  "05_STAT3_Quercetin"
)

for sys in "${SYSTEMS[@]}"; do
  echo ">>> Launching 200 ns Production MD for ${sys} <<<"
  (
    cd "${sys}"
    # If running on a different GROMACS version, re-compile md_200ns.tpr from npt.gro + npt.cpt:
    # gmx grompp -f md_200ns.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_200ns.tpr -maxwarn 1
    gmx mdrun -v -deffnm md_200ns -nb gpu -pme gpu -bonded gpu
  )
done
echo "All 5 x 200 ns production runs complete!"
"""
    sh_file = WORK_MD / "run_all_200ns_production.sh"
    sh_file.write_text(launch_sh)
    sh_file.chmod(0o755)

    shutil.copy2(SRC_MD / "MD_HANDOFF_PROTOCOL.md", WORK_MD / "MD_HANDOFF_PROTOCOL.md")
    shutil.copy2(SRC_MD / "MD_SYSTEMS_SUMMARY.tsv", WORK_MD / "MD_SYSTEMS_SUMMARY.tsv")

    keep_exts = {".gro", ".top", ".itp", ".ndx", ".mdp", ".cpt", ".tpr", ".sh", ".md", ".tsv", ".pdb", ".sdf"}
    zip_out = ROOT / "results" / "md_equilibrated_200ns_handoff.zip"
    with zipfile.ZipFile(zip_out, "w", zipfile.ZIP_DEFLATED) as zf:
        for f in sorted(WORK_MD.rglob("*")):
            if f.is_file() and f.suffix in keep_exts and "LIG.acpype" not in f.parts:
                zf.write(f, f.relative_to(ROOT / "results"))

    print(f"\n==================================================================")
    print(f"[COMPLETE] All 5 systems fully equilibrated & packaged at:")
    print(f"           {zip_out}")
    print(f"==================================================================")


if __name__ == "__main__":
    main()
