#!/usr/bin/env bash
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
