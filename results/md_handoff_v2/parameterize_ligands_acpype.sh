#!/usr/bin/env bash
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
