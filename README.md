# Multi-Cohort Transcriptomic Meta-Analysis, Macrophage Deconvolution, and Phytochemical Docking/MD Pipeline for Idiopathic Pulmonary Fibrosis (IPF)

## Reproducible Pipeline Order
1. **Differential Expression (`scripts/02_tximport_deseq2.R`)**: Imports Salmon `quant.sf` files via `tximport` (`data/reference/tx2gene.tsv`) and runs batch-corrected DESeq2 meta-analysis.
2. **Cohort QC & DEG Plots (`scripts/03_cohort_qc_and_deg_plots.R`)**: Generates PCA/UMAP batch-correction diagnostics, Volcano, MA, and top-100 DEG heatmaps (`Fig1`, `Fig2`, `SuppFigS1`).
3. **Macrophage Deconvolution (`scripts/04_macrophage_deconvolution.R`)**: Runs quanTIseq deconvolution and M1/M2/IPF-macrophage signature scoring (`Fig3`, `SuppFigS2`).
4. **IMPPAT Phytochemical Screen (`scripts/05_imppat_screen.py`)**: Filters IMPPAT phytochemicals for drug-likeness and M2-regulator activity.
5. **Molecular Docking (`scripts/06_prepare_docking.py`, `scripts/07_run_docking.py`)**: Performs AutoDock Vina docking and crystal ligand redocking validation across 5 IPF targets (`PPARG`, `STAT3`, `CCR2`, `HIF1A`, `NFKB1`).
6. **ADMET & Interaction Mapping (`scripts/08_admet_and_interactions.py`, `scripts/09_generate_docking_figures.R`, `scripts/12_fix_fig4_fig5.R`)**: Computes ADMET profiles and binding-pocket interactions (`Fig4`, `Fig5`).
7. **Automated Pre-MD Equilibration (`scripts/10_prepare_md_handoff.py`, `scripts/13_run_full_pre_md.py`)**: Repairs receptors at pH 7.4 (`PDBFixer`), assigns exact formal charges/bond orders via RDKit template mapping, parameterizes ligands with `acpype` (`GAFF2` + `AM1-BCC`), solvates in TIP3P + `0.15 M` NaCl (`amber99sb-ildn`), energy-minimizes, equilibrates (`100 ps` NVT + `100 ps` NPT at `310 K, 1 bar`), and compiles `md_200ns.tpr`.

## Regenerating Excluded Heavy MD Files
Heavy GROMACS binaries (`*.tpr`, `*.gro`, `*.cpt`, `*.zip`) are excluded from Git tracking and can be regenerated cleanly from `results/md_handoff_v2/`:

    conda activate md_prep
    python3 scripts/13_run_full_pre_md.py
