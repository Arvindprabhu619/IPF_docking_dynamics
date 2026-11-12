#!/usr/bin/env Rscript
# =============================================================================
#  04_macrophage_deconvolution.R (v2.0.0)
#  1. Macrophage polarization scoring (M1, M2, IPF_Mac with MERTK) on ComBat & limma
#  2. xCell cell-type deconvolution (M1 vs M2 macrophage enrichment) on Salmon TPM
#  3. Generates Fig 3A-D and Supp Fig S2A-C
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggpubr)
  library(rtracklayer)
  library(limma)
  library(reshape2)
})

set.seed(42)

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
root <- if (length(script_path) > 0) {
  normalizePath(file.path(dirname(script_path), ".."))
} else {
  getwd()
}

tab_dir <- file.path(root, "results/tables")
dec_dir <- file.path(root, "results/deconvolution")
fig_dir <- file.path(root, "results/figures")
dir.create(dec_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

samples <- read.delim(file.path(root, "data/metadata/samples.tsv"), stringsAsFactors = FALSE)
samples$condition <- factor(samples$condition, levels = c("Control", "IPF"))
samples$project   <- factor(samples$project)
rownames(samples) <- samples$sample_id

gtf <- import(file.path(root, "data/reference/Homo_sapiens.GRCh38.109.gtf"))
gtf_genes <- gtf[gtf$type == "gene"]
ensg2sym <- data.frame(
  ensg = sub("\\..*$", "", as.character(gtf_genes$gene_id)),
  sym  = as.character(gtf_genes$gene_name),
  stringsAsFactors = FALSE
)
ensg2sym <- ensg2sym[!duplicated(ensg2sym$ensg), ]

to_sym_mat <- function(dt, ensg2sym, sample_order) {
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- dt$gene
  mat <- mat[, sample_order, drop = FALSE]
  idx <- match(rownames(mat), ensg2sym$ensg)
  keep <- !is.na(idx) & !is.na(ensg2sym$sym[idx]) & ensg2sym$sym[idx] != ""
  mat <- mat[keep, , drop = FALSE]
  rownames(mat) <- ensg2sym$sym[idx[keep]]
  if (any(duplicated(rownames(mat)))) {
    grp <- rownames(mat)
    mat <- rowsum(mat, group = grp) / as.vector(table(grp))
  }
  mat
}

vst_raw_sym <- to_sym_mat(fread(file.path(tab_dir, "vst_counts.csv")), ensg2sym, samples$sample_id)
vst_cb_sym  <- to_sym_mat(fread(file.path(tab_dir, "vst_counts_batchcorrected.csv")), ensg2sym, samples$sample_id)
vst_limma   <- removeBatchEffect(vst_raw_sym, batch = samples$project)

# Curated macrophage signatures (MERTK replaces myofibroblast marker CTHRC1)
markers <- list(
  M1 = c("CXCL9", "CXCL10", "CXCL11", "IL12B", "IL23A", "NOS2", "TNF", "IL6", "CD80", "CD86"),
  M2 = c("CD163", "MRC1", "MSR1", "CD209", "IL10", "TGFB1", "CCL17", "CCL18", "CCL22", "MMP12"),
  IPF_Mac = c("SPP1", "CHIT1", "CHI3L1", "MERTK", "TREM2", "APOE", "GPNMB", "LPL")
)

score_sig <- function(mat, genes) colMeans(mat[intersect(genes, rownames(mat)), , drop = FALSE])

scores <- data.frame(
  sample_id     = samples$sample_id,
  condition     = samples$condition,
  project       = samples$project,
  M1            = score_sig(vst_cb_sym, markers$M1),
  M2            = score_sig(vst_cb_sym, markers$M2),
  IPF_Mac       = score_sig(vst_cb_sym, markers$IPF_Mac),
  M1_limma      = score_sig(vst_limma, markers$M1),
  M2_limma      = score_sig(vst_limma, markers$M2),
  IPF_Mac_limma = score_sig(vst_limma, markers$IPF_Mac),
  stringsAsFactors = FALSE
)
scores$M2_M1_diff <- scores$M2 - scores$M1

# xCell cell-type deconvolution on Salmon TPM
if (requireNamespace("immunedeconv", quietly = TRUE)) {
  message("[04] Running xCell deconvolution on gene-level TPM...")
  tpm_sym <- to_sym_mat(fread(file.path(tab_dir, "tpm_matrix.csv")), ensg2sym, samples$sample_id)
  xcell_res <- immunedeconv::deconvolute(tpm_sym, method = "xcell")
  xcell_df <- as.data.frame(xcell_res)
  rownames(xcell_df) <- xcell_df$cell_type
  fwrite(xcell_df, file.path(dec_dir, "xcell_cell_fractions.csv"))

  for (ct in c("Macrophages", "Macrophages M1", "Macrophages M2")) {
    if (ct %in% rownames(xcell_df)) {
      scores[[paste0("xCell_", gsub(" ", "_", ct))]] <- as.numeric(xcell_df[ct, samples$sample_id])
    }
  }
} else {
  message("[04] NOTE: immunedeconv not installed; skipping xCell deconvolution.")
}

fwrite(scores, file.path(dec_dir, "m1_m2_marker_scores_v2.csv"))

# Print Table 2 statistics
message("\n=== Table 2 Summary (ComBat VST & Unbiased limma Check) ===")
for (v in c("M1", "M2", "IPF_Mac", "M2_M1_diff", "M1_limma", "M2_limma", "IPF_Mac_limma")) {
  wt <- wilcox.test(scores[[v]] ~ scores$condition)
  message(sprintf("  %-15s Control=%.3f  IPF=%.3f  p=%.4e",
                  v, mean(scores[[v]][scores$condition == "Control"]),
                  mean(scores[[v]][scores$condition == "IPF"]), wt$p.value))
}
if ("xCell_Macrophages_M2" %in% colnames(scores)) {
  for (v in c("xCell_Macrophages_M1", "xCell_Macrophages_M2")) {
    wt <- wilcox.test(scores[[v]] ~ scores$condition)
    message(sprintf("  %-22s Control=%.4f  IPF=%.4f  p=%.4e",
                    v, mean(scores[[v]][scores$condition == "Control"]),
                    mean(scores[[v]][scores$condition == "IPF"]), wt$p.value))
  }
}

# Plot Fig 3A-D
plot_box <- function(df, y, title, out_base) {
  p <- ggplot(df, aes(condition, .data[[y]], fill = condition)) +
    geom_boxplot(alpha = 0.75, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 1.6, alpha = 0.65) +
    scale_fill_manual(values = c(Control = "#4393c3", IPF = "#d6604d")) +
    stat_compare_means(method = "wilcox.test", label.y.npc = 0.95, size = 3.5) +
    labs(title = title, x = NULL, y = "Signature score (VST)") +
    theme_bw(base_size = 12) + theme(plot.title = element_text(face = "bold"), legend.position = "none")
  ggsave(paste0(out_base, ".pdf"), p, width = 4.8, height = 4.8)
  ggsave(paste0(out_base, ".png"), p, width = 4.8, height = 4.8, dpi = 300)
}

plot_box(scores, "M1",         "A. M1 macrophage score",           file.path(fig_dir, "Fig3A_M1_score"))
plot_box(scores, "M2",         "B. M2 macrophage score",           file.path(fig_dir, "Fig3B_M2_score"))
plot_box(scores, "IPF_Mac",    "C. IPF-specific M2 macrophage score", file.path(fig_dir, "Fig3C_IPF_Mac_score"))
plot_box(scores, "M2_M1_diff", "D. M2 - M1 polarization balance",  file.path(fig_dir, "Fig3D_M2_M1_balance"))

# Supp Fig S2A-B (Scatter & Spearman Correlation)
p_sc <- ggplot(scores, aes(M1, M2, color = condition)) +
  geom_point(size = 2.4, alpha = 0.85) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c(Control = "#4393c3", IPF = "#d6604d")) +
  labs(title = "A. M2 vs M1 score per sample", x = "M1 score", y = "M2 score") +
  theme_bw(base_size = 11)
ggsave(file.path(fig_dir, "SuppFigS2A_M2_vs_M1_scatter.png"), p_sc, width = 5.5, height = 4.5, dpi = 300)

cormat <- cor(as.matrix(scores[, c("M1", "M2", "IPF_Mac", "M2_M1_diff")]), method = "spearman")
cm_long <- melt(cormat); colnames(cm_long) <- c("v1", "v2", "rho")
p_cor <- ggplot(cm_long, aes(v1, v2, fill = rho)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", rho)), fontface = "bold", size = 3.8) +
  scale_fill_gradient2(low = "#b2182b", mid = "#f7f7f7", high = "#2166ac", midpoint = 0, limits = c(-1, 1)) +
  labs(title = "B. Spearman correlation of macrophage scores", x = NULL, y = NULL) +
  theme_minimal(base_size = 11) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(fig_dir, "SuppFigS2B_macrophage_cor.png"), p_cor, width = 5.2, height = 4.8, dpi = 300)

message("[04] Complete. Fig 3 and Supp Fig S2 written to ", fig_dir)

# --- Additional Cell-Fraction Deconvolution via quantiseqr (M1 vs M2 fractions) ---
if (requireNamespace("quantiseqr", quietly = TRUE)) {
  message("[04] Running quanTIseq cell-fraction deconvolution...")
  tpm_sym <- to_sym_mat(fread(file.path(tab_dir, "tpm_matrix.csv")), ensg2sym, samples$sample_id)
  qres <- quantiseqr::run_quantiseq(
    expression_data = tpm_sym,
    signature_matrix = "TIL10",
    is_arraydata = FALSE,
    is_tumordata = FALSE,
    scale_mRNA = TRUE
  )
  fwrite(qres, file.path(dec_dir, "quantiseq_cell_fractions.csv"))
  scores$qTI_Macrophages_M1 <- qres[match(scores$sample_id, qres$Sample), "Macrophages.M1"]
  scores$qTI_Macrophages_M2 <- qres[match(scores$sample_id, qres$Sample), "Macrophages.M2"]
  fwrite(scores, file.path(dec_dir, "m1_m2_marker_scores_v2.csv"))

  message("\n=== quanTIseq Cell Fractions (Control vs IPF) ===")
  for (v in c("qTI_Macrophages_M1", "qTI_Macrophages_M2")) {
    wt <- wilcox.test(scores[[v]] ~ scores$condition)
    message(sprintf("  %-22s Control=%.4f  IPF=%.4f  p=%.4e",
                    v, mean(scores[[v]][scores$condition == "Control"], na.rm=TRUE),
                    mean(scores[[v]][scores$condition == "IPF"], na.rm=TRUE), wt$p.value))
  }
}
