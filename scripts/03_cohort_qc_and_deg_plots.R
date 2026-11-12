#!/usr/bin/env Rscript
# =============================================================================
#  03_cohort_qc_and_deg_plots.R (v2.0.0)
#  Generates:
#    - Fig 1A-D: PCA & UMAP before and after ComBat batch correction
#    - Supp Fig S1A: t-SNE after batch correction
#    - Fig 2A-C: Annotated Volcano plot, MA plot, and Top 50 Up/Down DEG Heatmap
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(umap)
  library(Rtsne)
  library(pheatmap)
  library(RColorBrewer)
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
fig_dir <- file.path(root, "results/figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

samples <- read.delim(file.path(root, "data/metadata/samples.tsv"), stringsAsFactors = FALSE)
samples$condition <- factor(samples$condition, levels = c("Control", "IPF"))
samples$project   <- factor(samples$project)

to_mat <- function(dt, sample_ids) {
  m <- as.matrix(dt[, -1, with = FALSE])
  rownames(m) <- dt$gene
  m[, sample_ids, drop = FALSE]
}

mat_before <- to_mat(fread(file.path(tab_dir, "vst_counts.csv")), samples$sample_id)
mat_after  <- to_mat(fread(file.path(tab_dir, "vst_counts_batchcorrected.csv")), samples$sample_id)

top_var <- function(m, n = 2000) {
  rv <- apply(m, 1, var)
  m[order(rv, decreasing = TRUE)[seq_len(min(n, length(rv)))], ]
}

# 1. PCA before & after (Fig 1A-B)
plot_pca <- function(m, title, out_base) {
  pca <- prcomp(t(top_var(m)), scale. = TRUE)
  pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                   condition = samples$condition, project = samples$project)
  p <- ggplot(df, aes(PC1, PC2, colour = project, shape = condition)) +
    geom_point(size = 3, alpha = 0.85) +
    labs(title = title, x = paste0("PC1 (", pct[1], "%)"), y = paste0("PC2 (", pct[2], "%)")) +
    theme_bw(base_size = 12) + theme(plot.title = element_text(face = "bold"))
  ggsave(paste0(out_base, ".pdf"), p, width = 6, height = 5)
  ggsave(paste0(out_base, ".png"), p, width = 6, height = 5, dpi = 300)
}

plot_pca(mat_before, "A. PCA before batch correction", file.path(fig_dir, "Fig1A_pca_before"))
plot_pca(mat_after,  "B. PCA after batch correction",  file.path(fig_dir, "Fig1B_pca_after"))

# 2. UMAP before & after (Fig 1C-D)
plot_umap <- function(m, title, out_base) {
  set.seed(42)
  u <- umap(t(top_var(m)), n_neighbors = 15, min_dist = 0.1)
  df <- data.frame(UMAP1 = u$layout[, 1], UMAP2 = u$layout[, 2],
                   condition = samples$condition, project = samples$project)
  p <- ggplot(df, aes(UMAP1, UMAP2, colour = project, shape = condition)) +
    geom_point(size = 3, alpha = 0.85) +
    labs(title = title, x = "UMAP1", y = "UMAP2") +
    theme_bw(base_size = 12) + theme(plot.title = element_text(face = "bold"))
  ggsave(paste0(out_base, ".pdf"), p, width = 6, height = 5)
  ggsave(paste0(out_base, ".png"), p, width = 6, height = 5, dpi = 300)
}

plot_umap(mat_before, "C. UMAP before batch correction", file.path(fig_dir, "Fig1C_umap_before"))
plot_umap(mat_after,  "D. UMAP after batch correction",  file.path(fig_dir, "Fig1D_umap_after"))

# 3. t-SNE after correction (Supp Fig S1A)
set.seed(42)
ts <- Rtsne(t(top_var(mat_after)), perplexity = 20, max_iter = 1000, check_duplicates = FALSE)
df_ts <- data.frame(tSNE1 = ts$Y[, 1], tSNE2 = ts$Y[, 2],
                    condition = samples$condition, project = samples$project)
p_ts <- ggplot(df_ts, aes(tSNE1, tSNE2, colour = project, shape = condition)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(title = "A. t-SNE after batch correction", x = "t-SNE 1", y = "t-SNE 2") +
  theme_bw(base_size = 12) + theme(plot.title = element_text(face = "bold"))
ggsave(file.path(fig_dir, "SuppFigS1A_tsne_after.pdf"), p_ts, width = 6, height = 5)
ggsave(file.path(fig_dir, "SuppFigS1A_tsne_after.png"), p_ts, width = 6, height = 5, dpi = 300)

# 4. Annotated Volcano (Fig 2A) & MA Plot (Fig 2B)
res <- fread(file.path(tab_dir, "deseq2_results_annotated.csv"))[!is.na(padj)]
res[, sig := "NS"]
res[log2FoldChange >=  1 & padj < 0.05, sig := "Up"]
res[log2FoldChange <= -1 & padj < 0.05, sig := "Down"]
res[, sig := factor(sig, levels = c("Down", "NS", "Up"))]

top15 <- res[order(padj)][1:15]
p_volc <- ggplot(res, aes(log2FoldChange, -log10(padj), colour = sig)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_colour_manual(values = c(Down = "steelblue", NS = "grey70", Up = "firebrick")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(data = top15, aes(label = symbol), size = 3, max.overlaps = 30) +
  labs(title = "A. Differential expression (IPF vs Control)",
       x = "log2 fold change", y = "-log10 adjusted p-value") +
  theme_bw(base_size = 12) + theme(legend.title = element_blank())
ggsave(file.path(fig_dir, "Fig2A_volcano.pdf"), p_volc, width = 7, height = 6)
ggsave(file.path(fig_dir, "Fig2A_volcano.png"), p_volc, width = 7, height = 6, dpi = 300)

p_ma <- ggplot(res, aes(log10(baseMean + 1), log2FoldChange, colour = sig)) +
  geom_point(alpha = 0.45, size = 0.9) +
  scale_colour_manual(values = c(Down = "steelblue", NS = "grey70", Up = "firebrick")) +
  geom_hline(yintercept = c(-1, 0, 1), linetype = c("dashed", "solid", "dashed")) +
  ylim(-6, 6) +
  labs(title = "B. MA plot (IPF vs Control)", x = "log10(mean normalized counts + 1)", y = "log2 fold change") +
  theme_bw(base_size = 12) + theme(legend.title = element_blank())
ggsave(file.path(fig_dir, "Fig2B_ma_plot.pdf"), p_ma, width = 7, height = 5.5)
ggsave(file.path(fig_dir, "Fig2B_ma_plot.png"), p_ma, width = 7, height = 5.5, dpi = 300)

# 5. Top 50 Up + 50 Down Heatmap (Fig 2C)
res$label <- make.unique(res$symbol)
top100 <- rbind(res[sig == "Up"][order(padj)][1:50], res[sig == "Down"][order(padj)][1:50])
mat_hm <- mat_after[top100$gene, samples$sample_id, drop = FALSE]
rownames(mat_hm) <- top100$label
mat_z <- t(scale(t(mat_hm)))
mat_z[mat_z > 2] <- 2; mat_z[mat_z < -2] <- -2

ann_col <- data.frame(Condition = samples$condition, Project = samples$project, row.names = samples$sample_id)
proj_lvls <- levels(samples$project)
ann_colors <- list(
  Condition = c(Control = "#4393c3", IPF = "#d6604d"),
  Project   = setNames(colorRampPalette(brewer.pal(8, "Set2"))(length(proj_lvls)), proj_lvls)
)

for (ext in c("png", "pdf")) {
  out_f <- file.path(fig_dir, paste0("Fig2C_top100_DEG_heatmap.", ext))
  if (ext == "png") png(out_f, width = 12, height = 11, units = "in", res = 300) else pdf(out_f, width = 12, height = 11)
  pheatmap(mat_z, color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           annotation_col = ann_col, annotation_colors = ann_colors,
           show_rownames = TRUE, show_colnames = FALSE, fontsize_row = 6, border_color = NA,
           main = "C. Top 50 upregulated and 50 downregulated DEGs")
  dev.off()
}

message("[03] Complete. Figures 1, 2, and Supp Fig S1A saved to ", fig_dir)
