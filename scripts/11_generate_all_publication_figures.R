#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(gridExtra)
  library(grid)
  library(umap)
  library(Rtsne)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
root <- if (length(script_path) > 0) normalizePath(file.path(dirname(script_path), "..")) else getwd()

meta_path <- file.path(root, "data/metadata/samples.tsv")
tab_dir   <- file.path(root, "results/tables")
fig_root  <- file.path(root, "results/figures")
pdf_dir   <- file.path(fig_root, "pdf")
png_dir   <- file.path(fig_root, "png")
tiff_dir  <- file.path(fig_root, "tiff")
for (d in c(pdf_dir, png_dir, tiff_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Export helper: PDF + 300 DPI PNG + 300 DPI LZW-compressed TIFF via pdftoppm
export_tri_format <- function(base_name, width_in, height_in, draw_fn) {
  pdf_path  <- file.path(pdf_dir,  paste0(base_name, ".pdf"))
  png_pref  <- file.path(png_dir,  base_name)
  tiff_pref <- file.path(tiff_dir, base_name)

  pdf(pdf_path, width = width_in, height = height_in, useDingbats = FALSE)
  draw_fn()
  dev.off()

  system2("pdftoppm", c("-png", "-r", "300", "-singlefile", pdf_path, png_pref))
  system2("pdftoppm", c("-tiff", "-tiffcompression", "lzw", "-r", "300", "-singlefile", pdf_path, tiff_pref))
  message("[11] Exported ", base_name, " (.pdf, .png, .tiff)")
}

pub_theme <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0, margin = margin(b = 8)),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(color = "grey15", size = base_size - 1),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size - 1),
      legend.key.height = unit(0.45, "cm"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey30", linewidth = 0.7),
      plot.margin = margin(12, 14, 12, 12)
    )
}

cohort_cols <- c(
  "PRJNA229054" = "#E69F00",
  "PRJNA358081" = "#009E73",
  "PRJNA388978" = "#0072B2",
  "PRJNA575200" = "#CC79A7"
)
cond_cols <- c("Control" = "#377EB8", "IPF" = "#E41A1C")

# Load metadata & VST matrices
meta <- fread(meta_path)
meta[, condition := factor(condition, levels = c("Control", "IPF"))]
meta[, project := factor(project, levels = sort(unique(project)))]

vst_b_dt <- fread(file.path(tab_dir, "vst_counts.csv"))
vst_before <- as.matrix(vst_b_dt[, -1, with = FALSE])
rownames(vst_before) <- vst_b_dt[[1]]
vst_before <- vst_before[, meta$sample_id]

vst_a_dt <- fread(file.path(tab_dir, "vst_counts_batchcorrected.csv"))
vst_after <- as.matrix(vst_a_dt[, -1, with = FALSE])
rownames(vst_after) <- vst_a_dt[[1]]
vst_after <- vst_after[, meta$sample_id]

# ==============================================================================
# FIGURE 1A–D & SUPP FIG S1A: PCA, UMAP & t-SNE
# ==============================================================================
plot_embedding <- function(df, xcol, ycol, xlab, ylab, title_str) {
  ggplot(df, aes(x = .data[[xcol]], y = .data[[ycol]], color = project, shape = condition)) +
    geom_point(size = 3.2, alpha = 0.88, stroke = 0.4) +
    scale_color_manual(name = "Cohort", values = cohort_cols) +
    scale_shape_manual(name = "Condition", values = c("Control" = 16, "IPF" = 17)) +
    guides(
      color = guide_legend(order = 1, override.aes = list(size = 3.8)),
      shape = guide_legend(order = 2, override.aes = list(size = 3.8))
    ) +
    labs(title = title_str, x = xlab, y = ylab) +
    pub_theme(11)
}

top_var_b <- head(order(apply(vst_before, 1, var), decreasing = TRUE), 2000)
pca_b <- prcomp(t(vst_before[top_var_b, ]), scale. = TRUE)
pve_b <- round(100 * (pca_b$sdev^2 / sum(pca_b$sdev^2)), 1)
df_pca_b <- cbind(meta, PC1 = pca_b$x[, 1], PC2 = pca_b$x[, 2])

top_var_a <- head(order(apply(vst_after, 1, var), decreasing = TRUE), 2000)
pca_a <- prcomp(t(vst_after[top_var_a, ]), scale. = TRUE)
pve_a <- round(100 * (pca_a$sdev^2 / sum(pca_a$sdev^2)), 1)
df_pca_a <- cbind(meta, PC1 = pca_a$x[, 1], PC2 = pca_a$x[, 2])

export_tri_format("Fig1A_pca_before", 6.8, 5.2, function() {
  print(plot_embedding(df_pca_b, "PC1", "PC2",
                       paste0("PC1 (", pve_b[1], "%)"), paste0("PC2 (", pve_b[2], "%)"),
                       "A. PCA before batch correction"))
})
export_tri_format("Fig1B_pca_after", 6.8, 5.2, function() {
  print(plot_embedding(df_pca_a, "PC1", "PC2",
                       paste0("PC1 (", pve_a[1], "%)"), paste0("PC2 (", pve_a[2], "%)"),
                       "B. PCA after batch correction"))
})

set.seed(42)
um_b <- umap(t(vst_before[top_var_b, ]))
df_um_b <- cbind(meta, UMAP1 = um_b$layout[, 1], UMAP2 = um_b$layout[, 2])
export_tri_format("Fig1C_umap_before", 6.8, 5.2, function() {
  print(plot_embedding(df_um_b, "UMAP1", "UMAP2", "UMAP 1", "UMAP 2",
                       "C. UMAP before batch correction"))
})

set.seed(42)
um_a <- umap(t(vst_after[top_var_a, ]))
df_um_a <- cbind(meta, UMAP1 = um_a$layout[, 1], UMAP2 = um_a$layout[, 2])
export_tri_format("Fig1D_umap_after", 6.8, 5.2, function() {
  print(plot_embedding(df_um_a, "UMAP1", "UMAP2", "UMAP 1", "UMAP 2",
                       "D. UMAP after batch correction"))
})

set.seed(42)
tsne_res <- Rtsne(t(vst_after[top_var_a, ]), perplexity = 20, check_duplicates = FALSE)
df_tsne <- cbind(meta, tSNE1 = tsne_res$Y[, 1], tSNE2 = tsne_res$Y[, 2])
export_tri_format("SuppFigS1A_tsne_after", 6.8, 5.2, function() {
  print(plot_embedding(df_tsne, "tSNE1", "tSNE2", "t-SNE 1", "t-SNE 2",
                       "A. t-SNE after batch correction"))
})

# ==============================================================================
# FIGURE 2A–C: VOLCANO, MA PLOT & TOP 100 DEG HEATMAP
# ==============================================================================
deg <- fread(file.path(tab_dir, "deseq2_results_annotated.csv"))
deg <- deg[!is.na(padj) & !is.na(log2FoldChange)]
deg[, label_gene := ifelse(!is.na(symbol) & symbol != "", symbol, gene)]
deg[, Status := ifelse(padj < 0.05 & log2FoldChange >= 1, "Upregulated",
                ifelse(padj < 0.05 & log2FoldChange <= -1, "Downregulated", "Not Significant"))]
deg[, Status := factor(Status, levels = c("Downregulated", "Not Significant", "Upregulated"))]
deg[, negLog10P := -log10(pmax(padj, 1e-300))]

top_up   <- deg[Status == "Upregulated"][order(padj)][1:min(12, .N)]
top_down <- deg[Status == "Downregulated"][order(padj)][1:min(6, .N)]
lab_deg  <- rbind(top_up, top_down)

p2a <- ggplot(deg, aes(x = log2FoldChange, y = negLog10P, color = Status)) +
  geom_point(alpha = 0.65, size = 1.6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey30", linewidth = 0.55) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey30", linewidth = 0.55) +
  geom_label_repel(
    data = lab_deg,
    aes(label = label_gene),
    size = 2.9,
    fontface = "bold",
    fill = alpha("white", 0.85),
    label.size = 0.15,
    box.padding = 0.45,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.color = "grey40",
    segment.size = 0.35,
    max.overlaps = 30,
    seed = 42,
    show.legend = FALSE
  ) +
  scale_color_manual(
    name = "Regulation",
    values = c("Downregulated" = "#2B6CB0", "Not Significant" = "grey75", "Upregulated" = "#C53030")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3.2, alpha = 1))) +
  scale_x_continuous(expand = expansion(mult = 0.08)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(title = "A. Differential expression (IPF vs Control)",
       x = expression(bold(log[2]~"fold change")),
       y = expression(bold(-log[10]~"adjusted p-value"))) +
  pub_theme(11)

export_tri_format("Fig2A_volcano", 7.4, 5.6, function() print(p2a))

p2b <- ggplot(deg, aes(x = log10(baseMean + 1), y = log2FoldChange, color = Status)) +
  geom_point(alpha = 0.65, size = 1.6) +
  geom_hline(yintercept = 0, color = "grey20", linewidth = 0.6) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "grey30", linewidth = 0.55) +
  scale_color_manual(
    name = "Regulation",
    values = c("Downregulated" = "#2B6CB0", "Not Significant" = "grey75", "Upregulated" = "#C53030")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3.2, alpha = 1))) +
  labs(title = "B. MA plot (IPF vs Control)",
       x = expression(bold(log[10]~"(mean normalized counts + 1)")),
       y = expression(bold(log[2]~"fold change"))) +
  pub_theme(11)

export_tri_format("Fig2B_ma_plot", 7.4, 5.6, function() print(p2b))

# Fig 2C: Top 50 Up + Top 50 Down Heatmap with clean gene symbols
sig_up_dt   <- deg[Status == "Upregulated"][order(padj)][1:min(50, .N)]
sig_down_dt <- deg[Status == "Downregulated"][order(padj)][1:min(50, .N)]
top100_dt   <- rbind(sig_up_dt, sig_down_dt)
top100_dt   <- top100_dt[gene %in% rownames(vst_after)]

ord_idx <- order(meta$condition, meta$project)
meta_ord <- meta[ord_idx]
mat_h <- vst_after[top100_dt$gene, meta_ord$sample_id]
rownames(mat_h) <- make.unique(top100_dt$label_gene)
mat_z <- t(scale(t(mat_h)))
mat_z[mat_z > 2.5] <- 2.5
mat_z[mat_z < -2.5] <- -2.5

ann_col <- data.frame(
  Condition = meta_ord$condition,
  Cohort = meta_ord$project,
  row.names = meta_ord$sample_id
)
ann_col_colors <- list(
  Condition = cond_cols,
  Cohort = cohort_cols
)

export_tri_format("Fig2C_top100_DEG_heatmap", 9.5, 13.5, function() {
  pheatmap(
    mat_z,
    color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    fontsize_row = 6.8,
    annotation_col = ann_col,
    annotation_colors = ann_col_colors,
    border_color = NA,
    main = "C. Top 50 upregulated and 50 downregulated DEGs"
  )
})

# ==============================================================================
# FIGURE 3A–D & SUPP FIG S2A–B: MACROPHAGE POLARIZATION SCORES
# ==============================================================================
# Build symbol-indexed VST matrix using deseq2 annotation mapping
sym_map <- unique(deg[!is.na(symbol) & symbol != "", .(gene, symbol)])
sym_map <- sym_map[gene %in% rownames(vst_after)]
vst_sym <- vst_after[sym_map$gene, , drop = FALSE]
rownames(vst_sym) <- make.unique(sym_map$symbol)

score_sig <- function(genes) {
  g_use <- intersect(genes, rownames(vst_sym))
  colMeans(vst_sym[g_use, , drop = FALSE], na.rm = TRUE)
}

m1_genes <- c("IL1B", "IL6", "TNF", "NOS2", "CXCL9", "CXCL10", "CXCL11", "CD80", "CD86", "CCR7", "IRF5", "STAT1")
m2_genes <- c("MRC1", "CD163", "MSR1", "ARG1", "IL10", "TGFB1", "CCL17", "CCL18", "CCL22", "PPARG", "STAT6", "IRF4", "CHI3L1", "MMP9")
ipf_mac_genes <- c("SPP1", "MERTK", "LGMN", "TREM2", "GPNMB", "FABP4", "CCL18", "MMP9", "CHIT1", "CD163", "MRC1", "FN1")

mac <- copy(meta)
mac[, M1 := score_sig(m1_genes)]
mac[, M2 := score_sig(m2_genes)]
mac[, IPF_Mac := score_sig(ipf_mac_genes)]
mac[, M2_M1_diff := M2 - M1]

plot_mac_box <- function(col_name, title_str, ylab_str) {
  y_vals <- mac[[col_name]]
  wt <- wilcox.test(y_vals ~ mac$condition)
  p_txt <- if (wt$p.value < 0.001) {
    paste0("Wilcoxon, p = ", formatC(wt$p.value, format = "e", digits = 1))
  } else {
    paste0("Wilcoxon, p = ", signif(wt$p.value, 2))
  }
  y_max <- max(y_vals, na.rm = TRUE)
  y_min <- min(y_vals, na.rm = TRUE)
  span  <- y_max - y_min
  y_bar <- y_max + 0.06 * span
  y_lbl <- y_max + 0.13 * span

  ggplot(mac, aes(x = condition, y = .data[[col_name]], fill = condition)) +
    geom_boxplot(width = 0.55, alpha = 0.85, outlier.shape = NA, color = "grey20", linewidth = 0.65) +
    geom_jitter(width = 0.14, size = 2.1, alpha = 0.72, color = "grey15") +
    annotate("segment", x = 1, xend = 2, y = y_bar, yend = y_bar, linewidth = 0.55, color = "grey20") +
    annotate("segment", x = 1, xend = 1, y = y_bar - 0.02 * span, yend = y_bar, linewidth = 0.55, color = "grey20") +
    annotate("segment", x = 2, xend = 2, y = y_bar - 0.02 * span, yend = y_bar, linewidth = 0.55, color = "grey20") +
    annotate("text", x = 1.5, y = y_lbl, label = p_txt, size = 3.8, fontface = "bold") +
    scale_fill_manual(values = cond_cols) +
    scale_y_continuous(limits = c(y_min - 0.04 * span, y_max + 0.21 * span)) +
    labs(title = title_str, x = NULL, y = ylab_str) +
    pub_theme(11) +
    theme(legend.position = "none", axis.text.x = element_text(face = "bold", size = 11))
}

export_tri_format("Fig3A_M1_score", 5.2, 5.2, function() {
  print(plot_mac_box("M1", "A. M1 macrophage score", "Signature score (VST)"))
})
export_tri_format("Fig3B_M2_score", 5.2, 5.2, function() {
  print(plot_mac_box("M2", "B. M2 macrophage score", "Signature score (VST)"))
})
export_tri_format("Fig3C_IPF_Mac_score", 5.2, 5.2, function() {
  print(plot_mac_box("IPF_Mac", "C. IPF-specific M2 macrophage score", "Signature score (VST)"))
})
export_tri_format("Fig3D_M2_M1_balance", 5.2, 5.2, function() {
  print(plot_mac_box("M2_M1_diff", "D. M2 - M1 polarization balance", "Polarization score (M2 - M1)"))
})

p_s2a <- ggplot(mac, aes(x = M1, y = M2, color = condition)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.6) +
  geom_point(size = 3.0, alpha = 0.85) +
  scale_color_manual(name = "Condition", values = cond_cols) +
  guides(color = guide_legend(override.aes = list(size = 3.6))) +
  labs(title = "A. M2 vs M1 signature score per sample", x = "M1 macrophage score", y = "M2 macrophage score") +
  pub_theme(11)

export_tri_format("SuppFigS2A_M2_vs_M1_scatter", 6.4, 5.2, function() print(p_s2a))

cor_mat <- cor(mac[, .(M1, M2, IPF_Mac, M2_M1_diff)], method = "spearman")
clean_names <- c("M1 Score", "M2 Score", "IPF-Mac Score", "M2 - M1 Balance")
dimnames(cor_mat) <- list(clean_names, clean_names)
cor_dt <- as.data.table(as.table(cor_mat))
setnames(cor_dt, c("Var1", "Var2", "rho"))
cor_dt[, Var1 := factor(Var1, levels = clean_names)]
cor_dt[, Var2 := factor(Var2, levels = clean_names)]

p_s2b <- ggplot(cor_dt, aes(x = Var1, y = Var2, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f", rho)), fontface = "bold", size = 4.2) +
  scale_fill_gradient2(name = expression(bold(rho)), low = "#B2182B", mid = "#F7F7F7", high = "#2166AC",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(title = "B. Spearman correlation of macrophage scores", x = NULL, y = NULL) +
  pub_theme(11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"))

export_tri_format("SuppFigS2B_macrophage_cor", 6.4, 5.4, function() print(p_s2b))

# ==============================================================================
# FIGURE 4: MULTI-TARGET DOCKING AFFINITY HEATMAP
# ==============================================================================
dock <- fread(file.path(tab_dir, "Table5_docking_affinities_v2.tsv"))
dock[Compound == "biochanin A", Compound := "Biochanin A"]
dock[, Role_Clean := fcase(
  Role == "Lead_Candidate",   "Lead Phytochemical",
  Role == "Clinical_Control", "Clinical Control",
  Role == "Native_Control",   "Native Co-Crystal Control"
)]

mat_d <- as.matrix(dock[, .(PPARG, STAT3, CCR2, NFKB1, HIF1A)])
rownames(mat_d) <- dock$Compound

ann_row <- data.frame(
  Category = factor(dock$Role_Clean, levels = c("Lead Phytochemical", "Clinical Control", "Native Co-Crystal Control")),
  row.names = dock$Compound
)
ann_colors <- list(
  Category = c(
    "Lead Phytochemical"        = "#009E73",
    "Clinical Control"          = "#D55E00",
    "Native Co-Crystal Control" = "#0072B2"
  )
)

export_tri_format("Fig4_docking_affinity_heatmap", 9.6, 6.8, function() {
  pheatmap(
    mat_d,
    color = colorRampPalette(c("#08306B", "#2171B5", "#6BAED6", "#DEEBF7", "#FEE0D2"))(100),
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    display_numbers = TRUE,
    number_format = "%.2f",
    fontsize_number = 9.5,
    fontsize_row = 10,
    fontsize_col = 10.5,
    angle_col = 0,
    gaps_row = 10,
    annotation_row = ann_row,
    annotation_colors = ann_colors,
    border_color = "white",
    main = "Multi-Target AutoDock Vina Binding Affinities (kcal/mol)"
  )
})

# ==============================================================================
# FIGURE 5: ADMET & BINDING EFFICIENCY LANDSCAPE
# ==============================================================================
admet <- fread(file.path(tab_dir, "Table6_ADMET_and_druglikeness_v2.tsv"))
admet[Compound == "biochanin A", Compound := "Biochanin A"]
admet[, Type := ifelse(Compound == "Pirfenidone", "Clinical Control", "Lead Phytochemical")]

p5a <- ggplot(admet, aes(x = QED, y = Mean_Vina_kcal_mol, color = Type)) +
  geom_hline(yintercept = -5.63, linetype = "dashed", color = "#D55E00", linewidth = 0.65) +
  geom_point(size = 3.8, alpha = 0.92) +
  geom_label_repel(
    aes(label = Compound),
    size = 3.0,
    fontface = "bold",
    fill = alpha("white", 0.85),
    label.size = 0.15,
    box.padding = 0.45,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.color = "grey45",
    max.overlaps = 25,
    seed = 42,
    show.legend = FALSE
  ) +
  scale_color_manual(name = "Category", values = c("Lead Phytochemical" = "#009E73", "Clinical Control" = "#D55E00")) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  scale_x_continuous(limits = c(0.38, 0.80)) +
  scale_y_continuous(limits = c(-7.35, -5.45)) +
  labs(title = "A. Drug-likeness (QED) vs. mean binding affinity",
       x = "Quantitative Estimate of Drug-Likeness (QED)",
       y = "Mean Vina affinity (kcal/mol)") +
  pub_theme(10.5) +
  theme(legend.position = "bottom")

p5b <- ggplot(admet, aes(x = LogP, y = TPSA_A2, color = Type)) +
  annotate("rect", xmin = -0.4, xmax = 5.88, ymin = 0, ymax = 131.6, alpha = 0.08, fill = "#009E73") +
  annotate("rect", xmin = 1.0, xmax = 4.5, ymin = 0, ymax = 90.0, alpha = 0.12, fill = "#0072B2") +
  geom_point(size = 3.8, alpha = 0.92) +
  geom_label_repel(
    aes(label = Compound),
    size = 3.0,
    fontface = "bold",
    fill = alpha("white", 0.85),
    label.size = 0.15,
    box.padding = 0.45,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.color = "grey45",
    max.overlaps = 25,
    seed = 42,
    show.legend = FALSE
  ) +
  scale_color_manual(name = "Category", values = c("Lead Phytochemical" = "#009E73", "Clinical Control" = "#D55E00")) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  scale_x_continuous(limits = c(-0.2, 5.2)) +
  scale_y_continuous(limits = c(0, 142)) +
  labs(title = "B. Gastrointestinal & BBB absorption space",
       x = "Wildman-Crippen LogP",
       y = expression(bold("Topological polar surface area ("*ring(A)^2*")"))) +
  pub_theme(10.5) +
  theme(legend.position = "bottom")

export_tri_format("Fig5_ADMET_landscape", 13.2, 5.6, function() {
  grid.arrange(p5a, p5b, ncol = 2)
})

message("[11] ALL PUBLICATION FIGURES GENERATED IN PDF, PNG (300 DPI), AND TIFF (300 DPI LZW)!")
