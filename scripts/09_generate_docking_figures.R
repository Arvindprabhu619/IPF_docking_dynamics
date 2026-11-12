#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
root <- if (length(script_path) > 0) normalizePath(file.path(dirname(script_path), "..")) else getwd()

tab_dir <- file.path(root, "results/tables")
fig_dir <- file.path(root, "results/figures")

# Figure 4: Multi-target Docking Affinity Heatmap
dock <- fread(file.path(tab_dir, "Table5_docking_affinities_v2.tsv"))
mat <- as.matrix(dock[, .(PPARG, STAT3, CCR2, NFKB1, HIF1A)])
rownames(mat) <- dock$Compound

ann_row <- data.frame(
  Category = factor(dock$Role, levels = c("Lead_Candidate", "Clinical_Control", "Native_Control")),
  row.names = dock$Compound
)
ann_colors <- list(
  Category = c(Lead_Candidate = "#1b9e77", Clinical_Control = "#d95f02", Native_Control = "#7570b3")
)

pdf(file.path(fig_dir, "Fig4_docking_affinity_heatmap.pdf"), width = 8.5, height = 6.5)
pheatmap(
  mat,
  color = colorRampPalette(c("#08306b", "#2171b5", "#6baed6", "#deebf7", "#fee0d2"))(100),
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  number_format = "%.2f",
  fontsize_number = 9,
  annotation_row = ann_row,
  annotation_colors = ann_colors,
  border_color = "white",
  main = "Multi-Target AutoDock Vina Binding Affinities (kcal/mol)"
)
dev.off()

# Figure 5: ADMET & Binding Efficiency Landscape (QED vs Mean Docking Affinity + Boiled-Egg TPSA vs LogP)
admet <- fread(file.path(tab_dir, "Table6_ADMET_and_druglikeness_v2.tsv"))
admet[, Type := ifelse(Compound == "Pirfenidone", "Clinical Control", "Lead Phytochemical")]

p5a <- ggplot(admet, aes(x = QED, y = Mean_Vina_kcal_mol, color = Type)) +
  geom_hline(yintercept = -5.63, linetype = "dashed", color = "#d95f02") +
  geom_point(size = 3.5, alpha = 0.9) +
  geom_text_repel(aes(label = Compound), size = 3.2, max.overlaps = 20) +
  scale_color_manual(values = c("Lead Phytochemical" = "#1b9e77", "Clinical Control" = "#d95f02")) +
  labs(title = "A. Quantitative Estimate of Drug-Likeness (QED) vs Mean Affinity",
       x = "QED Score (0-1)", y = "Mean Vina Binding Affinity (kcal/mol)") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

p5b <- ggplot(admet, aes(x = LogP, y = TPSA_A2, color = Type)) +
  annotate("rect", xmin = -0.4, xmax = 5.88, ymin = 0, ymax = 131.6, alpha = 0.08, fill = "#1b9e77") +
  annotate("rect", xmin = 1.0, xmax = 4.5, ymin = 0, ymax = 90.0, alpha = 0.12, fill = "#2171b5") +
  geom_point(size = 3.5, alpha = 0.9) +
  geom_text_repel(aes(label = Compound), size = 3.2, max.overlaps = 20) +
  scale_color_manual(values = c("Lead Phytochemical" = "#1b9e77", "Clinical Control" = "#d95f02")) +
  labs(title = "B. Physicochemical Absorption Space (LogP vs TPSA)",
       x = "Wildman-Crippen LogP", y = "Topological Polar Surface Area (A^2)") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

pdf(file.path(fig_dir, "Fig5_ADMET_landscape.pdf"), width = 11, height = 5.2)
gridExtra::grid.arrange(p5a, p5b, ncol = 2)
dev.off()

message("[09] Complete. Fig 4 and Fig 5 saved to ", fig_dir)
