#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
root <- if (length(script_path) > 0) normalizePath(file.path(dirname(script_path), "..")) else getwd()

tab_dir  <- file.path(root, "results/tables")
fig_root <- file.path(root, "results/figures")
pdf_dir  <- file.path(fig_root, "pdf")
png_dir  <- file.path(fig_root, "png")
tiff_dir <- file.path(fig_root, "tiff")
for (d in c(pdf_dir, png_dir, tiff_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

export_tri_format <- function(base_name, width_in, height_in, draw_fn) {
  pdf_path  <- file.path(pdf_dir,  paste0(base_name, ".pdf"))
  png_pref  <- file.path(png_dir,  base_name)
  tiff_pref <- file.path(tiff_dir, base_name)

  pdf(pdf_path, width = width_in, height = height_in, useDingbats = FALSE)
  draw_fn()
  dev.off()

  system2("pdftoppm", c("-png", "-r", "300", "-singlefile", pdf_path, png_pref))
  system2("pdftoppm", c("-tiff", "-tiffcompression", "lzw", "-r", "300", "-singlefile", pdf_path, tiff_pref))
  message("[12] Exported ", base_name, " (.pdf, .png, .tiff)")
}

# ==============================================================================
# FIGURE 4: CLEAN GGPLOT2 MULTI-TARGET DOCKING HEATMAP (ZERO CHOPPED LABELS)
# ==============================================================================
dock <- fread(file.path(tab_dir, "Table5_docking_affinities_v2.tsv"))
dock[Compound == "biochanin A", Compound := "Biochanin A"]
dock[, Group := fcase(
  Role == "Lead_Candidate",   "Top 10 Drug-Like Phytochemicals",
  Role == "Clinical_Control", "Clinical Control",
  Role == "Native_Control",   "Co-Crystallized Native Controls"
)]
dock[, Group := factor(Group, levels = c(
  "Top 10 Drug-Like Phytochemicals",
  "Clinical Control",
  "Co-Crystallized Native Controls"
))]

# Order compounds top-to-bottom by Mean_Affinity within each group
dock[, Compound := factor(Compound, levels = rev(dock$Compound))]

targets <- c("PPARG", "CCR2", "HIF1A", "STAT3", "NFKB1")
dock_long <- melt(
  dock,
  id.vars = c("Compound", "Group", "Mean_Affinity"),
  measure.vars = targets,
  variable.name = "Target",
  value.name = "Affinity"
)
dock_long[, Target := factor(Target, levels = targets)]
dock_long[, txt_col := ifelse(Affinity <= -7.8, "white", "grey10")]

p4 <- ggplot(dock_long, aes(x = Target, y = Compound, fill = Affinity)) +
  geom_tile(color = "white", linewidth = 0.9) +
  geom_text(aes(label = sprintf("%.2f", Affinity), color = txt_col),
            fontface = "bold", size = 3.6, show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_gradientn(
    name = "Binding Affinity (kcal/mol)",
    colours = c("#08306B", "#2171B5", "#6BAED6", "#DEEBF7", "#FEE0D2"),
    limits = c(-10.7, -4.2),
    breaks = seq(-10, -5, by = 1)
  ) +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
  scale_x_discrete(position = "top") +
  labs(
    title = "Multi-Target AutoDock Vina Binding Affinities (kcal/mol)",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0, margin = margin(b = 10)),
    axis.text.x.top = element_text(face = "bold", size = 11, color = "grey10"),
    axis.text.y = element_text(face = "bold", size = 10, color = "grey15"),
    strip.background = element_rect(fill = "#F0F4F8", color = "grey40", linewidth = 0.6),
    strip.text.y = element_text(face = "bold", size = 9.5, color = "grey10", angle = 0, hjust = 0),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "grey40", linewidth = 0.7),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 10, vjust = 0.8),
    legend.text = element_text(size = 9.5),
    legend.key.width = unit(1.6, "cm"),
    legend.key.height = unit(0.4, "cm"),
    plot.margin = margin(14, 18, 14, 14)
  )

export_tri_format("Fig4_docking_affinity_heatmap", 9.6, 6.8, function() print(p4))

# ==============================================================================
# FIGURE 5: MODERN TWO-PANEL BINDING POTENCY & ADMET PROFILE (ZERO CLUTTER)
# ==============================================================================
admet <- fread(file.path(tab_dir, "Table6_ADMET_and_druglikeness_v2.tsv"))
admet[Compound == "biochanin A", Compound := "Biochanin A"]

dock_sub <- dock[Role %in% c("Lead_Candidate", "Clinical_Control"),
                 .(Compound = as.character(Compound), Mean_Affinity, Best_Affinity)]
m5 <- merge(admet, dock_sub, by = "Compound")
m5[, Category := ifelse(Compound == "Pirfenidone", "Clinical Control", "Lead Phytochemical")]

# Sort by Mean_Affinity (strongest at top)
m5 <- m5[order(Mean_Affinity)]
m5[, Compound := factor(Compound, levels = rev(m5$Compound))]

# Panel A: Cleveland Dumbbell Plot (Mean vs Best Target Binding Affinity)
p5a <- ggplot(m5, aes(y = Compound)) +
  geom_vline(xintercept = -5.63, linetype = "dashed", color = "#D55E00", linewidth = 0.65) +
  geom_segment(aes(x = Best_Affinity, xend = Mean_Affinity, yend = Compound),
               color = "grey65", linewidth = 1.1) +
  geom_point(aes(x = Best_Affinity, shape = "Strongest Single Target", color = Category), size = 3.5) +
  geom_point(aes(x = Mean_Affinity, shape = "5-Target Mean Affinity", color = Category), size = 3.8) +
  geom_text(aes(x = Best_Affinity - 0.14, label = sprintf("%.2f", Best_Affinity)),
            size = 2.9, fontface = "bold", color = "grey20", hjust = 1) +
  scale_color_manual(
    name = "Compound Class",
    values = c("Lead Phytochemical" = "#009E73", "Clinical Control" = "#D55E00")
  ) +
  scale_shape_manual(
    name = "Affinity Metric",
    values = c("5-Target Mean Affinity" = 16, "Strongest Single Target" = 18)
  ) +
  scale_x_continuous(limits = c(-9.4, -5.2), breaks = seq(-9, -5, by = 1)) +
  labs(
    title = "A. Multi-target vs. peak binding potency",
    x = "AutoDock Vina affinity (kcal/mol)",
    y = NULL
  ) +
  theme_bw(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 11.5, hjust = 0, margin = margin(b = 8)),
    axis.title.x = element_text(face = "bold", size = 10.5),
    axis.text.y = element_text(face = "bold", size = 9.8, color = "grey15"),
    axis.text.x = element_text(size = 9.5, color = "grey15"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey92"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(face = "bold", size = 9.5),
    legend.text = element_text(size = 9),
    plot.margin = margin(12, 12, 12, 12)
  )

# Panel B: Clean Multi-Parameter ADMET & Drug-Likeness Dot/Tile Scorecard
admet_metrics <- rbind(
  m5[, .(Compound, Metric = "QED\n(0-1)",      Val_Text = sprintf("%.2f", QED),      Score = QED)],
  m5[, .(Compound, Metric = "MW\n(Da)",        Val_Text = sprintf("%.0f", MW),       Score = 1 - abs(MW - 300)/300)],
  m5[, .(Compound, Metric = "LogP",            Val_Text = sprintf("%.2f", LogP),     Score = 1 - abs(LogP - 2.8)/4)],
  m5[, .(Compound, Metric = "TPSA\n(A^2)",     Val_Text = sprintf("%.0f", TPSA_A2),  Score = 1 - TPSA_A2/160)],
  m5[, .(Compound, Metric = "HBD / HBA",       Val_Text = paste0(HBD, " / ", HBA),   Score = 0.85)],
  m5[, .(Compound, Metric = "GI\nAbsorp.",     Val_Text = GI_Absorption,             Score = ifelse(GI_Absorption == "High", 0.92, 0.4))],
  m5[, .(Compound, Metric = "PAINS\nAlerts",   Val_Text = ifelse(PAINS_Filter == "0 alerts", "Clean (0)", "Catechol"),
         Score = ifelse(PAINS_Filter == "0 alerts", 0.95, 0.45))]
)
admet_metrics[, Metric := factor(Metric, levels = c(
  "QED\n(0-1)", "MW\n(Da)", "LogP", "TPSA\n(A^2)", "HBD / HBA", "GI\nAbsorp.", "PAINS\nAlerts"
))]

p5b <- ggplot(admet_metrics, aes(x = Metric, y = Compound, fill = Score)) +
  geom_tile(color = "white", linewidth = 0.9) +
  geom_text(aes(label = Val_Text), size = 2.95, fontface = "bold", color = "grey10") +
  scale_fill_gradient(
    name = "Drug-Likeness Favorability",
    low = "#FDE0DD", high = "#A1D99B", limits = c(0.25, 1.0),
    breaks = c(0.35, 0.65, 0.95),
    labels = c("Moderate", "Favorable", "Optimal")
  ) +
  scale_x_discrete(position = "top") +
  labs(
    title = "B. Physicochemical & ADMET compliance scorecard",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 11.5, hjust = 0, margin = margin(b = 8)),
    axis.text.x.top = element_text(face = "bold", size = 9.2, color = "grey10"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 9.5),
    legend.text = element_text(size = 9),
    legend.key.width = unit(1.3, "cm"),
    legend.key.height = unit(0.38, "cm"),
    plot.margin = margin(12, 14, 12, 4)
  )

export_tri_format("Fig5_ADMET_landscape", 13.4, 6.2, function() {
  grid.arrange(p5a, p5b, ncol = 2, widths = c(1.05, 1.15))
})

message("[12] Figure 4 and Figure 5 regenerated cleanly in PDF, PNG (300 DPI), and TIFF (300 DPI LZW)!")
