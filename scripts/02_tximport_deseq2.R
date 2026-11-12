#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(sva)
  library(rtracklayer)
  library(data.table)
})

set.seed(42)

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
root <- if (length(script_path) > 0) {
  normalizePath(file.path(dirname(script_path), ".."))
} else {
  getwd()
}

meta_file <- file.path(root, "data/metadata/samples.tsv")
tx2gene_f <- file.path(root, "data/reference/tx2gene.tsv")
gtf_file  <- file.path(root, "data/reference/Homo_sapiens.GRCh38.109.gtf")
tab_dir   <- file.path(root, "results/tables")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

samples <- read.delim(meta_file, stringsAsFactors = FALSE)
samples$condition <- factor(samples$condition, levels = c("Control", "IPF"))
samples$project   <- factor(samples$project)
samples$batch     <- factor(samples$batch)

tx2gene <- read.delim(tx2gene_f, header = FALSE, stringsAsFactors = FALSE)
colnames(tx2gene) <- c("tx", "gene")

files <- file.path(root, samples$quant_path)
names(files) <- samples$sample_id
stopifnot(all(file.exists(files)))

message("[02] Running tximport on ", length(files), " samples...")
txi <- tximport(files, type = "salmon", tx2gene = tx2gene,
                ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")

fwrite(data.frame(gene = rownames(txi$abundance), txi$abundance, check.names = FALSE),
       file.path(tab_dir, "tpm_matrix.csv"))

coldata <- data.frame(
  condition = samples$condition,
  project   = samples$project,
  batch     = samples$batch,
  row.names = samples$sample_id
)

dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ project + condition)
keep <- rowSums(counts(dds) >= 10) >= min(table(coldata$condition))
dds  <- dds[keep, ]
message("[02] Retained ", nrow(dds), " genes after low-count filtering.")

dds <- DESeq(dds, parallel = FALSE)
res <- results(dds, contrast = c("condition", "IPF", "Control"), alpha = 0.05)
res <- res[order(res$padj), ]

message("[02] Mapping Ensembl gene IDs to HGNC symbols...")
gtf <- import(gtf_file)
gtf_genes <- gtf[gtf$type == "gene"]
gtf_df <- data.frame(
  ensg   = sub("\\..*$", "", as.character(gtf_genes$gene_id)),
  symbol = as.character(gtf_genes$gene_name),
  stringsAsFactors = FALSE
)
gtf_df <- gtf_df[!duplicated(gtf_df$ensg), ]

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- merge(res_df, gtf_df, by.x = "gene", by.y = "ensg", all.x = TRUE)
missing_sym <- is.na(res_df$symbol) | res_df$symbol == ""
res_df$symbol[missing_sym] <- res_df$gene[missing_sym]
res_df <- res_df[order(res_df$padj), c("gene", "symbol", setdiff(colnames(res_df), c("gene", "symbol")))]
fwrite(res_df, file.path(tab_dir, "deseq2_results_annotated.csv"))

sig_up   <- subset(res_df, !is.na(padj) & padj < 0.05 & log2FoldChange >=  1)
sig_down <- subset(res_df, !is.na(padj) & padj < 0.05 & log2FoldChange <= -1)
fwrite(sig_up[,   c("gene", "symbol", "log2FoldChange", "padj")], file.path(tab_dir, "DEGs_up_annotated.csv"))
fwrite(sig_down[, c("gene", "symbol", "log2FoldChange", "padj")], file.path(tab_dir, "DEGs_down_annotated.csv"))

vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)
mod <- model.matrix(~ condition, data = coldata)
vst_corrected <- ComBat(dat = vst_mat, batch = coldata$project, mod = mod,
                        par.prior = TRUE, prior.plots = FALSE)

fwrite(data.frame(gene = rownames(vst_mat), vst_mat, check.names = FALSE),
       file.path(tab_dir, "vst_counts.csv"))
fwrite(data.frame(gene = rownames(vst_corrected), vst_corrected, check.names = FALSE),
       file.path(tab_dir, "vst_counts_batchcorrected.csv"))

message("[02] Complete. Up DEGs: ", nrow(sig_up), " | Down DEGs: ", nrow(sig_down))
