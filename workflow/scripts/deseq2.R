# ─────────────────────────────────────────
# DESeq2 Differential Expression Script
# Called by Snakemake via the script directive
# Inputs, outputs, and params are passed
# automatically via the snakemake object
# ─────────────────────────────────────────

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(DESeq2)
library(EnhancedVolcano)
library(tidyverse)
library(pheatmap)

# ── Load count matrix ──────────────────────
# featureCounts output has 6 metadata columns
# before the actual counts begin
counts_raw <- read.table(
    snakemake@input[["counts"]],
    header = TRUE,
    skip = 1,          # skip the featureCounts comment line
    row.names = 1
)

# Drop metadata columns, keep only count columns
counts <- counts_raw %>%
    select(-Chr, -Start, -End, -Strand, -Length)

# Clean up sample names (featureCounts uses full BAM paths)
colnames(counts) <- colnames(counts) %>%
    basename() %>%
    str_remove("_Aligned.sortedByCoord.out.bam")

# ── Build sample metadata ──────────────────
# Read condition labels from config
sample_info <- data.frame(
    sample    = colnames(counts),
    condition = unlist(snakemake@config[["samples"]])
) %>% column_to_rownames("sample")

# ── Run DESeq2 ─────────────────────────────
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = sample_info,
    design    = ~ condition
)

# Filter lowly expressed genes
# Keep genes with at least 10 counts in at least 2 samples
keep <- rowSums(counts(dds) >= 10) >= 2
dds  <- dds[keep, ]

# Run the DESeq2 model
dds <- DESeq(dds)

# ── Extract results ────────────────────────
res <- results(
    dds,
    alpha    = snakemake@params[["fdr_threshold"]],
    lfcThreshold = snakemake@params[["lfc_threshold"]]
)

# Convert to dataframe and save
res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    arrange(padj)

write.csv(res_df,
    snakemake@output[["results"]],
    row.names = FALSE
)

# ── Normalized counts ──────────────────────
# Variance stabilizing transformation for
# visualization (not for DE testing)
vst_counts <- vst(dds, blind = FALSE)
norm_df <- assay(vst_counts) %>%
    as.data.frame() %>%
    rownames_to_column("gene_id")

write.csv(norm_df,
    snakemake@output[["norm_counts"]],
    row.names = FALSE
)

# ── Diagnostic plots ───────────────────────

# PCA plot: checks sample clustering by condition
pca_data <- plotPCA(vst_counts,
    intgroup = "condition",
    returnData = TRUE
)
pca_var <- round(100 * attr(pca_data, "percentVar"))

pdf(snakemake@output[["pca_plot"]])
ggplot(pca_data, aes(PC1, PC2, color = condition)) +
    geom_point(size = 4) +
    xlab(paste0("PC1: ", pca_var[1], "% variance")) +
    ylab(paste0("PC2: ", pca_var[2], "% variance")) +
    theme_minimal() +
    ggtitle("PCA of VST-normalized counts")
dev.off()

# Volcano plot: visualizes DE results
pdf(snakemake@output[["volcano"]])
EnhancedVolcano(res_df,
    lab    = res_df$gene_id,
    x      = "log2FoldChange",
    y      = "padj",
    pCutoff = snakemake@params[["fdr_threshold"]],
    FCcutoff = snakemake@params[["lfc_threshold"]],
    title  = "Treated vs Untreated"
)
dev.off()

# Heatmap: top 50 most variable genes
top_genes <- order(rowVars(assay(vst_counts)),
    decreasing = TRUE)[1:50]

pdf(snakemake@output[["heatmap"]])
pheatmap(
    assay(vst_counts)[top_genes, ],
    annotation_col = sample_info,
    show_rownames  = TRUE,
    scale          = "row"
)
dev.off()

sink()
sink(type = "message")