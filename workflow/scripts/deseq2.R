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
library(org.Dm.eg.db)

# Prevent AnnotationDbi from masking dplyr::select
select <- dplyr::select

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
    dplyr::select(-Chr, -Start, -End, -Strand, -Length)

# Clean up sample names (featureCounts uses full BAM paths)
colnames(counts) <- colnames(counts) %>%
    basename() %>%
    str_remove("_Aligned.sortedByCoord.out.bam")

# ── Build sample metadata ──────────────────
# Read condition labels from config
sample_info <- data.frame(
    condition = unlist(snakemake@config[["samples"]])
)
rownames(sample_info) <- colnames(counts)

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

# Convert to dataframe
res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    arrange(padj)

# Gene symbol mapping
res_df$gene_symbol <- mapIds(
    org.Dm.eg.db,
    keys      = res_df$gene_id,
    column    = "SYMBOL",
    keytype   = "FLYBASE",
    multiVals = "first"
)

res_df$label <- ifelse(
    is.na(res_df$gene_symbol),
    res_df$gene_id,
    res_df$gene_symbol
)

# Save results with gene symbols
write.csv(res_df,
    snakemake@output[["results"]],
    row.names = FALSE
)

# ── Normalized counts ──────────────────────
# Variance stabilizing transformation for
# visualization (not for DE testing)
vst_counts <- vst(dds, blind = FALSE)
norm_df <- assay(vst_counts) %>%
    as.data.frame()

# Clean column names before writing to CSV
colnames(norm_df) <- colnames(norm_df) %>%
    str_extract("GSM[0-9]+")

norm_df <- norm_df %>%
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
    geom_text_repel(aes(label = rownames(pca_data)),
                    size = 3, show.legend = FALSE,
                    box.padding = 0.5,
                    point.padding = 0.3) +
    xlab(paste0("PC1: ", pca_var[1], "% variance")) +
    ylab(paste0("PC2: ", pca_var[2], "% variance")) +
    theme_minimal() +
    ggtitle("PCA of VST-normalized counts") +
    labs(color = "Condition") +
    guides(color = guide_legend(override.aes = list(label = "")))
dev.off()

# Volcano plot: visualizes DE results

# Label only top 10 most significant genes with abs(log2FC) > 2
# Always include ps (pasilla) since it is the gene of interest
top_labels <- res_df %>%
    filter(!is.na(padj),
           padj < snakemake@params[["fdr_threshold"]],
           abs(log2FoldChange) > 2) %>%
    slice_min(padj, n = 10) %>%
    pull(label)

top_labels <- unique(c(top_labels, "ps"))

pdf(snakemake@output[["volcano"]])
EnhancedVolcano(res_df,
    lab             = res_df$label,
    selectLab       = top_labels,
    x               = "log2FoldChange",
    y               = "padj",
    pCutoff         = snakemake@params[["fdr_threshold"]],
    FCcutoff        = snakemake@params[["lfc_threshold"]],
    title           = "Treated vs Untreated",
    subtitle        = "Pasilla knockdown — Drosophila melanogaster",
    drawConnectors  = TRUE,
    widthConnectors = 0.5
)
dev.off()

# Heatmap: top 50 most variable genes
top_genes <- order(rowVars(assay(vst_counts)),
    decreasing = TRUE)[1:50]

heatmap_mat <- assay(vst_counts)[top_genes, ]

# Clean column names
colnames(heatmap_mat) <- colnames(heatmap_mat) %>%
    str_extract("GSM[0-9]+")

# Map FBgn row names to gene symbols
rownames(heatmap_mat) <- mapIds(
    org.Dm.eg.db,
    keys      = rownames(heatmap_mat),
    column    = "SYMBOL",
    keytype   = "FLYBASE",
    multiVals = "first"
) %>% ifelse(is.na(.), rownames(heatmap_mat), .)

pdf(snakemake@output[["heatmap"]])
pheatmap(
    heatmap_mat,
    annotation_col = sample_info,
    show_rownames  = TRUE,
    fontsize_row   = 8,
    scale          = "row",
    angle_col      = 45
)
dev.off()

sink()
sink(type = "message")