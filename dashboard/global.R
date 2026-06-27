# ─────────────────────────────────────────
# global.R
# Runs once at app startup. Loads libraries
# and data shared across the entire dashboard.
# ─────────────────────────────────────────

library(shiny)
library(tidyverse)
library(plotly)
library(pheatmap)
library(DT)
library(shinydashboard)
library(org.Dm.eg.db)

# ── Load pipeline outputs ──────────────────
# Pre-computed CSVs from the DESeq2 pipeline rule
# Stored in dashboard/data/ for local development
# and shinyapps.io deployment
deseq2_results <- read.csv(
    "data/deseq2_results.csv"
)

normalized_counts <- read.csv(
    "data/normalized_counts.csv",
    row.names = 1
)

# ── Derived objects ────────────────────────
# Significant genes at FDR < 0.05
sig_genes <- deseq2_results %>%
  filter(!is.na(padj), padj < 0.05)

# Top 50 most variable genes for heatmap
top_variable <- normalized_counts %>%
  dplyr::select(starts_with("GSM")) %>%
  mutate(variance = apply(across(everything()), 1, var)) %>%
  arrange(desc(variance)) %>%
  slice_head(n = 50) %>%
  dplyr::select(-variance)

# ── Gene symbol mapping ────────────────────
# Map FlyBase IDs to readable gene symbols
# Available app-wide for all visualizations
deseq2_results$gene_symbol <- mapIds(
  org.Dm.eg.db,
  keys      = deseq2_results$gene_id,
  column    = "SYMBOL",
  keytype   = "FLYBASE",
  multiVals = "first"
)

deseq2_results$label <- ifelse(
  is.na(deseq2_results$gene_symbol),
  deseq2_results$gene_id,
  deseq2_results$gene_symbol
)

# ── Sample metadata ────────────────────────
# Condition labels for heatmap annotation
# Update if using different samples/conditions
sample_info <- data.frame(
  condition = c("untreated", "untreated", "treated", "treated"),
  row.names = c("GSM461177", "GSM461178", "GSM461180", "GSM461181")
)

# ── Gene symbol mapping for normalized counts ──
# Map FBgn row names to gene symbols for heatmap
normalized_counts$gene_symbol <- mapIds(
  org.Dm.eg.db,
  keys      = rownames(normalized_counts),
  column    = "SYMBOL",
  keytype   = "FLYBASE",
  multiVals = "first"
)
