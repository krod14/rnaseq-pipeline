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
# These CSVs are produced by the DESeq2 rule
deseq2_results <- read.csv(
    "../results/diffexp/deseq2_results.csv"
)

normalized_counts <- read.csv(
    "../results/diffexp/normalized_counts.csv",
    row.names = 1
)

# ── Derived objects ────────────────────────
# Significant genes at FDR < 0.05
sig_genes <- deseq2_results %>%
    filter(!is.na(padj), padj < 0.05)

# Top 50 most variable genes for heatmap
top_variable <- normalized_counts %>%
    mutate(variance = apply(across(everything()), 1, var)) %>%
    arrange(desc(variance)) %>%
    slice_head(n = 50)