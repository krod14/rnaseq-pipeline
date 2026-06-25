# ─────────────────────────────────────────
# Differential Expression Rule
# DESeq2: statistical modeling of count
# data to identify differentially expressed
# genes between conditions
# ─────────────────────────────────────────

rule deseq2:
    """
    Run DESeq2 differential expression analysis
    on the featureCounts count matrix. Produces
    a results table, normalized counts, and
    diagnostic plots. Normalized counts feed
    directly into the Shiny dashboard.
    """
    input:
        counts  = config["paths"]["counts"] + "counts.txt",
        config  = "config/config.yaml"
    output:
        results    = "results/diffexp/deseq2_results.csv",
        norm_counts = "results/diffexp/normalized_counts.csv",
        pca_plot   = "results/diffexp/pca_plot.pdf",
        volcano    = "results/diffexp/volcano_plot.pdf",
        heatmap    = "results/diffexp/heatmap.pdf"
    conda:
        "../envs/diffexp.yaml"
    log:
        "logs/deseq2/deseq2.log"
    params:
        fdr_threshold = config["params"]["deseq2"]["fdr_threshold"],
        lfc_threshold = config["params"]["deseq2"]["lfc_threshold"]
    script:
        "../scripts/deseq2.R"