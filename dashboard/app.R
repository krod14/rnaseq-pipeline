# ─────────────────────────────────────────
# app.R
# Shiny dashboard for RNA-seq pipeline results
# Tabs: Overview | PCA | Volcano | Heatmap | Table
# ─────────────────────────────────────────

source("global.R")

# ── UI ─────────────────────────────────────
ui <- dashboardPage(
    dashboardHeader(title = "RNA-seq Results"),

    dashboardSidebar(
        sidebarMenu(
            menuItem("Overview",      tabName = "overview", icon = icon("chart-bar")),
            menuItem("PCA",           tabName = "pca",      icon = icon("circle-nodes")),
            menuItem("Volcano",       tabName = "volcano",  icon = icon("circle-dot")),
            menuItem("Heatmap",       tabName = "heatmap",  icon = icon("th")),
            menuItem("Results Table", tabName = "table",    icon = icon("table"))
        )
    ),

    dashboardBody(
        tabItems(

            # ── Overview Tab ───────────────
            tabItem(tabName = "overview",

                # ── Value Boxes ───────────────
                fluidRow(
                    valueBox(nrow(deseq2_results), "Genes Tested",
                        icon = icon("dna"), color = "purple", width = 3),
                    valueBox(nrow(sig_genes), "Sig. Genes (FDR < 0.05)",
                        icon = icon("star"), color = "green", width = 3),
                    valueBox(
                        sum(sig_genes$log2FoldChange > 0),
                        "Upregulated",
                        icon = icon("arrow-up"), color = "blue", width = 3),
                    valueBox(
                        sum(sig_genes$log2FoldChange < 0),
                        "Downregulated",
                        icon = icon("arrow-down"), color = "red", width = 3)
                ),

                # ── About This Dashboard ───────
                fluidRow(
                    box(
                        title = "About This Dashboard",
                        status = "primary",
                        width = 12,
                        p("This interactive dashboard explores the results of an end-to-end
                          bulk RNA-seq pipeline built with Snakemake, conda, and Docker.
                          It provides interactive visualizations of differential expression
                          results and normalized count data across four tabs: PCA, Volcano,
                          Heatmap, and Results Table."),
                        p("The pipeline and all source code are available on ",
                          a("GitHub", href = "https://github.com/krod14/rnaseq-pipeline",
                            target = "_blank"), ".")
                    )
                ),

                # ── About the Dataset ──────────
                fluidRow(
                    box(
                        title = "About the Dataset",
                        status = "success",
                        width = 12,
                        p("This dashboard displays results from the ",
                          a("Pasilla dataset (GSE18508)",
                            href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18508",
                            target = "_blank"),
                          " — a well-characterized ",
                          em("Drosophila melanogaster"),
                          " RNA-seq benchmark. The experiment compares RNAi knockdown
                          of the splicing factor ",
                          em("pasilla"),
                          " against untreated controls across 4 paired-end samples
                          (2 treated, 2 untreated)."),
                        p("Reference genome: dm6 (BDGP Release 6), Ensembl release 109.")
                    )
                ),

                # ── About the Pipeline ─────────
                fluidRow(
                    box(
                        title = "About the Pipeline",
                        status = "warning",
                        width = 12,
                        p("Raw FASTQ reads are processed through the following steps:"),
                        tags$ol(
                            tags$li("STAR genome index construction"),
                            tags$li("FastQC + MultiQC quality control"),
                            tags$li("Trimmomatic adapter trimming"),
                            tags$li("STAR splice-aware alignment"),
                            tags$li("featureCounts read quantification"),
                            tags$li("DESeq2 differential expression analysis")
                        )
                    )
                ),

                # ── How to Use ─────────────────
                fluidRow(
                    box(
                        title = "How to Use This Dashboard",
                        status = "danger",
                        width = 12,
                        p("Navigate between tabs using the sidebar:"),
                        tags$ul(
                            tags$li(strong("PCA"), " — sample clustering by condition; use to assess batch effects and replicate consistency"),
                            tags$li(strong("Volcano"), " — interactive plot of fold change vs significance; adjust FDR and fold change cutoffs using the sliders"),
                            tags$li(strong("Heatmap"), " — top N most variable genes; adjust the number of genes and scaling using the options panel"),
                            tags$li(strong("Results Table"), " — full DESeq2 results; sortable and filterable by any column")
                        )
                    )
                ),

                # ── Key Results ────────────────
                fluidRow(
                  box(
                    title = "Key Results",
                    status = "primary",
                    width = 12,
                    fluidRow(
                      style = "padding: 0px;",
                      column(4,
                             h4("Top 5 Upregulated Genes"),
                             tableOutput("top_up")
                      ),
                      column(4,
                             h4("Top 5 Downregulated Genes"),
                             tableOutput("top_down")
                      )
                    )
                  )
                )
            ),

            # ── PCA Tab ────────────────────
            tabItem(tabName = "pca",
                fluidRow(
                    box(
                        title = "PCA: Sample Clustering by Condition",
                        width = 12,
                        plotlyOutput("pca_plot", height = "500px")
                    )
                )
            ),

            # ── Volcano Tab ────────────────
            tabItem(tabName = "volcano",
                fluidRow(
                    box(width = 3,
                        title = "Filters",
                        sliderInput("fdr_cutoff",
                                    "FDR Cutoff",
                                    min = 0.001, max = 0.1,
                                    value = 0.05, step = 0.001),
                        p(em("Note: This dataset has a bimodal adjusted p-value distribution 
                              typical of small sample sizes (n=2 per condition). Most genes are 
                              either highly significant (padj < 0.001) or not significant (padj = 1). 
                              The FDR slider may show limited visual change as a result."),
                          style = "font-size: 11px; color: #888;"),
                        sliderInput("lfc_cutoff",
                            "Log2 Fold Change Cutoff",
                            min = 0, max = 4,
                            value = 1, step = 0.1)
                    ),
                    box(width = 9,
                        title = "Volcano Plot",
                        plotlyOutput("volcano_plot", height = "500px")
                    )
                )
            ),

            # ── Heatmap Tab ────────────────
            tabItem(tabName = "heatmap",
                fluidRow(
                    box(width = 3,
                        title = "Options",
                        sliderInput("n_genes",
                            "Number of Top Variable Genes",
                            min = 10, max = 100,
                            value = 50, step = 10),
                        selectInput("scale_by",
                            "Scale by",
                            choices = c("row", "column", "none"),
                            selected = "row")
                    ),
                    box(width = 9,
                        title = "Heatmap of Top Variable Genes",
                        plotOutput("heatmap_plot", height = "800px")
                    )
                )
            ),

            # ── Results Table Tab ──────────
            tabItem(tabName = "table",
                fluidRow(
                    box(width = 12,
                        title = "DESeq2 Results",
                        DTOutput("results_table")
                    )
                )
            )
        )
    )
)

# ── Server ─────────────────────────────────
server <- function(input, output, session) {

    # Reactive: filter results based on slider inputs
    filtered_results <- reactive({
        deseq2_results %>%
            filter(!is.na(padj),
                   padj < input$fdr_cutoff,
                   abs(log2FoldChange) > input$lfc_cutoff)
    })

    # ── Overview: Top genes tables ─────────
    output$top_up <- renderTable({
      deseq2_results %>%
        filter(!is.na(padj), padj < 0.05, log2FoldChange > 0) %>%
        arrange(padj, desc(log2FoldChange)) %>%
        slice_head(n = 5) %>%
        dplyr::select(label, log2FoldChange, padj) %>%
        rename("Gene"         = label,
               "Log2FC"       = log2FoldChange,
               "Adj. P-value" = padj) %>%
        mutate(
          Log2FC = round(Log2FC, 4),
          `Adj. P-value` = formatC(`Adj. P-value`, format = "e", digits = 2)
        )
    })
    
    output$top_down <- renderTable({
      deseq2_results %>%
        filter(!is.na(padj), padj < 0.05, log2FoldChange < 0) %>%
        arrange(padj, log2FoldChange) %>%
        slice_head(n = 5) %>%
        dplyr::select(label, log2FoldChange, padj) %>%
        rename("Gene"         = label,
               "Log2FC"       = log2FoldChange,
               "Adj. P-value" = padj) %>%
        mutate(
          Log2FC = round(Log2FC, 4),
          `Adj. P-value` = formatC(`Adj. P-value`, format = "e", digits = 2)
        )
    })

    # ── PCA plot ───────────────────────────
    output$pca_plot <- renderPlotly({
        pca_mat <- normalized_counts %>%
            dplyr::select(starts_with("GSM")) %>%
            t()

        pca_result <- prcomp(pca_mat, scale. = TRUE)
        pca_data <- as.data.frame(pca_result$x[, 1:2])
        pca_data$condition <- sample_info$condition
        pca_data$sample <- rownames(pca_data)
        pca_var <- round(100 * summary(pca_result)$importance[2, 1:2], 1)

        plot_ly(pca_data,
                x = ~PC1,
                y = ~PC2,
                color = ~condition,
                colors = c("treated" = "#F8766D", "untreated" = "#00BFC4"),
                text  = ~sample,
                type  = "scatter",
                mode  = "markers+text",
                textposition = "top center",
                marker = list(size = 10)
        ) %>%
          layout(
            xaxis = list(
              title = paste0("PC1: ", pca_var[1], "% variance"),
              zeroline = FALSE
            ),
            yaxis = list(
              title = paste0("PC2: ", pca_var[2], "% variance"),
              zeroline = FALSE
            ),
            title = "PCA of VST-normalized counts",
            legend = list(title = list(text = "Condition"))
          )
    })

    # ── Volcano plot ───────────────────────
    output$volcano_plot <- renderPlotly({
      df <- deseq2_results %>%
        filter(!is.na(padj), !is.na(log2FoldChange)) %>%
        mutate(
          significance = case_when(
            padj < input$fdr_cutoff & log2FoldChange >  input$lfc_cutoff ~ "Up",
            padj < input$fdr_cutoff & log2FoldChange < -input$lfc_cutoff ~ "Down",
            TRUE ~ "NS"
          ),
          neg_log10_padj = -log10(padj)
        )
      
      # Top genes to label
      top_genes <- df %>%
        filter(padj < input$fdr_cutoff, abs(log2FoldChange) > 2) %>%
        arrange(padj) %>%
        slice_head(n = 5)
      
      # Always include ps
      ps_gene <- df %>% filter(label == "ps")
      top_genes <- bind_rows(top_genes, ps_gene) %>% distinct()
      
      # Base plot
      p <- plot_ly(df,
                   x = ~log2FoldChange,
                   y = ~neg_log10_padj,
                   color = ~significance,
                   colors = c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey70"),
                   text  = ~label,
                   hoverinfo = "text",
                   type  = "scatter",
                   mode  = "markers",
                   marker = list(size = 8, opacity = 0.7)
      )
      
      # Add labels for top genes
      p <- p %>% add_annotations(
        x         = top_genes$log2FoldChange,
        y         = top_genes$neg_log10_padj,
        text      = top_genes$label,
        showarrow = FALSE,
        font      = list(size = 12),
        yshift    = 10
      )
      
      p %>% layout(
        xaxis = list(title = "Log2 Fold Change", zeroline = FALSE),
        yaxis = list(title = "-log10(adjusted p-value)", zeroline = FALSE),
        shapes = list(
          list(type = "line",
               x0 = -input$lfc_cutoff, x1 = -input$lfc_cutoff,
               y0 = 0, y1 = max(-log10(df$padj), na.rm = TRUE),
               line = list(dash = "dot", color = "black")),
          list(type = "line",
               x0 = input$lfc_cutoff, x1 = input$lfc_cutoff,
               y0 = 0, y1 = max(-log10(df$padj), na.rm = TRUE),
               line = list(dash = "dot", color = "black"))
        )
      )
    })

    # ── Heatmap ────────────────────────────
    output$heatmap_plot <- renderPlot({
      mat <- normalized_counts %>%
        dplyr::select(starts_with("GSM")) %>%
        mutate(variance = apply(across(everything()), 1, var)) %>%
        arrange(desc(variance)) %>%
        slice_head(n = input$n_genes) %>%
        dplyr::select(-variance) %>%
        as.matrix()
      
      rownames(mat) <- mapIds(
        org.Dm.eg.db,
        keys      = rownames(normalized_counts)[1:input$n_genes],
        column    = "SYMBOL",
        keytype   = "FLYBASE",
        multiVals = "first"
      ) %>% ifelse(is.na(.), rownames(normalized_counts)[1:input$n_genes], .)
      
      pheatmap(mat,
               annotation_col = sample_info,
               scale          = input$scale_by,
               show_rownames  = TRUE,
               fontsize_row   = 11,
               angle_col      = 45
      )
    })

    # ── Results table ──────────────────────
    output$results_table <- renderDT({
      deseq2_results %>%
        filter(!is.na(padj)) %>%
        dplyr::select(gene_id, label, log2FoldChange, padj, pvalue,
                      baseMean, lfcSE) %>%
        rename("Gene ID"       = gene_id,
               "Gene Symbol"   = label,
               "Log2FC"        = log2FoldChange,
               "Adj. P-value"  = padj,
               "P-value"       = pvalue,
               "Base Mean"     = baseMean,
               "LFC Std Error" = lfcSE) %>%
        mutate(
          Log2FC          = round(Log2FC, 4),
          `Base Mean`     = round(`Base Mean`, 2),
          `LFC Std Error` = round(`LFC Std Error`, 4),
          abs_lfc         = abs(Log2FC)
        ) %>%
        datatable(
          filter    = "top",
          options   = list(
            pageLength = 25,
            scrollX = TRUE,
            order = list(list(3, "asc"), list(7, "desc")),
            columnDefs = list(
              list(visible = FALSE, targets = 7),
              list(
                targets = c(3, 4),
                render = JS(
                  "function(data, type, row) {
                                if (type === 'display') {
                                    return parseFloat(data).toExponential(2);
                                }
                                return data;
                            }"
                )
              )
            )
          ),
          rownames = FALSE
        ) %>%
        formatStyle("Log2FC",
                    backgroundColor = styleInterval(
                      c(-1, 1),
                      c("#FFAAAA", "white", "#AEC6E8")
                    )
        )
    })
}

shinyApp(ui, server)