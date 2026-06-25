# ─────────────────────────────────────────
# app.R
# Shiny dashboard for RNA-seq pipeline results
# Tabs: Overview | Volcano | Heatmap | Table
# ─────────────────────────────────────────

source("global.R")

# ── UI ─────────────────────────────────────
ui <- dashboardPage(
    dashboardHeader(title = "RNA-seq Results"),

    dashboardSidebar(
        sidebarMenu(
            menuItem("Overview",  tabName = "overview",  icon = icon("chart-bar")),
            menuItem("Volcano",   tabName = "volcano",   icon = icon("circle-dot")),
            menuItem("Heatmap",   tabName = "heatmap",   icon = icon("th")),
            menuItem("Results Table", tabName = "table", icon = icon("table"))
        )
    ),

    dashboardBody(
        tabItems(

            # ── Overview Tab ───────────────
            tabItem(tabName = "overview",
                fluidRow(
                    valueBox(nrow(deseq2_results), "Genes Tested",
                        icon = icon("dna"), color = "blue"),
                    valueBox(nrow(sig_genes), "Significant Genes (FDR < 0.05)",
                        icon = icon("star"), color = "green"),
                    valueBox(
                        sum(sig_genes$log2FoldChange > 0),
                        "Upregulated",
                        icon = icon("arrow-up"), color = "red"),
                    valueBox(
                        sum(sig_genes$log2FoldChange < 0),
                        "Downregulated",
                        icon = icon("arrow-down"), color = "purple")
                ),
                fluidRow(
                    box(
                        title = "PCA: Sample Clustering by Condition",
                        width = 12,
                        plotlyOutput("pca_plot")
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
                        plotOutput("heatmap_plot", height = "600px")
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

    # PCA plot
    output$pca_plot <- renderPlotly({
        # Placeholder until real VST data is wired in
        plot_ly(type = "scatter", mode = "markers") %>%
            layout(title = "PCA plot will render with real data")
    })

    # Volcano plot
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

        plot_ly(df,
            x = ~log2FoldChange,
            y = ~neg_log10_padj,
            color = ~significance,
            colors = c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey70"),
            text  = ~gene_id,
            type  = "scatter",
            mode  = "markers",
            marker = list(size = 5, opacity = 0.7)
        ) %>%
        layout(
            xaxis = list(title = "Log2 Fold Change"),
            yaxis = list(title = "-log10(adjusted p-value)"),
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

    # Heatmap
    output$heatmap_plot <- renderPlot({
        mat <- normalized_counts %>%
            mutate(variance = apply(across(everything()), 1, var)) %>%
            arrange(desc(variance)) %>%
            slice_head(n = input$n_genes) %>%
            select(-variance) %>%
            as.matrix()

        pheatmap(mat,
            scale        = input$scale_by,
            show_rownames = input$n_genes <= 50,
            fontsize_row  = 8
        )
    })

    # Results table
    output$results_table <- renderDT({
        deseq2_results %>%
            filter(!is.na(padj)) %>%
            mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
            datatable(
                filter    = "top",
                options   = list(pageLength = 25, scrollX = TRUE),
                rownames  = FALSE
            ) %>%
            formatStyle("log2FoldChange",
                backgroundColor = styleInterval(
                    c(-1, 1),
                    c("#AEC6E8", "white", "#FFAAAA")
                )
            )
    })
}

shinyApp(ui, server)