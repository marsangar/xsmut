library(shiny)
library(xsmut)
library(data.table)

options(shiny.maxRequestSize = 2000 * 1024^2)  # COSMIC files are large

ui <- fluidPage(
  titlePanel("xsmut — species mutations vs human COSMIC"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("1. Genes"),
      textInput("genes", "Gene symbol(s), comma-separated", "TP53"),
      h4("2. Species"),
      textInput("species", "Species (any name Ensembl knows)", "macaque"),
      textInput("assembly", "Assembly (checked against Ensembl)", "Mmul_10"),
      selectInput("human_asm", "Human assembly (must match COSMIC file)", c("GRCh38", "GRCh37")),
      h4("3. COSMIC files (local)"),
      helpText("Download from cancer.sanger.ac.uk/cosmic/download (login required)."),
      fileInput("cosmic", "Mutation export (.tsv / .tsv.gz)"),
      fileInput("census", "Cancer Gene Census (optional)"),
      h4("4. Species mutations"),
      fileInput("muts", "CSV/TSV/RDS: xsmut columns (gene, aa_change or aa_pos, …) or dNdScv annotmuts"),
      h4("5. Display"),
      radioButtons("species_panel", "Show my mutations in the",
                   c("Top panel" = "top", "Bottom panel" = "bottom"), inline = TRUE),
      checkboxInput("compress", "Compress introns", TRUE),
      numericInput("min_label", "Label COSMIC hotspots with ≥ N samples", 5, min = 1),
      checkboxInput("highlight", "Highlight matches to top COSMIC hotspots", FALSE),
      conditionalPanel(
        "input.highlight",
        numericInput("hs_top", "Top N COSMIC hotspot residues", 10, min = 1),
        numericInput("hs_window", "Close match window (± aa; 0 = same residue)", 3, min = 0),
        numericInput("hs_min", "Min. COSMIC samples per hotspot", 2, min = 1),
        radioButtons("hs_col", "Halo colour", c("Gold" = "gold", "Red (slides)" = "red"), inline = TRUE),
        radioButtons("hs_style", "Highlight", c("Strong" = "strong", "Subtle" = "subtle"), inline = TRUE)
      ),
      selectInput("legend_pos", "Legend", c("right", "bottom", "none")),
      checkboxInput("pictos", "Species pictograms (PhyloPic; needs rphylopic + internet)", TRUE),
      actionButton("go", "Run", class = "btn-primary"),
      br(), br(),
      downloadButton("dl_plot", "Download plot (PDF)"),
      downloadButton("dl_table", "Download mapped mutations (CSV)"),
      downloadButton("dl_hits", "Download hotspot matches (CSV)")
    ),
    mainPanel(
      width = 9,
      uiOutput("gene_tabs")
    )
  )
)

server <- function(input, output, session) {

  cosmic_dt <- reactive({
    req(input$cosmic)
    genes <- trimws(strsplit(input$genes, ",")[[1]])
    withProgress(message = "Reading COSMIC…", cosmic_load(input$cosmic$datapath, genes = genes))
  })

  user_muts <- reactive({
    req(input$muts)
    f <- input$muts$datapath
    x <- if (grepl("\\.rds$", input$muts$name, ignore.case = TRUE)) readRDS(f) else fread(f)
    if (is_dndscv(x)) as_xsmut(x) else x
  })

  results <- eventReactive(input$go, {
    genes <- trimws(strsplit(input$genes, ",")[[1]])
    cos <- cosmic_dt(); mu <- user_muts()
    census <- if (!is.null(input$census)) input$census$datapath else NULL
    withProgress(message = "Querying Ensembl & mapping…", value = 0, {
      out <- list()
      for (g in genes) {
        incProgress(1 / length(genes), detail = g)
        out[[g]] <- tryCatch(
          xsmut_pipeline(g, mu, input$species, input$assembly, cos, census,
                         human_assembly = input$human_asm,
                         species_panel = input$species_panel,
                         highlight_hotspots = isTRUE(input$highlight),
                         hotspot_top_n = input$hs_top, hotspot_window = input$hs_window,
                         hotspot_min_samples = input$hs_min,
                         highlight_colour = input$hs_col, legend_position = input$legend_pos,
                         show_pictograms = isTRUE(input$pictos),
                         highlight_style = if (is.null(input$hs_style)) "strong" else input$hs_style,
                         compress_introns = input$compress, min_label_n = input$min_label),
          error = function(e) list(error = conditionMessage(e)))
      }
      out
    })
  })

  output$gene_tabs <- renderUI({
    res <- results()
    tabs <- lapply(names(res), function(g) {
      r <- res[[g]]
      if (!is.null(r$error)) return(tabPanel(g, div(class = "alert alert-danger", r$error)))
      chk <- r$gene_check
      orth <- attr(r$species_mutations, "ortholog")$target
      n_hits <- if (is.null(r$hotspot_matches)) NA_integer_ else nrow(r$hotspot_matches)
      tabPanel(g,
        div(class = if (isTRUE(chk$in_census)) "alert alert-success" else "alert alert-warning",
            sprintf("%s — %s. COSMIC mutations in file: %d. Ortholog in %s: %s (%s%% identity, %s).",
                    g,
                    if (isTRUE(chk$in_census)) paste0("in Cancer Gene Census", if (!is.na(chk$tier)) paste0(" (tier ", chk$tier, ")") else "")
                    else if (is.na(chk$in_census)) "census file not supplied" else "NOT in Cancer Gene Census",
                    chk$n_cosmic_mutations, r$species_name, orth$protein_id,
                    round(as.numeric(orth$perc_id), 1), orth$type)),
        if (!is.na(n_hits))
          div(class = if (n_hits > 0) "alert alert-info" else "alert alert-warning",
              sprintf("%d of your mutations match (±%d aa) the top %d COSMIC hotspot residues.",
                      n_hits, input$hs_window, input$hs_top)),
        plotOutput(paste0("plot_", g), height = "720px"),
        if (length(r$pictogram_credits))
          helpText("Silhouettes from PhyloPic: ", paste(r$pictogram_credits, collapse = "; ")),
        if (!is.na(n_hits)) tagList(h4("Matches to COSMIC hotspots"),
                                    DT::dataTableOutput(paste0("tbl_hit_", g))),
        h4("Species mutations mapped to human"),
        DT::dataTableOutput(paste0("tbl_sp_", g)),
        h4("Human COSMIC mutations"),
        DT::dataTableOutput(paste0("tbl_cos_", g))
      )
    })
    do.call(tabsetPanel, c(tabs, id = "tabs"))
  })

  observe({
    res <- results()
    for (g in names(res)) local({
      gg <- g; r <- res[[gg]]
      if (is.null(r$error)) {
        output[[paste0("plot_", gg)]] <- renderPlot(print(r$plot))
        if (!is.null(r$hotspot_matches))
          output[[paste0("tbl_hit_", gg)]] <- DT::renderDataTable(r$hotspot_matches,
                                                                  options = list(pageLength = 10, scrollX = TRUE))
        output[[paste0("tbl_sp_", gg)]] <- DT::renderDataTable(r$species_mutations, options = list(pageLength = 10, scrollX = TRUE))
        output[[paste0("tbl_cos_", gg)]] <- DT::renderDataTable(
          r$cosmic_mutations[, .(aa_change, cds_change, consequence, chrom, start, primary_site, sample_id)],
          options = list(pageLength = 10, scrollX = TRUE))
      }
    })
  })

  output$dl_plot <- downloadHandler(
    filename = function() paste0("xsmut_", gsub("[^A-Za-z0-9]", "_", input$genes), ".pdf"),
    content = function(file) {
      res <- results(); ok <- Filter(function(r) is.null(r$error), res)
      if (capabilities("cairo")) grDevices::cairo_pdf(file, width = 12, height = 8, onefile = TRUE) else pdf(file, width = 12, height = 8)
      for (r in ok) print(r$plot)
      dev.off()
    })

  output$dl_table <- downloadHandler(
    filename = function() "xsmut_mapped_mutations.csv",
    content = function(file) {
      res <- results(); ok <- Filter(function(r) is.null(r$error), res)
      fwrite(rbindlist(lapply(ok, `[[`, "species_mutations"), fill = TRUE), file)
    })

  output$dl_hits <- downloadHandler(
    filename = function() "xsmut_hotspot_matches.csv",
    content = function(file) {
      res <- results(); ok <- Filter(function(r) is.null(r$error) && !is.null(r$hotspot_matches), res)
      hits <- rbindlist(lapply(names(ok), function(g) {
        m <- ok[[g]]$hotspot_matches
        if (nrow(m)) cbind(gene = g, m) else NULL
      }), fill = TRUE)
      fwrite(hits, file)
    })
}

shinyApp(ui, server)
