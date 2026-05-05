library(shiny)
library(pROC)
library(caTools)
library(matrixTests)
library(ggplot2)
library(shinydashboard)
library(stringr)
library(reactable)
library(shinycustomloader)
library(cowplot)
library(shinyWidgets)
library(rlist)
library(funprog)
library(data.table)
library(Seurat)
library(RColorBrewer)
library(tools)
library(plotly)

options(shiny.maxRequestSize = 1024^3)

reactMeta <- read.csv("../../../data/fingerprint_prep_objects/human_reaction_meta.csv", stringsAsFactors = FALSE)
colnames(reactMeta)[1] <- "metabolite"

sc_object_files <- list.files("../../../data/sc_objects", pattern = "\\.rds$", full.names = TRUE)
if (length(sc_object_files) == 0) {
  stop("No .rds files found in sc_objects/")
}
sc_object_names <- file_path_sans_ext(basename(sc_object_files))
names(sc_object_files) <- sc_object_names

find_flux_file <- function(assay_name) {
  file.path("../../../data/sc_objects", paste0(assay_name, "_flux.csv"))
}

safe_read_flux <- function(assay_name) {
  flux_file <- find_flux_file(assay_name)
  if (!file.exists(flux_file)) {
    stop(sprintf("Expected flux file not found: %s", flux_file))
  }
  as.data.frame(as.matrix(fread(flux_file), rownames = 1))
}

get_choices_from_meta <- function(obj) {
  md <- obj@meta.data
  sample_choices <- if ("sample" %in% colnames(md)) sort(unique(as.character(md$sample))) else character(0)
  treatment_choices <- if ("treatment" %in% colnames(md)) sort(unique(as.character(md$treatment))) else character(0)
  cluster_choices <- sort(unique(as.character(Idents(obj))))
  list(
    sample = sample_choices,
    treatment = treatment_choices,
    cluster = cluster_choices
  )
}

subset_assay_object <- function(obj, samples = NULL, treatments = NULL, clusters = NULL) {
  if (is.null(obj)) return(NULL)
  md <- obj@meta.data
  keep <- rep(TRUE, nrow(md))
  names(keep) <- rownames(md)
  
  if (!is.null(samples) && length(samples) > 0 && "sample" %in% colnames(md)) {
    keep <- keep & md$sample %in% samples
  }
  if (!is.null(treatments) && length(treatments) > 0 && "treatment" %in% colnames(md)) {
    keep <- keep & md$treatment %in% treatments
  }
  
  kept_cells <- rownames(md)[keep]
  if (length(kept_cells) == 0) return(NULL)
  out <- subset(obj, cells = kept_cells)
  
  if (!is.null(clusters) && length(clusters) > 0) {
    if ("seurat_clusters" %in% colnames(out@meta.data)) {
      out <- subset(out, subset = seurat_clusters %in% clusters)
    } else {
      ident_keep <- as.character(Idents(out)) %in% as.character(clusters)
      kept_cells_2 <- colnames(out)[ident_keep]
      if (length(kept_cells_2) == 0) return(NULL)
      out <- subset(out, cells = kept_cells_2)
    }
  }
  
  if (ncol(out) == 0) return(NULL)
  out
}

build_assay_meta <- function(cell_names, label) {
  out <- data.frame(selected_assay = rep(label, length(cell_names)), stringsAsFactors = FALSE)
  rownames(out) <- cell_names
  out
}

run_marker_analysis <- function(flux_1, flux_2, names_1, names_2, auc_cutoff, fc_cutoff, fdr_cutoff, reactMeta) {
  validate(
    need(!is.null(flux_1), "Flux file for assay 1 could not be found."),
    need(!is.null(flux_2), "Flux file for assay 2 could not be found."),
    need(length(names_1) > 0, "Assay 1 has no selected cells."),
    need(length(names_2) > 0, "Assay 2 has no selected cells.")
  )
  
  flux <- rbind(flux_1, flux_2)
  flux <- flux[unique(rownames(flux)), , drop = FALSE]
  
  assay_meta_1 <- build_assay_meta(names_1, "assay_1")
  assay_meta_2 <- build_assay_meta(names_2, "assay_2")
  metaFinal <- rbind(assay_meta_1, assay_meta_2)
  
  common_cells <- intersect(rownames(metaFinal), rownames(flux))
  validate(need(length(common_cells) > 1, "Not enough overlapping cells between selected assays and flux matrices."))
  
  metaFinal <- metaFinal[common_cells, , drop = FALSE]
  flux_full <- cbind(metaFinal, flux[common_cells, , drop = FALSE])
  flux_full <- flux_full[complete.cases(flux_full[, 2, drop = TRUE]), , drop = FALSE]
  metaFinal <- metaFinal[rownames(flux_full), , drop = FALSE]
  
  data_filt_var <- flux_full
  outer <- "selected_assay"
  
  validate(need(nrow(subset(data_filt_var, selected_assay == "assay_1")) > 1, "Assay 1 needs at least two rows after filtering."))
  validate(need(nrow(subset(data_filt_var, selected_assay == "assay_2")) > 1, "Assay 2 needs at least two rows after filtering."))
  
  AUC <- colAUC(data_filt_var[, -1, drop = FALSE], factor(data_filt_var[, "selected_assay"]), plotROC = FALSE)
  AUC <- AUC[1, ]
  
  mean_1 <- apply(subset(data_filt_var, selected_assay == "assay_1")[, -1, drop = FALSE], 2, mean)
  mean_2 <- apply(subset(data_filt_var, selected_assay == "assay_2")[, -1, drop = FALSE], 2, mean)
  
  logFC <- log2(mean_1 / mean_2)
  logFC[!is.finite(logFC)] <- 0
  logFC[AUC < auc_cutoff] <- 0
  
  short <- as.data.frame(logFC[which(abs(logFC) > fc_cutoff)])
  validate(need(nrow(short) > 0, "No metabolites passed the AUC and Log2FC thresholds."))
  
  t <- col_t_welch(
    as.matrix(subset(data_filt_var, selected_assay == "assay_1")[, rownames(short), drop = FALSE]),
    as.matrix(subset(data_filt_var, selected_assay == "assay_2")[, rownames(short), drop = FALSE])
  )
  
  final <- cbind(t$mean.y, t$mean.x, short, t$statistic, t$pvalue)
  colnames(final) <- c("var_1_mean", "var_2_mean", "log2FC", "t_welch_statistic", "pval")
  final <- subset(final, pval < 0.05)
  validate(need(nrow(final) > 0, "No metabolites passed the nominal p-value threshold."))
  
  final$fdr <- p.adjust(final[, "pval"], method = "fdr")
  final <- subset(final, fdr < fdr_cutoff)
  validate(need(nrow(final) > 0, "No metabolites passed the FDR threshold."))
  
  metab_names <- unique(rownames(final))
  names_keep <- c(outer, metab_names)
  data_filt_var <- data_filt_var[, names_keep, drop = FALSE]
  
  t2 <- col_t_welch(
    as.matrix(subset(data_filt_var, data_filt_var[, outer] == "assay_1")[, metab_names, drop = FALSE]),
    as.matrix(subset(data_filt_var, data_filt_var[, outer] == "assay_2")[, metab_names, drop = FALSE])
  )
  
  final2 <- as.data.frame(cbind(metabolite = metab_names, as.numeric(t2$mean.y), as.numeric(t2$mean.x), as.numeric(t2$statistic), as.numeric(t2$pvalue)), stringsAsFactors = FALSE)
  colnames(final2) <- c("metabolite", "var_1_mean", "var_2_mean", "t_welch_statistic", "pval")
  final2$sign_t_log_pval <- ifelse(final2$t_welch_statistic > 0, -log10(as.numeric(final2$pval)), log10(as.numeric(final2$pval)))
  
  metabs_tab <- merge(final2, reactMeta, by = "metabolite", all.x = TRUE, all.y = FALSE)
  metabs_tab$sig <- ifelse(metabs_tab$sign_t_log_pval > 1.301, "up", ifelse(metabs_tab$sign_t_log_pval < -1.301, "down", "ns"))
  
  list(
    flux = flux,
    metaFinal = metaFinal,
    data_filt_var = data_filt_var,
    metabs_tab = metabs_tab,
    outer = outer,
    metab_names = metab_names,
    AUC = AUC
  )
}

header <- dashboardHeader(
  title = "RADAR | SCepter",
  titleWidth = 400
)

body <- dashboardBody(
  fluidRow(
    column(
      width = 12,
      box(
        width = 6,
        title = "Select assay 1",
        solidHeader = TRUE,
        status = "primary",
        selectizeInput("assay_file_1", "Assay", choices = sc_object_names, selected = sc_object_names[1])
      ),
      box(
        width = 6,
        title = "Select assay 2",
        solidHeader = TRUE,
        status = "primary",
        selectizeInput("assay_file_2", "Assay", choices = sc_object_names, selected = sc_object_names[min(2, length(sc_object_names))])
      )
    )
  ),
  fluidRow(
    column(
      width = 12,
      box(
        width = 6,
        title = textOutput("plot1Title"),
        solidHeader = TRUE,
        withLoader(plotOutput("plot1", click = "clickBar", height = 500), type = "html", loader = "dnaspin"),
        box(width = 12, noUiSliderInput("resolution_1", "Adjust Resolution", min = 0.1, max = 1.0, value = 0.2, tooltip = TRUE, step = 0.1)),
        box(width = 6, virtualSelectInput("sample_1", "Select Sample(s)", choices = character(0), showValueAsTags = TRUE, search = TRUE, multiple = TRUE, disableSelectAll = FALSE)),
        box(width = 6, virtualSelectInput("treatment_1", "Select Treatment(s)", choices = character(0), showValueAsTags = TRUE, search = TRUE, multiple = TRUE, disableSelectAll = FALSE)),
        box(width = 12, virtualSelectInput("cluster_1", "Select Cluster(s)", choices = character(0), showValueAsTags = TRUE, search = TRUE, multiple = TRUE, disableSelectAll = FALSE))
      ),
      box(
        width = 6,
        title = textOutput("plot2Title"),
        solidHeader = TRUE,
        withLoader(plotOutput("plot2", click = "clickBar", height = 500), type = "html", loader = "dnaspin"),
        box(width = 12, noUiSliderInput("resolution_2", "Adjust Resolution", min = 0.1, max = 1.0, value = 0.2, tooltip = TRUE, step = 0.1)),
        box(width = 6, virtualSelectInput("sample_2", "Select Sample(s)", choices = character(0), showValueAsTags = TRUE, search = TRUE, multiple = TRUE, disableSelectAll = FALSE)),
        box(width = 6, virtualSelectInput("treatment_2", "Select Treatment(s)", choices = character(0), showValueAsTags = TRUE, search = TRUE, multiple = TRUE, disableSelectAll = FALSE)),
        box(width = 12, virtualSelectInput("cluster_2", "Select Cluster(s)", choices = character(0), showValueAsTags = TRUE, search = TRUE, multiple = TRUE, disableSelectAll = FALSE))
      )
    )
  ),
  fluidRow(
    column(
      width = 12,
      box(
        width = 12,
        title = "RADAR | SCepter Pre-Analysis Summary",
        solidHeader = TRUE,
        status = "primary",
        box(width = 6, verbatimTextOutput("reportText")),
        box(
          width = 6,
          noUiSliderInput("customAUC", "AUC:", min = 0.6, max = 1.0, value = 0.65, tooltip = TRUE, step = 0.01),
          noUiSliderInput("customFC", "Log2FC:", min = 0, max = 5, value = 1, tooltip = TRUE, step = 0.01),
          noUiSliderInput("customFDR", "FDR:", min = 0, max = 1, value = 0.05, tooltip = TRUE, step = 0.01)
        ),
        fluidRow(column(width = 12, actionButton("submitReport", "Begin Analysis"))),
        column(width = 12, actionButton("submitReportDefault", "Begin Analysis (default parameters)"))
      )
    )
  ),
  fluidRow(
    column(
      width = 12,
      box(
        title = textOutput("plot1Title_metab"),
        width = 12, solidHeader = TRUE, status = "primary",
        plotlyOutput("plot_metab", height = 650),
        tags$div(
          style = "margin-top: 10px; padding: 10px 12px; border-radius: 6px;
               background-color: #e8ebff; border: 1px solid #d1d7ff;
               font-size: 13px; line-height: 1.4; color: #273043;",
          tags$ul(
            style = "padding-left: 18px; margin: 0;",
            tags$li("Use zoom/pan to inspect dense regions."),
            tags$li("Use box or lasso selection from the Plotly toolbar to select reactions."),
            tags$li("Click a point to make it the active reaction for detailed plots.")
          )
        ),
        downloadLink("dlPlot1", "Download Plot as PDF"),
        downloadLink("fingerprint", "Download Fingerprint Data")
      )
    )
  ),
  fluidRow(column(width = 8, box(title = textOutput("outerTitle"), solidHeader = TRUE, status = "primary", checkboxGroupInput("outerSelect", label = NULL, choices = c("Significant (+)" = "up", "Significant (-)" = "down"), selected = c("up", "down"))))),
  fluidRow(column(width = 12, box(width = 12, reactableOutput("info")))),
  fluidRow(column(width = 12, box(width = 12, withLoader(plotOutput("plot3_metab", height = 500), type = "html", loader = "dnaspin"), actionLink("invert", "Invert Boxplot | "), downloadLink("dlPlot3", "Download Plot as PDF")))),
  fluidRow(
    column(
      width = 12,
      box(
        width = 12,
        title = textOutput("plotFeatTitle"),
        solidHeader = TRUE,
        status = "primary",
        withLoader(plotOutput("plot_feat_combined", height = 550), type = "html", loader = "dnaspin")
      )
    )
  )
  )


ui <- dashboardPage(
  skin = "purple",
  header = header,
  sidebar = dashboardSidebar(disable = TRUE),
  body = body
)

server <- function(input, output, session) {
  selected_reaction <- reactiveVal(NULL)
  
  assay_1_full <- reactiveVal(NULL)
  assay_1 <- reactiveVal(NULL)
  assay_1_mid <- reactiveVal(NULL)
  assay_1_active <- reactiveVal(NULL)
  flux_1 <- reactiveVal(NULL)
  names_1 <- reactiveVal(character(0))
  
  assay_2_full <- reactiveVal(NULL)
  assay_2 <- reactiveVal(NULL)
  assay_2_mid <- reactiveVal(NULL)
  assay_2_active <- reactiveVal(NULL)
  flux_2 <- reactiveVal(NULL)
  names_2 <- reactiveVal(character(0))
  
  plot_max <- reactiveVal(NULL)
  plot_min <- reactiveVal(NULL)
  errors_global <- reactiveVal("NONE")
  analysis_results <- reactiveVal(NULL)
  selected_subsystem <- reactiveVal(NULL)
  invert_boxplot <- reactiveVal(FALSE)
  
  output$plot1Title <- renderText({
    paste0("DimPlot - ", input$assay_file_1)
  })
  
  output$plot2Title <- renderText({
    paste0("DimPlot - ", input$assay_file_2)
  })
  
  output$plotFeatTitle <- renderText({
    rxn <- selected_reaction()
    if (is.null(rxn)) {
      paste0("FeaturePlot - ", input$assay_file_1 %||% "Assay 1", " vs ", input$assay_file_2 %||% "Assay 2")
    } else {
      paste0("FeaturePlot - ", rxn, " (", input$assay_file_1, " vs ", input$assay_file_2, ")")
    }
  })
  
  load_assay_panel <- function(assay_name, resolution,
                               assay_full_val,   # new argument
                               assay_val, assay_mid_val, assay_active_val,
                               flux_val,
                               sample_id, treatment_id, cluster_id,
                               names_val) {
    req(assay_name)
    obj <- readRDS(sc_object_files[[assay_name]])
    obj <- FindClusters(obj, resolution = resolution)
    
    assay_full_val(obj)           # store full clustered object
    assay_val(obj)
    assay_mid_val(obj)
    assay_active_val(obj)
    flux_val(safe_read_flux(assay_name))
    names_val(rownames(obj@meta.data))
    
    ch <- get_choices_from_meta(obj)   # choices from FULL object
    updateVirtualSelect(sample_id,   choices = ch$sample,   selected = ch$sample)
    updateVirtualSelect(treatment_id, choices = ch$treatment, selected = ch$treatment)
    updateVirtualSelect(cluster_id,   choices = ch$cluster,   selected = ch$cluster)
  }
  
  observeEvent(list(input$assay_file_1, input$resolution_1), {
    load_assay_panel(
      assay_name      = input$assay_file_1,
      resolution      = input$resolution_1,
      assay_full_val  = assay_1_full,        # new
      assay_val       = assay_1,
      assay_mid_val   = assay_1_mid,
      assay_active_val= assay_1_active,
      flux_val        = flux_1,
      sample_id       = "sample_1",
      treatment_id    = "treatment_1",
      cluster_id      = "cluster_1",
      names_val       = names_1
    )
  }, ignoreNULL = FALSE)
  
  observeEvent(list(input$assay_file_2, input$resolution_2), {
    load_assay_panel(
      assay_name      = input$assay_file_2,
      resolution      = input$resolution_2,
      assay_full_val  = assay_2_full,        # new
      assay_val       = assay_2,
      assay_mid_val   = assay_2_mid,
      assay_active_val= assay_2_active,
      flux_val        = flux_2,
      sample_id       = "sample_2",
      treatment_id    = "treatment_2",
      cluster_id      = "cluster_2",
      names_val       = names_2
    )
  }, ignoreNULL = FALSE)
  
  observeEvent(list(input$sample_1, input$treatment_1, input$cluster_1, assay_1_full()), {
    req(assay_1_full())
    obj_sub <- subset_assay_object(assay_1_full(), input$sample_1, input$treatment_1, input$cluster_1)
    assay_1_mid(obj_sub)
    assay_1_active(obj_sub)
    if (!is.null(obj_sub)) {
      # DO NOT change cluster_1 choices here
      names_1(rownames(obj_sub@meta.data))
    } else {
      names_1(character(0))
    }
  }, ignoreInit = TRUE)
  
  observeEvent(list(input$sample_2, input$treatment_2, input$cluster_2, assay_2_full()), {
    req(assay_2_full())
    obj_sub <- subset_assay_object(assay_2_full(), input$sample_2, input$treatment_2, input$cluster_2)
    assay_2_mid(obj_sub)
    assay_2_active(obj_sub)
    if (!is.null(obj_sub)) {
      # do not touch cluster_2 choices
      names_2(rownames(obj_sub@meta.data))
    } else {
      names_2(character(0))
    }
  }, ignoreInit = TRUE)
  
  output$plot1 <- renderPlot({
    req(assay_1_active())
    DimPlot(assay_1_active())
  })
  
  output$plot2 <- renderPlot({
    req(assay_2_active())
    DimPlot(assay_2_active())
  })
  
  observeEvent(list(input$assay_file_1, input$assay_file_2, input$resolution_1, input$resolution_2, input$sample_1, input$sample_2, input$treatment_1, input$treatment_2, input$cluster_1, input$cluster_2, input$customAUC, input$customFC, input$customFDR), {
    errors <- c()
    if (is.null(assay_1_active()) || length(names_1()) == 0) errors <- c(errors, "Assay 1 selection returned no cells")
    if (is.null(assay_2_active()) || length(names_2()) == 0) errors <- c(errors, "Assay 2 selection returned no cells")
    if (is.null(flux_1())) errors <- c(errors, paste0("Flux file missing for assay 1: ", input$assay_file_1))
    if (is.null(flux_2())) errors <- c(errors, paste0("Flux file missing for assay 2: ", input$assay_file_2))
    
    error_text <- if (length(errors) == 0) "NONE" else paste(errors, collapse = " - ")
    errors_global(error_text)
    
    output$reportText <- renderText({
      paste(
        paste0("Assay 1: ", input$assay_file_1, "\n"),
        paste0("Assay 2: ", input$assay_file_2, "\n"),
        paste0("Minimum AUC: ", input$customAUC, "\n"),
        paste0("Minimum Log2FC: ", input$customFC, "\n"),
        paste0("Maximum FDR: ", input$customFDR, "\n"),
        paste0("Current errors: ", error_text, "\n"),
        paste0("Report Generated at: ", Sys.time())
      )
    })
  }, ignoreNULL = FALSE)
  
  perform_analysis <- function(fdr_cutoff) {
    validate(need(errors_global() == "NONE", errors_global()))
    res <- runMarkerWrapper(fdr_cutoff)
    analysis_results(res)
    selected_subsystem(NULL)
  }
  
  runMarkerWrapper <- function(fdr_cutoff) {
    run_marker_analysis(
      flux_1 = flux_1(),
      flux_2 = flux_2(),
      names_1 = names_1(),
      names_2 = names_2(),
      auc_cutoff = input$customAUC,
      fc_cutoff = input$customFC,
      fdr_cutoff = fdr_cutoff,
      reactMeta = reactMeta
    )
  }
  
  observeEvent(input$submitReportDefault, {
    res <- runMarkerWrapper(0.05)
    analysis_results(res)
    selected_subsystem(NULL)
  })
  
  observeEvent(input$submitReport, {
    res <- runMarkerWrapper(input$customFDR)
    analysis_results(res)
    selected_subsystem(NULL)
  })
  
  metabs_fin <- reactive({
    req(analysis_results())
    subset(analysis_results()$metabs_tab, sig %in% input$outerSelect)
  })
  
  output$outerTitle <- renderText({
    req(input$assay_file_1, input$assay_file_2)
    paste0(input$assay_file_1, " v. ", input$assay_file_2)
  })
  
  output$plot1Title_metab <- renderText({
    req(input$assay_file_1, input$assay_file_2)
    paste0("Significant Reactions Across All Subsystems - (", input$assay_file_1, " v. ", input$assay_file_2, ")")
  })
  
  plot_metab_gg <- reactive({
    df <- metabs_fin()
    req(df)
    
    ggplot(
      df,
      aes(
        x = reorder(subsystem, -as.numeric(t_welch_statistic)),
        y = as.numeric(t_welch_statistic),
        color = sig,
        fill = sig,
        key = metabolite
      )
    ) +
      geom_point(size = 1.8, alpha = 0.85) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      geom_hline(yintercept = 1.3, size = .05) +
      geom_hline(yintercept = -1.3, size = .05) +
      xlab("subsystem")
  })
  
  output$plot_metab <- renderPlotly({
    p <- plot_metab_gg()
    
    ggplotly(
      p,
      source = "metab_plot",
      tooltip = c("x", "y", "key")
    ) |>
      event_register("plotly_selected") |>
      event_register("plotly_click")
  })
  
  systems <- reactive({
    req(metabs_fin())
    ggbld <- ggplot_build(
      ggplot(metabs_fin(), aes(x = reorder(subsystem, -as.numeric(t_welch_statistic)), y = as.numeric(t_welch_statistic), color = sig, fill = sig)) +
        geom_point(size = 1, position = position_dodge(width = 0.5)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    )
    ggbld$layout$panel_params[[1]]$x$limits
  })
  
  
  
  selected_points_plot <- reactive({
    sel <- event_data("plotly_selected", source = "metab_plot")
    if (is.null(sel) || nrow(sel) == 0) return(NULL)
    
    df <- metabs_fin()
    as.data.frame(df[df$metabolite %in% sel$key, , drop = FALSE])
  })
  
  active_reaction <- reactive({
    clk <- event_data("plotly_click", source = "metab_plot")
    if (!is.null(clk) && nrow(clk) > 0) return(clk$key[[1]])
    
    sel <- selected_points_plot()
    if (!is.null(sel) && nrow(sel) > 0) return(sel$metabolite[[1]])
    
    NULL
  })
  
  
  
  output$info <- renderReactable({
    df <- selected_points_plot()
    req(df)
    
    # Optionally drop sig and subsystem if you want a cleaner table
    df_tab <- df
    # Example: remove sig and subsystem columns like old code removed 5,6
    drop_cols <- intersect(c("sig"), colnames(df_tab))
    if (length(drop_cols) > 0) {
      df_tab <- df_tab[, !colnames(df_tab) %in% drop_cols, drop = FALSE]
    }
    
    # Round numeric columns similar to original
    num_cols <- intersect(c("var_1_mean", "var_2_mean", "t_welch_statistic"), colnames(df_tab))
    for (cc in num_cols) {
      df_tab[[cc]] <- signif(as.numeric(df_tab[[cc]]), 4)
    }
    
    reactable(
      df_tab,
      striped = TRUE,
      highlight = TRUE,
      compact = TRUE,
      selection = "single",
      onClick = "select"  # clicking a row selects it
    )
  })
  
  observeEvent(getReactableState("info", "selected"), {
    sel <- getReactableState("info", "selected")
    df <- selected_points_plot()
    if (!is.null(sel) && length(sel) == 1 && !is.null(df) && nrow(df) >= sel) {
      selected_reaction(df$metabolite[sel])
    }
  })
  
  observeEvent(selected_points_plot(), {
    df <- selected_points_plot()
    if (!is.null(df) && nrow(df) > 0) {
      # Only set default if nothing is selected yet
      if (is.null(selected_reaction())) {
        selected_reaction(df$metabolite[1])
      }
    } else {
      selected_reaction(NULL)
    }
  }, ignoreNULL = FALSE)
  
  observeEvent(input$invert, {
    invert_boxplot(!isTRUE(invert_boxplot()))
  })
  
  output$plot3_metab <- renderPlot({
    req(analysis_results())
    req(selected_reaction())
    df <- analysis_results()$data_filt_var
    rxn <- selected_reaction()
    xvar <- if (isTRUE(invert_boxplot())) rxn else "selected_assay"
    yvar <- if (isTRUE(invert_boxplot())) "selected_assay" else rxn
    ggplot(df, aes_string(x = xvar, y = yvar)) +
      geom_boxplot() +
      geom_jitter() +
      ylab(yvar) +
      xlab(xvar)
  })
  
  output$plot_feat_combined <- renderPlot({
    req(analysis_results())
    req(selected_reaction())
    req(assay_1_active(), assay_2_active())
    
    flux <- analysis_results()$flux
    rxn <- selected_reaction()
    
    obj1 <- AddMetaData(
      assay_1_active(),
      flux[rownames(assay_1_active()@meta.data), rxn, drop = TRUE],
      col.name = rxn
    )
    
    obj2 <- AddMetaData(
      assay_2_active(),
      flux[rownames(assay_2_active()@meta.data), rxn, drop = TRUE],
      col.name = rxn
    )
    
    shared_min <- min(c(
      obj1@meta.data[, rxn],
      obj2@meta.data[, rxn]
    ), na.rm = TRUE)
    
    shared_max <- max(c(
      obj1@meta.data[, rxn],
      obj2@meta.data[, rxn]
    ), na.rm = TRUE)
    
    plot_min(shared_min)
    plot_max(shared_max)
    
    p1 <- FeaturePlot(obj1, rxn) +
      ggtitle(input$assay_file_1) +
      scale_color_gradient2(
        low = "blue2",
        mid = "lavenderblush",
        high = "firebrick1",
        midpoint = 0,
        limits = c(shared_min, shared_max)
      )
    
    p2 <- FeaturePlot(obj2, rxn) +
      ggtitle(input$assay_file_2) +
      scale_color_gradient2(
        low = "blue2",
        mid = "lavenderblush",
        high = "firebrick1",
        midpoint = 0,
        limits = c(shared_min, shared_max)
      )
    
    cowplot::plot_grid(p1, p2, ncol = 2, align = "hv")
  })
  
  output$fingerprint <- downloadHandler(
    filename = function() {
      paste0(paste0(input$outerSelect, collapse = "_"), "_", input$assay_file_1, "_vs_", input$assay_file_2, "_filtered_reactions.rds")
    },
    content = function(file) {
      req(analysis_results())
      data_to_save <- list(
        metabs_fin(),
        analysis_results()$metaFinal,
        data.frame(AUC = input$customAUC, Log2FC = input$customFC, FDR = input$customFDR),
        "single_cell"
      )
      saveRDS(data_to_save, file)
    },
    contentType = "application/rds"
  )
  outputOptions(output, "fingerprint", suspendWhenHidden = FALSE)
  
  output$dlPlot1 <- downloadHandler(
    filename = function() {
      "all_subsystems.pdf"
    },
    content = function(file) {
      req(metabs_fin())
      p <- ggplot(metabs_fin(), aes(x = as.numeric(t_welch_statistic), y = reorder(subsystem, -as.numeric(t_welch_statistic)), color = sig, fill = sig)) +
        geom_point(size = 1, position = position_dodge(width = 0.5)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        geom_vline(xintercept = 1.3, size = 0.05) +
        geom_vline(xintercept = -1.3, size = 0.05) +
        ylab("subsystem") +
        xlab("log(pval) with sign of t-statistic") +
        ggtitle(paste0("Significant Reactions Across All Subsystems - (", input$assay_file_1, " v. ", input$assay_file_2, ")"))
      ggsave(file, p, width = 11, height = 8.5, dpi = 400, units = "in")
    },
    contentType = "application/pdf"
  )
  outputOptions(output, "dlPlot1", suspendWhenHidden = FALSE)
  
  output$dlPlot2 <- downloadHandler(
    filename = function() {
      req(selected_subsystem())
      paste0(str_replace_all(selected_subsystem(), " ", "_"), "_", input$assay_file_1, "_vs_", input$assay_file_2, ".pdf")
    },
    content = function(file) {
      req(selected_subsystem(), metabs_fin())
      p <- ggplot(subset(metabs_fin(), subsystem == selected_subsystem()), aes(x = subsystem, y = as.numeric(t_welch_statistic), color = sig, fill = sig)) +
        geom_point(size = 2) +
        geom_line(aes(group = metabolite, x = subsystem), color = "grey") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        geom_hline(yintercept = 1.3, size = 0.05) +
        geom_hline(yintercept = -1.3, size = 0.05) +
        ggtitle(paste0("Significant Reactions - ", selected_subsystem(), " (", input$assay_file_1, " v. ", input$assay_file_2, ")"))
      ggsave(file, p, width = 8, height = 12, dpi = 400, units = "in")
    }
  )
  
  output$dlPlot3 <- downloadHandler(
    filename = function() {
      paste0(input$boxplot, "_boxplot.pdf")
    },
    content = function(file) {
      req(analysis_results(), input$boxplot != "")
      df <- analysis_results()$data_filt_var
      xvar <- if (isTRUE(invert_boxplot())) input$boxplot else "selected_assay"
      yvar <- if (isTRUE(invert_boxplot())) "selected_assay" else input$boxplot
      p <- ggplot(df, aes_string(x = xvar, y = yvar)) + geom_boxplot() + geom_jitter() + ylab(yvar) + xlab(xvar)
      ggsave(file, p, width = 12, height = 8, dpi = 400, units = "in")
    }
  )
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x

shinyApp(ui, server)
