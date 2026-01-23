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
library(MatrixGenerics)


reactMeta<-read.csv("../../../../data/RADAR_xCheck_cohort/human_reaction_meta.csv")
colnames(reactMeta)[1]<-"metabolite"

cohort_registry <- read.csv("../../../../data/RADAR_xCheck_cohort/xCheck_cohort_registry.csv", stringsAsFactors = FALSE)

# beataml_flux<-as.data.frame(as.matrix(fread("beataml2_flux.csv"), rownames = 1))
# beataml_gene<-as.data.frame(as.matrix(fread("beataml2_norm_filtered.csv"), rownames = 1))
# beataml_clinical<-read.csv("beataml2_clinical.csv", header = T, row.names = 1)
# beataml_meta_clinical<-read.csv("beataml2_meta_clinical.csv", header = T)
# 
# pdac_tcga_flux<-as.data.frame(as.matrix(fread("pdac_tcga_flux.csv"), rownames = 1))
# pdac_tcga_gene<-as.data.frame(as.matrix(fread("pdac_tcga_norm_filtered.csv", header = TRUE), rownames = 1))
# pdac_tcga_clinical<-read.csv("pdac_tcga_clinical.csv", header = T, row.names = 1)
# pdac_tcga_meta_clinical<-read.csv("pdac_tcga_meta_clinical.csv", header = T)
# 
# glio_flux<-as.data.frame(as.matrix(fread("glioblastoma_flux.csv"), rownames = 1))
# glio_gene<-as.data.frame(as.matrix(fread("glioblastoma_norm_filtered.csv", header = TRUE), rownames = 1))
# glio_clinical<-read.csv("glioblastoma_clinical.csv", header = T, row.names = 1)
# glio_meta_clinical<-read.csv("glioblastoma_meta_clinical.csv", header = T)

#Assign categories to be compared

header <- dashboardHeader(
  title = "RADAR | xCheck",
  titleWidth = 400 
)

body <- dashboardBody(
  
  fluidRow(
    column(
      width = 12,
      box(
        width = 12,
        title = "Select Cohort",
        solidHeader = TRUE, status = "primary",
        box(
          width = 4,
          selectInput(
            label   = "Select Cohort",
            inputId = "cohort",
            # value = folder_name, label = "cohort_name - disease"
            choices = stats::setNames(
              cohort_registry$folder_name,
              paste0(cohort_registry$cohort_name, " - ", cohort_registry$disease)
            ),
            selected = cohort_registry$folder_name[1]
          )
        )
      )
    )
  ),
  
  fluidRow(
    column(
      width = 12,
      box(
        width = 12,
        title = "Filter by Stratification Attributes",
        solidHeader = TRUE,
        status = "primary",
        
        # Stratification type
        box(
          width = 4,
          selectInput(
            label   = "Stratification type (optional)",
            inputId = "stratType",
            choices = c(
              "Three-Class Continuous (25th / 75th percentile)" = 1,
              "Two-Class Continuous (median split)"             = 2
            ),
            selected = 1
          )
        ),
        
        # Clinical stratification
        box(
          width = 4,
          virtualSelectInput(
            label   = "Clinical stratification (leave empty for no stratification)",
            inputId = "clinStrat",
            choices = NULL,
            showValueAsTags = TRUE,
            search = TRUE,
            multiple = TRUE,
            placeholder = "No clinical stratification (optional)"
          )
        ),
        
        # Gene stratification
        box(
          width = 4,
          virtualSelectInput(
            label   = "Gene expression stratification (leave empty for no stratification)",
            inputId = "geneStrat",
            choices = NULL,
            showValueAsTags = TRUE,
            search = TRUE,
            multiple = TRUE,
            placeholder = "No gene stratification (optional)"
          )
        ),
        
        # Info text + submit button
        fluidRow(
          column(
            width = 8,
            tags$div(
              style = "padding: 10px 15px; color: #555;",
              tags$strong("Note: "),
              "If you leave both stratification fields empty, the full cohort will be used (no stratification)."
            )
          ),
          column(
            width = 4,
            div(
              style = "text-align: right; padding-right: 15px;",
              actionButton("submitStrat", "Apply stratification")
            )
          )
        )
      ),
      tags$div(id = "stratification")
    )
  ),
  
  fluidRow(column(width=8,box(title=textOutput("plot1Title"),
                              width = 12, solidHeader = TRUE, status="primary",withLoader(plotOutput("plot1", click = "clickBar", height = 500), type = "html", loader = "dnaspin"),
                              downloadLink("dlPlot1", "Download Plot as PDF | "),downloadLink("reacts", "Download Fingerprint Data")),
  ), 
  column(width=4,box(title=textOutput("plot2Title"),
                     width = 12, solidHeader = TRUE, status="primary", withLoader(plotOutput("plot2", brush = "plot_brush2", height = 500), type = "html", loader = "dnaspin"),
                     downloadLink("dlPlot2", "Download Plot as PDF")))),
  fluidRow(column(width=8,box(title=textOutput("outerTitle"), solidHeader = TRUE,status="primary",checkboxGroupInput("outerSelect", label = NULL, choices = c("Significant (+)"="up","Not Signficant"="ns","Significant (-)"="down"), selected=c("up","ns","down"))),box(title=textOutput("innerTitle"), solidHeader = TRUE, status="primary",checkboxGroupInput("innerSelect", label = NULL, choices = c("high","baseline","low"))))),
  fluidRow(column(width=12,box(width=12,selectizeInput("boxplot", "Select Metabolic Flux to Analyze", choices = c(""))))),
  fluidRow(column(width=12,box(width=12,reactableOutput("info")))),
  fluidRow(column(width=12,box(width=12,withLoader(plotOutput("plot3", height = 500), type = "html", loader = "dnaspin"),actionLink("invert","Invert Boxplot | "), downloadLink("dlPlot3", "Download Plot as PDF"))))
  
)

ui<-dashboardPage(
  skin = "purple",
  header = header,
  sidebar = dashboardSidebar(disable = TRUE), # here, we only have one tab of our app, so we don't need a sidebar
  body = body
)

server <-function(input,output,session){
  cohort_name<-reactiveVal(NULL)
  cohort_flux<-reactiveVal(NULL)
  gene<-reactiveVal(NULL)
  clinical<-reactiveVal(NULL)
  meta_clinical<-reactiveVal(NULL)
  
  choices_gene_2_class<-reactiveVal(NULL)
  choices_gene_3_class<-reactiveVal(NULL)
  choices_clinical_2_class<-reactiveVal(NULL)
  choices_clinical_3_class<-reactiveVal(NULL)
  
  type_strat<-reactiveVal(NULL)
  clinical_strat<-reactiveVal(NULL)
  gene_strat<-reactiveVal(NULL)
  
  type_outer<-reactiveVal(NULL)
  outerVarFor<-reactiveVal(NULL)
  outerVarRev<-reactiveVal(NULL)
  
  type_inner<-reactiveVal(NULL)
  inner<-reactiveVal(NULL)
  
  metaFinal_out <- reactiveVal(NULL)
  
  errors_list<-reactiveVal(NULL)
  
  strat_ui_inserted  <- reactiveVal(FALSE)
  outer_ui_inserted  <- reactiveVal(FALSE)
  inner_ui_inserted  <- reactiveVal(FALSE)
  summary_ui_inserted <- reactiveVal(FALSE)
  
  observeEvent(input$cohort, {
    req(input$cohort)
    withProgress(message = "Loading cohort data", value = 0, {
      incProgress(0.2, detail = "Reading registry")
    row <- subset(cohort_registry, folder_name == input$cohort)
    
    cohort_name(paste0(row$cohort_name, " - ", row$disease))
    
    prefix <- row$folder_name
    
    flux_file          <- paste0("../../../../data/RADAR_xCheck_cohort/",file.path(prefix, paste0(prefix, "_flux.csv")))
    gene_file          <- paste0("../../../../data/RADAR_xCheck_cohort/",file.path(prefix, paste0(prefix, "_norm_filtered.csv")))
    clinical_file      <- paste0("../../../../data/RADAR_xCheck_cohort/",file.path(prefix, paste0(prefix, "_clinical.csv")))
    meta_clinical_file <- paste0("../../../../data/RADAR_xCheck_cohort/",file.path(prefix, paste0(prefix, "_meta_clinical.csv")))
    
    incProgress(0.3, detail = "Loading flux")
    flux <- as.data.frame(as.matrix(data.table::fread(flux_file), rownames = 1))
    
    incProgress(0.2, detail = "Loading gene expression")
    gene_mat <- as.data.frame(as.matrix(data.table::fread(gene_file, header = TRUE), rownames = 1))
    
    incProgress(0.2, detail = "Loading clinical data")
    clinical_df <- read.csv(clinical_file, header = TRUE, row.names = 1)
    meta_clin_df <- read.csv(meta_clinical_file, header = TRUE)
    
    cohort_flux(flux)
    gene(gene_mat)
    clinical(clinical_df)
    meta_clinical(meta_clin_df)
    
    incProgress(0.1, detail = "Finishing")
    })
  })
  
  observeEvent(list(input$cohort,input$stratType),{
    
    withProgress(message = "Updating stratifications", value = 0, {
      req(gene(), clinical(), meta_clinical())
      
      genes <- colnames(gene())
      n_genes <- length(genes)
      
      results_genes_2 <- vector("list", n_genes)
      names(results_genes_2) <- genes
      
      for (i in seq_along(genes)) {
        results_genes_2[[i]] <- c(
          paste0(genes[i], " - high"),
          paste0(genes[i], " - low")
        )
        if (i %% 100 == 0 || i == n_genes) {
          incProgress(0.2 * i / n_genes,
                      detail = sprintf("Preparing 2-class gene stratifications (%d/%d)", i, n_genes))
        }
      }
      choices_gene_2_class(results_genes_2)
      
      results_genes_3 <- vector("list", n_genes)
      names(results_genes_3) <- genes
      for (i in seq_along(genes)) {
        results_genes_3[[i]] <- c(
          paste0(genes[i], " - high"),
          paste0(genes[i], " - baseline"),
          paste0(genes[i], " - low")
        )
        if (i %% 100 == 0 || i == n_genes) {
          incProgress(0.2 + 0.2 * i / n_genes,
                      detail = sprintf("Preparing 3-class gene stratifications (%d/%d)", i, n_genes))
        }
      }
      choices_gene_3_class(results_genes_3)
      
      # Clinical stratifications
      clinicalVars <- colnames(clinical())
      n_clin <- length(clinicalVars)
      
      results_clinical_2 <- list()
      results_clinical_3 <- list()
      
      for (i in seq_along(clinicalVars)) {
        var <- clinicalVars[i]
        if (var %in% subset(meta_clinical(), type == "continuous")$variable) {
          results_clinical_2[[var]] <- c(paste0(var, " - high"),
                                         paste0(var, " - low"))
          results_clinical_3[[var]] <- c(paste0(var, " - high"),
                                         paste0(var, " - baseline"),
                                         paste0(var, " - low"))
        } else {
          vals <- na.omit(unique(clinical()[[var]]))
          labels <- paste0(var, " - ", vals)
          results_clinical_2[[var]] <- labels
          results_clinical_3[[var]] <- labels
        }
        if (i %% 10 == 0 || i == n_clin) {
          incProgress(0.4 + 0.2 * i / max(1, n_clin),
                      detail = sprintf("Building clinical stratifications (%d/%d)", i, n_clin))
        }
      }
      
      choices_clinical_2_class(results_clinical_2)
      choices_clinical_3_class(results_clinical_3)
      
      # Final UI updates (coarse progress to 1.0)
      if (as.numeric(input$stratType) == 1) {
        updateVirtualSelect("clinStrat", choices = choices_clinical_3_class())
        updateVirtualSelect("geneStrat", choices = choices_gene_3_class())
      } else {
        updateVirtualSelect("clinStrat", choices = choices_clinical_2_class())
        updateVirtualSelect("geneStrat", choices = choices_gene_2_class())
      }
      incProgress(0.2, detail = "Finishing stratification update")
    })
  })
  
  observeEvent(input$submitStrat, {
    req(input$submitStrat > 0)
    
    if (!isTRUE(strat_ui_inserted())){
      insertUI(
        selector = "#stratification",
        where = "afterEnd",
        immediate = TRUE,
        ui = fluidRow(column(width=6,box(width=12,title="Select Outer Comparison",solidHeader = TRUE, status="primary",
                                         box(width = 12, selectInput(label="Select Outer Comparison Type", inputId = "outerType", choices = c("Three-Class Continuous (from 25th and 75th percentiles)"=1,"Two-Class Continuous (from median)"=2, "Two-Class Binary (yes/no, positive/negative)"=3, "Categorical"=4),selected=1)),
                                         box(width = 12, selectizeInput(label="Select Variable Type", inputId = "outerVarType", choices = c("Clinical Data"=1,"Gene Expression"=2),selected=1)),
                                         box(width = 6, virtualSelectInput(label="Select Forward Variable(s)", inputId = "outerVarFor", choices = colnames(clinical),showValueAsTags = TRUE, search = TRUE, multiple = TRUE)),
                                         box(width = 6, virtualSelectInput(label="Select Inverted Variable(s)", inputId = "outerVarRev", choices = colnames(clinical),showValueAsTags = TRUE, search = TRUE, multiple = TRUE)),
                                         fluidRow(column(width=12,actionButton("submitOuter", "Submit"),tags$div(id = 'outerButton'))))),tags$div(id = 'outer')
                      
        )
      )
      strat_ui_inserted(TRUE)
    }
    type_strat(as.numeric(input$stratType))
    clinical_strat(list(input$clinStrat))
    gene_strat(list(input$geneStrat))
  })
  
  observeEvent(input$outerType,{
    req(input$outerType)
    if(input$outerType == 1 | input$outerType == 2){
      updateSelectizeInput(session, "outerVarType", choices = c("Clinical Data"=1,"Gene Expression"=2),selected=1)
    }else if(input$outerType == 3){
      updateSelectizeInput(session, "outerVarType", choices = c("Clinical Data"=1),selected=1)
    }else{
      updateSelectizeInput(session, "outerVarType", choices = c("Clinical Data"=1),selected=1)
    }
  })
  
  observeEvent(list(input$cohort,input$outerType,input$outerVarType),{
    req(input$outerType)
    if(input$outerType == 1){
      if(input$outerVarType == 1){
        updateVirtualSelect("outerVarFor", label = "Select Upper Variable(s)", choices = choices_clinical_3_class()[subset(meta_clinical(), type == "continuous")$variable])
        updateVirtualSelect("outerVarRev", label = "Select Lower Variable(s)", choices = choices_clinical_3_class()[subset(meta_clinical(), type == "continuous")$variable])
      }else{
        updateVirtualSelect("outerVarFor", label = "Select Upper Variable(s)", choices = choices_gene_3_class())
        updateVirtualSelect("outerVarRev", label = "Select Lower Variable(s)", choices = choices_gene_3_class())
      }
    }else if(input$outerType == 2){
      if(input$outerVarType == 1){
        updateVirtualSelect("outerVarFor", label = "Select Upper Variable(s)", choices = choices_clinical_2_class()[subset(meta_clinical(), type == "continuous")$variable])
        updateVirtualSelect("outerVarRev", label = "Select Lower Variable(s)", choices = choices_clinical_2_class()[subset(meta_clinical(), type == "continuous")$variable])
      }else{
        updateVirtualSelect("outerVarFor", label = "Select Upper Variable(s)", choices = choices_gene_2_class())
        updateVirtualSelect("outerVarRev", label = "Select Lower Variable(s)", choices = choices_gene_2_class())
      }
    }else if(input$outerType == 3){
      updateVirtualSelect("outerVarFor", label = "Select Upper Variable(s)", choices = choices_clinical_3_class()[subset(meta_clinical(), type == "binary")$variable])
      updateVirtualSelect("outerVarRev", label = "Select Lower Variable(s)", choices = choices_clinical_3_class()[subset(meta_clinical(), type == "binary")$variable])
      
    }else{
      updateVirtualSelect("outerVarFor", label = "Select Upper Variable","outerVarType", choices = choices_clinical_3_class()[subset(meta_clinical(), type == "categorical")$variable])
      updateVirtualSelect("outerVarRev", label = "Select Lower Variable","outerVarType", choices = choices_clinical_3_class()[subset(meta_clinical(), type == "categorical")$variable])
    }
    
  })
  
  observeEvent(input$submitOuter, {
    req(input$submitOuter > 0)
    
    if (is.null(input$outerVarFor) && is.null(input$outerVarRev)) {
      showModal(
        modalDialog(
          title = "Missing outer comparison",
          "At least one outer comparison needs to be selected for upper/lower.",
          easyClose = TRUE,
          footer = modalButton("OK")
        )
      )
      return()  # do nothing else
    }
    
    if (!isTRUE(outer_ui_inserted())){
      insertUI(
        selector  = "#outer",
        where     = "beforeEnd",
        immediate = TRUE,
        ui = column(
          width = 6,
          box(
            width = 12,
            title = "Select Inner Comparison (optional)",
            solidHeader = TRUE,
            status = "primary",
            
            # Inner comparison type, with an explicit "None" option
            box(
              width = 12,
              selectInput(
                label   = "Inner comparison type",
                inputId = "innerType",
                choices = c(
                  "None (no inner comparison)"                            = 0,
                  "Three-Class Continuous (from 25th and 75th percentiles)" = 1,
                  "Two-Class Continuous (from median)"                      = 2,
                  "Two-Class Binary (yes/no, positive/negative)"           = 3,
                  "Categorical"                                            = 4
                ),
                selected = 0
              )
            ),
            
            # Variable type – disabled / ignored when innerType = 0
            box(
              width = 12,
              selectInput(
                label   = "Inner variable type",
                inputId = "innerVarType",
                choices = c("Clinical Data" = 1, "Gene Expression" = 2),
                selected = 1
              )
            ),
            
            # Inner variable – optional, with a clear 'none' placeholder
            box(
              width = 12,
              virtualSelectInput(
                label       = "Inner variable (leave empty for none)",
                inputId     = "innerVar",
                choices     = choices_clinical_3_class(),
                showValueAsTags = TRUE,
                search      = TRUE,
                multiple    = FALSE,
                placeholder = "No inner comparison (default)"
              )
            ),
            
            fluidRow(
              column(
                width = 12,
                actionButton("submitInner", "Submit inner settings")
              )
            )
          )
        ),
        tags$div(id = "inner")
      )
      outer_ui_inserted(TRUE)
    }
    type_outer(input$outerType)
    outerVarFor(list(input$outerVarFor))
    outerVarRev(list(input$outerVarRev))
  })
  
  observeEvent(input$innerType,{
    req(input$innerType)
    if(input$innerType == 1 | input$innerType == 2){
      updateSelectizeInput(session, "innerVarType", choices = c("Clinical Data"=1,"Gene Expression"=2),selected=1)
    }else if(input$innerType == 3){
      updateSelectizeInput(session, "innerVarType", choices = c("Clinical Data"=1),selected=1)
    }else{
      updateSelectizeInput(session, "innerVarType", choices = c("Clinical Data"=1),selected=1)
    }
  })
  
  observeEvent(list(input$cohort,input$innerType,input$innerVarType),{
    req(input$innerType)
    if(input$innerType == 1 | input$innerType == 2){
      if(input$innerVarType == 1){
        updateVirtualSelect("innerVar", label = "Select Variable", choices = subset(meta_clinical(), type == "continuous")$variable)
      }else{
        updateVirtualSelect("innerVar", label = "Select Variable", choices = colnames(gene()))
      }
    }else if(input$innerType == 3){
      updateVirtualSelect("innerVar", label = "Select Variable", choices = subset(meta_clinical(), type == "binary")$variable)
    }else{
      updateVirtualSelect("innerVar", label = "Select Variable", choices = subset(meta_clinical(), type == "categorical")$variable)
    }
  })
  
  observeEvent(input$submitInner, {
    req(input$submitInner > 0)
    if (!isTRUE(inner_ui_inserted())){
      insertUI(
        selector = "#outer",
        where = "afterEnd",
        immediate = TRUE,
        ui = fluidRow(column(width=12,box(width=12,title="RADAR | xCheck Pre-Analysis Summary",solidHeader = TRUE, status="primary", 
                                          box(width = 6, verbatimTextOutput('reportText')),
                                          box(width = 6, noUiSliderInput(inputId = "customAUC", label = "AUC:", min = 0.6, max = 1.0, value = 0.65, tooltip = TRUE, step = 0.01),
                                              noUiSliderInput(inputId = "customFC", label = "Log2FC:", min = 0, max = 5, value = 1, tooltip = TRUE, step = 0.01),
                                              noUiSliderInput(inputId = "customFDR", label = "FDR:", min = 0, max = 1, value = 0.05, tooltip = TRUE, step = 0.01)),
                                          fluidRow(column(width=12,actionButton("submitReport", "Begin Analysis (set parameters)")),column(width=12,actionButton("submitReportDefault", "Begin Analysis (default parameters)")))))
                      
        )
      )
      inner_ui_inserted(TRUE)
    }
      type_inner(input$innerType)
      inner(input$innerVar)

  })
  
  observeEvent(
    list(input$submitStrat, input$submitOuter, input$submitInner),
    {
      # Require that stratification, outer and inner settings exist
      req(type_strat(), type_outer(), type_inner())
      
      warnings <- character(0)
      errors   <- character(0)
      
      ########################
      # 1. STRATIFIED COHORT #
      ########################
      
      if (is.null(clinical_strat()[[1]]) && is.null(gene_strat()[[1]])) {
        
        # No stratification: keep all samples
        strat_names <- rownames(gene())
        
      } else if (!is.null(clinical_strat()[[1]]) && is.null(gene_strat()[[1]])) {
        
        # Clinical-only stratification
        clinical_keys <- lapply(clinical_strat(), function(x) {
          split_elements <- unlist(strsplit(x, " - "))
          split_elements[seq_along(split_elements) %% 2 == 1]
        })
        clinical_values <- lapply(clinical_strat(), function(x) {
          split_elements <- unlist(strsplit(x, " - "))
          split_elements[seq_along(split_elements) %% 2 == 0]
        })
        
        dictionary <- vector("list", length(unique(clinical_keys[[1]])))
        names(dictionary) <- unique(clinical_keys[[1]])
        for (i in seq_along(clinical_keys[[1]])) {
          key   <- clinical_keys[[1]][i]
          value <- clinical_values[[1]][i]
          dictionary[[key]] <- c(dictionary[[key]], value)
        }
        
        strat_df <- as.data.frame(clinical()[, unique(clinical_keys[[1]])])
        colnames(strat_df) <- unique(clinical_keys[[1]])
        rownames(strat_df) <- rownames(clinical())
        
        for (i in seq_len(ncol(strat_df))) {
          if (colnames(strat_df)[i] %in% subset(meta_clinical(), type == "continuous")$variable) {
            if (type_strat() == 1) {
              q25 <- quantile(strat_df[, i], 0.25, na.rm = TRUE)
              q75 <- quantile(strat_df[, i], 0.75, na.rm = TRUE)
              strat_df[, i] <- ifelse(
                strat_df[, i] > q75, "high",
                ifelse(strat_df[, i] >= q25, "baseline", "low")
              )
            } else {
              q50 <- quantile(strat_df[, i], 0.5, na.rm = TRUE)
              strat_df[, i] <- ifelse(strat_df[, i] > q50, "high", "low")
            }
          }
        }
        
        for (i in seq_along(clinical_keys[[1]])) {
          strat_df <- subset(
            strat_df,
            subset = strat_df[, clinical_keys[[1]][i]] %in% dictionary[[clinical_keys[[1]][i]]]
          )
        }
        strat_names <- rownames(strat_df)
        
      } else if (is.null(clinical_strat()[[1]]) && !is.null(gene_strat()[[1]])) {
        
        # Gene-only stratification
        gene_keys <- lapply(gene_strat(), function(x) {
          split_elements <- unlist(strsplit(x, " - "))
          split_elements[seq_along(split_elements) %% 2 == 1]
        })
        gene_values <- lapply(gene_strat(), function(x) {
          split_elements <- unlist(strsplit(x, " - "))
          split_elements[seq_along(split_elements) %% 2 == 0]
        })
        
        dictionary <- vector("list", length(unique(gene_keys[[1]])))
        names(dictionary) <- unique(gene_keys[[1]])
        for (i in seq_along(gene_keys[[1]])) {
          key   <- gene_keys[[1]][i]
          value <- gene_values[[1]][i]
          dictionary[[key]] <- c(dictionary[[key]], value)
        }
        
        strat_df <- as.data.frame(gene()[, unique(gene_keys[[1]])])
        rownames(strat_df) <- rownames(gene())
        colnames(strat_df) <- unique(gene_keys[[1]])
        
        for (i in seq_len(ncol(strat_df))) {
          if (type_strat() == 1) {
            q25 <- quantile(strat_df[, i], 0.25, na.rm = TRUE)
            q75 <- quantile(strat_df[, i], 0.75, na.rm = TRUE)
            strat_df[, i] <- ifelse(
              strat_df[, i] > q75, "high",
              ifelse(strat_df[, i] >= q25, "baseline", "low")
            )
          } else {
            q50 <- quantile(strat_df[, i], 0.5, na.rm = TRUE)
            strat_df[, i] <- ifelse(strat_df[, i] > q50, "high", "low")
          }
        }
        
        for (i in seq_along(gene_keys[[1]])) {
          strat_df <- subset(
            strat_df,
            subset = strat_df[, gene_keys[[1]][i]] %in% dictionary[[gene_keys[[1]][i]]]
          )
        }
        strat_names <- rownames(strat_df)
        
      } else {
        
        # Clinical + gene stratification
        clinical_keys <- lapply(clinical_strat(), function(x) {
          split_elements <- unlist(strsplit(x, " - "))
          split_elements[seq_along(split_elements) %% 2 == 1]
        })
        clinical_values <- lapply(clinical_strat(), function(x) {
          split_elements <- unlist(strsplit(x, " - "))
          split_elements[seq_along(split_elements) %% 2 == 0]
        })
        gene_keys <- lapply(gene_strat(), function(x) {
          split_elements <- unlist(strsplit(x, " - "))
          split_elements[seq_along(split_elements) %% 2 == 1]
        })
        gene_values <- lapply(gene_strat(), function(x) {
          split_elements <- unlist(strsplit(x, " - "))
          split_elements[seq_along(split_elements) %% 2 == 0]
        })
        
        dictionary_clinical <- vector("list", length(unique(clinical_keys[[1]])))
        names(dictionary_clinical) <- unique(clinical_keys[[1]])
        for (i in seq_along(clinical_keys[[1]])) {
          key   <- clinical_keys[[1]][i]
          value <- clinical_values[[1]][i]
          dictionary_clinical[[key]] <- c(dictionary_clinical[[key]], value)
        }
        
        dictionary_gene <- vector("list", length(unique(gene_keys[[1]])))
        names(dictionary_gene) <- unique(gene_keys[[1]])
        for (i in seq_along(gene_keys[[1]])) {
          key   <- gene_keys[[1]][i]
          value <- gene_values[[1]][i]
          dictionary_gene[[key]] <- c(dictionary_gene[[key]], value)
        }
        
        strat_df_clinical <- as.data.frame(clinical()[, unique(clinical_keys[[1]])])
        rownames(strat_df_clinical) <- rownames(clinical())
        colnames(strat_df_clinical) <- unique(clinical_keys[[1]])
        
        for (i in seq_len(ncol(strat_df_clinical))) {
          if (colnames(strat_df_clinical)[i] %in% subset(meta_clinical(), type == "continuous")$variable) {
            if (type_strat() == 1) {
              q25 <- quantile(strat_df_clinical[, i], 0.25, na.rm = TRUE)
              q75 <- quantile(strat_df_clinical[, i], 0.75, na.rm = TRUE)
              strat_df_clinical[, i] <- ifelse(
                strat_df_clinical[, i] > q75, "high",
                ifelse(strat_df_clinical[, i] >= q25, "baseline", "low")
              )
            } else {
              q50 <- quantile(strat_df_clinical[, i], 0.5, na.rm = TRUE)
              strat_df_clinical[, i] <- ifelse(strat_df_clinical[, i] > q50, "high", "low")
            }
          }
        }
        
        strat_df_gene <- as.data.frame(gene()[, unique(gene_keys[[1]])])
        rownames(strat_df_gene) <- rownames(gene())
        colnames(strat_df_gene) <- unique(gene_keys[[1]])
        
        for (i in seq_len(ncol(strat_df_gene))) {
          if (type_strat() == 1) {
            q25 <- quantile(strat_df_gene[, i], 0.25, na.rm = TRUE)
            q75 <- quantile(strat_df_gene[, i], 0.75, na.rm = TRUE)
            strat_df_gene[, i] <- ifelse(
              strat_df_gene[, i] > q75, "high",
              ifelse(strat_df_gene[, i] >= q25, "baseline", "low")
            )
          } else {
            q50 <- quantile(strat_df_gene[, i], 0.5, na.rm = TRUE)
            strat_df_gene[, i] <- ifelse(strat_df_gene[, i] > q50, "high", "low")
          }
        }
        
        strat_df <- cbind(strat_df_clinical, strat_df_gene)
        dictionary <- merge(dictionary_clinical, dictionary_gene, all = TRUE)
        
        for (i in seq_along(clinical_keys[[1]])) {
          strat_df <- subset(
            strat_df,
            subset = strat_df[, clinical_keys[[1]][i]] %in% dictionary[[clinical_keys[[1]][i]]]
          )
        }
        for (i in seq_along(gene_keys[[1]])) {
          strat_df <- subset(
            strat_df,
            subset = strat_df[, gene_keys[[1]][i]] %in% dictionary[[gene_keys[[1]][i]]]
          )
        }
        strat_names <- rownames(strat_df)
      }
      
      if (length(strat_names) < 2) {
        errors <- c(errors, "Stratified Cohort Must Contain at Least Two Samples")
        metaFinal_out(NULL)
        errors_list(errors)
        output$reportText <- renderText(paste(
          "Cohort:", cohort_name(), "\n",
          "Errors:", paste(errors, collapse = "; "), "\n"
        ))
        return()
      }
      
      ##############################
      # 2. OUTER DF AND OUTER FLAG #
      ##############################
      
      outer_df <- cbind(clinical(), gene())[strat_names, , drop = FALSE]
      
      upper_keys <- lapply(outerVarFor(), function(x) {
        split_elements <- unlist(strsplit(x, " - "))
        split_elements[seq_along(split_elements) %% 2 == 1]
      })
      upper_values <- lapply(outerVarFor(), function(x) {
        split_elements <- unlist(strsplit(x, " - "))
        split_elements[seq_along(split_elements) %% 2 == 0]
      })
      lower_keys <- lapply(outerVarRev(), function(x) {
        split_elements <- unlist(strsplit(x, " - "))
        split_elements[seq_along(split_elements) %% 2 == 1]
      })
      lower_values <- lapply(outerVarRev(), function(x) {
        split_elements <- unlist(strsplit(x, " - "))
        split_elements[seq_along(split_elements) %% 2 == 0]
      })
      
      total_keys <- unique(c(upper_keys[[1]], lower_keys[[1]]))
      outer_df   <- as.data.frame(outer_df[, total_keys, drop = FALSE])
      rownames(outer_df) <- strat_names
      colnames(outer_df) <- total_keys
      
      dictionary_upper <- vector("list", length(unique(upper_keys[[1]])))
      names(dictionary_upper) <- unique(upper_keys[[1]])
      for (i in seq_along(upper_keys[[1]])) {
        key   <- upper_keys[[1]][i]
        value <- upper_values[[1]][i]
        dictionary_upper[[key]] <- c(dictionary_upper[[key]], value)
      }
      
      dictionary_lower <- vector("list", length(unique(lower_keys[[1]])))
      names(dictionary_lower) <- unique(lower_keys[[1]])
      for (i in seq_along(lower_keys[[1]])) {
        key   <- lower_keys[[1]][i]
        value <- lower_values[[1]][i]
        dictionary_lower[[key]] <- c(dictionary_lower[[key]], value)
      }
      
      # Define continuous clinical variables once
      continuous_clinical <- if ("continuous" %in% meta_clinical()[, "type"]) {
        subset(meta_clinical(), type == "continuous")$variable
      } else {
        character(0)
      }
      
      # Discretize outer_df variables and map to high/baseline/low or high/low
      for (i in seq_len(ncol(outer_df))) {
        var_name <- colnames(outer_df)[i]
        
        if (var_name %in% continuous_clinical || var_name %in% colnames(gene())) {
          
          x_num <- suppressWarnings(as.numeric(outer_df[, i]))
          if (all(is.na(x_num))) {
            next
          }
          
          if (type_outer() == 1) {
            q25 <- quantile(x_num, 0.25, na.rm = TRUE)
            q75 <- quantile(x_num, 0.75, na.rm = TRUE)
            outer_df[, i] <- ifelse(
              x_num > q75, "high",
              ifelse(x_num >= q25, "baseline", "low")
            )
          } else {
            q50 <- quantile(x_num, 0.5, na.rm = TRUE)
            outer_df[, i] <- ifelse(x_num > q50, "high", "low")
          }
        }
      }
      
      # Map to upper/lower using dictionaries
      for (i in seq_along(total_keys)) {
        key <- total_keys[i]
        if (!is.null(dictionary_upper[[key]])) {
          outer_df[, i] <- ifelse(outer_df[, i] %in% dictionary_upper[[key]], "upper", outer_df[, i])
        }
        if (!is.null(dictionary_lower[[key]])) {
          outer_df[, i] <- ifelse(outer_df[, i] %in% dictionary_lower[[key]], "lower", outer_df[, i])
        }
      }
      
      outer_df <- na.omit(outer_df)
      if (nrow(outer_df) == 0) {
        errors <- c(errors, "No samples remain after outer discretization")
        metaFinal_out(NULL)
        errors_list(errors)
        output$reportText <- renderText(paste(
          "Cohort:", cohort_name(), "\n",
          "Errors:", paste(errors, collapse = "; "), "\n"
        ))
        return()
      }
      
      # Determine rows that are consistently upper or lower
      v1 <- do.call(pmin, c(outer_df, na.rm = TRUE))
      v2 <- do.call(pmax, c(outer_df, na.rm = TRUE))
      outer_label <- ifelse(v1 == v2, v1, NA_character_)
      
      valid_rows <- which(outer_label %in% c("upper", "lower"))
      if (length(valid_rows) == 0) {
        errors <- c(errors, "No samples remain after applying outer comparison filters")
        metaFinal_out(NULL)
        errors_list(errors)
        output$reportText <- renderText(paste(
          "Cohort:", cohort_name(), "\n",
          "Errors:", paste(errors, collapse = "; "), "\n"
        ))
        return()
      }
      
      outer_df    <- outer_df[valid_rows, , drop = FALSE]
      outer_label <- outer_label[valid_rows]
      
      #############################
      # 3. INNER DF FROM SCRATCH  #
      #############################
      
      inner_df <- data.frame(
        outer = outer_label,
        row.names = rownames(outer_df),
        stringsAsFactors = FALSE
      )
      
      inner_choices <- cbind(clinical(), gene())[rownames(inner_df), , drop = FALSE]
      
      if (inner() == "" || type_inner() == 0) {
        inner_df$inner <- "All"
      } else {
        if (!inner() %in% colnames(inner_choices)) {
          errors <- c(errors, paste0("Inner variable '", inner(), "' not found in current cohort"))
          metaFinal_out(NULL)
          errors_list(errors)
          output$reportText <- renderText(paste(
            "Cohort:", cohort_name(), "\n",
            "Errors:", paste(errors, collapse = "; "), "\n"
          ))
          return()
        }
        inner_df$inner <- inner_choices[, inner(), drop = TRUE]
        
        if (inner() %in% continuous_clinical || inner() %in% colnames(gene())) {
          if (type_inner() == 1) {
            q25 <- quantile(as.numeric(inner_df$inner), 0.25, na.rm = TRUE)
            q75 <- quantile(as.numeric(inner_df$inner), 0.75, na.rm = TRUE)
            inner_df$inner <- ifelse(
              as.numeric(inner_df$inner) > q75, "high",
              ifelse(as.numeric(inner_df$inner) >= q25, "baseline", "low")
            )
          } else if (type_inner() == 2) {
            q50 <- quantile(as.numeric(inner_df$inner), 0.5, na.rm = TRUE)
            inner_df$inner <- ifelse(as.numeric(inner_df$inner) > q50, "high", "low")
          }
          # type_inner() 3/4: leave categorical/binary as-is
        }
      }
      
      inner_df <- na.omit(inner_df)
      if (nrow(inner_df) == 0) {
        errors <- c(errors, "No samples remain after applying inner comparison filters")
        metaFinal_out(NULL)
        errors_list(errors)
        output$reportText <- renderText(paste(
          "Cohort:", cohort_name(), "\n",
          "Errors:", paste(errors, collapse = "; "), "\n"
        ))
        return()
      }
      
      ##################################
      # 4. STORE META AND UPDATE TEXT  #
      ##################################
      
      metaFinal_out(inner_df)
      errors_list(errors)
      
      clinical_strat_text <- paste(clinical_strat(), sep = ", ")
      gene_strat_text     <- paste(gene_strat(),     sep = ", ")
      outer_upper_text    <- paste(outerVarFor(),    sep = ", ")
      outer_lower_text    <- paste(outerVarRev(),    sep = ", ")
      
      output$reportText <- renderText({
        paste(
          paste0("Cohort: ", cohort_name(), "\n"),
          paste0("Patient Stratifications (Clinical): ", clinical_strat_text, "\n"),
          paste0("Patient Stratifications (Gene): ", gene_strat_text, "\n"),
          paste0("Outer Comparisons (Upper): ", outer_upper_text, "\n"),
          paste0("Outer Comparisons (Lower): ", outer_lower_text, "\n"),
          paste0("Inner Comparisons: ", inner(), "\n"),
          paste0("Minimum AUC: ", input$customAUC, "\n"),
          paste0("Minimum Log2FC: ", input$customFC, "\n"),
          paste0("Maximum FDR: ", input$customFDR, "\n"),
          paste0("Errors: ", paste(errors, collapse = "; "), "\n"),
          paste0("Report Generated at: ", Sys.time())
        )
      })
    }
  )
  
  observeEvent(input$submitReport,{
    if (length(errors_list()) > 0) {
      showModal(
        modalDialog(
          title = "Cannot run analysis",
          paste("Errors present:", paste(errors_list(), collapse = "; ")),
          easyClose = TRUE
        )
      )
      return()
    }
    
    metaFinal <- metaFinal_out()
    if (is.null(metaFinal) || nrow(metaFinal) == 0) {
      showModal(
        modalDialog(
          title = "No filtered cohort",
          "Please configure stratification, outer and inner comparisons, then run the pre-analysis summary before starting the analysis.",
          easyClose = TRUE
        )
      )
      return()
    }
    
    flux_full <- cbind(metaFinal, cohort_flux()[rownames(metaFinal), ])
    flux_comp <- flux_full
    
    outer <- colnames(metaFinal)[1]
    inner <- colnames(metaFinal)[2]
    
    outer_for_vec <- if (!is.null(outerVarFor())) outerVarFor()[[1]] else character(0)
    outer_rev_vec <- if (!is.null(outerVarRev())) outerVarRev()[[1]] else character(0)
    outer_all     <- c(outer_for_vec, outer_rev_vec)
    
    parse_outer <- function(x) {
      parts <- strsplit(x, " - ", fixed = TRUE)[[1]]
      if (length(parts) >= 2) list(var = parts[1], level = parts[2]) else list(var = x, level = x)
    }
    
    # Base comparison name (e.g. "FLT3.ITD")
    outer_comp_name <- if (length(outer_all) > 0) {
      parse_outer(outer_all[1])$var
    } else {
      "Outer comparison"
    }
    
    # Human-readable upper / lower labels
    upper_label <- if (length(outer_for_vec) > 0) parse_outer(outer_for_vec[1])$level else "upper"
    lower_label <- if (length(outer_rev_vec) > 0) parse_outer(outer_rev_vec[1])$level else "lower"
    
      #   if(input$outerType != 1){
      outer_select_1<-c("upper", "lower")
      inner_opts<-unique(metaFinal[,"inner"])
      metab_names<-c()
      data_filt<-flux_comp
      
      for (i in 1:length(inner_opts)){
        data_filt_var<-subset(data_filt, data_filt[,"inner"] == inner_opts[[i]])
        AUC = colAUC(data_filt_var[,-c(1:dim(metaFinal)[2])],factor(data_filt_var[,outer]),plotROC=FALSE)
        
        AUC = AUC[1,]
        #plot(density(AUC),main=paste0("AUC for ", inner_opts[i]," features"))
        logFC = log2(apply(subset(data_filt_var, subset = data_filt_var[,"outer"] == "upper")[,-c(1:dim(metaFinal)[2])],2,mean) / apply(subset(data_filt_var, subset = data_filt_var[,"outer"] == "lower")[,-c(1:dim(metaFinal)[2])],2,mean))
        logFC[AUC<input$customAUC]=0
        short<-as.data.frame(logFC[which(abs(logFC)>input$customFC)])
        t<-col_t_welch(as.matrix(subset(data_filt_var, subset = data_filt_var[,"outer"] == "upper")[,rownames(short)]),as.matrix(subset(data_filt_var, subset = data_filt_var[,"outer"] == "lower")[,rownames(short)]))
        final<-cbind(t$mean.y, t$mean.x, short, t$statistic, t$pvalue)
        colnames(final)<-c("var_1_mean","var_2_mean","log2FC","t_welch_statistic","pval")
        
        #final<-subset(final, subset = pval < 0.05)
        final$fdr = p.adjust(final[,"pval"],method="fdr")
        final<-subset(final, subset = fdr < input$customFDR)
        metab_names<-c(metab_names, rownames(final))
      }
      metab_names<-unique(metab_names)
      
      test_list<-data.frame()
      for (i in 1:length(inner_opts)){
        names<-c(outer, metab_names)
        inner_comparison<-inner_opts[i]
        data_filt_var<-subset(data_filt, subset = data_filt[,"inner"] == inner_opts[i])
        data_filt_var<-data_filt_var[,names]
        t<-col_t_welch(as.matrix(subset(data_filt_var, subset = data_filt_var[,"outer"] == "upper")[,metab_names]),as.matrix(subset(data_filt_var, subset = data_filt_var[,"outer"] == "lower")[,metab_names]))
        final<-as.data.frame(cbind(metabolite=metab_names,inner_comparison=inner_comparison,as.numeric(t$mean.y), as.numeric(t$mean.x),as.numeric(t$statistic), as.numeric(t$pvalue)))
        colnames(final)<-c("metabolite","inner_comparison","var_1_mean","var_2_mean","t_welch_statistic","pval")
        final$sign_t_log_pval = ifelse(final$t_welch_statistic>0, -log10(as.numeric(final$pval)), log10(as.numeric(final$pval)))
        test_list<-as.data.frame(rbind(test_list,final))
      }
      
      metabs_tab <- merge(test_list, reactMeta, by = "metabolite", all.x = TRUE, all.y = FALSE)
      # If there was no inner comparison, ensure the column exists and is "All"
      
      metabs_tab$sig <- ifelse(
        metabs_tab$sign_t_log_pval > 1.301, "up",
        ifelse(metabs_tab$sign_t_log_pval < -1.301, "down", "ns")
      )
      
      # Re-level inner_comparison depending on innerType
      metabs_tab$inner_comparison <- switch(
        as.numeric(input$innerType),
        # 1: three-class continuous
        factor(metabs_tab$inner_comparison, levels = c("high", "baseline", "low")),
        # 2: two-class continuous
        factor(metabs_tab$inner_comparison, levels = c("high", "low")),
        # 3: binary
        factor(metabs_tab$inner_comparison, levels = unique(metabs_tab$inner_comparison)),
        # 4: categorical
        factor(metabs_tab$inner_comparison, levels = unique(metabs_tab$inner_comparison))
      )
      
      if (inner() == "") {
        metabs_tab$inner_comparison <- "All"
      }
      
      updateCheckboxGroupInput(session,"innerSelect", choices = unique(data_filt[,2]), selected = unique(data_filt[,2]))
      
      output$outerTitle<-renderText(outer)
      output$innerTitle<-renderText(inner)
      
      output$plot1Title <- renderText(
        paste0(
          "Significant Reactions Across All Subsystems - ",
          outer_comp_name,
          " (", upper_label, " v. ", lower_label, ")"
        )
      )
      
      metabs_fin <- reactive({
        req(exists("metabs_tab"), !is.null(metabs_tab))
        
        if (!"inner_comparison" %in% colnames(metabs_tab)) {
          return(metabs_tab)
        }
        
        if (nrow(metabs_tab) == 0) {
          return(metabs_tab)
        }
        
        inner_levels <- unique(metabs_tab[,"inner_comparison", drop = TRUE])
        
        if (length(inner_levels) == 1) {
          subset(metabs_tab, sig %in% input$outerSelect)
        } else {
          subset(
            metabs_tab,
            inner_comparison %in% input$innerSelect & sig %in% input$outerSelect
          )
        }
      })
      
      plot1 <- reactive({
        if(length(unique(metabs_fin()[,"inner_comparison"])) == 1){
          ggplot(metabs_fin(), aes(x=reorder(subsystem, -sign_t_log_pval), y=sign_t_log_pval, color=sig, fill = sig)) + geom_point(size = 1, position = position_dodge(width=0.5))+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)+xlab("subsystem")
        }else{
          ggplot(metabs_fin(), aes(x=reorder(subsystem, -sign_t_log_pval), y=sign_t_log_pval, color=inner_comparison, fill = inner_comparison)) + geom_point(size = 1, position = position_dodge(width=0.5))+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)+xlab("subsystem")
        }
        
      })
      
      output$plot1<- renderPlot({
        plot1()
      })
      
      plot1pdf <- reactive({
        ggplot(metabs_fin(), aes(x=as.numeric(t_welch_statistic), y=reorder(subsystem, -as.numeric(t_welch_statistic)), color=sig, fill = sig)) + geom_point(size = 1, position = position_dodge(width=0.5))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_vline(xintercept=1.3, size=.05)+geom_vline(xintercept=-1.3, size=.05)+ylab("subsystem")+xlab("log(pval) with sign of t-statistic")+ggtitle(paste0("Significant Reactions Across All Subsystems - (assay_1 v. assay_2)"))
      })
      
      output$reacts <- downloadHandler(
        filename = function() {
          paste(paste0(input$innerSelect, collapse = "_"),"_",inner, "_",paste0(input$outerSelect, collapse = "_"), outer, '_filtered_reactions.rds', sep='')
        },
        content = function(file) {
          saveRDS(list(metabs_fin(),metaFinal,data.frame(AUC=input$customAUC,Log2FC=input$customFC,FDR=input$customFDR),"cohort"),file)
        }
      )
      
      # output$reacts_meta <- downloadHandler(
      #   filename = function() {
      #     paste(paste0(input$innerSelect, collapse = "_"),"_",inner, "_",paste0(input$outerSelect, collapse = "_"), outer, '_meta_data.csv', sep='')
      #   },
      #   content = function(file) {
      #     write.csv(metaFinal,file)
      #   }
      # )
      
      systems<-reactive({
        if(length(unique(metabs_fin()[,"inner_comparison"])) == 1){
          ggbld <- ggplot_build(ggplot(metabs_fin(), aes(x=reorder(subsystem, -sign_t_log_pval), y=sign_t_log_pval, color=sig, fill = sig)) + geom_point(size = 1, position = position_dodge(width=0.5))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)+xlab("subsystem"))
          ggbld$layout$panel_params[[1]]$x$limits
        }else{
          ggbld <- ggplot_build(ggplot(metabs_fin(), aes(x=reorder(subsystem, -sign_t_log_pval), y=sign_t_log_pval, color=inner_comparison, fill = inner_comparison)) + geom_point(size = 1, position = position_dodge(width=0.5))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)+xlab("subsystem"))
          ggbld$layout$panel_params[[1]]$x$limits
        }
        
      })
      
      
      output$dlPlot1 <- downloadHandler(
        filename = function() {
          paste('all_subsystems.pdf', sep='')
        },
        content = function(file) {
          withProgress(message = 'Generating PDF', value = 0, {
            incProgress(0.3, detail = "Preparing data")
            
            pdf(file, width = 11, height = 8.5)  # Standard letter size
            
            # Plot or write your content to the PDF here
            # For example:
            plot(plot1pdf())
            # Add more plots or text as needed
            
            incProgress(0.3, detail = "Creating PDF")
            
            # Close the PDF device
            dev.off()
            
            incProgress(0.4, detail = "Finalizing")
          })
        },
        contentType = "application/pdf"
      )
      
      observeEvent(input$clickBar, {
        groupId <- round(input$clickBar$x)
        output$plot2Title <- renderText(
          paste0(
            "Significant Reactions ", outer_comp_name, " - ", systems()[groupId],
            " (", upper_label, " v. ", lower_label, ")"
          )
        )
        plot2 <- reactive({
          if(length(unique(metabs_fin()[,"inner_comparison"])) == 1){
            ggplot(subset(metabs_fin(), subset = subsystem == systems()[groupId]), aes(x=inner_comparison, y=sign_t_log_pval, color=sig, fill = sig)) + geom_point(size = 2)+geom_line(aes(group = metabolite, x=inner_comparison),color="grey")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)
          }else{
            ggplot(subset(metabs_fin(), subset = subsystem == systems()[groupId]), aes(x=inner_comparison, y=sign_t_log_pval, color=inner_comparison, fill = inner_comparison)) + geom_point(size = 2)+geom_line(aes(group = metabolite, x=inner_comparison),color="grey")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)
          }
          
        })
        output$plot2<- renderPlot({
          plot2()
        })
        plot2pdf <- reactive({
          if(length(unique(metabs_fin()[,"inner_comparison"])) == 1){
            ggplot(subset(metabs_fin(), subset = subsystem == systems()[groupId]), aes(x=inner_comparison, y=sign_t_log_pval, color=sig, fill = sig)) + geom_point(size = 2)+geom_line(aes(group = metabolite, x=inner_comparison),color="grey")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)+ggtitle(paste0("Significant Reactions ",outer, " - ",systems[groupId], " (upper v. lower)"))
          }else{
            ggplot(subset(metabs_fin(), subset = subsystem == systems()[groupId]), aes(x=inner_comparison, y=sign_t_log_pval, color=inner_comparison, fill = inner_comparison)) + geom_point(size = 2)+geom_line(aes(group = metabolite, x=inner_comparison),color="grey")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)+ggtitle(paste0("Significant Reactions ",outer, " - ",systems[groupId], " (upper v. lower)"))
          }
          
        })
        output$dlPlot2 <- downloadHandler(
          filename = function() {
            paste(str_replace(systems()[groupId]," ","_"),'_',outer,'.pdf', sep='')
          },
          content = function(file) {
            ggsave(file, plot2pdf(), width = 8, height = 12, dpi = 400, units = "in")
          }
        )
        output$info<-renderReactable({
          metabs_fin<-metabs_fin()
          metabs_fin$var_1_mean<-round(as.numeric(metabs_fin$var_1_mean), digits=4)
          metabs_fin$var_2_mean<-round(as.numeric(metabs_fin$var_2_mean), digits=4)
          metabs_fin$sign_t_log_pval<-round(as.numeric(metabs_fin$sign_t_log_pval), digits=4)
          reactable(brushedPoints(subset(metabs_fin, subset = subsystem == systems()[groupId]), input$plot_brush2, xvar = "inner_comparison", yvar="sign_t_log_pval"))
        })
        data <- reactive({
          as.data.frame(brushedPoints(subset(metabs_fin, subset = subsystem == systems()[groupId]), input$plot_brush2, xvar = "inner_comparison", yvar="sign_t_log_pval"))
        })
        output$data <- downloadHandler(
          filename = function() {
            paste(str_replace(systems()[groupId]," ","_"),'_selected_reactions_',outer, '.csv', sep='')
          },
          content = function(file) {
            write.csv(data(),file)
          }
        )
        
        observe({
          updateSelectizeInput(session, "boxplot",
                               choices = brushedPoints(subset(metabs_fin(), subset = subsystem == systems()[groupId]), input$plot_brush2, xvar = "inner_comparison", yvar="sign_t_log_pval")$metabolite
          )})
      })
 
      plot3 <- reactive({
        if (input$boxplot == "") {
          return(ggplot())
        }
        
        metaFinal <- metaFinal_out()
        req(metaFinal)
        
        outer_var <- colnames(metaFinal)[1]
        inner_var <- colnames(metaFinal)[2]
        
        outer_type <- as.numeric(input$outerType)
        has_inner  <- !(is.null(inner()) || inner() == "" || type_inner() == 0)
        inner_name <- if (has_inner) inner() else "None"
        
        # Base data
        df <- data_filt[, c(colnames(metaFinal), input$boxplot), drop = FALSE]
        colnames(df)[1] <- "outer_group"
        colnames(df)[2] <- "inner"
        
        # ---------- Parse outer selections for labels ----------
        
        outer_for_vec <- if (!is.null(outerVarFor())) outerVarFor()[[1]] else character(0)
        outer_rev_vec <- if (!is.null(outerVarRev())) outerVarRev()[[1]] else character(0)
        outer_all     <- c(outer_for_vec, outer_rev_vec)
        
        parse_outer <- function(x) {
          parts <- strsplit(x, " - ", fixed = TRUE)[[1]]
          if (length(parts) >= 2) list(var = parts[1], level = parts[2]) else list(var = x, level = x)
        }
        
        # Default labels
        cont_name       <- NULL
        outer_level_map <- NULL
        outer_levels    <- sort(unique(df$outer_group))
        
        if (length(outer_all) > 0) {
          parsed <- lapply(outer_all, parse_outer)
          cont_name <- parsed[[1]]$var                    # e.g. "FLT3.ITD"
          for_levels <- vapply(outer_for_vec, function(z) parse_outer(z)$level, character(1))
          rev_levels <- vapply(outer_rev_vec, function(z) parse_outer(z)$level, character(1))
          outer_levels <- unique(c(for_levels, rev_levels))
          
          # Map original outer_group values to pretty labels:
          # assume outer_group encodes upper/lower in your metaFinal,
          # but here we simply relabel based on order of levels we derived
          if (length(outer_levels) >= 2) {
            outer_level_map <- setNames(outer_levels, c("upper", "lower")[seq_along(outer_levels)])
          }
        }
        
        # Apply outer label mapping to df$outer_group (if available)
        df$outer_group <- factor(as.character(df$outer_group))
        if (!is.null(outer_level_map)) {
          # Map factor levels "upper"/"lower" (or whatever you used) to user-friendly labels
          df$outer_group <- factor(
            outer_level_map[as.character(df$outer_group)],
            levels = outer_levels
          )
        }
        
        # ---------- BOX PLOTS (original behavior, now with pretty x labels) ----------
        flux_col <- input$boxplot
        
        pval_to_stars <- function(p) {
          if (is.na(p)) return("ns")
          if (p < 0.001) return("***")
          if (p < 0.01)  return("**")
          if (p < 0.05)  return("*")
          "ns"
        }
        
        if (!has_inner) {
          # global test between the two outer groups
          sub <- df[!is.na(df$outer_group) & !is.na(df[[flux_col]]), ]
          if (nrow(sub) >= 2 && nlevels(sub$outer_group) == 2) {
            p_val <- tryCatch(
              {
                t.test(df[[flux_col]] ~ df$outer_group)$p.value
              },
              error = function(e) NA_real_
            )
          } else {
            p_val <- NA_real_
          }
          star_label <- pval_to_stars(p_val)
          
          # y position for the star: slightly above max
          y_max <- max(sub[[flux_col]], na.rm = TRUE)
          y_star <- y_max + 0.05 * diff(range(sub[[flux_col]], na.rm = TRUE))
          
          p_box <- ggplot(sub, aes(x = outer_group, y = .data[[flux_col]])) +
            geom_boxplot() +
            geom_jitter(width = 0.1, alpha = 0.6) +
            # add one star spanning both boxes
            geom_segment(aes(
              x = 1, xend = 2,
              y = y_star, yend = y_star
            )) +
            geom_text(
              aes(x = 1.5, y = y_star, label = star_label),
              vjust = -0.3
            ) +
            xlab(if (is.null(cont_name)) "Outer group" else cont_name) +
            ylab(flux_col)
          
        } else {
          # per-inner tests
          inner_levels <- sort(unique(df$inner))
          # compute a small data frame with positions and stars
          star_df <- do.call(rbind, lapply(inner_levels, function(lv) {
            sub <- df[df$inner == lv & !is.na(df$outer_group) & !is.na(df[[flux_col]]), ]
            if (nrow(sub) < 2 || nlevels(sub$outer_group) < 2) {
              return(NULL)
            }
            p_val <- tryCatch(
              {
                t.test(sub[[flux_col]] ~ sub$outer_group)$p.value
              },
              error = function(e) NA_real_
            )
            star_label <- pval_to_stars(p_val)
            y_max <- max(sub[[flux_col]], na.rm = TRUE)
            y_star <- y_max + 0.05 * diff(range(sub[[flux_col]], na.rm = TRUE))
            data.frame(
              inner = lv,
              x_start = 1,
              x_end   = 2,
              x_mid   = 1.5,
              y_star  = y_star,
              label   = star_label,
              stringsAsFactors = FALSE
            )
          }))
          
          p_box <- ggplot(df, aes(x = outer_group, y = .data[[flux_col]])) +
            geom_boxplot() +
            geom_jitter(width = 0.1, alpha = 0.6) +
            facet_grid(. ~ inner) +
            xlab(if (is.null(cont_name)) "Outer group" else cont_name) +
            ylab(flux_col)
          
          if (!is.null(star_df) && nrow(star_df) > 0) {
            p_box <- p_box +
              geom_segment(
                data = star_df,
                aes(x = x_start, xend = x_end, y = y_star, yend = y_star),
                inherit.aes = FALSE
              ) +
              geom_text(
                data = star_df,
                aes(x = x_mid, y = y_star, label = label),
                vjust = -0.3,
                inherit.aes = FALSE
              )
          }
        }
        
        # ---------- IF OUTER NOT CONTINUOUS-DERIVED: ONLY BOX PLOTS ----------
        
        if (!(outer_type %in% c(1, 2))) {
          return(p_box)
        }
        
        # ---------- RECONSTRUCT RAW CONTINUOUS OUTER VARIABLE ----------
        
        if (is.null(cont_name)) {
          return(p_box)
        }
        
        all_covars <- cbind(clinical(), gene())
        if (!cont_name %in% colnames(all_covars)) {
          return(p_box)
        }
        
        vals  <- all_covars[rownames(metaFinal), cont_name]
        x_num <- suppressWarnings(as.numeric(vals))
        if (all(is.na(x_num))) {
          return(p_box)
        }
        
        df$outer_cont <- x_num[match(rownames(df), rownames(metaFinal))]
        
        # ---------- CORRELATION PLOTS (outer_cont vs flux, r annotated) ----------
        
        flux_col <- input$boxplot
        
        if (!has_inner) {
          sub <- df[!is.na(df$outer_cont) & !is.na(df[[flux_col]]), ]
          if (nrow(sub) == 0) {
            return(p_box)
          }
          r_val <- suppressWarnings(cor(sub$outer_cont, sub[[flux_col]], use = "complete.obs"))
          p_corr <- ggplot(sub, aes(x = outer_cont, y = .data[[flux_col]])) +
            geom_point(alpha = 0.6) +
            geom_smooth(method = "lm", se = FALSE, color = "red") +
            xlab(cont_name) +
            ylab(flux_col) +
            ggtitle(
              paste0("Correlation across full cohort (r = ",
                     sprintf("%.2f", r_val), ")")
            )
        } else {
          inner_levels <- sort(unique(df$inner))
          corr_plots <- lapply(inner_levels, function(lv) {
            sub <- df[df$inner == lv & !is.na(df$outer_cont) & !is.na(df[[flux_col]]), ]
            if (nrow(sub) == 0) {
              return(ggplot() + ggtitle(paste("No data for", lv)))
            }
            r_val <- suppressWarnings(cor(sub$outer_cont, sub[[flux_col]], use = "complete.obs"))
            ggplot(sub, aes(x = outer_cont, y = .data[[flux_col]])) +
              geom_point(alpha = 0.6) +
              geom_smooth(method = "lm", se = FALSE, color = "red") +
              xlab(cont_name) +
              ylab(flux_col) +
              ggtitle(
                paste0(inner_name, " = ", lv, " (r = ", sprintf("%.2f", r_val), ")")
              )
          })
          p_corr <- cowplot::plot_grid(plotlist = corr_plots, ncol = 1)
        }
        
        cowplot::plot_grid(p_box, p_corr, ncol = 2, rel_widths = c(2, 1))
      })
      
      output$plot3<- renderPlot({
        plot3()
      })
      
      output$dlPlot3 <- downloadHandler(
        filename = function() {
          paste(input$boxplot,'_',label, '.pdf', sep='')
        },
        content = function(file) {
          ggsave(file, plot3(), width = 12, height = 8, dpi = 400, units = "in")
        }
      )
      
      #     }
      #   }
  })
  
  observeEvent(input$submitReportDefault,{
    if (length(errors_list()) > 0) {
      showModal(
        modalDialog(
          title = "Cannot run analysis",
          paste("Errors present:", paste(errors_list(), collapse = "; ")),
          easyClose = TRUE
        )
      )
      return()
    }
    
    metaFinal <- metaFinal_out()
    if (is.null(metaFinal) || nrow(metaFinal) == 0) {
      showModal(
        modalDialog(
          title = "No filtered cohort",
          "Please configure stratification, outer and inner comparisons, then run the pre-analysis summary before starting the analysis.",
          easyClose = TRUE
        )
      )
      return()
    }
    
    flux_full <- cbind(metaFinal, cohort_flux()[rownames(metaFinal), ])
    flux_comp <- flux_full
    
    outer <- colnames(metaFinal)[1]
    inner <- colnames(metaFinal)[2]
    
    outer_for_vec <- if (!is.null(outerVarFor())) outerVarFor()[[1]] else character(0)
    outer_rev_vec <- if (!is.null(outerVarRev())) outerVarRev()[[1]] else character(0)
    outer_all     <- c(outer_for_vec, outer_rev_vec)
    
    parse_outer <- function(x) {
      parts <- strsplit(x, " - ", fixed = TRUE)[[1]]
      if (length(parts) >= 2) list(var = parts[1], level = parts[2]) else list(var = x, level = x)
    }
    
    # Base comparison name (e.g. "FLT3.ITD")
    outer_comp_name <- if (length(outer_all) > 0) {
      parse_outer(outer_all[1])$var
    } else {
      "Outer comparison"
    }
    
    # Human-readable upper / lower labels
    upper_label <- if (length(outer_for_vec) > 0) parse_outer(outer_for_vec[1])$level else "upper"
    lower_label <- if (length(outer_rev_vec) > 0) parse_outer(outer_rev_vec[1])$level else "lower"
    
    #   if(input$outerType != 1){
    outer_select_1<-c("upper", "lower")
    inner_opts<-unique(metaFinal[,"inner"])
    metab_names<-c()
    data_filt<-flux_comp
    
    for (i in 1:length(inner_opts)){
      data_filt_var<-subset(data_filt, data_filt[,"inner"] == inner_opts[[i]])
      AUC = colAUC(data_filt_var[,-c(1:dim(metaFinal)[2])],factor(data_filt_var[,outer]),plotROC=FALSE)
      
      AUC = AUC[1,]
      #plot(density(AUC),main=paste0("AUC for ", inner_opts[i]," features"))
      logFC = log2(apply(subset(data_filt_var, subset = data_filt_var[,"outer"] == "upper")[,-c(1:dim(metaFinal)[2])],2,mean) / apply(subset(data_filt_var, subset = data_filt_var[,"outer"] == "lower")[,-c(1:dim(metaFinal)[2])],2,mean))
      short<-as.data.frame(logFC)
      t<-col_t_welch(as.matrix(subset(data_filt_var, subset = data_filt_var[,"outer"] == "upper")[,rownames(short)]),as.matrix(subset(data_filt_var, subset = data_filt_var[,"outer"] == "lower")[,rownames(short)]))
      final<-cbind(t$mean.y, t$mean.x, short, t$statistic, t$pvalue)
      colnames(final)<-c("var_1_mean","var_2_mean","log2FC","t_welch_statistic","pval")
      
      #final<-subset(final, subset = pval < 0.05)
      final$fdr = p.adjust(final[,"pval"],method="fdr")
      final<-subset(final, subset = pval < 0.05)
      metab_names<-c(metab_names, rownames(final))
    }
    metab_names<-unique(metab_names)
    
    test_list<-data.frame()
    for (i in 1:length(inner_opts)){
      names<-c(outer, metab_names)
      inner_comparison<-inner_opts[i]
      data_filt_var<-subset(data_filt, subset = data_filt[,"inner"] == inner_opts[i])
      data_filt_var<-data_filt_var[,names]
      t<-col_t_welch(as.matrix(subset(data_filt_var, subset = data_filt_var[,"outer"] == "upper")[,metab_names]),as.matrix(subset(data_filt_var, subset = data_filt_var[,"outer"] == "lower")[,metab_names]))
      final<-as.data.frame(cbind(metabolite=metab_names,inner_comparison=inner_comparison,as.numeric(t$mean.y), as.numeric(t$mean.x),as.numeric(t$statistic), as.numeric(t$pvalue)))
      colnames(final)<-c("metabolite","inner_comparison","var_1_mean","var_2_mean","t_welch_statistic","pval")
      final$sign_t_log_pval = ifelse(final$t_welch_statistic>0, -log10(as.numeric(final$pval)), log10(as.numeric(final$pval)))
      test_list<-as.data.frame(rbind(test_list,final))
    }
    
    metabs_tab <- merge(test_list, reactMeta, by = "metabolite", all.x = TRUE, all.y = FALSE)
    # If there was no inner comparison, ensure the column exists and is "All"
    
    metabs_tab$sig <- ifelse(
      metabs_tab$sign_t_log_pval > 1.301, "up",
      ifelse(metabs_tab$sign_t_log_pval < -1.301, "down", "ns")
    )
    
    # Re-level inner_comparison depending on innerType
    metabs_tab$inner_comparison <- switch(
      as.numeric(input$innerType),
      # 1: three-class continuous
      factor(metabs_tab$inner_comparison, levels = c("high", "baseline", "low")),
      # 2: two-class continuous
      factor(metabs_tab$inner_comparison, levels = c("high", "low")),
      # 3: binary
      factor(metabs_tab$inner_comparison, levels = unique(metabs_tab$inner_comparison)),
      # 4: categorical
      factor(metabs_tab$inner_comparison, levels = unique(metabs_tab$inner_comparison))
    )
    
    if (inner() == "") {
      metabs_tab$inner_comparison <- "All"
    }
    
    updateCheckboxGroupInput(session,"innerSelect", choices = unique(data_filt[,2]), selected = unique(data_filt[,2]))
    
    output$outerTitle<-renderText(outer)
    output$innerTitle<-renderText(inner)
    
    output$plot1Title <- renderText(
      paste0(
        "Significant Reactions Across All Subsystems - ",
        outer_comp_name,
        " (", upper_label, " v. ", lower_label, ")"
      )
    )
    
    metabs_fin <- reactive({
      req(exists("metabs_tab"), !is.null(metabs_tab))
      
      if (!"inner_comparison" %in% colnames(metabs_tab)) {
        return(metabs_tab)
      }
      
      if (nrow(metabs_tab) == 0) {
        return(metabs_tab)
      }
      
      inner_levels <- unique(metabs_tab[,"inner_comparison", drop = TRUE])
      
      if (length(inner_levels) == 1) {
        subset(metabs_tab, sig %in% input$outerSelect)
      } else {
        subset(
          metabs_tab,
          inner_comparison %in% input$innerSelect & sig %in% input$outerSelect
        )
      }
    })
    
    plot1 <- reactive({
      if(length(unique(metabs_fin()[,"inner_comparison"])) == 1){
        ggplot(metabs_fin(), aes(x=reorder(subsystem, -sign_t_log_pval), y=sign_t_log_pval, color=sig, fill = sig)) + geom_point(size = 1, position = position_dodge(width=0.5))+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)+xlab("subsystem")
      }else{
        ggplot(metabs_fin(), aes(x=reorder(subsystem, -sign_t_log_pval), y=sign_t_log_pval, color=inner_comparison, fill = inner_comparison)) + geom_point(size = 1, position = position_dodge(width=0.5))+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)+xlab("subsystem")
      }
      
    })
    
    output$plot1<- renderPlot({
      plot1()
    })
    
    plot1pdf <- reactive({
      ggplot(metabs_fin(), aes(x=as.numeric(t_welch_statistic), y=reorder(subsystem, -as.numeric(t_welch_statistic)), color=sig, fill = sig)) + geom_point(size = 1, position = position_dodge(width=0.5))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_vline(xintercept=1.3, size=.05)+geom_vline(xintercept=-1.3, size=.05)+ylab("subsystem")+xlab("log(pval) with sign of t-statistic")+ggtitle(paste0("Significant Reactions Across All Subsystems - (assay_1 v. assay_2)"))
    })
    
    output$reacts <- downloadHandler(
      filename = function() {
        paste(paste0(input$innerSelect, collapse = "_"),"_",inner, "_",paste0(input$outerSelect, collapse = "_"), outer, '_filtered_reactions.rds', sep='')
      },
      content = function(file) {
        saveRDS(list(metabs_fin(),metaFinal,data.frame(AUC=input$customAUC,Log2FC=input$customFC,FDR=input$customFDR),"cohort"),file)
      }
    )
    
    # output$reacts_meta <- downloadHandler(
    #   filename = function() {
    #     paste(paste0(input$innerSelect, collapse = "_"),"_",inner, "_",paste0(input$outerSelect, collapse = "_"), outer, '_meta_data.csv', sep='')
    #   },
    #   content = function(file) {
    #     write.csv(metaFinal,file)
    #   }
    # )
    
    systems<-reactive({
      if(length(unique(metabs_fin()[,"inner_comparison"])) == 1){
        ggbld <- ggplot_build(ggplot(metabs_fin(), aes(x=reorder(subsystem, -sign_t_log_pval), y=sign_t_log_pval, color=sig, fill = sig)) + geom_point(size = 1, position = position_dodge(width=0.5))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)+xlab("subsystem"))
        ggbld$layout$panel_params[[1]]$x$limits
      }else{
        ggbld <- ggplot_build(ggplot(metabs_fin(), aes(x=reorder(subsystem, -sign_t_log_pval), y=sign_t_log_pval, color=inner_comparison, fill = inner_comparison)) + geom_point(size = 1, position = position_dodge(width=0.5))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)+xlab("subsystem"))
        ggbld$layout$panel_params[[1]]$x$limits
      }
      
    })
    
    
    output$dlPlot1 <- downloadHandler(
      filename = function() {
        paste('all_subsystems.pdf', sep='')
      },
      content = function(file) {
        withProgress(message = 'Generating PDF', value = 0, {
          incProgress(0.3, detail = "Preparing data")
          
          pdf(file, width = 11, height = 8.5)  # Standard letter size
          
          # Plot or write your content to the PDF here
          # For example:
          plot(plot1pdf())
          # Add more plots or text as needed
          
          incProgress(0.3, detail = "Creating PDF")
          
          # Close the PDF device
          dev.off()
          
          incProgress(0.4, detail = "Finalizing")
        })
      },
      contentType = "application/pdf"
    )
    
    observeEvent(input$clickBar, {
      groupId <- round(input$clickBar$x)
      output$plot2Title <- renderText(
        paste0(
          "Significant Reactions ", outer_comp_name, " - ", systems()[groupId],
          " (", upper_label, " v. ", lower_label, ")"
        )
      )
      plot2 <- reactive({
        if(length(unique(metabs_fin()[,"inner_comparison"])) == 1){
          ggplot(subset(metabs_fin(), subset = subsystem == systems()[groupId]), aes(x=inner_comparison, y=sign_t_log_pval, color=sig, fill = sig)) + geom_point(size = 2)+geom_line(aes(group = metabolite, x=inner_comparison),color="grey")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)
        }else{
          ggplot(subset(metabs_fin(), subset = subsystem == systems()[groupId]), aes(x=inner_comparison, y=sign_t_log_pval, color=inner_comparison, fill = inner_comparison)) + geom_point(size = 2)+geom_line(aes(group = metabolite, x=inner_comparison),color="grey")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)
        }
        
      })
      output$plot2<- renderPlot({
        plot2()
      })
      plot2pdf <- reactive({
        if(length(unique(metabs_fin()[,"inner_comparison"])) == 1){
          ggplot(subset(metabs_fin(), subset = subsystem == systems()[groupId]), aes(x=inner_comparison, y=sign_t_log_pval, color=sig, fill = sig)) + geom_point(size = 2)+geom_line(aes(group = metabolite, x=inner_comparison),color="grey")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)+ggtitle(paste0("Significant Reactions ",outer, " - ",systems[groupId], " (upper v. lower)"))
        }else{
          ggplot(subset(metabs_fin(), subset = subsystem == systems()[groupId]), aes(x=inner_comparison, y=sign_t_log_pval, color=inner_comparison, fill = inner_comparison)) + geom_point(size = 2)+geom_line(aes(group = metabolite, x=inner_comparison),color="grey")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)+ggtitle(paste0("Significant Reactions ",outer, " - ",systems[groupId], " (upper v. lower)"))
        }
        
      })
      output$dlPlot2 <- downloadHandler(
        filename = function() {
          paste(str_replace(systems()[groupId]," ","_"),'_',outer,'.pdf', sep='')
        },
        content = function(file) {
          ggsave(file, plot2pdf(), width = 8, height = 12, dpi = 400, units = "in")
        }
      )
      output$info<-renderReactable({
        metabs_fin<-metabs_fin()
        metabs_fin$var_1_mean<-round(as.numeric(metabs_fin$var_1_mean), digits=4)
        metabs_fin$var_2_mean<-round(as.numeric(metabs_fin$var_2_mean), digits=4)
        metabs_fin$sign_t_log_pval<-round(as.numeric(metabs_fin$sign_t_log_pval), digits=4)
        reactable(brushedPoints(subset(metabs_fin, subset = subsystem == systems()[groupId]), input$plot_brush2, xvar = "inner_comparison", yvar="sign_t_log_pval"))
      })
      data <- reactive({
        as.data.frame(brushedPoints(subset(metabs_fin, subset = subsystem == systems()[groupId]), input$plot_brush2, xvar = "inner_comparison", yvar="sign_t_log_pval"))
      })
      output$data <- downloadHandler(
        filename = function() {
          paste(str_replace(systems()[groupId]," ","_"),'_selected_reactions_',outer, '.csv', sep='')
        },
        content = function(file) {
          write.csv(data(),file)
        }
      )
      
      observe({
        updateSelectizeInput(session, "boxplot",
                             choices = brushedPoints(subset(metabs_fin(), subset = subsystem == systems()[groupId]), input$plot_brush2, xvar = "inner_comparison", yvar="sign_t_log_pval")$metabolite
        )})
    })
    
    plot3 <- reactive({
      if (input$boxplot == "") {
        return(ggplot())
      }
      
      metaFinal <- metaFinal_out()
      req(metaFinal)
      
      outer_var <- colnames(metaFinal)[1]
      inner_var <- colnames(metaFinal)[2]
      
      outer_type <- as.numeric(input$outerType)
      has_inner  <- !(is.null(inner()) || inner() == "" || type_inner() == 0)
      inner_name <- if (has_inner) inner() else "None"
      
      # Base data
      df <- data_filt[, c(colnames(metaFinal), input$boxplot), drop = FALSE]
      colnames(df)[1] <- "outer_group"
      colnames(df)[2] <- "inner"
      
      # ---------- Parse outer selections for labels ----------
      
      outer_for_vec <- if (!is.null(outerVarFor())) outerVarFor()[[1]] else character(0)
      outer_rev_vec <- if (!is.null(outerVarRev())) outerVarRev()[[1]] else character(0)
      outer_all     <- c(outer_for_vec, outer_rev_vec)
      
      parse_outer <- function(x) {
        parts <- strsplit(x, " - ", fixed = TRUE)[[1]]
        if (length(parts) >= 2) list(var = parts[1], level = parts[2]) else list(var = x, level = x)
      }
      
      # Default labels
      cont_name       <- NULL
      outer_level_map <- NULL
      outer_levels    <- sort(unique(df$outer_group))
      
      if (length(outer_all) > 0) {
        parsed <- lapply(outer_all, parse_outer)
        cont_name <- parsed[[1]]$var                    # e.g. "FLT3.ITD"
        for_levels <- vapply(outer_for_vec, function(z) parse_outer(z)$level, character(1))
        rev_levels <- vapply(outer_rev_vec, function(z) parse_outer(z)$level, character(1))
        outer_levels <- unique(c(for_levels, rev_levels))
        
        # Map original outer_group values to pretty labels:
        # assume outer_group encodes upper/lower in your metaFinal,
        # but here we simply relabel based on order of levels we derived
        if (length(outer_levels) >= 2) {
          outer_level_map <- setNames(outer_levels, c("upper", "lower")[seq_along(outer_levels)])
        }
      }
      
      # Apply outer label mapping to df$outer_group (if available)
      df$outer_group <- factor(as.character(df$outer_group))
      if (!is.null(outer_level_map)) {
        # Map factor levels "upper"/"lower" (or whatever you used) to user-friendly labels
        df$outer_group <- factor(
          outer_level_map[as.character(df$outer_group)],
          levels = outer_levels
        )
      }
      
      # ---------- BOX PLOTS (original behavior, now with pretty x labels) ----------
      flux_col <- input$boxplot
      
      pval_to_stars <- function(p) {
        if (is.na(p)) return("ns")
        if (p < 0.001) return("***")
        if (p < 0.01)  return("**")
        if (p < 0.05)  return("*")
        "ns"
      }
      
      if (!has_inner) {
        # global test between the two outer groups
        sub <- df[!is.na(df$outer_group) & !is.na(df[[flux_col]]), ]
        if (nrow(sub) >= 2 && nlevels(sub$outer_group) == 2) {
          p_val <- tryCatch(
            {
              t.test(df[[flux_col]] ~ df$outer_group)$p.value
            },
            error = function(e) NA_real_
          )
        } else {
          p_val <- NA_real_
        }
        star_label <- pval_to_stars(p_val)
        
        # y position for the star: slightly above max
        y_max <- max(sub[[flux_col]], na.rm = TRUE)
        y_star <- y_max + 0.05 * diff(range(sub[[flux_col]], na.rm = TRUE))
        
        p_box <- ggplot(sub, aes(x = outer_group, y = .data[[flux_col]])) +
          geom_boxplot() +
          geom_jitter(width = 0.1, alpha = 0.6) +
          # add one star spanning both boxes
          geom_segment(aes(
            x = 1, xend = 2,
            y = y_star, yend = y_star
          )) +
          geom_text(
            aes(x = 1.5, y = y_star, label = star_label),
            vjust = -0.3
          ) +
          xlab(if (is.null(cont_name)) "Outer group" else cont_name) +
          ylab(flux_col)
        
      } else {
        # per-inner tests
        inner_levels <- sort(unique(df$inner))
        # compute a small data frame with positions and stars
        star_df <- do.call(rbind, lapply(inner_levels, function(lv) {
          sub <- df[df$inner == lv & !is.na(df$outer_group) & !is.na(df[[flux_col]]), ]
          if (nrow(sub) < 2 || nlevels(sub$outer_group) < 2) {
            return(NULL)
          }
          p_val <- tryCatch(
            {
              t.test(sub[[flux_col]] ~ sub$outer_group)$p.value
            },
            error = function(e) NA_real_
          )
          star_label <- pval_to_stars(p_val)
          y_max <- max(sub[[flux_col]], na.rm = TRUE)
          y_star <- y_max + 0.05 * diff(range(sub[[flux_col]], na.rm = TRUE))
          data.frame(
            inner = lv,
            x_start = 1,
            x_end   = 2,
            x_mid   = 1.5,
            y_star  = y_star,
            label   = star_label,
            stringsAsFactors = FALSE
          )
        }))
        
        p_box <- ggplot(df, aes(x = outer_group, y = .data[[flux_col]])) +
          geom_boxplot() +
          geom_jitter(width = 0.1, alpha = 0.6) +
          facet_grid(. ~ inner) +
          xlab(if (is.null(cont_name)) "Outer group" else cont_name) +
          ylab(flux_col)
        
        if (!is.null(star_df) && nrow(star_df) > 0) {
          p_box <- p_box +
            geom_segment(
              data = star_df,
              aes(x = x_start, xend = x_end, y = y_star, yend = y_star),
              inherit.aes = FALSE
            ) +
            geom_text(
              data = star_df,
              aes(x = x_mid, y = y_star, label = label),
              vjust = -0.3,
              inherit.aes = FALSE
            )
        }
      }
      
      # ---------- IF OUTER NOT CONTINUOUS-DERIVED: ONLY BOX PLOTS ----------
      
      if (!(outer_type %in% c(1, 2))) {
        return(p_box)
      }
      
      # ---------- RECONSTRUCT RAW CONTINUOUS OUTER VARIABLE ----------
      
      if (is.null(cont_name)) {
        return(p_box)
      }
      
      all_covars <- cbind(clinical(), gene())
      if (!cont_name %in% colnames(all_covars)) {
        return(p_box)
      }
      
      vals  <- all_covars[rownames(metaFinal), cont_name]
      x_num <- suppressWarnings(as.numeric(vals))
      if (all(is.na(x_num))) {
        return(p_box)
      }
      
      df$outer_cont <- x_num[match(rownames(df), rownames(metaFinal))]
      
      # ---------- CORRELATION PLOTS (outer_cont vs flux, r annotated) ----------
      
      flux_col <- input$boxplot
      
      if (!has_inner) {
        sub <- df[!is.na(df$outer_cont) & !is.na(df[[flux_col]]), ]
        if (nrow(sub) == 0) {
          return(p_box)
        }
        r_val <- suppressWarnings(cor(sub$outer_cont, sub[[flux_col]], use = "complete.obs"))
        p_corr <- ggplot(sub, aes(x = outer_cont, y = .data[[flux_col]])) +
          geom_point(alpha = 0.6) +
          geom_smooth(method = "lm", se = FALSE, color = "red") +
          xlab(cont_name) +
          ylab(flux_col) +
          ggtitle(
            paste0("Correlation across full cohort (r = ",
                   sprintf("%.2f", r_val), ")")
          )
      } else {
        inner_levels <- sort(unique(df$inner))
        corr_plots <- lapply(inner_levels, function(lv) {
          sub <- df[df$inner == lv & !is.na(df$outer_cont) & !is.na(df[[flux_col]]), ]
          if (nrow(sub) == 0) {
            return(ggplot() + ggtitle(paste("No data for", lv)))
          }
          r_val <- suppressWarnings(cor(sub$outer_cont, sub[[flux_col]], use = "complete.obs"))
          ggplot(sub, aes(x = outer_cont, y = .data[[flux_col]])) +
            geom_point(alpha = 0.6) +
            geom_smooth(method = "lm", se = FALSE, color = "red") +
            xlab(cont_name) +
            ylab(flux_col) +
            ggtitle(
              paste0(inner_name, " = ", lv, " (r = ", sprintf("%.2f", r_val), ")")
            )
        })
        p_corr <- cowplot::plot_grid(plotlist = corr_plots, ncol = 1)
      }
      
      cowplot::plot_grid(p_box, p_corr, ncol = 2, rel_widths = c(2, 1))
    })
    
    output$plot3<- renderPlot({
      plot3()
    })
    
    output$dlPlot3 <- downloadHandler(
      filename = function() {
        paste(input$boxplot,'_',label, '.pdf', sep='')
      },
      content = function(file) {
        ggsave(file, plot3(), width = 12, height = 8, dpi = 400, units = "in")
      }
    )
    
    #     }
    #   }
  })
}

shinyApp(ui,server)
