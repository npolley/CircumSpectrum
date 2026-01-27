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
library(readr)
library(data.table)
library(caret)
library(dplyr)
library(shinycssloaders)
library(SHAPforxgboost)
library(xgboost)
library(fastshap)
library(shapviz)
library(umap)
library(inflection)

cohort_registry <- readr::read_csv("../../../../data/RADAR_xCheck_experimental_assay/xCheck_experiment_registry.csv")

dataset_registry <- cohort_registry %>%
  mutate(
    id        = folder_name,  # internal key used everywhere
    label     = paste(experiment_name, "-", treatment),
    flux_file = paste0("../../../../data/RADAR_xCheck_experimental_assay/",folder_name,"/",folder_name, "_flux",".csv"),
    meta_file = paste0("../../../../data/RADAR_xCheck_experimental_assay/",folder_name,"/",folder_name, "_meta",".csv")
  )

#Read in flux table, metadata, and organism reaction information

load_flux <- function(path) {
  as.data.frame(as.matrix(data.table::fread(path), rownames = 1))
}

load_meta <- function(path) {
  read.csv(path, header = TRUE, row.names = 1)
}

flux_list <- setNames(
  lapply(dataset_registry$flux_file, load_flux),
  dataset_registry$id
)

meta_list <- setNames(
  lapply(dataset_registry$meta_file, load_meta),
  dataset_registry$id
)
#data<-data[rownames(meta),]

reactMeta<-read.csv(paste0("../../../../data/RADAR_xCheck_experimental_assay/human_reaction_meta.csv"))
colnames(reactMeta)[1]<-"metabolite"


#Assign categories to be compared

header <- dashboardHeader(
  title = "RADAR | xCheck",
  titleWidth = 400 
)

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Fingerprint Designer", tabName = "main", icon = icon("dashboard"))
  )
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "main",
  h2("Metabolic Fingerprint Designer for Experiments"),
  fluidRow(column(width=12,box(width=12,title="Select Experiment and Stratification Attributes",solidHeader = TRUE, status="primary", 
                               box(width = 6, selectInput(
                                 label   = "Select Experiment",
                                 inputId = "experiment",
                                 choices = setNames(dataset_registry$id, dataset_registry$label),
                                 selected = dataset_registry$id[1]
                               )),
                               box(width = 6, virtualSelectInput(label="Select Stratification", inputId = "strat", choices = c(), showValueAsTags = TRUE, search = TRUE, multiple = TRUE)),
                               fluidRow(column(width=12,actionButton("submitStrat", "Submit")))),tags$div(id = 'stratification')),
                               ),

  fluidRow(column(width=8,box(title=textOutput("plot1Title"),
                              width = 12, solidHeader = TRUE, status="primary",withLoader(plotOutput("plot1", click = "clickBar", height = 500), type = "html", loader = "dnaspin"),
                              downloadLink("dlPlot1", "Download Plot as PDF | "), downloadLink("dlSvg", "Download Image as SVG | "), downloadLink("reacts", "Download Fingerprint Data")),
  ), 
  column(width=4,box(title=textOutput("plot2Title"),
                     width = 12, solidHeader = TRUE, status="primary", withLoader(plotOutput("plot2", brush = "plot_brush2", height = 500), type = "html", loader = "dnaspin"),
                     downloadLink("dlPlot2", "Download Plot as PDF")))),
  fluidRow(column(width=8,box(title=textOutput("outerTitle"), solidHeader = TRUE,status="primary",checkboxGroupInput("outerSelect", label = NULL, choices = c("Significant (+)"="up","Not Signficant"="ns","Significant (-)"="down"), selected=c("up","ns","down"))),box(title=textOutput("innerTitle"), solidHeader = TRUE, status="primary",checkboxGroupInput("innerSelect", label = NULL, choices = c("high","baseline","low"))))),
  fluidRow(column(width=12,box(width=12,selectizeInput("boxplot", "Select Metabolic Flux to Analyze", choices = c(""))))),
  fluidRow(column(width=12,box(width=12,reactableOutput("info")))),
  fluidRow(column(width=12,box(width=12,withLoader(plotOutput("plot3", height = 500), type = "html", loader = "dnaspin"),actionLink("invert","Invert Boxplot | "), downloadLink("dlPlot3", "Download Plot as PDF"))))
  )
  )
)

ui<-dashboardPage(
  skin = "purple",
  header = header,
  sidebar = sidebar, 
  body = body
)

server <-function(input,output,session){
  options(shiny.maxRequestSize=100*1024^2)
  
  experiment_name      <- reactiveVal(NULL)
  experiment_flux      <- reactiveVal(NULL)
  experiment_meta      <- reactiveVal(NULL)
  experiment_meta_filt <- reactiveVal(NULL)
  
  outer     <- reactiveVal(NULL)
  inner     <- reactiveVal(NULL)
  top       <- reactiveVal(NULL)
  bottom    <- reactiveVal(NULL)
  innerVars <- reactiveVal(NULL)
  
  error     <- reactiveVal(NULL)
  
  customAUC <- reactiveVal(NULL)
  customFC  <- reactiveVal(NULL)
  customFDR <- reactiveVal(NULL)

  rds_data <- reactiveVal(NULL)
  load_success <- reactiveVal(FALSE)
  fingerprint_flux<-reactiveVal(FALSE)
  react_list<-reactiveVal(vector("list", 100))
  
  observeEvent(input$load_rds, {
    req(input$rds_upload)
    
    tryCatch({
      loaded_data <- readRDS(input$rds_upload$datapath)
      rds_data(loaded_data)
      output$rds_load_status <- renderText("Fingerprint file successfully loaded.")
      flux_list<-list("scotland"=metformin_flux,"ustinova_healthy"=metformin_healthy_flux, "ustinova_diabetes"=metformin_diabetes_flux, "kulkarni_bulk"=kulkarni_flux, "kulkarni_2_bulk"=kulkarni_2_flux, "beataml_response"=beataml_flux, "glioblastoma_hr"=glio_flux,"sabatier_GILT"=sabatier_gilt_flux, "TUH84"=TUH84_flux, "beataml"=beataml_flux,"cervical"=cervical_flux,"breast"=breast_flux,"melanoma"=melanoma_flux,"PDAC"=PDAC_flux)
      fingerprint_flux(flux_list[[rds_data()$All$assay]])
      updateVirtualSelect('selected_features_subsystem', choices = list(rds_data()$All$fingerprint_final$coefnames), selected = rds_data()$All$fingerprint_final$coefnames)
      load_success(TRUE)
    }, error = function(e) {
      output$rds_load_status <- renderText(paste("Error loading fingerprint file:", e$message))
      load_success(FALSE)
    })
    
    
  })
  
  shap_data <- reactive({
    req(rds_data())
    
    classifier_final <- rds_data()
    vals_subsystem <- classifier_final$All$self_score_primary
    shap_values <- classifier_final$All$shaps
    shap_values<-data.frame(colMeans(abs(shap_values)))
    colnames(shap_values)[1]<-"shap"
    shap_values<-data.frame(shap_values, spacer=0)
    shap_values <- shap_values[order(shap_values[["shap"]], decreasing = TRUE), ]
    filt_shap <- subset(shap_values, shap > 0)
    
    return(list(
      shap_long_frag = shap_values,
      filt_shap = filt_shap
    ))
  })
  
  output$shap_plot_box <- renderUI({
    req(load_success())
    box(
      title = "SHAP Plot (by Subsystem)",
      width = NULL,
      solidHeader = TRUE,
      status = "primary",
      shinycssloaders::withSpinner(plotOutput("shap_plot", height = "655px")),
      downloadButton("download_shap_plot", "Download SHAP Plot")
    )
  })
  
  observe({
    classifier_final <- rds_data()
    selected_features <- classifier_final$All$fingerprint_final$coefnames
    
    updateSelectInput(session, "feature_selector",
                      choices = selected_features,
                      selected = selected_features)
  })
  
  output$umap_plot_box <- renderUI({
    req(load_success())
    classifier_final <- rds_data()
    selected_features <- input$selected_features_subsystem
    
    box(
      title = "UMAP Plot (by Subsystem)",
      width = NULL,
      solidHeader = TRUE,
      status = "primary",
      shinycssloaders::withSpinner(plotOutput("umap_plot", height = "655px")),
      downloadButton("download_umap_plot", "Download UMAP Plot")
    )
  })
  
  output$shap_density <- renderUI({
    req(input$selected_react, rds_data())
    box(
      title = paste0(input$selected_react," Density Plot"),
      width = NULL,
      solidHeader = TRUE,
      status = "primary",
      shinycssloaders::withSpinner(plotOutput("density_plot", height = "600px")),
      downloadButton("download_density_plot", "Download Density Plot")
    )
  })
  
  output$shap_boxplot <- renderUI({
    req(input$selected_react, rds_data())
    box(
      title = paste0(input$selected_react," Boxplot"),
      width = NULL,
      solidHeader = TRUE,
      status = "primary",
      shinycssloaders::withSpinner(plotOutput("shap_box_plot", height = "600px")),
      downloadButton("download_boxplot", "Download Boxplot")
    )
  })
  
  selected_reaction_data <- reactive({
    req(input$selected_react, rds_data())
    
    flux<-flux_list[[rds_data()$All$assay]]
    meta<-meta_list[[rds_data()$All$assay]]
    
    flux<-flux[rownames(meta),]
    
    # Filter the data for the selected reaction
    reaction_data <- na.omit(data.frame(outer=meta[,"IDH2_MUT"],flux[, input$selected_react]))
    
    return(reaction_data)
  })
  
  # Render the density plot
  output$density_plot <- renderPlot({
    req(selected_reaction_data())
    
    data <- selected_reaction_data()
    View(data)
    data[,"outer"]<-factor(data[,"outer"])
    data[,2]<-as.numeric(data[,2])
    colnames(data)[2]<-"flux"
    
    ggplot(data, aes(x = flux, fill = outer, alpha=0.7)) +
      geom_density() + 
      labs(title = paste("Density Plot for", input$selected_react),
           x = "Value",
           y = "Density") +
      theme(plot.title = element_text(hjust = 0.5))
  })
  
  output$shap_box_plot <- renderPlot({
    req(selected_reaction_data())
    
    data <- selected_reaction_data()
    data[,"outer"]<-factor(data[,"outer"])
    data[,2]<-as.numeric(data[,2])
    colnames(data)[2]<-"flux"
    
    ggplot(data, aes(x = flux, y=outer, fill = outer, alpha=0.7)) +
      geom_boxplot() + 
      geom_jitter() + 
      labs(title = paste("Boxplot for", input$selected_react),
           x = "Value",
           y = "Density") +
      theme(plot.title = element_text(hjust = 0.5))
  })
  
  # Download handler for the density plot
  output$download_density_plot <- downloadHandler(
    filename = function() {
      paste0("density_plot_", input$selected_react, "_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      data <- selected_reaction_data()
      data[,"outer"]<-factor(data[,"outer"])
      data[,2]<-as.numeric(data[,2])
      colnames(data)[2]<-"flux"
      
      plot <- ggplot(data, aes(x = flux, fill = outer)) +
        geom_density(alpha=0.7) + 
        labs(title = paste("Density Plot for", input$selected_react),
             x = "Value",
             y = "Density") +
        theme(plot.title = element_text(hjust = 0.5))
      
      ggsave(file, plot = plot, device = "pdf", width = 10, height = 6)
    }
  )
  
  output$download_box_plot <- downloadHandler(
    filename = function() {
      paste0("boxplot_", input$selected_react, "_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      data <- selected_reaction_data()
      data[,"outer"]<-factor(data[,"outer"])
      data[,2]<-as.numeric(data[,2])
      colnames(data)[2]<-"flux"
      
      plot <- ggplot(data, aes(x = flux, fill = outer)) +
        geom_boxplot(alpha=0.7) + geom_jitter()+
        labs(title = paste("Boxplot for", input$selected_react),
             x = "Value",
             y = "Density") +
        theme(plot.title = element_text(hjust = 0.5))
      
      ggsave(file, plot = plot, device = "pdf", width = 10, height = 6)
    }
  )
  
  output$react_info<-renderReactable({
    req(input$selected_react)
    reactable(subset(reactMeta,metabolite==input$selected_react))
  })
  
  output$shap_dropdown_box_2 <- renderUI({
    
    # req(react_list())
    # 
    # models_list <- rds_data()$All$fingerprints_primary
    # vals<-rds_data()$All$self_score_primary
    # 
    # model_n <- which(names(vals) == input$selected_shap_var)-1
    # rownames <- models_list[[model_n]]$feature_names
    
    box(
      title = "Select Metabolic Reaction",
      width = 12,
      solidHeader = TRUE,
      status = "primary",
      selectInput("selected_react", "Select Reaction", 
                  choices = react_list(), 
                  selected = react_list()[1])
    )
  })
  
  output$shap_dropdown_box <- renderUI({
    req(load_success(), shap_data())
    
    rownames <- rownames(shap_data()$filt_shap)
    
    box(
      title = "Select Metabolic Subsystem",
      width = 12,
      solidHeader = TRUE,
      status = "primary",
      selectInput("selected_shap_var", "Select Subsystem", 
                  choices = rownames, 
                  selected = rownames[1])
    )
  })
  
  subsystem_shap_plot_reactive <- reactive({
    req(input$selected_shap_var, rds_data())
    
    classifier_final <- rds_data()
    auc <- rds_data()$All$auc_primary
    models_list <- rds_data()$All$fingerprints_primary
    vals<-rds_data()$All$self_score_primary
    
    # Find the model number corresponding to the selected subsystem
    model_n <- which(names(vals) == input$selected_shap_var)-1
    shap_long_frag <- shap.prep(xgb_model = models_list[[model_n]]$finalModel, 
                                X_train = as.matrix(classifier_final$All$train_primary[, models_list[[model_n]]$coefnames]))
    filt_shap_1 <- subset(shap_long_frag, mean_value > 0)
    
    react_list(unique(filt_shap_1$variable))
    
    shap_plot <- shap.plot.summary.wrap1(models_list[[model_n]]$finalModel, 
                                         X = as.matrix(classifier_final$All$train_primary[, models_list[[model_n]]$coefnames]), 
                                         top_n = length(unique(filt_shap_1$variable)))
    return(shap_plot)
  })
  
  output$subsystem_shap_plot <- renderPlot({
    subsystem_shap_plot_reactive()
  })
  
  output$subsystem_shap_plot_box <- renderUI({
    req(input$selected_shap_var)
    if(length(react_list())<15){
      box_height<-300
    }else{
      box_height<-length(react_list())*20
    }
    
    box(
      title = paste("SHAP Plot for", input$selected_shap_var),
      width = NULL,
      solidHeader = TRUE,
      status = "primary",
      shinycssloaders::withSpinner(plotOutput("subsystem_shap_plot", height = paste0(box_height,"px"))),
      downloadButton("download_subsystem_shap_plot", "Download Subsystem SHAP Plot")
    )
  })
  
  output$download_subsystem_shap_plot <- downloadHandler(
    filename = function() {
      paste("subsystem_shap_plot_", input$selected_shap_var, "_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = subsystem_shap_plot_reactive(), device = "pdf", width = 12, height = 8)
    }
  )
  
  umap_plot_reactive <- reactive({
    req(rds_data())
    
    classifier_final <- rds_data()
    vals_subsystem <- rds_data()$All$self_score_primary
    
    selected_features<-input$selected_features_subsystem
    
    if(length(selected_features)>1){
      umap_result<-umap(vals_subsystem[,selected_features])
      umap_df <- data.frame(
        UMAP1 = umap_result$layout[,1],
        UMAP2 = umap_result$layout[,2],
        Class = factor(vals_subsystem$outer)
      )
      umap_plot<-ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Class)) +
        geom_point() +
        theme_minimal() +
        labs(title = "UMAP Plot for Kulkarni Bulk",
             x = "UMAP1",
             y = "UMAP2")
      
      return(umap_plot)
    }else{
      no_plot<-plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
      text(0.5, 0.5, "UMAP Plot Unavailable for Individual Subsystems\nPlease Consult SHAP Plot for Individual Subsystem Effect.", cex = 1.2, col = "blue")
      return(no_plot)
    }
    
  })
  
  
  output$umap_plot <- renderPlot({
    umap_plot_reactive()
  })
  
  output$download_umap_plot <- downloadHandler(
    filename = function() {
      paste("umap_plot_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = umap_plot_reactive(), device = "pdf", width = 12, height = 8)
    }
  )
  
  shap_plot_reactive <- reactive({
    req(rds_data())
    
    classifier_final <- rds_data()
    vals_subsystem <- rds_data()$All$self_score_primary
    
    pred_wrapper <- function(object, newdata) {
      # Replace this with the appropriate prediction method for your model
      predict(object, newdata = newdata, type = "prob")[,1]
    }
    
    shap_values <- classifier_final$All$shaps
    sv<-shapviz(shap_values, X = vals_subsystem[,classifier_final$All$fingerprint_final$coefnames])
    shap_plot<-sv_importance(sv, kind = "bee", max_display = 30)
    return(shap_plot)
  })
  
  output$shap_plot <- renderPlot({
    shap_plot_reactive()
  })
  
  output$download_shap_plot <- downloadHandler(
    filename = function() {
      paste("shap_plot_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = shap_plot_reactive(), device = "pdf", width = 12, height = 8)
    }
  )
  
  output$performance_plot_box <- renderUI({
    req(load_success())
    fluidRow(
      box(
        title = "Performance Plot",
        width = 12,
        solidHeader = TRUE,
        status = "primary",
        shinycssloaders::withSpinner(plotOutput("performance_plot", height = "600px", width = "1500px")),
        downloadButton("download_performance_plot", "Download Performance Plot")
      )
    )
  })
  
  performance_plot_reactive <- reactive({
    req(rds_data())
    
    # Assuming the RDS file contains the necessary objects
    classifier_final <- rds_data()
    vals_subsystem <- rds_data()$All$self_score_primary
    perf_subsystem <- rds_data()$All$auc_primary
    auc_value <- rds_data()$All$auc_final
    assay_name <- rds_data()$All$assay
    
    shap_values <- classifier_final$All$shaps
    shap_values<-data.frame(colMeans(abs(shap_values)))
    colnames(shap_values)[1]<-"shap"
    shap_values<-data.frame(shap_values, spacer=0)
    shap_values <- shap_values[order(shap_values[["shap"]], decreasing = TRUE), ]
    
    imp<-varImp(classifier_final$All$fingerprint_final, scale = TRUE)[[1]][,"Class1"]
    mean_imp<-imp/sum(imp)
    
    shaps <- data.frame(subsystem=rownames(shap_values), score=shap_values$shap)
    rownames(shaps) <- rownames(shap_values)
    stats_fin <- data.frame(index=rownames(perf_subsystem), subsystem=perf_subsystem$subsystem, loocv_test_auc=as.numeric(perf_subsystem$cv_auc))
    stats_fin <- data.frame(stats_fin, shaps=shaps[colnames(vals_subsystem)[as.numeric(stats_fin$index)+1],"score"])
    stats_fin$shaps <- ifelse(is.na(stats_fin$shaps), min(stats_fin$shaps, na.rm = TRUE), stats_fin$shaps)
    stats_fin <- rbind(data.frame(index=0, subsystem="FINAL MODEL", loocv_test_auc=auc_value, shaps=NA), stats_fin)
    stats_fin$subsystem <- factor(stats_fin$subsystem, levels = stats_fin$subsystem[order(stats_fin$loocv_test_auc, decreasing = TRUE)])
    length_x <- dim(stats_fin)[1]-10
    
    ggplot(stats_fin, aes(x=subsystem, y = loocv_test_auc, alpha = shaps)) +
      geom_bar(stat="identity", fill = "steelblue") +
      coord_cartesian(ylim=c(0.4,1)) +
      theme(axis.text.x = element_text(size = 10,angle = 90, vjust = 0.5, hjust = 1)) +
      geom_text(aes(y = 0.91, label = "Excellent"), x = length_x, hjust = 0, size = 3) +
      geom_segment(aes(x=2,y = 0.9, xend=dim(stats_fin)[1]-3), linewidth = 0.05, linetype = "dashed" ) +
      geom_text(aes(y = 0.81, label = "Very Good"), x = length_x, hjust = 0, size = 3) +
      geom_segment(aes(x=2,y = 0.8, xend=dim(stats_fin)[1]-3), linewidth = 0.05, linetype = "dashed" ) +
      geom_text(aes(y = 0.71, label = "Good"), x = length_x, hjust = 0, size = 3) +
      geom_segment(aes(x=2,y = 0.7, xend=dim(stats_fin)[1]-3), linewidth = 0.05, linetype = "dashed" ) +
      geom_text(aes(y = 0.61, label = "Adequate"), x = length_x, hjust = 0, size = 3) +
      geom_segment(aes(x=2,y = 0.6, xend=dim(stats_fin)[1]-3), linewidth = 0.05, linetype = "dashed" ) +
      geom_text(aes(y = 0.5, label = "Irrelevant"), x = length_x, hjust = 0, size = 3) +
      labs(title = paste0("Metabolic Fingerprint Performance - ",assay_name," (internal prediction)"), 
           y = paste0("Test AUC - LOOCV (",length(rds_data()$All$self_score_final)," Folds)"), 
           alpha="SHAP Importance in Final Model") + 
      theme(plot.margin = unit(c(1, 1, 1, 4), "lines"))
  })
  
  performance_plot_reactive_pdf <- reactive({
    req(rds_data())
    
    # Assuming the RDS file contains the necessary objects
    classifier_final <- rds_data()
    vals_subsystem <- rds_data()$All$self_score_primary
    perf_subsystem <- rds_data()$All$auc_primary
    auc_value <- rds_data()$All$auc_final
    assay_name <- rds_data()$All$assay
    
    shap_values <- classifier_final$All$shaps
    shap_values<-data.frame(colMeans(abs(shap_values)))
    colnames(shap_values)[1]<-"shap"
    shap_values<-data.frame(shap_values, spacer=0)
    shap_values <- shap_values[order(shap_values[["shap"]], decreasing = TRUE), ]
    
    imp<-varImp(classifier_final$All$fingerprint_final, scale = TRUE)[[1]][,"Class1"]
    mean_imp<-imp/sum(imp)
    
    shaps <- data.frame(subsystem=rownames(shap_values), score=shap_values$shap)
    rownames(shaps) <- rownames(shap_values)
    stats_fin <- data.frame(index=rownames(perf_subsystem), subsystem=perf_subsystem$subsystem, loocv_test_auc=as.numeric(perf_subsystem$loocv_test_auc))
    stats_fin <- data.frame(stats_fin, shaps=shaps[colnames(vals_subsystem)[as.numeric(stats_fin$index)+1],"score"])
    stats_fin$shaps <- ifelse(is.na(stats_fin$shaps), min(stats_fin$shaps, na.rm = TRUE), stats_fin$shaps)
    stats_fin <- rbind(data.frame(index=0, subsystem="FINAL MODEL", loocv_test_auc=auc_value, shaps=NA), stats_fin)
    stats_fin$subsystem <- factor(stats_fin$subsystem, levels = stats_fin$subsystem[order(stats_fin$loocv_test_auc, decreasing = FALSE)])
    
    inflection<-findiplist(stats_fin$loocv_test_auc, 1:length(stats_fin$loocv_test_auc), 0)[2,3]
    stats_fin<-subset(stats_fin, loocv_test_auc > inflection)
    
    length_x <- dim(stats_fin)[1]-10
    
    if(length(rds_data()$All$self_score_final) > 100){
      folds<-100
    }else{
      folds<-length(rds_data()$All$self_score_final)
    }
    
    ggplot(stats_fin, aes(x = loocv_test_auc, y = subsystem, alpha = shaps)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_cartesian(xlim = c(0.4, 1)) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 6, hjust = 1),
        plot.margin = unit(c(1, 1, 1, 1), "lines"),
        aspect.ratio = 16/9  # Adjust this ratio to make the plot less wide
      ) +
      labs(
        title = paste0("Metabolic Fingerprint Performance - ", assay_name, " (internal prediction)"),
        x = paste0("Test AUC - LOOCV (", folds, " Folds)"),
        alpha = "SHAP Importance in Final Model"
      ) +
      scale_x_continuous(breaks = seq(0.4, 1, by = 0.1))+
      geom_vline(xintercept = c(0.9, 0.8, 0.7, 0.6), linetype = "dashed", color = "gray50", size = 0.5)
  })
  
  output$performance_plot <- renderPlot({
    performance_plot_reactive()
  })
  
  output$download_performance_plot <- downloadHandler(
    filename = function() {
      paste("performance_plot_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = performance_plot_reactive_pdf(), device = "pdf", width = 18, height = 8)
    }
  )
  
  observeEvent(input$experiment, {
    id  <- input$experiment
    row <- dataset_registry[dataset_registry$id == id, ]
    
    experiment_flux(flux_list[[id]])
    experiment_meta(meta_list[[id]])
    experiment_name(row$label)
    
    # reset downstream state
    experiment_meta_filt(NULL)
    outer(NULL); inner(NULL)
    top(NULL); bottom(NULL)
    innerVars(NULL)
    error(NULL)
    
    # clear inputs
    updateVirtualSelect("strat", choices = list(), selected = NULL)
    updateSelectizeInput(session, "outerType", choices = character(0), selected = NULL)
    updateVirtualSelect("outerVarTop",    choices = list(), selected = NULL)
    updateVirtualSelect("outerVarBottom", choices = list(), selected = NULL)
    updateSelectizeInput(session, "innerType", choices = character(0), selected = NULL)
    updateVirtualSelect("innerVar", choices = list(), selected = NULL)
    
    # repopulate stratification for the new experiment
    result <- lapply(experiment_meta(), function(col) {
      list(unique(col))
    })
    updateVirtualSelect(
      "strat",
      choices  = result,
      selected = unique(unlist(experiment_meta()))
    )
  })

  observeEvent(input$submitStrat, {
    # reset downstream state
    outer(NULL); inner(NULL)
    top(NULL); bottom(NULL)
    innerVars(NULL)
    error(NULL)
    
    # filter meta
    matching_rows <- apply(experiment_meta(), 1, function(row) all(row %in% input$strat))
    filtered_df   <- experiment_meta()[matching_rows, ]
    experiment_meta_filt(filtered_df)
    
    # first time only: create outer box + anchor
    if (input$submitStrat == 1) {
      insertUI(
        selector  = "#stratification",
        where     = "afterEnd",
        immediate = TRUE,
        ui = fluidRow(
          column(
            width = 6,
            box(
              width = 12,
              title = "Select Outer Comparison",
              solidHeader = TRUE,
              status = "primary",
              box(
                width = 12,
                selectizeInput("outerType", "Select Outer Comparison", choices = c())
              ),
              box(
                width = 6,
                virtualSelectInput(
                  "outerVarTop", "Select Top Variable",
                  choices = c(), showValueAsTags = TRUE,
                  search = TRUE, multiple = TRUE
                )
              ),
              box(
                width = 6,
                virtualSelectInput(
                  "outerVarBottom", "Select Bottom Variable",
                  choices = c(), showValueAsTags = TRUE,
                  search = TRUE, multiple = TRUE
                )
              ),
              fluidRow(
                column(width = 12, actionButton("submitOuter", "Submit"))
              )
            )
          ),
          tags$div(id = "outer")  # anchor for inner box
        )
      )
    }
    
    # always: update outerType choices
    updateSelectizeInput(
      session,
      "outerType",
      choices  = colnames(experiment_meta_filt())[colnames(experiment_meta_filt()) != "none"],
      selected = NULL
    )
  })
  
  observeEvent(input$outerType, {
    req(input$outerType)
    updateVirtualSelect('outerVarTop', choices = unique(experiment_meta_filt()[,input$outerType]))
    updateVirtualSelect('outerVarBottom', choices = unique(experiment_meta_filt()[,input$outerType]))
  })

  observeEvent(input$submitOuter, {
    req(input$submitOuter)
    
    # reset inner-level state
    inner(NULL)
    innerVars(NULL)
    error(NULL)
    
    outer(input$outerType)
    top(input$outerVarTop)
    bottom(input$outerVarBottom)
    
    # first time only: create inner box under #outer
    if (input$submitOuter == 1) {
      insertUI(
        selector  = "#outer",
        where     = "beforeEnd",
        immediate = TRUE,
        ui = column(
          width = 6,
          box(
            width = 12,
            title = "Select Inner Comparison",
            solidHeader = TRUE,
            status = "primary",
            box(
              width = 12,
              selectizeInput("innerType", "Select Inner Comparison", choices = c())
            ),
            box(
              width = 12,
              virtualSelectInput(
                "innerVar", "Select Variable(s)",
                choices = c(), multiple = TRUE,
                showValueAsTags = TRUE, search = TRUE
              )
            ),
            fluidRow(
              column(width = 12, actionButton("submitInner", "Submit"))
            )
          ),
          tags$div(id = "inner")
        )
      )
    }
    
    # always: update innerType choices for new outer selection
    updateSelectizeInput(
      session,
      "innerType",
      choices  = colnames(experiment_meta_filt())[colnames(experiment_meta_filt()) != outer()],
      selected = NULL
    )
  })

  observeEvent(list(input$submitOuter,input$innerType),{
    req(input$submitOuter)
    req(input$innerType)
    
    updateVirtualSelect(
      "innerVar",
      choices  = unique(experiment_meta_filt()[, input$innerType]),
      selected = unique(experiment_meta_filt()[, input$innerType])
    )
  })

  observeEvent(input$submitInner, {
    req(input$submitInner)
    
    inner(input$innerType)
    innerVars(input$innerVar)
    error(NULL)
    
    # first time only: create summary box
    if (input$submitInner == 1) {
      insertUI(
        selector  = "#outer",
        where     = "afterEnd",
        immediate = TRUE,
        ui = fluidRow(
          column(
            width = 12,
            box(
              width = 12,
              title = "RADAR | xCheck Pre-Analysis Summary",
              solidHeader = TRUE,
              status = "primary",
              box(width = 6, verbatimTextOutput("reportText")),
              box(
                width = 6,
                noUiSliderInput("customAUC", "AUC:",   min = 0.6, max = 1.0, value = 0.65, step = 0.01),
                noUiSliderInput("customFC",  "Log2FC:",min = 0,   max = 5,   value = 1,    step = 0.01),
                noUiSliderInput("customFDR", "FDR:",   min = 0,   max = 1,   value = 0.05, step = 0.01)
              ),
              fluidRow(
                column(width = 12, actionButton("submitReport",        "Begin Analysis (set parameters)")),
                column(width = 12, actionButton("submitReportDefault", "Begin Analysis (default parameters)"))
              )
            )
          )
        )
      )
    }
  })

  observeEvent(list(input$submitStrat, input$submitOuter,input$submitInner), {
    req(input$submitInner)
    warnings <- c()
    error(c())
    
    if (is.null(top()))      { error(c(error(), "At least one top outer variable must be chosen")) }
    if (is.null(bottom()))   { error(c(error(), "At least one bottom outer variable must be chosen")) }
    if (is.null(innerVars())){ error(c(error(), "At least one inner variable must be chosen")) }
    
    if (is.null(warnings)) {
      warnings <- "NONE"
    } else {
      warnings <- paste0(warnings, collapse = " - ")
    }
    
    if (is.null(error())) {
      error("NONE")
    } else {
      error(paste0(error(), collapse = " - "))
    }
    
    output$reportText <- renderText({
      paste(
        paste0("Experiment: ", experiment_name(), "\n"),
        paste0("Outer Category: ", outer(), "\n"),
        paste0("Outer Variables: ", paste(top(), collapse = ', '), " (top) - ",
               paste(bottom(), collapse = ', '), " (bottom)", "\n"),
        paste0("Inner Category: ", inner(), "\n"),
        paste0("Inner Variables: ", paste(innerVars(), collapse = ', '), "\n"),
        paste0("Minimum AUC: ", input$customAUC, "\n"),
        paste0("Minimum Log2FC: ", input$customFC, "\n"),
        paste0("Maximum FDR: ", input$customFDR, "\n"),
        paste0("Warnings: ", warnings, "\n"),
        paste0("Errors: ", error(), "\n"),
        paste0("Report Generated at: ", Sys.time())
      )
    })
  })
  
    observeEvent(input$submitReport,{
      customAUC(input$customAUC)
      customFC(input$customFC)
      customFDR(input$customFDR)
      
      if(error()=="NONE"){
        metaFinal<-experiment_meta_filt()
          outer<-outer()
          inner<-inner()
            quant1_out<-top()
            quant2_out<-bottom()
            metaFinal[,outer]<-ifelse(metaFinal[,outer] %in% quant1_out, paste(quant1_out,collapse = ', '),paste(quant2_out,collapse = ', '))
            inner_opts<-innerVars()
            metab_names<-c()
            flux_full<-cbind(metaFinal,experiment_flux()[rownames(metaFinal),])
            data_filt<-flux_full

            metaFinal_condensed<-metaFinal[,c(outer, inner)]
            colnames(metaFinal_condensed)<-c("outer","inner")
            metaFinal_condensed$outer<-ifelse(metaFinal_condensed$outer == quant1_out, "upper", "lower")

            for (i in 1:length(inner_opts)){
              data_filt_var<-subset(data_filt, data_filt[,inner] == inner_opts[[i]])
              AUC = colAUC(data_filt_var[,-c(1:dim(metaFinal)[2])],factor(data_filt_var[,outer]),plotROC=FALSE)
              AUC = AUC[1,]
              plot(density(AUC),main=paste0("AUC for ", inner_opts[i]," features"))
              
              logFC = log2(apply(subset(data_filt_var, subset = data_filt_var[,outer] == paste(quant1_out,collapse = ', '))[,-c(1:dim(metaFinal)[2])],2,mean) / apply(subset(data_filt_var, subset = data_filt_var[,outer] == paste(quant2_out,collapse = ', '))[,-c(1:dim(metaFinal)[2])],2,mean))
              logFC[AUC<customAUC()]=0
              short<-as.data.frame(logFC[which(abs(logFC)>customFC())])
              t<-col_t_welch(as.matrix(subset(data_filt_var, subset = data_filt_var[,outer] == paste(quant1_out,collapse = ', '))[,rownames(short)]),as.matrix(subset(data_filt_var, subset = data_filt_var[,outer] == paste(quant2_out,collapse = ', '))[,rownames(short)]))
              final<-cbind(t$mean.y, t$mean.x, short, t$statistic, t$pvalue)
              colnames(final)<-c("var_1_mean","var_2_mean","log2FC","t_welch_statistic","pval")
              #final<-subset(final, subset = pval < 0.05)
              final$fdr = p.adjust(final[,"pval"],method="fdr")
              final<-subset(final, subset = fdr < customFDR())
              metab_names<-c(metab_names, rownames(final))
            }

            metab_names<-unique(metab_names)
            test_list<-data.frame()
            
            for (i in 1:length(inner_opts)){
              names<-c(outer, metab_names)
              inner_comparison<-inner_opts[i]
              data_filt_var<-subset(data_filt, subset = data_filt[,inner] == inner_opts[i])
              data_filt_var<-data_filt_var[,names]
              t<-col_t_welch(as.matrix(subset(data_filt_var, subset = data_filt_var[,outer] == paste(quant1_out,collapse = ', '))[,metab_names]),as.matrix(subset(data_filt_var, subset = data_filt_var[,outer] == paste(quant2_out,collapse = ', '))[,metab_names]))
              final<-as.data.frame(cbind(metabolite=metab_names,inner_comparison=inner_comparison,as.numeric(t$mean.y), as.numeric(t$mean.x),as.numeric(t$statistic), as.numeric(t$pvalue)))
              colnames(final)<-c("metabolite","inner_comparison","var_1_mean","var_2_mean","t_welch_statistic","pval")
              final$sign_t_log_pval = ifelse(final$t_welch_statistic>0, -log10(as.numeric(final$pval)), log10(as.numeric(final$pval)))
              test_list<-as.data.frame(rbind(test_list,final))

            }
           
            metabs_tab<-merge(test_list, reactMeta, by = "metabolite", all.x=TRUE, all.y=FALSE)
            metabs_tab$sig<-factor(ifelse(metabs_tab$sign_t_log_pval > 1.301, "up", ifelse(metabs_tab$sign_t_log_pval < -1.301, "down", "ns")), levels = c("up","ns","down"))
            
            updateCheckboxGroupInput(session,"innerSelect", choices = innerVars(), selected = innerVars())
            
            metabs_tab$inner_comparison<-factor(metabs_tab$inner_comparison, levels = unique(metabs_tab$inner_comparison))
            
            output$outerTitle<-renderText(outer)
            output$innerTitle<-renderText(inner)
            
            output$plot1Title<-renderText(paste0("Significant Reactions Across All Subsystems - ",outer, " (",quant1_out, " v. ",quant2_out,")"))
            
            metabs_fin<-reactive({
              subset(metabs_tab, inner_comparison %in% input$innerSelect & sig %in% input$outerSelect)
            })

            plot1 <- reactive({
              ggplot(metabs_fin(), aes(x=reorder(subsystem, -sign_t_log_pval), y=sign_t_log_pval, color=sig, fill = sig)) + geom_point(size = 1, position = position_dodge(width=0.5))+
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)+xlab("subsystem")
            })
            
            output$plot1<- renderPlot({
              plot1()
            })
            
            plot1pdf <- reactive({
              ggplot(metabs_fin(), aes(x=reorder(subsystem, -sign_t_log_pval), y=sign_t_log_pval, color=inner_comparison, fill = inner_comparison)) + geom_point(size = 1, position = position_dodge(width=0.5))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)+xlab("subsystem")+ggtitle(paste0("Significant Reactions Across All Subsystems - ",outer, " (",quant1_out, " v. ",quant2_out,")"))
            })
            
            output$reacts <- downloadHandler(
              filename = function() {
                paste(paste0(input$innerSelect, collapse = "_"),"_",inner, "_",paste0(input$outerSelect, collapse = "_"), outer, '_filtered_reactions_.rds', sep='')
              },
              content = function(file) {
                saveRDS(list(metabs_fin(),metaFinal_condensed,data.frame(AUC=input$customAUC,Log2FC=input$customFC,FDR=input$customFDR),"experiment"),file)
              }
            )
            
            # output$reacts_meta <- downloadHandler(
            #   filename = function() {
            #     paste(paste0(input$innerSelect, collapse = "_"),"_",inner, "_",paste0(input$outerSelect, collapse = "_"), outer, '_meta_data.csv', sep='')
            #   },
            #   content = function(file) {
            #     write.csv(metaFinal_condensed,file)
            #   }
            # )
            
            systems<-reactive({
              ggbld <- ggplot_build(ggplot(metabs_fin(), aes(x=reorder(subsystem, -sign_t_log_pval), y=sign_t_log_pval, color=inner_comparison, fill = inner_comparison)) + geom_point(size = 1, position = position_dodge(width=0.5))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)+xlab("subsystem"))
              ggbld$layout$panel_params[[1]]$x$limits
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
            
            output$dlSvg <- downloadHandler(
              filename = function() {
                paste('all_subsystems_',label,'.svg', sep='')
              },
              content = function(file) {
                ggsave(file, plot1Svg(), width = 24, height = 12, dpi = 400, units = "in")
              }
            )
            
            observeEvent(input$clickBar, {
              groupId <- round(input$clickBar$x)
              output$plot2Title<-renderText(paste0("Significant Reactions ",outer, " - ",systems()[groupId], " (",paste(quant1_out,collapse = ', '), " v. ",paste(quant2_out,collapse = ', '),")"))
              plot2 <- reactive({
                ggplot(subset(metabs_fin(), subset = subsystem == systems()[groupId]), aes(x=inner_comparison, y=sign_t_log_pval, color=inner_comparison, fill = inner_comparison)) + geom_point(size = 2)+geom_line(aes(group = metabolite, x=inner_comparison),color="grey")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)
              })
              output$plot2<- renderPlot({
                plot2()
              })
              plot2pdf <- reactive({
                ggplot(subset(metabs_fin(), subset = subsystem == systems()[groupId]), aes(x=inner_comparison, y=sign_t_log_pval, color=inner_comparison, fill = inner_comparison)) + geom_point(size = 2)+geom_line(aes(group = metabolite, x=inner_comparison),color="grey")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)+ggtitle(paste0("Significant Reactions ",outer, " - ",systems[groupId], " (",paste(quant1_out,collapse = ', '), " v. ",paste(quant2_out,collapse = ', '),")"))
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
                metabs_fin<-metabs_fin()[,-c(5,6)]
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
              if(input$boxplot != ""){

                  forBoxplot<-data_filt[,c(input$outerType,input$innerType,input$boxplot)]
                  colnames(forBoxplot)[1]<-outer
                  colnames(forBoxplot)[2]<-"inner"
                  forBoxplot$inner<-factor(forBoxplot$inner, unique(forBoxplot$inner))
                  if(input$invert %% 2 == 0){
                    ggplot(forBoxplot, aes_string(x=outer, y=input$boxplot))+geom_boxplot()+facet_grid(.~inner)+geom_jitter()+ylab(input$boxplot)
                  }else{
                    ggplot(forBoxplot, aes_string(x="inner", y=input$boxplot))+geom_boxplot()+facet_grid(outer)+geom_jitter()+ylab(input$boxplot)+xlab(inner)
                  }
                  
                
              }else{
                ggplot()
              }
              
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
        
      }
    
  

    
  })
    observeEvent(input$submitReportDefault,{
      if(error()=="NONE"){
        metaFinal<-experiment_meta_filt()
        outer<-outer()
        inner<-inner()
        quant1_out<-top()
        quant2_out<-bottom()
        metaFinal[,outer]<-ifelse(metaFinal[,outer] %in% quant1_out, paste(quant1_out,collapse = ', '),paste(quant2_out,collapse = ', '))
        inner_opts<-innerVars()
        metab_names<-c()
        flux_full<-cbind(metaFinal,experiment_flux()[rownames(metaFinal),])
        data_filt<-flux_full
        
        metaFinal_condensed<-metaFinal[,c(outer, inner)]
        colnames(metaFinal_condensed)<-c("outer","inner")
        metaFinal_condensed$outer<-ifelse(metaFinal_condensed$outer == quant1_out, "upper", "lower")
        
        for (i in 1:length(inner_opts)){
          data_filt_var<-subset(data_filt, data_filt[,inner] == inner_opts[[i]])
          AUC = colAUC(data_filt_var[,-c(1:dim(metaFinal)[2])],factor(data_filt_var[,outer]),plotROC=FALSE)
          AUC = AUC[1,]
          
          logFC = log2(apply(subset(data_filt_var, subset = data_filt_var[,outer] == paste(quant1_out,collapse = ', '))[,-c(1:dim(metaFinal)[2])],2,mean) / apply(subset(data_filt_var, subset = data_filt_var[,outer] == paste(quant2_out,collapse = ', '))[,-c(1:dim(metaFinal)[2])],2,mean))
          short<-as.data.frame(logFC)
          t<-col_t_welch(as.matrix(subset(data_filt_var, subset = data_filt_var[,outer] == paste(quant1_out,collapse = ', '))[,rownames(short)]),as.matrix(subset(data_filt_var, subset = data_filt_var[,outer] == paste(quant2_out,collapse = ', '))[,rownames(short)]))
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
          data_filt_var<-subset(data_filt, subset = data_filt[,inner] == inner_opts[i])
          data_filt_var<-data_filt_var[,names]
          t<-col_t_welch(as.matrix(subset(data_filt_var, subset = data_filt_var[,outer] == paste(quant1_out,collapse = ', '))[,metab_names]),as.matrix(subset(data_filt_var, subset = data_filt_var[,outer] == paste(quant2_out,collapse = ', '))[,metab_names]))
          final<-as.data.frame(cbind(metabolite=metab_names,inner_comparison=inner_comparison,as.numeric(t$mean.y), as.numeric(t$mean.x),as.numeric(t$statistic), as.numeric(t$pvalue)))
          colnames(final)<-c("metabolite","inner_comparison","var_1_mean","var_2_mean","t_welch_statistic","pval")
          final$sign_t_log_pval = ifelse(final$t_welch_statistic>0, -log10(as.numeric(final$pval)), log10(as.numeric(final$pval)))
          test_list<-as.data.frame(rbind(test_list,final))
          
        }
        
        metabs_tab<-merge(test_list, reactMeta, by = "metabolite", all.x=TRUE, all.y=FALSE)
        metabs_tab$sig<-factor(ifelse(metabs_tab$sign_t_log_pval > 1.301, "up", ifelse(metabs_tab$sign_t_log_pval < -1.301, "down", "ns")), levels = c("up","ns","down"))
        
        updateCheckboxGroupInput(session,"innerSelect", choices = innerVars(), selected = innerVars())
        
        metabs_tab$inner_comparison<-factor(metabs_tab$inner_comparison, levels = unique(metabs_tab$inner_comparison))
        
        output$outerTitle<-renderText(outer)
        output$innerTitle<-renderText(inner)
        
        output$plot1Title<-renderText(paste0("Significant Reactions Across All Subsystems - ",outer, " (",quant1_out, " v. ",quant2_out,")"))
        
        metabs_fin<-reactive({
          subset(metabs_tab, inner_comparison %in% input$innerSelect & sig %in% input$outerSelect)
        })
        
        plot1 <- reactive({
          p<-ggplot(metabs_fin(), aes(x=reorder(subsystem, -sign_t_log_pval), y=sign_t_log_pval, color=sig, fill = sig)) + geom_point(size = 1, position = position_dodge(width=0.5))+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)+xlab("subsystem")
          p
        })
        
        output$plot1<- renderPlot({
          plot1()
        })
        
        plot1pdf <- reactive({
          ggplot(metabs_fin(), aes(x=as.numeric(t_welch_statistic), y=reorder(subsystem, -as.numeric(t_welch_statistic)), color=sig, fill = sig)) + geom_point(size = 1, position = position_dodge(width=0.5))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_vline(xintercept=1.3, size=.05)+geom_vline(xintercept=-1.3, size=.05)+ylab("subsystem")+xlab("log(pval) with sign of t-statistic")+ggtitle(paste0("Significant Reactions Across All Subsystems - (assay_1 v. assay_2)"))
        })
        
        output$reacts <- downloadHandler(
          filename = function() {
            paste(paste0(input$innerSelect, collapse = "_"),"_",inner, "_",paste0(input$outerSelect, collapse = "_"), outer, '_filtered_reactions_.rds', sep='')
          },
          content = function(file) {
            saveRDS(list(metabs_fin(),metaFinal_condensed,data.frame(AUC=input$customAUC,Log2FC=input$customFC,FDR=input$customFDR),"experiment"),file)
          }
        )
        
        # output$reacts_meta <- downloadHandler(
        #   filename = function() {
        #     paste(paste0(input$innerSelect, collapse = "_"),"_",inner, "_",paste0(input$outerSelect, collapse = "_"), outer, '_meta_data.csv', sep='')
        #   },
        #   content = function(file) {
        #     write.csv(metaFinal_condensed,file)
        #   }
        # )
        
        systems<-reactive({
          ggbld <- ggplot_build(ggplot(metabs_fin(), aes(x=reorder(subsystem, -sign_t_log_pval), y=sign_t_log_pval, color=inner_comparison, fill = inner_comparison)) + geom_point(size = 1, position = position_dodge(width=0.5))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)+xlab("subsystem"))
          ggbld$layout$panel_params[[1]]$x$limits
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
        
        output$dlSvg <- downloadHandler(
          filename = function() {
            paste('all_subsystems_',label,'.svg', sep='')
          },
          content = function(file) {
            ggsave(file, plot1Svg(), width = 24, height = 12, dpi = 400, units = "in")
          }
        )
        
        observeEvent(input$clickBar, {
          groupId <- round(input$clickBar$x)
          output$plot2Title<-renderText(paste0("Significant Reactions ",outer, " - ",systems()[groupId], " (",paste(quant1_out,collapse = ', '), " v. ",paste(quant2_out,collapse = ', '),")"))
          plot2 <- reactive({
            ggplot(subset(metabs_fin(), subset = subsystem == systems()[groupId]), aes(x=inner_comparison, y=sign_t_log_pval, color=inner_comparison, fill = inner_comparison)) + geom_point(size = 2)+geom_line(aes(group = metabolite, x=inner_comparison),color="grey")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)
          })
          output$plot2<- renderPlot({
            plot2()
          })
          plot2pdf <- reactive({
            ggplot(subset(metabs_fin(), subset = subsystem == systems()[groupId]), aes(x=inner_comparison, y=sign_t_log_pval, color=inner_comparison, fill = inner_comparison)) + geom_point(size = 2)+geom_line(aes(group = metabolite, x=inner_comparison),color="grey")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)+ggtitle(paste0("Significant Reactions ",outer, " - ",systems[groupId], " (",paste(quant1_out,collapse = ', '), " v. ",paste(quant2_out,collapse = ', '),")"))
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
            metabs_fin<-metabs_fin()[,-c(5,6)]
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
          if(input$boxplot != ""){
            
            forBoxplot<-data_filt[,c(input$outerType,input$innerType,input$boxplot)]
            colnames(forBoxplot)[1]<-outer
            colnames(forBoxplot)[2]<-"inner"
            forBoxplot$inner<-factor(forBoxplot$inner, unique(forBoxplot$inner))
            if(input$invert %% 2 == 0){
              ggplot(forBoxplot, aes_string(x=outer, y=input$boxplot))+geom_boxplot()+facet_grid(.~inner)+geom_jitter()+ylab(input$boxplot)
            }else{
              ggplot(forBoxplot, aes_string(x="inner", y=input$boxplot))+geom_boxplot()+facet_grid(outer)+geom_jitter()+ylab(input$boxplot)+xlab(inner)
            }
            
            
          }else{
            ggplot()
          }
          
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
        
      }
      
      
      
      
    })
  
  # observeEvent(list(input$assay,input$comp),{
  #     if(as.numeric(input$assay)==1){
  # 
  #     }
  #   if(as.numeric(input$assay)==2){
  #     label<-"Combination PIK3CA/PIK3CG (threefold)"
  #     metabs_fin<-switch(as.numeric(input$comp),combination_1_high_v_low, combination_1_high_v_baseline, combination_1_baseline_v_low)
  #     data_filt<-flux_full
  #     data_filt<-data_filt[complete.cases(data_filt[,"combination_1"]),]
  #     label2<-switch(as.numeric(input$comp), "(High v. Low)", "(High v. Baseline)","(Baseline v. Low)")
  #     data_filt$combination_1<-factor(data_filt$combination_1, levels=c("low","baseline","high"))
  #     metabs_fin$variable<-factor(metabs_fin$variable, levels=c("PIK3CA-PIK3CG"))
  #   }
  #   if(as.numeric(input$assay)==3){
  #     label<-"Combination PIK3CA/PIK3CG (by median)"
  #     metabs_fin<-combination_2_high_v_low
  #     data_filt<-flux_full
  #     data_filt<-data_filt[complete.cases(data_filt[,"combination_2"]),]
  #     label2<-"(High v. Low)"
  #     data_filt$combination_2<-factor(data_filt$combination_2, levels=c("low","high"))
  #     metabs_fin$variable<-factor(metabs_fin$variable, levels=c("PIK3CA-PIK3CG"))
  #   }
  #   
  #   output$plot1Title<-renderText(paste0("Significant Reactions Across All Subsystems - ",label, " ", label2))
  #   

  #   
  #   observeEvent(input$clickBar, {
  #     groupId <- round(input$clickBar$x)
  #     output$plot2Title<-renderText(paste0("Significant Reactions ",label, " - ",systems[groupId], " ",label2))
  #     plot2 <- reactive({
  #       ggplot(subset(metabs_fin, subset = subsystem == systems[groupId]), aes(x=variable, y=sign_t_log_pval, color=variable, fill = variable)) + geom_point(size = 2)+geom_line(aes(group = metabolite, x=variable),color="grey")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)
  #     })
  #     output$plot2<- renderPlot({
  #       plot2()
  #     })
  #     plot2pdf <- reactive({
  #       ggplot(subset(metabs_fin, subset = subsystem == systems[groupId]), aes(x=variable, y=sign_t_log_pval, color=variable, fill = variable)) + geom_point(size = 2)+geom_line(aes(group = metabolite, x=variable),color="grey")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)+ggtitle(paste0("Significant Reactions ",label, " - ",systems[groupId], " ",label2))
  #     })
  #     output$dlPlot2 <- downloadHandler(
  #       filename = function() {
  #         paste(str_replace(systems[groupId]," ","_"),'_',label,'.pdf', sep='')
  #       },
  #       content = function(file) {
  #         ggsave(file, plot2pdf(), width = 8, height = 12, dpi = 400, units = "in")
  #       }
  #     )
  #     output$info<-renderReactable({
  #       metabs_fin<-metabs_fin[,-c(5,6)]
  #       metabs_fin$var_1_mean<-round(metabs_fin$var_1_mean, digits=4)
  #       metabs_fin$var_2_mean<-round(metabs_fin$var_2_mean, digits=4)
  #       metabs_fin$sign_t_log_pval<-round(metabs_fin$sign_t_log_pval, digits=4)
  #       reactable(brushedPoints(subset(metabs_fin, subset = subsystem == systems[groupId]), input$plot_brush2, xvar = "variable", yvar="sign_t_log_pval"))
  #     })
  #     data <- reactive({
  #       as.data.frame(brushedPoints(subset(metabs_fin, subset = subsystem == systems[groupId]), input$plot_brush2, xvar = "variable", yvar="sign_t_log_pval"))
  #     })
  #     output$data <- downloadHandler(
  #       filename = function() {
  #         paste(str_replace(systems[groupId]," ","_"),'_selected_reactions_',label, '.csv', sep='')
  #       },
  #       content = function(file) {
  #         write.csv(data(),file)
  #       }
  #     )
  #     
  #     observe({
  #       updateSelectizeInput(session, "boxplot",
  #                            choices = brushedPoints(subset(metabs_fin, subset = subsystem == systems[groupId]), input$plot_brush2, xvar = "variable", yvar="sign_t_log_pval")$metabolite
  #       )})
  #   })
  # 
  #   plot3 <- reactive({
  #     if(input$boxplot != ""){
  #       if(as.numeric(input$assay)==1){
  #         plot1<-ggplot(data_filt[,c(colnames(meta),input$boxplot)], aes_string(x="PIK3CA", y=input$boxplot))+geom_boxplot()+geom_jitter()+ylab(input$boxplot)
  #         plot2<-ggplot(data_filt[,c(colnames(meta),input$boxplot)], aes_string(x="PIK3CG", y=input$boxplot))+geom_boxplot()+geom_jitter()+ylab(input$boxplot)
  #         plot_grid(plot1, plot2, labels = "AUTO")
  #       }
  #       else if(as.numeric(input$assay)==2){
  #         ggplot(data_filt[,c(colnames(meta),input$boxplot)], aes_string(x="combination_1", y=input$boxplot))+geom_boxplot()+geom_jitter()+ylab(input$boxplot)
  #       }
  #       else{
  #         ggplot(data_filt[,c(colnames(meta),input$boxplot)], aes_string(x="combination_2", y=input$boxplot))+geom_boxplot()+geom_jitter()+ylab(input$boxplot)
  #       }
  #     }else{
  #       ggplot()
  #     }
  # 
  #   })
  #   
  #   output$plot3<- renderPlot({
  #     plot3()
  #   })
  #   
  #   output$dlPlot3 <- downloadHandler(
  #     filename = function() {
  #       paste(input$boxplot,'_',label, '.pdf', sep='')
  #     },
  #     content = function(file) {
  #       ggsave(file, plot3(), width = 12, height = 8, dpi = 400, units = "in")
  #     }
  #   )
  # })
}

shinyApp(ui,server)
