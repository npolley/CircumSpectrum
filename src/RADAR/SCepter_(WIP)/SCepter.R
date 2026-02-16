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

bosc<-readRDS("bosc_TUH07_TUH69.rds")
bosc_flux<-as.data.frame(as.matrix(fread("bosc_flux.csv"), rownames = 1))
sample_choices_bosc<-c("Patient-TUH07","Patient-TUH69")
treatment_choices_bosc<-c("Diagnosis","Cytarabine", "Venetoclax", "Cytarabine-Venetoclax")

sabatier<-readRDS("sabatier_TUH110_TUH93_TUH84.rds")
sabatier_flux<-as.data.frame(as.matrix(fread("sabatier_flux.csv"), rownames = 1))
sample_choices_sabatier<-c("TUH110","TUH93", "TUH84")
treatment_choices_sabatier<-c("CTL_D0","GILT_D7", "GILT_REL_D37", "GILT_D14", "APR_D7", "GILT_APR_D7")

boet<-readRDS("boet_TUH109.rds")
boet_flux<-as.data.frame(as.matrix(fread("boet_flux.csv"), rownames = 1))
sample_choices_boet<-c("Patient-TUH109")
treatment_choices_boet<-c("CTL_D0","ARAC_D1", "ARAC_D2", "ARAC_D5", "ARAC_D8", "ARAC_RLP_D29")

kulkarni<-readRDS("kulkarni_sc.rds")
kulkarni_flux<-as.data.frame(as.matrix(fread("kulkarni_sc_flux.csv"), rownames = 1))
sample_choices_kulkarni<-c("Muscle")
treatment_choices_kulkarni<-c("CTL","MET")

condition_choices<-c("Acute Myeloid Leukemia (AML)","Pancreatic Ductal Adenocarcinoma (PDAC)", "Metformin - black6")
aml_studies<-c("Cytarabine and Venetoclax - Bosc et al. 2021","Gilteritinib FLT3/CEBPA - Sabatier et al. 2022","CD39 - Boet et al. (Pending)")
pdac_studies<-c("VSP34 - Shalhoub et al. (Pending)")
metformin_studies<-c("Metformin - Kulkarni et al. 2019")

reactMeta<-read.csv("human_reaction_meta.csv")
colnames(reactMeta)[1]<-"metabolite"

header <- dashboardHeader(
  title = "RADAR | SCepter",
  titleWidth = 400 
)

body <- dashboardBody(
  
  fluidRow(column(width=12,box(width=6,title="Select assay 1",solidHeader = TRUE, status="primary", 
                               box(width = 6, selectizeInput(label="Select Condition", inputId = "condition_1", choices = c("Acute Myeloid Leukemia (AML)"=1,"Pancreatic Ductal Adenocarcinoma (PDAC)"=2,"Metformin"=3),selected=1)),
                               box(width = 6, selectizeInput(label="Select Study", inputId = "study_1", choices = c("Cytarabine and Venetoclax - Bosc et al. 2021"=1,"Gilteritinib FLT3/CEBPA - Sabatier et al. 2022"=2,"CD39 - Boet et al. (Pending)"=3,"Metformin - Kulkarni et al. 2019"=4),selected=2))),
                  box(width=6,title="Select assay 2",solidHeader = TRUE, status="primary", 
                      box(width = 6, selectizeInput(label="Select Condition", inputId = "condition_2", choices = c("Acute Myeloid Leukemia (AML)"=1,"Pancreatic Ductal Adenocarcinoma (PDAC)"=2,"Metformin"=3),selected=1)),
                      box(width = 6, selectizeInput(label="Select Study", inputId = "study_2", choices = c("Cytarabine and Venetoclax - Bosc et al. 2021"=1,"Gilteritinib FLT3/CEBPA - Sabatier et al. 2022"=2,"CD39 - Boet et al. (Pending)"=3,"Metformin - Kulkarni et al. 2019"=4),selected=2))),
                               ),
  ),
  fluidRow(column(width=12,box(width = 6, title=textOutput("plot1Title"),
                               solidHeader = TRUE, withLoader(plotOutput("plot1", click = "clickBar", height = 500), type = "html", loader = "dnaspin"),
                               box(width = 12, noUiSliderInput(label = "Adjust Resolution", inputId = "resolution_1", min = 0.1, max = 1.0, value = 0.2, tooltip = TRUE, step = 0.1)),
                               box(width = 6, virtualSelectInput(label="Select Sample(s)", inputId = "sample_1", choices = sample_choices_sabatier, showValueAsTags = TRUE, search = TRUE, multiple = TRUE,disableSelectAll=FALSE)),
                               box(width = 6, virtualSelectInput(label="Select Treatment(s)", inputId = "treatment_1", choices = treatment_choices_sabatier, showValueAsTags = TRUE, search = TRUE, multiple = TRUE,disableSelectAll=FALSE)),
                               box(width = 12, virtualSelectInput(label="Select Cluster(s)", inputId = "cluster_1", choices = c(0:9), showValueAsTags = TRUE, search = TRUE, multiple = TRUE, disableSelectAll = FALSE))),
                  box(width = 6, title=textOutput("plot2Title"),
                      solidHeader = TRUE, withLoader(plotOutput("plot2", click = "clickBar", height = 500), type = "html", loader = "dnaspin"),
                      box(width = 12, noUiSliderInput(label = "Adjust Resolution", inputId = "resolution_2", min = 0.1, max = 1.0, value = 0.2, tooltip = TRUE, step = 0.1)),
                      box(width = 6, virtualSelectInput(label="Select Sample(s)", inputId = "sample_2", choices = sample_choices_sabatier, showValueAsTags = TRUE, search = TRUE, multiple = TRUE,disableSelectAll=FALSE)),
                      box(width = 6, virtualSelectInput(label="Select Treatment(s)", inputId = "treatment_2", choices = treatment_choices_sabatier, showValueAsTags = TRUE, search = TRUE, multiple = TRUE,disableSelectAll=FALSE)),
                      box(width = 12, virtualSelectInput(label="Select Cluster(s)", inputId = "cluster_2", choices = c(0:9), showValueAsTags = TRUE, search = TRUE, multiple = TRUE, disableSelectAll = FALSE)))),
                               
           ),
  fluidRow(column(width=12,box(width=12,title="RADAR | SCepter Pre-Analysis Summary",solidHeader = TRUE, status="primary", 
                               box(width = 6, verbatimTextOutput('reportText')),
                               box(width = 6, noUiSliderInput(inputId = "customAUC", label = "AUC:", min = 0.6, max = 1.0, value = 0.65, tooltip = TRUE, step = 0.01),
                                   noUiSliderInput(inputId = "customFC", label = "Log2FC:", min = 0, max = 5, value = 1, tooltip = TRUE, step = 0.01),
                                   noUiSliderInput(inputId = "customFDR", label = "FDR:", min = 0, max = 1, value = 0.05, tooltip = TRUE, step = 0.01)),
                               fluidRow(column(width=12,actionButton("submitReport", "Begin Analysis"))),column(width=12,actionButton("submitReportDefault", "Begin Analysis (default parameters)"))))),
  
  fluidRow(column(width=8,box(title=textOutput("plot1Title_metab"),
                              width = 12, solidHeader = TRUE, status="primary",withLoader(plotOutput("plot1_metab", click = "clickBar", height = 500), type = "html", loader = "dnaspin"),
                              downloadLink("dlPlot1", "Download Plot as PDF | "), downloadLink("fingerprint", "Download Fingerprint Data")),
  ), 
  column(width=4,box(title=textOutput("plot2Title_metab"),
                     width = 12, solidHeader = TRUE, status="primary", withLoader(plotOutput("plot2_metab", brush = "plot_brush2", height = 500), type = "html", loader = "dnaspin"),
                     downloadLink("dlPlot2", "Download Plot as PDF")))),
  fluidRow(column(width=8,box(title=textOutput("outerTitle"), solidHeader = TRUE,status="primary",checkboxGroupInput("outerSelect", label = NULL, choices = c("Significant (+)"="up", "Significant (-)"="down"), selected=c("up","down"))))),
  fluidRow(column(width=12,box(width=12,selectizeInput("boxplot", "Select Metabolic Flux to Analyze", choices = c(""))))),
  fluidRow(column(width=12,box(width=12,reactableOutput("info")))),
  fluidRow(column(width=12,box(width=12,withLoader(plotOutput("plot3_metab", height = 500), type = "html", loader = "dnaspin"),actionLink("invert","Invert Boxplot | "), downloadLink("dlPlot3", "Download Plot as PDF")))),
  
  fluidRow(column(width=12,box(width = 6, title="Assay_1",
                               solidHeader = TRUE, status="primary", withLoader(plotOutput("plot1_feat", click = "clickBar", height = 500), type = "html", loader = "dnaspin")),
                  box(width = 6, title="Assay_2",
                      solidHeader = TRUE, status="primary", withLoader(plotOutput("plot3_feat", click = "clickBar", height = 500), type = "html", loader = "dnaspin")),
                  ),

  ),
  
  
)

ui<-dashboardPage(
  skin = "purple",
  header = header,
  sidebar = dashboardSidebar(disable = TRUE),
  body = body
)

server <-function(input,output,session){
  assay_merged <- reactiveVal(NULL)
  
  assay_1 <- reactiveVal(NULL)
  assay_1_mid <- reactiveVal(NULL)
  assay_1_active<-reactiveVal(NULL)
  
  plot_max<-reactiveVal(NULL)
  plot_min<-reactiveVal(NULL)
  
  flux_1<-reactiveVal(NULL)
  names_1<-reactiveVal(NULL)
  
  errors_global<-reactiveVal(NULL)
  
  observeEvent(list(input$condition_1, input$study_1, input$resolution_1),{
    if(input$condition_1 == 1){
      updateSelectizeInput(session, "study_1", choices = c("Cytarabine and Venetoclax - Bosc et al. 2021"=1,"Gilteritinib FLT3/CEBPA - Sabatier et al. 2022"=2,"CD39 - Boet et al. (Pending)"=3),selected=2)
      if(input$study_1 == 1){
        assay_1(bosc)
        flux_1(bosc_flux)
        assay_1(FindClusters(assay_1(), resolution = input$resolution_1))
        updateVirtualSelect('cluster_1', choices = c(as.data.frame(table(Idents(assay_1())))[[1]]), selected = c(as.data.frame(table(Idents(assay_1())))[[1]]))
        updateVirtualSelect('sample_1', choices = sample_choices_bosc, selected = sample_choices_bosc)
        updateVirtualSelect('treatment_1', choices = treatment_choices_bosc, selected = treatment_choices_bosc)
        assay_1_mid(assay_1())
        assay_1_active(assay_1())
      }else if(input$study_1 == 2){
        assay_1(sabatier)
        flux_1(sabatier_flux)
        assay_1(FindClusters(assay_1(), resolution = input$resolution_1))
        updateVirtualSelect('cluster_1', choices = c(as.data.frame(table(Idents(assay_1())))[[1]]), selected = c(as.data.frame(table(Idents(assay_1())))[[1]]))
        updateVirtualSelect('sample_1', choices = sample_choices_sabatier, selected = sample_choices_sabatier)
        updateVirtualSelect('treatment_1', choices = treatment_choices_sabatier, selected = treatment_choices_sabatier)
        assay_1_mid(assay_1())
        assay_1_active(assay_1())
      }else if(input$study_1 == 3){
        
      }else{
        
      }
    }else if(input$condition_1 == 3){
      assay_1(kulkarni)
      flux_1(kulkarni_flux)
      assay_1(FindClusters(assay_1(), resolution = input$resolution_1))
      updateVirtualSelect('cluster_1', choices = c(as.data.frame(table(Idents(assay_1())))[[1]]), selected = c(as.data.frame(table(Idents(assay_1())))[[1]]))
      updateVirtualSelect('sample_1', choices = sample_choices_kulkarni, selected = sample_choices_kulkarni)
      updateVirtualSelect('treatment_1', choices = treatment_choices_kulkarni, selected = treatment_choices_kulkarni)
      assay_1_mid(assay_1())
      assay_1_active(assay_1())
    }
    else{
      updateSelectizeInput(session, "study_1", choices = c("VSP34 - Shalhoub et al. (Pending)"=1),selected=1)
      if(input$study_1 == 1){
        
      }else{
        
      }
    }
  })
  
  observeEvent(list(input$treatment_1, input$sample_1),{
    req(input$treatment_1)
    req(input$sample_1)
    assay_1_mid(subset(assay_1(), treatment %in% input$treatment_1 & sample %in% input$sample_1))
    updateVirtualSelect('cluster_1', choices = c(as.data.frame(table(Idents(assay_1_mid())))[[1]]), selected = c(as.data.frame(table(Idents(assay_1_mid())))[[1]]))
    assay_1_active(assay_1_mid())
    names_1(rownames(assay_1_active()@meta.data))
  })
  
  observeEvent(input$cluster_1,{
    req(input$cluster_1)
    if(is.null(assay_1_mid())){
      assay_1_active(subset(assay_1(), seurat_clusters %in% input$cluster_1))
      names_1(rownames(assay_1_active()@meta.data))
    }else{
      assay_1_active(subset(assay_1_mid(), seurat_clusters %in% input$cluster_1))
      names_1(rownames(assay_1_active()@meta.data))
    }
  })
  
  plot1 <- reactive({
    if(is.null(assay_1_active())){
      DimPlot(assay_1_mid())
    }else{
      DimPlot(assay_1_active())
    }
  })
  
  output$plot1<- renderPlot({
    plot1()
  })
  
  assay_2 <- reactiveVal(NULL)
  assay_2_mid <- reactiveVal(NULL)
  assay_2_active<-reactiveVal(NULL)
  
  flux_2<-reactiveVal(NULL)
  names_2<-reactiveVal(NULL)
  
  observeEvent(list(input$condition_2, input$study_2, input$resolution_2),{
    if(input$condition_2 == 1){
      updateSelectizeInput(session, "study_2", choices = c("Cytarabine and Venetoclax - Bosc et al. 2021"=1,"Gilteritinib FLT3/CEBPA - Sabatier et al. 2022"=2,"CD39 - Boet et al. (Pending)"=3),selected=2)
      if(input$study_2 == 1){
        assay_2(bosc)
        flux_2(bosc_flux)
        assay_2(FindClusters(assay_2(), resolution = input$resolution_2))
        updateVirtualSelect('cluster_2', choices = c(as.data.frame(table(Idents(assay_2())))[[1]]), selected = c(as.data.frame(table(Idents(assay_2())))[[1]]))
        updateVirtualSelect('sample_2', choices = sample_choices_bosc, selected = sample_choices_bosc)
        updateVirtualSelect('treatment_2', choices = treatment_choices_bosc, selected = treatment_choices_bosc)
      }else if(input$study_2 == 2){
        assay_2(sabatier)
        flux_2(sabatier_flux)
        assay_2(FindClusters(assay_2(), resolution = input$resolution_2))
        updateVirtualSelect('cluster_2', choices = c(as.data.frame(table(Idents(assay_2())))[[1]]), selected = c(as.data.frame(table(Idents(assay_2())))[[1]]))
        updateVirtualSelect('sample_2', choices = sample_choices_sabatier, selected = sample_choices_sabatier)
        updateVirtualSelect('treatment_2', choices = treatment_choices_sabatier, selected = treatment_choices_sabatier)
      }else if(input$study_2 == 3){
        assay_2(boet)
        flux_2(boet_flux)
        assay_2(FindClusters(assay_2(), resolution = input$resolution_2))
        updateVirtualSelect('cluster_2', choices = c(as.data.frame(table(Idents(assay_2())))[[1]]), selected = c(as.data.frame(table(Idents(assay_2())))[[1]]))
        updateVirtualSelect('sample_2', choices = sample_choices_boet, selected = sample_choices_boet)
        updateVirtualSelect('treatment_2', choices = treatment_choices_boet, selected = treatment_choices_boet)
      }
    }else if(input$condition_2 == 3){
      assay_2(kulkarni)
      flux_2(kulkarni_flux)
      assay_2(FindClusters(assay_2(), resolution = input$resolution_2))
      updateVirtualSelect('cluster_2', choices = c(as.data.frame(table(Idents(assay_2())))[[1]]), selected = c(as.data.frame(table(Idents(assay_2())))[[1]]))
      updateVirtualSelect('sample_2', choices = sample_choices_kulkarni, selected = sample_choices_kulkarni)
      updateVirtualSelect('treatment_2', choices = treatment_choices_kulkarni, selected = treatment_choices_kulkarni)
      assay_2_mid(assay_2())
      assay_2_active(assay_2())
    }else{
      updateSelectizeInput(session, "study_2", choices = c("VSP34 - Shalhoub et al. (Pending)"=1),selected=1)
      if(input$study_2 == 1){
        
      }else{
        
      }
    }
  })
  
  observeEvent(list(input$treatment_2, input$sample_2),{
    req(input$treatment_2)
    req(input$sample_2)
    assay_2_mid(subset(assay_2(), treatment %in% input$treatment_2 & sample %in% input$sample_2))
    updateVirtualSelect('cluster_2', choices = c(as.data.frame(table(Idents(assay_2_mid())))[[1]]), selected = c(as.data.frame(table(Idents(assay_2_mid())))[[1]]))
    assay_2_active(assay_2_mid())
    names_2(rownames(assay_2_active()@meta.data))
  })
  
  observeEvent(input$cluster_2,{
    req(input$cluster_2)
    if(is.null(assay_2_mid())){
      assay_2_active(subset(assay_2(), seurat_clusters %in% input$cluster_2))
      names_2(rownames(assay_2_active()@meta.data))
    }else{
      assay_2_active(subset(assay_2_mid(), seurat_clusters %in% input$cluster_2))
      names_2(rownames(assay_2_active()@meta.data))
      
    }

  })
  
  plot2 <- reactive({
    if(is.null(assay_2_active())){
      DimPlot(assay_2_mid())
    }else{
      DimPlot(assay_2_active())
    }
  })
  
  output$plot2<- renderPlot({
    plot2()
  })
  
  observeEvent(list(input$condition_1, input$study_1, input$resolution_1,input$treatment_1, input$sample_1, input$cluster_1,input$condition_2, input$study_2, input$resolution_2,input$treatment_2, input$sample_2, input$cluster_2), {
    
    warnings<-c()
    errors<-c()
    
    metaFinal_out <- reactiveVal(NULL)
    
    if(is.null(warnings)){
      warnings<-"NONE"
    }else{
      warnings<-paste0(warnings, collapse = " - ")
    }
    
    if(is.null(errors)){
      errors<-"NONE"
    }else{
      errors<-paste0(errors, collapse = " - ")
    }
    
    errors_global(errors)
    output$reportText <- renderText({paste(paste0("Assay 1: ",paste0(input$condition_1," - ",input$study_1),"\n"), 
                                           paste0("Assay 2: ",paste0(input$condition_2," - ",input$study_2),"\n"),
                                           paste0("Minimum AUC: ",input$customAUC,"\n"),
                                           paste0("Minimum Log2FC: ",input$customFC,"\n"),
                                           paste0("Maximum FDR: ",input$customFDR,"\n"),
                                           paste0("Report Generated at: ", Sys.time()))})
    

  })
  observeEvent(input$submitReportDefault,{
    if(errors_global()=="NONE"){
      print(1)
      flux<-rbind(flux_1(), flux_2())
      flux<-flux[unique(rownames(flux)),]
      
      assay_meta_1<-as.data.frame(matrix(nrow=length(names_1()),ncol=1))
      colnames(assay_meta_1)<-c("selected_assay")
      rownames(assay_meta_1)<-names_1()
      assay_meta_1$selected_assay<-"assay_1"
      
      assay_meta_2<-as.data.frame(matrix(nrow=length(names_2()),ncol=1))
      colnames(assay_meta_2)<-c("selected_assay")
      rownames(assay_meta_2)<-names_2()
      assay_meta_2$selected_assay<-"assay_2"
      
      metaFinal<-rbind(assay_meta_1, assay_meta_2)
      
      
      flux_full<-cbind(metaFinal,flux[rownames(metaFinal),])
      flux_full<-flux_full[complete.cases(flux_full[,2]),]
      metaFinal<-as.data.frame(metaFinal[rownames(flux_full),])
      colnames(metaFinal)[1]<-"selected_assay"
      rownames(metaFinal)<-rownames(flux_full)
      
      flux_comp<-flux_full
      print(3)
      outer<-colnames(metaFinal)[1]
      
      outer_select_1<-c("assay_1", "assay_2")
      metab_names<-c()
      data_filt_var<-flux_comp
      
      AUC = colAUC(data_filt_var[,-1],factor(data_filt_var[,"selected_assay"]),plotROC=FALSE)
      AUC = AUC[1,]
      plot(density(AUC),main=paste0("AUC for features"))
      print(3)
      logFC = log2(apply(subset(data_filt_var, subset = selected_assay == "assay_1")[,-1],2,mean) / apply(subset(data_filt_var, subset = selected_assay == "assay_2")[,-1],2,mean))
      logFC[AUC<input$customAUC]=0
      print(4)
      short<-as.data.frame(logFC[which(abs(logFC)>input$customFC)])
      print(5)
      t<-col_t_welch(as.matrix(subset(data_filt_var, subset = selected_assay == "assay_1")[,rownames(short)]),as.matrix(subset(data_filt_var, subset = selected_assay == "assay_2")[,rownames(short)]))
      final<-cbind(t$mean.y, t$mean.x, short, t$statistic, t$pvalue)
      colnames(final)<-c("var_1_mean","var_2_mean","log2FC","t_welch_statistic","pval")
      print(6)
      final<-subset(final, subset = pval < 0.05)
      final$fdr = p.adjust(final[,"pval"],method="fdr")
      final<-subset(final, subset = fdr < 0.05)
      metab_names<-c(metab_names, rownames(final))
      
      metab_names<-unique(metab_names)
      
      test_list<-data.frame()
      
      names<-c(outer, metab_names)
      
      
      data_filt_var<-data_filt_var[,names]
      t<-col_t_welch(as.matrix(subset(data_filt_var, subset = data_filt_var[,outer] == "assay_1")[,metab_names]),as.matrix(subset(data_filt_var, subset = data_filt_var[,outer] == "assay_2")[,metab_names]))
      final<-as.data.frame(cbind(metabolite=metab_names,as.numeric(t$mean.y), as.numeric(t$mean.x),as.numeric(t$statistic), as.numeric(t$pvalue)))
      colnames(final)<-c("metabolite","var_1_mean","var_2_mean","t_welch_statistic","pval")
      final$sign_t_log_pval = ifelse(final$t_welch_statistic>0, -log10(as.numeric(final$pval)), log10(as.numeric(final$pval)))
      test_list<-as.data.frame(rbind(test_list,final))
      
      
      metabs_tab<-merge(test_list, reactMeta, by = "metabolite", all.x=TRUE, all.y=FALSE)
      
      metabs_tab$sig<-ifelse(metabs_tab$sign_t_log_pval > 1.301, "up", ifelse(metabs_tab$sign_t_log_pval < -1.301, "down", "ns"))
      
      output$outerTitle<-renderText("assay_1 v. assay_2")
      
      
      output$plot1Title_metab<-renderText(paste0("Significant Reactions Across All Subsystems - (assay_1 v. assay_2)"))
      
      metabs_fin<-reactive({
        subset(metabs_tab, sig %in% input$outerSelect)
      })
      
      plot1_metab <- reactive({
        ggplot(metabs_fin(), aes(x=reorder(subsystem, -as.numeric(t_welch_statistic)), y=as.numeric(t_welch_statistic), color=sig, fill = sig)) + geom_point(size = 1, position = position_dodge(width=0.5))+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)+xlab("subsystem")
      })
      
      output$plot1_metab<- renderPlot({
        plot1_metab()
      })
      
      plot1pdf <- reactive({
        ggplot(metabs_fin(), aes(x=as.numeric(t_welch_statistic), y=reorder(subsystem, -as.numeric(t_welch_statistic)), color=sig, fill = sig)) + geom_point(size = 1, position = position_dodge(width=0.5))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_vline(xintercept=1.3, size=.05)+geom_vline(xintercept=-1.3, size=.05)+ylab("subsystem")+xlab("log(pval) with sign of t-statistic")+ggtitle(paste0("Significant Reactions Across All Subsystems - (assay_1 v. assay_2)"))
      })
      
      output$fingerprint <- downloadHandler(
        filename = function() {
          paste0(paste0(input$outerSelect, collapse = "_"), outer, '_filtered_reactions_.rds')
        },
        content = function(file) {
          withProgress(message = 'Generating file', value = 0, {
            incProgress(0.3, detail = "Preparing data")
            
            # Prepare the data to be saved
            data_to_save <- list(
              metabs_fin(),
              metaFinal,
              data.frame(AUC = input$customAUC, Log2FC = input$customFC, FDR = input$customFDR),
              "single_cell"
            )
            
            incProgress(0.3, detail = "Saving file")
            
            # Save the data
            saveRDS(data_to_save, file)
            
            incProgress(0.4, detail = "Finalizing")
          })
        },
        contentType = "application/rds"
      )
      outputOptions(output, "fingerprint", suspendWhenHidden = FALSE)
      
      systems<-reactive({
        ggbld <- ggplot_build(ggplot(metabs_fin(), aes(x=reorder(subsystem, -as.numeric(t_welch_statistic)), y=as.numeric(t_welch_statistic), color=sig, fill = sig)) + geom_point(size = 1, position = position_dodge(width=0.5))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)+xlab("subsystem"))
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
      outputOptions(output, "dlPlot1", suspendWhenHidden = FALSE)
      observeEvent(input$clickBar, {
        groupId <- round(input$clickBar$x)
        output$plot2Title_metab<-renderText(paste0("Significant Reactions - ",systems()[groupId], " (assay_1 v. assay_2)"))
        plot2_metab <- reactive({
          ggplot(subset(metabs_fin(), subset = subsystem == systems()[groupId]), aes(x=subsystem, y=as.numeric(t_welch_statistic), color=sig, fill = sig)) + geom_point(size = 2)+geom_line(aes(group = metabolite, x=subsystem),color="grey")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)
        })
        output$plot2_metab<- renderPlot({
          plot2_metab()
        })
        plot2pdf <- reactive({
          ggplot(subset(metabs_fin(), subset = subsystem == systems()[groupId]), aes(x=subsystem, y=as.numeric(t_welch_statistic), color=sig, fill = sig)) + geom_point(size = 2)+geom_line(aes(group = metabolite, x=subsystem),color="grey")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)+ggtitle(paste0("Significant Reactions  - ",systems[groupId], " (assay_1 v. assay_2)"))
        })
        output$dlPlot2 <- downloadHandler(
          filename = function() {
            paste(str_replace(systems()[groupId]," ","_"),'_',outer,'.pdf', sep='')
          },
          content = function(file) {
            ggsave(file, plot2pdf(), width = 8, height = 12, dpi = 400, units = "in")
          }
        )
        
        View(metabs_fin())
        output$info<-renderReactable({
          metabs_fin<-metabs_fin()[,-c(5,6)]
          metabs_fin$var_1_mean<-round(as.numeric(metabs_fin$var_1_mean), digits=4)
          metabs_fin$var_2_mean<-round(as.numeric(metabs_fin$var_2_mean), digits=4)
          metabs_fin$t_welch_statistic<-round(as.numeric(metabs_fin$t_welch_statistic), digits=4)
          reactable(brushedPoints(subset(metabs_fin, subset = subsystem == systems()[groupId]), input$plot_brush2, xvar = "subsystem", yvar="t_welch_statistic"))
        })
        data <- reactive({
          as.data.frame(brushedPoints(subset(metabs_fin, subset = subsystem == systems()[groupId]), input$plot_brush2, xvar = "subsystem", yvar="t_welch_statistic"))
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
                               choices = brushedPoints(subset(metabs_fin(), subset = subsystem == systems()[groupId]), input$plot_brush2, xvar = "subsystem", yvar="t_welch_statistic")$metabolite
          )})
      })
      plot3_metab <- reactive({
        if(input$boxplot != ""){
          ggplot(data_filt_var, aes_string(x="selected_assay", y=input$boxplot))+geom_boxplot()+geom_jitter()+ylab(input$boxplot)
        }else{
          ggplot()
        }
      })
      
      output$plot3_metab<- renderPlot({
        plot3_metab()
      })
      
      output$dlPlot3 <- downloadHandler(
        filename = function() {
          paste(input$boxplot,'_',label, '.pdf', sep='')
        },
        content = function(file) {
          ggsave(file, plot3_metab(), width = 12, height = 8, dpi = 400, units = "in")
        }
      )
      
      plot1_feat <- reactive({
        if(input$boxplot != ""){
          assay_1_active(AddMetaData(assay_1_active(),flux[rownames(assay_1_active()@meta.data),input$boxplot],col.name=input$boxplot))
          assay_2_active(AddMetaData(assay_2_active(),flux[rownames(assay_2_active()@meta.data),input$boxplot],col.name=input$boxplot))
          
          plot_max(max(c(assay_1_active()@meta.data[,input$boxplot],assay_2_active()@meta.data[,input$boxplot])))
          plot_min(min(c(assay_1_active()@meta.data[,input$boxplot],assay_2_active()@meta.data[,input$boxplot])))
          FeaturePlot(assay_1_active(), input$boxplot)+ scale_color_gradient2(low = "blue2",
                                                                              mid = "lavenderblush",
                                                                              high = "firebrick1", 
                                                                              midpoint = 0, limits = c(plot_min(),plot_max()))
        }else{
          ggplot()
        }
      })
      
      output$plot1_feat<- renderPlot({
        plot1_feat()
      })
      
      plot3_feat <- reactive({
        if(input$boxplot != ""){
          FeaturePlot(assay_2_active(), input$boxplot)+ scale_color_gradient2(low = "blue2",
                                                                              mid = "lavenderblush",
                                                                              high = "firebrick1", 
                                                                              midpoint = 0, limits = c(plot_min(),plot_max()))
        }else{
          ggplot()
        }
      })
      
      output$plot3_feat<- renderPlot({
        plot3_feat()
      })
      
    }
  })
  
  observeEvent(input$submitReport,{
    if(errors_global()=="NONE"){
      print(1)
      flux<-rbind(flux_1(), flux_2())
      flux<-flux[unique(rownames(flux)),]
      
      assay_meta_1<-as.data.frame(matrix(nrow=length(names_1()),ncol=1))
      colnames(assay_meta_1)<-c("selected_assay")
      rownames(assay_meta_1)<-names_1()
      assay_meta_1$selected_assay<-"assay_1"
      
      assay_meta_2<-as.data.frame(matrix(nrow=length(names_2()),ncol=1))
      colnames(assay_meta_2)<-c("selected_assay")
      rownames(assay_meta_2)<-names_2()
      assay_meta_2$selected_assay<-"assay_2"
      
      metaFinal<-rbind(assay_meta_1, assay_meta_2)
      
      
      flux_full<-cbind(metaFinal,flux[rownames(metaFinal),])
      flux_full<-flux_full[complete.cases(flux_full[,2]),]
      metaFinal<-as.data.frame(metaFinal[rownames(flux_full),])
      colnames(metaFinal)[1]<-"selected_assay"
      rownames(metaFinal)<-rownames(flux_full)
      
      flux_comp<-flux_full
      print(3)
      outer<-colnames(metaFinal)[1]
      
      outer_select_1<-c("assay_1", "assay_2")
      metab_names<-c()
      data_filt_var<-flux_comp
      
      AUC = colAUC(data_filt_var[,-1],factor(data_filt_var[,"selected_assay"]),plotROC=FALSE)
      AUC = AUC[1,]
      plot(density(AUC),main=paste0("AUC for features"))
      print(3)
      logFC = log2(apply(subset(data_filt_var, subset = selected_assay == "assay_1")[,-1],2,mean) / apply(subset(data_filt_var, subset = selected_assay == "assay_2")[,-1],2,mean))
      logFC[AUC<input$customAUC]=0
      print(4)
      short<-as.data.frame(logFC[which(abs(logFC)>input$customFC)])
      print(5)
      t<-col_t_welch(as.matrix(subset(data_filt_var, subset = selected_assay == "assay_1")[,rownames(short)]),as.matrix(subset(data_filt_var, subset = selected_assay == "assay_2")[,rownames(short)]))
      final<-cbind(t$mean.y, t$mean.x, short, t$statistic, t$pvalue)
      colnames(final)<-c("var_1_mean","var_2_mean","log2FC","t_welch_statistic","pval")
      print(6)
      final<-subset(final, subset = pval < 0.05)
      final$fdr = p.adjust(final[,"pval"],method="fdr")
      final<-subset(final, subset = fdr < input$customFDR)
      metab_names<-c(metab_names, rownames(final))
      
      metab_names<-unique(metab_names)
      
      test_list<-data.frame()
      
      names<-c(outer, metab_names)
      
      
      data_filt_var<-data_filt_var[,names]
      t<-col_t_welch(as.matrix(subset(data_filt_var, subset = data_filt_var[,outer] == "assay_1")[,metab_names]),as.matrix(subset(data_filt_var, subset = data_filt_var[,outer] == "assay_2")[,metab_names]))
      final<-as.data.frame(cbind(metabolite=metab_names,as.numeric(t$mean.y), as.numeric(t$mean.x),as.numeric(t$statistic), as.numeric(t$pvalue)))
      colnames(final)<-c("metabolite","var_1_mean","var_2_mean","t_welch_statistic","pval")
      final$sign_t_log_pval = ifelse(final$t_welch_statistic>0, -log10(as.numeric(final$pval)), log10(as.numeric(final$pval)))
      test_list<-as.data.frame(rbind(test_list,final))
      
      
      metabs_tab<-merge(test_list, reactMeta, by = "metabolite", all.x=TRUE, all.y=FALSE)
      
      metabs_tab$sig<-ifelse(metabs_tab$sign_t_log_pval > 1.301, "up", ifelse(metabs_tab$sign_t_log_pval < -1.301, "down", "ns"))
      
      output$outerTitle<-renderText("assay_1 v. assay_2")
      
      
      output$plot1Title_metab<-renderText(paste0("Significant Reactions Across All Subsystems - (assay_1 v. assay_2)"))
      
      metabs_fin<-reactive({
        subset(metabs_tab, sig %in% input$outerSelect)
      })
      
      plot1_metab <- reactive({
        ggplot(metabs_fin(), aes(x=reorder(subsystem, -as.numeric(t_welch_statistic)), y=as.numeric(t_welch_statistic), color=sig, fill = sig)) + geom_point(size = 1, position = position_dodge(width=0.5))+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)+xlab("subsystem")
      })
      
      output$plot1_metab<- renderPlot({
        plot1_metab()
      })
      
      plot1pdf <- reactive({
        ggplot(metabs_fin(), aes(x=reorder(subsystem, -as.numeric(t_welch_statistic)), y=as.numeric(t_welch_statistic), color=sig, fill = sig)) + geom_point(size = 1, position = position_dodge(width=0.5))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)+xlab("subsystem")+ggtitle(paste0("Significant Reactions Across All Subsystems - (assay_1 v. assay_2)"))
      })

      
      output$fingerprint <- downloadHandler(
        filename = function() {
          paste(paste0(input$outerSelect, collapse = "_"), outer, '_filtered_reactions_.rds', sep='')
        },
        content = function(file) {
          withProgress(message = 'Generating file', value = 0, {
            incProgress(0.3, detail = "Preparing data")
            
            # Prepare the data to be saved
            data_to_save <- list(
              metabs_fin(),
              metaFinal_condensed,
              data.frame(AUC = input$customAUC, Log2FC = input$customFC, FDR = input$customFDR),
              "single_cell"
            )
            
            incProgress(0.3, detail = "Saving file")
            
            # Save the data
            saveRDS(data_to_save, file)
            
            incProgress(0.4, detail = "Finalizing")
          })
        }
      )
      outputOptions(output, "fingerprint", suspendWhenHidden = FALSE)

            systems<-reactive({
              ggbld <- ggplot_build(ggplot(metabs_fin(), aes(x=reorder(subsystem, -as.numeric(t_welch_statistic)), y=as.numeric(t_welch_statistic), color=sig, fill = sig)) + geom_point(size = 1, position = position_dodge(width=0.5))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)+xlab("subsystem"))
              ggbld$layout$panel_params[[1]]$x$limits
            })


            output$dlPlot1 <- downloadHandler(
              filename = function() {
                paste('all_subsystems_.pdf', sep='')
              },
              content = function(file) {
                ggsave(file, plot1pdf(), width = 24, height = 12, dpi = 400, units = "in")
              }
            )
            observeEvent(input$clickBar, {
              groupId <- round(input$clickBar$x)
              output$plot2Title_metab<-renderText(paste0("Significant Reactions - ",systems()[groupId], " (assay_1 v. assay_2)"))
              plot2_metab <- reactive({
                ggplot(subset(metabs_fin(), subset = subsystem == systems()[groupId]), aes(x=subsystem, y=as.numeric(t_welch_statistic), color=sig, fill = sig)) + geom_point(size = 2)+geom_line(aes(group = metabolite, x=subsystem),color="grey")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)
              })
              output$plot2_metab<- renderPlot({
                plot2_metab()
              })
              plot2pdf <- reactive({
                ggplot(subset(metabs_fin(), subset = subsystem == systems()[groupId]), aes(x=subsystem, y=as.numeric(t_welch_statistic), color=sig, fill = sig)) + geom_point(size = 2)+geom_line(aes(group = metabolite, x=subsystem),color="grey")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)+ggtitle(paste0("Significant Reactions  - ",systems[groupId], " (assay_1 v. assay_2)"))
              })
              output$dlPlot2 <- downloadHandler(
                filename = function() {
                  paste(str_replace(systems()[groupId]," ","_"),'_',outer,'.pdf', sep='')
                },
                content = function(file) {
                  ggsave(file, plot2pdf(), width = 8, height = 12, dpi = 400, units = "in")
                }
              )
              
              View(metabs_fin())
              output$info<-renderReactable({
                metabs_fin<-metabs_fin()[,-c(5,6)]
                metabs_fin$var_1_mean<-round(as.numeric(metabs_fin$var_1_mean), digits=4)
                metabs_fin$var_2_mean<-round(as.numeric(metabs_fin$var_2_mean), digits=4)
                metabs_fin$t_welch_statistic<-round(as.numeric(metabs_fin$t_welch_statistic), digits=4)
                reactable(brushedPoints(subset(metabs_fin, subset = subsystem == systems()[groupId]), input$plot_brush2, xvar = "subsystem", yvar="t_welch_statistic"))
              })
              data <- reactive({
                as.data.frame(brushedPoints(subset(metabs_fin, subset = subsystem == systems()[groupId]), input$plot_brush2, xvar = "subsystem", yvar="t_welch_statistic"))
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
                                     choices = brushedPoints(subset(metabs_fin(), subset = subsystem == systems()[groupId]), input$plot_brush2, xvar = "subsystem", yvar="t_welch_statistic")$metabolite
                )})
            })
             plot3_metab <- reactive({
               if(input$boxplot != ""){
                  ggplot(data_filt_var, aes_string(x="selected_assay", y=input$boxplot))+geom_boxplot()+geom_jitter()+ylab(input$boxplot)
               }else{
                 ggplot()
               }
            })

            output$plot3_metab<- renderPlot({
              plot3_metab()
            })

            output$dlPlot3 <- downloadHandler(
              filename = function() {
                paste(input$boxplot,'_',label, '.pdf', sep='')
              },
              content = function(file) {
                ggsave(file, plot3_metab(), width = 12, height = 8, dpi = 400, units = "in")
              }
            )
            
            plot1_feat <- reactive({
              if(input$boxplot != ""){
                assay_1_active(AddMetaData(assay_1_active(),flux[rownames(assay_1_active()@meta.data),input$boxplot],col.name=input$boxplot))
                assay_2_active(AddMetaData(assay_2_active(),flux[rownames(assay_2_active()@meta.data),input$boxplot],col.name=input$boxplot))

                plot_max(max(c(assay_1_active()@meta.data[,input$boxplot],assay_2_active()@meta.data[,input$boxplot])))
                plot_min(min(c(assay_1_active()@meta.data[,input$boxplot],assay_2_active()@meta.data[,input$boxplot])))
                FeaturePlot(assay_1_active(), input$boxplot)+ scale_color_gradient2(low = "blue2",
                                                                                    mid = "lavenderblush",
                                                                                    high = "firebrick1", 
                                                                                    midpoint = 0, limits = c(plot_min(),plot_max()))
              }else{
                ggplot()
              }
            })
            
            output$plot1_feat<- renderPlot({
              plot1_feat()
            })
            
            plot3_feat <- reactive({
              if(input$boxplot != ""){
                FeaturePlot(assay_2_active(), input$boxplot)+ scale_color_gradient2(low = "blue2",
                                                                                    mid = "lavenderblush",
                                                                                    high = "firebrick1", 
                                                                                    midpoint = 0, limits = c(plot_min(),plot_max()))
              }else{
                ggplot()
              }
            })
            
            output$plot3_feat<- renderPlot({
              plot3_feat()
            })

    }
  })
}

shinyApp(ui,server)
