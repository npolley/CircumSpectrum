library(data.table)

#CORALIE Tier 1
assay_metabolites_tot<-read.csv("beataml2_flux.csv", header = T, row.names = 1)
assay_metabolites<-predict(preProcess(assay_metabolites_tot,method = c("center", "scale")), assay_metabolites_tot)
assay_metabolites<- na.omit(assay_metabolites)

beataml_clin<-as.data.frame(as.matrix(fread("beataml2_clinical.csv"), rownames = 1))
beataml_flt3_pos_assay<-assay_metabolites[rownames(subset(beataml_clin, FLT3.ITD == "positive")),]
beataml_flt3_neg_assay<-assay_metabolites[rownames(subset(beataml_clin, FLT3.ITD == "negative")),]

beataml_gilt_res<-assay_metabolites[rownames(subset(beataml_clin, GILTERITINIB_RESPONSE == "resistant" & FLT3.ITD == "positive")),]
beataml_gilt_sens<-assay_metabolites[rownames(subset(beataml_clin, GILTERITINIB_RESPONSE == "sensitive" & FLT3.ITD == "positive")),]

QUIZ_MOLM14_1<-readRDS("QUIZ_MOLM14_1_.rds")
QUIZ_MOLM14_2<-readRDS("QUIZ_MOLM14_2_.rds")
QUIZ_MV411<-readRDS("QUIZ_MV411_.rds")
QUIZ_PDX<-readRDS("QUIZ_PDX_.rds")
GILT_PDX<-readRDS("GILT_PDX_.rds")

BEATAML_FLT3_MUT<-readRDS("BEATAML_FLT3_MUT_.rds")
BEATAML_CEBPA_MUT<-readRDS("BEATAML_CEBPA_MUT_.rds")

BEATAML_SCD<-readRDS("BEATAML_SCD_.rds")
BEATAML_FASN<-readRDS("BEATAML_FASN_.rds")
BEATAML_FADS1<-readRDS("BEATAML_FADS1_.rds")
BEATAML_FADS2<-readRDS("BEATAML_FADS2_.rds")
BEATAML_GILT_RESPONSE<-readRDS("BEATAML_GILT_RESPONSE_.rds")

fingerprints<-list(QUIZ_MOLM14_1, QUIZ_MOLM14_2, QUIZ_MV411, QUIZ_PDX, GILT_PDX, BEATAML_FLT3_MUT, BEATAML_CEBPA_MUT, BEATAML_SCD, BEATAML_FASN, BEATAML_FADS1, BEATAML_FADS2,BEATAML_GILT_RESPONSE)
fingerprint_names<-list("QUIZ_MOLM14_1", "QUIZ_MOLM14_2", "QUIZ_MV411", "QUIZ_PDX", "GILT_PDX", "BEATAML_FLT3_MUT", "BEATAML_CEBPA_MUT", "BEATAML_SCD","BEATAML_FASN", "BEATAML_FADS1", "BEATAML_FADS2","BEATAML_GILT_RESPONSE")

total_scores<-list()
subsystem_scores<-list()

for(i in 1:length(fingerprints)){
  fingerprint<-fingerprints[[i]]
  fingerprint_name<-fingerprint_names[[i]]
  print(fingerprint_name)
  
  auc_list<-fingerprint$All$auc_primary
  
  index<-rownames(auc_list, subsystem)
  system_names<-auc_list[index,]$subsystem
  
  final_classifiers<-fingerprint$All$fingerprints_primary
  
  ex_subsystems_no_norm<-list()
  
  #Generate prediction values for each individual XGBoost subsystem model on reference assay
  for(x in 1:length(index)){
    newDat<-as.matrix(beataml_gilt_sens[,final_classifiers[[as.numeric(index[x])]]$coefnames])
    colnames(newDat)<-final_classifiers[[as.numeric(index[x])]]$coefnames
    pred_tot<-predict(final_classifiers[[as.numeric(index[x])]], newdata = newDat,type = "prob", na.action = na.pass)[,1]
    ex_subsystems_no_norm[[system_names[x]]] <- pred_tot
  }
  
  ex_subsystems_no_norm<-data.frame(ex_subsystems_no_norm)
  
  score<-rowMeans(ex_subsystems_no_norm)
  score<-(score - min(score, na.rm = TRUE)) / (max(score, na.rm = TRUE) - min(score, na.rm = TRUE))
  
  total_scores[[fingerprint_name]]<-score
  subsystem_scores[[fingerprint_name]]<-ex_subsystems_no_norm
}

total_scores<-as.data.frame(total_scores)

cor_table<-cor(total_scores, method = "spearman")

ref_assay_name<-"BEATAML2 - Patients at Diagnosis (all patients)"

corrplot(cor_table, method = 'ellipse', type = 'upper', title = paste0("Correlations of Fingerprint Predictions\nReference Assay: ",ref_assay_name), tl.cex = 1, addCoef.col = "grey55",mar=c(0,0,2.5,0))

assay<-readRDS("sabatier_TUH110_TUH93_TUH84.rds")
assay <- FindClusters(assay, resolution = 0.4)

assay_metabolites_tot<-read.csv("sabatier_flux.csv", header = T, row.names = 1)

TUH84<-subset(assay, sample == "TUH84" & (treatment == "CTL_D0" | treatment == "GILT_D7"))

DimPlot(TUH84, split.by = "treatment")

pre_therapy <- subset(TUH84, subset = treatment == "CTL_D0")
post_therapy <- subset(TUH84, subset = treatment == "GILT_D7")

pre_proportions <- table(Idents(pre_therapy)) / ncol(pre_therapy)
post_proportions <- table(Idents(post_therapy)) / ncol(post_therapy)

# Create cluster_changes dataframe
cluster_changes <- data.frame(
  cluster = names(pre_proportions),
  CTL_D0 = as.numeric(pre_proportions),
  GILT_D7 = as.numeric(post_proportions),
  fold_change = as.numeric(post_proportions) / as.numeric(pre_proportions)
)

# Function to perform prop.test and return p-value
prop_test <- function(pre_prop, post_prop, pre_total, post_total) {
  x <- c(pre_prop * pre_total, post_prop * post_total)
  n <- c(pre_total, post_total)
  test_result <- prop.test(x, n, alternative = "two.sided", correct = FALSE)
  return(test_result$p.value)
}

# Calculate p-values
cluster_changes$p_value <- mapply(prop_test, 
                                  cluster_changes$CTL_D0, 
                                  cluster_changes$GILT_D7, 
                                  MoreArgs = list(pre_total = ncol(pre_therapy), 
                                                  post_total = ncol(post_therapy)))

# Add significance stars
cluster_changes$significance <- cut(cluster_changes$p_value, 
                                    breaks = c(0, 0.001, 0.01, 0.05, 1), 
                                    labels = c("***", "**", "*", ""))

# Categorize clusters
cluster_changes$category <- ifelse(cluster_changes$p_value < 0.05,
                                   ifelse(cluster_changes$fold_change > 1, "Resistant", "Sensitive"),
                                   "Persistent")

cluster_changes <- cluster_changes %>%
  mutate(cluster_num = as.numeric(str_extract(cluster, "\\d+"))) %>%
  arrange(cluster_num)

# Set the cluster factor levels in the correct order
cluster_changes$cluster <- factor(cluster_changes$cluster, levels = cluster_changes$cluster)

plot_data <- cluster_changes %>%
  select(cluster, CTL_D0, GILT_D7, significance, category) %>%
  pivot_longer(cols = c(CTL_D0, GILT_D7), names_to = "condition", values_to = "proportion")

# Create the plot
ggplot(plot_data, aes(x = cluster, y = proportion, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
  geom_text(data = cluster_changes, 
            aes(x = cluster, y = pmax(CTL_D0,GILT_D7)+0.01, label = significance),
            inherit.aes = FALSE) +
  geom_text(data = cluster_changes,
            aes(x = cluster, y = -0.05, label = category),
            angle = 90, hjust = 1, vjust = 0.5, inherit.aes = FALSE) +
  labs(x = "Cluster", y = "Proportion of Cells", fill = "Condition") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  coord_cartesian(ylim = c(-0.1, max(plot_data$proportion) * 1.1))

assay_metabolites<-assay_metabolites_tot[rownames(TUH84@meta.data),]
assay_metabolites<-predict(preProcess(assay_metabolites,method = c("center", "scale")), assay_metabolites)
assay_metabolites <- na.omit(assay_metabolites)

TUH84 <- subset(TUH84, cells = rownames(assay_metabolites))

TUH84$GILT_status<-ifelse(TUH84$seurat_clusters %in% subset(cluster_changes, category == "Sensitive")$cluster, "sensitive",ifelse(TUH84$seurat_clusters %in% subset(cluster_changes, category == "Persistent")$cluster,"persistent","resistant"))
TUH84$GILT_status<-factor(TUH84$GILT_status, levels=c("resistant","persistent","sensitive"))

DimPlot(TUH84, group.by = "GILT_status", split.by = "treatment")

TUH84$GILT_status<-factor(TUH84$GILT_status, levels=c("sensitive","persistent","resistant"))

assay_metabolites<-assay_metabolites[rownames(TUH84@meta.data),]

scale_type = "binary"

total_scores<-list()
subsystem_scores<-list()

pdf(paste0("TUH84_by_status.pdf"),width = 11.8, height = 4.4)
sample_name<-"TUH84"
for(i in 1:length(fingerprints)){
  fingerprint<-fingerprints[[i]]
  fingerprint_name<-fingerprint_names[[i]]
  print(fingerprint_name)
  
  auc_list<-fingerprint$All$auc_primary
  
  index<-rownames(auc_list, subsystem)
  system_names<-auc_list[index,]$subsystem
  
  final_classifiers<-fingerprint$All$fingerprints_primary
  
  ex_subsystems_no_norm<-list()
  
  #Generate prediction values for each individual XGBoost subsystem model on reference assay
  for(x in 1:length(index)){
    newDat<-as.matrix(assay_metabolites[,final_classifiers[[as.numeric(index[x])]]$coefnames])
    colnames(newDat)<-final_classifiers[[as.numeric(index[x])]]$coefnames
    pred_tot<-predict(final_classifiers[[as.numeric(index[x])]], newdata = newDat,type = "prob", na.action = na.pass)[,1]
    ex_subsystems_no_norm[[system_names[x]]] <- pred_tot
  }
  
  ex_subsystems_no_norm<-data.frame(ex_subsystems_no_norm)
  
  score<-rowMeans(ex_subsystems_no_norm)
  score<-as.data.frame((score - min(score, na.rm = TRUE)) / (max(score, na.rm = TRUE) - min(score, na.rm = TRUE)))
  
  total_scores[[fingerprint_name]]<-score
  subsystem_scores[[fingerprint_name]]<-ex_subsystems_no_norm
  
  colnames(score)[1]<-fingerprint_name
  rownames(score)<-rownames(assay_metabolites)
  
  if(scale_type == "binary"){
    TUH84@meta.data[,fingerprint_name]<-score
    
    max <- 1
    min <- 0
    mid<-0.5
  }else{
    TUH84@meta.data[,fingerprint_name]<-scale(score)
    
    max <- max(scale(score), na.rm = T)
    min <- min(scale(score), na.rm = T) 
    mid<-0
  }
  
  sub_1<-subset(TUH84, subset = treatment == "CTL_D0" & GILT_status == "sensitive")
  p1 <- FeaturePlot(sub_1, features = fingerprint_name)+ scale_color_gradient2(low = "blue2",mid = "lavenderblush",high = "firebrick1", midpoint = mid, limits = c(min,max))+ ggtitle(paste0(sample_name," - Sensitive D0"))
  
  sub_2<-subset(TUH84, subset = treatment == "GILT_D7" & GILT_status == "persistent")
  p2 <- FeaturePlot(sub_2, features = fingerprint_name)+ scale_color_gradient2(low = "blue2",mid = "lavenderblush",high = "firebrick1", midpoint = mid, limits = c(min,max))+ ggtitle(paste0(sample_name," - Persistent D14"))
  
  sub_3<-subset(TUH84, subset = treatment == "GILT_D7" & GILT_status == "resistant")
  p3 <- FeaturePlot(sub_3, features = fingerprint_name)+ scale_color_gradient2(low = "blue2",mid = "lavenderblush",high = "firebrick1", midpoint = mid, limits = c(min,max))+ ggtitle(paste0(sample_name," - Resistant D14"))
  
  metadata <- data.frame(
    GILT_status = TUH84$GILT_status,
    treatment = TUH84$treatment,
    fingerprint_name = TUH84@meta.data[,fingerprint_name]
  )
  
  # Create the boxplot
  box_1 <- ggplot(metadata, aes(x = GILT_status, y = fingerprint_name, fill = GILT_status)) +
    geom_violin(width = 0.5) +
    stat_compare_means(method = "wilcox.test", 
                       comparisons = list(c("resistant", "sensitive"),c("persistent","sensitive")),
                       label = "p.signif") +
    labs(
      x = "GILT Status",
      y = "Fingerprint Value"
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 23,
      size = 4,
      fill = "white",
      color = "black"
    ) +
    stat_summary(
      fun = mean,
      geom = "line",
      aes(group = 1),
      color = "blue",
      linetype = "dashed"
    ) +
    scale_fill_manual(values = c("sensitive" = "#619CFF", persistent="lightgreen",resistant = "#F8766D")) +
    theme_minimal() +
    theme(legend.position = "none")
  
  box_2<-ggplot(metadata, aes(x = treatment, y = fingerprint_name, fill = treatment)) +
    geom_boxplot(width = 0.5) +
    stat_compare_means(method = "wilcox.test", 
                       comparisons = list(c("CTL_D0", "GILT_D7")),
                       label = "p.signif")+
    labs(
      x = "treatment",
      y = paste0(fingerprint_name,
                 " (Fingerprint Score)")
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      
      shape = 23,
      size = 4,
      fill = "white",
      color = "black"
    ) +stat_summary(
      fun = mean,
      geom = "line",
      aes(group = 1),
      color = "blue",
      linetype = "dashed"
    )+
    theme_minimal() +
    theme(legend.position = "none")
  
  plot_row <- plot_grid(p1,p2,p3,box_1, nrow =1)
  title <- ggdraw() + draw_label(fingerprint_name)
  print(plot_grid(title, plot_row, ncol = 1, rel_heights = c(0.1,1)))
}

dev.off()

total_scores<-as.data.frame(total_scores)
colnames(total_scores)<-fingerprint_names
rownames(total_scores)<-rownames(TUH84@meta.data)

total_scores_res<-total_scores[rownames(subset(TUH84@meta.data, GILT_status == "resistant")),]

cor_table<-cor(total_scores_res, method = "spearman")

ref_assay_name<-"Gilteritinib SC-TUH84 (resistant clusters)"

corrplot(cor_table, method = 'ellipse', type = 'upper', title = paste0("Correlations of Fingerprint Predictions\nReference Assay: ",ref_assay_name), tl.cex = 1, addCoef.col = "grey55",mar=c(0,0,2.5,0))

