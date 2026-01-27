library(Seurat)
library(ggplot2)
library(ggpubr)
library(cowplot)

#AUC Boxplot
auc_all<-read.csv("auc_all.csv", header = T, row.names = 1)

ggplot(auc_all, aes(x = reorder(model_type, AUC, FUN = mean), y = AUC, fill = model_type)) +
  geom_boxplot(width = 0.2) +
  
  labs(title = "Performance Distribution by Model Type",
       subtitle = "Boxplot of AUC Values (n = 98 LOOCV Iterations)",
       x = "Model Type",
       y = "Area Under the Curve (AUC)") +
  theme(
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.position = "none",
    plot.margin = unit(c(0.1, 0.5, 0.5, 0.5), "cm") # Adjusts plot margins
  ) +
  theme(aspect.ratio = 1.2)

#Single-Cell
kulkarni_sc<-readRDS("MET_CTL_strat_1.rds")

kulkarni_sc_flux<-read.csv("MET_CTL_1_flux.csv", header = T, row.names = 1)
kulkarni_sc_flux<-kulkarni_sc_flux[rownames(kulkarni_sc@meta.data),]
kulkarni_sc_flux<-predict(preProcess(kulkarni_sc_flux,method = c("center", "scale")), kulkarni_sc_flux)
kulkarni_sc_flux<- na.omit(kulkarni_sc_flux)

kulkarni_sc<-subset(kulkarni_sc, cells = rownames(kulkarni_sc_flux))
kulkarni_sc$Fingerprint_Stratifications<-ifelse(kulkarni_sc$Fingerprint_Stratifications == "CTL (sensitive)", "CTL (declining)", "MET (increasing)")

kulkarni_young_sc_gene<-read.csv("kulkarni_young_strat_gene.csv", header = T, row.names = 1)
kulkarni_young_sc_gene_norm<-predict(preProcess(kulkarni_young_sc_gene,method = c("center", "scale")), kulkarni_young_sc_gene)

kulkarni_young_sc_kegg<-read.csv("kulkarni_young_strat_kegg.csv", header = T, row.names = 1)


flux_fingerprint_model<-readRDS("kulkarni_young_flux_fingerprint_.rds")
KEGG_fingerprint_model<-readRDS("kulkarni_young_KEGG_fingerprint_.rds")
xgboost_flux_model<-readRDS("kulkarni_young_xgboost_flux.rds")
xgboost_kegg_model<-readRDS("kulkarni_young_xgboost_kegg.rds")
xgboost_gene_model<-readRDS()
top_genes<-read.csv()
bottom_genes<-read.csv()

#iMAT Flux Fingerprint

fingerprint<-flux_fingerprint_model$All
final_classifiers<-fingerprint$fingerprints_primary

auc_list<-fingerprint$auc_primary

systems<-list()
sub_model_perf<-list()
hidden_units<-list()

for(x in auc_list$index){
  units<-final_classifiers[[as.numeric(x)]]$finalModel$snnsObject$getUnitDefinitions()
  hidden_units[[auc_list[x,]$subsystem]]<-length(units$unitNo[units$type == "UNIT_HIDDEN"])
}

auc_list<-cbind(auc_list, data.frame(hidden_units=t(as.data.frame(hidden_units))[make.names(auc_list$subsystem),]))

rownames(auc_list)<-auc_list$index

#auc_list<-subset(auc_list, hidden_units > quantile(auc_list$hidden_units)[2])

index<-auc_list$index
system_names<-auc_list[index,]$subsystem

ex_subsystems_no_norm<-list()

#Generate prediction values for each individual XGBoost subsystem model on reference assay
for(x in 1:length(index)){
  newDat<-as.matrix(kulkarni_sc_flux[,final_classifiers[[as.numeric(index[x])]]$coefnames])
  colnames(newDat)<-final_classifiers[[as.numeric(index[x])]]$coefnames
  pred_tot<-predict(final_classifiers[[as.numeric(index[x])]], newdata = newDat,type = "prob", na.action = na.pass)[,1]
  ex_subsystems_no_norm[[system_names[x]]] <- pred_tot
  sub_model_perf[[system_names[x]]]<-round(auc(kulkarni_sc$Fingerprint_Stratifications, pred_tot),2)
}

ex_subsystems_no_norm<-data.frame(ex_subsystems_no_norm)

normalize01 <- function(v) {
  rng <- range(v, na.rm = TRUE)
  (v - rng[1]) / (rng[2] - rng[1])
}

ex_subsystems_no_norm <- as.data.frame(lapply(ex_subsystems_no_norm, normalize01))

score<-rowMeans(ex_subsystems_no_norm)
score<-(score - min(score, na.rm = TRUE)) / (max(score, na.rm = TRUE) - min(score, na.rm = TRUE))

kulkarni_sc$iMAT_Flux_Fingerprint<-score
kulkarni_sc$iMAT_Flux_Fingerprint<-as.numeric(kulkarni_sc$iMAT_Flux_Fingerprint)

sub_1<-subset(kulkarni_sc, Fingerprint_Stratifications == "CTL (declining)")
p1 <- FeaturePlot(sub_1, features = "iMAT_Flux_Fingerprint")+ scale_color_gradient2(low = "blue2",mid = "lavenderblush",high = "firebrick1", midpoint = 0.5, limits = c(0,1))+ ggtitle("CTL (declining)")+NoLegend()

sub_2<-subset(kulkarni_sc, Fingerprint_Stratifications == "MET (increasing)")
p2 <- FeaturePlot(sub_2, features = "iMAT_Flux_Fingerprint")+ scale_color_gradient2(low = "blue2",mid = "lavenderblush",high = "firebrick1", midpoint = 0.5, limits = c(0,1))+ ggtitle("MET (increasing)")

metadata <- data.frame(
  Fingerprint_Stratifications = kulkarni_sc$Fingerprint_Stratifications,
  iMAT_Flux_Fingerprint = kulkarni_sc@meta.data[,"iMAT_Flux_Fingerprint"]
)


# Create the boxplot
box_1 <- ggplot(metadata, aes(x = Fingerprint_Stratifications, y = iMAT_Flux_Fingerprint, fill = Fingerprint_Stratifications)) +
  geom_violin(width = 0.5) +
  # Add shaded regions and dashed line
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", size = 0.3) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = Inf, alpha = 0.1, fill = "red") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.5, alpha = 0.1, fill = "blue") +
  stat_compare_means(method = "t.test", 
                     comparisons = list(c("MET (increasing)", "CTL (declining)")),
                     label = "p.signif") +
  labs(
    x = "Metformin Response",
    y = "Prediction Value",
    title = paste0("AUC = ", round(auc(kulkarni_sc$Fingerprint_Stratifications, kulkarni_sc$iMAT_Flux_Fingerprint),2))
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
  scale_fill_manual(values = c("CTL (declining)" = "#619CFF", "MET (increasing)" = "#F8766D")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold"),
    axis.title.y = element_text(size = 11),
    plot.margin = margin(5, 60, 5, 5)
  ) +
  annotate("text", x = 2.5, y = 0.375, label = "CTL-prediction →", size = 3,
           color = "#619CFF", fontface = "bold", hjust = 0, angle = 270) +
  annotate("text", x = 2.5, y = 0.9, label = "← MET-prediction", size = 3,
           color = "#F8766D", fontface = "bold", hjust = 0, angle = 270) +
  coord_cartesian(ylim = c(0,1.05),clip = "off")


plot_row <- plot_grid(p1,p2,box_1, nrow =1)
title <- ggdraw() + draw_label("SSGSEA Signature (one-sided)")
print(plot_grid(title, plot_row, ncol = 1, rel_heights = c(0.1,1)))

#KEGG Pathway Fingerprint

fingerprint<-KEGG_fingerprint_model$All
fingerprint$auc_primary$subsystem<-make.names(fingerprint$auc_primary$subsystem)

index<-rownames(subset(fingerprint$auc_primary, subsystem %in% fingerprint$fingerprint_final$coefnames))
system_names<-fingerprint$auc_primary[index,]$subsystem

final_classifiers<-fingerprint$fingerprints_primary

ex_subsystems<-list()

for(x in 1:length(index)){
  pred_tot<-predict(final_classifiers[[as.numeric(index[x])]], newdata = as.matrix(kulkarni_young_sc_gene_norm[rownames(kulkarni_young_sc_gene_norm),final_classifiers[[as.numeric(index[x])]]$coefnames]),type = "prob", na.action = na.pass)[,1]
  ex_subsystems[[system_names[x]]] <- pred_tot
}

ex_subsystems<-data.frame(ex_subsystems)
rownames(ex_subsystems)<-rownames(kulkarni_young_sc_gene_norm)
ex_subsystems<-predict(preProcess(ex_subsystems,method = c("center", "scale")), ex_subsystems)

score<-as.data.frame(predict(fingerprint$fingerprint_final$finalModel, as.matrix(ex_subsystems[,fingerprint$fingerprint_final$finalModel$xNames])),type = "prob", na.action = na.pass)[1]

kulkarni_young_sc$KEGG_Fingerprint<-score

sub_1<-subset(kulkarni_young_sc, Fingerprint_Stratifications == "YOUNG_CTL (sensitive)")
p1 <- FeaturePlot(sub_1, features = "KEGG_Fingerprint")+ scale_color_gradient2(low = "blue2",mid = "lavenderblush",high = "firebrick1", midpoint = 0.5, limits = c(0,1))+ ggtitle("YOUNG_CTL (sensitive)")

sub_2<-subset(kulkarni_young_sc, Fingerprint_Stratifications == "MET (resistant)")
p2 <- FeaturePlot(sub_2, features = "KEGG_Fingerprint")+ scale_color_gradient2(low = "blue2",mid = "lavenderblush",high = "firebrick1", midpoint = 0.5, limits = c(0,1))+ ggtitle("MET (resistant)")

metadata <- data.frame(
  Fingerprint_Stratifications = kulkarni_young_sc$Fingerprint_Stratifications,
  KEGG_Fingerprint = kulkarni_young_sc@meta.data[,"KEGG_Fingerprint"]
)


# Create the boxplot
box_1 <- ggplot(metadata, aes(x = Fingerprint_Stratifications, y = KEGG_Fingerprint, fill = Fingerprint_Stratifications)) +
  geom_violin(width = 0.5) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("MET (resistant)","YOUNG_CTL (sensitive)")),
                     label = "p.signif") +
  labs(
    x = "Metformin Response",
    y = "Prediction Value"
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
  scale_fill_manual(values = c("YOUNG_CTL (sensitive)" = "#619CFF", "MET (resistant)" = "#F8766D")) +
  theme_minimal() +
  theme(legend.position = "none")

plot_row <- plot_grid(p1,p2,box_1, nrow =1)
title <- ggdraw() + draw_label("SSGSEA_KEGG_Fingerprint")
print(plot_grid(title, plot_row, ncol = 1, rel_heights = c(0.1,1)))

#XGBoost KEGG

score<-as.data.frame(predict(xgboost_kegg_model, as.matrix(kulkarni_young_sc_kegg[,xgboost_kegg_model$coefnames]),type = "prob"), na.action = na.pass)[1]
rownames(score)<-rownames(kulkarni_young_sc_kegg)

kulkarni_young_sc$XGBoost_kegg<-score

sub_1<-subset(kulkarni_young_sc, Fingerprint_Stratifications == "YOUNG_CTL (sensitive)")
p1 <- FeaturePlot(sub_1, features = "XGBoost_kegg")+ scale_color_gradient2(low = "blue2",mid = "lavenderblush",high = "firebrick1", midpoint = 0.5, limits = c(0,1))+ ggtitle("YOUNG_CTL (sensitive)")

sub_2<-subset(kulkarni_young_sc, Fingerprint_Stratifications == "MET (resistant)")
p2 <- FeaturePlot(sub_2, features = "XGBoost_kegg")+ scale_color_gradient2(low = "blue2",mid = "lavenderblush",high = "firebrick1", midpoint = 0.5, limits = c(0,1))+ ggtitle("MET (resistant)")

metadata <- data.frame(
  Fingerprint_Stratifications = kulkarni_sc$Fingerprint_Stratifications,
  XGBoost_flux = kulkarni_sc@meta.data[,"gene_signature"]
)


# Create the boxplot
box_1 <- ggplot(metadata, aes(x = Fingerprint_Stratifications, y = XGBoost_flux, fill = Fingerprint_Stratifications)) +
  geom_violin(width = 0.5) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("MET (resistant)","CTL (sensitive)")),
                     label = "p.signif") +
  labs(
    x = "Metformin Response",
    y = "Prediction Value"
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
  scale_fill_manual(values = c("YOUNG_CTL (sensitive)" = "#619CFF", "MET (resistant)" = "#F8766D")) +
  theme_minimal() +
  theme(legend.position = "none")

plot_row <- plot_grid(p1,p2,box_1, nrow =1)
title <- ggdraw() + draw_label("XGBoost_SSGSEA_KEGG")
print(plot_grid(title, plot_row, ncol = 1, rel_heights = c(0.1,1)))