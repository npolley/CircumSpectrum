library(ggplot2)
library(patchwork)
library(corrplot)
library(heatmaply)
library(caret)
library(inflection)
library(data.table)

fingerprint<-readRDS("day_28_.rds")
fingerprint_data<-readRDS("day_28.rds")

features<-make.names(fingerprint$All$subsystems$subsystem)

classifier_final <- fingerprint
#vals_subsystem <- fingerprint$All$self_score_primary
perf_subsystem <- fingerprint$All$auc_primary
perf_subsystem$names<-make.names(perf_subsystem$subsystem)
perf_subsystem<-subset(perf_subsystem, names %in% features)
auc_value <- fingerprint$All$auc_final
assay_name <- fingerprint$All$assay

model_fingerprints<-fingerprint$All$fingerprints_primary[as.numeric(perf_subsystem$index)]

reactions<-c()

for(i in 1: length(model_fingerprints)){
  reactions<-c(reactions, model_fingerprints[[i]]$coefnames)
}

subsystems_final<-fingerprint_data[[1]]
subsystems_final<-subset(subsystems_final, metabolite %in% reactions)

stats_fin <- data.frame(index=perf_subsystem$index, subsystem=perf_subsystem$subsystem, loocv_test_auc=as.numeric(perf_subsystem$loocv_test_auc))
stats_fin <- subset(stats_fin, loocv_test_auc > 0.5)
#stats_fin <- data.frame(stats_fin, shaps=shaps[colnames(vals_subsystem)[as.numeric(stats_fin$index)+1],"score"])
#stats_fin$shaps <- ifelse(is.na(stats_fin$shaps), min(stats_fin$shaps, na.rm = TRUE), stats_fin$shaps)
#stats_fin <- rbind(data.frame(index=0, subsystem="FINAL MODEL", loocv_test_auc=auc_value, shaps=NA), stats_fin)
rownames(stats_fin)<-stats_fin$subsystem

subsystems_final<-subset(subsystems_final, subsystem %in% stats_fin$subsystem)

plot1 <- ggplot(subsystems_final, aes(x = as.numeric(t_welch_statistic), y = reorder(subsystem, -as.numeric(t_welch_statistic)), color = sig, fill = sig)) +
  geom_point(size = 1, position = position_dodge(width = 0.5)) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(face = "bold"),
    legend.position = "none"  # Remove the legend
  ) +
  geom_vline(xintercept = 1.3, size = .05) +
  geom_vline(xintercept = -1.3, size = .05) +
  ylab("subsystem") +
  xlab(expression(atop(paste("← Deceased", "           ", "Survived →"), 
                       "log(pval) with sign of t-statistic"))) +
  ggtitle(paste0("Significant Reactions Across All Subsystems")) +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "mm"),
    plot.title = element_text(size = 12)  # Set title font size here
  )

gb <- ggplot_build(plot1)

# Extract y-axis labels (top to bottom)
y_labels <- gb$layout$panel_params[[1]]$y$get_labels()
stats_fin<-stats_fin[y_labels,]
stats_fin$subsystem<-factor(stats_fin$subsystem, y_labels)
#shap_values <- classifier_final$All$shaps_exp
#shap_values<-data.frame(colMeans(abs(shap_values)))
#colnames(shap_values)[1]<-"shap"
#shap_values<-data.frame(shap_values, spacer=0)
#shap_values <- shap_values[order(shap_values[["shap"]], decreasing = TRUE), ]

#imp<-varImp(classifier_final$All$fingerprint_exp, scale = TRUE)[[1]][,"Class1"]
#mean_imp<-imp/sum(imp)

#shaps <- data.frame(subsystem=rownames(shap_values), score=shap_values$shap)
#rownames(shaps) <- rownames(shap_values)


#inflection<-findiplist(stats_fin$loocv_test_auc, 1:length(stats_fin$loocv_test_auc), 0)[2,3]
#stats_fin<-subset(stats_fin, loocv_test_auc > inflection)

length_x <- dim(stats_fin)[1]-10

folds<-dim(fingerprint$All$train_primary)[1]

plot2 <- ggplot(stats_fin, aes(x = loocv_test_auc, y = subsystem, alpha = loocv_test_auc)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_cartesian(xlim = c(0.4, 1)) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "lines"),
    aspect.ratio = 16/9,
    legend.position = "bottom"
  ) +
  labs(
    title = paste0("Subsystem Model Performance"),
    x = paste0("LOOCV AUC (", folds, " Folds)")
  ) +
  scale_x_continuous(breaks = seq(0.4, 1, by = 0.1)) +
  geom_vline(xintercept = c(0.9, 0.8, 0.7, 0.6), linetype = "dashed", color = "gray50", size = 0.5) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"),plot.title = element_text(size = 12) ) +
  scale_alpha(guide = "none") 

plot1 + plot2 + plot_layout(ncol = 2)

write.csv(stats_fin, "healthy_stats_fin.csv")
write.csv(subsystems_final, "healthy_subsystems_final.csv")

#CORALIE Tier 1

kulkarni_sc<-readRDS("MET_CTL_strat_1.rds")

assay_metabolites_tot<-read.csv("beataml2_flux.csv", header = T, row.names = 1)
assay_metabolites<-predict(preProcess(assay_metabolites_tot,method = c("center", "scale")), assay_metabolites_tot)
assay_metabolites<- na.omit(assay_metabolites)

beataml_clin<-as.data.frame(as.matrix(fread("beataml2_clinical.csv"), rownames = 1))
beataml_flt3_pos_assay<-assay_metabolites[rownames(subset(beataml_clin, FLT3.ITD == "positive")),]
beataml_flt3_neg_assay<-assay_metabolites[rownames(subset(beataml_clin, FLT3.ITD == "negative")),]

kulkarni_sc_met<-subset(kulkarni_sc, Fingerprint_Stratifications == "MET (resistant)")

assay_metabolites_tot<-read.csv("MET_CTL_1_flux.csv", header = T, row.names = 1)
assay_metabolites_tot<-assay_metabolites_tot[rownames(kulkarni_sc_met@meta.data),]
assay_metabolites<-predict(preProcess(assay_metabolites_tot,method = c("center", "scale")), assay_metabolites_tot)
assay_metabolites<- na.omit(assay_metabolites)

kulkarni_2020<-readRDS("kulkarni_2020_.rds")
kulkarni_2017<-readRDS("kulkarni_2017_.rds")
ustinova_healthy_7_days<-readRDS("ustinova_healthy_7_days_.rds")
ustinova_healthy_10_hours<-readRDS("ustinova_healthy_10_hours_.rds")
ustinova_diabetes<-readRDS("ustinova_diabetes_.rds")

MOLM14<-readRDS("MET_MOLM14_.rds")
HL60<-readRDS("MET_HL60_.rds")
U937<-readRDS("MET_U937_.rds")
KG1a<-readRDS("MET_KG1a_.rds")
breast<-readRDS("MET_breast_.rds")
PDAC<-readRDS("MET_PDAC_.rds")
melanoma_20mM<-readRDS("MET_melanoma_20mM_.rds")
melanoma_40mM<-readRDS("MET_melanoma_40mM_.rds")
cervical_11_days<-readRDS("MET_cervical_11_days_.rds")
cervical_24_days<-readRDS("MET_cervical_24_days_.rds")


fingerprints<-list(kulkarni_2020, kulkarni_2017, ustinova_healthy_7_days, ustinova_healthy_10_hours, ustinova_diabetes)
fingerprint_names<-list("MET_kulkarni_2020","MET_kulkarni_2017", "MET_ustinova_healthy_7_days","MET_ustinova_healthy_10_hours","MET_ustinova_diabetes")

fingerprints<-list(MOLM14, HL60, U937, KG1a, breast, PDAC, melanoma_20mM, melanoma_40mM, cervical_11_days,cervical_24_days)
fingerprint_names<-list("MET_AML_MOLM14","MET_AML_HL60", "MET_AML_U937","MET_AML_KG1a","MET_breast","MET_PDAC","MET_melanoma_20mM","MET_melanoma_40mM","MET_cervical_11_days","MET_cervical_24_days")


total_scores<-list()
subsystem_scores<-list()

for(i in 1:length(fingerprints)){
  fingerprint<-fingerprints[[i]]
  fingerprint_name<-fingerprint_names[[i]]
  print(fingerprint_name)
  
  auc_list<-fingerprint$All$auc_primary
  rownames(auc_list)<-auc_list$index
  
  index<-auc_list$index
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
  score<-(score - min(score, na.rm = TRUE)) / (max(score, na.rm = TRUE) - min(score, na.rm = TRUE))
  
  total_scores[[fingerprint_name]]<-score
  subsystem_scores[[fingerprint_name]]<-ex_subsystems_no_norm
}

total_scores<-as.data.frame(total_scores)

cor_table<-cor(total_scores, method = "pearson")

ref_assay_name<-"Kulkarni 2020 - Metformin Murine Single-Cell"

corrplot(cor_table, method = 'ellipse', type = 'upper', title = paste0("Correlations of Fingerprint Predictions\nReference Assay: ",ref_assay_name), tl.cex = 0.9, addCoef.col = "grey",mar=c(0,0,2.5,0))

#CORALIE Tier 2

fingerprint_1<-"MET_MOLM14"
fingerprint_2<-"MET_U937"

cor_table_subsystem<-cor(subsystem_scores[[fingerprint_1]],subsystem_scores[[fingerprint_2]])

cor_table_subsystem<- cor_table_subsystem[rowSums(!is.na(cor_table_subsystem)) > 0, ]
cor_table_subsystem <- cor_table_subsystem[, colSums(!is.na(cor_table_subsystem)) > 0]

out_1<-ifelse(total_scores[[fingerprint_1]]>0.5,1,0)
sub_scores_1<-data.frame(out=out_1,subsystem_scores[[fingerprint_1]])

t_test_1<-data.frame(sapply(names(subset(sub_scores_1,out==1)[,-1]), function(col) t.test(subset(sub_scores_1,out==1)[,-1][[col]], subset(sub_scores_1,out==0)[,-1][[col]])$p.value))
colnames(t_test_1)[1]<-"pval"
t_test_1$pval<- p.adjust(t_test_1$pval, method = "bonferroni")

fingerprint_features_1<-rownames(subset(t_test_1, pval < 0.05))

out_2<-ifelse(total_scores[[fingerprint_2]]>0.5,1,0)
sub_scores_2<-data.frame(out=out_2,subsystem_scores[[fingerprint_2]])

t_test_2<-data.frame(sapply(names(subset(sub_scores_2,out==1)[,-1]), function(col) t.test(subset(sub_scores_2,out==1)[,-1][[col]], subset(sub_scores_2,out==0)[,-1][[col]])$p.value))
colnames(t_test_2)[1]<-"pval"
t_test_2$pval<- p.adjust(t_test_2$pval, method = "bonferroni")

fingerprint_features_2<-rownames(subset(t_test_2, pval < 0.05))

cor_table_subsystem<-cor_table_subsystem[fingerprint_features_1,fingerprint_features_2]

if(dim(cor_table_subsystem)[1]>dim(cor_table_subsystem)[2]){
  cor_table_subsystem<-t(cor_table_subsystem)
  temp<-fingerprint_2
  fingerprint_2<-fingerprint_1
  fingerprint_1<-temp
}

heatmap<-heatmaply(cor_table_subsystem, scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "red", high = "blue", midpoint = 0, limits = c(min(cor_table_subsystem), max(cor_table_subsystem))), main = paste0("Correlation of Effects Among Fingerprint Features - ", fingerprint_1, " and ", fingerprint_2),xlab = fingerprint_1, ylab = fingerprint_2)
rows<-rev(heatmap$x$layout$yaxis2$ticktext)
columns<-heatmap$x$layout$xaxis$ticktext

heatmaply(cor_table_subsystem[rows,columns], scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "red", high = "blue", midpoint = 0, limits = c(min(cor_table_subsystem), max(cor_table_subsystem))), main = paste0("Correlation of Effects Among Fingerprint Features","\n", fingerprint_1, " and ", fingerprint_2, "\n", "Reference Assay: ", ref_assay_name),xlab = fingerprint_1, ylab = fingerprint_2,Rowv = FALSE, Colv = FALSE, margins = c(50, 50, 130, 50))

rowSum_max<-as.data.frame(apply(cor_table_subsystem, 1, max))
colnames(rowSum_max)[1]<-"cor_max"
colSum_max<-as.data.frame(apply(cor_table_subsystem, 2, max))
colnames(colSum_max)[1]<-"cor_max"

rowSum_min<-as.data.frame(rowSums(cor_table_subsystem))
colnames(rowSum_min)[1]<-"cor_min"
colSum_min<-as.data.frame(colSums(cor_table_subsystem))
colnames(colSum_min)[1]<-"cor_min"

#shared subsystems (two fingerprints)

rowNames<-rownames(subset(rowSum_max, cor_max > quantile(rowSum_max$cor_max, 0.75)))
colNames<-rownames(subset(colSum_max, cor_max > quantile(rowSum_max$cor_max, 0.75)))

common_cor_table<-cor_table_subsystem[rowNames,colNames]

# 1. Get absolute correlations for optimal pairing
abs_cor <- abs(common_cor_table[rows, columns])

# 2. Find optimal column order to maximize diagonal correlations
assignment <- solve_LSAP(abs_cor, maximum = TRUE)
column_order <- as.integer(assignment)

# 3. Reorder matrix columns using optimal pairing
reordered_cor <- common_cor_table[rows, columns][, column_order]

heatmaply(
  reordered_cor,
  scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
    low = "red", 
    high = "blue",
    midpoint = 0,
    limits = c(min(cor_table_subsystem), max(cor_table_subsystem))
  ),
  main = paste0(
    "Correlation of Effects Among Common Fingerprint Features",
    "\n", fingerprint_1, " vs. ", fingerprint_2, 
    "\nReference Assay: ", ref_assay_name
  ),
  xlab = fingerprint_1,
  ylab = fingerprint_2,
  Rowv = FALSE, 
  Colv = FALSE,
  margins = c(50, 50, 130, 50),
  rect_gp = grid::gpar(col = "white", lwd = 2)  # White grid lines
)

#shared subsystems (all fingerprints)

df_pairs <- combn(names(subsystem_scores), 2, simplify = FALSE)

correlation_list <- list()

for(pair in df_pairs) {
  df1 <- subsystem_scores[[pair[1]]]
  df2 <- subsystem_scores[[pair[2]]]
  
  # Calculate correlation matrix
  cor_matrix <- cor(df1, df2, method = "pearson")
  
  # Create descriptive name and store
  list_name <- paste0(pair[1], "_vs_", pair[2])
  correlation_list[[list_name]] <- cor_matrix
}

# Get all unique pairs of row and column names
all_pairs <- expand.grid(
  row = unique(unlist(lapply(correlation_list, rownames))),
  col = unique(unlist(lapply(correlation_list, colnames))),
  stringsAsFactors = FALSE
)

# Function to extract all values for a given pair across the list
get_pair_values <- function(pair, df_list) {
  values <- sapply(df_list, function(df) {
    if (pair$row %in% rownames(df) && pair$col %in% colnames(df)) {
      val <- df[pair$row, pair$col]
      if (!is.na(val)) return(val)
    }
    return(NA)
  })
  values <- values[!is.na(values)]
  return(values)
}

# For each pair, extract values and compute stats
pair_stats <- apply(all_pairs, 1, function(pair) {
  values <- get_pair_values(list(row = pair[1], col = pair[2]), correlation_list)
  avg <- if (length(values) > 0) mean(values) else NA
  count <- length(values)
  stddev <- if (length(values) > 1) sd(values) else NA
  data.frame(row = pair[1], col = pair[2], average = avg, count = count, std_dev = stddev)
})

# Combine results into a single data frame
result_df <- do.call(rbind, pair_stats)

diagonal_df <- result_df[result_df$row == result_df$col, ]
diagonal_df <- diagonal_df[order(-diagonal_df$average), ]
diagonal_df$row<-factor(diagonal_df$row,diagonal_df$row)

diagonal_df<-subset(diagonal_df, average > 0.55)

ggplot(diagonal_df, aes(x = row, y = average)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, face = "bold")
  ) +
  ggtitle("Average Correlation of Fingerprint Subsystem Features Across All Metformin Assays")

#diverging subsystems

rowNames<-rownames(subset(rowSum_min, cor_min < quantile(rowSum_min$cor_min, 0.1)))
colNames<-rownames(subset(colSum_min, cor_min < quantile(colSum_min$cor_min, 0.1)))

union<-union(rowNames, colNames)
inter<-intersect(intersect(union,rownames(rowSum_min)),intersect(union,rownames(colSum_min)))

div_cor_table<-cor_table_subsystem[inter,inter]

neg_diag <- diag(div_cor_table) < 0

div_cor_table <- div_cor_table[neg_diag, neg_diag]

heatmap<-heatmaply(div_cor_table, scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "red", high = "blue", midpoint = 0, limits = c(min(cor_table_subsystem), max(cor_table_subsystem))), main = paste0("Correlation of Effects Among Common Fingerprint Features - ", fingerprint_1, " and ", fingerprint_2),xlab = fingerprint_1, ylab = fingerprint_2)
rows<-rev(heatmap$x$layout$yaxis2$ticktext)
columns<-heatmap$x$layout$xaxis$ticktext

heatmaply(div_cor_table[rows,columns], scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "red", high = "blue", midpoint = 0, limits = c(min(cor_table_subsystem), max(cor_table_subsystem))), main = paste0("Correlation of Effects Among Diverging Fingerprint Features","\n", fingerprint_1, " and ", fingerprint_2, "\n", "Reference Assay: ", ref_assay_name),xlab = fingerprint_1, ylab = fingerprint_2,Rowv = FALSE, Colv = FALSE, margins = c(50, 50, 130, 50))


#CORALIE Tier 3 (all fingerprints)
reactMeta<-read.csv("human_reaction_meta.csv", header = T, row.names = 1)
reactMetaMouse<-read.csv("mouse_reaction_meta.csv", header = T, row.names = 1)
reactMeta_filt<-reactMeta[intersect(rownames(reactMeta), rownames(reactMetaMouse)),]

subsystem_name<-"Bile.acid.biosynthesis"

react_list<-c()
react_cors_list<-list()

for(i in 1:length(fingerprints)){
  fingerprint<-fingerprints[[i]]
  auc_list<-fingerprint$All$auc_primary
  auc_list$subsystem<-make.names(auc_list$subsystem)
  if(subsystem_name %in% auc_list$subsystem){
  model_num<-as.numeric(unlist(subset(auc_list, subsystem == subsystem_name)$index))
  model<-fingerprint$All$fingerprints_primary[[model_num]]
  features<-model$coefnames
  react_list<-c(react_list, features)}
}

for(i in 1:dim(total_scores)[2]){
  if(subsystem_name %in% colnames(subsystem_scores[[i]])){
  react_cors_list[[fingerprint_names[[i]]]]<-cor(assay_metabolites[,unique(react_list)], subsystem_scores[[i]][,subsystem_name], method = "spearman")
  }
}

react_cors_df<-as.data.frame(react_cors_list)
react_cors_df<-na.omit(react_cors_df)

rownames(react_cors_df)<-reactMeta_filt[rownames(react_cors_df),"RECON3D"]

heatmap<-heatmaply(react_cors_df, scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "red", high = "blue", midpoint = 0, limits = c(min(react_cors_df), max(react_cors_df))), main = paste0("Correlation of Effects Among Common Fingerprint Features - "))

rows<-rev(heatmap$x$layout$yaxis2$ticktext)
columns<-heatmap$x$layout$xaxis$ticktext

p<-heatmaply(
  react_cors_df[rows, columns],
  scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
    low = "red", high = "blue", midpoint = 0, limits = c(min(react_cors_df), max(react_cors_df))
  ),
  main = paste0(
    "<span style='font-size:12px;'>Correlation of Subsystem Reaction Activity and Fingerprint Prediction - ",
    subsystem_name, "<br>",
    "Reference Assay: ", ref_assay_name, "</span>"
  ),
  Rowv = FALSE,
  Colv = FALSE,
  margins = c(50, 50, 130, 50)
)

p <- plotly::layout(p, yaxis = list(tickfont = list(size = 9, color = "black"))) 

p


react_cors_df_means<-as.data.frame(rowMeans(react_cors_df))

reactMeta_select<-reactMeta_filt[reactMeta_filt$RECON3D %in% rownames(react_cors_df),]
reactMeta_select<-data.frame(correlation=react_cors_df_means$`rowMeans(react_cors_df)`,reactMeta_select)

write.csv(reactMeta_select, "tca.csv")

final_reacts<-subset()
#write.csv(subset(subsystems_final, subsystem == "Androgen metabolism"),"androgen_metabolism.csv")
