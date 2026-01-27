library(tidyr)
library(ggplot2)
library(ggpubr)
library(caret)
library(data.table)
library(corrplot)

assay_metabolites_tot<-read.csv("beataml2_flux.csv", header = T, row.names = 1)
assay_metabolites<-predict(preProcess(assay_metabolites_tot,method = c("center", "scale")), assay_metabolites_tot)
assay_metabolites<- na.omit(assay_metabolites)

beataml_clin<-as.data.frame(as.matrix(fread("beataml2_clinical.csv"), rownames = 1))
beataml_ALT_up<-assay_metabolites[rownames(subset(beataml_clin, as.numeric(beataml_clin$ALT) > median(as.numeric(beataml_clin$ALT), na.rm = TRUE))),]
beataml_ALT_down<-assay_metabolites[rownames(subset(beataml_clin, as.numeric(beataml_clin$ALT) < median(as.numeric(beataml_clin$ALT), na.rm = TRUE))),]

beataml_met_high_response<-assay_metabolites[rownames(subset(exvivo, auc < median(auc))),]
beataml_met_low_response<-assay_metabolites[rownames(subset(exvivo, auc > median(auc))),]

assay_metabolites_tot<-read.csv("metabTable_kulkarni_2020.csv", header = T, row.names = 1)
assay_metabolites<-predict(preProcess(assay_metabolites_tot,method = c("center", "scale")), assay_metabolites_tot)
assay_metabolites<- na.omit(assay_metabolites)

kulkarni_meta<-as.data.frame(read.csv("meta_kulkarni_2020.csv", header = T, row.names = 1))
kulkarni_met<-assay_metabolites[rownames(subset(kulkarni_meta, treatment == "MET")),]
kulkarni_ctl<-assay_metabolites[rownames(subset(kulkarni_meta, treatment == "CTL")),]

ustinova_healthy<-read.csv("metabTable_metformin_healthy.csv", header = T, row.names = 1)
ustinova_healthy<-predict(preProcess(ustinova_healthy,method = c("center", "scale")), ustinova_healthy)
ustinova_healthy<- na.omit(ustinova_healthy)

ustinova_meta<-as.data.frame(read.csv("meta_metformin_healthy.csv", header = T, row.names = 1))
ustinova_healthy<-ustinova_healthy[rownames(subset(ustinova_meta, treatment != "MET_10_HOURS")),]

ustinova_diabetes<-read.csv("metabTable_ustinova_diabetes.csv", header = T, row.names = 1)
ustinova_diabetes<-predict(preProcess(ustinova_diabetes,method = c("center", "scale")), ustinova_diabetes)
ustinova_diabetes<- na.omit(ustinova_diabetes)

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

fingerprints<-list(MOLM14, HL60, U937, KG1a, breast, PDAC, melanoma_20mM, melanoma_40mM, cervical_11_days,cervical_24_days)
fingerprint_names<-list("MET_AML_MOLM14","MET_AML_HL60", "MET_AML_U937","MET_AML_KG1a","MET_breast","MET_PDAC","MET_melanoma_20mM","MET_melanoma_40mM","MET_cervical_11_days","MET_cervical_24_days")

assay<-beataml_met_low_response
ref_assay_name<-"BEATAML2 - Metformin Low Responder"

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
    newDat<-as.matrix(assay[,final_classifiers[[as.numeric(index[x])]]$coefnames])
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

corrplot(cor_table, method = 'ellipse', type = 'upper', title = paste0("Correlations of Fingerprint Predictions\nReference Assay: ",ref_assay_name), tl.cex = 0.9, addCoef.col = "grey",mar=c(0,0,2.5,0))

exp_name<-"BEATAML2_MET_HIGH_RESPONSE"
ctl_name<-"BEATAML2_MET_LOW_RESPONSE"

cor_df <- as.data.frame(as.table(cor_table_exp))

# Remove self-correlations and duplicate pairs
cor_df <- cor_df[cor_df$Var1 != cor_df$Var2, ]
cor_df <- cor_df[!duplicated(t(apply(cor_df[,1:2], 1, sort))), ]

# Convert to list
cor_list_exp <- split(cor_df$Freq, paste(cor_df$Var1, cor_df$Var2, sep = "-"))
cor_list_exp<-t(as.data.frame(cor_list_exp))
colnames(cor_list_exp)[1]<-exp_name

cor_df <- as.data.frame(as.table(cor_table_ctl))

# Remove self-correlations and duplicate pairs
cor_df <- cor_df[cor_df$Var1 != cor_df$Var2, ]
cor_df <- cor_df[!duplicated(t(apply(cor_df[,1:2], 1, sort))), ]

# Convert to list
cor_list_ctl <- split(cor_df$Freq, paste(cor_df$Var1, cor_df$Var2, sep = "-"))
cor_list_ctl<-t(as.data.frame(cor_list_ctl))
colnames(cor_list_ctl)[1]<-ctl_name

cor_list_tot<-data.frame(cor_list_exp,cor_list_ctl)

cor_list_tot$row_id <- seq_len(nrow(cor_list_tot))

df_long <- pivot_longer(
  cor_list_tot,
  cols = c(exp_name, ctl_name),
  names_to = "Stratification_Group",
  values_to = "Pearson_Correlation"
)

df_long$Stratification_Group<-factor(df_long$Stratification_Group, levels = c(ctl_name, exp_name))

ggplot(df_long, aes(x = Stratification_Group, y = Pearson_Correlation)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_point(aes(group = row_id), position = position_jitter(width = 0.1), size = 2) +
  geom_line(aes(group = row_id), color = "blue", alpha = 0.3) +
  stat_compare_means(method = "t.test", 
                     comparisons = list(c(exp_name, ctl_name)),
                     label = "p.signif")+ggtitle("Correlations of Cancer Metformin Fingerprint Models \nBetween BEATAML2 High and Low Exvivo Response to Metformin")


df_pairs <- combn(names(subsystem_scores_low), 2, simplify = FALSE)

correlation_list <- list()

for(pair in df_pairs) {
  df1 <- subsystem_scores_low[[pair[1]]]
  df2 <- subsystem_scores_low[[pair[2]]]
  
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

#mean_by_name<-mean_by_name[1:inflection_point_val,]

ggplot(diagonal_df, aes(x = row, y = average)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))+ggtitle("Average Maximum Correlation of Fingerprint Subsystem Features Across All Metformin Assays")

merged_df <- merge(diagonal_df_high, diagonal_df_low, by = "row", suffixes = c("_high_response", "_low_response"))

merged_df$t_stat <- (merged_df$average_high_response - merged_df$average_low_response) /
  sqrt((merged_df$std_dev_high_response^2 / merged_df$count_high_response) + (merged_df$std_dev_low_response^2 / merged_df$count_low_response))

numerator <- (merged_df$std_dev_high_response^2 / merged_df$count_high_response + merged_df$std_dev_low_response^2 / merged_df$count_low_response)^2
denominator <- ((merged_df$std_dev_high_response^2 / merged_df$count_high_response)^2 / (merged_df$count_high_response - 1)) +
  ((merged_df$std_dev_low_response^2 / merged_df$count_low_response)^2 / (merged_df$count_low_response - 1))

merged_df$df <- numerator / denominator

merged_df$p_value <- 2 * pt(-abs(merged_df$t_stat), df = merged_df$df)

merged_df$abs_diff <- abs(merged_df$average_high_response - merged_df$average_low_response)

result_high_response<-subset(merged_df, p_value < 0.05 & t_stat > 0 & average_high_response > 0.4 & abs_diff > 0.1)
result_low_response<-subset(merged_df, p_value < 0.05 & t_stat < 0 & average_low_response > 0.4 & abs_diff > 0.1)

similar_high<-subset(merged_df, abs_diff < 0.05 & average_low_response > 0.4 & p_value > 0.05)

metabolic_data <- result_low_response %>%
  mutate(diff = average_high_response - average_low_response) %>%
  arrange(desc(diff)) %>% 
  mutate(row = factor(row, levels = unique(row)))

plot_data <- metabolic_data %>%
  pivot_longer(cols = c(average_high_response, average_low_response),
               names_to = "response_type",
               values_to = "average_response") %>%
  mutate(std_dev = ifelse(response_type == "average_high_response",
                          std_dev_high_response,
                          std_dev_low_response))

ggplot(plot_data, aes(y = row, x = average_response, fill = response_type)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D"),
                    labels = c("High Response", "Low Response")) +
  labs(title = "Metabolic Subsystem Model Correlations Significant in BEATAML2 High Response to Metformin (Ex Vivo)",
       x = "Metabolic Subsystem Model",
       y = "Average Pearson Correlation",
       fill = "Response Type") +
  theme_minimal() +
  theme(axis.text.y = element_text(hjust = 1, face = "bold", size = 11),
        plot.title = element_text(hjust = 0.5, face = "bold"))

reactMeta<-read.csv("human_reaction_meta.csv", header = T, row.names = 1)
reactMetaMouse<-read.csv("mouse_reaction_meta.csv", header = T, row.names = 1)
reactMeta_filt<-reactMeta[intersect(rownames(reactMeta), rownames(reactMetaMouse)),]

subsystem_name<-"Vitamin.E.metabolism"

ref_assay_name_1<-"BEATAML2 - Metformin Sensitive Response"
ref_assay_name_2<-"BEATAML2 - Metformin Resistant Response"

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
  if(subsystem_name %in% colnames(subsystem_scores_high[[i]])){
    react_cors_list[[fingerprint_names[[i]]]]<-cor(beataml_met_low_response[,unique(react_list)], subsystem_scores_low[[i]][,subsystem_name], method = "pearson")
  }
}

react_cors_df_1<-as.data.frame(react_cors_list)
react_cors_df_1<-na.omit(react_cors_df_1)

rownames(react_cors_df_1)<-reactMeta_filt[rownames(react_cors_df_1),"RECON3D"]

for(i in 1:dim(total_scores)[2]){
  if(subsystem_name %in% colnames(subsystem_scores_low[[i]])){
    react_cors_list[[fingerprint_names[[i]]]]<-cor(beataml_met_high_response[,unique(react_list)], subsystem_scores_high[[i]][,subsystem_name], method = "pearson")
  }
}

react_cors_df_2<-as.data.frame(react_cors_list)
react_cors_df_2<-na.omit(react_cors_df_2)

rownames(react_cors_df_2)<-reactMeta_filt[rownames(react_cors_df_2),"RECON3D"]

react_cors_df_1<-react_cors_df_1[intersect(rownames(react_cors_df_1),rownames(react_cors_df_2)),]
react_cors_df_2<-react_cors_df_2[intersect(rownames(react_cors_df_1),rownames(react_cors_df_2)),]

react_cors_df_delta<-react_cors_df_2-react_cors_df_1

heatmap<-heatmaply(react_cors_df_1, scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "red", high = "blue", midpoint = 0, limits = c(min(react_cors_df_1), max(react_cors_df))), main = paste0("Correlation of Effects Among Common Fingerprint Features - "))

rows<-rev(heatmap$x$layout$yaxis2$ticktext)
columns<-heatmap$x$layout$xaxis$ticktext

p <- heatmaply(
  react_cors_df_1[rows, columns],
  scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
    low = "red", high = "blue", midpoint = 0, limits = c(-1, 1)
  ),
  main = paste0(
    "<span style='font-size:12px;'>Correlation of Subsystem Reaction Activity and Fingerprint Prediction - ",
    subsystem_name, "<br>",
    "Reference Assay: ", ref_assay_name_2, "</span>"
  ),
  Rowv = FALSE,
  Colv = FALSE,
  margins = c(50, 50, 130, 50)
)

# Modify the Y-axis tick label font size
p1 <- plotly::layout(p, yaxis = list(tickfont = list(size = 10,color = "black")))  # 👈 Set Y-axis label font size here

p1

p <- heatmaply(
  react_cors_df_2[rows, columns],
  scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
    low = "red", high = "blue", midpoint = 0, limits = c(-1, 1)
  ),
  main = paste0(
    "<span style='font-size:12px;'>Correlation of Subsystem Reaction Activity and Fingerprint Prediction - ",
    subsystem_name, "<br>",
    "Reference Assay: ", ref_assay_name_1, "</span>"
  ),
  Rowv = FALSE,
  Colv = FALSE,
  margins = c(50, 50, 130, 50)
)

# Modify the Y-axis tick label font size
p2 <- plotly::layout(p, yaxis = list(tickfont = list(size = 10, color = "black")))  # 👈 Set Y-axis label font size here

p2

gc()

p <- heatmaply(
  react_cors_df_delta[rows, columns],
  scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
    low = "red", high = "blue", midpoint = 0, limits = c(min(react_cors_df_delta), max(react_cors_df_delta))
  ),
  main = paste0(
    "<span style='font-size:12px;'>Correlation of Subsystem Reaction Activity and Fingerprint Prediction - ",
    subsystem_name, "<br>",
    "Reference Assay: Δ Metformin Sensitive - Metformin Resistant", "</span>"
  ),
  Rowv = FALSE,
  Colv = FALSE,
  margins = c(50, 50, 130, 50)
)

p3 <- plotly::layout(p, yaxis = list(tickfont = list(size = 7, color = "black")))  # 👈 Set Y-axis label font size here

p3

# Modify the Y-axis tick label font size

gc()

all_values_flt3<-as.data.frame(all_values)
all_values_flt3$subsystem<-names(all_values)
colnames(all_values_flt3)[1]<-"avg_correlation"
all_values_flt3$stratification<-"BEATAML2_FLT3_ITD_POS"

all_values_wt<-as.data.frame(all_values)
all_values_wt$subsystem<-names(all_values)
colnames(all_values_wt)[1]<-"avg_correlation"
all_values_wt$stratification<-"BEATAML2_FLT3_ITD_WT"

all_values_pos_wt<-rbind(all_values_flt3,all_values_wt)

t_test_results <- all_values_pos_wt %>%
  group_by(subsystem) %>%
  summarise(
    pos_values = list(avg_correlation[stratification == "BEATAML2_FLT3_ITD_POS"]),
    wt_values = list(avg_correlation[stratification == "BEATAML2_FLT3_ITD_WT"])
  ) %>%
  rowwise() %>%
  mutate(
    t_test = list(t.test(pos_values, wt_values, var.equal = FALSE)),
    t_statistic = t_test$statistic,
    p_value = t_test$p.value,
    pos_mean = mean(unlist(pos_values)),
    wt_mean = mean(unlist(wt_values)),
    n_pos = length(unlist(pos_values)),
    n_wt = length(unlist(wt_values))
  ) %>%
  select(subsystem, t_statistic, p_value, pos_mean, wt_mean, n_pos, n_wt) %>%
  arrange(p_value)

sig_t_test_results_up<-subset(t_test_results, t_statistic > 0 & p_value < 0.05)
sig_t_test_results_down<-subset(t_test_results, t_statistic < 0 & p_value < 0.05)

plot_data_up<-subset(all_values_pos_wt, subsystem %in% sig_t_test_results_up$subsystem)
plot_data_down<-subset(all_values_pos_wt, subsystem %in% sig_t_test_results_down$subsystem)

plot_data <- plot_data_up %>%
  group_by(subsystem, stratification) %>%
  summarise(mean_corr = mean(avg_correlation, na.rm = TRUE), .groups = 'drop')

plot_data$subsystem <- factor(plot_data$subsystem, levels = sig_t_test_results_up$subsystem)

ggplot(plot_data, aes(x = subsystem, y = mean_corr, fill = stratification)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("BEATAML2_FLT3_ITD_POS" = "blue", "BEATAML2_FLT3_ITD_WT" = "red")) +
  labs(
    title = "Average Correlation by Subsystem Significant in BEATAML2 FLT3-ITD POS Patients",
    x = "Subsystem",
    y = "Average Correlation",
    fill = "Stratification"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

plot_data <- plot_data_down %>%
  group_by(subsystem, stratification) %>%
  summarise(mean_corr = mean(avg_correlation, na.rm = TRUE), .groups = 'drop')

plot_data$subsystem <- factor(plot_data$subsystem, levels = sig_t_test_results_down$subsystem)

ggplot(plot_data, aes(x = subsystem, y = mean_corr, fill = stratification)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("BEATAML2_FLT3_ITD_POS" = "blue", "BEATAML2_FLT3_ITD_WT" = "red")) +
  labs(
    title = "Average Correlation by Subsystem Significant in BEATAML2 FLT3-ITD WT Patients",
    x = "Subsystem",
    y = "Average Correlation",
    fill = "Stratification"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

