library(caret)
library(Boruta)
library(pROC)
library(ROCR)
library(singscore)
library(Biobase)
library(DESeq2)
library(limma)
library(dplyr)
library(ggpubr)

fingerprint<-readRDS("ustinova_healthy_10_hours_.rds")

flux<-read.csv("kulkarni_sc_flux.csv", header = T, row.names = 1)
assay_metabolites_tot<-flux[rownames(kulkarni_sc@meta.data),]
assay_metabolites_tot<-predict(preProcess(assay_metabolites_tot,method = c("center", "scale")), assay_metabolites_tot)
assay_metabolites_tot<-na.omit(assay_metabolites_tot)
cells_use <- intersect(rownames(assay_metabolites_tot), colnames(kulkarni_sc))
cells_use <- cells_use[match(rownames(assay_metabolites_tot), cells_use, nomatch = 0)]
kulkarni_sc<-subset(kulkarni_sc, cells = cells_use)


test_vals<-ifelse(grepl("CTL_",rownames(assay_metabolites)),0,1)
test_vals<-factor(test_vals, levels = c(1,0))

subsystem_auc_train<-list()

for(i in 2:dim(train)[2]){
  # roc_obj <- roc(response = test_vals, predictor = ex_subsystems_no_norm[,i])
  # optimal_cutoff <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "youden")$threshold
  # subsystem_auc[[colnames(ex_subsystems_no_norm)[i]]]<-confusionMatrix(factor(test_vals, levels=c(0,1)),factor(ifelse(ex_subsystems_no_norm[,i]>median(ex_subsystems_no_norm[,i]),1,0)))$byClass[["Balanced Accuracy"]]
  subsystem_auc_train[[colnames(train)[i]]]<-t.test(subset(train, outer == "Class1")[,i],subset(train, outer == "Class0")[,i], alternative = "two.sided")[["p.value"]]
}



#Calculate AUC with LOOCV
test_auc_list<-list()

counter<-0

sampled_rows_list<-list()
test_vals_list<-list()

for(i in 1:50){
  rows<-sample(nrow(assay_metabolites), 200)
  sampled_rows_list[[i]] <- rows
  test_vals_list[[i]] <- rows
}


for(i in 1:length(sampled_rows_list)){
  
  print(counter)
  
  auc_list<-fingerprint$All$auc_primary
  
  index<-rownames(auc_list, subsystem)
  system_names<-auc_list[index,]$subsystem
  
  final_classifiers<-fingerprint$All$fingerprints_primary
  
  ex_subsystems_no_norm<-list()
  
  #Generate prediction values for each individual XGBoost subsystem model on reference assay
  for(x in 1:length(index)){
    newDat<-as.matrix(assay_metabolites[sampled_rows_list[[i]],final_classifiers[[as.numeric(index[x])]]$coefnames])
    colnames(newDat)<-final_classifiers[[as.numeric(index[x])]]$coefnames
    pred_tot<-predict(final_classifiers[[as.numeric(index[x])]], newdata = newDat,type = "prob", na.action = na.pass)[,1]
    ex_subsystems_no_norm[[system_names[x]]] <- pred_tot
  }
  
  ex_subsystems_no_norm<-data.frame(ex_subsystems_no_norm)
  
  score<-rowMeans(ex_subsystems_no_norm)
  score<-(score - min(score, na.rm = TRUE)) / (max(score, na.rm = TRUE) - min(score, na.rm = TRUE))
  
  test_auc<-pROC::auc(test_vals[test_vals_list[[i]]],unlist(score))
  print(test_auc)
  test_auc_list[[i]]<-test_auc
  
  counter<-counter+1
}

auc_imat_fingerprint<-data.frame(model_type="iMAT Flux Fingerprint", AUC=unlist(test_auc_list))

assay_genes<-data.frame(read.csv("MET_CTL_strat_1_gene.csv", header = T, row.names =1))
assay_genes<-data.frame(outer=assay_genes[,1],predict(preProcess(assay_genes,method = c("center", "scale")), assay_genes)[,-1])

test_auc_list<-list()

counter<-0

sampled_rows_list<-list()
test_vals_list<-list()

for(i in 1:50){
  rows<-sample(nrow(assay_genes), 200)
  sampled_rows_list[[i]] <- rows
  test_vals_list[[i]] <- rows
}


for(i in 1:length(sampled_rows_list)){
  
  print(counter)
  
  index<-rownames(auc_list, subsystem)
  system_names<-auc_list[index,]$subsystem
  
  final_classifiers<-models_list
  
  ex_subsystems_no_norm<-list()
  
  #Generate prediction values for each individual XGBoost subsystem model on reference assay
  for(x in 1:length(index)){
    pred_tot<-predict(final_classifiers[[as.numeric(index[x])]], newdata = assay_genes[sampled_rows_list[[i]],final_classifiers[[as.numeric(index[x])]]$coefnames],type = "prob", na.action = na.pass)[,1]
    ex_subsystems_no_norm[[system_names[x]]] <- pred_tot
  }
  
  ex_subsystems_no_norm<-data.frame(ex_subsystems_no_norm)
  
  score<-rowMeans(ex_subsystems_no_norm)
  score<-(score - min(score, na.rm = TRUE)) / (max(score, na.rm = TRUE) - min(score, na.rm = TRUE))
  
  test_auc<-pROC::auc(test_vals[test_vals_list[[i]]],unlist(score))
  print(test_auc)
  test_auc_list[[i]]<-test_auc
  
  counter<-counter+1
}

auc_kegg_fingerprint<-data.frame(model_type="KEGG Pathway Fingerprint", AUC=unlist(test_auc_list))

ex_subsystems_no_norm[is.na(ex_subsystems_no_norm)] <- 0

subsystem_auc<-list()

for(i in 1:dim(ex_subsystems)[2]){
  roc_obj <- pROC::roc(response = test_vals, predictor = ex_subsystems[,i])
  subsystem_auc[[colnames(ex_subsystems)[i]]]<-roc_obj$auc
}

subsystem_auc<-t(data.frame(subsystem_auc))
colnames(subsystem_auc)[1]<-"AUC"

metrics<-auc_list
rownames(metrics)<-auc_list[,1]
metrics<-metrics[,-1]
metrics_all<-data.frame(subsystem_auc, metrics[rownames(subsystem_auc),])

train_imp<-readRDS("All.rds")[["data"]]

reactions<-c()
model_num<-rownames(subset(auc_list))

ex_subsystems_no_norm<-list()
system_names<-auc_list[model_num,]$subsystem
for(x in 1:length(model_num)){
  newDat<-as.matrix(assay_metabolites[rownames(assay_metabolites),models_list[[as.numeric(model_num[x])]]$coefnames])
  colnames(newDat)<-models_list[[as.numeric(model_num[x])]]$coefnames
  pred_tot<-predict(models_list[[as.numeric(model_num[x])]], newdata = newDat,type = "prob", na.action = na.pass)[,1]
  ex_subsystems_no_norm[[system_names[x]]] <- pred_tot
}

ex_subsystems_no_norm<-data.frame(ex_subsystems_no_norm)
rownames(ex_subsystems_no_norm)<-rownames(assay_metabolites)
ex_subsystems<-predict(preProcess(ex_subsystems_no_norm,method = c("center", "scale")), ex_subsystems_no_norm)
ex_subsystems[is.na(ex_subsystems)] <- 0

boruta_res <- Boruta(y = factor(test_vals), x = ex_subsystems, pValue=0.001, doTrace = 0)
selected_features <- getSelectedAttributes(boruta_res, withTentative = FALSE)

subsystem_cors<-cor(ex_subsystems, method = "spearman")
diag(subsystem_cors) <- NA

cor_means<-colMeans(subsystem_cors, na.rm = TRUE)

threshold <- 0.4
high_corr_idx <- which(abs(subsystem_cors) > threshold, arr.ind = TRUE)

results <- data.frame(
  Var1 = rownames(subsystem_cors)[high_corr_idx[,1]],
  Var2 = colnames(subsystem_cors)[high_corr_idx[,2]],
  Correlation = subsystem_cors[high_corr_idx]
)

results <- results[as.numeric(high_corr_idx[,1]) < as.numeric(high_corr_idx[,2]), ]

systems<-unique(c(results$Var1,results$Var2))

model_num<-rownames(subset(auc_list))

for(i in 1:length(model_num)){
  reactions<-c(reactions,models_list[[as.numeric(model_num[i])]]$coefnames)
}

train_imp<-train_imp[,c("outer",reactions)]

# class1_data <- subset(train_imp, outer == 1)[,-1]
# class1_numeric <- as.data.frame(lapply(class1_data, as.numeric))
# class1 <- as.matrix(class1_numeric)
# 
# class0_data <- subset(train_imp, outer == 0)[,-1]
# class0_numeric <- as.data.frame(lapply(class0_data, as.numeric))
# class0 <- as.matrix(class0_numeric)
# 
# t_test<-col_t_welch(class1, class0)
# 
# reactions<-rownames(subset(t_test, statistic > 0))
# train_imp<-train_imp[,c("outer",reactions)]

train_imp$outer <- factor(ifelse(train_imp$outer==1, "Class1","Class0"), levels = c("Class1", "Class0"))

control <- trainControl(method = "cv",
                        number = 5,
                        summaryFunction = defaultSummary,
                        savePredictions = "final",
                        classProbs = TRUE)

classifier_exp <- caret::train(outer ~ ., data = train_imp,
                               method = "svmRadialSigma",
                               trControl = control,
                               metric = "AUC")

boruta_res <- Boruta(y = factor(train_imp[,1]), x = train_imp[,-1], pValue=0.001, doTrace = 0)
selected_features <- getSelectedAttributes(boruta_res, withTentative = FALSE)

train_imp<-as.data.frame(train_imp[,c("outer",selected_features)])

control <- trainControl(method = "cv",
                        number = 5,
                        summaryFunction = defaultSummary,
                        savePredictions = "final",
                        classProbs = TRUE)

classifier_final <- caret::train(outer ~ ., data = train_imp,
                                 method = "mlpML",
                                 trControl = control,
                                 metric = "AUC")

fingerprint<-list(fingerprint_final=classifier_final, fingerprint_exp=classifier_exp)

score_final<-as.data.frame(predict(fingerprint$fingerprint_final, newdata= as.data.frame(assay_metabolites[,fingerprint$fingerprint_final$coefnames]), type = "prob"))[,1]
score_exp<-as.data.frame(predict(fingerprint$fingerprint_exp, newdata= as.data.frame(assay_metabolites[,fingerprint$fingerprint_exp$coefnames]), type = "prob"))[,1]

pROC::roc(test_vals,as.numeric(score_final))
pROC::roc(test_vals,as.numeric(score_exp))


pdf(paste0("MET_CTL_2.pdf"),width = 11.8, height = 8.4)
for(i in 1:dim(ex_subsystems)[2]){
  name<-colnames(ex_subsystems)[i]
  
  score<-ex_subsystems[,name]
  
  kulkarni_sc$iMAT_Flux_Fingerprint<-score
  
  sub_1<-subset(kulkarni_sc, Fingerprint_Stratifications == "CTL (sensitive)")
  p1 <- FeaturePlot(sub_1, features = "iMAT_Flux_Fingerprint")+ scale_color_gradient2(low = "blue2",mid = "lavenderblush",high = "firebrick1", midpoint = median(score), limits = c(min(score),max(score)))+ ggtitle("CTL (sensitive)")
  
  sub_2<-subset(kulkarni_sc, Fingerprint_Stratifications == "MET (resistant)")
  p2 <- FeaturePlot(sub_2, features = "iMAT_Flux_Fingerprint")+ scale_color_gradient2(low = "blue2",mid = "lavenderblush",high = "firebrick1", midpoint = median(score), limits = c(min(score),max(score)))+ ggtitle("MET (resistant)")
  
  metadata <- data.frame(
    Fingerprint_Stratifications = kulkarni_sc$Fingerprint_Stratifications,
    iMAT_Flux_Fingerprint = kulkarni_sc@meta.data[,"iMAT_Flux_Fingerprint"]
  )
  
  
  # Create the boxplot
  box_1 <- ggplot(metadata, aes(x = Fingerprint_Stratifications, y = iMAT_Flux_Fingerprint, fill = Fingerprint_Stratifications)) +
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
    scale_fill_manual(values = c("CTL (sensitive)" = "#619CFF", "MET (resistant)" = "#F8766D")) +
    theme_minimal() +
    theme(legend.position = "none")
  
  plot_row <- plot_grid(p1,p2,box_1, nrow =1)
  title <- ggdraw() + draw_label(paste0(name," (t = ",t_table[i,"t.score"],") (AUC = ", auc(test_vals,unlist(score)),")"))
  print(plot_grid(title, plot_row, ncol = 1, rel_heights = c(0.1,1)))
}
dev.off()
#preds <- sapply(train_auc_list, function(x) x$pred)
#train_actuals <- sapply(train_auc_list, function(x) ifelse(x$actual=="Class0",0,1))

#train_auc <- auc(unlist(train_actuals), unlist(preds))
#print(paste0("Self AUC: ",train_auc))
print(paste0("External AUC: ",mean(unlist(test_auc_list))))

#self_out<-data.frame(model_type="iMAT_Flux_Fingerprint_SELF",AUC=train_auc)
test_out<-data.frame(model_type="iMAT_Flux_Fingerprint",AUC=unlist(test_auc_list))

#write.csv(self_out, "auc_score_flux_fingerprint_SELF.csv")
write.csv(test_out, "auc_score_flux_fingerprint_TEST.csv")

train<-data.frame(read.csv("ustinova_healthy_7_days_gene.csv", header = T, row.names =1))
meta_train<-as.data.frame(train$outer)
colnames(meta_train)[1]<-"outer"
meta_train$outer<-ifelse(meta_train$outer == 1, "MET", "CTL")
counts_train<-train[,-1]

assay_genes_tot<-read.csv("MET_CTL_strat_1_gene.csv", header = T, row.names = 1)
meta_test<-as.data.frame(assay_genes_tot$outer)
colnames(meta_test)[1]<-"outer"
meta_test$outer<-ifelse(meta_test$outer == 1, "MET", "CTL")
counts_test<-assay_genes_tot[,-1]

Train <- counts_train
Train_meta <- as.data.frame(meta_train)
colnames(Train_meta)[1]<-"outer"

design <- model.matrix(~ 0 + factor(Train_meta$outer))
colnames(design) <- c("outer_1", "outer_0")
fit <- lmFit(t(Train), design)
contrast <- makeContrasts(outer_1 - outer_0, levels=design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

top_genes <- topTable(fit2, number=Inf, adjust.method="bonferroni")

significant_genes_up <- top_genes %>%
  filter(P.Value < 0.01 & t > 0)

gene_signature_up <- rownames(significant_genes_up)

train_gene_expression_matrix<-rankGenes(Train)

train_scores <- simpleScore(t(train_gene_expression_matrix), upSet = gene_signature_up)

train_control <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)

model <- train(x = as.matrix(data.frame(score=scale(train_scores[1]),space=0)), y = factor(Train_meta[,1]),
               method = "glmnet",
               trControl = train_control,
               metric = "ROC")

#Calculate AUC with LOOCV
auc_list_test<-list()

sampled_rows_list<-list()
test_vals_list<-list()

for(i in 1:50){
  rows<-sample(nrow(assay_genes_tot), 200)
  sampled_rows_list[[i]] <- rows
  test_vals_list[[i]] <- rows
}

counter<-0

for(i in 1:length(sampled_rows_list)){
  print(counter)
  
  test_gene_expression_matrix<-rankGenes(counts_test[sampled_rows_list[[i]],])
  test_scores <- simpleScore(t(test_gene_expression_matrix), upSet = gene_signature_up)
  
  score_label<-as.numeric(predict(model, data.frame(score=test_scores[1],space=0), type = "prob")[,"MET"])
  print(pROC::auc(meta_test$outer[test_vals_list[[i]]], score_label))
  auc_list_test[[i]]<-pROC::auc(meta_test$outer[test_vals_list[[i]]], score_label)
  
  counter<-counter+1
}

auc_signature_one_sided<-data.frame(model_type="SSGSEA Signature (one-sided)", AUC=unlist(auc_list_test))

significant_genes_down <- top_genes %>%
  filter(P.Value < 0.01 & t < 0)

gene_signature_down <- rownames(significant_genes_down)

train_scores_down <- simpleScore(t(train_gene_expression_matrix), upSet = gene_signature_down)
train_control <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)
model <- train(x = as.matrix(data.frame(score=scale(train_scores[1]),down_score=scale(train_scores_down[1]))), y = factor(Train_meta[,1]),
               method = "glmnet",
               trControl = train_control,
               metric = "ROC")

counter<-0

auc_list_test_2<-list()

for(i in 1:length(sampled_rows_list)){
  print(counter)
  
  test_gene_expression_matrix<-rankGenes(counts_test[sampled_rows_list[[i]],])
  test_scores <- simpleScore(t(test_gene_expression_matrix), upSet = gene_signature_up)
  test_scores_down <- simpleScore(t(test_gene_expression_matrix), upSet = gene_signature_down)
  
  score_label<-as.numeric(predict(model, data.frame(score=scale(test_scores[1]),down_score=scale(test_scores_down[1])), type = "prob")[,"MET"])
  print(pROC::auc(meta_test$outer[test_vals_list[[i]]], score_label))
  auc_list_test_2[[i]]<-pROC::auc(meta_test$outer[test_vals_list[[i]]], score_label)
  
  counter<-counter+1
}

auc_signature_two_sided<-data.frame(model_type="SSGSEA Signature (two-sided)", AUC=unlist(auc_list_test_2))

train<-data.frame(read.csv("kulkarni_flux.csv", header = T, row.names =1))
train<-train[,-dim(train)[2]]
meta_train<-as.data.frame(train$outer)
colnames(meta_train)[1]<-"outer"
meta_train$outer<-ifelse(meta_train$outer == 1, "MET", "CTL")
flux_train<-train[,-1]

assay_flux_tot<-assay_metabolites
meta_test<-as.data.frame(test_vals)
colnames(meta_test)[1]<-"outer"
meta_test$outer<-ifelse(meta_test$outer == 1, "MET", "CTL")
flux_test<-assay_flux_tot[,-1]

Train <- flux_train
Train_meta <- as.data.frame(meta_train)
colnames(Train_meta)[1]<-"outer"

design <- model.matrix(~ 0 + factor(Train_meta$outer))
colnames(design) <- c("outer_1", "outer_0")
fit <- lmFit(t(Train), design)
contrast <- makeContrasts(outer_1 - outer_0, levels=design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

top_flux <- topTable(fit2, number=Inf, adjust.method="BH")

significant_flux <- top_flux %>%
  filter(P.Value < 0.05)

Train_flux<-flux_train[,rownames(significant_flux)]
# 
boruta_res <- Boruta(factor(Train_meta$outer) ~ ., data = Train_flux, pValue=0.01, doTrace = 0)
selected_features <- getSelectedAttributes(boruta_res, withTentative = TRUE)

Train_flux<-Train_flux[,selected_features]

model <- train(x = as.matrix(Train_flux), y = factor(Train_meta[,1]),
               method = "rbfDDA",
               trControl = train_control,
               metric = "ROC")

counter<-0

auc_list_test_knn_flux<-list()

for(i in 1:length(sampled_rows_list)){
  print(counter)
  
  score_label<-as.numeric(predict(model, data.frame(flux_test[sampled_rows_list[[i]],selected_features]), type = "prob")[,"MET"])
  print(pROC::auc(meta_test$outer[test_vals_list[[i]]], score_label))
  auc_list_test_knn_flux[[i]]<-pROC::auc(meta_test$outer[test_vals_list[[i]]], score_label)
  
  counter<-counter+1
}

auc_knn_flux<-data.frame(model_type="iMAT Flux (rbfDDA)", AUC=unlist(auc_list_test_knn_flux))

Train_gene <- counts_train

significant_genes <- top_genes %>%
  filter(P.Value < 0.05)

Train_gene<-Train_gene[,rownames(significant_genes)]
# 
boruta_res <- Boruta(factor(Train_meta$outer) ~ ., data = Train_gene, pValue=0.01, doTrace = 0)
selected_features <- getSelectedAttributes(boruta_res, withTentative = TRUE)

Train_gene<-Train_gene[,selected_features]

model <- train(x = as.matrix(Train_gene), y = factor(Train_meta[,1]),
               method = "bayesglm",
               trControl = train_control,
               metric = "ROC")

counter<-0

auc_list_test_bayes_flux<-list()

for(i in 1:length(sampled_rows_list)){
  print(counter)
  
  score_label<-as.numeric(predict(model, data.frame(counts_test[sampled_rows_list[[i]],]), type = "prob")[,"MET"])
  print(pROC::auc(meta_test$outer[test_vals_list[[i]]], score_label))
  auc_list_test_bayes_gene[[i]]<-pROC::auc(meta_test$outer[test_vals_list[[i]]], score_label)
  
  counter<-counter+1
}

auc_bayes_gene<-data.frame(model_type="Gene Expression (bayesglm)", AUC=unlist(auc_list_test_bayes_gene))

Train_kegg<- read.csv("kulkarni_kegg_1.csv", header = T, row.names = 1)
Train_kegg<-data.frame(outer=Train_kegg[,1],predict(preProcess(Train_kegg,method = c("center", "scale")), Train_kegg)[,-1])
Train_kegg$outer<-ifelse(Train_kegg$outer==1,"MET","CTL")

Test_kegg<-read.csv("MET_CTL_strat_kegg_1.csv", header = T, row.names = 1)
Test_kegg<-predict(preProcess(Test_kegg,method = c("center", "scale")), Test_kegg)

model <- train(x = as.matrix(Train_kegg[,-1]), y = factor(Train_kegg[,1]),
               method = "glmnet",
               trControl = train_control,
               metric = "ROC")

auc_list_test_glmnet_kegg<-list()

for(i in 1:length(sampled_rows_list)){
  print(counter)
  
  score_label<-as.numeric(predict(model, data.frame(Test_kegg[sampled_rows_list[[i]],]), type = "prob")[,"MET"])
  print(pROC::auc(meta_test$outer[test_vals_list[[i]]], score_label))
  auc_list_test_glmnet_kegg[[i]]<-pROC::auc(meta_test$outer[test_vals_list[[i]]], score_label)
  
  counter<-counter+1
}

auc_glmnet_kegg<-data.frame(model_type="KEGG Signatures (glmnet)", AUC=unlist(auc_list_test_glm_kegg))

auc_all<-rbind(auc_imat_fingerprint, auc_kegg_fingerprint, auc_signature_one_sided,auc_signature_two_sided,auc_bayes_gene,auc_knn_flux,auc_glmnet_kegg)

colnames(auc_all)[1]<-"Model Type"
auc_all$`Model Type`<-factor(auc_all$`Model Type`, levels = c("iMAT Flux Fingerprint", "KEGG Pathway Fingerprint","SSGSEA Signature (two-sided)","SSGSEA Signature (one-sided)","iMAT Flux (rbfDDA)","Gene Expression (bayesglm)", "KEGG Signatures (glmnet)"))

first_type <- unique(auc_all$`Model Type`)[1]
auc_all$highlight <- ifelse(auc_all$`Model Type` == first_type, "highlight", "normal")
auc_all$highlight <- ifelse(auc_all$`Model Type` == "SSGSEA Signature (one-sided)", "highlight_2", auc_all$highlight)

# Step 2: Plot with custom fill
ggplot(auc_all, aes(x = `Model Type`, y = AUC, fill = highlight)) +
  geom_boxplot(width = 0.3) +  
  geom_jitter(alpha = 0.1, width = 0.2) +
  scale_fill_manual(values = c("highlight" = "#FFDD57", highlight_2 = "pink","normal" = "#CCCCCC")) +
  scale_color_manual(values = c("highlight" = "#FFDD57", highlight_2 = "pink","normal" = "#CCCCCC")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 10),
    plot.caption = element_text(hjust = 0.5)
  ) +
  labs(
    title = "Comparison of Model AUCs:\nPredicting Metformin Response in Mouse Single Cell Clusters\nUsing Patient-Trained Data",
    caption = "(n = 50 bootstrap iterations)"
  ) +
  stat_compare_means(
    comparisons = list(c("iMAT Flux Fingerprint", "SSGSEA Signature (one-sided)")),
    method = "wilcox.test",
    label = "p.signif"
  )

colnames(auc_all)[1]<-"Model Type"
auc_all$`Model Type`<-factor(auc_all$`Model Type`, levels = c("iMAT Flux Fingerprint", "SSGSEA Signature (one-sided)"))

auc_all$highlight <- ifelse(auc_all$`Model Type` == "iMAT Flux Fingerprint", "highlight", "normal")
auc_all$highlight <- ifelse(auc_all$`Model Type` == "SSGSEA Signature (one-sided)", "highlight_2", auc_all$highlight)

ggplot(auc_all, aes(x = `Model Type`, y = AUC, fill = highlight)) +
  geom_boxplot(width = 0.3) +  
  geom_jitter(alpha = 0.1, width = 0.2) +
  scale_fill_manual(values = c("highlight" = "#FFDD57", highlight_2 = "pink","normal" = "#CCCCCC")) +
  scale_color_manual(values = c("highlight" = "#FFDD57", highlight_2 = "pink","normal" = "#CCCCCC")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 10),
    plot.caption = element_text(hjust = 0.5)
  ) +
  labs(
    title = "Comparison of Model AUCs:\nPredicting Metformin Response in Mouse Single Cell Clusters\nUsing Patient-Trained Data",
    caption = "(n = 50 bootstrap iterations)"
  ) +
  stat_compare_means(
    comparisons = list(c("iMAT Flux Fingerprint", "SSGSEA Signature (one-sided)")),
    method = "wilcox.test",
    label = "p.signif"
  )
