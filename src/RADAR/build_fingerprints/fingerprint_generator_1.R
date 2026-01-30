library(data.table)
library(caret)
library(xgboost)
library(SHAPforxgboost)
library(matrixTests)
library(pROC)
library(zoo)
library(FNN)
library(smotefamily)
library(Boruta)
library(nnet)
library(MLmetrics)
library(randomForest)

args<-commandArgs(trailingOnly = TRUE)
x<-as.integer(args[1])
fingerprint_name<-args[2]

reactMeta<-read.csv("../../../data/fingerprint_prep_objects/human_reaction_meta.csv", header = T, row.names = 1)
reactMetaMouse<-read.csv("../../../data/fingerprint_prep_objects/mouse_reaction_meta.csv", header = T, row.names = 1)

# Base directory that contains the fingerprint_name folder
base_dir <- getwd()   # or a fixed path if needed

fp_dir <- file.path(base_dir, fingerprint_name)

# List immediate subdirectories inside fingerprint_name (each inner)
inner_dirs <- list.dirs(path = fp_dir, full.names = TRUE, recursive = FALSE)

xg_subset <- list()

for (inner_dir in inner_dirs) {
  inner_name <- basename(inner_dir)  # e.g. "MET_cervical_11_days"
  
  rds_path <- file.path(inner_dir, paste0(inner_name, ".rds"))
  
  if (!file.exists(rds_path)) {
    warning("RDS not found for inner: ", inner_name, " at ", rds_path)
    next
  }
  
  df <- readRDS(rds_path)
  
  # Use inner_name as the list key
  xg_subset[[inner_name]] <- df
}

for (i in 1:length(xg_subset)){
  sub<-xg_subset[[i]][["data"]]
  
  Test_tot<-factor(sub[,1])
  final_eval<-list()
  train_logloss<-list()
  shaps<-list()
  pred_tot<-list()
  
  final_classifiers<-list()
  subsystem_vals_list <- list()
  
  reactMeta_filt<-reactMeta[intersect(unique(colnames(sub)),intersect(rownames(reactMeta), rownames(reactMetaMouse))),]
  subsystems<-unique(reactMeta_filt$subsystem)
  
  selected_features<-rownames(subset(reactMeta_filt, subsystem == subsystems[x]))
  
  train<-data.frame(sub[,c("outer",selected_features)])
  train$outer <- factor(ifelse(train$outer == 1, "Class1", "Class0"), levels = c("Class1", "Class0"))
  
  # class1_data <- subset(train, outer == "Class1")[,-1]
  # class1_numeric <- as.data.frame(lapply(class1_data, as.numeric))
  # class1 <- as.matrix(class1_numeric)
  # 
  # class0_data <- subset(train, outer == "Class0")[,-1]
  # class0_numeric <- as.data.frame(lapply(class0_data, as.numeric))
  # class0 <- as.matrix(class0_numeric)
  # 
  # t_test<-col_t_welch(class1, class0)
  
  #selected_features<-rownames(subset(t_test, pvalue <	0.05)) 
  #boruta_res <- Boruta(outer ~ ., data = train[,c("outer",rownames(subset(t_test, pvalue < 0.05)))], pValue=0.001, doTrace = 0)
  #selected_features <- getSelectedAttributes(boruta_res, withTentative = FALSE)
  
  # rf_model <- randomForest(outer ~ ., data = train, importance = TRUE)
  # feature_importance <- importance(rf_model, type = 2)
  # cutoff <- quantile(feature_importance[, "MeanDecreaseGini"], 0.9)
  # selected_features <- rownames(feature_importance)[feature_importance[, "MeanDecreaseGini"] > cutoff]
  
  train<-train[,c("outer",selected_features)]
  train$outer<-factor(train$outer, levels = c("Class1","Class0"))
  
  customSummary <- function(data, lev = NULL, model = NULL) {
    obs <- data$obs
    pred <- data$pred
    probs <- data[, lev[1]]
    
    # AUC
    suppressMessages({
    auc <- pROC::auc(obs, probs)
    })
    # Confusion Matrix
    cm <- confusionMatrix(pred, obs, positive = lev[1])
    
    # F1 Score
    f1 <- cm$byClass["F1"]
    
    # Balanced Accuracy
    ba <- cm$byClass["Balanced Accuracy"]
    
    c(AUC = auc, F1 = f1, BalancedAccuracy = ba)
  }
  

  control <- trainControl(method = "repeatedcv",
                          number = 5,
                          repeats = 3,
                          summaryFunction = customSummary,
                          verboseIter = FALSE,
                          classProbs = TRUE)
  
  
  pre_classifier <- caret::train(outer ~ ., data = train,
                             method = "rbfDDA",
                             trControl = control,
                             verbosity = 0,
                             metric = "AUC")
  
  avg_auc <- mean(pre_classifier$results$AUC)
  avg_f1 <- mean(pre_classifier$results$F1.F1)
  avg_balanced_accuracy <- mean(pre_classifier$results$`BalancedAccuracy.Balanced Accuracy`)
  
  selected_features<-rownames(subset(reactMeta_filt, subsystem == subsystems[x]))
  
  train<-data.frame(sub[,c("outer",selected_features)])
  
  train$outer <- factor(ifelse(train$outer == 1, "Class1", "Class0"), levels = c("Class1", "Class0"))
  
  for(j in 1:dim(sub)[1]){
    print(j)
    selected_features<-rownames(subset(reactMeta_filt, subsystem == subsystems[x]))
    
    train<-data.frame(sub[-j,c("outer",selected_features)])
    test<-sub[j, c("outer",selected_features)]
    
    class1_data <- subset(train, outer == 1)[,-1]
    class1_numeric <- as.data.frame(lapply(class1_data, as.numeric))
    class1 <- as.matrix(class1_numeric)
    
    class0_data <- subset(train, outer == 0)[,-1]
    class0_numeric <- as.data.frame(lapply(class0_data, as.numeric))
    class0 <- as.matrix(class0_numeric)
    
    t_test<-col_t_welch(class1, class0)
    #selected_features<-rownames(subset(t_test, pvalue < 0.05))
#    boruta_res <- Boruta(outer ~ ., data = train[,c("outer",rownames(subset(t_test, pvalue < 0.05)))], pValue=0.001, doTrace = 0)
#    selected_features <- getSelectedAttributes(boruta_res, withTentative = FALSE)
    
    # rf_model <- randomForest(outer ~ ., data = train, importance = TRUE)
    # feature_importance <- importance(rf_model, type = 2)
    # cutoff <- quantile(feature_importance[, "MeanDecreaseGini"], 0.9)
    # selected_features <- rownames(feature_importance)[feature_importance[, "MeanDecreaseGini"] > cutoff]
    
    train<-train[,c("outer",selected_features)]
    test<-test[,c("outer",selected_features)]
    
    train$outer <- factor(ifelse(train$outer == 1, "Class1", "Class0"), levels = c("Class1", "Class0"))
    test$outer <- factor(ifelse(test$outer == 1, "Class1", "Class0"), levels = c("Class1", "Class0"))
    
    control <- trainControl(method = "cv",
                            number = 5,
                            summaryFunction = defaultSummary,
                            verboseIter = FALSE,
                            classProbs = TRUE)
    
    grid <- pre_classifier$bestTune

   classifier <- caret::train(outer ~ ., data = train,
                             method = "rbfDDA",
                             trControl = control,
                             tuneGrid=grid)

   
    pred <- predict(classifier, newdata = test, type = "prob")
    pred_tot[[j]]<-list(pred = pred[[1]], actual = test$outer)
  }
  
  preds <- sapply(pred_tot, function(x) x$pred)
  preds <- (preds - min(preds)) / (max(preds) - min(preds))
  test_actuals <- sapply(pred_tot, function(x) ifelse(x$actual=="Class0",0,1))
  
  # Calculate AUC
  auc_value <- auc(test_actuals, preds)
  print(subsystems[x])
  print(auc_value)
  
  subsystem_vals_list[[subsystems[x]]] <- pred_tot
  
  calculate_logloss <- function(actuals, predicted_probs) {
    predicted_probs<-ifelse(predicted_probs == 1, 0.999,predicted_probs)
    predicted_probs<-ifelse(predicted_probs == 0, 0.001,predicted_probs)
    logloss <- -mean(actuals * log(predicted_probs) + (1 - actuals) * log(1-predicted_probs))
    return(logloss)
  }
  
  test_logloss<-calculate_logloss(test_actuals, preds)
  
  train_pred <- predict(classifier$finalModel, newdata = as.matrix(train[,-1]))[,1]
  train_pred <- (train_pred - min(train_pred)) / (max(train_pred) - min(train_pred))
  actuals <- ifelse(train$outer == "Class0", 0, 1)
  
  train_logloss<-calculate_logloss(actuals, train_pred)
  train_auc<-auc(actuals, train_pred)
  
  pred_wrapper <- function(model, newdata) {
    predict(model, newdata, type = "prob")[,2]
  }
  
  print(paste0("CV AUC: ", avg_auc))
  print(paste0("CV F1: ", avg_f1))
  print(paste0("CV Balanced Accuracy: ", avg_f1))
  
  auc_stats<-data.frame(subsystem=subsystems[[x]], loocv_test_auc=auc_value, train_auc=train_auc, auc_diff=auc_value-train_auc, loocv_test_logloss=test_logloss, train_logloss=train_logloss, logloss_diff=test_logloss-train_logloss ,cv_auc=avg_auc, cv_f1=avg_f1, cv_balanced_accuracy=avg_balanced_accuracy)
                          
  sys<-subsystems[[x]]
  subsystem_vals <- data.frame(setNames(data.frame(preds), sys))
  
  out_1<-list(assay=fingerprint_name, inner= names(xg_subset)[[i]], type=xg_subset[[i]][["type"]], auc_stats=auc_stats, model=pre_classifier, subsystem_vals=subsystem_vals)
  saveRDS(out_1, paste0(fingerprint_name,"/",gsub(".rds","",names(xg_subset)[i]),"/temp/",fingerprint_name,"_",names(xg_subset)[i],"_",colnames(subsystem_vals)[1],"_classifier.rds"))
}


  
    
    




