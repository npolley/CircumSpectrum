library(caret)
library(pROC)
library(ggplot2)
library(SHAPforxgboost)
library(caTools)
library(matrixTests)
library(Boruta)

args<-commandArgs(trailingOnly = TRUE)
fingerprint_name<-args[1]

fingerprint_object<-readRDS(paste0("../../../data/fingerprint_prep_objects/RADAR_objects/",fingerprint_name,".rds"))
inner_var<-fingerprint_object[[6]]

#Read in Files from "temp" folder produced by fingerprint_generator_1
for(k in 1:length(inner_var)){
  temp_dir <- file.path(
    getwd(),           # or a fixed base directory if needed
    fingerprint_name,
    inner_var[[k]],
    "temp"
  )
  
  file_list <- list.files(
    path       = temp_dir,
    pattern    = "\\.rds$",
    full.names = TRUE
  )
  
  if (length(file_list) == 0) {
    warning("No .rds files found in ", temp_dir)
    next
  }
  
  perf_list   <- list()
  models_list <- list()
  vals_list   <- list()
  
  data       <- readRDS(file_list[[1]])
  assay      <- data[[1]]
  assay_type <- data[[3]]
  model_prim <- data[[5]]
  inner      <- gsub("\\.rds$", ".csv", data[[2]])

i<-1
for (file in file_list) {
  # Read the .rds file
  data <- readRDS(file)
  
  perf <- data[[4]]
  model<-data[[5]]
  val<-data[[6]]
  #inner_var<-data[[7]]
  
  perf_list <- append(perf_list, list(perf))
  models_list[[i]] <- model
  vals_list <- append(vals_list, list(val))
  
  i<-i+1
}

vals_subsystem<-do.call(cbind, vals_list)
vals_subsystem<-data.frame(outer=model_prim$trainingData[,1], vals_subsystem)

inverted<-ifelse(colMeans(subset(vals_subsystem, outer == "Class1")[,-1])<colMeans(subset(vals_subsystem, outer == "Class0")[,-1]), 1, 0)

#Filter Subsystem Models by Performance
perf_subsystem_unfilt <- do.call(rbind, perf_list)
perf_subsystem_unfilt$inverted<-inverted
perf_subsystem_unfilt$auc_diff<-perf_subsystem_unfilt$auc_diff*-1
n_features<-c()
for(i in 1:length(models_list)){
  n_features<-c(n_features, length(models_list[[i]]$coefnames))
}

perf_subsystem_unfilt$n_features<-n_features
perf_subsystem_unfilt$index<-rownames(perf_subsystem_unfilt)

vals_subsystem<-predict(preProcess(vals_subsystem,method = c("center", "scale")), vals_subsystem)

models_list<-models_list#[-c(rm)]

total_fingerprints<-list()
if(data[[3]] != "small_experiment"){
  total_fingerprints[["All"]]<-list(type=assay_type, assay=assay, inner=inner , fingerprints_primary=models_list, auc_primary=perf_subsystem_unfilt, subsystems=perf_subsystem_unfilt,
                                    self_score_primary=vals_subsystem, train_primary=read.csv(paste0(temp_dir,"/",inner,".csv"),header=T,row.names=1)
  )
}else{
  total_fingerprints[["All"]]<-list(type=assay_type, assay=assay, inner=inner , fingerprints_primary=models_list, auc_primary=perf_subsystem_unfilt, subsystems=perf_subsystem_unfilt, 
                                    self_score_primary=vals_subsystem,  train_primary=read.csv(paste0(temp_dir,"/",inner,".csv"),header=T,row.names=1))
}

saveRDS(total_fingerprints, paste0(assay,"_.rds"))
}
