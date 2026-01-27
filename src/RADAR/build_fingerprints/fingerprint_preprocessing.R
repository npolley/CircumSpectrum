library(caret)
library(FNN)
library(smotefamily)
library(data.table)
library(dplyr)

inner_var<-list("All")
outer_var<-"day_637"

fingerprint_object<-readRDS(paste0(outer_var,".rds"))
metabs<-fingerprint_object[[1]]
flux<-as.data.frame(as.matrix(fread("beataml2_flux.csv"), rownames = 1))

#For RADAR-SCepter single-cell implementation (WIP)
if(fingerprint_object[[4]] == "single_cell"){
  meta<-data.frame(outer=ifelse(fingerprint_object[[2]][[1]]=="assay_1","upper", "lower"), inner = "All")
  rownames(meta)<-rownames(fingerprint_object[[2]])
  microclusters<-read.csv("microclusters.csv", header = T, row.names = 1)
  microclusters<-microclusters[rownames(meta),]
  merged<-data.frame(cells=rownames(meta),outer=meta$outer,microclusters=microclusters$microclusters)
  reduced_merged <- merged %>%
    group_by(microclusters) %>%
    summarize(
      outer = ifelse(sum(outer == "upper") > sum(outer == "lower"), "upper", "lower"),
      .groups = 'drop'  # Drop grouping after summarizing
    )
  rownames(reduced_merged)<-reduced_merged$microclusters
  meta<-data.frame(outer=reduced_merged$outer, inner="All")
  rownames(meta)<-rownames(reduced_merged)
}else{
  meta<-fingerprint_object[[2]]
}

reactMeta<-read.csv("human_reaction_meta.csv", header = T, row.names = 1)
reactMetaMouse<-read.csv("mouse_reaction_meta.csv", header = T, row.names = 1)

flux<-flux[rownames(meta),intersect(unique(colnames(flux)),intersect(rownames(reactMeta), rownames(reactMetaMouse)))]

xg<-cbind(meta, flux)
xg_subset<-split(xg, xg[["inner"]])

for (i in 1:length(xg_subset)){
  new_dir<-paste0(getwd(),"/",outer_var,"/",gsub(".rds","",inner_var[[i]]))
  
  if (!dir.exists(new_dir)) {
    dir.create(new_dir, recursive = TRUE)
  }
  
  if(dim(xg_subset[[i]])[1] < 10){
    sub<-data.frame(xg_subset[[i]][,-2])
    class_0 <- sub[sub$outer == "lower", ]
    class_0 <- class_0[sample(nrow(class_0), 30, replace = TRUE), ]
    for (col in names(class_0)) {
      if (is.numeric(class_0[[col]])) {
        # Add small random variations (e.g., within 5% of the original value)
        class_0[[col]] <- class_0[[col]] * runif(30, 0.95, 1.05)
      }
    }
    class_1 <- sub[sub$outer == "upper", ]
    class_1 <- class_1[sample(nrow(class_1), 30, replace = TRUE), ]
    for (col in names(class_1)) {
      if (is.numeric(class_1[[col]])) {
        # Add small random variations (e.g., within 5% of the original value)
        class_1[[col]] <- class_1[[col]] * runif(30, 0.95, 1.05)
      }
    }
    sub<-data.frame(rbind(class_0,class_1))
    sub$outer<-factor(ifelse(sub$outer == "upper",1,0))
    sub[,2:dim(sub)[2]]<-predict(preProcess(sub[,2:dim(sub)[2]],method = c("center", "scale")), sub[,2:dim(sub)[2]])
    data_type<-"small_experiment"
  }else{
    sub<-xg_subset[[i]][,-2]
    sub$outer<-factor(ifelse(sub$outer == "upper",1,0))
    sub<-data.frame(upSample(x=sub, y = factor(sub$outer)))
    sub[2:dim(sub)[2]] <- data.frame(apply(sub[2:dim(sub)[2]], 2, function(x) as.numeric(as.character(x))))
    sub[,2:dim(sub)[2]]<-predict(preProcess(sub[,2:dim(sub)[2]],method = c("center", "scale")), sub[,2:dim(sub)[2]])
    data_type<-"large_cohort"
  }
  
  out_1<-list(data=sub, type=data_type)
  saveRDS(out_1, file.path(new_dir,paste0(inner_var,".rds")))
  dir.create(file.path(new_dir,"/temp"))
  write.csv(sub,file.path(new_dir,paste0("/temp/",inner_var[[i]],".csv")))
}

