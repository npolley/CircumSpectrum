library("AnnotationDbi")
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library(Seurat)
library(parallel)

#1 Read in files from args

args <- commandArgs(trailingOnly = TRUE)

sc_object_name<-args[1]
model_name  <- args[2]
res <- as.numeric(args[3])

sc_object<-readRDS(paste0("../../../../../data/sc_objects/",sc_object_name,".rds"))
model <- readRDS(paste0("../../../../../data/gem_models/", model_name,".rds"))

#2 Microclustering

sc_object <- FindNeighbors(object = sc_object, 
                       dims = 1:40)

sc_object <- FindClusters(object = sc_object,
                      resolution = res)

sc_object$microclusters<-paste0("mc_",sc_object$seurat_clusters)

microclusters<-unique(sc_object$microclusters)

expr<-t(as.matrix(sc_object@assays$SCT@data))
expr<-data.frame(microclusters=sc_object@meta.data$microclusters,expr)

expr_means <- aggregate(. ~ microclusters, data = expr, FUN = mean, na.rm = TRUE)
rownames(expr_means)<-expr_means[,1]
expr_means<-expr_means[,-1]

expr<-as.data.frame(expr_means)

#2 Determine default types of gene names used in input assay

guess_keytype <- function(ids, db) {
  kts <- keytypes(db)
  props <- vapply(kts, function(kt) {
    valid <- keys(db, keytype = kt)
    mean(ids %in% valid)
  }, numeric(1))
  
  if (all(props == 0)) {
    stop("Invalid gene names. Please transpose expression table.")
  }
  
  data.frame(
    keytype    = kts,
    match_prop = props,
    stringsAsFactors = FALSE
  )[order(-props), ]
}

#3 Map human genes to ENSEMBL IDs and mouse genes to SYMBOL IDs (per model specification)

map_to_ensembl_human <- function(ids, db, multiVals = "first",
                                 min_prop = 0) {
  
  kt_scores <- guess_keytype(ids, db)
  
  # keep only keytypes with some matches
  kt_scores <- subset(kt_scores, match_prop > min_prop)
  
  message("Using keytypes: ",
          paste0(kt_scores$keytype, " (",
                 round(kt_scores$match_prop * 100, 1), "%)", collapse = ", "))
  
  # initialise result table
  res <- data.frame(input_id = ids, stringsAsFactors = FALSE)
  
  # one ENSEMBL column per nonzero keytype
  for (kt in kt_scores$keytype) {
    colname <- paste0("ENSEMBL_", kt)
    if (kt == "ENSEMBL") {
      # IDs already ENSEMBL
      res[[colname]] <- ids
    } else {
      res[[colname]] <- mapIds(db,
                               keys     = ids,
                               column   = "ENSEMBL",
                               keytype  = kt,
                               multiVals = multiVals)
    }
  }
  
  # consensus ENSEMBL: majority vote across non‑NA values
  ens_cols <- grep("^ENSEMBL_", names(res), value = TRUE)
  
  res$ensid_auto <- apply(res[ens_cols], 1, function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) return(NA_character_)
    tab <- sort(table(x), decreasing = TRUE)
    names(tab)[1]         # most frequent ENSEMBL ID
  })
  res
}

map_to_ensembl_mouse <- function(ids, db, multiVals = "first",
                                 min_prop = 0) {
  
  kt_scores <- guess_keytype(ids, db)
  
  # keep only keytypes with some matches
  kt_scores <- subset(kt_scores, match_prop > min_prop)
  
  message("Using keytypes: ",
          paste0(kt_scores$keytype, " (",
                 round(kt_scores$match_prop * 100, 1), "%)", collapse = ", "))
  
  # initialise result table
  res <- data.frame(input_id = ids, stringsAsFactors = FALSE)
  
  # one ENSEMBL column per nonzero keytype
  for (kt in kt_scores$keytype) {
    colname <- paste0("SYMBOL_", kt)
    if (kt == "SYMBOL") {
      # IDs already ENSEMBL
      res[[colname]] <- ids
    } else {
      res[[colname]] <- mapIds(db,
                               keys     = ids,
                               column   = "SYMBOL",
                               keytype  = kt,
                               multiVals = multiVals)
    }
  }
  
  # consensus ENSEMBL: majority vote across non‑NA values
  sym_cols <- grep("^SYMBOL_", names(res), value = TRUE)
  
  res$ensid_auto <- apply(res[sym_cols], 1, function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) return(NA_character_)
    tab <- sort(table(x), decreasing = TRUE)
    names(tab)[1]         # most frequent ENSEMBL ID
  })
  res
}

convert <- data.frame(symbol = colnames(expr), 
                      stringsAsFactors = FALSE)

if(model_name == "humanGEM"){
  tmp <- map_to_ensembl_human(convert$symbol, org.Hs.eg.db)
  convert <- cbind(convert, tmp[ , setdiff(names(tmp), "input_id")])
}else if(model_name == "mouseGEM"){
  tmp <- map_to_ensembl_mouse(convert$symbol, org.Mm.eg.db)
  convert <- cbind(convert, tmp[ , setdiff(names(tmp), "input_id")])
}else{
  stop("Unsupported Model")
}

message(sprintf(paste0("%.2f%% of assay genes in ", model_name, " model."),
                sum(model$genes %in% convert$ensid_auto) / length(model$genes) * 100))


colnames(expr)<-convert$ensid_auto

#4 Log normalization (comment out if previously SCT assay)

#expr_log<-log(expr+1)

#5 Descretize Expression Values

expr_log<-t(expr)

exp_int<-as.data.frame(matrix(nrow=length(model$genes), ncol = length(colnames(expr_log))))
rownames(exp_int)<-model$genes
colnames(exp_int)<-colnames(expr_log)

common_rows <- intersect(row.names(expr_log), row.names(exp_int))
exp_int[common_rows, ] <- expr_log[common_rows, ]

descretize<-function(row){
  if(any(is.na(row))){
    desc_row<-rep(0, length(row))
  }else{
    row<-as.numeric(row)
    row_quants<-quantile(row[row != 0])
    desc_row<-ifelse(row == 0 | row < row_quants[2], -1, ifelse(row > row_quants[4], 1, 0))
  }
  return(desc_row)
}

#6 Write out int and names files for "imat_prep_RNAseq" folder

out_int<-t(apply(exp_int, 1, descretize))
colnames(out_int)<-colnames(exp_int)

write.csv(out_int, paste0(sc_object_name,"_int.csv"))
write.csv(colnames(out_int), paste0(sc_object_name,"_names.csv"))
write.csv(data.frame(cells=rownames(sc_object@meta.data),microclusters=sc_object$microclusters),paste0(sc_object_name, "_microclusters.csv"))

