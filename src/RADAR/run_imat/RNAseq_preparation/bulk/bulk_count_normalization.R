library(DESeq2)
library(org.Hs.eg.db)

args <- commandArgs(trailingOnly = TRUE)

assay_name<-args[1]

countData <- read.csv(paste0("../../../../../data/count_RNAseq/",assay_name,"_counts.csv"), header = T, row.names = 1)
meta<-read.csv(paste0("../../../../../data/count_RNAseq/",assay_name,"_meta.csv"), header = T, row.names = 1)

countData<-countData[rownames(meta),]
countData<-t(countData)
dds <- DESeqDataSetFromMatrix(countData = countData, colData = meta, design = ~ treatment)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

normalized_counts<-t(normalized_counts)

ensembl_ids <- colnames(normalized_counts)

mapping <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = ensembl_ids,
  columns = "SYMBOL",
  keytype = "ENSEMBL"
)