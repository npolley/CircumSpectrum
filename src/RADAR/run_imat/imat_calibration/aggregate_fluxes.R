library(stringr)

rds_combo <- list.files(pattern = "*.rds", full.names = TRUE )

metabTable<-data.frame(matrix("", ncol = length(rds_combo), nrow = 13070))
for (i in 1:length(rds_combo)){
  sample<-readRDS(rds_combo[i])
  fluxes <- rowMeans(sample$result.model$sample$pnts)
  metabTable[,i]<-fluxes
  print(rds_combo[i])
}

metabTable<-t(metabTable)

rownames(metabTable)<-str_remove(str_remove(rds_combo, "./"), ".rds")
colnames(metabTable)<-sample$result.model$rxns

write.csv(metabTable, "metabTable.csv")

