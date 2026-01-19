library(gembox)

args <- commandArgs(trailingOnly = TRUE)

array_id<-args[1]
assay  <- args[2]
medium <- args[3]
model_name  <- args[4]
seed<-args[5]
n_sims <-if (length(args) >= 6) as.integer(args[6]) else 100
cores  <- if (length(args) >= 7) as.integer(args[7]) else 10L

file<-read.csv(paste0("../../../../data/imat_prep_RNAseq/",assay,"_names.csv"), header = T, row.names = 1)
names<-file[,1]
media <- readRDS(paste0("../../../../data/starting_metabolites/",medium,".rds"))
expr_int<-read.csv(paste0("../../../../data/imat_prep_RNAseq/",assay,"_int.csv"), row.names = 1)

i <- as.integer(array_id)

model <- readRDS(paste0("../../../../data/gem_models/", model_name,".rds"))
model <- set.medium(model, media, set.all=TRUE)

bm <- get.biomass.idx(model, rgx = "MAR04413")
model <- set.rxn.bounds(model, bm, lbs=0, relative=TRUE)

set.seed(seed)

patient<-as.integer(expr_int[,i])
print(paste0("Running sample ", names[i]))
samples<-imat(model, patient, samp.pars=list(nc=cores, n.sample=n_sims))
saveRDS(samples, paste0(names[i], ".rds"))

