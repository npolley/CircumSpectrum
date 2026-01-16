library(gembox)

args <- commandArgs(trailingOnly = TRUE)

assay  <- args[2]
medium <- args[3]
model_name  <- args[4]
n_sims <-if (length(args) >= 5) as.integer(args[5]) else 100
cores  <- if (length(args) >= 6) as.integer(args[6]) else 10L



file<-read.csv(paste0(assay,"_names.csv"), header = T, row.names = 1)
names<-file[,1]
media <- readRDS(paste0("../../../../data/starting_metabolites/",medium,".rds"))
expr_int<-read.csv(paste0(assay,"_int.csv"), row.names = 1)

i <- as.integer(args[1])

model <- readRDS(paste0("../../../../data/GEM_models/", model_name,".rds"))
model <- set.medium(model, media, set.all=TRUE)

bm <- get.biomass.idx(model, rgx = "MAR04413")
model <- set.rxn.bounds(model, bm, lbs=0, relative=TRUE)

set.seed(3742109)

patient<-as.integer(expr_int[,i])
print(paste0("Running sample ", names[i]))
samples<-imat(model, patient, samp.pars=list(nc=cores, n.sample=n_sims))
saveRDS(samples, paste0(names[i], ".rds"))

