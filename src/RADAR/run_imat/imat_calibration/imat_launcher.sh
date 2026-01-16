#!/bin/bash

ASSAY="$1"
MEDIUM="$2"
MODEL="$3"
N_SIMS="$4"
CORES="$5"

module load statistics/R/4.3.0

N=$(Rscript -e "dim(read.csv('../../../../data/imat_prep_RNAseq/${ASSAY}_names.csv'))[1]")

sbatch --array=1-"$N" \
  --cpus-per-task="${CORES}" \
  --job-name="${ASSAY}" \
  --output="${ASSAY}.out" \
  imat_slurm.sh "$ASSAY" "$MEDIUM" "$MODEL" "$N_SIMS" "$CORES" 
