#!/bin/bash

ASSAY="$1"
MEDIUM="$2"
MODEL="$3"
N_SIMS="$4"
CORES="$5"

module load statistics/R/4.3.0

# Number of array tasks from the *_names.csv file
N=$(Rscript -e "cat(dim(read.csv('../../../../data/imat_prep_RNAseq/${ASSAY}_names.csv'))[1])")

echo "Array size N = $N"
# Submit the main array job and capture its ID
ARRAY_JOBID=$(sbatch --parsable \
  --array=1-"$N" \
  --cpus-per-task="${CORES}" \
  --job-name="${ASSAY}" \
  --output=logs/%x_%A_%a.out \
  --error=logs/%x_%A_%a.err \
  imat_slurm.sh "$ASSAY" "$MEDIUM" "$MODEL" "$N_SIMS" "$CORES")

echo "Submitted array job: ${ARRAY_JOBID}"

# OPTIONAL: follow-up job that runs ONLY if all array tasks succeed
# Uncomment and adapt if needed.
FOLLOW_JOBID=$(sbatch --parsable \
   --dependency=afterany:${ARRAY_JOBID} \
   aggregate_fluxes_launcher.sh "$ASSAY" "$MODEL")
echo "Submitted follow-up job (afterany): ${FOLLOW_JOBID}"
