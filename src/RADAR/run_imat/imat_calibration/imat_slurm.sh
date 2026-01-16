#!/bin/bash
#SBATCH -J launchRscript
#SBATCH -t 0-6:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8        # default, can be overridden by sbatch
#SBATCH --mem=35G
#SBATCH --array=1-1              # placeholder; overridden by sbatch

module load statistics/R/4.3.0

ASSAY="$1"
MEDIUM="$2"
MODEL="$3"
N_SIMS="$4"
CORES="$5"


Rscript imat.R "$SLURM_ARRAY_TASK_ID" "$ASSAY" "$MEDIUM" "$MODEL" "$N_SIMS" "$CORES" 
