#!/bin/bash
#SBATCH -J fingerprint_pipeline
#SBATCH -t 0-10:00:00
#SBATCH -o fingerprint_pipeline.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=15G
#SBATCH -o logs/fingerprint_pipeline_%j.out
#SBATCH -e logs/fingerprint_pipeline_%j.err

OBJ="$1"
FLUX="$2"

module load statistics/R/4.3.0

# 1) preprocessing
Rscript fingerprint_preprocessing.R "$OBJ" "$FLUX"

# 2) submit array job for generator 1, capture job ID
GEN1_JOBID=$(sbatch --export=ALL,OBJ="$OBJ" fingerprint_generator_1.job | awk '{print $4}')

# 3) submit generator 2, after array finishes
sbatch --dependency=afterany:${GEN1_JOBID} --export=ALL,OBJ="$OBJ" fingerprint_generator_2.job


