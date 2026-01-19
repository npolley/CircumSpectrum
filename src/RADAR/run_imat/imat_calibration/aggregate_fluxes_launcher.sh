#!/bin/bash
#SBATCH -J aggregate_fluxes
#SBATCH --nodes=1
#SBATCH --mem=100G

ASSAY="$1"
MODEL="$2"

module load statistics/R/4.3.0

Rscript aggregate_fluxes.R "ASSAY" "MODEL"

