#!/bin/bash
#SBATCH -J aggregate_fluxes
#SBATCH --nodes=1
#SBATCH --mem=100G

module load statistics/R/4.3.0

Rscript aggregate_fluxes.R

