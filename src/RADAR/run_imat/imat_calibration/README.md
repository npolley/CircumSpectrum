# iMAT calibration

This folder contains scripts and job launchers to run the iMAT algorithm on discretized RNA-seq data using GEMbox on an HPC cluster, and to aggregate the resulting flux samples into reaction-level flux tables.

## Files

- `imat.R`  
  R script that runs the iMAT algorithm from GEMbox for a single sample (array task). It reads discretized expression (`<ASSAY>_int.csv`) and sample names (`<ASSAY>_names.csv`) from `data/imat_prep_RNAseq/`, loads the appropriate metabolic model and starting medium, configures biomass and reaction bounds, and performs flux sampling with user-specified `n_sims` and number of cores. Results are saved as `<SAMPLE>.rds`.

- `imat_slurm.sh`  
  SLURM batch script that executes `imat.R` for one array task. It sets resource requests, derives task-specific seeds, implements simple retry logic with up to three requeues on failure, and logs output and errors per array index.

- `imat_launcher.sh`  
  Launcher script that determines the number of samples from `<ASSAY>_names.csv`, submits an array job via `imat_slurm.sh` with user-specified assay, medium, model, number of simulations, and cores, and then submits a dependent follow-up job to aggregate flux results once the array completes.

- `imat_tui.sh`  
  Terminal UI wrapper that scans `data/imat_prep_RNAseq/` for available `<ASSAY>_int.csv` files, allows interactive selection of assay, medium, model, `N_SIMS`, and `CORES` via `whiptail`, and calls `imat_launcher.sh` with the chosen parameters.

- `aggregate_fluxes.R`  
  R script that reads all `.rds` files produced by iMAT for a given assay and model, extracts mean fluxes per reaction from each sample, and writes a combined `<ASSAY>_metabTable.csv` with samples as rows and reactions as columns.

- `aggregate_fluxes_launcher.sh`  
  SLURM launcher script that requests resources and calls `aggregate_fluxes.R` for a specified assay and model.

## Usage

1. Ensure discretized input files (`<ASSAY>_int.csv` and `<ASSAY>_names.csv`) and corresponding model and medium `.rds` files are available in the `data/` directory.  
2. Run `imat_tui.sh` to configure and submit the iMAT SLURM array job.  
3. After all array tasks finish, the aggregation step produces `<ASSAY>_metabTable.csv`, which can be used for downstream analyses.
