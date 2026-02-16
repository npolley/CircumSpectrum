# run_imat

This folder contains the scripts and configuration needed to run the iMAT workflow.

## Subdirectories

- `RNAseq_preparation/`  
  Executes normalization of bulk or single-cell RNA-seq gene count tables and formats the discretized input matrices required for iMAT calibration (run locally). 

- `imat_calibration/`  
  Uses the iMAT algorithm from the GEMbox library to compute flux scores from the discretized RNA-seq data, intended to be run on the cluster.
