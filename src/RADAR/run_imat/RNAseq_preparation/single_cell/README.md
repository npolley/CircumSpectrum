# single-cell RNA-seq preparation

This folder contains scripts to aggregate single-cell RNA-seq data into microclusters and generate iMAT-ready discretized expression matrices.

## Files

- `expr_2_imat_sc.R`  
  R script that reads a Seurat object from `data/sc_objects/`, performs microclustering using Seurat (neighbor graph and clustering at a user-specified resolution), computes mean expression per microcluster, maps gene identifiers to the appropriate metabolic model (`humanGEM` or `mouseGEM`), and discretizes expression into iMAT-compatible integer states. It writes out `<ASSAY>_int.csv` and `<ASSAY>_names.csv` for iMAT calibration, as well as `<ASSAY>_microclusters.csv` linking individual cells to their assigned microclusters.

- `expr_2_imat_sc_tui.sh`  
  Terminal UI wrapper that scans `data/sc_objects/` for available `.rds` Seurat objects, lets you select an assay, choose a model, and set the microclustering resolution (default 40) via `whiptail`, then launches `expr_2_imat_sc.R` with the chosen parameters.

- `expr_2_imat_launcher_sc.sh`  
  Lightweight launcher script used by the TUI wrapper to call the single-cell iMAT preparation pipeline.

## Usage

1. Place a preprocessed Seurat object (`<ASSAY>.rds`) in `data/sc_objects/`.  
2. Run `expr_2_imat_sc_tui.sh` and select the assay, model, and desired microclustering resolution.  
3. Use the generated `<ASSAY>_int.csv`, `<ASSAY>_names.csv`, and `<ASSAY>_microclusters.csv` as inputs for downstream iMAT calibration.
