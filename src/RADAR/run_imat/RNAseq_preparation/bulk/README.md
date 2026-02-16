# bulk RNA-seq normalization

This folder contains scripts to normalize bulk RNA-seq gene count tables and prepare them for downstream iMAT discretization.

## Files

- `bulk_count_normalization.R`  
  R script that reads `*_counts.csv` and corresponding `*_meta.csv` files from `data/count_RNAseq/`, performs size-factor normalization with DESeq2, and writes normalized expression matrices as `*_norm.csv` into `data/norm_RNAseq/`.

- `bulk_count_normalization_tui.sh`  
  Terminal UI wrapper that scans `data/count_RNAseq/` for available `*_counts.csv` assays, lets you select an assay via `whiptail`, and launches `bulk_count_normalization.R` with the chosen assay name.

- `expr_2_imat.R`  
  R script that takes a normalized assay (`*_norm.csv`) and a metabolic model (`humanGEM` or `mouseGEM`), maps gene identifiers to the appropriate model IDs, log-transforms expression, and discretizes gene-level expression into iMAT-compatible integer states, writing `*_int.csv` and `*_names.csv` for use in `imat_prep_RNAseq/`.

- `expr_2_imat_tui.sh`  
  Terminal UI wrapper that scans `data/norm_RNAseq/` for `*_norm.csv` files, lets you select an assay and model via `whiptail`, and launches `expr_2_imat.R` with the chosen options.

## Usage

1. Ensure input count and metadata files are present in `data/count_RNAseq/` as `<ASSAY>_counts.csv` and `<ASSAY>_meta.csv`.  
2. Run `bulk_count_normalization_tui.sh` to generate the corresponding `<ASSAY>_norm.csv` in `data/norm_RNAseq/`.  
3. Run `expr_2_imat_tui.sh` to select a normalized assay and model and produce the discretized `*_int.csv` and `*_names.csv` files for iMAT calibration.
