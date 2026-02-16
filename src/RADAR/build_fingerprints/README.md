# build_fingerprints

This folder contains the scripts and SLURM launchers used to construct subsystem-level metabolic fingerprints from RADAR flux data. The pipeline:

1. Reads a fingerprint preparation object (RADAR_xCheck-derived metadata and design).  
2. Reads the corresponding flux table for a chosen cohort or experiment.  
3. Preprocesses and stratifies the data into inner comparisons.  
4. Trains subsystem-level classifiers in parallel (generator 1).  
5. Aggregates models and performance into final fingerprint objects (generator 2).

## Overview of the pipeline

The main entry point is:

- `generate_fingerprint_tui.sh` → interactive selection of:
  - A fingerprint prep object from `data/fingerprint_prep_objects/RADAR_objects/`.  
  - A flux assay (cohort or experimental) from `data/RADAR_xCheck_cohort/` or `data/RADAR_xCheck_experimental_assay/`.  
- `generate_fingerprint.sh` → launches the full SLURM-based fingerprint pipeline:
  1. `fingerprint_preprocessing.R`  
  2. `fingerprint_generator_1.job` → `fingerprint_generator_1.R`  
  3. `fingerprint_generator_2.job` → `fingerprint_generator_2.R` 

---

## Files

### Launchers and TUI

- `generate_fingerprint_tui.sh`  
  - Terminal UI wrapper.  
  - Lists available fingerprint prep objects in `data/fingerprint_prep_objects/RADAR_objects/` (one `.rds` per fingerprint design).  
  - Lists available flux assays (directories) from both `data/RADAR_xCheck_cohort/` and `data/RADAR_xCheck_experimental_assay/`.  
  - After selection, runs:
    - `bash generate_fingerprint.sh "<OBJ>" "<FLUX>"`.

- `generate_fingerprint.sh`  
  - SLURM job script that orchestrates the full pipeline for a given fingerprint object (`OBJ`) and flux assay (`FLUX`).  
  - Steps:
    1. Calls `Rscript fingerprint_preprocessing.R "$OBJ" "$FLUX"` to prepare data.  
    2. Submits `fingerprint_generator_1.job` as a SLURM array job (per-subsystem or per-inner unit), passing `OBJ` as an environment variable.  
    3. Submits `fingerprint_generator_2.job` with `--dependency=afterany:${GEN1_JOBID}` to aggregate all generator 1 outputs into a final fingerprint object.

- `fingerprint_generator_1.job`  
  - SLURM array job configuration for `fingerprint_generator_1.R`.  
  - Specifies walltime, resources, and `--array=1-142` (number of subsystems to process).  
  - Loads R and runs:
    - `Rscript fingerprint_generator_1.R "$SLURM_ARRAY_TASK_ID" "$OBJ"`.

- `fingerprint_generator_2.job`  
  - SLURM job for the final aggregation stage.  
  - Loads R and runs:
    - `Rscript fingerprint_generator_2.R "$OBJ"`.

- `mlp_fingerprint_generator_1.sh`  
  - Alternative SLURM launcher for `fingerprint_generator_1.R` (without explicit `OBJ`, array index only), kept for manual or experimental runs.

---

### R scripts

- `fingerprint_preprocessing.R`  
  - Arguments: `fingerprint_name` (`OBJ`), `flux_table_name` (`FLUX`).  
  - Reads the fingerprint prep object from  
    - `data/fingerprint_prep_objects/RADAR_objects/<fingerprint_name>.rds`.  
  - Loads the corresponding flux table based on object type:  
    - If type is `"cohort"` → `data/RADAR_xCheck_cohort/<FLUX>/<FLUX>_flux.csv`.  
    - If type is `"experiment"` → `data/RADAR_xCheck_experimental_assay/<FLUX>/<FLUX>_flux.csv`.  
  - Aligns flux columns with reaction metadata (`human_reaction_meta.csv` and `mouse_reaction_meta.csv`).  
  - Merges flux with metadata (outer/inner labels) from the fingerprint object.  
  - Splits data into **inner** comparisons (one folder per inner label).  
  - For each inner comparison:
    - Creates a directory: `<pwd>/<fingerprint_name>/<inner_name>/`.  
    - Performs:
      - Synthetic variation and upsampling for small assays (`data_type = "small_experiment"`), or  
      - Class-balanced upsampling and scaling for larger cohorts (`data_type = "large_cohort"`).  
    - Writes a preprocessed `.rds` file and a CSV in a `temp` subfolder:
      - `<fingerprint_name>/<inner_name>/<inner_name>.rds` (list with `data` and `type`).  
      - `<fingerprint_name>/<inner_name>/temp/<inner_name>.csv`.

- `fingerprint_generator_1.R`  
  - Arguments: subsystem index (`x` = array task ID), `fingerprint_name` (`OBJ`).  
  - Reads reaction metadata for human and mouse.  
  - For the given fingerprint:
    - Enumerates per-inner directories under `<pwd>/<fingerprint_name>/`.  
    - Loads the preprocessed data for each inner comparison.  
  - For each inner comparison:
    - Filters reactions belonging to subsystem `x` using reaction metadata.  
    - Builds training data (`outer` vs features for that subsystem).  
    - Performs model selection and training using `caret` with `rbfDDA`:
      - Repeated cross-validation to choose hyperparameters.  
      - LOOCV-like evaluation across samples, storing predicted probabilities.  
    - Computes performance metrics (AUC, F1, balanced accuracy).  
    - Saves per-subsystem, per-inner intermediate results (models, predictions, SHAP/valuation info) to `.rds` in the corresponding `temp` folders.  
  - In effect, this stage trains a **local classifier per subsystem per inner comparison**, capturing how each subsystem’s flux pattern separates the outer groups.

- `fingerprint_generator_2.R`  
  - Argument: `fingerprint_name` (`OBJ`).  
  - Reads the fingerprint prep object.  
  - For each inner comparison:
    - Locates the `temp` directory created by generator 1:
      - `<pwd>/<fingerprint_name>/<inner_name>/temp/`.  
    - Loads all `.rds` files produced by generator 1 and extracts:
      - Per-subsystem performance objects.  
      - Trained models.  
      - Subsystem-level valuation vectors.  
    - Combines subsystem-level valuation outputs into a matrix and inspects directionality (e.g. whether higher subsystem scores correspond to “upper” vs “lower” group).  
  - Aggregates across all inner comparisons to construct the final fingerprint structure:
    - Overall assay name and type.  
    - Vector of subsystem models, performances, and indices.  
    - Any metadata needed by RADAR-xCheck and HALA (e.g. `auc_primary`, `auc_final`, inner labels).  
  - Writes the final fingerprint object to `data/fingerprints/` and/or updates the corresponding prep object under `data/fingerprint_prep_objects/RADAR_objects/` (depending on your configuration).

---

## Typical usage

1. Ensure you have:
   - A fingerprint prep object at `data/fingerprint_prep_objects/RADAR_objects/<OBJ>.rds`.  
   - A matching flux assay directory `<FLUX>` under either:
     - `data/RADAR_xCheck_cohort/<FLUX>/`, or  
     - `data/RADAR_xCheck_experimental_assay/<FLUX>/`,  
     containing `<FLUX>_flux.csv`.  

2. From the `build_fingerprints` folder (or appropriate working directory), run the TUI:
   - `bash generate_fingerprint_tui.sh`  
   - Select the **fingerprint object** (`OBJ`) and **flux assay** (`FLUX`).  
   - Confirm to launch `generate_fingerprint.sh`.

3. Monitor the SLURM jobs:
   - `fingerprint_pipeline` (top-level job from `generate_fingerprint.sh`).  
   - `launchRscript` array (`fingerprint_generator_1.job`).  
   - `fingerprint_generator_2` (aggregation stage).  

4. After completion, the final fingerprint(s) and updated prep objects will be available under:
   - `data/fingerprints/`  
   - `data/fingerprint_prep_objects/RADAR_objects/`  

These fingerprints can then be used for interactive analysis in SOURIS and CORALIE along with structural inspection with the HA/LA Graph scripts.
