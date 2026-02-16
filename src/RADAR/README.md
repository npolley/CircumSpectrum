# RADAR analysis pipeline overview

This README summarizes the main analysis flow in the RADAR (Rapid Assessment of Differentially-Active Reactions) project, from RNA-seq–based flux inference all the way to fingerprint anatomy with High-Activity / Low-Activity (HA/LA) graphs:

> **run_imat → RADAR_xCheck → build_fingerprints → HA/LA graphs**

---

## 1. From RNA-seq to fluxes: `run_imat/`

The `run_imat` folder converts RNA-seq measurements into iMAT-based flux distributions on the cluster.

### 1.1 RNAseq_preparation

- Normalize bulk or single-cell RNA-seq data.  
- Map gene identifiers to model-specific IDs (e.g. humanGEM or mouseGEM).  
- Discretize expression into iMAT-compatible integer states (−1, 0, +1).  
- Write discrete inputs to `data/imat_prep_RNAseq/`:
  - `<ASSAY>_int.csv`  
  - `<ASSAY>_names.csv`  

### 1.2 imat_calibration

- Use `imat_tui.sh` to:
  - Select an assay (`<ASSAY>_int.csv` / `<ASSAY>_names.csv`),  
  - Choose medium and model,  
  - Set `N_SIMS` and `CORES`.  
- Submit a SLURM array job via `imat_launcher.sh` / `imat_slurm.sh` that runs `imat.R` for each sample, producing per-sample flux samples as `<SAMPLE>.rds`.  
- After the array completes, run `aggregate_fluxes.R` (via `aggregate_fluxes_launcher.sh`) to combine samples into:
  - `<ASSAY>_metabTable.csv` (samples × reactions).  

**Output of this stage:** aggregated flux tables and per-sample `.rds` files for each assay, ready for downstream analysis.

---

## 2. Interactive exploration of flux fingerprints: `RADAR_xCheck/`

The `RADAR_xCheck` folder hosts the RADAR-xCheck Shiny application for interactive, stratified analysis of flux-derived fingerprints.

### 2.1 Application inputs

RADAR-xCheck uses:

- Reaction-level metadata and flux-derived features (e.g. from `data/RADAR_xCheck_cohort/` and `data/RADAR_xCheck_experimental_assay/`).  
- Cohort and experiment registries that reference flux/fingerprint-derived metrics for observational and experimental datasets.  

### 2.2 Cohort (observational) mode

- Select a RADAR cohort from the cohort registry.  
- Optionally define:
  - Clinical stratifications (e.g. response groups, risk categories).  
  - Gene-expression stratifications.  
- Choose default or custom cutoffs for AUC, log2FC, and FDR.  
- Run the xCheck analysis to identify subsystems and reactions whose flux fingerprints differ between groups.  

### 2.3 Experimental mode

- Select an experimental assay (e.g. metformin studies).  
- Define “outer” comparisons (e.g. treatment vs control) and inner comparisons.  
- Apply default/custom thresholds as above.  
- Run the xCheck analysis to detect intervention-associated fingerprint changes.

### 2.4 Key visual outputs

- Subsystem-level volcano-style plots of signed log p-values.  
- Reaction-level boxplots across selected groups.  
- Filterable tables of significant features.  
- A textual analysis report summarizing cohorts/assays, stratifications, inner/outer comparisons, thresholds, and run time.

**Output of this stage:** stratified views and exports of flux-derived comparisons and features, which then feed into fingerprint building.

---

## 3. From xCheck-ready flux features to fingerprints: `build_fingerprints/`

The `build_fingerprints` folder constructs formal subsystem-level **fingerprints** from the flux and feature outputs that have been prepared and explored via RADAR-xCheck.

- Take fingerprint preparation objects produced by RADAR-xCheck.  
- Fit subsystem-level models (e.g. classifiers) that summarize metabolic behavior as fingerprint scores.  
- For each assay, write:
  - Fingerprint object: `data/fingerprints/<FINGERPRINT_NAME>_.rds`.  
  - Preparation/metadata object: `data/fingerprint_prep_objects/RADAR_objects/<FINGERPRINT_NAME>.rds`.  

These objects contain:

- Subsystem-level performance metrics (e.g. LOOCV AUC per subsystem).  
- Reaction-level coefficients for each subsystem fingerprint model.  
- Labels describing which groups the fingerprint separates (e.g. high vs low activity).

**Output of this stage:** reusable fingerprint objects and prep metadata for each assay, suitable for phenotypic analysis and downstream visualization.

---

## 4. Fingerprint anatomy: HA/LA graphs (`HALA.R`)

Once fingerprints are built, their internal feature anatomy can be explored using the **HA/LA (High-Activity / Low-Activity) graph** scripts in the RADAR folder.

### 4.1 Files

- `HALA.R`  
  - Reads:
    - Fingerprint object: `data/fingerprints/<FINGERPRINT_NAME>_.rds`.  
    - Prep object: `data/fingerprint_prep_objects/RADAR_objects/<FINGERPRINT_NAME>.rds`.  
  - Derives:
    - Subsystem-level performance data (`loocv_test_auc`) and final AUC.  
    - Reaction-level coefficients contributing to each subsystem fingerprint.  
  - Produces:
    - A subsystem plot of signed –log(p) values with t-statistic sign, labeled with directional group labels (left/right low- vs high-activity).  
    - A barplot of subsystem LOOCV AUCs (if mean performance < 0.95), with dashed reference lines for AUC thresholds.  
  - Writes:
    - `<FINGERPRINT_NAME>_HA_LA.svg` – High-Activity / Low-Activity fingerprint anatomy graph.  
    - `<FINGERPRINT_NAME>_stats.csv` – subsystem-level performance summary.  
    - `<FINGERPRINT_NAME>_subsystems.csv` – reaction-level statistics for the fingerprint.  

- `HALA_tui.sh`  
  - Scans `data/fingerprint_prep_objects/RADAR_objects/` for available fingerprint `.rds` files.  
  - Presents a menu to select a fingerprint name.  
  - Calls `Rscript HALA.R "<FINGERPRINT_NAME>"` to generate HA/LA outputs.

### 4.2 Typical HA/LA usage

1. Confirm the fingerprint you want to inspect has been built and saved under:
   - `data/fingerprints/`  
   - `data/fingerprint_prep_objects/RADAR_objects/`  
2. From the RADAR folder, run:
   - `bash HALA_tui.sh`  
3. Choose the fingerprint from the interactive menu.  
4. Examine:
   - `<FINGERPRINT_NAME>_HA_LA.svg` to see which subsystems are high- or low-activity and how strong their models are.  
   - The accompanying `.csv` files for numerical and reaction-level detail.

---

## End-to-end flow

1. **run_imat/**  
   - RNA-seq → discretized expression → iMAT flux samples → aggregated flux tables.  

2. **RADAR_xCheck/**  
   - Use those flux-derived features to define cohorts/assays, stratify groups, and identify key subsystems/reactions.  

3. **build_fingerprints/**  
   - Turn xCheck-ready flux data into formal subsystem-level fingerprints for each assay.  

4. **HA/LA graphs (HALA scripts)**  
   - Dissect each fingerprint’s anatomy, highlighting high- and low-activity subsystems and their driving reactions.
