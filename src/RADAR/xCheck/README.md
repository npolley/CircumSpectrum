 <p align="center">
   <img src="../../../img/radar_xcheck.png" alt="xCheck logo" width="550">
 </p>

 <h1 align="center">RADAR-xCheck</h1>

<p align="center">
  RADAR-xCheck is an interactive interface for exploring metabolic flux distributions by subsystem and associated clinical or experimental metadata from the RADAR project. It also facilitates the export of selected stratification comparisons as fingerprint preparation objects.
  It is implemented as a Shiny dashboard (`xCheck.R`) with a small shell wrapper (`xCheck_gui.sh`) for convenient launch.</p>


# RADAR-xCheck


## Overview

The app supports two main analysis modes:

- **Cohort (observational study) mode** – stratified analysis of iMAT fluxes across clinical and/or gene-expression–defined subgroups within predefined RADAR cohorts.  
- **Experimental mode** – analysis of intervention datasets (e.g. metformin assays) using flexible “outer” and “inner” comparison definitions.

In both modes, the interface provides interactive plots, tables, and a textual report summarizing the stratification and statistics used.

## Input data

`xCheck.R` expects the following data under `data/`:

- `RADAR_xCheck_cohort/`
  - `human_reaction_meta.csv`: reaction-level metadata (e.g. metabolite IDs, subsystems, annotations).
  - `xCheck_cohort_registry.csv`: registry of observational cohorts, including folder names, cohort names, and disease labels.

- `RADAR_xCheck_experimental_assay/`
  - `xCheck_experiment_registry.csv`: registry of experimental datasets (e.g. metformin assays), including folder names, study names, GEO accessions, and treatments.

Flux fingerprint matrices and any additional per-cohort or per-experiment metadata are expected in the corresponding data subfolders referenced by the registries.

## Files in this folder

- `xCheck.R`  
  Main Shiny application script. It:
  - Loads reaction metadata and the cohort/experiment registries.
  - Defines the RADAR-xCheck dashboard layout (header, sidebar, tabbed body).
  - Provides controls for selecting cohorts or experiments, defining stratification schemes, and configuring statistical thresholds.
  - Executes the xCheck analysis, produces interactive plots and tables, and generates a textual report of the analysis settings and results.

- `xCheck_gui.sh`  
  Convenience launcher script that starts the RADAR-xCheck Shiny app from the command line (e.g. on a workstation or interactive node), opening the web interface in your browser.

## Core functionality

### Cohort (observational study) mode

In the “observational” tab, the app:

- Lets you select a RADAR cohort from `xCheck_cohort_registry.csv`.  
- Offers optional stratification controls:
  - **Stratification type** for continuous variables (e.g. tertiles, median split).
  - **Clinical stratification** (select one or more clinical variables).
  - **Gene expression stratification** (select one or more genes to define groups).
- Allows you to:
  - Use **default thresholds** (e.g. AUC = 0, log2FC = 0, FDR = 0.05).
  - Or specify **custom** cutoffs for AUC, log2FC, and FDR.
- Runs the xCheck analysis to identify flux features (metabolites/reactions) that differ across the defined groups, labeling them as significantly increased, decreased, or nonsignificant.

### Experimental mode

In the “experimental” tab, the app:

- Lets you select an experimental assay from `xCheck_experiment_registry.csv`.  
- Requires you to define:
  - An **outer comparison** parameter (e.g. treatment vs control).
  - The corresponding **upper** and/or **lower** outer comparison groups.
- Offers similar AUC/log2FC/FDR controls and “default vs custom” threshold options as in cohort mode.
- Runs an analogous xCheck analysis for experimental conditions, highlighting flux changes associated with the intervention.

## Interactive outputs

After an analysis is run (in either mode), the app provides:

- **Subsystem plot**  
  - X-axis: metabolic subsystems.  
  - Y-axis: signed log p-value (direction + significance).  
  - Points: metabolites/flux features, colored by:
    - Significance category (e.g. increased, decreased, nonsignificant), or
    - Inner comparison (for multi-inner designs, with x-offsets per inner comparison).  
  - Horizontal lines at fixed |log p| thresholds to mark significance cutoffs.  
  - Hover and selection support via Plotly.

- **Linked selection**  
  - Brushing/selection on the main plot passes the selected metabolites into detailed tables or secondary views for closer inspection.

- **Filterable results table**  
  - A table of significant (and optionally all) features with:
    - Inner comparison labels.
    - Significance category.
    - Statistics (AUC, log2FC, p-values, FDR, etc., as defined in the analysis).
  - Filters for:
    - Specific inner comparisons.
    - Selected significance categories (e.g. only “Significant – increased”, “Significant – decreased”).

- **Textual report**  
  - A text block summarizing the analysis:
    - Cohort or experiment name.
    - Clinical and gene-based stratifications used.
    - Outer and inner comparisons.
    - AUC/log2FC/FDR thresholds.
    - Errors or warnings (if any).
    - Timestamp of report generation.
   
- **Reaction-level boxplots**  
  - Displays boxplots of flux values for selected reactions across user-defined groups (e.g. clinical strata or experimental conditions).
    - Helps interpret and validate significant reactions from xCheck analyses at a single-reaction level.
    - Correlation plots are provided for continuous outer variables.

## User experience and theming

The dashboard includes:

- A custom RADAR-themed header, typography, and color palette.  
- Styled boxes, buttons, and inputs tuned for dense analytical interfaces.  
- A full-screen **loading overlay** with:
  - Animated logo.
  - Dynamic status text (“Reading cohort registry…”, “Running xCheck analysis…”, etc.).
  - Smooth fade-in transitions when results become available.

These elements are implemented via custom CSS and JavaScript hooks within `xCheck.R`.

## Running RADAR-xCheck

1. Ensure the required R packages are installed (e.g. `shiny`, `shinydashboard`, `shinyWidgets`, `ggplot2`, `plotly`, `reactable`, `data.table`, `MatrixGenerics`, and dependencies).  
2. Confirm that the RADAR data folders are populated under `data/`:
   - `data/RADAR_xCheck_cohort/`
   - `data/RADAR_xCheck_experimental_assay/`
3. From this directory, start the app using either:
   - `bash xCheck_gui.sh`, or  
   - `Rscript xCheck.R` (or `shiny::runApp("xCheck.R")` from an R session).
4. Open the displayed URL in a web browser if it does not open automatically.  
5. Use the “observational” and “experimental” tabs to configure stratifications, run analyses, inspect plots/tables, and export figures or tables as needed.
6. Click "Download Fingerprint Preparation Object" – Selected cohort stratifications, associated sample metadata, and their corresponding flux profiles are bundled into self‑contained fingerprint preparation objects and written to the fingerprint_prep_objects/ folder, where they can be used as standardized inputs for downstream fingerprint training and evaluation.
