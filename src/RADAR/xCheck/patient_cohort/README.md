<p align="center">
  <img src="../../../../img/radar_xcheck.png" alt="RADAR | xCheck logo" width="420">
</p>

<h1 align="center">RADAR | xCheck

---

## Overview

RADAR | xCheck is an R Shiny dashboard for interactive exploration of iMAT flux output across multiple patient cohorts. The app implements a cohort registry, stratification-aware outer and inner comparisons.

The dashboard loads flux, gene expression, and clinical data for a selected cohort from a centralized registry (`xCheck_cohort_registry.csv`), then guides the user through:

1. Cohort selection and data loading.  
2. Optional stratification of patients by clinical and/or gene expression features.  
3. Definition of outer (upper vs lower) and inner (subgroup) comparisons.  
4. Pre-analysis summary and parameter tuning (AUC, log2FC, FDR).  
5. Full RADAR analysis with interactive plots and exportable fingerprint objects.  

---

## Cohort registry and data inputs

- Cohorts are defined in `xCheck_cohort_registry.csv` with `folder_name`, `cohort_name`, and `disease`.  
- For a given `folder_name = PREFIX`, the app expects the following files under `../../../../data/RADAR_xCheck_cohort/`:
  - `PREFIX_flux.csv` – reaction-level iMAT flux matrix.  
  - `PREFIX_norm_filtered.csv` – gene expression matrix used for stratifications.  
  - `PREFIX_clinical.csv` – clinical annotation table.  
  - `PREFIX_meta_clinical.csv` – metadata for clinical variables, including variable `type` (`continuous`, `binary`, `categorical`).  

All paths are constructed from `folder_name`, so adding a new cohort only requires registering it in the CSV and providing the corresponding files.

---

## Stratification workflow

After selecting a cohort, the app builds stratification options dynamically:

- Global stratification mode:
  - Three-class continuous: high / baseline / low via 25th and 75th percentiles.  
  - Two-class continuous: high / low via median split.  
- Clinical stratification:
  - Continuous variables are binned using the chosen scheme.  
  - Binary and categorical variables use their existing levels.  
- Gene stratification:
  - For each gene, 2- or 3-class labels (e.g. `GENE - high`, `GENE - baseline`, `GENE - low`) are precomputed and exposed via virtual select controls.  

Leaving both clinical and gene stratification empty uses the full cohort without pre-filtering.

---

## Outer and inner comparisons

### Outer comparison (required)

- Outer comparison type:
  - Three-class continuous (derived from 25th/75th percentiles).  
  - Two-class continuous (median split).  
  - Two-class binary.  
  - Categorical.  
- Variable source:
  - Clinical data or gene expression, depending on outer type.  
- Users specify:
  - “Upper” variable-level combinations (e.g. `RISK_SCORE - high`, `FLT3.ITD - high`).  
  - “Lower” variable-level combinations (e.g. `RISK_SCORE - low`).  

The app discretizes continuous variables when needed, maps selected levels to `upper` and `lower`, then derives a consistent outer label per patient; samples that cannot be assigned cleanly are removed from the analysis.

### Inner comparison (optional)

- Inner comparison types:
  - None (no inner comparison; entire stratified cohort treated as one group).  
  - Three-class continuous (high / baseline / low).  
  - Two-class continuous (high / low).  
  - Two-class binary.  
  - Categorical.  
- Variable source:
  - Clinical or gene expression, with options constrained by the selected inner type.  

If no inner comparison is chosen, all patients share an `inner = "All"` label and analysis proceeds as a two-group comparison only on the outer dimension.

---

## Pre-analysis summary

Once stratification, outer, and (optional) inner settings are submitted, the app:

- Constructs the final meta table (`metaFinal_out`) containing `outer` and `inner` labels for all retained samples.  
- Reports any errors (e.g. fewer than two samples, empty outer/inner groups) and halts analysis when necessary.  
- Displays a textual pre-analysis summary including:
  - Chosen clinical and gene stratifications.  
  - Outer upper and lower selections.  
  - Inner comparison variable.  
  - Chosen thresholds for minimum AUC, minimum log2FC, and maximum FDR.  
  - Any warnings or errors encountered.  

Two buttons are provided to start the main analysis:

- “Begin Analysis (set parameters)” – uses the user-specified AUC, log2FC, and FDR cutoffs.  
- “Begin Analysis (p < 0.05 only)” – runs a simplified analysis with AUC = 0, log2FC = 0, FDR = 0.05.  

---

## Analysis and statistics

For each inner subgroup and flux variable, the analysis pipeline:

- Combines the meta data with flux for the stratified cohort.  
- Computes ROC AUC for the outer comparison using `caTools::colAUC`.  
- Computes log2 fold change between `upper` and `lower` mean flux.  
- Applies AUC and log2FC thresholds to identify candidate reactions.  
- Performs Welch’s t-tests per reaction with `matrixTests::col_t_welch`.  
- Computes FDR-adjusted p-values, retaining only reactions below the specified FDR cutoff.  
- Encodes the signed enrichment metric as `sign_t_log_pval` = signed \(-\log_{10}(p)\) and categorizes reactions as “upper”, “lower”, or “ns” using a threshold of \(|\log_{10}(p)| > 1.301\).  

The result is merged with reaction metadata (`human_reaction_meta.csv`) to attach subsystem and annotation information and stored as `metabs_tab` for downstream visualization.

---

## Visual design and loading behavior

The app includes an updated UI theme:

- Custom `Poppins`-based typography and adjusted box styling (rounded corners, soft shadows).  
- Purple–blue gradient header with a flat matching logo tile for the Shiny dashboard.  
- Global loading overlay featuring:
  - A gradient “RADAR-style” square with an animated SVG hex-spiral.  
  - Sparkle effects rendered via CSS keyframes.  
  - Dynamic status text updated via Shiny custom message handlers (e.g. “Reading cohort registry…”, “Preparing stratification options…”, “Running xCheck analysis…”).  

Panels fade in with subtle animations when new outputs are generated.

---

## Interactive plots

### Plot 1 – Significant reactions across subsystems

- Displays `sign_t_log_pval` versus metabolic subsystem, using:
  - Single-inner mode: points colored by significance category (upper / lower / ns).  
  - Multi-inner mode: offsets on the x-axis per inner level, joined by thin lines per reaction to show trajectories across inner subgroups.  
- Horizontal reference lines at \(\pm 1.3\) mark the \(|\log_{10}(p)| = 1.3\) threshold.  
- Interactive features:
  - Plotly zoom/pan controls.  
  - Box/lasso selection to choose subsets of reactions.  
  - A custom handler to highlight the selected subsystem label.  

A PDF version of this plot is generated and downloadable as `all_subsystems.pdf`.

### Reaction table and selection

- Brushed points from Plot 1 populate a `reactable` table containing:
  - Means in upper and lower groups.  
  - log2FC, Welch’s t-statistic, p-value, FDR, signed \(-\log_{10}(p)\), and metadata columns.  
- Selecting a row in this table sets the current `flux_name` and updates the boxplot variable selector.  
- The brushed subset can be exported as CSV via the “selected reactions” download handler.  

### Plot 3 – Boxplots and correlation view

Plot 3 focuses on a single selected reaction:

- Boxplot view:
  - Compares flux distributions between outer groups, optionally faceted across inner categories.  
  - Draws significance bars and “star” annotations based on t-test p-values.  
  - Uses a consistent “radar_theme” for minimal, legible styling.  
- Correlation view (for continuous outer variables):
  - Reconstructs the underlying continuous outer variable used to define upper/lower.  
  - Plots flux versus this variable, with linear regression and Pearson \(r\) annotated:
    - Single inner: one scatter/regression panel.  
    - Multi-inner: one panel per inner level, combined via `cowplot::plot_grid`.  

Boxplots and correlation plots are combined side-by-side when applicable and can be exported as high-resolution PDFs.

---

## Outputs and exports

The app provides multiple export options:

- Fingerprint preparation object (`.rds`):
  - Contains the filtered reaction table, final meta table, analysis thresholds, cohort label, upper/lower labels, and inner comparison levels.  
- Plot exports:
  - PDF for Plot 1 (`all_subsystems.pdf`).  
  - PDF for Plot 3 with file naming based on the selected flux variable.  
- CSV of selected reactions from Plot 1 for further downstream processing.  

---

## Dependencies

The application relies on:

- `shiny`, `shinydashboard`, `shinyWidgets`  
- `caTools`, `matrixTests`, `MatrixGenerics`  
- `ggplot2`, `plotly`, `cowplot`  
- `reactable`, `data.table`  

The UI is built with:

```r
dashboardPage(
  skin   = "purple",
  header = header,
  sidebar = dashboardSidebar(disable = TRUE),
  body   = body
)
