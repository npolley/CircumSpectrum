# SOURIS

SOURIS is an interactive Shiny application to explore metabolic fingerprint classifiers on bulk and single‑cell transcriptomics data, compute subsystem‑ and reaction‑level metrics, and visualize them with linked plots.

The app provides:
- A “Bulk samples” dashboard for running fingerprint models on reference assays.
- A “Single‑cell” dashboard to project fingerprint scores onto Seurat objects and generate FeaturePlots.

---

## Features

- **Modern** dashboard UI with a custom SOURIS blue theme, animated loading overlay, and smooth fade‑in panels.
- **Bulk tab**
  - Load one fingerprint object (`.rds`).
  - Select two reference assays (`.csv`).
  - Visualize aggregated fingerprint scores by reference condition (violin + boxplots with Wilcoxon p‑values and significance stars).
  - Inspect subsystem‑level AUC and drill down to reaction‑level AUC and raw values.
- **Single‑cell tab**
  - Reuse the same fingerprints and reference assays.
  - Load Seurat objects (`.rds`) with matching sample IDs.
  - Visualize per‑cell fingerprint, subsystem, and reaction scores via `FeaturePlot`, split by reference group.
- **Statistics**
  - AUC computation via Mann–Whitney U on predicted scores.
  - Per‑subsystem and per‑reaction Wilcoxon tests with multiple visual encodings of significance.

---

## Launching the app

The preferred way to launch SOURIS is through the provided TUI shell wrapper:

```bash
#!/usr/bin/env bash

Rscript -e "shiny::runApp('SOURIS.R', launch.browser = TRUE)"

From the repository root:

bash
chmod +x run_souris.sh      # if you name the script run_souris.sh
./run_souris.sh

This starts the Shiny app defined in SOURIS.R and opens it in your default web browser.

You can also launch it directly from R:

r
shiny::runApp("SOURIS.R", launch.browser = TRUE)

Installation
1. Clone the repository

bash
git clone https://github.com/<user>/<repo>.git
cd <repo>

2. Install R dependencies

SOURIS uses the following R packages:

    shiny, shinydashboard, shinyjs, shinyWidgets, shinyjqui

    plotly, ggplot2, ggpubr

    dplyr, tidyr, uwot

    Seurat (single‑cell functionality)

Install them in R:

r
pkgs <- c(
  "shiny", "shinydashboard", "shinyjs", "shinyWidgets", "shinyjqui",
  "plotly", "ggplot2", "ggpubr",
  "dplyr", "tidyr", "uwot",
  "Seurat"
)

install.packages(setdiff(pkgs, rownames(installed.packages())))

If you use renv or another environment manager, you can instead restore from the project lockfile, if provided.
Data layout

By default, SOURIS expects the following data structure (relative to the project root):

text
data/
  fingerprints/
    <fingerprint_set_1>.rds
    <fingerprint_set_2>.rds
    ...
  referenceassays/
    <reference_assay_1>.csv
    <reference_assay_2>.csv
    ...
  scobjects/
    <single_cell_object_1>.rds
    <single_cell_object_2>.rds
    ...

These paths correspond to the variables:

r
fingerprints_dir <- "./data/fingerprints"
reference_dir    <- "./data/referenc_eassays"
sc_dir           <- "./data/sc_objects"

You can change them in SOURIS.R if you prefer a different layout.
Expected input formats
Fingerprints (data/fingerprints/*.rds)

    R objects containing:

        An All component with:

            fingerprints$primary (a list of trained classifiers).

            aucprimary (AUC metadata per subsystem).

    Classifiers must support predict(model, newdata = ..., type = "prob") with the positive‑class probability in column 1.

Reference assays (data/referenceassays/*.csv)

    CSV files with:

        First column: sample IDs (used as row names).

        Remaining columns: reaction features.

    Column names should match the internal reaction IDs used by the fingerprint models, or be mappable via the reactmetafilt metadata (e.g. to RECON3D IDs).

Single‑cell objects (data/scobjects/*.rds)

    Seurat objects.

    rownames(obj@meta.data) must contain at least some of the sample IDs from the first column of the reference .csv files.

    SOURIS attaches fields:

        SOURISfingerprintscore

        SOURISreferencegroup

        SOURISsubsystemscore and SOURISreactionscore for specific analyses.

If no IDs match between reference assays and a Seurat object, the app shows a modal (“No matching samples”) and does not keep the object in memory.
How to use
Bulk tab

    Load a fingerprint

        Use the “Fingerprint selection” picker to choose one .rds file.

        Click Select fingerprint (enabled only when exactly one is selected).

    Load reference assays

        After a fingerprint is loaded, a “Reference assay selection” panel appears.

        Select exactly two .csv files and click Select reference assays.

    Explore results

        The app computes predictions for each subsystem and each sample, aggregates scores, and computes AUC between the two references.

        Plots:

            Violin/boxplot of aggregated fingerprint scores per reference.

            Subsystem AUC bar plot (click any bar to open a detail panel with subsystem and reaction‑level plots).

Single‑cell tab

    Load a fingerprint

        Select one fingerprint in the single‑cell “Fingerprint selection” box and click Select fingerprint.

    Load reference assays

        Select exactly two .csv files and click Select reference assays to compute single‑cell‑aligned scores.

    Load a Seurat object

        In “Single‑cell object selection”, choose a .rds Seurat object and click Load single-cell object.

        The app aligns scores to cells, then inserts:

            A violin/boxplot panel for scores by reference group.

            Subsystem AUC plots.

            A FeaturePlot of fingerprint scores.

    Drill‑down

        Click a subsystem bar to generate subsystem‑specific FeaturePlots.

        Click a reaction entry in the reaction‑level AUC panel to generate a reaction‑specific FeaturePlot.

Troubleshooting

Common validation messages and what they mean:

    “Please select exactly one fingerprint.” — The picker must have one, not zero or multiple entries selected.

    “Please select exactly two reference assays before proceeding.” — Two reference .csv files are required for AUC comparisons.

    “Fingerprint object has no All component.” — The fingerprint .rds structure does not match the expected layout.

    “No primary fingerprint classifiers found.” — All$fingerprints$primary is missing or empty.

    “Reference assays must be loaded before loading the single-cell object.” — Load reference assays on the single‑cell tab before SC object.

    “No IDs from the first column of the reference assays were found in the single-cell object metadata rownames.” — Sample IDs do not match between reference assays and Seurat object.

    “No reactions for subsystem X in reference assay Y.” — Columns for that subsystem’s reactions are missing from the reference matrix.

Inspect the R console logs for full stack traces and internal validation checks if a panel fails to appear.
