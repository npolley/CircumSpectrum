> <p align="center">
>   <img src="../../img/souris.png" alt="SOURIS logo" width="160">
> </p>
>
> <h1 align="center">SOURIS</h1>
>
> <p align="center">
>   SOURIS (Superimposition of Universal Reactions In Silico) is a Shiny-based visualization module of the CircumSpectrum framework that projects an individual metabolic fingerprint onto larger bulk and single-cell reference assays, and provides multi-tier, reaction-level exploration of metabolic subsystems.
> </p>
>
> <p align="center">
>   It offers linked bulk and single-cell dashboards, allowing you to move from cohort-level fingerprint scores down to subsystem and reaction AUCs, with Seurat FeaturePlots for spatial and cluster-aware visualization.
> </p>
>
> ---
>
> ## Features
>
> - Modern SOURIS visual theme with animated loading overlay and smooth fade‑in panels.
> - **Bulk tab**
>   - Load one fingerprint object (`.rds`) and two reference assays (`.csv`).
>   - Visualize aggregated fingerprint scores per reference assay (violin + boxplots with Wilcoxon p‑values and significance stars).
>   - Drill down to subsystem‑level and reaction‑level AUC and raw values.
> - **Single‑cell tab**
>   - Reuse the same fingerprints and reference assays.
>   - Load Seurat objects and align fingerprint scores to cells.
>   - Generate subsystem and reaction FeaturePlots split by reference group.
>
> ---
>
> ## R dependencies
>
> SOURIS uses the following R packages:
>
> - `shiny`, `shinydashboard`, `shinyjs`, `shinyWidgets`, `shinyjqui`
> - `plotly`, `ggplot2`, `ggpubr`
> - `dplyr`, `tidyr`, `uwot`
> - `Seurat`
>
> In the CircumSpectrum repository, these are installed via the shared helper:
>
> ```r
> source("src/install_local_dependencies.R")
> ```
>
> You can also install manually:
>
> ```r
> pkgs <- c(
>   "shiny", "shinydashboard", "shinyjs", "shinyWidgets", "shinyjqui",
>   "plotly", "ggplot2", "ggpubr",
>   "dplyr", "tidyr", "uwot",
>   "Seurat"
> )
> install.packages(setdiff(pkgs, rownames(installed.packages())))
> ```
>
> ---
>
> ## Data layout
>
> By default, SOURIS expects the following data directories relative to the project root (adjust in `SOURIS.R` if needed):
>
> ```r
> fingerprints_dir <- "./data/fingerprints"
> reference_dir    <- "./data/reference_assays"
> sc_dir           <- "./data/sc_objects"
> ```
>
> - `data/fingerprints/`: fingerprint objects (`.rds`) with an `All` component containing `fingerprints_primary` and `auc_primary`.
> - `data/reference_assays/`: reference matrices (`.csv`), samples in rows, iMAT reactions/features in columns (MAR...); first column is sample ID.
> - `data/sc_objects/`: single‑cell Seurat objects (`.rds`) with `rownames(meta.data)` matching reference sample IDs.
>
> ---
>
> ## Launching SOURIS
>
> You can launch SOURIS via a Bash GUI wrapper (e.g. `SOURIS_gui.sh`) on Unix‑like systems or directly from R on any platform.
>
> ### Unix / Linux / macOS
>
> From the project root:
>
> ```bash
> chmod +x SOURIS_gui.sh    # one time, if not already executable
> ./SOURIS_gui.sh
> ```
>
> A typical launcher runs:
>
> ```bash
> Rscript -e "shiny::runApp('SOURIS.R', launch.browser = TRUE)"
> ```
>
> ### Any platform (R console / RStudio)
>
> From an R session with the working directory set to the project root:
>
> ```r
> source("SOURIS.R")
> ```
>
> or:
>
> ```r
> shiny::runApp("SOURIS.R", launch.browser = TRUE)
> ```
>
> ---
>
> ## Basic workflow
>
> ### Bulk tab
>
> 1. **Load a fingerprint**
>    - Use the **Fingerprint selection** picker to choose one `.rds` fingerprint.
>    - Click **Select fingerprint** (enabled only when exactly one is selected).
>
> 2. **Load reference assays**
>    - After a fingerprint is loaded, a **Reference assay selection** panel appears.
>    - Select exactly two `.csv` reference assays and click **Select reference assays**.
>
> 3. **Explore results**
>    - SOURIS computes predictions for each subsystem and sample, aggregates scores, and computes AUC between the two references.
>    - Visualizations include:
>      - Violin/boxplots of aggregated fingerprint scores per reference assay, with Wilcoxon p‑values and significance stars.
>      - Subsystem AUC bar plot; clicking a bar opens subsystem and reaction‑level detail panels.
>
> ### Single‑cell tab
>
> 1. **Load a fingerprint and references** as above on the single‑cell tab.
> 2. **Load a Seurat object**
>    - Select a `.rds` object from `data/scobjects/` and click **Load single-cell object**.
>    - The app intersects reference sample IDs with `rownames(meta.data)` and attaches fingerprint scores to cells.
> 3. **Explore single‑cell views**
>    - Violin/boxplots of aggregated fingerprint scores by reference group.
>    - Subsystem AUC plots; clicking a subsystem opens reaction‑level detail.
>    - Seurat `FeaturePlot` panels showing fingerprint, subsystem, or reaction scores split by reference group.
>
> ---
>
> ## Troubleshooting
>
> Typical validation messages and their meaning:
>
> - “Please select exactly one fingerprint.” – choose a single `.rds` fingerprint file.
> - “Please select exactly two reference assays before proceeding.” – two reference `.csv` files are required for AUC comparisons.
> - “Fingerprint object has no All component.” – the fingerprint object structure does not match the expected layout (`All$fingerprints_primary`, `All$auc_primary`).
> - “No primary fingerprint classifiers found.” – `All$fingerprints_primary` is missing or empty.
> - “No matching samples between reaction values and single-cell object.” – sample IDs in reference assays do not match the Seurat object’s metadata rownames.
>
> Check the R console for full error messages and internal validation notes if a panel fails to appear or a plot is empty.
