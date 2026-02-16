 <p align="center">
   <img src="../../img/coralie.png" alt="CORALIE logo" width="400">
 </p>

 <h1 align="center">CORALIE</h1>

 <p align="center">
    CORALIE (Correlation of Reaction Activities Localized in Experiments) is an interactive Shiny application for multi-tier exploration of metabolic fingerprints and related features. It is organized around a single **Analysis Interface** that exposes three tiers of comparison: cohort-level fingerprints, pairwise subsystem correlations, and reaction-level drill-downs.
 </p>

---

# Using CORALIE

## Starting the application

1. Ensure the required R packages are installed (`shiny`, `shinydashboard`, `shinyjs`, `shinyWidgets`, `shinyjqui`, `dplyr`, `tidyr`, `ggplot2`, `ggpubr`, `plotly`, and dependencies).  
2. Prepare the input fingerprint and metadata files expected by CORALIE and confirm that file paths inside `CORALIE.R` point to the correct locations.  
3. From the CORALIE application directory, start the app via:
   - `bash CORALIE_gui.sh`, or  
   - `Rscript CORALIE.R` (or `shiny::runApp("CORALIE.R")` from R).  
4. Open the URL printed in the console if your browser does not open automatically.  
5. In the sidebar, select **Analysis Interface** (Tier 1 by default), then move through tiers as needed.

## Tier 1 – Compare fingerprints to a reference

Use Tier 1 to compare multiple fingerprints (e.g. cohorts, conditions, timepoints) against a common reference assay.

Typical workflow:

1. **Select reference assay**  
   - Choose the fingerprint that will serve as the reference (e.g. a control cohort or baseline condition).  
2. **Select comparison fingerprints**  
   - Pick two or more additional fingerprints to compare against the reference.  
3. **Configure display options**  
   - Choose whether to visualize differences at the subsystem level, aggregate metrics, or other summary views exposed in the UI.  
4. **Generate plots**  
   - CORALIE displays how subsystem-level fingerprint scores shift relative to the reference, highlighting global changes and outliers across fingerprints.  

Use Tier 1 to obtain a high-level overview of how different cohorts or conditions diverge from a chosen baseline in terms of their metabolic fingerprint. High positive correlations indicate high phenotypic similarity.
High negative correlations indicate diverging mechanistic difference, and no correlation is indicative of no relationship or inconclusive experimental results.

## Tier 2 – Correlate subsystem features between two fingerprints

Tier 2 focuses on pairwise comparison of subsystem-level fingerprint features for two selected fingerprints.

Typical workflow:

1. **Choose two fingerprints**  
   - Select an intersecting pair of fingerprints from the interactive Tier 1 correlation plot for detailed comparison.  
2. **Select display features**  
   - Optionally filter to specific subsystem models significant to phenotypic output in reference assay. Select hierarchical clustering or diagonal optimisation.  
3. **Compute and visualize correlations**  
   - CORALIE calculates correlations between subsystem-level fingerprint features of the two fingerprints and plots them.  
4. **Interpret correlations**  
   - Identify subsystems that co-vary strongly between fingerprints (suggesting shared responses) and those that diverge (indicating condition-specific metabolic changes).  

Use Tier 2 when you want to understand how subsystem signatures relate between two particular fingerprints beyond the global view provided in Tier 1.

## Tier 3 – Reaction–subsystem correlation drill-down

Tier 3 provides a fine-grained view that links reaction-level fluxes to subsystem-level fingerprint scores across the reference assay.

Typical workflow:

1. **Select reference assay and subsystem(s)**  
   - Click the subsystems of interest provided in the Tier-1 bar graph (e.g. those that changed significantly).  
2. **Compute correlations**  
   - CORALIE computes correlations between each reaction’s flux and the overall subsystem fingerprint scores across samples of the reference assay.  
3. **Visualize and interpret**  
   - View reaction-level correlation plots and/or tables that show which reactions are most strongly associated with subsystem-level changes.  
   - Use this information to connect subsystem-level signatures back to concrete reactions that may drive or report the observed patterns.  

Use Tier 3 as a mechanistic drill-down, translating subsystem fingerprints into interpretable reaction-level hypotheses.

## General tips

- Progress and heavier computations are accompanied by a CORALIE-branded loading overlay with a logo and progress bar. Wait for the overlay to disappear before interacting with updated plots and controls.  
- You can typically move back and forth between tiers without restarting the app; previous selections often serve as context for deeper tiers.  
- When building figures for reports or presentations, use Tier 1 for global overviews, Tier 2 for concise pairwise comparisons, and Tier 3 for mechanistic detail at the reaction level.
