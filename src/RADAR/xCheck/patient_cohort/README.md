RADAR | xCheck

RADAR | xCheck is an R Shiny dashboard for interactive exploration of iMAT flux assays across multiple patient cohorts. The app lets you:

    Select a cohort (e.g. BEATAML2.0 AML, TCGA PDAC, glioblastoma) from a registry file.

Optionally stratify patients by clinical variables or gene expression into two- or three-class groups.

​

Define an outer comparison (upper vs lower groups) and an inner comparison (subgroups) to compute AUC, log2FC, and differential flux statistics.

​

Visualize significant reactions across metabolic subsystems and drill down to reaction-level and flux-distribution views.

​

Export filtered reaction fingerprints and plots for downstream analysis and reporting.

    ​

Features

    Cohort registry–driven configuration

        Cohorts are not hard-coded; they are defined in a CSV registry (xCheck_cohort_registry.csv) with folder_name, cohort_name, and disease.

The app constructs all file paths from folder_name, enabling easy addition of new cohorts without changing code.

    ​

Dynamic stratification UI

    Clinical and gene expression stratifications are auto-generated based on the selected cohort.

​

Two modes for continuous variables:

    Three-class (high / baseline / low) via 25th and 75th percentiles.

    Two-class (high / low) via median split.

    ​

Non-continuous clinical variables are treated as categorical levels.

​

Stratification is optional; leaving both clinical and gene stratifications empty uses the full cohort.

    ​

Two-level comparison design

    Outer comparison: defines “upper” vs “lower” groups (e.g. high vs low risk), based on selected variables and dichotomization rules.

​

Inner comparison: defines subgroups within the stratified cohort (e.g. high/baseline/low gene expression or clinical bins).

    ​

Statistical workflow

    For each inner subgroup and flux variable, the app computes:

        ROC AUC for outer (upper vs lower).

        log2 fold change between upper and lower mean flux.

        Welch’s t-test and FDR-adjusted p-values.

    ​

User-defined thresholds for minimum AUC, minimum log2FC, and maximum FDR control which reactions are retained.

    ​

Visualizations

    Plot 1: “Significant Reactions Across All Subsystems” – scatter of signed −log⁡10(p)−log10(p) by metabolic subsystem, colored by significance or inner comparison.

​

Plot 2: Reaction trajectories within a selected subsystem, across inner comparison levels.

​

Plot 3: Boxplot view of selected flux variables across outer and inner groups, with invert option.

​

Interactive brushing and clicking link the plots and update detail views.

    ​

Outputs and export

    Download filtered reaction fingerprints (reaction-level stats and metadata) as an RDS object.

​

Download publication-ready PDFs of subsystem and detailed plots.

​

Tabular summary of selected reactions and their annotations via reactable
