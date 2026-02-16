# install_local_dependencies.R
# Install all required CRAN packages for local CircumSpectrum applications

required_packages <- c(
  "shiny",
  "shinydashboard",
  "shinyjs",
  "shinyWidgets",
  "shinyjqui",
  "plotly",
  "ggplot2",
  "dplyr",
  "tidyr",
  "uwot",
  "ggpubr",
  "Seurat",
  "caTools",
  "matrixTests",
  "reactable",
  "cowplot",
  "data.table",
  "MatrixGenerics"
)

installed <- rownames(installed.packages())
to_install <- setdiff(required_packages, installed)

if (length(to_install) > 0) {
  install.packages(to_install)
} else {
  message("All required packages are already installed.")
}
