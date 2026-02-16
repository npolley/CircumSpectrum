# install_cluster_dependencies.R
# Install all required CRAN packages for CircumSpectrum cluster applications ("run_imat" and "build_fingerprints" from RADAR)

required_packages <- c(
  "stringr",
  "gembox",
  "data.table",
  "caret",
  "matrixTests",
  "pROC",
  "zoo",
  "FNN",
  "smotefamily",
  "nnet",
  "MLmetrics",
  "randomForest",
  "ggplot2",
  "caTools"
)

installed <- rownames(installed.packages())
to_install <- setdiff(required_packages, installed)

if (length(to_install) > 0) {
  install.packages(to_install)
} else {
  message("All required packages are already installed.")
}
