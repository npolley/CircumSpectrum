# install_cplex_rcplex.R
# Helper script to install the Rcplex R interface to IBM CPLEX.
# NOTE: This assumes CPLEX Optimization Studio is already installed
# and that the required environment variables / paths are set by your cluster.

message("Checking for CPLEX environment...")

cplex_env <- Sys.getenv(c(
  "CPLEX_STUDIO_DIR",
  "CPLEX_STUDIO_LIB",
  "CPLEX_PATH"
))

print(cplex_env)

if (!any(nzchar(cplex_env))) {
  warning(
    paste(
      "CPLEX environment variables (e.g. CPLEX_STUDIO_DIR, CPLEX_STUDIO_LIB, CPLEX_PATH)",
      "are not set. Please configure them in your module or shell before installing Rcplex."
    )
  )
} else {
  message("CPLEX environment variables detected (at least one non‑empty).")
}

# Install dependencies required by Rcplex
dep_pkgs <- c("slam", "Matrix")
to_install_dep <- setdiff(dep_pkgs, rownames(installed.packages()))
if (length(to_install_dep) > 0) {
  install.packages(to_install_dep)
}

# Try to install Rcplex from CRAN (or the CRAN archive / local tarball on the cluster)
# On some systems Rcplex is only available as a source tarball.
message("Attempting to install Rcplex from source...")

tryCatch(
  {
    install.packages("Rcplex", type = "source")
    message("Rcplex installation attempt finished.")
  },
  error = function(e) {
    message("Rcplex installation failed:")
    message(conditionMessage(e))
    message(
      paste(
        "If this is a cluster, you may need to:",
        "\n  1) Load the CPLEX module (e.g. `module load cplex`),",
        "\n  2) Ensure CPLEX include and lib directories are in your compiler/linker paths,",
        "\n  3) Use a site‑provided Rcplex tarball and install with:",
        "\n     install.packages('Rcplex_0.3-6.tar.gz', repos = NULL, type = 'source')"
      )
    )
  }
)
