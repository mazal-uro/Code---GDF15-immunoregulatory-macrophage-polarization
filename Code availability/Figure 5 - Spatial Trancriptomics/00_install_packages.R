# 00_install_packages.R
# One-time setup for the Takehara et al. prostate Xenium pipeline.
# Tested with R 4.3+ and Seurat 5.1+.

# ---- CRAN ----
cran_pkgs <- c(
  "Seurat",        # core single-cell + Xenium loader (>= 5.1)
  "SeuratObject",
  "Matrix",
  "data.table",
  "arrow",         # read Xenium .parquet transcript files
  "R.utils",       # gunzip the .parquet.gz files in 01_load_xenium.R
  "ggplot2",
  "patchwork",
  "dplyr",
  "tidyr",
  "RColorBrewer",
  "viridis",
  "future",        # parallelization
  "future.apply",
  "harmony",       # integration across samples
  "presto"         # fast Wilcoxon for FindMarkers
)
new_cran <- cran_pkgs[!cran_pkgs %in% installed.packages()[, "Package"]]
if (length(new_cran)) install.packages(new_cran)

# ---- Bioconductor ----
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
bioc_pkgs <- c(
  "SpatialExperiment",
  "SingleCellExperiment",
  "scran",
  "scater",
  "BiocParallel",
  "rhdf5"            # Xenium cell_feature_matrix.h5
)
new_bioc <- bioc_pkgs[!bioc_pkgs %in% installed.packages()[, "Package"]]
if (length(new_bioc)) BiocManager::install(new_bioc, update = FALSE, ask = FALSE)

# ---- Optional: Banksy for spatial-aware clustering ----
# Banksy uses both gene expression and local neighborhood composition
# to detect spatially coherent clusters — useful for niche identification.
if (!"Banksy" %in% installed.packages()[, "Package"]) {
  BiocManager::install("Banksy", update = FALSE, ask = FALSE)
}

# ---- Sanity check ----
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})
cat("Seurat version: ", as.character(packageVersion("Seurat")), "\n")
stopifnot(packageVersion("Seurat") >= "5.1.0")
cat("Setup complete.\n")
