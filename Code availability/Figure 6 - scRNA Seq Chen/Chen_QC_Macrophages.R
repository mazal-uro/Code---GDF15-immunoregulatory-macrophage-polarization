# =====================================================================
# Chen et al. (GSE141445) prostate cancer scRNA-seq
# Robust QC  ->  cluster whole dataset  ->  isolate macrophages  ->  GDF15
# ---------------------------------------------------------------------
# Fresh rebuild from the GSM4203181 data matrix. Key fixes vs the old script:
#   * percent.mt now uses the correct human pattern "^MT-" (was "^mt-",
#     which matched nothing, so the mito filter did nothing).
#   * Matrix is loaded with gene symbols set as row names, and gene/cell names
#     are de-duplicated (the old as.matrix() could coerce to character, and
#     duplicate symbols crashed CreateSeuratObject).
#   * QC is done TWO ways (MAD-based, data-driven + fixed thresholds) so you
#     can compare; MAD is the default.
#   * Macrophages are identified by a marker MODULE SCORE per cluster, not a
#     hard-coded cluster number (your old comments disagreed: 7 vs 9).
#
# CONFIRMED ABOUT THIS DATASET (checked in-session):
#   - No mitochondrial genes present -> percent.mt is 0, so MT-based QC is
#     unavailable; QC relies on nFeature/nCount.
#   - The matrix is ALREADY NORMALIZED (nFeature > nCount for some cells;
#     non-integer values), so STEP 5 does NOT run NormalizeData.
#
# NOTE: I could not run this (no R / data here). Run top-to-bottom and check
# the diagnostic plots + printed messages at each QC / annotation step.
# Written for Seurat v5 (also works on v4).
# =====================================================================

# -----------------------------------------
# Libraries
# -----------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
library(ggplot2)
library(cowplot)
library(SCpubr)

# =====================================================================
# STEP 1: Load the raw Chen count matrix (robustly)
# =====================================================================
file_path <- "~/Documents/scRNA/Chen et al/GSE141445_RAW/GSM4203181_data.matrix.txt.gz"

# fread is fast; genes are in the first column, cells across the header.
counts_dt <- data.table::fread(file_path, sep = "\t", header = TRUE)
counts_df <- as.data.frame(counts_dt)
rownames(counts_df) <- counts_df[[1]]   # first column = gene symbols
counts_df[[1]] <- NULL                  # drop that column
counts_mat <- as.matrix(counts_df)      # now purely numeric
storage.mode(counts_mat) <- "numeric"
counts_mat <- as(counts_mat, "dgCMatrix")  # sparse -> memory-efficient
rm(counts_dt, counts_df)                # free memory

# Gene symbols: mirror Seurat's "_" -> "-" replacement and force uniqueness.
# Without this, two symbols that differ only by _ vs - collapse to the same
# name and CreateSeuratObject fails with:
#   invalid class "LogMap" object: Duplicate rownames not allowed
n_dup <- sum(duplicated(gsub("_", "-", rownames(counts_mat))))
cat("Duplicate gene symbols after '_'->'-':", n_dup, "(made unique)\n")
rownames(counts_mat) <- make.unique(gsub("_", "-", rownames(counts_mat)))
colnames(counts_mat) <- make.unique(colnames(counts_mat))  # dedupe cell barcodes too

# Hard stop if anything is still duplicated (otherwise CreateSeuratObject throws
# the LogMap "Duplicate rownames" error).
stopifnot(!any(duplicated(rownames(counts_mat))),
          !any(duplicated(colnames(counts_mat))))

cat("Matrix loaded:", nrow(counts_mat), "genes x", ncol(counts_mat), "cells\n")

# =====================================================================
# STEP 2: Create the Seurat object
# =====================================================================
chen <- CreateSeuratObject(
  counts       = counts_mat,
  project      = "Chen_Data",
  min.cells    = 3,      # gene must be seen in >=3 cells
  min.features = 200     # cell must have >=200 genes
)
chen

# =====================================================================
# STEP 3: QC metrics
# =====================================================================
# Human mitochondrial genes are "MT-..." (UPPERCASE). Verify they exist —
# if this prints 0, the matrix uses a different naming and the mito filter
# would silently do nothing (exactly the bug in the old script).
n_mt <- sum(grepl("^MT-", rownames(chen)))
cat("Mitochondrial genes found (^MT-):", n_mt, "\n")

chen[["percent.mt"]]   <- PercentageFeatureSet(chen, pattern = "^MT-")
chen[["percent.ribo"]] <- PercentageFeatureSet(chen, pattern = "^RPS|^RPL")
chen[["percent.hb"]]   <- PercentageFeatureSet(chen, pattern = "^HB[AB]")  # hemoglobin (RBC contamination)

# --- Sanity check: are these RAW integer counts, or already normalized? ---
# For real UMI counts, every cell must have nFeature <= nCount, and values
# should be integers. If either check below fails, this matrix is already
# processed -> do NOT run NormalizeData again (see STEP 5 note).
cat("All nFeature <= nCount (expect TRUE for raw counts):",
    all(chen$nFeature_RNA <= chen$nCount_RNA), "\n")
cat("Non-integer values present (TRUE => already normalized):",
    any(GetAssayData(chen, layer = "counts")@x %% 1 != 0), "\n")

# Diagnostic plots (look here BEFORE trusting any threshold)
VlnPlot(chen, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
        ncol = 4, pt.size = 0)
FeatureScatter(chen, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(chen, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# =====================================================================
# STEP 4: QC filtering — MAD-based (default) AND fixed, for comparison
# =====================================================================
# Choose which filter to actually apply downstream:
qc_method <- "MAD"     # "MAD" (data-driven, recommended) or "fixed"

# ---- (a) MAD-based outlier detection (scran/scater style) ----
# A cell is an outlier if a metric is > nmads MADs from the median.
# nFeature & nCount: flag both tails (on log scale). percent.mt: upper tail only.
qc_nmads <- 3
mad_outlier <- function(x, nmads = 3, direction = c("both", "higher", "lower"), use_log = FALSE) {
  direction <- match.arg(direction)
  v   <- if (use_log) log1p(x) else x
  med <- median(v, na.rm = TRUE)
  dev <- mad(v, na.rm = TRUE)
  lower <- med - nmads * dev
  upper <- med + nmads * dev
  switch(direction,
         both   = v < lower | v > upper,
         higher = v > upper,
         lower  = v < lower)
}

# Guard: only apply the mito MAD if percent.mt actually varies (if most cells
# are 0, mad = 0 and a naive rule would flag every MT-containing cell).
mt_fail <- if (mad(chen$percent.mt, na.rm = TRUE) > 0) {
  mad_outlier(chen$percent.mt, qc_nmads, "higher", use_log = FALSE)
} else rep(FALSE, ncol(chen))

mad_fail <-
  mad_outlier(chen$nFeature_RNA, qc_nmads, "both", use_log = TRUE) |
  mad_outlier(chen$nCount_RNA,   qc_nmads, "both", use_log = TRUE) |
  mt_fail

# ---- (b) Fixed thresholds (tune to the diagnostic plots) ----
min_features <- 200
max_features <- 6000
max_counts   <- 40000
max_mt       <- 10     # percent
fixed_fail <-
  chen$nFeature_RNA < min_features |
  chen$nFeature_RNA > max_features |
  chen$nCount_RNA   > max_counts   |
  chen$percent.mt   > max_mt

# ---- Report how many each removes ----
cat(sprintf("MAD filter  : removes %d of %d cells (keeps %d)\n",
            sum(mad_fail),   length(mad_fail),   sum(!mad_fail)))
cat(sprintf("Fixed filter: removes %d of %d cells (keeps %d)\n",
            sum(fixed_fail), length(fixed_fail), sum(!fixed_fail)))

# ---- Apply the chosen filter ----
chen$qc_keep <- if (qc_method == "MAD") !mad_fail else !fixed_fail
chen_sub <- subset(chen, subset = qc_keep)
cat("Kept", ncol(chen_sub), "cells after", qc_method, "QC\n")

# =====================================================================
# STEP 5: (matrix is ALREADY NORMALIZED — do NOT re-normalize)
# =====================================================================
# The two checks above confirmed this is a processed matrix, not raw counts.
# So we DON'T call NormalizeData. Instead we place the provided values into the
# "data" layer directly. If they turn out to be linear-normalized (not logged),
# we log1p them; if already log-scaled, we leave them as-is.
rng <- range(GetAssayData(chen_sub, layer = "counts")@x)
cat("Value range:", round(rng[1], 3), "to", round(rng[2], 3), "\n")

dat <- GetAssayData(chen_sub, layer = "counts")   # the authors' normalized values
if (rng[2] > 30) {
  dat <- log1p(dat)
  cat("Applied log1p (values looked linear).\n")
} else {
  cat("Left as-is (values already appear log-scaled).\n")
}
chen_sub <- SetAssayData(chen_sub, layer = "data", new.data = dat)

# HVGs + scaling + PCA on the normalized data (no NormalizeData step).
chen_sub <- FindVariableFeatures(chen_sub, selection.method = "vst", nfeatures = 3000)
chen_sub <- ScaleData(chen_sub, features = rownames(chen_sub))
chen_sub <- RunPCA(chen_sub, features = VariableFeatures(chen_sub), verbose = FALSE)

# How many PCs? Look at the elbow.
ElbowPlot(chen_sub, ndims = 40)
n_dims <- 20   # adjust to where the elbow flattens

# =====================================================================
# STEP 6: Cluster the whole dataset (resolution high enough to split cell types)
# =====================================================================
# 0.5 typically resolves the major lineages (epithelial, T/NK, B, myeloid,
# endothelial, fibroblast) in a whole-tumor sample. Raise it if macrophages
# don't separate from other myeloid cells.
chen_sub <- FindNeighbors(chen_sub, dims = 1:n_dims)
chen_sub <- FindClusters(chen_sub, resolution = 0.5)
chen_sub <- RunUMAP(chen_sub, dims = 1:n_dims)
DimPlot(chen_sub, reduction = "umap", label = TRUE) + NoLegend()

# =====================================================================
# STEP 7: Identify the macrophage cluster(s) by marker MODULE SCORE
# =====================================================================
# Canonical macrophage / monocyte-macrophage markers. Only genes present in
# the object are used (AddModuleScore would otherwise complain).
mac_markers <- c("CD68", "CD14", "LYZ", "CSF1R", "AIF1", "TYROBP", "FCER1G",
                 "C1QA", "C1QB", "C1QC", "FCGR3A", "MRC1", "CD163", "ITGAM")
mac_markers <- intersect(mac_markers, rownames(chen_sub))
chen_sub <- AddModuleScore(chen_sub, features = list(mac_markers), name = "MacScore")
# AddModuleScore appends "MacScore1"
mac_score_by_cluster <- sort(
  tapply(chen_sub$MacScore1, chen_sub$seurat_clusters, mean),
  decreasing = TRUE
)
cat("Mean macrophage score per cluster:\n"); print(round(mac_score_by_cluster, 3))

# Auto-select cluster(s) whose mean score clearly stands out
# (mean + 1 SD across clusters). Inspect the confirmation plots below.
score_cutoff <- mean(mac_score_by_cluster) + sd(mac_score_by_cluster)
macrophage_clusters <- names(mac_score_by_cluster)[mac_score_by_cluster > score_cutoff]
cat("Auto-selected macrophage cluster(s):", paste(macrophage_clusters, collapse = ", "), "\n")
# If this looks wrong after the plots below, set it manually, e.g.:
# macrophage_clusters <- c("7")

# ---- Confirm visually before subsetting ----
FeaturePlot(chen_sub, features = c("CD68", "CD14", "LYZ", "CSF1R", "C1QA", "MRC1"), ncol = 3)
DotPlot(chen_sub, features = mac_markers) + RotatedAxis()
# Also show non-macrophage lineages so you can rule out contamination:
FeaturePlot(chen_sub, features = c("CD3E", "MS4A1", "KRT19", "AR", "PECAM1", "GDF15"), ncol = 3)

# Optional: full marker table for manual sanity check
# all_markers <- FindAllMarkers(chen_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# all_markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10) %>% print(n = Inf)

# =====================================================================
# STEP 8: Subset macrophages and re-analyze
# =====================================================================
Idents(chen_sub) <- "seurat_clusters"
macrophages <- subset(chen_sub, idents = macrophage_clusters)
cat("Macrophage cells:", ncol(macrophages), "\n")

# No NormalizeData: the subset already carries the normalized "data" layer.
macrophages <- FindVariableFeatures(macrophages, selection.method = "vst", nfeatures = 2000)
macrophages <- ScaleData(macrophages, features = rownames(macrophages))
macrophages <- RunPCA(macrophages, verbose = FALSE)
ElbowPlot(macrophages, ndims = 30)
mac_dims <- 15   # adjust to the macrophage elbow

macrophages <- FindNeighbors(macrophages, dims = 1:mac_dims)
macrophages <- FindClusters(macrophages, resolution = 0.2)
macrophages <- RunUMAP(macrophages, dims = 1:mac_dims)
DimPlot(macrophages, reduction = "umap", label = TRUE)

# =====================================================================
# STEP 9: GDF15 in macrophages
# =====================================================================
FeaturePlot(macrophages,
            features = c("IL1B", "TNF", "CD86", "CD163", "MRC1", "IL10", "GDF15"),
            ncol = 4)
VlnPlot(macrophages, features = "GDF15", pt.size = 0.1)
DotPlot(macrophages, features = c("LYZ", "CD68", "CD14", "CD163", "MRC1", "IL1B", "TNF", "GDF15")) +
  RotatedAxis()

# Macrophage subcluster markers (for annotation / GDF15-high subcluster)
macrophage_markers <- FindAllMarkers(macrophages, only.pos = TRUE,
                                     min.pct = 0.25, logfc.threshold = 0.25)
top5_mac <- macrophage_markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5)
DoHeatmap(macrophages, features = top5_mac$gene) + NoLegend()

# =====================================================================
# STEP 10: Publication-style plots
# =====================================================================
p1 <- SCpubr::do_DimPlot(macrophages, reduction = "umap", group.by = "seurat_clusters",
                         pt.size = 1.0, label = TRUE,
                         plot.title = "Chen macrophage subclusters")
p2 <- SCpubr::do_FeaturePlot(macrophages, features = "GDF15", pt.size = 1.0,
                             plot.title = "GDF15 expression")
cowplot::plot_grid(p1, p2, labels = c("A", "B"), ncol = 2, align = "hv")

# =====================================================================
# STEP 11: Save
# =====================================================================
setwd("/Users/alekhyamazumdar/Documents/scRNA/Chen et al")
saveRDS(macrophages, "Seurat_ChenEtAlmacrophages.rds")
saveRDS(chen_sub,    "Seurat_ChenEtAl_QCfiltered_full.rds")
