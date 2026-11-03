# =====================================================================
# Tabula Sapiens - Prostate  (TS_Prostate.h5ad)
# Robust QC  ->  cluster  ->  isolate macrophages  ->  GDF15
# ---------------------------------------------------------------------
# Counterpart to Chen_QC_Macrophages.R, adapted for the Tabula Sapiens h5ad.
# IMPORTANT difference from Chen:
#   * This h5ad exposes a "raw_counts" assay -> these ARE raw integer counts,
#     so NormalizeData IS correct here (do NOT skip it).
#   * MT genes ARE present (~37), so MT-based QC works.
# (Chen was pre-normalized with no MT genes; that mismatch is exactly why MT/
#  batch effects appeared when the two were integrated.)
#
# Fixes vs the old Sapiens script:
#   * install.packages/BiocManager lines moved to a one-off, commented block.
#   * Removed the redundant double subset() in QC (the 2nd dropped the MT rule).
#   * QC done TWO ways (MAD default + fixed) with a data-type sanity check.
#   * Macrophages identified by a marker MODULE SCORE, not hard-coded cluster 3
#     (old comments disagreed: "Monocytic - No 3" vs "M2 Macrophage - No 6").
#   * scVI UMAP preserved as its own reduction ("scvi_umap") so RunUMAP doesn't
#     overwrite it.
#   * SCpubr FeaturePlot max.cutoff clamped so it can't error (GDF15 max < 3).
#
# NOTE: I could not run this (no R / data here). Run top-to-bottom and check the
# printed messages + diagnostic plots. Written for Seurat v5 (also works on v4).
# =====================================================================

# -----------------------------------------
# STEP 0 (ONE-OFF): install packages — run once, then leave commented
# -----------------------------------------
# install.packages("BiocManager")
# BiocManager::install(c("zellkonverter", "SingleCellExperiment"))

# -----------------------------------------
# Libraries
# -----------------------------------------
library(zellkonverter)
library(SingleCellExperiment)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(data.table)
library(cowplot)
library(SCpubr)

# Helper: SCpubr::do_FeaturePlot errors if max.cutoff exceeds a feature's max.
safe_cutoff <- function(obj, feature, desired, layer = "data") {
  vals <- FetchData(obj, vars = feature, layer = layer)[[1]]
  min(desired, max(vals, na.rm = TRUE))
}

# =====================================================================
# STEP 1: Unzip + read the h5ad
# =====================================================================
base_dir <- "/Users/alekhyamazumdar/Documents/scRNA/Tabula Sapiens - Prostate"
# Unzip only if not already extracted
h5ad_path <- file.path(base_dir, "TS_Prostate.h5ad")
if (!file.exists(h5ad_path)) {
  unzip(file.path(base_dir, "TS_Prostate.h5ad.zip"), exdir = base_dir)
}
adata <- readH5AD(h5ad_path)

# Inspect what's inside so we grab the right assay / reduction (don't guess).
cat("Assays:\n");        print(assayNames(adata))
cat("Reductions:\n");    print(reducedDimNames(adata))

# =====================================================================
# STEP 2: Build the Seurat object from RAW counts
# =====================================================================
# Confirm the raw-counts assay name from the printout above; adjust if needed.
counts_assay <- if ("raw_counts" %in% assayNames(adata)) "raw_counts" else assayNames(adata)[1]
cat("Using assay:", counts_assay, "\n")

counts_matrix <- assays(adata)[[counts_assay]]
metadata      <- as.data.frame(colData(adata))

sapiens <- CreateSeuratObject(counts = counts_matrix, meta.data = metadata,
                              project = "Sapiens_Prostate",
                              min.cells = 3, min.features = 200)

# Preserve the scVI UMAP as its OWN reduction (so RunUMAP later doesn't clobber it)
if ("X_scvi_umap" %in% reducedDimNames(adata)) {
  scvi <- as.matrix(reducedDims(adata)[["X_scvi_umap"]])
  rownames(scvi) <- colnames(adata)                 # ensure cell names are set
  scvi <- scvi[colnames(sapiens), , drop = FALSE]   # align to cells kept by CreateSeuratObject
  colnames(scvi) <- paste0("scviUMAP_", seq_len(ncol(scvi)))
  sapiens[["scvi_umap"]] <- CreateDimReducObject(
    embeddings = scvi, key = "scviUMAP_", assay = DefaultAssay(sapiens)
  )
}

# =====================================================================
# STEP 3: QC metrics
# =====================================================================
n_mt <- sum(grepl("^MT-", rownames(sapiens)))
cat("Mitochondrial genes found (^MT-):", n_mt, "\n")

sapiens[["percent.mt"]]   <- PercentageFeatureSet(sapiens, pattern = "^MT-")
sapiens[["percent.ribo"]] <- PercentageFeatureSet(sapiens, pattern = "^RPS|^RPL")
sapiens[["percent.hb"]]   <- PercentageFeatureSet(sapiens, pattern = "^HB[AB]")

# Sanity check: expect RAW counts here (unlike Chen). Both should say "raw".
cat("All nFeature <= nCount (expect TRUE for raw counts):",
    all(sapiens$nFeature_RNA <= sapiens$nCount_RNA), "\n")
cat("Non-integer values present (expect FALSE for raw counts):",
    any(GetAssayData(sapiens, layer = "counts")@x %% 1 != 0), "\n")

VlnPlot(sapiens, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
        ncol = 4, pt.size = 0)
FeatureScatter(sapiens, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(sapiens, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# =====================================================================
# STEP 4: QC filtering — MAD-based (default) AND fixed, for comparison
# =====================================================================
qc_method <- "MAD"     # "MAD" (data-driven, recommended) or "fixed"

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

# percent.mt has a real distribution here, so the mito MAD is meaningful.
mt_fail <- if (mad(sapiens$percent.mt, na.rm = TRUE) > 0) {
  mad_outlier(sapiens$percent.mt, qc_nmads, "higher", use_log = FALSE)
} else rep(FALSE, ncol(sapiens))

mad_fail <-
  mad_outlier(sapiens$nFeature_RNA, qc_nmads, "both", use_log = TRUE) |
  mad_outlier(sapiens$nCount_RNA,   qc_nmads, "both", use_log = TRUE) |
  mt_fail

# Fixed thresholds (from your original script; tune to the diagnostic plots)
min_features <- 200
max_features <- 5000
max_counts   <- 40000
max_mt       <- 20     # percent
fixed_fail <-
  sapiens$nFeature_RNA < min_features |
  sapiens$nFeature_RNA > max_features |
  sapiens$nCount_RNA   > max_counts   |
  sapiens$percent.mt   > max_mt

cat(sprintf("MAD filter  : removes %d of %d cells (keeps %d)\n",
            sum(mad_fail),   length(mad_fail),   sum(!mad_fail)))
cat(sprintf("Fixed filter: removes %d of %d cells (keeps %d)\n",
            sum(fixed_fail), length(fixed_fail), sum(!fixed_fail)))

sapiens$qc_keep <- if (qc_method == "MAD") !mad_fail else !fixed_fail
sapiens_sub <- subset(sapiens, subset = qc_keep)
cat("Kept", ncol(sapiens_sub), "cells after", qc_method, "QC\n")

# =====================================================================
# STEP 5: Normalize, variable features, scale, PCA  (raw counts -> normalize)
# =====================================================================
sapiens_sub <- NormalizeData(sapiens_sub, normalization.method = "LogNormalize", scale.factor = 1e4)
sapiens_sub <- FindVariableFeatures(sapiens_sub, selection.method = "vst", nfeatures = 3000)
sapiens_sub <- ScaleData(sapiens_sub, features = rownames(sapiens_sub))
sapiens_sub <- RunPCA(sapiens_sub, features = VariableFeatures(sapiens_sub), verbose = FALSE)

ElbowPlot(sapiens_sub, ndims = 40)
n_dims <- 15   # adjust to where the elbow flattens

# =====================================================================
# STEP 6: Cluster (resolution high enough to split cell types)
# =====================================================================
sapiens_sub <- FindNeighbors(sapiens_sub, dims = 1:n_dims)
sapiens_sub <- FindClusters(sapiens_sub, resolution = 0.5)
sapiens_sub <- RunUMAP(sapiens_sub, dims = 1:n_dims)   # writes "umap"; scVI UMAP kept as "scvi_umap"
DimPlot(sapiens_sub, reduction = "umap", label = TRUE) + NoLegend()

# =====================================================================
# STEP 7: Identify the macrophage cluster(s) by marker MODULE SCORE
# =====================================================================
mac_markers <- c("CD68", "CD14", "LYZ", "CSF1R", "AIF1", "TYROBP", "FCER1G",
                 "C1QA", "C1QB", "C1QC", "FCGR3A", "MRC1", "CD163", "ITGAM")
mac_markers <- intersect(mac_markers, rownames(sapiens_sub))
sapiens_sub <- AddModuleScore(sapiens_sub, features = list(mac_markers), name = "MacScore")
mac_score_by_cluster <- sort(
  tapply(sapiens_sub$MacScore1, sapiens_sub$seurat_clusters, mean), decreasing = TRUE
)
cat("Mean macrophage score per cluster:\n"); print(round(mac_score_by_cluster, 3))

score_cutoff <- mean(mac_score_by_cluster) + sd(mac_score_by_cluster)
macrophage_clusters <- names(mac_score_by_cluster)[mac_score_by_cluster > score_cutoff]
cat("Auto-selected macrophage cluster(s):", paste(macrophage_clusters, collapse = ", "), "\n")
# If wrong after the plots below, set manually, e.g.: macrophage_clusters <- c("3")

# Confirm visually before subsetting
FeaturePlot(sapiens_sub, features = c("CD68", "CD14", "LYZ", "CSF1R", "C1QA", "MRC1"), ncol = 3)
DotPlot(sapiens_sub, features = mac_markers) + RotatedAxis()
FeaturePlot(sapiens_sub, features = c("AR", "KRT19", "KRT18", "CD3E", "PECAM1", "GDF15"), ncol = 3)

# =====================================================================
# STEP 8: Subset macrophages and re-analyze
# =====================================================================
Idents(sapiens_sub) <- "seurat_clusters"
macrophages <- subset(sapiens_sub, idents = macrophage_clusters)
cat("Macrophage cells:", ncol(macrophages), "\n")

# No NormalizeData: the subset already carries the normalized "data" layer.
macrophages <- FindVariableFeatures(macrophages, selection.method = "vst", nfeatures = 2000)
macrophages <- ScaleData(macrophages, features = rownames(macrophages))
macrophages <- RunPCA(macrophages, verbose = FALSE)
ElbowPlot(macrophages, ndims = 30)
mac_dims <- 10

macrophages <- FindNeighbors(macrophages, dims = 1:mac_dims)
macrophages <- FindClusters(macrophages, resolution = 0.2)
macrophages <- RunUMAP(macrophages, dims = 1:mac_dims)
DimPlot(macrophages, reduction = "umap", label = TRUE)

# =====================================================================
# STEP 9: GDF15 in macrophages
# =====================================================================
FeaturePlot(macrophages,
            features = c("IL1B", "TNF", "CD86", "CD163", "MRC1", "IL10", "GDF15"), ncol = 4)
VlnPlot(macrophages, features = "GDF15", pt.size = 0.1, layer = "data")
DotPlot(macrophages, features = c("LYZ", "CD68", "CD14", "CD163", "MRC1", "IL1B", "TNF", "GDF15")) +
  RotatedAxis()

macrophage_markers <- FindAllMarkers(macrophages, only.pos = TRUE,
                                     min.pct = 0.25, logfc.threshold = 0.25)
top5_mac <- macrophage_markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5)
DoHeatmap(macrophages, features = top5_mac$gene) + NoLegend()

# =====================================================================
# STEP 10: Publication-style plots
# =====================================================================
p1 <- SCpubr::do_DimPlot(macrophages, reduction = "umap", group.by = "seurat_clusters",
                         pt.size = 1.0, label = TRUE,
                         plot.title = "Tabula Sapiens macrophage subclusters")
p2 <- SCpubr::do_FeaturePlot(macrophages, features = "GDF15", pt.size = 1.0,
                             min.cutoff = 0, max.cutoff = safe_cutoff(macrophages, "GDF15", 3.0),
                             plot.title = "GDF15 expression")
cowplot::plot_grid(p1, p2, labels = c("A", "B"), ncol = 2, align = "hv")

# =====================================================================
# STEP 11: Save
# =====================================================================
setwd(base_dir)
saveRDS(macrophages, "SapiensProstateMacrophages.rds")
saveRDS(sapiens_sub, "SapiensProstate_QCfiltered_full.rds")
