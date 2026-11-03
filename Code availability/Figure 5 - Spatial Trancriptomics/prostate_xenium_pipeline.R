# =============================================================================
# Prostate Xenium Spatial Analysis — Takehara et al. (Immunity 2026)
# Single-file pipeline. Open in RStudio and step through section by section.
# Use Code -> Show Document Outline (Ctrl/Cmd+Shift+O) to navigate sections.
# =============================================================================
#
# Data:    /Users/alekhyamazumdar/Documents/Xenium/GSE296275_RAW   (4 samples, GEO GSE296275)
# Output:  /Users/alekhyamazumdar/Documents/Xenium/Output          (RDS, figures, CSVs)
#
# Sections:
#   00. Install packages (run once)
#   01. Load Xenium data (build Seurat object from .mtx + transcripts.parquet)
#   02. QC, normalization, clustering, cell-type annotation
#   03. Spatial niche identification
#   04. CD8+ T cell heterogeneity
#   05. Cytokine niche mapping (IL-15 / IL-7 / TGF-beta)
#   06. Neighborhood / cell-type co-localization
#   07. Macrophage analysis
#   08. Macrophage deep dive: DE heatmap, GO enrichment, proximity pseudobulk
#
# Each section saves an .rds checkpoint into OUT_DIR so you can resume mid-pipeline.
# =============================================================================


# 00. INSTALL PACKAGES ---------------------------------------------------------
# Run this section once. It is idempotent — re-running it is harmless.

if (FALSE) {  # set to TRUE the first time
  cran_pkgs <- c(
    "Seurat", "SeuratObject", "Matrix", "data.table", "arrow",
    "R.utils", "ggplot2", "patchwork", "dplyr", "tidyr",
    "RColorBrewer", "viridis", "future", "future.apply", "harmony",
    "presto", "FNN"
  )
  new_cran <- cran_pkgs[!cran_pkgs %in% installed.packages()[, "Package"]]
  if (length(new_cran)) install.packages(new_cran)

  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  bioc_pkgs <- c("SpatialExperiment", "SingleCellExperiment", "scran",
                 "scater", "BiocParallel", "rhdf5",
                 # Section 08 dependencies (macrophage deep dive)
                 "DESeq2", "clusterProfiler", "org.Mm.eg.db",
                 "ComplexHeatmap")
  new_bioc <- bioc_pkgs[!bioc_pkgs %in% installed.packages()[, "Package"]]
  if (length(new_bioc)) BiocManager::install(new_bioc, update = FALSE, ask = FALSE)

  # CRAN extras for section 08
  extra_cran <- c("pheatmap", "ggrepel")
  new_extra <- extra_cran[!extra_cran %in% installed.packages()[, "Package"]]
  if (length(new_extra)) install.packages(new_extra)
}


# GLOBAL SETUP -----------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(arrow)
  library(data.table)
  library(R.utils)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(harmony)
  library(future)
  library(FNN)
})

# Run sequentially — Xenium objects are large enough that future's parallel
# data shipping blows past the default globals budget.
plan("sequential")
options(future.globals.maxSize = 16 * 1024^3)  # 16 GB

DATA_DIR <- "/Users/alekhyamazumdar/Documents/Xenium/GSE296275_RAW"
OUT_DIR  <- "/Users/alekhyamazumdar/Documents/Xenium/Output"
TMP_DIR  <- file.path(tempdir(), "xenium_parquet")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TMP_DIR, showWarnings = FALSE, recursive = TRUE)

# Helper: print a plot to RStudio's Plots pane AND save it to OUT_DIR.
# Use this everywhere instead of bare ggsave().
save_plot <- function(p, filename, width = 8, height = 6) {
  print(p)
  ggsave(file.path(OUT_DIR, filename), p, width = width, height = height)
  invisible(p)
}

# Helper: auto-annotate clusters by picking, for each cluster, the signature
# whose mean module score is highest. Returns a named character vector
# (cluster id -> label), and writes a CSV + heatmap so you can review.
auto_annotate_clusters <- function(meta, cluster_col, score_prefix = "score_",
                                   csv_name = "cluster_signature_scores.csv",
                                   heatmap_name = "cluster_signature_heatmap.png",
                                   strip_prefix = "score_") {
  score_cols <- grep(paste0("^", score_prefix), colnames(meta), value = TRUE)
  # Filter to numeric columns only — guards against character columns like
  # 'mac_subset' that may share the prefix from a previous pipeline run.
  if (length(score_cols) > 0) {
    score_cols <- score_cols[
      vapply(meta[, score_cols, drop = FALSE], is.numeric, logical(1))
    ]
  }
  if (length(score_cols) == 0) {
    stop("No numeric score columns matching '", score_prefix, "' found in metadata.")
  }
  mat <- meta %>%
    dplyr::select(dplyr::all_of(c(cluster_col, score_cols))) %>%
    dplyr::group_by(.data[[cluster_col]]) %>%
    dplyr::summarise(across(everything(), mean), .groups = "drop") %>%
    as.data.frame()
  rownames(mat) <- as.character(mat[[cluster_col]])
  mat[[cluster_col]] <- NULL
  write.csv(mat, file.path(OUT_DIR, csv_name), row.names = TRUE)

  # AddModuleScore appends a digit to whatever 'name' you give it, so a column
  # asked to be 'score_CD8_Tcell' becomes 'score_CD8_Tcell1'. Strip the prefix
  # AND the trailing digits to get a clean label.
  clean_label <- function(x) {
    x <- sub(paste0("^", strip_prefix), "", x)
    x <- sub("\\d+$", "", x)
    x
  }
  df_heat <- mat %>%
    tibble::rownames_to_column("cluster") %>%
    tidyr::pivot_longer(-cluster, names_to = "signature", values_to = "score") %>%
    dplyr::mutate(signature = clean_label(signature),
                  cluster   = factor(cluster,
                                     levels = as.character(sort(
                                       suppressWarnings(as.numeric(unique(cluster)))))))
  p_heat <- ggplot(df_heat, aes(x = signature, y = cluster, fill = score)) +
    geom_tile() +
    scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick",
                         midpoint = 0, name = "Mean\nscore") +
    labs(x = "Signature", y = "Cluster") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank())
  save_plot(p_heat, heatmap_name, width = max(7, ncol(mat) * 0.5),
            height = max(5, nrow(mat) * 0.35))

  auto <- apply(mat, 1, function(row) clean_label(names(row)[which.max(row)]))
  message("Auto-assigned cluster labels (review the heatmap and override below if needed):")
  print(auto)
  auto
}

sample_meta <- tibble::tribble(
  ~gsm,         ~sample_id,        ~timepoint, ~replicate,
  "GSM8967906", "D7_Prostate_1",   "7dpi",     "rep1",
  "GSM8967907", "D7_Prostate_2",   "7dpi",     "rep2",
  "GSM8967908", "D30_Prostate_1",  "30dpi",    "rep1",
  "GSM8967909", "D30_Prostate_2",  "30dpi",    "rep2"
)


# 01. LOAD XENIUM DATA ---------------------------------------------------------
# GEO deposits Xenium output flattened (no cells.csv.gz), so we build the
# Seurat object manually:
#   * read cell x gene matrix from .mtx + barcodes + features
#   * derive centroids from transcripts.parquet (median x/y per cell_id)
#   * wrap in a Seurat v5 object with one Centroids FOV per sample

.fname <- function(gsm, sample, suffix) {
  file.path(DATA_DIR, paste0(gsm, "_", sample, "_", suffix))
}

read_xenium_matrix <- function(gsm, sample) {
  mtx_path <- .fname(gsm, sample, "cell_feature_matrix_matrix.mtx.gz")
  bc_path  <- .fname(gsm, sample, "cell_feature_matrix_barcodes.tsv.gz")
  ft_path  <- .fname(gsm, sample, "cell_feature_matrix_features.tsv.gz")
  message("  reading matrix: ", basename(mtx_path))
  m  <- readMM(gzfile(mtx_path))
  bc <- readLines(gzfile(bc_path))
  ft <- read.table(gzfile(ft_path), sep = "\t", header = FALSE,
                   stringsAsFactors = FALSE)
  rownames(m) <- make.unique(ft$V2)
  colnames(m) <- bc
  m <- as(m, "CsparseMatrix")
  keep <- ft$V3 == "Gene Expression"
  m[keep, , drop = FALSE]
}

derive_centroids <- function(gsm, sample, valid_cells) {
  pq_gz <- .fname(gsm, sample, "transcripts.parquet.gz")
  pq    <- file.path(TMP_DIR, paste0(sample, "_transcripts.parquet"))
  if (!file.exists(pq)) {
    message("  decompressing parquet: ", basename(pq_gz))
    gunzip(pq_gz, destname = pq, remove = FALSE, overwrite = TRUE)
  }
  message("  reading parquet schema")
  ds <- open_dataset(pq, format = "parquet")
  print(ds$schema)
  cell_col <- "cell_id"; x_col <- "x_location"; y_col <- "y_location"
  message("  computing centroids")
  cols <- c(cell_col, x_col, y_col)
  if ("is_gene" %in% names(ds$schema)) {
    txs <- ds |>
      dplyr::filter(is_gene == TRUE,
                    .data[[cell_col]] %in% valid_cells) |>
      dplyr::select(dplyr::all_of(cols)) |>
      dplyr::collect()
  } else {
    txs <- ds |>
      dplyr::filter(.data[[cell_col]] %in% valid_cells) |>
      dplyr::select(dplyr::all_of(cols)) |>
      dplyr::collect()
  }
  setDT(txs)
  setnames(txs, c(cell_col, x_col, y_col), c("cell", "x", "y"))
  centroids <- txs[, .(x = median(x, na.rm = TRUE),
                       y = median(y, na.rm = TRUE)), by = cell]
  as.data.frame(centroids)
}

build_one <- function(gsm, sample_id, timepoint, replicate) {
  message("=== ", sample_id, " (", gsm, ") ===")
  m <- read_xenium_matrix(gsm, sample_id)
  cents <- derive_centroids(gsm, sample_id, valid_cells = colnames(m))
  keep_cells <- intersect(colnames(m), cents$cell)
  m     <- m[, keep_cells, drop = FALSE]
  cents <- cents[match(keep_cells, cents$cell), ]
  obj <- CreateSeuratObject(counts = m, assay = "Xenium", project = sample_id)
  obj$sample_id <- sample_id
  obj$timepoint <- timepoint
  obj$replicate <- replicate
  cents_obj <- CreateCentroids(
    coords = data.frame(x = cents$x, y = cents$y, cell = cents$cell)
  )
  fov <- CreateFOV(coords = list(centroids = cents_obj),
                   type = "centroids", assay = "Xenium",
                   key = paste0(make.names(sample_id), "_"))
  obj[[sample_id]] <- fov
  obj
}

obj_list <- mapply(build_one, sample_meta$gsm, sample_meta$sample_id,
                   sample_meta$timepoint, sample_meta$replicate,
                   SIMPLIFY = FALSE)
names(obj_list) <- sample_meta$sample_id

xen <- merge(obj_list[[1]], y = obj_list[-1], add.cell.ids = names(obj_list))
DefaultAssay(xen) <- "Xenium"

print(xen)
cat("Total cells:    ", ncol(xen), "\n")
cat("Genes in panel: ", nrow(xen), "\n")
print(table(xen$sample_id, xen$timepoint))

saveRDS(xen, file.path(OUT_DIR, "01_xenium_raw.rds"))
rm(obj_list); gc()


# 02. QC, NORMALIZATION, CLUSTERING --------------------------------------------
# Resume point: xen <- readRDS(file.path(OUT_DIR, "01_xenium_raw.rds"))

# Join per-sample counts layers into one (Seurat v5 keeps them split after merge)
xen <- JoinLayers(xen, assay = "Xenium")

# QC violins (nCount/nFeature already populated by CreateSeuratObject)
p_qc <- VlnPlot(xen, features = c("nCount_Xenium", "nFeature_Xenium"),
                group.by = "sample_id", pt.size = 0) +
  plot_layout(ncol = 1)
save_plot(p_qc, "qc_violins.png", 9, 6)

# Filter very-low-count cells (segmentation artifacts)
xen <- subset(xen, subset = nCount_Xenium >= 10 & nFeature_Xenium >= 5)
cat("Cells after QC: ", ncol(xen), "\n")

# Normalize + scale (use all 541 panel genes — no variable feature subsetting)
xen <- NormalizeData(xen, assay = "Xenium", verbose = FALSE)
xen <- ScaleData(xen, assay = "Xenium", features = rownames(xen),
                 verbose = FALSE)

# PCA on all genes, then Harmony to integrate across the 4 samples
xen <- RunPCA(xen, assay = "Xenium", features = rownames(xen),
              npcs = 30, verbose = FALSE)
xen <- RunHarmony(xen, group.by.vars = "sample_id", reduction = "pca",
                  reduction.save = "harmony", assay.use = "Xenium",
                  verbose = FALSE)

# UMAP + Leiden clustering (algorithm = 4)
xen <- RunUMAP(xen, reduction = "harmony", dims = 1:30, verbose = FALSE)
xen <- FindNeighbors(xen, reduction = "harmony", dims = 1:30, verbose = FALSE)
xen <- FindClusters(xen, resolution = 0.5, algorithm = 4, verbose = FALSE)

p_umap <- DimPlot(xen, group.by = "seurat_clusters", label = TRUE, raster = FALSE) +
  DimPlot(xen, group.by = "sample_id", raster = FALSE)
save_plot(p_umap, "umap_clusters.png", 12, 5)

# Marker-based annotation. Lists are pre-filtered to the 480-gene panel of
# this dataset (Takehara et al.). Genes commonly assumed for these cell types
# but absent from the panel (e.g. Cd68, Krt14, Nkx3-1, Pbsn, Trp63, Mrc1,
# Cd163, C1qa/b, Mpo, S100a9, Tigit, Tox, Eomes, Hobit) have been removed.
marker_panels <- list(
  Epithelial         = c("Krt8", "Krt5", "Krt18", "Cdh1", "Epcam"),
  Luminal_epithelial = c("Krt8", "Krt18", "Ar", "Adgrg2"),
  Basal_epithelial   = c("Krt5", "Prom1"),
  Stem_epithelial    = c("Lgr5", "Lrig1", "Prom1", "Olfm4"),
  Smooth_muscle      = c("Acta2", "Myh11", "Tagln2", "Cnn2"),
  Fibroblast         = c("Col1a1", "Col3a1", "Pdgfra", "Dpt", "Tcf21"),
  Endothelial        = c("Pecam1", "Cdh5", "Cd34", "Vwf"),
  CD8_Tcell          = c("Cd8a", "Cd8b1", "Cd3e"),
  CD4_Tcell          = c("Cd4", "Cd3e", "Foxp3"),
  B_cell             = c("Cd19", "Cd79a"),
  Macrophage         = c("Adgre1", "Csf1r", "Lyz2", "Aif1", "Mertk"),
  Dendritic          = c("Itgax", "Xcr1", "Clec9a", "H2-Ab1"),
  NK_cell            = c("Ncr1", "Gzma", "Nkg7"),
  Neutrophil         = c("S100a8", "Ly6g"),
  Trm_signature      = c("Cd69", "Itgae", "Runx3", "Prdm1", "Notch1")
)

# Defensive: drop any genes that aren't in this object (handles future panels
# with different content). Warn so you know what was dropped.
.report_missing <- function(panels, available, label) {
  for (n in names(panels)) {
    miss <- setdiff(panels[[n]], available)
    if (length(miss) > 0) {
      message(label, " :: '", n, "' missing from data: ",
              paste(miss, collapse = ", "))
    }
  }
}
.filter_panels <- function(panels, available) {
  out <- lapply(panels, function(g) intersect(g, available))
  out[sapply(out, length) > 0]   # drop empty panels
}
.report_missing(marker_panels, rownames(xen), "marker_panels")
marker_panels <- .filter_panels(marker_panels, rownames(xen))

DefaultAssay(xen) <- "Xenium"
# AddModuleScore defaults (nbin=24, ctrl=100) are tuned for whole-transcriptome
# scRNA-seq with ~20k genes. On a 480-gene panel those defaults fail with
# "cannot take a sample larger than the population". Use small panel-friendly
# parameters instead.
for (n in names(marker_panels)) {
  feats <- intersect(marker_panels[[n]], rownames(xen))
  if (length(feats) >= 1) {
    xen <- AddModuleScore(xen, features = list(feats),
                          name = paste0("score_", n), assay = "Xenium",
                          nbin = 10, ctrl = 10)
  }
}

all_markers <- unique(unlist(marker_panels))
all_markers <- intersect(all_markers, rownames(xen))
p_dot <- DotPlot(xen, features = all_markers, group.by = "seurat_clusters") +
  RotatedAxis()
save_plot(p_dot, "marker_dotplot.png", 14, 6)

cluster_markers <- FindAllMarkers(xen, only.pos = TRUE, min.pct = 0.25,
                                  logfc.threshold = 0.5)
write.csv(cluster_markers, file.path(OUT_DIR, "cluster_markers.csv"),
          row.names = FALSE)

# Auto-annotate clusters by signature score: each cluster gets the label of
# the marker_panels signature whose mean module score is highest in that
# cluster. Writes:
#   - cluster_signature_scores.csv : full cluster x signature mean-score matrix
#   - cluster_signature_heatmap.png : visual heatmap of the same
# Open the heatmap to see whether the auto-assignment matches your biology.
cluster_annotation <- auto_annotate_clusters(
  xen@meta.data, cluster_col = "seurat_clusters",
  csv_name     = "cluster_signature_scores.csv",
  heatmap_name = "cluster_signature_heatmap.png"
)

# >>> OPTIONAL: review cluster_signature_heatmap.png + marker_dotplot.png +
# >>> cluster_markers.csv, then override individual cluster labels here.
# >>> Examples (uncomment and edit as needed):
#   cluster_annotation["3"]  <- "Smooth_muscle"
#   cluster_annotation["7"]  <- "CD4_Tcell"
#   cluster_annotation["10"] <- "B_cell"

xen$cell_type <- ifelse(
  as.character(xen$seurat_clusters) %in% names(cluster_annotation),
  cluster_annotation[as.character(xen$seurat_clusters)],
  paste0("cluster_", xen$seurat_clusters)
)

p_celltype <- DimPlot(xen, group.by = "cell_type", label = TRUE, raster = FALSE)
save_plot(p_celltype, "umap_celltype.png", 9, 7)

for (s in unique(xen$sample_id)) {
  p_sp <- ImageDimPlot(xen, fov = s, group.by = "cell_type",
                       size = 0.5, dark.background = FALSE)
  save_plot(p_sp, paste0("spatial_celltype_", s, ".png"), 8, 8)
}

saveRDS(xen, file.path(OUT_DIR, "02_xenium_annotated.rds"))


# 03. NICHE IDENTIFICATION -----------------------------------------------------
# Resume: xen <- readRDS(file.path(OUT_DIR, "02_xenium_annotated.rds"))

niche_results <- list()
for (s in unique(xen$sample_id)) {
  message("Niche assay for ", s)
  sub <- subset(xen, subset = sample_id == s)
  sub <- BuildNicheAssay(object = sub, fov = s, group.by = "cell_type",
                         niches.k = 5, neighbors.k = 30)
  niche_results[[s]] <- sub@meta.data[, c("niches"), drop = FALSE]
}
niche_df <- do.call(rbind, lapply(niche_results, function(x) {
  data.frame(cell = rownames(x), niche = x$niches)
}))
xen$niche <- niche_df$niche[match(colnames(xen), niche_df$cell)]

# Niche composition heatmap
comp <- xen@meta.data %>%
  dplyr::filter(!is.na(niche)) %>%
  dplyr::count(niche, cell_type) %>%
  dplyr::group_by(niche) %>%
  dplyr::mutate(frac = n / sum(n)) %>%
  dplyr::ungroup()
p_comp <- ggplot(comp, aes(x = factor(niche), y = cell_type, fill = frac)) +
  geom_tile() + scale_fill_viridis_c() +
  labs(x = "Niche", y = "Cell type", fill = "Fraction\nin niche") +
  theme_minimal(base_size = 11) + theme(panel.grid = element_blank())
save_plot(p_comp, "niche_composition.png", 7, 6)

for (s in unique(xen$sample_id)) {
  p_sp <- ImageDimPlot(xen, fov = s, group.by = "niche", size = 0.5,
                       dark.background = FALSE) + ggtitle(s)
  save_plot(p_sp, paste0("spatial_niche_", s, ".png"), 8, 8)
}

# Auto-name each niche by its dominant cell type (highest fraction).
niche_top_celltype <- comp %>%
  dplyr::group_by(niche) %>%
  dplyr::slice_max(frac, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()
niche_names <- setNames(
  paste0("niche_", niche_top_celltype$niche, "_",
         tolower(niche_top_celltype$cell_type)),
  as.character(niche_top_celltype$niche)
)
message("Auto-named niches:")
print(niche_names)

# >>> OPTIONAL: override individual niche names after viewing the heatmap.
# >>> Examples:
#   niche_names["1"] <- "epithelial_luminal"
#   niche_names["3"] <- "stromal_fibroblast"

xen$niche_label <- ifelse(
  !is.na(xen$niche) & as.character(xen$niche) %in% names(niche_names),
  niche_names[as.character(xen$niche)], NA_character_
)

saveRDS(xen, file.path(OUT_DIR, "03_xenium_niches.rds"))

# 1) Re-source the GLOBAL SETUP block (lines 47-138 of the script).
#    This redefines save_plot() and auto_annotate_clusters() with the fix.
#    In RStudio: select those lines and hit Cmd/Ctrl+Enter, OR use:
source("/Users/alekhyamazumdar/Library/Application Support/Claude/local-agent-mode-sessions/2ac14bbb-1733-464b-a511-2c93036d083a/cbccebac-ddf4-4f09-a413-9819c2fafc66/local_75a6f55d-bde8-4ca8-9214-98e3cbc9bdd8/outputs/prostate_xenium_pipeline.R", echo = FALSE)
# (sourcing the whole file is safe; it'll re-define functions but the heavy
#  data-loading/clustering code only runs if you also reach those lines —
#  if you don't want that, just select-and-run lines ~47-138 in the editor.)

# 2) Confirm the fix is live:
body(auto_annotate_clusters)
# You should see "clean_label" defined inside the function body. If you don't,
# the new version isn't loaded yet.

# 3) Now re-run just the auto-annotation call (you don't need to redo PCA/UMAP):
cluster_annotation <- auto_annotate_clusters(
  xen@meta.data, cluster_col = "seurat_clusters",
  csv_name     = "cluster_signature_scores.csv",
  heatmap_name = "cluster_signature_heatmap.png"
)
print(cluster_annotation)
# Labels should now be clean: "Epithelial", "CD8_Tcell", "Macrophage", etc.
# (No trailing 1's.)

# 4) Re-apply to the object and re-save:
xen$cell_type <- ifelse(
  as.character(xen$seurat_clusters) %in% names(cluster_annotation),
  cluster_annotation[as.character(xen$seurat_clusters)],
  paste0("cluster_", xen$seurat_clusters)
)
saveRDS(xen, file.path(OUT_DIR, "02_xenium_annotated.rds"))
table(xen$cell_type)   # sanity check


# 04. CD8+ T CELL HETEROGENEITY ------------------------------------------------
# Resume: xen <- readRDS(file.path(OUT_DIR, "03_xenium_niches.rds"))

tc <- subset(xen, subset = cell_type == "CD8_Tcell")
cat("CD8+ T cells: ", ncol(tc), "\n")

DefaultAssay(tc) <- "Xenium"
tc <- ScaleData(tc, assay = "Xenium", features = rownames(tc), verbose = FALSE)
tc <- RunPCA(tc, assay = "Xenium", features = rownames(tc), npcs = 20, verbose = FALSE)
tc <- RunUMAP(tc, dims = 1:20, verbose = FALSE)
tc <- FindNeighbors(tc, dims = 1:20, verbose = FALSE)
tc <- FindClusters(tc, resolution = 0.4, algorithm = 4, verbose = FALSE)

# Trm signatures, pre-filtered to the 480-gene panel.
# Removed (not on panel): Hobit, Rgs1, Gzmk, Prf1, Lef1, Bach2, Havcr2, Tox,
#   Tigit, Smad7, Skil, Bcl2, Mcl1, Eomes, Bhlhe40.
# IL-15 responsiveness is now indexed by the IL-2/IL-15 receptor beta chain
# (Il2rb, shared receptor) plus Klrg1 (effector marker amplified by IL-15)
# and Il15ra itself.
trm_signatures <- list(
  Trm_core       = c("Cd69", "Itgae", "Runx3", "Prdm1", "Notch1"),
  Trm_cytotoxic  = c("Gzmb", "Gzma", "Ifng", "Klrg1", "Nkg7"),
  Trm_progenitor = c("Tcf7", "Sell", "Ccr7"),
  Trm_exhausted  = c("Pdcd1", "Lag3"),
  Trm_tgfb_resp  = c("Itgae", "Itgb7"),
  Trm_il15_resp  = c("Il15ra", "Il2rb", "Klrg1")
)
.report_missing(trm_signatures, rownames(tc), "trm_signatures")
trm_signatures <- .filter_panels(trm_signatures, rownames(tc))
for (n in names(trm_signatures)) {
  feats <- intersect(trm_signatures[[n]], rownames(tc))
  if (length(feats) > 0) {
    tc <- AddModuleScore(tc, features = list(feats),
                         name = paste0("score_", n), assay = "Xenium",
                         nbin = 10, ctrl = 10)
  }
}

# Use the columns we just added (avoids picking up upstream score_* columns
# carried over from section 02, and also avoids any character columns sharing
# the prefix from a previous run).
score_cols <- paste0("score_", names(trm_signatures), "1")
score_cols <- intersect(score_cols, colnames(tc@meta.data))
score_cols <- score_cols[
  vapply(tc@meta.data[, score_cols, drop = FALSE], is.numeric, logical(1))
]
p_scores <- FeaturePlot(tc, features = score_cols, ncol = 3, raster = FALSE) &
  scale_color_viridis_c()
save_plot(p_scores, "tcell_signature_scores.png", 12, 8)

tc_markers <- FindAllMarkers(tc, only.pos = TRUE, min.pct = 0.2,
                             logfc.threshold = 0.4)
write.csv(tc_markers, file.path(OUT_DIR, "tcell_subcluster_markers.csv"),
          row.names = FALSE)

# Auto-annotate T cell subclusters by Trm signature score.
tc_annotation <- auto_annotate_clusters(
  tc@meta.data, cluster_col = "seurat_clusters",
  csv_name     = "tcell_signature_scores_by_cluster.csv",
  heatmap_name = "tcell_signature_heatmap.png"
)

# >>> OPTIONAL: override individual T cell subset labels after viewing the heatmap.
# >>> Examples:
#   tc_annotation["1"] <- "Trm_cytotoxic"
#   tc_annotation["3"] <- "Trm_progenitor"

tc$tcell_subset <- ifelse(
  as.character(tc$seurat_clusters) %in% names(tc_annotation),
  tc_annotation[as.character(tc$seurat_clusters)],
  paste0("tcell_cluster_", tc$seurat_clusters)
)

xen$tcell_subset <- NA_character_
xen$tcell_subset[colnames(xen) %in% colnames(tc)] <- tc$tcell_subset[match(
  colnames(xen)[colnames(xen) %in% colnames(tc)], colnames(tc))]

df <- xen@meta.data %>%
  dplyr::filter(!is.na(tcell_subset), !is.na(niche_label)) %>%
  dplyr::count(tcell_subset, niche_label) %>%
  dplyr::group_by(niche_label) %>%
  dplyr::mutate(frac = n / sum(n))
p_dist <- ggplot(df, aes(x = niche_label, y = frac, fill = tcell_subset)) +
  geom_col(position = "stack") +
  labs(x = "Niche", y = "Fraction of T cells", fill = "T cell subset") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p_dist, "tcell_subset_by_niche.png", 8, 6)

for (s in unique(xen$sample_id)) {
  p_sp <- ImageDimPlot(
    xen, fov = s, group.by = "tcell_subset",
    cells = WhichCells(xen, expression = !is.na(tcell_subset)),
    size = 0.7, dark.background = FALSE
  ) + ggtitle(paste0("CD8+ T cell subsets — ", s))
  save_plot(p_sp, paste0("spatial_tcell_subset_", s, ".png"), 8, 8)
}

df_time <- xen@meta.data %>%
  dplyr::filter(!is.na(tcell_subset)) %>%
  dplyr::count(timepoint, tcell_subset) %>%
  dplyr::group_by(timepoint) %>%
  dplyr::mutate(frac = n / sum(n))
p_time <- ggplot(df_time,
                 aes(x = timepoint, y = frac, fill = tcell_subset)) +
  geom_col(position = "stack") +
  labs(x = "Timepoint", y = "Fraction of CD8+ T cells", fill = "Subset") +
  theme_minimal(base_size = 11)
save_plot(p_time, "tcell_subset_by_timepoint.png", 6, 6)

saveRDS(tc,  file.path(OUT_DIR, "04_tcell_subset.rds"))
saveRDS(xen, file.path(OUT_DIR, "04_xenium_with_tcells.rds"))


# 05. CYTOKINE NICHE MAPPING ---------------------------------------------------
# Resume: xen <- readRDS(file.path(OUT_DIR, "04_xenium_with_tcells.rds"))

cytokines <- list(
  IL15 = c("Il15", "Il15ra"),
  IL7  = c("Il7",  "Il7r"),
  TGFb = c("Tgfb1", "Tgfb2", "Tgfb3", "Tgfbr1", "Tgfbr2")
)
cytokines <- lapply(cytokines, function(g) intersect(g, rownames(xen)))
print(cytokines)

DefaultAssay(xen) <- "Xenium"
for (s in unique(xen$sample_id)) {
  for (cyt in names(cytokines)) {
    feats <- cytokines[[cyt]]
    if (length(feats) == 0) next
    p_sp <- ImageFeaturePlot(xen, fov = s, features = feats,
                             cols = c("lightgrey", "red"), size = 0.5) +
      plot_annotation(title = paste0(cyt, " — ", s))
    save_plot(p_sp, paste0("spatial_", cyt, "_", s, ".png"),
              4 * length(feats), 4)
  }
}

# Top-decile cytokine "source" cells per cytokine
source_summary <- list()
for (cyt in names(cytokines)) {
  feats <- cytokines[[cyt]]
  if (length(feats) == 0) next
  expr <- GetAssayData(xen, assay = "Xenium", layer = "data")[feats, , drop = FALSE]
  agg <- if (length(feats) > 1) Matrix::colMeans(expr) else expr[1, ]
  thr <- quantile(agg, 0.90, na.rm = TRUE)
  is_src <- agg >= thr
  xen[[paste0("source_", cyt)]] <- is_src
  src_celltypes <- xen$cell_type[is_src]
  source_summary[[cyt]] <- as.data.frame(table(src_celltypes)) %>%
    dplyr::rename(cell_type = src_celltypes, n = Freq) %>%
    dplyr::mutate(frac = n / sum(n), cytokine = cyt)
}
src_df <- bind_rows(source_summary)
p_src <- ggplot(src_df, aes(x = cytokine, y = frac, fill = cell_type)) +
  geom_col(position = "stack") +
  labs(x = "Cytokine", y = "Fraction of source cells", fill = "Cell type") +
  theme_minimal(base_size = 11)
save_plot(p_src, "cytokine_source_celltypes.png", 7, 6)

# Distance from each CD8+ T cell to nearest cytokine source per sample
distance_records <- list()
for (s in unique(xen$sample_id)) {
  sub <- subset(xen, subset = sample_id == s)
  fov <- sub@images[[s]]
  if (is.null(fov)) next
  coords <- GetTissueCoordinates(fov)
  rownames(coords) <- coords$cell

  tc_idx   <- sub$cell_type == "CD8_Tcell"
  tc_cells <- colnames(sub)[tc_idx]
  if (length(tc_cells) == 0) next
  tc_xy <- coords[tc_cells, c("x", "y"), drop = FALSE]

  for (cyt in names(cytokines)) {
    src_col   <- paste0("source_", cyt)
    src_cells <- colnames(sub)[sub[[src_col]][, 1]]
    if (length(src_cells) < 5) next
    src_xy <- coords[src_cells, c("x", "y"), drop = FALSE]
    nn <- get.knnx(as.matrix(src_xy), as.matrix(tc_xy), k = 1)
    distance_records[[length(distance_records) + 1]] <- data.frame(
      cell         = tc_cells,
      sample_id    = s,
      timepoint    = sub$timepoint[1],
      tcell_subset = sub$tcell_subset[tc_idx],
      cytokine     = cyt,
      distance_um  = nn$nn.dist[, 1]
    )
  }
}
dist_df <- bind_rows(distance_records)
write.csv(dist_df, file.path(OUT_DIR, "tcell_cytokine_distances.csv"),
          row.names = FALSE)

p_dist2 <- ggplot(dist_df,
                  aes(x = tcell_subset, y = distance_um, fill = cytokine)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ cytokine, scales = "free_y") +
  scale_y_log10() +
  labs(x = "T cell subset",
       y = "Distance to nearest source cell (μm, log10)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none")
save_plot(p_dist2, "tcell_cytokine_distances.png", 10, 5)

niche_cyt <- xen@meta.data %>%
  dplyr::filter(!is.na(niche_label)) %>%
  dplyr::select(niche_label, starts_with("score_")) %>%
  dplyr::group_by(niche_label) %>%
  dplyr::summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")
write.csv(niche_cyt, file.path(OUT_DIR, "niche_cytokine_scores.csv"),
          row.names = FALSE)

saveRDS(xen, file.path(OUT_DIR, "05_xenium_cytokine.rds"))


# 06. NEIGHBORHOOD CO-LOCALIZATION ---------------------------------------------
# Resume: xen <- readRDS(file.path(OUT_DIR, "05_xenium_cytokine.rds"))
# Permutation-based neighborhood enrichment z-scores between cell-type pairs.

K <- 30
nbhd_records <- list()
for (s in unique(xen$sample_id)) {
  sub <- subset(xen, subset = sample_id == s)
  fov <- sub@images[[s]]
  if (is.null(fov)) next
  coords <- GetTissueCoordinates(fov)
  rownames(coords) <- coords$cell
  xy <- coords[colnames(sub), c("x", "y"), drop = FALSE]

  ct  <- sub$cell_type
  nn  <- get.knn(as.matrix(xy), k = K)$nn.index
  nbr_ct <- matrix(ct[nn], nrow = nrow(nn))

  ct_levels <- sort(unique(ct))
  obs <- matrix(0, nrow = length(ct_levels), ncol = length(ct_levels),
                dimnames = list(ct_levels, ct_levels))
  for (i in seq_along(ct)) {
    nbr <- nbr_ct[i, ]
    obs[ct[i], ] <- obs[ct[i], ] +
      tabulate(match(nbr, ct_levels), nbins = length(ct_levels))
  }
  set.seed(42); n_perm <- 100
  perm <- array(0, dim = c(length(ct_levels), length(ct_levels), n_perm))
  for (p in seq_len(n_perm)) {
    ctp <- sample(ct)
    nbrp <- matrix(ctp[nn], nrow = nrow(nn))
    op <- matrix(0, nrow = length(ct_levels), ncol = length(ct_levels),
                 dimnames = list(ct_levels, ct_levels))
    for (i in seq_along(ctp)) {
      op[ctp[i], ] <- op[ctp[i], ] +
        tabulate(match(nbrp[i, ], ct_levels), nbins = length(ct_levels))
    }
    perm[, , p] <- op
  }
  mu  <- apply(perm, c(1, 2), mean)
  sd_ <- apply(perm, c(1, 2), sd)
  z   <- (obs - mu) / pmax(sd_, 1e-9)
  long <- as.data.frame(as.table(z)) %>%
    setNames(c("cell_type_a", "cell_type_b", "z")) %>%
    dplyr::mutate(sample_id = s)
  nbhd_records[[s]] <- long
}
nbhd_df <- bind_rows(nbhd_records)
write.csv(nbhd_df, file.path(OUT_DIR, "neighborhood_zscores.csv"),
          row.names = FALSE)

nbhd_avg <- nbhd_df %>%
  dplyr::left_join(unique(xen@meta.data[, c("sample_id", "timepoint")]),
                   by = "sample_id") %>%
  dplyr::group_by(timepoint, cell_type_a, cell_type_b) %>%
  dplyr::summarise(z = mean(z, na.rm = TRUE), .groups = "drop")

p_nbhd <- ggplot(nbhd_avg, aes(x = cell_type_a, y = cell_type_b, fill = z)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       name = "Enrichment\nz-score") +
  facet_wrap(~ timepoint) + labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())
save_plot(p_nbhd, "neighborhood_enrichment.png", 11, 5)

cd8_partners <- nbhd_df %>%
  dplyr::filter(cell_type_a == "CD8_Tcell") %>%
  dplyr::left_join(unique(xen@meta.data[, c("sample_id", "timepoint")]),
                   by = "sample_id") %>%
  dplyr::group_by(timepoint, cell_type_b) %>%
  dplyr::summarise(z_mean = mean(z), z_sd = sd(z), .groups = "drop") %>%
  dplyr::arrange(desc(z_mean))
write.csv(cd8_partners, file.path(OUT_DIR, "cd8_neighborhood_partners.csv"),
          row.names = FALSE)
p_cd8 <- ggplot(cd8_partners,
                aes(x = reorder(cell_type_b, z_mean), y = z_mean, fill = timepoint)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = z_mean - z_sd, ymax = z_mean + z_sd),
                position = position_dodge(width = 0.8), width = 0.3) +
  coord_flip() + geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Neighbor cell type", y = "Enrichment z-score (vs CD8 T cell)") +
  theme_minimal(base_size = 11)
save_plot(p_cd8, "cd8_neighborhood_partners.png", 8, 6)

saveRDS(xen, file.path(OUT_DIR, "06_xenium_neighborhood.rds"))


# 07. MACROPHAGE ANALYSIS ------------------------------------------------------
# Resume: xen <- readRDS(file.path(OUT_DIR, "06_xenium_neighborhood.rds"))

mac <- subset(xen, subset = cell_type == "Macrophage")
cat("Macrophages: ", ncol(mac), "\n")

DefaultAssay(mac) <- "Xenium"
mac <- ScaleData(mac, assay = "Xenium", features = rownames(mac), verbose = FALSE)
mac <- RunPCA(mac, assay = "Xenium", features = rownames(mac), npcs = 20, verbose = FALSE)
mac <- RunUMAP(mac, dims = 1:20, verbose = FALSE)
mac <- FindNeighbors(mac, dims = 1:20, verbose = FALSE)
mac <- FindClusters(mac, resolution = 0.4, algorithm = 4, verbose = FALSE)

# Macrophage subset signatures, pre-filtered to the 480-gene panel.
# Dropped entirely (no markers on panel): Phagocytic, IFN_responsive.
# Sparse markers removed: Mrc1, Cd163, Folr2, Cbr2, Timd4, Ccr2, Ifitm6, Nos2,
#   Tnf, Stat1, Chil3, Trem2, Lpl, Gpnmb, Fabp5, Axl, Tyro3, Gas6, Isg15,
#   Ifit1, Ifit3, Irf7, Mx1.
mac_sigs <- list(
  TissueResident = c("Lyve1", "Mertk", "Mgl2"),
  Infiltrating   = c("Ly6c2", "Plac8"),
  M1_like        = c("Cd86", "Il1b", "Cxcl9", "Cxcl10"),
  M2_like        = c("Arg1", "Retnla", "Mgl2"),
  Lipid_assoc    = c("Spp1", "Apoe")
)
.report_missing(mac_sigs, rownames(mac), "mac_sigs")
mac_sigs <- .filter_panels(mac_sigs, rownames(mac))
for (n in names(mac_sigs)) {
  feats <- intersect(mac_sigs[[n]], rownames(mac))
  if (length(feats) > 0) {
    mac <- AddModuleScore(mac, features = list(feats),
                          name = paste0("mac_", n), assay = "Xenium",
                          nbin = 10, ctrl = 10)
  }
}

# Use the columns we just added (skipping mac_subset which is character and
# would have been carried over from a previous pipeline run).
mac_score_cols <- paste0("mac_", names(mac_sigs), "1")
mac_score_cols <- intersect(mac_score_cols, colnames(mac@meta.data))
mac_score_cols <- mac_score_cols[
  vapply(mac@meta.data[, mac_score_cols, drop = FALSE], is.numeric, logical(1))
]
p_mscores <- FeaturePlot(mac, features = mac_score_cols,
                         ncol = 3, raster = FALSE) & scale_color_viridis_c()
save_plot(p_mscores, "macrophage_signature_scores.png", 12, 8)

mac_markers <- FindAllMarkers(mac, only.pos = TRUE, min.pct = 0.2,
                              logfc.threshold = 0.4)
write.csv(mac_markers, file.path(OUT_DIR, "macrophage_subcluster_markers.csv"),
          row.names = FALSE)

# Auto-annotate macrophage subclusters by macrophage signature score.
mac_annotation <- auto_annotate_clusters(
  mac@meta.data, cluster_col = "seurat_clusters",
  score_prefix = "mac_", strip_prefix = "mac_",
  csv_name     = "macrophage_signature_scores_by_cluster.csv",
  heatmap_name = "macrophage_signature_heatmap.png"
)

# >>> OPTIONAL: override individual macrophage subset labels here.
# >>> Examples:
#   mac_annotation["1"] <- "TissueResident"
#   mac_annotation["3"] <- "M1_like"

mac$mac_subset <- ifelse(
  as.character(mac$seurat_clusters) %in% names(mac_annotation),
  mac_annotation[as.character(mac$seurat_clusters)],
  paste0("mac_cluster_", mac$seurat_clusters)
)

xen$mac_subset <- NA_character_
xen$mac_subset[colnames(xen) %in% colnames(mac)] <- mac$mac_subset[match(
  colnames(xen)[colnames(xen) %in% colnames(mac)], colnames(mac))]

df_mac <- xen@meta.data %>%
  dplyr::filter(!is.na(mac_subset), !is.na(niche_label)) %>%
  dplyr::count(timepoint, niche_label, mac_subset) %>%
  dplyr::group_by(timepoint, niche_label) %>%
  dplyr::mutate(frac = n / sum(n))
p_mn <- ggplot(df_mac, aes(x = niche_label, y = frac, fill = mac_subset)) +
  geom_col(position = "stack") + facet_wrap(~ timepoint) +
  labs(x = "Niche", y = "Fraction of macrophages",
       fill = "Macrophage subset") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p_mn, "macrophage_by_niche_time.png", 11, 5)

for (s in unique(xen$sample_id)) {
  p_sp <- ImageDimPlot(
    xen, fov = s, group.by = "mac_subset",
    cells = WhichCells(xen, expression = !is.na(mac_subset)),
    size = 0.7, dark.background = FALSE
  ) + ggtitle(paste0("Macrophage subsets — ", s))
  save_plot(p_sp, paste0("spatial_macrophage_", s, ".png"), 8, 8)
}

# Macrophage <-> CD8 spatial coupling
mac_subsets_vec <- na.omit(unique(xen$mac_subset))
records <- list()
for (s in unique(xen$sample_id)) {
  sub <- subset(xen, subset = sample_id == s)
  fov <- sub@images[[s]]
  if (is.null(fov)) next
  coords <- GetTissueCoordinates(fov)
  rownames(coords) <- coords$cell
  tc_cells <- colnames(sub)[sub$cell_type == "CD8_Tcell"]
  if (length(tc_cells) == 0) next
  tc_xy <- coords[tc_cells, c("x", "y"), drop = FALSE]
  for (m in mac_subsets_vec) {
    src_cells <- colnames(sub)[which(sub$mac_subset == m)]
    if (length(src_cells) < 5) next
    src_xy <- coords[src_cells, c("x", "y"), drop = FALSE]
    nn <- get.knnx(as.matrix(src_xy), as.matrix(tc_xy), k = 1)
    records[[length(records) + 1]] <- data.frame(
      cell         = tc_cells,
      sample_id    = s,
      timepoint    = sub$timepoint[1],
      tcell_subset = sub$tcell_subset[match(tc_cells, colnames(sub))],
      mac_subset   = m,
      distance_um  = nn$nn.dist[, 1]
    )
  }
}
mac_dist <- bind_rows(records)
write.csv(mac_dist, file.path(OUT_DIR, "tcell_to_macrophage_distances.csv"),
          row.names = FALSE)

p_md <- ggplot(mac_dist,
               aes(x = mac_subset, y = distance_um, fill = tcell_subset)) +
  geom_boxplot(outlier.size = 0.3) + facet_wrap(~ timepoint) +
  scale_y_log10() +
  labs(x = "Macrophage subset",
       y = "CD8+ T cell distance to nearest mac (μm, log10)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p_md, "tcell_macrophage_distances.png", 11, 6)

saveRDS(mac, file.path(OUT_DIR, "07_macrophage_subset.rds"))
saveRDS(xen, file.path(OUT_DIR, "07_xenium_full.rds"))


# 08. MACROPHAGE DEEP DIVE -----------------------------------------------------
# Goals:
#   8a. Heatmap of genes differentially expressed between macrophage subsets
#   8b. GO term enrichment per subset, using the 480-gene panel as the background
#       universe (critical for targeted panels — avoids false significance from
#       a non-genome-wide gene set).
#   8c. M2-like deep dive: M2_like vs all other macrophages
#   8d. Spatial proximity to epithelium: classify each macrophage as
#       epithelium-proximal vs epithelium-distal based on distance.
#   8e. Pseudobulk DE (DESeq2): proximal vs distal macrophages, both globally
#       and within the M2_like subset.
#   8f. GO enrichment of proximity-DE genes.
#
# Resume: 
mac <- readRDS(file.path(OUT_DIR, "07_macrophage_subset.rds"))
xen <- readRDS(file.path(OUT_DIR, "07_xenium_full.rds"))

suppressPackageStartupMessages({
  library(DESeq2)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(pheatmap)
  library(ggrepel)
})


# 8a. DE genes between macrophage subsets --------------------------------------

Idents(mac) <- "mac_subset"
mac_subset_markers <- FindAllMarkers(
  mac, only.pos = TRUE,
  min.pct = 0.10, logfc.threshold = 0.25
)
write.csv(mac_subset_markers,
          file.path(OUT_DIR, "macrophage_subset_DE_markers.csv"),
          row.names = FALSE)

# Top 10 markers per subset for the heatmap
top_per_subset <- mac_subset_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(avg_log2FC, n = 10, with_ties = FALSE) %>%
  dplyr::ungroup()
top_genes <- unique(top_per_subset$gene)

# Average normalized expression per subset for a clean heatmap
avg_expr <- AverageExpression(mac, assays = "Xenium",
                              group.by = "mac_subset",
                              layer = "data")[["Xenium"]]
avg_expr <- avg_expr[rownames(avg_expr) %in% top_genes, , drop = FALSE]
# Scale rows (z-score per gene)
mat_z <- t(scale(t(as.matrix(avg_expr))))
mat_z[!is.finite(mat_z)] <- 0

p_mac_heat <- pheatmap::pheatmap(
  mat_z,
  cluster_rows = TRUE, cluster_cols = TRUE,
  fontsize_row = 7, fontsize_col = 10,
  color = colorRampPalette(c("navy", "white", "firebrick"))(100),
  main = "Top markers per macrophage subset (z-scored mean expression)",
  filename = file.path(OUT_DIR, "macrophage_subset_DE_heatmap.png"),
  width = 7, height = max(8, 0.18 * nrow(mat_z))
)
# pheatmap's filename argument auto-saves; also push to active device:
print(p_mac_heat)


# 8b. GO enrichment per macrophage subset --------------------------------------
# Use the 480-gene panel as the background universe so we don't get spurious
# enrichment from genes that simply weren't on the panel.

panel_genes <- rownames(mac)  # 480 panel genes that survived QC
panel_entrez <- suppressMessages(
  bitr(panel_genes, fromType = "SYMBOL",
       toType = "ENTREZID", OrgDb = org.Mm.eg.db)
) %>%
  dplyr::filter(!is.na(ENTREZID), nzchar(ENTREZID)) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)
universe_entrez <- unique(as.character(panel_entrez$ENTREZID))
cat("Panel universe size for GO: ", length(universe_entrez), " Entrez IDs\n")

# Wrapper that catches every flavor of clusterProfiler error and returns NULL.
# Notes:
#  - Use suppressWarnings() rather than a tryCatch warning handler. A warning
#    handler in tryCatch re-invokes the call and the error then escapes.
#  - With a 480-gene universe, q-value estimation can fail (empirical p-value
#    distribution too small). qvalueCutoff = 1 bypasses the q-value filter.
#  - Drop rows with NA qvalue from the @result table — clusterProfiler's
#    dotplot() crashes internally on NA q-values.
safe_enrichGO <- function(genes_entrez, universe_entrez, ont = "BP",
                          p_cut = 0.05) {
  if (length(genes_entrez) < 3) return(NULL)
  res <- tryCatch(
    suppressWarnings(enrichGO(
      gene          = genes_entrez,
      universe      = universe_entrez,
      OrgDb         = org.Mm.eg.db,
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = p_cut,
      qvalueCutoff  = 1,
      readable      = TRUE
    )),
    error = function(e) {
      message("  enrichGO error: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(res) || is.null(res@result) || nrow(res@result) == 0) {
    return(NULL)
  }
  # Strip rows that would crash dotplot()
  res@result <- res@result[!is.na(res@result$qvalue) &
                           !is.na(res@result$p.adjust), , drop = FALSE]
  if (nrow(res@result) == 0) return(NULL)
  res
}

# Wrap dotplot too — even with cleaned q-values it can still throw on edge cases.
safe_dotplot <- function(ego, title = "", showCategory = 15) {
  tryCatch(
    suppressWarnings(dotplot(ego, showCategory = showCategory)) +
      ggtitle(title),
    error = function(e) {
      message("  dotplot error (skipping plot): ", conditionMessage(e))
      NULL
    }
  )
}

go_results_per_subset <- list()
for (sub_name in unique(mac_subset_markers$cluster)) {
  genes <- mac_subset_markers %>%
    dplyr::filter(cluster == sub_name, p_val_adj < 0.05) %>%
    dplyr::pull(gene) %>% unique()
  if (length(genes) < 3) {
    message("Skipping ", sub_name, " — too few DE genes (", length(genes), ").")
    next
  }
  ent <- panel_entrez$ENTREZID[match(genes, panel_entrez$SYMBOL)]
  ent <- unique(as.character(ent[!is.na(ent) & nzchar(ent)]))
  cat(sub_name, ": ", length(ent), " mappable Entrez IDs out of ",
      length(genes), " DE genes\n", sep = "")

  ego <- safe_enrichGO(ent, universe_entrez, ont = "BP")
  go_results_per_subset[[sub_name]] <- ego

  if (!is.null(ego)) {
    write.csv(ego@result,
              file.path(OUT_DIR, paste0("GO_BP_", sub_name, ".csv")),
              row.names = FALSE)
    p_go <- safe_dotplot(ego, title = paste0("GO BP — ", sub_name))
    if (!is.null(p_go)) {
      save_plot(p_go, paste0("GO_BP_", sub_name, ".png"), 9, 6)
    }
  } else {
    message("  No significant GO BP terms for ", sub_name,
            " against the 480-gene panel universe.")
  }
}


# 8c. M2-like deep dive: M2_like vs all other macrophages ----------------------

m2_label <- grep("M2", unique(mac$mac_subset), value = TRUE)[1]
if (is.na(m2_label)) {
  warning("No mac_subset label containing 'M2' was found — skipping M2 deep dive.")
} else {
  message("Using '", m2_label, "' as the M2 label.")
  mac$is_M2 <- mac$mac_subset == m2_label

  Idents(mac) <- "is_M2"
  m2_vs_rest <- FindMarkers(mac, ident.1 = TRUE, ident.2 = FALSE,
                            logfc.threshold = 0.1, min.pct = 0.05) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::arrange(p_val_adj)
  write.csv(m2_vs_rest,
            file.path(OUT_DIR, "M2_vs_other_macrophages_DE.csv"),
            row.names = FALSE)

  # Volcano plot
  m2_vs_rest$log10p <- -log10(pmax(m2_vs_rest$p_val_adj, 1e-300))
  m2_vs_rest$signif <- m2_vs_rest$p_val_adj < 0.05 &
                       abs(m2_vs_rest$avg_log2FC) > 0.25
  top_label <- m2_vs_rest %>%
    dplyr::arrange(desc(log10p * sign(avg_log2FC))) %>%
    dplyr::slice_head(n = 20)
  p_volc <- ggplot(m2_vs_rest,
                   aes(x = avg_log2FC, y = log10p, color = signif)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = c("grey60", "firebrick")) +
    ggrepel::geom_text_repel(data = top_label, aes(label = gene),
                             size = 3, max.overlaps = 30) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed") +
    labs(x = "log2 FC (M2 vs other mac)", y = "-log10 adjusted p",
         title = "M2-like vs other macrophages") +
    theme_minimal(base_size = 11) + theme(legend.position = "none")
  save_plot(p_volc, "M2_vs_other_volcano.png", 8, 6)

  # GO of M2-up genes
  m2_up <- m2_vs_rest %>%
    dplyr::filter(p_val_adj < 0.05, avg_log2FC > 0.25) %>%
    dplyr::pull(gene)
  ent_m2 <- panel_entrez$ENTREZID[match(m2_up, panel_entrez$SYMBOL)]
  ent_m2 <- unique(as.character(ent_m2[!is.na(ent_m2) & nzchar(ent_m2)]))
  ego_m2 <- safe_enrichGO(ent_m2, universe_entrez, ont = "BP")
  if (!is.null(ego_m2)) {
    write.csv(ego_m2@result,
              file.path(OUT_DIR, "GO_BP_M2_up.csv"), row.names = FALSE)
    p_go_m2 <- safe_dotplot(ego_m2, title = "GO BP — M2-up genes",
                            showCategory = 20)
    if (!is.null(p_go_m2)) save_plot(p_go_m2, "GO_BP_M2_up.png", 9, 7)
  } else {
    message("  No significant GO BP terms for M2-up genes.")
  }
}


# 8d. Distance to nearest epithelial cell --------------------------------------
# Compute per-macrophage distance to the closest epithelial-class cell
# within the SAME sample (FOV), then store on the xen object.

# Treat any 'epithelial' label as epithelium (catches Luminal_epithelial,
# Basal_epithelial, Stem_epithelial, etc.)
epi_pattern <- "epithelial|Epithelial"
xen$is_epithelial <- grepl(epi_pattern, xen$cell_type)
cat("Epithelial cells: ", sum(xen$is_epithelial), "\n")

dist_to_epi <- rep(NA_real_, ncol(xen))
names(dist_to_epi) <- colnames(xen)

for (s in unique(xen$sample_id)) {
  sub <- subset(xen, subset = sample_id == s)
  fov <- sub@images[[s]]
  if (is.null(fov)) next
  coords <- GetTissueCoordinates(fov)
  rownames(coords) <- coords$cell

  # All cells in this sample with valid coords
  all_cells <- intersect(colnames(sub), rownames(coords))
  xy_all   <- as.matrix(coords[all_cells, c("x", "y"), drop = FALSE])

  epi_cells <- all_cells[sub$is_epithelial[all_cells]]
  if (length(epi_cells) < 5) {
    message("Skipping ", s, " — fewer than 5 epithelial cells.")
    next
  }
  xy_epi <- as.matrix(coords[epi_cells, c("x", "y"), drop = FALSE])

  nn <- get.knnx(xy_epi, xy_all, k = 1)
  dist_to_epi[all_cells] <- nn$nn.dist[, 1]
}
xen$dist_to_epithelium_um <- dist_to_epi

# Mirror this onto the macrophage object
mac$dist_to_epithelium_um <- xen$dist_to_epithelium_um[colnames(mac)]

# Histogram of macrophage distances
p_dhist <- ggplot(data.frame(d = mac$dist_to_epithelium_um,
                              subset = mac$mac_subset),
                  aes(x = d, fill = subset)) +
  geom_histogram(bins = 60, alpha = 0.6, position = "identity") +
  scale_x_log10() +
  geom_vline(xintercept = 50, linetype = "dashed") +
  labs(x = "Distance from nearest epithelial cell (μm, log10)",
       y = "Number of macrophages",
       title = "Macrophage distance from epithelium") +
  theme_minimal(base_size = 11)
save_plot(p_dhist, "macrophage_distance_histogram.png", 9, 5)


# 8e. Pseudobulk DE: epithelium-proximal vs epithelium-distal ------------------
# Cutoff = 50 μm (roughly 2 epithelial-cell diameters). Adjust to taste.
PROX_CUTOFF_UM <- 50

mac$proximity <- ifelse(
  is.na(mac$dist_to_epithelium_um), NA_character_,
  ifelse(mac$dist_to_epithelium_um <= PROX_CUTOFF_UM, "proximal", "distal")
)

# Composition table — how do the subsets distribute across proximal/distal?
comp_tbl <- mac@meta.data %>%
  dplyr::filter(!is.na(proximity)) %>%
  dplyr::count(mac_subset, proximity) %>%
  dplyr::group_by(mac_subset) %>%
  dplyr::mutate(frac = n / sum(n))
write.csv(comp_tbl,
          file.path(OUT_DIR, "macrophage_proximity_composition.csv"),
          row.names = FALSE)
p_comp_prox <- ggplot(comp_tbl,
                      aes(x = mac_subset, y = frac, fill = proximity)) +
  geom_col() + labs(x = "Macrophage subset", y = "Fraction") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p_comp_prox, "macrophage_proximity_composition.png", 7, 5)

# Helper: pseudobulk a Seurat object by (sample_id, proximity) and run DESeq2
run_pseudobulk_DE <- function(seu, label = "all_mac") {
  obj <- subset(seu, cells = colnames(seu)[!is.na(seu$proximity)])
  if (ncol(obj) < 50) {
    message("Too few cells for ", label, " (", ncol(obj), ") — skipping.")
    return(NULL)
  }
  obj$pb_group <- paste(obj$sample_id, obj$proximity, sep = "__")

  # Aggregate raw counts across cells in each pseudobulk group
  pb <- AggregateExpression(obj, assays = "Xenium", group.by = "pb_group",
                            return.seurat = FALSE)[["Xenium"]]
  pb <- as.matrix(pb)
  pb <- round(pb)  # ensure integers for DESeq2

  # Build coldata
  coldata <- data.frame(
    pb_group  = colnames(pb),
    sample_id = sub("__.*$", "",  colnames(pb)),
    proximity = sub("^.*__", "",  colnames(pb)),
    row.names = colnames(pb),
    stringsAsFactors = FALSE
  )
  coldata$proximity <- factor(coldata$proximity, levels = c("distal", "proximal"))
  coldata$sample_id <- factor(coldata$sample_id)

  # Need both proximity levels and >= 2 samples per level
  ok <- table(coldata$proximity)
  if (any(ok < 2)) {
    message("Not enough pseudobulk reps for ", label, " (",
            paste(names(ok), ok, sep = "=", collapse = ", "), ") — skipping.")
    return(NULL)
  }

  # Drop genes that are zero across all groups
  keep_g <- rowSums(pb) > 0
  pb <- pb[keep_g, , drop = FALSE]

  dds <- DESeqDataSetFromMatrix(countData = pb, colData = coldata,
                                design = ~ sample_id + proximity)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("proximity", "proximal", "distal"))
  df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::arrange(padj)
  write.csv(df, file.path(OUT_DIR,
            paste0("pseudobulk_DE_proximal_vs_distal_", label, ".csv")),
            row.names = FALSE)

  # Volcano
  df$log10p <- -log10(pmax(df$padj, 1e-300))
  df$signif <- !is.na(df$padj) & df$padj < 0.05 & abs(df$log2FoldChange) > 0.5
  top_lab <- df %>% dplyr::filter(!is.na(padj)) %>%
    dplyr::arrange(desc(log10p * sign(log2FoldChange))) %>%
    dplyr::slice_head(n = 20)
  p_v <- ggplot(df, aes(x = log2FoldChange, y = log10p, color = signif)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = c("grey60", "firebrick")) +
    ggrepel::geom_text_repel(data = top_lab, aes(label = gene),
                             size = 3, max.overlaps = 30) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
    labs(x = "log2 FC (proximal vs distal)",
         y = "-log10 adjusted p",
         title = paste0("Pseudobulk DE — ", label,
                        " (proximal vs distal epithelium)")) +
    theme_minimal(base_size = 11) + theme(legend.position = "none")
  save_plot(p_v, paste0("pseudobulk_volcano_", label, ".png"), 9, 6)

  # GO enrichment of proximity-up and proximity-down genes
  for (direction in c("up", "down")) {
    sig <- if (direction == "up") {
      df %>% dplyr::filter(!is.na(padj), padj < 0.05, log2FoldChange > 0.5)
    } else {
      df %>% dplyr::filter(!is.na(padj), padj < 0.05, log2FoldChange < -0.5)
    }
    if (nrow(sig) < 3) next
    ent <- panel_entrez$ENTREZID[match(sig$gene, panel_entrez$SYMBOL)]
    ent <- unique(as.character(ent[!is.na(ent) & nzchar(ent)]))
    ego <- safe_enrichGO(ent, universe_entrez, ont = "BP")
    if (!is.null(ego)) {
      write.csv(ego@result, file.path(OUT_DIR,
                paste0("GO_BP_", label, "_proximal_", direction, ".csv")),
                row.names = FALSE)
      p_g <- safe_dotplot(ego,
                          title = paste0("GO BP — ", label,
                                         " (proximal-", direction, ")"))
      if (!is.null(p_g)) {
        save_plot(p_g, paste0("GO_BP_", label, "_proximal_",
                              direction, ".png"), 9, 6)
      }
    } else {
      message("  No significant GO BP terms for ", label, "/proximal-", direction, ".")
    }
  }
  df
}

# Run for all macrophages, then specifically for M2-like
de_all <- run_pseudobulk_DE(mac, label = "all_mac")
if (exists("m2_label") && !is.na(m2_label)) {
  mac_m2 <- subset(mac, subset = mac_subset == m2_label)
  de_m2  <- run_pseudobulk_DE(mac_m2, label = "M2_like")
}

saveRDS(mac, file.path(OUT_DIR, "08_macrophage_with_proximity.rds"))
saveRDS(xen, file.path(OUT_DIR, "08_xenium_with_dist_to_epi.rds"))

message("Pipeline complete. Outputs in: ", OUT_DIR)
