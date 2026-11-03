# 02_qc_clustering.R
# QC, normalization, dimensionality reduction, clustering, and marker-based
# cell-type annotation for the Xenium prostate dataset.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(harmony)
  library(future)
})

# For a 541-gene panel x ~550k cells, parallel data shipping by future blows
# past the default 500MB globals budget. Run sequentially and lift the cap.
plan("sequential")
options(future.globals.maxSize = 16 * 1024^3)  # 16 GB

OUT_DIR <- "/Users/alekhyamazumdar/Documents/Xenium/Output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
xen <- readRDS(file.path(OUT_DIR, "01_xenium_raw.rds"))
DefaultAssay(xen) <- "Xenium"

# Seurat v5 keeps a separate counts layer per merged sample (counts.1, counts.2, ...).
# Join them into a single 'counts' layer so downstream functions that expect one
# layer (FindAllMarkers, GetAssayData, etc.) work without error.
xen <- JoinLayers(xen, assay = "Xenium")

# ---------- 1. QC ----------
# nCount_Xenium and nFeature_Xenium are populated automatically by
# CreateSeuratObject() at load time, so we don't need to recompute them.
# (If for some reason they're missing, uncomment the two lines below.)
# xen$nCount_Xenium   <- Matrix::colSums(LayerData(xen, layer = "counts"))
# xen$nFeature_Xenium <- Matrix::colSums(LayerData(xen, layer = "counts") > 0)

# Per-sample distribution of transcripts per cell
p_qc <- VlnPlot(
  xen,
  features = c("nCount_Xenium", "nFeature_Xenium"),
  group.by = "sample_id",
  pt.size = 0
) + plot_layout(ncol = 1)
ggsave(file.path(OUT_DIR, "qc_violins.png"), p_qc, width = 9, height = 6)

# Filter very low-count cells. The Xenium pipeline already removes empties,
# but cells with <10 transcripts are usually segmentation artifacts.
xen <- subset(xen, subset = nCount_Xenium >= 10 & nFeature_Xenium >= 5)
cat("Cells after QC: ", ncol(xen), "\n")

# ---------- 2. Normalize ----------
# For a targeted 541-gene panel SCTransform is overkill (and memory-heavy).
# Standard log-normalization + scaling is appropriate and far cheaper.
# We don't subset to "variable features" because all 541 genes were chosen
# deliberately for the panel — use them all.
xen <- NormalizeData(xen, assay = "Xenium", verbose = FALSE)
xen <- ScaleData(xen, assay = "Xenium", features = rownames(xen),
                 verbose = FALSE)

# ---------- 3. PCA + Harmony integration across samples ----------
# RunPCA on all 541 features (no VariableFeatures step needed for a panel).
xen <- RunPCA(xen, assay = "Xenium", features = rownames(xen),
              npcs = 30, verbose = FALSE)

# Integrate so clustering is driven by biology, not by replicate/batch
xen <- RunHarmony(xen, group.by.vars = "sample_id", reduction = "pca",
                  reduction.save = "harmony", assay.use = "Xenium", verbose = FALSE)

# ---------- 4. UMAP + Leiden clustering ----------
xen <- RunUMAP(xen, reduction = "harmony", dims = 1:30, verbose = FALSE)
xen <- FindNeighbors(xen, reduction = "harmony", dims = 1:30, verbose = FALSE)
xen <- FindClusters(xen, resolution = 0.5, algorithm = 4, verbose = FALSE)  # 4 = Leiden

p_umap <- DimPlot(xen, group.by = "seurat_clusters", label = TRUE, raster = FALSE) +
  DimPlot(xen, group.by = "sample_id", raster = FALSE)
ggsave(file.path(OUT_DIR, "umap_clusters.png"), p_umap, width = 12, height = 5)

# ---------- 5. Marker-based cell type annotation ----------
# These markers are present on Takehara's 480-gene Xenium panel. Adjust if your
# panel differs. Genes with no expression in your data are silently skipped.
marker_panels <- list(
  Epithelial         = c("Krt8", "Krt5", "Krt18", "Cdh1", "Epcam", "Krt14"),
  Luminal_epithelial = c("Krt8", "Krt18", "Nkx3-1", "Pbsn", "Ar"),
  Basal_epithelial   = c("Krt5", "Krt14", "Trp63"),
  Smooth_muscle      = c("Acta2", "Myh11", "Tagln", "Cnn1"),
  Fibroblast         = c("Col1a1", "Col1a2", "Pdgfra", "Dcn", "Lum"),
  Endothelial        = c("Pecam1", "Cdh5", "Cd34", "Vwf"),
  CD8_Tcell          = c("Cd8a", "Cd8b1", "Cd3e", "Cd3d"),
  CD4_Tcell          = c("Cd4", "Cd3e", "Foxp3"),
  B_cell             = c("Cd19", "Ms4a1", "Cd79a", "Cd79b"),
  Macrophage         = c("Cd68", "Adgre1", "C1qa", "C1qb", "Csf1r", "Lyz2"),
  Dendritic          = c("Itgax", "Xcr1", "Clec9a", "H2-Aa"),
  NK_cell            = c("Ncr1", "Klrb1c", "Gzma"),
  Neutrophil         = c("S100a8", "S100a9", "Ly6g", "Mpo"),
  Trm_signature      = c("Cd69", "Itgae", "Runx3", "Prdm1", "Hobit", "Notch1")
)

# Score each cluster against the marker sets
DefaultAssay(xen) <- "Xenium"
for (n in names(marker_panels)) {
  feats <- intersect(marker_panels[[n]], rownames(xen))
  if (length(feats) >= 1) {
    xen <- AddModuleScore(xen, features = list(feats), name = paste0("score_", n),
                          assay = "Xenium")
  }
}

# DotPlot of canonical markers per cluster
all_markers <- unique(unlist(marker_panels))
all_markers <- intersect(all_markers, rownames(xen))
p_dot <- DotPlot(xen, features = all_markers, group.by = "seurat_clusters") +
  RotatedAxis()
ggsave(file.path(OUT_DIR, "marker_dotplot.png"), p_dot, width = 14, height = 6)

# Fast top-marker discovery per cluster (presto-accelerated Wilcoxon)
cluster_markers <- FindAllMarkers(
  xen, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5
)
write.csv(cluster_markers, file.path(OUT_DIR, "cluster_markers.csv"), row.names = FALSE)

# ---------- 6. Manual cell-type annotation ----------
# Inspect cluster_markers.csv and the dotplot, then fill this lookup.
# The Idents are the cluster IDs as character.
cluster_annotation <- c(
  "1"  = "Luminal_epithelial",
  "2"  = "Basal_epithelial",
  "3"  = "Smooth_muscle",
  "4"  = "Fibroblast",
  "5"  = "Endothelial",
  "6"  = "CD8_Tcell",
  "7"  = "CD4_Tcell",
  "8"  = "Macrophage",
  "9"  = "Dendritic",
  "10" = "B_cell",
  "11" = "NK_cell",
  "12" = "Neutrophil"
  # Extend / adjust based on your clustering
)
xen$cell_type <- ifelse(
  as.character(xen$seurat_clusters) %in% names(cluster_annotation),
  cluster_annotation[as.character(xen$seurat_clusters)],
  paste0("cluster_", xen$seurat_clusters)
)

# Visualize
p_celltype <- DimPlot(xen, group.by = "cell_type", label = TRUE, raster = FALSE)
ggsave(file.path(OUT_DIR, "umap_celltype.png"), p_celltype, width = 9, height = 7)

# Spatial plot per sample
for (s in unique(xen$sample_id)) {
  p_sp <- ImageDimPlot(xen, fov = s, group.by = "cell_type", size = 0.5, dark.background = FALSE)
  ggsave(file.path(OUT_DIR, paste0("spatial_celltype_", s, ".png")),
         p_sp, width = 8, height = 8)
}

saveRDS(xen, file.path(OUT_DIR, "02_xenium_annotated.rds"))
message("Saved: ", file.path(OUT_DIR, "02_xenium_annotated.rds"))
