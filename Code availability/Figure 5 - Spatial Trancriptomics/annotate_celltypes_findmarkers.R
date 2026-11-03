# =============================================================================
# annotate_celltypes_findmarkers.R
# Cell-type annotation of the 10x Xenium Prime Human Prostate dataset using
# DATA-DRIVEN cluster markers (FindAllMarkers) instead of module-score
# auto-annotation.
#
# This is a FRESH run from the converted MTX — it does NOT read any .rds.
#
# Workflow:
#   1. Load converted Xenium data (MTX + cells.csv)  — from scratch
#   2. QC filtering
#   3. Normalize (CLR), PCA, cluster, UMAP
#   4. FindAllMarkers on the unsupervised clusters  ← the key step
#   5. Marker tables + heatmap + dot plots (curated panel shown only as a guide)
#   6. Manual cluster -> cell-type mapping (YOU fill it in), then apply,
#      run a per-cell-type FindMarkers, and save the annotated object
#
# PRE-REQUISITE: 00_convert_xenium_zarr.py must already have produced:
#     ~/Documents/Xenium/Xenium_10x_converted/cell_feature_matrix/
#     ~/Documents/Xenium/Xenium_10x_converted/cells.csv.gz
#
# Required packages: Seurat (>=5.0), ggplot2, dplyr, patchwork, viridis, Matrix
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(viridis)
library(Matrix)

set.seed(42)   # reproducible clustering / UMAP

# ── Paths ─────────────────────────────────────────────────────────────────────
converted_path <- "/Users/alekhyamazumdar/Documents/Xenium/Xenium_10x_converted"
mtx_path       <- file.path(converted_path, "cell_feature_matrix")
cells_csv      <- file.path(converted_path, "cells.csv.gz")
output_dir     <- "/Users/alekhyamazumdar/Documents/Xenium/Output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

if (!dir.exists(mtx_path)) {
  stop(paste0(
    "\nConverted data not found at: ", mtx_path,
    "\nPlease run  00_convert_xenium_zarr.py  first."
  ))
}

# ── Analysis parameters (change here, in one place) ───────────────────────────
CLUSTER_RES   <- 0.4    # clustering resolution (match script 08 by default)
PCS_USE       <- 20     # PCs for neighbours / UMAP
N_TOP_MARKERS <- 10     # top markers per cluster to display

# ── Plotting theme ─────────────────────────────────────────────────────────────
theme_pub <- theme_classic(base_size = 12) +
  theme(
    plot.title   = element_text(face = "bold", size = 13),
    axis.text    = element_text(color = "black"),
    legend.title = element_text(face = "bold")
  )

# ── Curated prostate panel — REFERENCE ONLY ───────────────────────────────────
# NOTE: this is used purely as a visual cross-check (a reference dot plot). The
# actual cell-type calls come from FindAllMarkers below + your manual mapping.
cell_type_markers <- list(
  Epithelial_Luminal  = c("KRT8", "KRT18", "AR", "FOXA1", "NKX3-1"),
  Epithelial_Basal    = c("KRT5", "KRT14", "TP63", "ACTA2"),
  Epithelial_Club     = c("SCGB1A1", "PIGR", "KRT7"),
  Fibroblast          = c("VIM", "DCN", "COL1A1", "COL3A1", "PDGFRA"),
  Smooth_Muscle       = c("ACTA2", "TAGLN", "MYH11", "CNN1"),
  Endothelial         = c("PECAM1", "CDH5", "VWF", "CLDN5"),
  Macrophage          = c("CD68", "CD163", "MRC1", "CSF1R", "ITGAM"),
  Mast_Cell           = c("KIT", "TPSAB1", "CPA3"),
  T_Cell              = c("CD3E", "CD3D", "CD8A", "CD4", "FOXP3"),
  B_NK_Cell           = c("MS4A1", "CD79A", "NKG7", "GNLY"),
  Dendritic_Cell      = c("FCER1A", "CD1C", "CLEC9A", "LILRA4"),
  Plasma_Cell         = c("IGKC", "IGLC2", "SDC1", "JCHAIN")
)


# =============================================================================
# SECTION 1 — LOAD XENIUM DATA (from converted MTX)
# =============================================================================
message("[1] Loading Xenium data from converted MTX ...")

counts <- Read10X(data.dir = mtx_path)
xenium <- CreateSeuratObject(
  counts  = counts,
  assay   = "Xenium",
  project = "XeniumProstatePrime5K"
)

# Spatial coordinates
cells_meta   <- read.csv(cells_csv)
common_cells <- intersect(colnames(xenium), cells_meta$cell_id)
cells_meta   <- cells_meta[match(common_cells, cells_meta$cell_id), ]
xenium       <- xenium[, common_cells]

xenium <- AddMetaData(xenium, metadata = cells_meta$x_centroid, col.name = "x_centroid")
xenium <- AddMetaData(xenium, metadata = cells_meta$y_centroid, col.name = "y_centroid")
xenium <- AddMetaData(xenium, metadata = cells_meta$cell_area,  col.name = "cell_area")

coords_df <- data.frame(row.names = common_cells,
                        x = cells_meta$x_centroid,
                        y = cells_meta$y_centroid)
xenium[["fov"]] <- CreateFOV(coords = coords_df, type = "centroids", assay = "Xenium")

cat(sprintf("Cells loaded  : %d\n", ncol(xenium)))
cat(sprintf("Genes in panel: %d\n", nrow(xenium)))


# =============================================================================
# SECTION 2 — QC FILTERING
# =============================================================================
message("[2] QC filtering ...")

p_qc <- VlnPlot(xenium, features = "nCount_Xenium", pt.size = 0) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
  scale_y_log10() + ggtitle("Transcript counts per cell") + theme_pub
ggsave(file.path(output_dir, "annot_qc_transcript_counts.png"), p_qc, width = 6, height = 5)

p_area <- ggplot(xenium@meta.data, aes(x = cell_area)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "white") +
  scale_x_log10() +
  geom_vline(xintercept = c(10, 5000), linetype = "dashed", color = "red") +
  labs(title = "Cell area distribution", x = "Cell area (µm²)", y = "Count") + theme_pub
ggsave(file.path(output_dir, "annot_qc_cell_area.png"), p_area, width = 6, height = 4)

xenium <- subset(xenium, subset = nCount_Xenium >= 10 & cell_area >= 10 & cell_area <= 5000)
cat(sprintf("Cells after QC: %d\n", ncol(xenium)))


# =============================================================================
# SECTION 3 — NORMALIZE, PCA, CLUSTER, UMAP
# =============================================================================
message("[3] Normalizing, clustering, UMAP ...")

xenium <- NormalizeData(xenium, normalization.method = "CLR", margin = 2)
xenium <- FindVariableFeatures(xenium, nfeatures = 2000)
xenium <- ScaleData(xenium)
xenium <- RunPCA(xenium, npcs = 30, verbose = FALSE)

p_elbow <- ElbowPlot(xenium, ndims = 30) + theme_pub
ggsave(file.path(output_dir, "annot_pca_elbow.png"), p_elbow, width = 6, height = 4)

xenium <- FindNeighbors(xenium, reduction = "pca", dims = 1:PCS_USE)
xenium <- FindClusters(xenium, resolution = CLUSTER_RES, verbose = FALSE)
xenium <- RunUMAP(xenium, dims = 1:PCS_USE, verbose = FALSE)

cat(sprintf("Clusters found (res %.2f): %d\n",
            CLUSTER_RES, nlevels(xenium$seurat_clusters)))

p_umap <- DimPlot(xenium, label = TRUE, label.size = 4) +
  ggtitle(sprintf("Unsupervised clusters (res %.2f)", CLUSTER_RES)) + theme_pub
ggsave(file.path(output_dir, "annot_umap_clusters.png"), p_umap, width = 7, height = 6)
print(p_umap)

# ── Spatial map: every cell at its tissue coordinate, coloured by cluster ──────
# (x/y swapped as in script 08's overlay so the tissue reads in the tall
#  orientation; drop the swap if you prefer the raw x = x_centroid layout.)
n_cl_sp    <- nlevels(xenium$seurat_clusters)
clust_cols <- setNames(DiscretePalette(n_cl_sp, palette = "polychrome"),
                       levels(xenium$seurat_clusters))

spatial_df <- data.frame(
  x       = xenium$y_centroid,          # swapped for tall orientation
  y       = xenium$x_centroid,
  cluster = xenium$seurat_clusters
)
p_spatial_clusters <- ggplot(spatial_df, aes(x = x, y = y, color = cluster)) +
  geom_point(size = 0.25, stroke = 0, alpha = 0.9) +
  scale_color_manual(values = clust_cols, name = "Cluster") +
  coord_fixed() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1)) +
  labs(title = sprintf("Spatial position of all clusters (res %.2f)", CLUSTER_RES)) +
  theme_pub +
  theme(axis.line  = element_blank(), axis.text  = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank())
# Printed to the plot pane only — not saved to disk.
# (To export later: ggsave(file.path(output_dir, "annot_spatial_all_clusters.png"),
#  p_spatial_clusters, width = 9, height = 11, dpi = 300, bg = "white"))
print(p_spatial_clusters)


# =============================================================================
# SECTION 3b — SUGGESTED CELL-TYPE ANNOTATION + annotated spatial map
# =============================================================================
# Cluster -> cell-type calls suggested from the FindAllMarkers table (Section 4)
# cross-checked against the curated panel. EDIT any label you disagree with.
# Confidence is noted in the comments. Clusters 7 & 12 had NO significant
# markers (likely low-quality / mixed) and are left "Unknown".
cluster_annotation <- c(
  "0"  = "Epithelial_Luminal",   # AMACR, KLK4, SIM2, GDF15 — tumour luminal   (high)
  "1"  = "Epithelial_Luminal",   # weak markers (SPON2, ANPEP) low-distinction  (low)
  "2"  = "Smooth_Muscle",        # SYNPO2, ADRA1A, KCNMB1, FLNC                 (high)
  "3"  = "Fibroblast",           # PDGFRA, PI16, COL5A2, ELN                    (high)
  "4"  = "Endothelial",          # CLDN5, TEK/TIE1, PLVAP, SOX17                (very high)
  "5"  = "Epithelial_Basal",     # TP63, COL17A1, ITGB4, LAMB3                  (high)
  "6"  = "Macrophage",           # CD163, CSF1R, MRC1, CYBB                     (very high)
  "7"  = "Unknown",              # no significant markers — likely low-quality
  "8"  = "T_NK_Cell",            # CD3E, CD8A, GZMK, KLRB1, CTLA4               (high)
  "9"  = "Epithelial_Other",     # MME, CYP3A5, CD177 — uncertain subtype       (low)
  "10" = "Mast_Cell",            # MS4A2, HDC, KIT, IL1RL1                      (very high)
  "11" = "Schwann_Neural",       # SOX10, MPZ, MAL, L1CAM — peripheral nerve    (high)
  "12" = "Unknown"               # no significant markers — likely low-quality
)

# Map only if every existing cluster has an entry (guards against a resolution
# change producing different cluster IDs).
stopifnot(all(levels(xenium$seurat_clusters) %in% names(cluster_annotation)))
# unname() is essential: a *named* vector makes Seurat match the names to cell
# barcodes ("No cell overlap ..." error). Unnamed -> assigned per cell in order.
ct_vec <- unname(cluster_annotation[as.character(xenium$seurat_clusters)])
xenium$cell_type <- factor(ct_vec, levels = unique(cluster_annotation))

cat("\n── Cells per assigned cell type ──\n")
print(table(xenium$cell_type))

# Annotated spatial map — same tall orientation & style as the cluster map above
n_ct          <- nlevels(xenium$cell_type)
celltype_cols <- setNames(DiscretePalette(n_ct, palette = "polychrome"),
                          levels(xenium$cell_type))
spatial_ct_df <- data.frame(x = xenium$y_centroid,   # swapped for tall orientation
                            y = xenium$x_centroid,
                            cell_type = xenium$cell_type)
p_spatial_ct <- ggplot(spatial_ct_df, aes(x = x, y = y, color = cell_type)) +
  geom_point(size = 0.25, stroke = 0, alpha = 0.9) +
  scale_color_manual(values = celltype_cols, name = "Cell type") +
  coord_fixed() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1)) +
  labs(title = "Spatial position of annotated cell types") +
  theme_pub +
  theme(axis.line  = element_blank(), axis.text  = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank())
print(p_spatial_ct)   # printed to plot pane only
# To export: ggsave(file.path(output_dir, "annot_spatial_celltypes.png"),
#   p_spatial_ct, width = 9, height = 11, dpi = 300, bg = "white")

# Annotated UMAP (printed)
p_umap_ct <- DimPlot(xenium, group.by = "cell_type", label = TRUE,
                     label.size = 3, repel = TRUE) +
  ggtitle("Annotated cell types") + theme_pub
print(p_umap_ct)


# =============================================================================
# SECTION 3c — GDF15 EXPRESSION PER CELL TYPE  (uses the annotated cell_type)
# =============================================================================
# % of cells expressing GDF15 in each annotated cell type, coloured by the
# number of GDF15+ cells. Same style as script 08's per-cell-type version.
# (Clusters sharing a label — e.g. 0 & 1 -> Epithelial_Luminal — are merged.)
if ("GDF15" %in% rownames(xenium)) {

  gdf15_vec <- as.numeric(
    GetAssayData(xenium, assay = "Xenium", layer = "data")["GDF15", ])

  gdf15_by_celltype <- data.frame(
    celltype = as.character(xenium$cell_type),
    pos      = gdf15_vec > 0
  ) %>%
    group_by(celltype) %>%
    summarise(n_cells        = n(),
              n_gdf15_pos    = sum(pos),
              pct_expressing = 100 * mean(pos),
              .groups = "drop")

  # order bars so the highest % sits at the top (matches the reference figure)
  gdf15_by_celltype$celltype_f <- factor(
    gdf15_by_celltype$celltype,
    levels = gdf15_by_celltype$celltype[order(gdf15_by_celltype$pct_expressing)])

  cat("\n── GDF15 expression per cell type ──\n")
  print(gdf15_by_celltype[order(-gdf15_by_celltype$pct_expressing),
                          c("celltype", "n_cells", "n_gdf15_pos", "pct_expressing")],
        row.names = FALSE, digits = 3)

  p_gdf15_celltype <- ggplot(gdf15_by_celltype,
      aes(x = pct_expressing, y = celltype_f, fill = n_gdf15_pos)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", pct_expressing)),
              hjust = -0.2, size = 3.5) +
    scale_fill_viridis(option = "mako", name = "# GDF15+\ncells") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = "% of each cell type expressing GDF15",
         x = "% cells with GDF15 > 0", y = NULL) +
    theme_pub
  print(p_gdf15_celltype)   # plot pane only
  # To export: ggsave(file.path(output_dir, "annot_GDF15_pct_by_celltype.png"),
  #   p_gdf15_celltype, width = 8, height = 6, dpi = 300, bg = "white")

} else {
  message("GDF15 not in panel — skipping per-cell-type GDF15 plot.")
}


# =============================================================================
# SECTION 4 — FIND ALL MARKERS  (data-driven cluster characterisation)
# =============================================================================
message("[4] FindAllMarkers on unsupervised clusters ...")

Idents(xenium) <- "seurat_clusters"

# Positive markers per cluster (Wilcoxon). only.pos = TRUE keeps up-genes.
all_markers <- FindAllMarkers(
  xenium,
  assay           = "Xenium",
  only.pos        = TRUE,
  min.pct         = 0.10,
  logfc.threshold = 0.25,
  test.use        = "wilcox"
)

write.csv(all_markers,
          file.path(output_dir, "annot_FindAllMarkers_all.csv"),
          row.names = FALSE)
cat(sprintf("FindAllMarkers: %d marker rows across %d clusters\n",
            nrow(all_markers), nlevels(xenium$seurat_clusters)))

# Top-N markers per cluster (ranked by avg_log2FC among padj-significant genes)
top_markers <- all_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = N_TOP_MARKERS, with_ties = FALSE) %>%
  ungroup()

write.csv(top_markers,
          file.path(output_dir, "annot_FindAllMarkers_topN.csv"),
          row.names = FALSE)

# Print top markers per cluster to the console for quick inspection
cat("\n── Top markers per cluster (avg_log2FC, padj<0.05) ───────────────────────\n")
for (cl in levels(xenium$seurat_clusters)) {
  g <- top_markers$gene[top_markers$cluster == cl]
  cat(sprintf("Cluster %-3s: %s\n", cl, paste(g, collapse = ", ")))
}
cat("\n")


# =============================================================================
# SECTION 5 — MARKER VISUALS (help you map clusters -> cell types)
# =============================================================================
message("[5] Marker heatmap + dot plots ...")

top_genes <- unique(top_markers$gene)

# Heatmap needs these genes in scale.data; scale just the marker genes (cheap),
# and downsample cells per cluster so the heatmap stays legible.
xenium <- ScaleData(xenium, features = top_genes, verbose = FALSE)

p_heat <- DoHeatmap(
  subset(xenium, downsample = 100),
  features = top_genes,
  size     = 3
) +
  scale_fill_viridis(option = "magma") +
  ggtitle(sprintf("Top %d FindAllMarkers genes per cluster", N_TOP_MARKERS))
ggsave(file.path(output_dir, "annot_marker_heatmap.png"),
       p_heat, width = 12, height = max(8, 0.18 * length(top_genes)),
       dpi = 300, bg = "white", limitsize = FALSE)
print(p_heat)

# Dot plot of the data-driven top markers (top 5 per cluster to keep readable)
top5 <- all_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE) %>%
  ungroup()
top5_genes <- unique(top5$gene)

p_dot_data <- DotPlot(xenium, features = top5_genes, group.by = "seurat_clusters",
                      dot.scale = 5) +
  scale_color_viridis(option = "plasma") +
  ggtitle("Data-driven top-5 markers per cluster") +
  xlab("Gene") + ylab("Cluster") + theme_pub +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7))
ggsave(file.path(output_dir, "annot_dotplot_datadriven_markers.png"),
       p_dot_data, width = max(14, 0.20 * length(top5_genes)), height = 7.5,
       dpi = 300, bg = "white", limitsize = FALSE)
print(p_dot_data)

# Reference dot plot of the curated prostate panel (GUIDE ONLY — not the call)
panel_genes <- intersect(unique(unlist(cell_type_markers)), rownames(xenium))
p_dot_ref <- DotPlot(xenium, features = panel_genes, group.by = "seurat_clusters",
                     dot.scale = 5) +
  scale_color_viridis(option = "viridis") +
  ggtitle("Reference: curated prostate marker panel by cluster") +
  xlab("Gene") + ylab("Cluster") + theme_pub +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7))
ggsave(file.path(output_dir, "annot_dotplot_curated_reference.png"),
       p_dot_ref, width = max(14, 0.22 * length(panel_genes)), height = 7.5,
       dpi = 300, bg = "white", limitsize = FALSE)
print(p_dot_ref)

# ── Feature plots: key macrophage markers on the whole-tissue UMAP ─────────────
# CD68 (pan-macrophage), MRC1 (= CD206, M2/TAM), CD200 — shown across all cells
# so you can see which cluster(s) carry the macrophage programme.
fp_genes <- intersect(c("CD68", "CD200", "MRC1"), rownames(xenium))
if (length(fp_genes) > 0) {
  p_feat <- FeaturePlot(xenium, features = fp_genes, order = TRUE,
                        ncol = length(fp_genes)) &
    scale_color_viridis(option = "inferno") &
    theme_pub &
    theme(axis.text = element_blank(), axis.ticks = element_blank())
  ggsave(file.path(output_dir, "annot_featureplot_CD68_CD200_MRC1.png"),
         p_feat, width = 5 * length(fp_genes), height = 5, dpi = 300, bg = "white")
  print(p_feat)

  fp_missing <- setdiff(c("CD68", "CD200", "MRC1"), fp_genes)
  if (length(fp_missing) > 0)
    message("  Feature plot: not in panel, skipped: ",
            paste(fp_missing, collapse = ", "))
}


# ── Save the annotated object for downstream scripts (e.g. neighbourhood) ─────
# Carries xenium$cell_type (from Section 3b) so later scripts can reuse the new
# annotation without re-clustering.
saveRDS(xenium, file.path(output_dir, "annot_xenium_celltype_annotated.rds"))
message("Annotated object saved -> annot_xenium_celltype_annotated.rds")

message("Done. Clusters printed; marker tables & plots written to: ", output_dir)
