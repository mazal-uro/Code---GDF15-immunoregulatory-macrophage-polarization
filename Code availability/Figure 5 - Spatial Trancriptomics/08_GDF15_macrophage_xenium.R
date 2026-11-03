# =============================================================================
# 08_GDF15_macrophage_xenium.R
# GDF15 Expression in Macrophages — 10x Xenium Prime Human Prostate Adenocarcinoma
#
# Analyses:
#   1. Load & QC the Xenium dataset
#   2. Normalize, cluster, and annotate cell types
#   3. Subset macrophages and subcluster
#   4. Score macrophage phenotypes (Tissue Resident, Infiltrating, M1, M2, LAM)
#   5. GDF15 expression by subcluster (violin + dot plot)
#   6. Spatial distribution of GDF15+ macrophages
#   7. Genome-wide Spearman correlation with GDF15 in macrophages
#   8. Marker scatter plots for top correlated genes
#
# PRE-REQUISITE: Run 00_convert_xenium_zarr.py first to convert the zarr
#   bundle to standard MTX + cells.csv format. That script produces:
#     ~/Documents/Xenium/Xenium_10x_converted/cell_feature_matrix/
#     ~/Documents/Xenium/Xenium_10x_converted/cells.csv.gz
#
# Required packages: Seurat (>=5.0), ggplot2, dplyr, patchwork, viridis,
#                    ggrepel, corrplot, Matrix
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(viridis)
library(ggrepel)
library(corrplot)
library(Matrix)
library(cowplot)   # loaded up front so theme_cowplot() is available everywhere
                   # (previously only loaded in Section 7c, after Section 6 used it)

# ── Paths ─────────────────────────────────────────────────────────────────────
# Point to the converted output (no spaces in path — avoids data.table::fread bug)
converted_path <- "/Users/alekhyamazumdar/Documents/Xenium/Xenium_10x_converted"
mtx_path       <- file.path(converted_path, "cell_feature_matrix")
cells_csv      <- file.path(converted_path, "cells.csv.gz")
output_dir     <- "/Users/alekhyamazumdar/Documents/Xenium/Output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ── Resume toggle ─────────────────────────────────────────────────────────────
# TRUE  = load the Seurat objects saved by a previous run (Section 14 outputs)
#         and SKIP the heavy rebuild in Sections 1–6 (load, QC, cluster, annotate,
#         subset, subcluster). Everything from Section 4b onward still runs.
# FALSE = run the full pipeline from the converted MTX (original behaviour).
RESUME_FROM_RDS <- TRUE
xenium_rds <- file.path(output_dir, "08_xenium_full_annotated.rds")
macro_rds  <- file.path(output_dir, "08_macrophage_GDF15_xenium.rds")

# Sanity check — fail early if conversion hasn't been run
# (only relevant on the full-pipeline path; the RESUME path never touches the MTX)
if (!RESUME_FROM_RDS && !dir.exists(mtx_path)) {
  stop(paste0(
    "\nConverted data not found at: ", mtx_path,
    "\nPlease run  00_convert_xenium_zarr.py  first:\n",
    "  pip install zarr scipy\n",
    "  python3 /Users/alekhyamazumdar/Documents/Xenium/00_convert_xenium_zarr.py\n"
  ))
}

# ── Plotting theme ─────────────────────────────────────────────────────────────
theme_pub <- theme_classic(base_size = 12) +
  theme(
    plot.title   = element_text(face = "bold", size = 13),
    axis.text    = element_text(color = "black"),
    legend.title = element_text(face = "bold")
  )


# ── Constants used by downstream sections ─────────────────────────────────────
# Defined here (not only inside Sections 4/6) so the RESUME path — which skips
# those sections — still has them. In the full run they are simply re-defined
# identically inside their sections, which is harmless.

# Tableau 10 — muted, high-contrast, publication-standard cluster palette
# (used in Sections 7 and 8)
tableau10 <- c(
  "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
  "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"
)

# Broad marker panel for human prostate tissue
# (used for auto-annotation in Section 4 and for Supplementary Table S1 at the end)
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
# RESUME PATH — load saved objects instead of rebuilding (Sections 1–6)
# =============================================================================
if (RESUME_FROM_RDS) {
  message("RESUME_FROM_RDS = TRUE — loading saved Seurat objects; skipping Sections 1–6.")

  if (!file.exists(xenium_rds) || !file.exists(macro_rds)) {
    stop(paste0(
      "\nRESUME_FROM_RDS is TRUE but the saved objects were not found:\n",
      "  ", xenium_rds, "\n",
      "  ", macro_rds, "\n",
      "Set RESUME_FROM_RDS <- FALSE to rebuild them from the converted MTX,\n",
      "or check that Section 14 of a previous run saved them."
    ))
  }

  xenium <- readRDS(xenium_rds)
  macro  <- readRDS(macro_rds)

  # The downstream code expects these metadata columns to be present in the
  # saved objects — fail loudly now rather than deep inside a later section.
  if (!"auto_celltype"   %in% colnames(xenium@meta.data))
    stop("Loaded 'xenium' is missing the 'auto_celltype' column (was it saved after Section 4?).")
  if (!"seurat_clusters" %in% colnames(macro@meta.data))
    stop("Loaded 'macro' is missing 'seurat_clusters' (was it saved after Section 6?).")

  cat(sprintf("Loaded xenium : %d cells x %d genes\n", ncol(xenium), nrow(xenium)))
  cat(sprintf("Loaded macro  : %d cells x %d genes\n", ncol(macro),  nrow(macro)))
  cat(sprintf("GDF15 in panel: %s\n", "GDF15" %in% rownames(xenium)))
  message(
    "Note: figures produced *inside* the skipped build steps are NOT regenerated ",
    "on the RESUME path — specifically the Section 4 cluster/annotation panels ",
    "(08_umap_clusters, 08_umap_celltypes, 08_spatial_clusters, 08_spatial_celltypes, ",
    "08_dotplot_broad_markers) and the Section 6 macrophage UMAP + phenotype heatmap ",
    "(08_mac_umap_subclusters, 08_mac_phenotype_heatmap). Those already exist in Output/. ",
    "Everything from Section 4b onward (GDF15, neighbourhood, pathway, zoom figures) is regenerated. ",
    "Set RESUME_FROM_RDS <- FALSE for a full rebuild that also refreshes those figures."
  )
}


# =============================================================================
# SECTIONS 1–6 — BUILD OBJECTS (skipped when RESUME_FROM_RDS = TRUE)
# =============================================================================
if (!RESUME_FROM_RDS) {

# =============================================================================
# SECTION 1 — LOAD XENIUM DATA
# =============================================================================
message("Loading Xenium data from converted MTX ...")

# ── 1a. Expression matrix ─────────────────────────────────────────────────────
counts <- Read10X(data.dir = mtx_path)
xenium <- CreateSeuratObject(
  counts = counts,
  assay  = "Xenium",
  project = "XeniumProstatePrime5K"
)

# ── 1b. Cell spatial coordinates ─────────────────────────────────────────────
cells_meta <- read.csv(cells_csv)

# Align barcodes (Read10X uses the barcode as cell name)
common_cells <- intersect(colnames(xenium), cells_meta$cell_id)
cells_meta   <- cells_meta[match(common_cells, cells_meta$cell_id), ]
xenium       <- xenium[, common_cells]

# Add spatial metadata
xenium <- AddMetaData(xenium, metadata = cells_meta$x_centroid,  col.name = "x_centroid")
xenium <- AddMetaData(xenium, metadata = cells_meta$y_centroid,  col.name = "y_centroid")
xenium <- AddMetaData(xenium, metadata = cells_meta$cell_area,   col.name = "cell_area")

# Create a FOV (Field of View) for spatial plots
# Seurat v5 expects a FOV stored as a data.frame with columns x, y
coords_df <- data.frame(
  row.names   = common_cells,
  x           = cells_meta$x_centroid,
  y           = cells_meta$y_centroid
)
xenium[["fov"]] <- CreateFOV(
  coords      = coords_df,
  type        = "centroids",
  assay       = "Xenium"
)

# Quick sanity check
cat(sprintf("Cells loaded : %d\n", ncol(xenium)))
cat(sprintf("Genes in panel: %d\n", nrow(xenium)))
cat(sprintf("GDF15 in panel: %s\n", "GDF15" %in% rownames(xenium)))


# =============================================================================
# SECTION 2 — QC FILTERING
# =============================================================================
message("Running QC ...")

# Transcript count distribution
# Note: CreateSeuratObject names the count column "nCount_Xenium" automatically
p_qc <- VlnPlot(xenium, features = "nCount_Xenium", pt.size = 0) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
  scale_y_log10() +
  ggtitle("Transcript counts per cell (Xenium)") +
  theme_pub
ggsave(file.path(output_dir, "08_qc_transcript_counts.png"), p_qc, width = 6, height = 5)
print(p_qc)

# Also plot cell area distribution
p_area <- ggplot(xenium@meta.data, aes(x = cell_area)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "white") +
  scale_x_log10() +
  geom_vline(xintercept = c(10, 5000), linetype = "dashed", color = "red") +
  labs(title = "Cell area distribution", x = "Cell area (µm²)", y = "Count") +
  theme_pub
ggsave(file.path(output_dir, "08_qc_cell_area.png"), p_area, width = 6, height = 4)
print(p_area)

# Filter: >= 10 transcripts and reasonable cell area
xenium <- subset(xenium, subset = nCount_Xenium >= 10 & cell_area >= 10 & cell_area <= 5000)
cat(sprintf("Cells after QC: %d\n", ncol(xenium)))


# =============================================================================
# SECTION 3 — NORMALIZE, REDUCE, CLUSTER
# =============================================================================
message("Normalizing and clustering ...")

# CLR normalization is recommended for Xenium panel data
xenium <- NormalizeData(xenium, normalization.method = "CLR", margin = 2)
xenium <- FindVariableFeatures(xenium, nfeatures = 2000)
xenium <- ScaleData(xenium)
xenium <- RunPCA(xenium, npcs = 30, verbose = FALSE)

# Elbow plot to guide choice of PCs
p_elbow <- ElbowPlot(xenium, ndims = 30) + theme_pub
ggsave(file.path(output_dir, "08_pca_elbow.png"), p_elbow, width = 6, height = 4)
print(p_elbow)

xenium <- FindNeighbors(xenium, reduction = "pca", dims = 1:20)
xenium <- FindClusters(xenium,  resolution = 0.4, verbose = FALSE)
xenium <- RunUMAP(xenium, dims = 1:20, verbose = FALSE)


# =============================================================================
# SECTION 4 — BROAD CELL TYPE ANNOTATION
# =============================================================================
message("Annotating broad cell types ...")

# Broad marker panel for human prostate tissue
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

# Score all cell types
for (ct in names(cell_type_markers)) {
  genes_present <- intersect(cell_type_markers[[ct]], rownames(xenium))
  if (length(genes_present) > 0) {
    xenium <- AddModuleScore(xenium, features = list(genes_present), name = ct)
  }
}

# UMAP colored by clusters
p_umap <- DimPlot(xenium, label = TRUE, label.size = 3) +
  ggtitle("Seurat clusters (all cells)") +
  theme_pub
ggsave(file.path(output_dir, "08_umap_clusters.png"), p_umap, width = 7, height = 6)
print(p_umap)

# Dot plot of cell type markers across clusters — helps assign identities manually
p_dot_broad <- DotPlot(
  xenium,
  features = c("KRT8", "KRT18", "AR",        # Luminal epithelial
               "KRT5", "TP63",                 # Basal epithelial
               "DCN", "COL1A1",                # Fibroblast
               "ACTA2", "MYH11",               # Smooth muscle
               "PECAM1", "CDH5",               # Endothelial
               "CD68", "CD163", "MRC1", "CSF1R", # Macrophage
               "CD3E", "CD8A", "CD4",          # T cell
               "MS4A1", "CD79A",               # B cell
               "KIT", "TPSAB1")               # Mast cell
) +
  RotatedAxis() +
  ggtitle("Marker expression by cluster (broad annotation)") +
  theme_pub
ggsave(file.path(output_dir, "08_dotplot_broad_markers.png"), p_dot_broad,
       width = 16, height = 8)
print(p_dot_broad)

# Spatial overview: cell clusters on tissue
p_spatial_clusters <- ImageDimPlot(xenium, fov = "fov", size = 0.2, alpha = 0.7) +
  ggtitle("Seurat clusters — tissue coordinates") +
  theme_pub
ggsave(file.path(output_dir, "08_spatial_clusters.png"), p_spatial_clusters,
       width = 10, height = 10)
print(p_spatial_clusters)

# ── Manual annotation step ────────────────────────────────────────────────────
# Inspect the dot plot (08_dotplot_broad_markers.png) and the UMAP above,
# then fill in the cluster-to-celltype mapping below.
#
# Example (replace with your actual cluster IDs after inspection):
#   new_labels <- c("0"  = "Luminal_Epithelial",
#                   "1"  = "Macrophage",
#                   "2"  = "Fibroblast", ...)
#
# For now we use a marker-score-based auto-annotation as a starting point:

score_cols <- paste0(names(cell_type_markers), "1")
score_cols <- score_cols[score_cols %in% colnames(xenium@meta.data)]

xenium$auto_celltype <- apply(
  xenium@meta.data[, score_cols, drop = FALSE], 1,
  function(x) sub("1$", "", names(which.max(x)))
)

p_umap_ct <- DimPlot(xenium, group.by = "auto_celltype", label = TRUE,
                      label.size = 3, repel = TRUE) +
  ggtitle("Auto-annotated cell types") +
  theme_pub
ggsave(file.path(output_dir, "08_umap_celltypes.png"), p_umap_ct, width = 9, height = 7)
print(p_umap_ct)

p_spatial_ct <- ImageDimPlot(xenium, fov = "fov", group.by = "auto_celltype",
                               size = 0.2, alpha = 0.7) +
  ggtitle("Cell types — tissue coordinates") +
  theme_pub
ggsave(file.path(output_dir, "08_spatial_celltypes.png"), p_spatial_ct,
       width = 10, height = 10)
print(p_spatial_ct)

}  # ── end SECTIONS 1–4 build block (skipped when RESUME_FROM_RDS = TRUE) ──────


# =============================================================================
# SECTION 4b — GDF15 EXPRESSION ACROSS ALL CELL TYPES
# =============================================================================
message("Profiling GDF15 expression across all cell types ...")

# Pull GDF15 expression for every cell
all_gdf15 <- as.numeric(
  GetAssayData(xenium, assay = "Xenium", layer = "data")["GDF15", ]
)
names(all_gdf15) <- colnames(xenium)

# Tag each cell as GDF15+ (any detectable expression) or GDF15-
xenium$GDF15_expressed <- all_gdf15 > 0
xenium$GDF15_expr_val  <- all_gdf15

n_pos   <- sum(xenium$GDF15_expressed)
pct_pos <- 100 * mean(xenium$GDF15_expressed)
cat(sprintf("GDF15+ cells (expr > 0): %d / %d  (%.1f%% of all cells)\n",
            n_pos, ncol(xenium), pct_pos))

# ── Build summary table ────────────────────────────────────────────────────────
gdf15_meta <- data.frame(
  celltype     = xenium$auto_celltype,
  gdf15_val    = all_gdf15,
  gdf15_pos    = all_gdf15 > 0,
  row.names    = colnames(xenium)
)

# Per cell type: total cells, GDF15+ count, % expressing, mean expression in +ve cells
gdf15_summary <- gdf15_meta %>%
  group_by(celltype) %>%
  summarise(
    n_cells      = n(),
    n_gdf15_pos  = sum(gdf15_pos),
    pct_expressing = 100 * mean(gdf15_pos),
    mean_expr_all  = mean(gdf15_val),
    mean_expr_pos  = ifelse(any(gdf15_pos), mean(gdf15_val[gdf15_pos]), NA_real_),
    .groups = "drop"
  ) %>%
  arrange(desc(n_gdf15_pos))

# Print ranking to console
cat("\n── GDF15+ cell counts by type (ranked) ─────────────────────────────────\n")
print(as.data.frame(gdf15_summary), digits = 3)
cat("\n")

write.csv(gdf15_summary,
          file.path(output_dir, "04b_GDF15_celltype_summary.csv"),
          row.names = FALSE)

# ── 4b-i. Bar chart: absolute GDF15+ cell counts per type (ranked) ────────────
ct_order_count <- gdf15_summary$celltype[order(gdf15_summary$n_gdf15_pos)]
gdf15_summary$celltype_f <- factor(gdf15_summary$celltype, levels = ct_order_count)

p_gdf15_counts <- ggplot(gdf15_summary,
  aes(x = n_gdf15_pos, y = celltype_f, fill = pct_expressing)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_gdf15_pos), hjust = -0.2, size = 3.5) +
  scale_fill_viridis(option = "plasma", name = "% expressing") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "GDF15+ cell counts by cell type",
    x     = "Number of GDF15+ cells",
    y     = NULL
  ) +
  theme_pub
ggsave(file.path(output_dir, "04b_GDF15_counts_by_celltype.png"),
       p_gdf15_counts, width = 9, height = 6)
print(p_gdf15_counts)

# ── 4b-ii. Bar chart: % of each cell type that expresses GDF15 ────────────────
ct_order_pct <- gdf15_summary$celltype[order(gdf15_summary$pct_expressing)]
gdf15_summary$celltype_pct_f <- factor(gdf15_summary$celltype, levels = ct_order_pct)

p_gdf15_pct <- ggplot(gdf15_summary,
  aes(x = pct_expressing, y = celltype_pct_f, fill = n_gdf15_pos)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", pct_expressing)), hjust = -0.2, size = 3.5) +
  scale_fill_viridis(option = "mako", name = "# GDF15+\ncells") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "% of each cell type expressing GDF15",
    x     = "% cells with GDF15 > 0",
    y     = NULL
  ) +
  theme_pub
ggsave(file.path(output_dir, "04b_GDF15_pct_by_celltype.png"),
       p_gdf15_pct, width = 9, height = 6)
print(p_gdf15_pct)

# ── 4b-iii. Violin: GDF15 expression level across cell types (GDF15+ only) ─────
gdf15_pos_df <- gdf15_meta[gdf15_meta$gdf15_pos, ]
ct_order_vln <- gdf15_summary$celltype[order(gdf15_summary$mean_expr_pos, na.last = FALSE)]
gdf15_pos_df$celltype_f <- factor(gdf15_pos_df$celltype, levels = ct_order_vln)

p_gdf15_vln <- ggplot(gdf15_pos_df,
  aes(x = gdf15_val, y = celltype_f, fill = celltype_f)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.85) +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA, color = "black") +
  scale_fill_viridis(discrete = TRUE, option = "turbo", guide = "none") +
  labs(
    title = "GDF15 expression level by cell type\n(GDF15+ cells only)",
    x     = "GDF15 normalized expression",
    y     = NULL
  ) +
  theme_pub
ggsave(file.path(output_dir, "04b_GDF15_expression_violin_by_celltype.png"),
       p_gdf15_vln, width = 9, height = 7)
print(p_gdf15_vln)

# ── 4b-iv. Stacked bar: GDF15+ vs GDF15- composition within each cell type ─────
gdf15_stacked <- gdf15_meta %>%
  group_by(celltype) %>%
  summarise(
    pos = sum(gdf15_pos),
    neg = sum(!gdf15_pos),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(cols = c("pos", "neg"),
                      names_to = "status", values_to = "count") %>%
  mutate(status = factor(status, levels = c("pos", "neg"),
                          labels = c("GDF15+", "GDF15-")),
         celltype = factor(celltype, levels = ct_order_pct))

p_gdf15_stack <- ggplot(gdf15_stacked,
  aes(x = count, y = celltype, fill = status)) +
  geom_col(position = "fill", width = 0.7) +
  scale_fill_manual(values = c("GDF15+" = "#D62728", "GDF15-" = "#C8C8C8")) +
  scale_x_continuous(labels = scales::percent_format()) +
  labs(
    title = "GDF15 detection rate within each cell type",
    x     = "Proportion of cells",
    y     = NULL,
    fill  = NULL
  ) +
  theme_pub
ggsave(file.path(output_dir, "04b_GDF15_stacked_proportion.png"),
       p_gdf15_stack, width = 9, height = 6)
print(p_gdf15_stack)

# ── 4b-v. Spatial map: all GDF15+ cells colored by cell type ──────────────────
gdf15_pos_barcodes <- names(all_gdf15)[all_gdf15 > 0]

coords_gdf15 <- xenium@meta.data[, c("x_centroid", "y_centroid", "auto_celltype",
                                       "GDF15_expressed")]
# Plot background (GDF15-) first, then GDF15+ cells on top
coords_bg  <- coords_gdf15[!coords_gdf15$GDF15_expressed, ]
coords_pos <- coords_gdf15[ coords_gdf15$GDF15_expressed, ]

p_gdf15_spatial <- ggplot() +
  geom_point(data = coords_bg,
             aes(x = x_centroid, y = y_centroid),
             color = "#D8D8D8", size = 0.1, alpha = 0.25, stroke = 0) +
  geom_point(data = coords_pos,
             aes(x = x_centroid, y = y_centroid, color = auto_celltype),
             size = 0.8, alpha = 0.85, stroke = 0) +
  scale_color_viridis(discrete = TRUE, option = "turbo", name = "Cell type") +
  coord_fixed() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title = "Spatial distribution of GDF15+ cells by cell type") +
  theme_pub +
  theme(axis.line  = element_blank(), axis.text  = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank())
ggsave(file.path(output_dir, "04b_GDF15_spatial_all_celltypes.png"),
       p_gdf15_spatial, width = 10, height = 12)
print(p_gdf15_spatial)


# =============================================================================
# SECTIONS 5–6 — SUBSET + SUBCLUSTER MACROPHAGES
#   (skipped when RESUME_FROM_RDS = TRUE; `macro` is loaded from RDS instead)
# =============================================================================
if (!RESUME_FROM_RDS) {

# =============================================================================
# SECTION 5 — SUBSET MACROPHAGES
# =============================================================================
message("Subsetting macrophages ...")

# Select macrophage cells based on auto-annotation
# (adjust the string if your annotation uses a different name)
mac_cells <- WhichCells(xenium, expression = auto_celltype == "Macrophage")

# Fallback: if auto-annotation misses some, also grab any cluster where
# both CD68 and CSF1R score are high
if (length(mac_cells) < 50) {
  message("Auto-annotation found few macrophages; falling back to marker-based selection.")
  cd68_expr  <- GetAssayData(xenium, assay = "Xenium", layer = "data")["CD68", ]
  csf1r_expr <- GetAssayData(xenium, assay = "Xenium", layer = "data")["CSF1R", ]
  mac_cells  <- names(which(cd68_expr > 0 & csf1r_expr > 0))
}

cat(sprintf("Macrophages identified: %d\n", length(mac_cells)))
macro <- subset(xenium, cells = mac_cells)


# =============================================================================
# SECTION 6 — MACROPHAGE SUBCLUSTERING
# =============================================================================
message("Subclustering macrophages ...")

macro <- NormalizeData(macro, normalization.method = "CLR", margin = 2)
macro <- FindVariableFeatures(macro, nfeatures = 500)
macro <- ScaleData(macro)
macro <- RunPCA(macro, npcs = 20, verbose = FALSE)
macro <- FindNeighbors(macro, dims = 1:15)
macro <- FindClusters(macro, resolution = 0.2, verbose = FALSE)
macro <- RunUMAP(macro, dims = 1:15, verbose = FALSE)

cat(sprintf("Macrophage subclusters: %d\n", nlevels(Idents(macro))))

umap_cl_df          <- as.data.frame(Embeddings(macro, reduction = "umap"))
colnames(umap_cl_df) <- c("UMAP1", "UMAP2")
umap_cl_df$cluster  <- factor(macro$seurat_clusters)
centroids_cl        <- aggregate(cbind(UMAP1, UMAP2) ~ cluster,
                                  data = umap_cl_df, FUN = mean)

n_cl <- nlevels(umap_cl_df$cluster)
# Tableau 10 — muted, high-contrast, publication standard
tableau10 <- c(
  "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
  "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"
)
cl_palette <- tableau10[seq_len(n_cl)]

p_mac_umap <- ggplot(umap_cl_df, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size = 1.4, alpha = 0.80) +
  geom_label(data = centroids_cl,
             mapping = aes(x = UMAP1, y = UMAP2, label = cluster),
             color = "black", fill = "white", size = 3.5,
             label.size = NA, alpha = 0.80, fontface = "bold",
             inherit.aes = FALSE) +
  scale_color_manual(values = cl_palette, name = "Cluster") +
  coord_equal() +
  labs(x = NULL, y = NULL, title = "Macrophage Subclusters") +
  theme_cowplot(font_size = 11) +
  theme(
    axis.line  = element_blank(),
    axis.ticks = element_blank(),
    axis.text  = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "right"
  )
ggsave(file.path(output_dir, "08_mac_umap_subclusters.png"), p_mac_umap,
       width = 7, height = 6, dpi = 300, bg = "white")
print(p_mac_umap)

# ── Phenotype signature scoring ───────────────────────────────────────────────
tissue_resident_genes <- intersect(
  c("LYVE1", "TIMD4", "FOLR2", "CX3CR1", "MRC1", "MARCO"), rownames(macro))
infiltrating_genes    <- intersect(
  c("CD74", "CSF1R", "VIM", "TGFBI", "CCR2", "SELL"),       rownames(macro))
m1_like_genes         <- intersect(
  c("CD86", "IL1B", "TNF", "CXCL9", "CXCL10", "NOS2"),      rownames(macro))
m2_like_genes         <- intersect(
  c("CD163", "MRC1", "IL10", "TGFB1", "CCL18", "CCL22"),    rownames(macro))
lam_genes             <- intersect(
  c("TREM2", "APOE", "LPL", "FABP4", "GPNMB", "SPP1"),      rownames(macro))

macro <- AddModuleScore(macro, features = list(tissue_resident_genes), name = "TissueResident")
macro <- AddModuleScore(macro, features = list(infiltrating_genes),    name = "Infiltrating")
macro <- AddModuleScore(macro, features = list(m1_like_genes),         name = "M1_like")
macro <- AddModuleScore(macro, features = list(m2_like_genes),         name = "M2_like")
macro <- AddModuleScore(macro, features = list(lam_genes),             name = "LAM")

# Heatmap of phenotype scores across subclusters
sig_scores <- c("TissueResident1", "Infiltrating1", "M1_like1", "M2_like1", "LAM1")
avg_scores <- AggregateExpression(macro, features = NULL,
                                   group.by = "seurat_clusters",
                                   return.seurat = FALSE)
# Score heatmap using AverageExpression on metadata
meta_means <- macro@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(sig_scores), mean), .groups = "drop")

score_mat <- as.matrix(meta_means[, sig_scores])
rownames(score_mat) <- paste0("Cluster ", meta_means$seurat_clusters)
colnames(score_mat) <- c("Tissue Resident", "Infiltrating", "M1-like", "M2-like", "LAM")

png(file.path(output_dir, "08_mac_phenotype_heatmap.png"), width = 800, height = 600)
corrplot(t(score_mat), method = "color", is.corr = FALSE, tl.cex = 1,
         title = "Macrophage phenotype scores by subcluster", mar = c(0,0,2,0))
dev.off()
corrplot(t(score_mat), method = "color", is.corr = FALSE, tl.cex = 1,
         title = "Macrophage phenotype scores by subcluster", mar = c(0,0,2,0))

}  # ── end SECTIONS 5–6 build block (skipped when RESUME_FROM_RDS = TRUE) ──────


# =============================================================================
# SECTION 6b — QC: WHICH WHOLE-TISSUE CLUSTERS EXPRESS CD68?
#   Sanity check on the macrophage subset. Section 5 selected macrophages via
#   auto_celltype == "Macrophage" (a per-cell which.max over module scores, NOT
#   a manual cluster assignment). This block reports, for each ORIGINAL
#   whole-tissue Seurat cluster: macrophage-marker expression, the dominant
#   auto_celltype call, and how much of the cluster actually ended up in `macro`
#   — so a CD68+ cluster that was mis-labelled or missed becomes obvious.
#   Uses the annotated `xenium` + the `macro` subset, so it runs on both paths.
# =============================================================================
message("Section 6b: CD68 / macrophage-marker QC across whole-tissue clusters ...")

mac_markers_qc <- intersect(c("CD68", "CD163", "MRC1", "CSF1R", "ITGAM"),
                            rownames(xenium))
cat(sprintf("Macrophage QC markers present in panel: %s\n",
            paste(mac_markers_qc, collapse = ", ")))

xen_data     <- GetAssayData(xenium, assay = "Xenium", layer = "data")
cluster_chr  <- as.character(xenium$seurat_clusters)
clusters_all <- if (!is.null(levels(xenium$seurat_clusters)))
                  levels(xenium$seurat_clusters) else sort(unique(cluster_chr))
macro_bc     <- Cells(macro)

pct_expr <- function(g, idx) 100 * mean(as.numeric(xen_data[g, idx]) > 0)

qc_tab <- do.call(rbind, lapply(clusters_all, function(cl) {
  idx <- which(cluster_chr == cl)
  bc  <- colnames(xenium)[idx]
  row <- data.frame(
    cluster          = cl,
    n_cells          = length(idx),
    CD68_pct         = pct_expr("CD68", idx),
    CD68_mean        = mean(as.numeric(xen_data["CD68", idx])),
    top_celltype     = names(sort(table(xenium$auto_celltype[idx]),
                                  decreasing = TRUE))[1],
    pct_called_macro = 100 * mean(xenium$auto_celltype[idx] == "Macrophage"),
    pct_in_subset    = 100 * mean(bc %in% macro_bc),
    stringsAsFactors = FALSE
  )
  for (g in setdiff(mac_markers_qc, "CD68"))
    row[[paste0(g, "_pct")]] <- pct_expr(g, idx)
  row
}))
qc_tab <- qc_tab[order(-qc_tab$CD68_pct), ]

cat("\n── Macrophage-marker expression by whole-tissue cluster (ranked by CD68 %) ──\n")
print(qc_tab, row.names = FALSE, digits = 3)
write.csv(qc_tab,
          file.path(output_dir, "06b_CD68_macrophage_QC_by_cluster.csv"),
          row.names = FALSE)

# Flag clusters that look macrophage-like by CD68 but were largely NOT captured.
# Heuristic thresholds — CD68 detected in >= 25% of cells, yet < 50% of the
# cluster ended up in the macrophage subset. Tune to taste.
CD68_PCT_FLAG   <- 25
SUBSET_PCT_FLAG <- 50
missed <- qc_tab[qc_tab$CD68_pct >= CD68_PCT_FLAG &
                 qc_tab$pct_in_subset < SUBSET_PCT_FLAG, ]
if (nrow(missed) > 0) {
  cat(sprintf("\n[FLAG] %d cluster(s) are CD68+ (>=%.0f%%) but <%.0f%% captured in the macrophage subset:\n",
              nrow(missed), CD68_PCT_FLAG, SUBSET_PCT_FLAG))
  print(missed[, c("cluster", "n_cells", "CD68_pct", "top_celltype",
                   "pct_called_macro", "pct_in_subset")],
        row.names = FALSE, digits = 3)
  cat("   -> review these: they may be macrophages that auto_celltype assigned elsewhere.\n")
} else {
  cat(sprintf("\n[OK] No cluster with CD68 >= %.0f%% fell below %.0f%% capture — subset looks complete.\n",
              CD68_PCT_FLAG, SUBSET_PCT_FLAG))
}

# Cross-tabulation: whole-tissue cluster x auto_celltype call
cat("\n── Cross-tabulation: seurat_clusters (rows) x auto_celltype (cols) ──\n")
print(table(cluster = xenium$seurat_clusters, celltype = xenium$auto_celltype))

# ── Visual 1: macrophage markers across all whole-tissue clusters (dot plot) ──
p_qc_dot <- DotPlot(xenium, features = mac_markers_qc,
                    group.by = "seurat_clusters", dot.scale = 6) +
  scale_color_viridis(option = "plasma") +
  RotatedAxis() +
  ggtitle("Macrophage markers across whole-tissue clusters") +
  xlab("Gene") + ylab("Seurat cluster") +
  theme_pub
ggsave(file.path(output_dir, "06b_CD68_macrophage_markers_dotplot.png"),
       p_qc_dot, width = 8, height = 6, dpi = 300, bg = "white")
print(p_qc_dot)

# ── Visual 2: CD68 feature UMAP alongside the whole-tissue cluster UMAP ────────
p_qc_feat <- FeaturePlot(xenium, features = "CD68", order = TRUE) +
  scale_color_viridis(option = "inferno") +
  ggtitle("CD68 expression (all cells)") + theme_pub
p_qc_clus <- DimPlot(xenium, group.by = "seurat_clusters",
                     label = TRUE, label.size = 3) +
  ggtitle("Whole-tissue clusters") + theme_pub
p_qc_umap <- p_qc_feat + p_qc_clus
ggsave(file.path(output_dir, "06b_CD68_featureplot_vs_clusters.png"),
       p_qc_umap, width = 13, height = 6, dpi = 300, bg = "white")
print(p_qc_umap)


# =============================================================================
# SECTION 7 — GDF15 EXPRESSION BY MACROPHAGE SUBCLUSTER
# =============================================================================
message("Plotting GDF15 by subcluster ...")

# ── 7a. Violin plot ───────────────────────────────────────────────────────────
p_vln <- VlnPlot(
  macro, features = "GDF15",
  group.by = "seurat_clusters",
  pt.size  = 0,
  cols     = tableau10[seq_len(nlevels(Idents(macro)))]
) +
  stat_summary(fun = median, geom = "crossbar", width = 0.4,
               color = "black", linewidth = 0.6) +
  ggtitle("GDF15 expression by macrophage subcluster") +
  xlab("Macrophage subcluster") +
  ylab("Normalized expression") +
  theme_pub
ggsave(file.path(output_dir, "08_GDF15_violin_subclusters.png"), p_vln,
       width = 10, height = 5)
print(p_vln)

# ── 7b. Dot plot: GDF15 + key macrophage markers ─────────────────────────────
mac_dot_genes <- c(
  "GDF15",
  # Pan-macrophage
  "CD68", "CSF1R", "MRC1",
  # TAM / immunosuppressive
  "CD163", "TGFB1", "IL10", "VEGFA",
  # LAM
  "TREM2", "APOE", "SPP1", "GPNMB",
  # M1 / inflammatory
  "CD86", "IL1B", "CXCL9",
  # Other
  "HIF1A", "TP53", "CD274"
)
mac_dot_genes <- intersect(mac_dot_genes, rownames(macro))

p_dot <- DotPlot(macro, features = mac_dot_genes, group.by = "seurat_clusters",
                  dot.scale = 6) +
  scale_color_viridis(option = "plasma") +
  RotatedAxis() +
  ggtitle("GDF15 and macrophage marker expression by subcluster") +
  theme_pub
ggsave(file.path(output_dir, "08_GDF15_dotplot_subclusters.png"), p_dot,
       width = 14, height = 6)
print(p_dot)

# ── 7c. Feature plot on UMAP (cowplot publication style) ─────────────────────
# Open circles + black outline, white → teal → navy colour scale.
library(cowplot)

# ── Build data frame ──────────────────────────────────────────────────────────
umap_df <- as.data.frame(Embeddings(macro, reduction = "umap"))
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$GDF15   <- as.numeric(
  GetAssayData(macro, assay = "Xenium", layer = "data")["GDF15", ]
)
umap_df$cluster <- as.character(macro$seurat_clusters)
umap_df <- umap_df[order(umap_df$GDF15), ]   # zeros rendered first

# ── Colour scale ──────────────────────────────────────────────────────────────
gdf15_max  <- max(umap_df$GDF15)
scale_vals <- scales::rescale(c(0, 0.02, 1, 2, gdf15_max), from = c(0, gdf15_max))

# ── Feature plot ──────────────────────────────────────────────────────────────
p_feat_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, fill = GDF15)) +
  geom_point(shape = 21, size = 1.9, stroke = 0.28, color = "black") +
  scale_fill_gradientn(
    colors = c("white", "#CDEAF5", "#6BB8D4", "#1D7DB5", "#0B3060"),
    values = scale_vals,
    limits = c(0, gdf15_max),
    breaks = c(0, 1, 2, 3),
    name   = "GDF15 RNA\nExpression"
  ) +
  coord_equal() +
  theme_cowplot(font_size = 11) +
  theme(
    axis.line  = element_blank(),
    axis.ticks = element_blank(),
    axis.text  = element_blank(),
    axis.title = element_blank(),
    legend.position   = "bottom",
    legend.title      = element_text(size = 9, face = "bold", hjust = 0.5),
    legend.text       = element_text(size = 8.5),
    legend.key.height = unit(0.38, "cm"),
    legend.key.width  = unit(1.6, "cm"),
    plot.background   = element_rect(fill = "white", color = NA),
    plot.margin       = margin(10, 10, 10, 10)
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.hjust    = 0.5,
    barwidth       = unit(4.5, "cm"),
    barheight      = unit(0.4, "cm"),
    ticks.colour   = "grey30",
    frame.colour   = "grey30"
  ))

ggsave(file.path(output_dir, "08_GDF15_featureplot_umap.png"),
       p_feat_umap, width = 5, height = 6, dpi = 300, bg = "white")
print(p_feat_umap)


# ── 7d. Feature plots on UMAP for additional macrophage markers ───────────────
# Same cowplot publication style as the GDF15 feature plot (7c), applied to a
# set of macrophage markers, one PNG per gene.
#
# NOTE: CD206 is a protein/CD alias of the gene MRC1 — they are the SAME feature
# in the Xenium panel. A single MRC1 plot (labelled "MRC1 (CD206)") covers both.
# Panel check confirmed present: MRC1, CD163, CD200.
marker_umap_genes  <- c("MRC1", "CD163", "CD200")   # CD206 == MRC1
marker_umap_titles <- c(
  MRC1  = "MRC1 (CD206) RNA\nExpression",
  CD163 = "CD163 RNA\nExpression",
  CD200 = "CD200 RNA\nExpression"
)

make_marker_feature_umap <- function(gene, title_lab) {
  df <- as.data.frame(Embeddings(macro, reduction = "umap"))
  colnames(df) <- c("UMAP1", "UMAP2")
  df$expr <- as.numeric(
    GetAssayData(macro, assay = "Xenium", layer = "data")[gene, ]
  )
  df <- df[order(df$expr), ]                       # zeros rendered first
  gmax <- max(df$expr)
  if (gmax < 1e-9) gmax <- 1                        # guard against all-zero genes

  ggplot(df, aes(x = UMAP1, y = UMAP2, fill = expr)) +
    geom_point(shape = 21, size = 1.9, stroke = 0.28, color = "black") +
    scale_fill_gradientn(
      colors = c("white", "#CDEAF5", "#6BB8D4", "#1D7DB5", "#0B3060"),
      # low-end emphasis (mirrors 7c) but expressed as fractions of the range so
      # it stays valid for any per-gene maximum, large or small
      values = c(0, 0.01, 0.5, 0.8, 1),
      limits = c(0, gmax),
      breaks = pretty(c(0, gmax), n = 3),
      name   = title_lab
    ) +
    coord_equal() +
    theme_cowplot(font_size = 11) +
    theme(
      axis.line  = element_blank(),
      axis.ticks = element_blank(),
      axis.text  = element_blank(),
      axis.title = element_blank(),
      legend.position   = "bottom",
      legend.title      = element_text(size = 9, face = "bold", hjust = 0.5),
      legend.text       = element_text(size = 8.5),
      legend.key.height = unit(0.38, "cm"),
      legend.key.width  = unit(1.6, "cm"),
      plot.background   = element_rect(fill = "white", color = NA),
      plot.margin       = margin(10, 10, 10, 10)
    ) +
    guides(fill = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = unit(4.5, "cm"),
      barheight      = unit(0.4, "cm"),
      ticks.colour   = "grey30",
      frame.colour   = "grey30"
    ))
}

for (g in marker_umap_genes) {
  if (!g %in% rownames(macro)) {
    message(sprintf("  [7d] '%s' not in the panel / macrophage object — skipping.", g))
    next
  }
  p_g   <- make_marker_feature_umap(g, marker_umap_titles[[g]])
  out_g <- file.path(output_dir, sprintf("08_%s_featureplot_umap.png", g))
  ggsave(out_g, p_g, width = 5, height = 6, dpi = 300, bg = "white")
  print(p_g)
  cat(sprintf("  [7d] saved %s\n", basename(out_g)))
}


# =============================================================================
# SECTION 8 — SPATIAL DISTRIBUTION OF GDF15+ MACROPHAGES
# =============================================================================
message("Plotting spatial distribution of GDF15+ macrophages ...")

# ── 8a. GDF15 expression in macrophages on tissue ─────────────────────────────
p_spatial_mac_gdf15 <- ImageFeaturePlot(
  macro, fov = "fov",
  features = "GDF15",
  max.cutoff = "q95",
  size = 1.2
) +
  scale_fill_viridis(option = "inferno", name = "GDF15") +
  ggtitle("GDF15 expression — macrophages on tissue") +
  theme_pub
ggsave(file.path(output_dir, "08_spatial_GDF15_macrophages.png"),
       p_spatial_mac_gdf15, width = 10, height = 10)
print(p_spatial_mac_gdf15)

# ── 8b. GDF15-high vs GDF15-low macrophages on tissue ────────────────────────
gdf15_expr <- GetAssayData(macro, assay = "Xenium", layer = "data")["GDF15", ]
gdf15_threshold <- quantile(gdf15_expr[gdf15_expr > 0], 0.5)  # median of expressing cells

macro$GDF15_status <- ifelse(
  gdf15_expr > gdf15_threshold, "GDF15-high", "GDF15-low/neg"
)
cat(sprintf("GDF15-high macrophages: %d (%.1f%%)\n",
            sum(macro$GDF15_status == "GDF15-high"),
            100 * mean(macro$GDF15_status == "GDF15-high")))

p_spatial_gdf15_cat <- ImageDimPlot(
  macro, fov = "fov",
  group.by = "GDF15_status",
  cols = c("GDF15-high" = "#D62728", "GDF15-low/neg" = "#AEC7E8"),
  size = 1.2, alpha = 0.9
) +
  ggtitle("GDF15-high vs GDF15-low macrophages on tissue") +
  theme_pub
ggsave(file.path(output_dir, "08_spatial_GDF15_high_low_macrophages.png"),
       p_spatial_gdf15_cat, width = 10, height = 10)
print(p_spatial_gdf15_cat)

# ── 8b_overlay. Full tissue (cell types, semi-transparent) + macrophage GDF15 ─
# Build a single coordinate frame: background cells keep their auto_celltype
# colour; macrophage rows are overwritten with their GDF15 status label so they
# render on top with distinct colours.
coords_ov <- xenium@meta.data[, c("x_centroid", "y_centroid", "auto_celltype")]
coords_ov$display_group <- as.character(coords_ov$auto_celltype)

# Substitute macrophage rows with GDF15 status (already computed in 8b above)
mac_cells_ov <- Cells(macro)
mac_gdf15_ov <- setNames(macro$GDF15_status, mac_cells_ov)
mac_rows_ov  <- rownames(coords_ov) %in% mac_cells_ov
coords_ov$display_group[mac_rows_ov] <-
  ifelse(mac_gdf15_ov[rownames(coords_ov)[mac_rows_ov]] == "GDF15-high",
         "Macrophage (GDF15-high)", "Macrophage (GDF15-low/neg)")

# Render order: all background cell types first, then GDF15-low macs, GDF15-high last
mac_ov_groups <- c("Macrophage (GDF15-low/neg)", "Macrophage (GDF15-high)")
bg_ov_groups  <- sort(setdiff(unique(coords_ov$display_group), mac_ov_groups))
level_ov      <- c(bg_ov_groups, mac_ov_groups)
coords_ov$display_group <- factor(coords_ov$display_group, levels = level_ov)
coords_ov <- coords_ov[order(coords_ov$display_group), ]

# Colours: background cell types share the same hue-wheel palette as the
# original cell-type spatial plot; macrophages get strong distinct colours.
n_bg_ov    <- length(bg_ov_groups)
bg_ov_cols <- setNames(scales::hue_pal()(n_bg_ov), bg_ov_groups)
mac_ov_cols <- c(
  "Macrophage (GDF15-low/neg)" = "#4E9AC7",
  "Macrophage (GDF15-high)"    = "#D62728"
)
all_ov_cols  <- c(bg_ov_cols, mac_ov_cols)

# Size / alpha: background tiny; GDF15-low macs small+faded; GDF15-high macs prominent
all_ov_sizes <- setNames(
  c(rep(0.15, n_bg_ov), 0.4, 2.2),
  level_ov
)
all_ov_alpha <- setNames(
  c(rep(0.90, n_bg_ov), 0.40, 0.95),
  level_ov
)

p_spatial_overlay <- ggplot(
  coords_ov,
  aes(x = y_centroid, y = x_centroid,
      color = display_group,
      size  = display_group,
      alpha = display_group)
) +
  geom_point(stroke = 0) +
  scale_color_manual(values = all_ov_cols, name = NULL) +
  scale_size_manual(values  = all_ov_sizes,  guide = "none") +
  scale_alpha_manual(values = all_ov_alpha,  guide = "none") +
  coord_fixed() +
  guides(color = guide_legend(
    override.aes = list(size = 3.5, alpha = 1),
    ncol = 1
  )) +
  labs(title = "Macrophage GDF15 status in full tissue context") +
  theme_pub +
  theme(
    axis.line  = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )
ggsave(file.path(output_dir, "08_spatial_GDF15_overlay_context.png"),
       p_spatial_overlay, width = 10, height = 12, dpi = 300, bg = "white")
print(p_spatial_overlay)

# ── 8c. Macrophage subclusters on tissue (cell-type background overlay) ────────
# Reuse the bg coordinate frame from 8b_overlay; substitute macrophage rows
# with their subcluster label so they render on top in Tableau colours.
coords_8c <- xenium@meta.data[, c("x_centroid", "y_centroid", "auto_celltype")]
coords_8c$display_group <- as.character(coords_8c$auto_celltype)

mac_clusters_8c <- setNames(as.character(macro$seurat_clusters), Cells(macro))
mac_rows_8c <- rownames(coords_8c) %in% Cells(macro)
coords_8c$display_group[mac_rows_8c] <-
  paste0("Cluster ", mac_clusters_8c[rownames(coords_8c)[mac_rows_8c]])

cl_labels_8c  <- sort(unique(coords_8c$display_group[mac_rows_8c]))
bg_labels_8c  <- sort(setdiff(unique(coords_8c$display_group), cl_labels_8c))
level_8c      <- c(bg_labels_8c, cl_labels_8c)
coords_8c$display_group <- factor(coords_8c$display_group, levels = level_8c)
coords_8c <- coords_8c[order(coords_8c$display_group), ]

n_bg_8c    <- length(bg_labels_8c)
n_cl_8c    <- length(cl_labels_8c)
bg_cols_8c <- setNames(scales::hue_pal()(n_bg_8c), bg_labels_8c)
cl_cols_8c <- setNames(tableau10[seq_len(n_cl_8c)], cl_labels_8c)
all_cols_8c  <- c(bg_cols_8c, cl_cols_8c)
all_sizes_8c <- setNames(c(rep(0.15, n_bg_8c), rep(1.5, n_cl_8c)), level_8c)
all_alpha_8c <- setNames(c(rep(0.90, n_bg_8c), rep(0.90, n_cl_8c)), level_8c)

p_spatial_mac_subclusters <- ggplot(
  coords_8c,
  aes(x = y_centroid, y = x_centroid,
      color = display_group,
      size  = display_group,
      alpha = display_group)
) +
  geom_point(stroke = 0) +
  scale_color_manual(values = all_cols_8c, name = NULL) +
  scale_size_manual(values  = all_sizes_8c, guide = "none") +
  scale_alpha_manual(values = all_alpha_8c, guide = "none") +
  coord_fixed() +
  guides(color = guide_legend(override.aes = list(size = 3.5, alpha = 1), ncol = 1)) +
  labs(title = "Macrophage subclusters on tissue") +
  theme_pub +
  theme(
    axis.line  = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )
ggsave(file.path(output_dir, "08_spatial_mac_subclusters.png"),
       p_spatial_mac_subclusters, width = 10, height = 12, dpi = 300, bg = "white")
print(p_spatial_mac_subclusters)

# ── 8d. GDF15-high macrophages highlighted in full tissue context ──────────────
# Same cell-type background; GDF15-low macs small/faded, GDF15-high macs prominent.
coords_8d <- xenium@meta.data[, c("x_centroid", "y_centroid", "auto_celltype")]
coords_8d$display_group <- as.character(coords_8d$auto_celltype)

mac_gdf15_8d <- setNames(macro$GDF15_status, Cells(macro))
mac_rows_8d  <- rownames(coords_8d) %in% Cells(macro)
coords_8d$display_group[mac_rows_8d] <-
  ifelse(mac_gdf15_8d[rownames(coords_8d)[mac_rows_8d]] == "GDF15-high",
         "Macrophage (GDF15-high)", "Macrophage (GDF15-low/neg)")

mac_grps_8d <- c("Macrophage (GDF15-low/neg)", "Macrophage (GDF15-high)")
bg_grps_8d  <- sort(setdiff(unique(coords_8d$display_group), mac_grps_8d))
level_8d    <- c(bg_grps_8d, mac_grps_8d)
coords_8d$display_group <- factor(coords_8d$display_group, levels = level_8d)
coords_8d <- coords_8d[order(coords_8d$display_group), ]

n_bg_8d     <- length(bg_grps_8d)
bg_cols_8d  <- setNames(rep("#C8C8C8", n_bg_8d), bg_grps_8d)   # all background = grey
mac_cols_8d <- c("Macrophage (GDF15-low/neg)" = "#4E9AC7",
                 "Macrophage (GDF15-high)"    = "#D62728")
all_cols_8d  <- c(bg_cols_8d, mac_cols_8d)
# Background: tiny grey. GDF15-low/neg: clearly smaller than red but solid & visible.
# GDF15-high: large red. Order in level_8d is c(background..., low/neg, high).
all_sizes_8d <- setNames(c(rep(0.15, n_bg_8d), 1.3,  2.2),  level_8d)
all_alpha_8d <- setNames(c(rep(0.90, n_bg_8d), 0.90, 0.95), level_8d)

p_spatial_context <- ggplot(
  coords_8d,
  aes(x = y_centroid, y = x_centroid,
      color = display_group,
      size  = display_group,
      alpha = display_group)
) +
  geom_point(stroke = 0) +
  scale_color_manual(values = all_cols_8d, breaks = mac_grps_8d, name = NULL) +
  scale_size_manual(values  = all_sizes_8d, guide = "none") +
  scale_alpha_manual(values = all_alpha_8d, guide = "none") +
  coord_fixed() +
  guides(color = guide_legend(override.aes = list(size = 3.5, alpha = 1), ncol = 1)) +
  labs(title = "GDF15-high macrophages in full tissue context") +
  theme_pub +
  theme(
    axis.line  = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )
ggsave(file.path(output_dir, "08_spatial_GDF15_tissue_context.png"),
       p_spatial_context, width = 10, height = 12, dpi = 300, bg = "white")
print(p_spatial_context)


# =============================================================================
# SECTION 10 — SPATIAL NEIGHBORHOOD ANALYSIS
# GDF15-high vs GDF15-low macrophage microenvironments
# =============================================================================
message("Spatial neighborhood analysis: GDF15-high vs GDF15-low macrophages ...")

RADIUS_UM <- 100  # search radius in micrometers (adjust if tissue scale differs)

# Coordinates and cell types for the full post-QC xenium object
all_coords <- data.frame(
  x = xenium@meta.data$x_centroid,
  y = xenium@meta.data$y_centroid,
  row.names = colnames(xenium)
)
all_celltypes <- setNames(xenium$auto_celltype, colnames(xenium))

# GDF15-high and GDF15-low macrophage barcodes
gdf15_high_cells <- WhichCells(macro, expression = GDF15_status == "GDF15-high")
gdf15_low_cells  <- WhichCells(macro, expression = GDF15_status == "GDF15-low/neg")

# Match sample sizes — use all GDF15-high cells, sample equal from GDF15-low
n_compare <- min(length(gdf15_high_cells), length(gdf15_low_cells))
set.seed(42)
gdf15_high_sample <- sample(gdf15_high_cells, n_compare)
gdf15_low_sample  <- sample(gdf15_low_cells,  n_compare)
cat(sprintf("Neighborhood comparison: %d GDF15-high vs %d GDF15-low macrophages\n",
            n_compare, n_compare))

# Return cell types of every non-focal cell within radius_um of any focal cell
get_neighborhood_types <- function(focal_barcodes, all_coords, all_celltypes, radius_um) {
  fx <- all_coords[focal_barcodes, "x"]
  fy <- all_coords[focal_barcodes, "y"]
  ax <- all_coords$x
  ay <- all_coords$y
  r2 <- radius_um^2
  nbr_types <- lapply(seq_along(focal_barcodes), function(i) {
    d2 <- (ax - fx[i])^2 + (ay - fy[i])^2
    all_celltypes[d2 > 0 & d2 <= r2]
  })
  unlist(nbr_types)
}

nbr_high <- get_neighborhood_types(gdf15_high_sample, all_coords, all_celltypes, RADIUS_UM)
nbr_low  <- get_neighborhood_types(gdf15_low_sample,  all_coords, all_celltypes, RADIUS_UM)

prop_high <- prop.table(table(nbr_high))
prop_low  <- prop.table(table(nbr_low))

all_types <- union(names(prop_high), names(prop_low))
nbr_df <- data.frame(
  celltype   = all_types,
  GDF15_high = as.numeric(prop_high[all_types]),
  GDF15_low  = as.numeric(prop_low[all_types])
)
nbr_df[is.na(nbr_df)] <- 0
nbr_df$diff <- nbr_df$GDF15_high - nbr_df$GDF15_low

write.csv(nbr_df, file.path(output_dir, "10_neighborhood_composition.csv"), row.names = FALSE)

# ── 10a. Lollipop: Δ proportion (GDF15-high minus GDF15-low) ──────────────────
nbr_df$celltype_f <- factor(nbr_df$celltype, levels = nbr_df$celltype[order(nbr_df$diff)])

p_nbr_diff <- ggplot(nbr_df, aes(x = diff, y = celltype_f, color = diff > 0)) +
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = "dashed") +
  geom_segment(aes(x = 0, xend = diff, yend = celltype_f), linewidth = 1.2) +
  geom_point(size = 4) +
  scale_color_manual(
    values = c("TRUE" = "#D62728", "FALSE" = "#4E9AC7"),
    labels = c("TRUE" = "Enriched near GDF15-high", "FALSE" = "Enriched near GDF15-low")
  ) +
  labs(
    title = sprintf("Neighborhood enrichment: GDF15-high vs GDF15-low\n(within %d µm, n = %d per group)",
                    RADIUS_UM, n_compare),
    x     = "Δ proportion  (GDF15-high − GDF15-low)",
    y     = NULL,
    color = NULL
  ) +
  theme_pub
ggsave(file.path(output_dir, "10_neighborhood_enrichment_lollipop.png"), p_nbr_diff,
       width = 9, height = 6, dpi = 300, bg = "white")
print(p_nbr_diff)


# =============================================================================
# SECTION 10b — TGFB1 / TGFB2 / TGFB3 IN LUMINAL EPITHELIAL CELLS
# Epithelial cells within 100 µm of GDF15-high vs GDF15-low macrophages
# =============================================================================
message("TGFB1/2/3 expression in luminal epithelial cells by macrophage neighborhood ...")

luminal_cells <- WhichCells(xenium, expression = auto_celltype == "Epithelial_Luminal")
cat(sprintf("Luminal epithelial cells: %d\n", length(luminal_cells)))
lum_coords <- all_coords[luminal_cells, ]

get_epi_neighbors <- function(focal_barcodes, lum_coords, all_coords, radius_um) {
  r2 <- radius_um^2
  lx <- lum_coords$x;  ly <- lum_coords$y
  unique(unlist(lapply(focal_barcodes, function(bc) {
    d2 <- (lx - all_coords[bc, "x"])^2 + (ly - all_coords[bc, "y"])^2
    rownames(lum_coords)[d2 > 0 & d2 <= r2]
  })))
}

epi_near_high <- get_epi_neighbors(gdf15_high_sample, lum_coords, all_coords, RADIUS_UM)
epi_near_low  <- get_epi_neighbors(gdf15_low_sample,  lum_coords, all_coords, RADIUS_UM)
epi_overlap       <- intersect(epi_near_high, epi_near_low)
epi_near_high_unq <- setdiff(epi_near_high, epi_overlap)
epi_near_low_unq  <- setdiff(epi_near_low,  epi_overlap)

cat(sprintf("Luminal epithelial near GDF15-high only: %d\n", length(epi_near_high_unq)))
cat(sprintf("Luminal epithelial near GDF15-low only : %d\n", length(epi_near_low_unq)))
cat(sprintf("Near both (excluded)                   : %d\n", length(epi_overlap)))

xenium_expr <- GetAssayData(xenium, assay = "Xenium", layer = "data")
tgfb_genes  <- intersect(c("TGFB1", "TGFB2", "TGFB3"), rownames(xenium_expr))
cat(sprintf("TGFB genes in panel: %s\n", paste(tgfb_genes, collapse = ", ")))

if (length(tgfb_genes) > 0) {
  # Long-format data frame: all TGFB genes × both groups
  tgfb_long <- do.call(rbind, lapply(tgfb_genes, function(g) {
    rbind(
      data.frame(gene  = g,
                 expr  = as.numeric(xenium_expr[g, epi_near_high_unq]),
                 group = sprintf("Near GDF15-high\n(n=%d)", length(epi_near_high_unq))),
      data.frame(gene  = g,
                 expr  = as.numeric(xenium_expr[g, epi_near_low_unq]),
                 group = sprintf("Near GDF15-low\n(n=%d)",  length(epi_near_low_unq)))
    )
  }))
  tgfb_long$gene  <- factor(tgfb_long$gene,  levels = tgfb_genes)
  tgfb_long$group <- factor(tgfb_long$group, levels = unique(tgfb_long$group))

  # Wilcoxon test + summary per gene
  wt_list <- lapply(tgfb_genes, function(g) {
    h <- as.numeric(xenium_expr[g, epi_near_high_unq])
    l <- as.numeric(xenium_expr[g, epi_near_low_unq])
    data.frame(gene    = g,
               p_value = wilcox.test(h, l)$p.value,
               mean_high = mean(h), mean_low = mean(l),
               pct_high  = 100 * mean(h > 0),
               pct_low   = 100 * mean(l > 0))
  })
  wt_df <- do.call(rbind, wt_list)
  print(wt_df)
  write.csv(wt_df, file.path(output_dir, "10b_TGFB_expression_stats.csv"),
            row.names = FALSE)

  # Clip to 99th percentile of expressing cells per gene (removes extreme outliers)
  tgfb_plot <- do.call(rbind, lapply(tgfb_genes, function(g) {
    sub <- tgfb_long[tgfb_long$gene == g, ]
    cap <- quantile(sub$expr[sub$expr > 0], 0.99, na.rm = TRUE)
    if (is.na(cap)) cap <- max(sub$expr)
    sub$expr <- pmin(sub$expr, cap)
    sub
  }))

  # Power transformation: bottom 20% of the y-range occupies 50% of visual space.
  # Exponent p = log(0.5) / log(0.2) ≈ 0.431
  # Works per-facet with scales = "free_y" because the transform is scale-invariant.
  SPLIT_P <- log(0.5) / log(0.2)   # ≈ 0.431
  split_trans <- scales::trans_new(
    name      = "split_y",
    transform = function(x) x^SPLIT_P,
    inverse   = function(x) x^(1 / SPLIT_P),
    breaks    = function(x) {
      rng <- range(x, na.rm = TRUE)
      b   <- scales::pretty_breaks(n = 5)(c(0, rng[2]))
      b[b >= 0]
    },
    format = scales::label_number(accuracy = 0.001, trim = TRUE)
  )

  # p-value labels at top; % expressing at the foot of each facet
  y_max_df <- aggregate(expr ~ gene, data = tgfb_plot, FUN = max)
  ann_df   <- merge(wt_df, y_max_df, by = "gene")
  ann_df$label <- sprintf("p = %.3g", ann_df$p_value)
  ann_df$x     <- 1.5
  ann_df$y     <- ann_df$expr

  pct_ann <- rbind(
    data.frame(gene  = wt_df$gene, x = 1,
               label = sprintf("%.1f%% expr", wt_df$pct_high),
               group = levels(tgfb_long$group)[1]),
    data.frame(gene  = wt_df$gene, x = 2,
               label = sprintf("%.1f%% expr", wt_df$pct_low),
               group = levels(tgfb_long$group)[2])
  )
  pct_ann$gene <- factor(pct_ann$gene, levels = tgfb_genes)

  p_tgfb <- ggplot(tgfb_plot, aes(x = group, y = expr, fill = group)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.85) +
    geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA, color = "grey20") +
    geom_text(data = ann_df, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, size = 3.5, color = "grey25",
              vjust = -0.5) +
    geom_text(data = pct_ann, aes(x = x, y = 0, label = label),
              inherit.aes = FALSE, size = 2.9, color = "grey40",
              vjust = 1.6) +
    scale_fill_manual(values = c("#D62728", "#4E9AC7")) +
    scale_y_continuous(
      trans   = split_trans,
      expand  = expansion(mult = c(0.08, 0.18))
    ) +
    facet_wrap(~ gene, scales = "free_y", nrow = 1) +
    labs(
      title    = "TGFB1 / TGFB2 / TGFB3 in luminal epithelial cells",
      subtitle = sprintf(
        "All cells within %d µm of each macrophage group  •  y-axis expanded: bottom 20%% of range = 50%% of visual space  •  %% expr = fraction of cells > 0",
        RADIUS_UM),
      x = NULL, y = "Normalized expression"
    ) +
    theme_pub +
    theme(legend.position = "none",
          strip.text = element_text(face = "bold", size = 11))

  ggsave(file.path(output_dir, "10b_TGFB_epithelial_violin.png"),
         p_tgfb, width = max(8, 4 * length(tgfb_genes)), height = 6,
         dpi = 300, bg = "white")
  print(p_tgfb)
}


# =============================================================================
# SECTION 10c — HALLMARK TGFB SIGNALLING GENES IN EPITHELIAL NEIGHBOURHOODS
# Dot plot + spatial expression maps for Hallmark TGF-Beta Signalling genes
# that are present in the Xenium panel.
# Comparison: epithelial cells exclusively near GDF15-high vs GDF15-low macs.
# =============================================================================
message("Section 10c: Hallmark TGFB signalling genes in epithelial neighbourhoods ...")

library(msigdbr)

# ── Get Hallmark TGF-Beta Signalling gene set ─────────────────────────────────
tgfb_hallmark_raw <- msigdbr(species = "Homo sapiens", category = "H")
tgfb_hallmark_raw <- tgfb_hallmark_raw[
  tgfb_hallmark_raw$gs_name == "HALLMARK_TGF_BETA_SIGNALING", ]
tgfb_hall_genes   <- unique(tgfb_hallmark_raw$gene_symbol)

# Intersect with genes present in the Xenium panel
tgfb_in_panel <- intersect(tgfb_hall_genes, rownames(xenium_expr))
cat(sprintf("Hallmark TGF-Beta Signalling genes in panel: %d / %d\n",
            length(tgfb_in_panel), length(tgfb_hall_genes)))
cat("Genes found:", paste(sort(tgfb_in_panel), collapse = ", "), "\n")

if (length(tgfb_in_panel) == 0) {
  message("No Hallmark TGFB signalling genes found in panel — skipping Section 10c.")
} else {

  # ── Build per-cell expression matrix for the two epithelial groups ───────────
  # epi_near_high_unq and epi_near_low_unq carried forward from Section 10b
  expr_high_epi <- as.matrix(xenium_expr[tgfb_in_panel, epi_near_high_unq, drop = FALSE])
  expr_low_epi  <- as.matrix(xenium_expr[tgfb_in_panel, epi_near_low_unq,  drop = FALSE])

  # ── Dot plot data: % expressing + mean expression per group per gene ──────────
  make_dot_df <- function(expr_mat, group_label) {
    data.frame(
      gene       = rownames(expr_mat),
      group      = group_label,
      pct_expr   = 100 * rowMeans(expr_mat > 0),
      mean_expr  = rowMeans(expr_mat)
    )
  }

  dot_df <- rbind(
    make_dot_df(expr_high_epi,
                sprintf("Near GDF15-high\n(n=%d)", length(epi_near_high_unq))),
    make_dot_df(expr_low_epi,
                sprintf("Near GDF15-low\n(n=%d)",  length(epi_near_low_unq)))
  )

  # Order genes by descending fold-change (mean_high / mean_low)
  fc_df   <- merge(
    dot_df[dot_df$group == dot_df$group[1], c("gene","mean_expr")],
    dot_df[dot_df$group == dot_df$group[nrow(dot_df)], c("gene","mean_expr")],
    by = "gene", suffixes = c("_high", "_low")
  )
  fc_df$logfc    <- log2((fc_df$mean_expr_high + 1e-6) / (fc_df$mean_expr_low + 1e-6))
  gene_order     <- fc_df$gene[order(fc_df$logfc)]
  dot_df$gene    <- factor(dot_df$gene, levels = gene_order)
  dot_df$group   <- factor(dot_df$group, levels = unique(dot_df$group))
  expr_max       <- max(dot_df$mean_expr)

  p_tgfb_dot <- ggplot(dot_df,
    aes(x = group, y = gene,
        size = pct_expr, color = mean_expr)) +
    geom_point() +
    scale_size_continuous(
      name   = "% expressing",
      range  = c(1, 10),
      breaks = c(5, 15, 30, 50),
      limits = c(0, max(dot_df$pct_expr))
    ) +
    scale_color_gradientn(
      colors = c("grey92", "#CDEAF5", "#6BB8D4", "#1D7DB5", "#0B3060"),
      limits = c(0, expr_max),
      name   = "Mean\nexpression"
    ) +
    labs(
      title    = "Hallmark TGF-β Signalling genes in luminal epithelial cells",
      subtitle = sprintf(
        "Cells exclusively within %d µm of each macrophage group  •  ordered by fold-change (high / low)",
        RADIUS_UM),
      x = NULL, y = NULL
    ) +
    theme_pub +
    theme(axis.text.x  = element_text(size = 9.5, hjust = 0.5),
          axis.text.y  = element_text(size = 9),
          panel.grid.major = element_line(color = "grey91", linewidth = 0.35),
          plot.subtitle    = element_text(size = 9, color = "grey42"))

  ggsave(file.path(output_dir, "10c_TGFB_hallmark_dotplot.png"),
         p_tgfb_dot,
         width = max(7, 2.5 + 0.9 * length(tgfb_in_panel)),
         height = 6, dpi = 300, bg = "white")
  print(p_tgfb_dot)

  # ── Compute Hallmark TGF-β module score on the full xenium object ─────────────
  # AddModuleScore aggregates all available pathway genes into a single per-cell
  # score (subtracting a background of similarly expressed control genes).
  # This captures coordinated subtle shifts invisible in individual gene plots.
  message("  Computing Hallmark TGF-β module score ...")
  xenium <- AddModuleScore(xenium,
                           features = list(tgfb_in_panel),
                           name     = "TGFb_score",
                           assay    = "Xenium")
  xenium$TGFb_score  <- xenium$TGFb_score1   # Seurat appends "1"
  xenium$TGFb_score1 <- NULL

  # ── Violin: module score in epithelial cells near GDF15-high vs GDF15-low ─────
  score_high_epi <- xenium$TGFb_score[epi_near_high_unq]
  score_low_epi  <- xenium$TGFb_score[epi_near_low_unq]
  wt_mod <- wilcox.test(score_high_epi, score_low_epi)
  cat(sprintf("TGF-β module score Wilcoxon p = %.4g\n", wt_mod$p.value))
  cat(sprintf("Mean near GDF15-high: %.4f  |  near GDF15-low: %.4f\n",
              mean(score_high_epi), mean(score_low_epi)))

  mod_vln_df <- rbind(
    data.frame(score = score_high_epi,
               group = sprintf("Near GDF15-high\n(n=%d)", length(score_high_epi))),
    data.frame(score = score_low_epi,
               group = sprintf("Near GDF15-low\n(n=%d)",  length(score_low_epi)))
  )
  mod_vln_df$group <- factor(mod_vln_df$group, levels = unique(mod_vln_df$group))
  y_ann <- max(mod_vln_df$score) + diff(range(mod_vln_df$score)) * 0.06

  p_mod_vln <- ggplot(mod_vln_df, aes(x = group, y = score, fill = group)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.85) +
    geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA,
                 color = "grey20", linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey50") +
    annotate("text", x = 1.5, y = y_ann,
             label = sprintf("Wilcoxon p = %.3g", wt_mod$p.value),
             size = 3.8, color = "grey25") +
    scale_fill_manual(values = c("#D62728", "#4E9AC7")) +
    labs(
      title    = "Hallmark TGF-β module score — luminal epithelial cells",
      subtitle = sprintf(
        "Cells exclusively within %d µm of each macrophage group\nPositive score = coordinated TGF-β pathway upregulation",
        RADIUS_UM),
      x = NULL, y = "TGF-β module score"
    ) +
    theme_pub +
    theme(legend.position = "none",
          plot.subtitle = element_text(size = 8.5, color = "grey42"))

  ggsave(file.path(output_dir, "10c_TGFb_module_score_violin.png"),
         p_mod_vln, width = 6, height = 6, dpi = 300, bg = "white")
  print(p_mod_vln)

  # ── Spatial expression maps — one row per gene, two columns (high | low) ─────
  # Use one representative focal macrophage from each group (random, clock-seeded)
  set.seed(as.integer(Sys.time()) %% 100000L)
  focal_10c_high <- sample(gdf15_high_sample, 1)
  focal_10c_low  <- sample(gdf15_low_sample,  1)
  cat(sprintf("10c focal GDF15-high : %s  (x=%.1f, y=%.1f)\n",
              focal_10c_high,
              all_coords[focal_10c_high, "x"],
              all_coords[focal_10c_high, "y"]))
  cat(sprintf("10c focal GDF15-low  : %s  (x=%.1f, y=%.1f)\n",
              focal_10c_low,
              all_coords[focal_10c_low, "x"],
              all_coords[focal_10c_low, "y"]))

  PAD_10C <- 150   # µm window half-width

  # Collect cells in each window
  get_win_10c <- function(focal_bc, pad) {
    cx <- all_coords[focal_bc, "x"];  cy <- all_coords[focal_bc, "y"]
    meta <- xenium@meta.data
    meta$x <- all_coords[rownames(meta), "x"]
    meta$y <- all_coords[rownames(meta), "y"]
    idx <- abs(meta$x - cx) <= pad & abs(meta$y - cy) <= pad
    df  <- meta[idx, c("x", "y")]
    df[, tgfb_in_panel] <- t(as.matrix(
      xenium_expr[tgfb_in_panel, rownames(df), drop = FALSE]))
    df$cx <- cx;  df$cy <- cy
    df
  }

  win_10c_high <- get_win_10c(focal_10c_high, PAD_10C)
  win_10c_low  <- get_win_10c(focal_10c_low,  PAD_10C)

  # Shared colour limit per gene across both windows
  theme_10c <- theme_cowplot(font_size = 9) +
    theme(plot.title      = element_text(size = 9, face = "bold", hjust = 0.5),
          axis.title      = element_text(size = 8),
          axis.text       = element_text(size = 7),
          legend.title    = element_text(size = 7.5, face = "bold"),
          legend.text     = element_text(size = 7),
          plot.background = element_rect(fill = "white", color = NA))

  make_gene_panel <- function(win, focal_bc, gene, shared_lim, focal_fill, pad) {
    cx    <- win$cx[1];  cy <- win$cy[1]
    theta <- seq(0, 2 * pi, length.out = 400)
    circ  <- data.frame(x = cx + RADIUS_UM * cos(theta),
                         y = cy + RADIUS_UM * sin(theta))
    focal_df <- data.frame(x = cx, y = cy)
    win_s <- win[order(win[[gene]]), ]

    ggplot(win_s, aes(x = x, y = y, color = .data[[gene]])) +
      geom_point(size = 1.4, alpha = 0.88, stroke = 0) +
      geom_path(data = circ, aes(x = x, y = y), inherit.aes = FALSE,
                color = "red", linetype = "dashed", linewidth = 0.9) +
      geom_point(data = focal_df, aes(x = x, y = y), inherit.aes = FALSE,
                 shape = 23, size = 5, fill = focal_fill,
                 color = "black", stroke = 1.3) +
      scale_color_gradientn(
        colors = c("grey92", "#CDEAF5", "#6BB8D4", "#1D7DB5", "#0B3060"),
        limits = shared_lim, name = gene,
        guide  = guide_colorbar(barwidth  = unit(2.8, "cm"),
                                barheight = unit(0.32, "cm"),
                                title.position = "top", title.hjust = 0.5)) +
      coord_fixed(xlim = c(cx - pad, cx + pad),
                  ylim = c(cy - pad, cy + pad)) +
      labs(x = "x (µm)", y = "y (µm)") +
      theme_10c + theme(legend.position = "bottom")
  }

  # Build grid: genes × 2 columns
  # Limit to 8 genes max (highest mean expression across groups) to keep readable
  top_genes_10c <- if (length(tgfb_in_panel) > 8) {
    mean_expr_all <- rowMeans(cbind(rowMeans(expr_high_epi), rowMeans(expr_low_epi)))
    names(sort(mean_expr_all, decreasing = TRUE))[1:8]
  } else tgfb_in_panel

  panel_list <- lapply(top_genes_10c, function(g) {
    lim <- c(0, max(c(win_10c_high[[g]], win_10c_low[[g]]), na.rm = TRUE))
    if (lim[2] < 1e-6) lim[2] <- 1
    p_h <- make_gene_panel(win_10c_high, focal_10c_high, g, lim, "#D62728", PAD_10C) +
      labs(title = if (g == top_genes_10c[1])
                    sprintf("GDF15-high focal [%s]", focal_10c_high) else NULL)
    p_l <- make_gene_panel(win_10c_low,  focal_10c_low,  g, lim, "#4E9AC7", PAD_10C) +
      labs(title = if (g == top_genes_10c[1])
                    sprintf("GDF15-low focal [%s]", focal_10c_low) else NULL)
    # Row tag on left
    wrap_plots(p_h, p_l, ncol = 2) +
      plot_annotation(tag_levels   = NULL,
                      title        = g,
                      theme        = theme(
                        plot.title = element_text(face = "bold.italic", size = 10,
                                                  hjust = 0, margin = margin(b = 2))))
  })

  fig_10c <- wrap_plots(panel_list, ncol = 1) +
    plot_annotation(
      title    = sprintf("Hallmark TGF-β Signalling genes — neighbourhood expression (%d µm)",
                         RADIUS_UM),
      subtitle = "Left column: GDF15-high focal macrophage (red ◆)  •  Right column: GDF15-low focal macrophage (blue ◆)\nColour scales shared between left and right for each gene",
      caption  = "Only genes present in the Xenium panel shown",
      theme    = theme(
        plot.title    = element_text(face = "bold", size = 12, hjust = 0),
        plot.subtitle = element_text(size = 8.5, color = "grey42", hjust = 0),
        plot.caption  = element_text(size = 7.5, color = "grey55", hjust = 1)
      )
    )

  fig_10c_height <- 5.5 * length(top_genes_10c)
  ggsave(file.path(output_dir, "10c_TGFB_hallmark_spatial.png"),
         fig_10c, width = 12, height = fig_10c_height, dpi = 300, bg = "white")
  print(fig_10c)

}   # end if (length(tgfb_in_panel) > 0)


# =============================================================================
# SECTION 11 — EPITHELIAL PSEUDO-BULK: PATHWAY LOLLIPOP PLOTS
# Luminal epithelial cells near GDF15-high vs GDF15-absent macrophages
# =============================================================================
message("Pseudo-bulk signaling: epithelial cells near GDF15-high vs GDF15-absent macrophages ...")

library(limma); library(edgeR); library(fgsea); library(msigdbr)

RADIUS_PB     <- 100   # µm
MIN_EPI_CELLS <- 3

# GDF15-high (> 0) vs GDF15-absent (== 0) — clean binary contrast
macro_gdf15_raw <- as.numeric(
  GetAssayData(macro, assay = "Xenium", layer = "data")["GDF15", ])
names(macro_gdf15_raw) <- colnames(macro)
mac_high   <- names(macro_gdf15_raw)[macro_gdf15_raw  > 0]
mac_absent <- names(macro_gdf15_raw)[macro_gdf15_raw == 0]
cat(sprintf("GDF15-expressing macrophages : %d\nGDF15-absent macrophages     : %d\n",
            length(mac_high), length(mac_absent)))

raw_counts <- GetAssayData(xenium, assay = "Xenium", layer = "counts")

build_pseudobulk <- function(focal_mac_bcs, lum_coords, all_coords,
                              raw_counts, radius_um, min_cells) {
  r2 <- radius_um^2;  lx <- lum_coords$x;  ly <- lum_coords$y
  pb_list <- lapply(focal_mac_bcs, function(bc) {
    d2  <- (lx - all_coords[bc, "x"])^2 + (ly - all_coords[bc, "y"])^2
    epi <- rownames(lum_coords)[d2 > 0 & d2 <= r2]
    if (length(epi) < min_cells) return(NULL)
    Matrix::rowSums(raw_counts[, epi, drop = FALSE])
  })
  pb_list <- Filter(Negate(is.null), pb_list)
  if (length(pb_list) < 3)
    stop(sprintf("Only %d valid neighbourhoods — need ≥ 3.", length(pb_list)))
  mat <- do.call(cbind, pb_list)
  colnames(mat) <- paste0("nb_", seq_len(ncol(mat)))
  as.matrix(mat)
}

pb_high   <- build_pseudobulk(mac_high,   lum_coords, all_coords,
                               raw_counts, RADIUS_PB, MIN_EPI_CELLS)
pb_absent <- build_pseudobulk(mac_absent, lum_coords, all_coords,
                               raw_counts, RADIUS_PB, MIN_EPI_CELLS)
cat(sprintf("Pseudo-bulk samples — GDF15-high: %d  |  GDF15-absent: %d\n",
            ncol(pb_high), ncol(pb_absent)))

pb_mat <- cbind(pb_high, pb_absent)
group  <- factor(c(rep("GDF15_high", ncol(pb_high)),
                   rep("GDF15_absent", ncol(pb_absent))),
                 levels = c("GDF15_absent", "GDF15_high"))
dge  <- DGEList(counts = pb_mat, group = group)
keep <- filterByExpr(dge, min.count = 5, min.total.count = 15)
dge  <- dge[keep, , keep.lib.sizes = FALSE]
dge  <- calcNormFactors(dge, method = "TMM")
design <- model.matrix(~ group)
v      <- voom(dge, design, plot = FALSE)
fit    <- lmFit(v, design)
fit    <- eBayes(fit)
de_res      <- topTable(fit, coef = "groupGDF15_high", number = Inf, sort.by = "P")
de_res$gene <- rownames(de_res)
write.csv(de_res, file.path(output_dir, "11_pseudobulk_DE_results.csv"), row.names = FALSE)
cat("\n── Top 20 up (near GDF15-high) ─────\n")
print(head(de_res[de_res$logFC > 0, c("gene","logFC","P.Value","adj.P.Val")], 20))
cat("\n── Top 20 down ─────────────────────\n")
print(head(de_res[de_res$logFC < 0, c("gene","logFC","P.Value","adj.P.Val")], 20))

rank_epi <- setNames(sign(de_res$logFC) * (-log10(de_res$P.Value + 1e-300)),
                     de_res$gene)
rank_epi <- sort(rank_epi, decreasing = TRUE)

# ── Publication theme + lollipop helper (used in Sections 11 & 12) ────────────
theme_pub12 <- theme_classic(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 13, hjust = 0, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 9.5, color = "grey42", hjust = 0, margin = margin(b = 10)),
    axis.title    = element_text(face = "bold", size = 10.5, color = "grey15"),
    axis.text     = element_text(size = 9.5, color = "grey20"),
    axis.line     = element_line(color = "grey30", linewidth = 0.45),
    axis.ticks    = element_line(color = "grey30", linewidth = 0.35),
    legend.title  = element_text(face = "bold", size = 9),
    legend.text   = element_text(size = 8.5),
    legend.key    = element_blank(), legend.background = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    plot.margin      = margin(14, 18, 12, 14)
  )

make_pub_lollipop <- function(fgsea_df, title, subtitle = NULL,
                               padj_cutoff = 1.0, sig_cutoff = 0.05, n_show = 15,
                               col_up = "#B2182B", col_down = "#2166AC",
                               out_png, width = 12, height = 9) {
  sig <- if (padj_cutoff < 1) {
    s <- fgsea_df[!is.na(fgsea_df$padj) & fgsea_df$padj < padj_cutoff, ]
    if (nrow(s) == 0) { message("No terms at padj < ", padj_cutoff); return(invisible(NULL)) }
    s
  } else fgsea_df[!is.na(fgsea_df$NES), ]

  sig$label_tc <- tools::toTitleCase(tolower(sig$pathway_clean))
  sig <- sig[order(sig$NES), ]
  n   <- min(n_show, floor(nrow(sig) / 2))
  top <- unique(rbind(head(sig, n), tail(sig, n)))
  top$direction <- factor(
    ifelse(top$NES > 0, "GDF15-high enriched", "GDF15-absent enriched"),
    levels = c("GDF15-high enriched", "GDF15-absent enriched"))

  # ── Significance flag (BH-adjusted p) + display labels ───────────────────────
  # The plot still shows the top-N terms by NES in each direction (nothing is
  # dropped); terms passing FDR < sig_cutoff are additionally marked with a
  # trailing " *", a black ring on the bubble, and a bold NES value.
  top$significant <- !is.na(top$padj) & top$padj < sig_cutoff
  top$disp <- ifelse(top$significant, paste0(top$label_tc, " *"), top$label_tc)
  top$disp <- factor(top$disp, levels = top$disp[order(top$NES)])
  n_sig    <- sum(top$significant)
  sig_pts  <- top[top$significant, , drop = FALSE]   # ringed overlay layer

  p <- ggplot(top, aes(x = NES, y = disp)) +
    geom_vline(xintercept = 0, linewidth = 0.55, color = "grey35") +
    geom_segment(aes(x = 0, xend = NES, yend = disp, color = direction),
                 linewidth = 1.15, alpha = 0.55, lineend = "round") +
    geom_point(aes(color = direction, size = size, alpha = -log10(pval + 1e-10))) +
    # black ring marks terms significant at FDR < sig_cutoff
    geom_point(data = sig_pts, aes(size = size),
               shape = 21, fill = NA, color = "black", stroke = 1.15,
               show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.2f", NES),
                  hjust = ifelse(NES > 0, -0.22, 1.22)),
              size = 2.8, color = "grey25",
              fontface = ifelse(top$significant, "bold", "plain")) +
    scale_color_manual(values = c("GDF15-high enriched" = col_up,
                                   "GDF15-absent enriched" = col_down),
                        guide = guide_legend(override.aes = list(size = 4, alpha = 1))) +
    scale_size_continuous(name = "Gene set\nsize", range = c(3, 9),
                           guide = guide_legend(override.aes = list(alpha = 0.8))) +
    scale_alpha_continuous(name = expression(-log[10]~italic(p)), range = c(0.35, 1.0)) +
    scale_x_continuous(expand = expansion(mult = c(0.20, 0.20))) +
    labs(title = title, subtitle = subtitle,
         x = "Normalised Enrichment Score (NES)", y = NULL, color = NULL,
         caption = sprintf(
           "Black ring + * = significant at FDR (BH-adjusted p) < %.3g   |   %d of %d shown terms significant",
           sig_cutoff, n_sig, nrow(top))) +
    theme_pub12 +
    theme(axis.text.y = element_text(size = 9.5, color = "grey15"),
          panel.grid.major.x = element_line(color = "grey91", linewidth = 0.35),
          plot.caption = element_text(size = 8, color = "grey40", hjust = 1,
                                       margin = margin(t = 8)),
          legend.position = "right")
  ggsave(out_png, p, width = width, height = height, dpi = 300, bg = "white")
  print(p); invisible(p)
}

# Shared helper: run fgsea, save CSV, call make_pub_lollipop
run_lollipop <- function(gene_sets, rank_vec, set_label, title,
                          out_csv, out_png, n_show = 15, width = 13, height = 10,
                          sig_cutoff = 0.05) {
  set.seed(42)
  res <- fgsea(pathways = gene_sets, stats = rank_vec,
               minSize = 5, maxSize = 500, nPermSimple = 10000)
  res$pathway_clean <- gsub(paste0("^", set_label, "_"), "", res$pathway)
  res$pathway_clean <- gsub("_", " ", res$pathway_clean)
  res_df <- as.data.frame(res)[, c("pathway","pathway_clean","pval","padj","NES","size")]
  write.csv(res_df[order(res_df$padj), ], out_csv, row.names = FALSE)
  n_sig_total <- sum(res_df$padj < sig_cutoff, na.rm = TRUE)
  cat(sprintf("\n── %s: %d term(s) significant at FDR < %.3g;  top 10 by |NES| ──\n",
              set_label, n_sig_total, sig_cutoff))
  print(head(res_df[order(-abs(res_df$NES)), c("pathway_clean","NES","padj")], 10))
  make_pub_lollipop(res_df, title = title,
    subtitle = "Top 15 terms by NES in each direction  •  bubble size = gene set size  •  opacity = −log₁₀(p)  •  ring/* = FDR<0.05",
    padj_cutoff = 1.0, sig_cutoff = sig_cutoff,
    n_show = n_show, out_png = out_png, width = width, height = height)
  invisible(res_df)
}

# Load all gene sets once (reused in Section 12)
message("Loading MSigDB gene sets ...")
hallmark_sets <- split(
  msigdbr(species = "Homo sapiens", category = "H")$gene_symbol,
  msigdbr(species = "Homo sapiens", category = "H")$gs_name)

go_bp_raw  <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
go_bp_sets <- split(go_bp_raw$gene_symbol, go_bp_raw$gs_name)

go_mf_raw  <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF")
go_mf_sets <- split(go_mf_raw$gene_symbol, go_mf_raw$gs_name)

reactome_raw  <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
reactome_sets <- split(reactome_raw$gene_symbol, reactome_raw$gs_name)

# ── 11a. GO Biological Process ────────────────────────────────────────────────
message("  fgsea: GO BP (epithelial) ...")
fgsea_epi_gobp <- run_lollipop(
  go_bp_sets, rank_epi, "GOBP",
  title   = "GO Biological Process — Epithelial Cells Near GDF15-high vs GDF15-absent Macrophages",
  out_csv = file.path(output_dir, "11_epi_fgsea_GOBP.csv"),
  out_png = file.path(output_dir, "11_epi_GOBP_lollipop.png"),
  n_show = 15, width = 13, height = 10)

# ── 11b. GO Molecular Function ────────────────────────────────────────────────
message("  fgsea: GO MF (epithelial) ...")
fgsea_epi_gomf <- run_lollipop(
  go_mf_sets, rank_epi, "GOMF",
  title   = "GO Molecular Function — Epithelial Cells Near GDF15-high vs GDF15-absent Macrophages",
  out_csv = file.path(output_dir, "11_epi_fgsea_GOMF.csv"),
  out_png = file.path(output_dir, "11_epi_GOMF_lollipop.png"),
  n_show = 15, width = 13, height = 9)

# ── 11c. Hallmark pathways ────────────────────────────────────────────────────
message("  fgsea: Hallmark (epithelial) ...")
fgsea_epi_hallmark <- run_lollipop(
  hallmark_sets, rank_epi, "HALLMARK",
  title   = "Hallmark Pathways — Epithelial Cells Near GDF15-high vs GDF15-absent Macrophages",
  out_csv = file.path(output_dir, "11_epi_fgsea_hallmark.csv"),
  out_png = file.path(output_dir, "11_epi_hallmark_lollipop.png"),
  n_show = 15, width = 12, height = 9)

# ── 11d. Reactome pathways ────────────────────────────────────────────────────
message("  fgsea: Reactome (epithelial) ...")
fgsea_epi_reactome <- run_lollipop(
  reactome_sets, rank_epi, "REACTOME",
  title   = "Reactome Pathways — Epithelial Cells Near GDF15-high vs GDF15-absent Macrophages",
  out_csv = file.path(output_dir, "11_epi_fgsea_reactome.csv"),
  out_png = file.path(output_dir, "11_epi_reactome_lollipop.png"),
  n_show = 15, width = 14, height = 10)


# =============================================================================
# SECTION 12 — MACROPHAGE PSEUDO-BULK: PATHWAY LOLLIPOP PLOTS
# GDF15-high vs GDF15-absent macrophages
# =============================================================================
message("Macrophage pseudo-bulk DE + GO/pathway analysis ...")

MIN_MAC_CELLS <- 3
mac_raw_expr <- as.numeric(
  GetAssayData(macro, assay = "Xenium", layer = "data")["GDF15", ])
names(mac_raw_expr) <- colnames(macro)
mac_pb_high   <- names(mac_raw_expr)[mac_raw_expr  > 0]
mac_pb_absent <- names(mac_raw_expr)[mac_raw_expr == 0]
cat(sprintf("GDF15-expressing macrophages : %d  (%.1f%%)\nGDF15-absent macrophages     : %d  (%.1f%%)\n",
            length(mac_pb_high),   100 * length(mac_pb_high)   / ncol(macro),
            length(mac_pb_absent), 100 * length(mac_pb_absent) / ncol(macro)))

macro_raw_counts <- GetAssayData(macro, assay = "Xenium", layer = "counts")
clusters <- as.character(sort(unique(macro$seurat_clusters)))
pb_list <- list(); group_labels <- character(0); cluster_labels <- character(0)

for (cl in clusters) {
  cells_cl  <- WhichCells(macro, idents = cl)
  high_cl   <- intersect(mac_pb_high,   cells_cl)
  absent_cl <- intersect(mac_pb_absent, cells_cl)
  if (length(high_cl) >= MIN_MAC_CELLS) {
    pb_list[[paste0("cl", cl, "_high")]]   <-
      Matrix::rowSums(macro_raw_counts[, high_cl,   drop = FALSE])
    group_labels   <- c(group_labels,   "GDF15_high")
    cluster_labels <- c(cluster_labels, paste0("cl", cl))
  }
  if (length(absent_cl) >= MIN_MAC_CELLS) {
    pb_list[[paste0("cl", cl, "_absent")]] <-
      Matrix::rowSums(macro_raw_counts[, absent_cl, drop = FALSE])
    group_labels   <- c(group_labels,   "GDF15_absent")
    cluster_labels <- c(cluster_labels, paste0("cl", cl))
  }
}
if (length(pb_list) < 4)
  stop(sprintf("Only %d pseudo-bulk samples — need ≥ 4.", length(pb_list)))

pb_mac_mat <- do.call(cbind, lapply(pb_list, as.numeric))
rownames(pb_mac_mat) <- rownames(macro_raw_counts)
colnames(pb_mac_mat) <- names(pb_list)
cat(sprintf("Pseudo-bulk samples: %d  (%d GDF15-high | %d GDF15-absent)\n",
            length(pb_list), sum(group_labels == "GDF15_high"),
            sum(group_labels == "GDF15_absent")))

group_mac   <- factor(group_labels,   levels = c("GDF15_absent", "GDF15_high"))
cluster_fac <- factor(cluster_labels)
dge_mac  <- DGEList(counts = pb_mac_mat, group = group_mac)
keep_mac <- filterByExpr(dge_mac, min.count = 5, min.total.count = 10)
dge_mac  <- dge_mac[keep_mac, , keep.lib.sizes = FALSE]
dge_mac  <- calcNormFactors(dge_mac, method = "TMM")
design_mac <- model.matrix(~ cluster_fac + group_mac)
v_mac      <- voom(dge_mac, design_mac, plot = FALSE)
fit_mac    <- lmFit(v_mac, design_mac)
fit_mac    <- eBayes(fit_mac)
de_mac      <- topTable(fit_mac, coef = "group_macGDF15_high", number = Inf, sort.by = "P")
de_mac$gene <- rownames(de_mac)
write.csv(de_mac, file.path(output_dir, "12_mac_pseudobulk_DE.csv"), row.names = FALSE)
cat("\n── Top 20 up (GDF15-high macrophages) ─────\n")
print(head(de_mac[de_mac$logFC > 0, c("gene","logFC","P.Value","adj.P.Val")], 20))
cat("\n── Top 20 down ─────────────────────────────\n")
print(head(de_mac[de_mac$logFC < 0, c("gene","logFC","P.Value","adj.P.Val")], 20))

rank_mac <- setNames(sign(de_mac$logFC) * (-log10(de_mac$P.Value + 1e-300)), de_mac$gene)
rank_mac <- sort(rank_mac, decreasing = TRUE)

# ── 12a. GO Biological Process ───────────────────────────────────────────────
message("  fgsea: GO BP (macrophage) ...")
fgsea_mac_gobp <- run_lollipop(
  go_bp_sets, rank_mac, "GOBP",
  title   = "GO Biological Process — GDF15-high vs GDF15-absent Macrophages",
  out_csv = file.path(output_dir, "12_mac_fgsea_GOBP.csv"),
  out_png = file.path(output_dir, "12_mac_GOBP_lollipop.png"),
  n_show = 15, width = 13, height = 10)

# ── 12b. GO Molecular Function ───────────────────────────────────────────────
message("  fgsea: GO MF (macrophage) ...")
fgsea_mac_gomf <- run_lollipop(
  go_mf_sets, rank_mac, "GOMF",
  title   = "GO Molecular Function — GDF15-high vs GDF15-absent Macrophages",
  out_csv = file.path(output_dir, "12_mac_fgsea_GOMF.csv"),
  out_png = file.path(output_dir, "12_mac_GOMF_lollipop.png"),
  n_show = 15, width = 13, height = 8)

# ── 12c. Hallmark pathways ────────────────────────────────────────────────────
message("  fgsea: Hallmark (macrophage) ...")
fgsea_mac_hallmark <- run_lollipop(
  hallmark_sets, rank_mac, "HALLMARK",
  title   = "Hallmark Pathways — GDF15-high vs GDF15-absent Macrophages",
  out_csv = file.path(output_dir, "12_mac_fgsea_hallmark.csv"),
  out_png = file.path(output_dir, "12_mac_hallmark_lollipop.png"),
  n_show = 15, width = 12, height = 9)

# ── 12d. Reactome pathways ────────────────────────────────────────────────────
message("  fgsea: Reactome (macrophage) ...")
fgsea_mac_reactome <- run_lollipop(
  reactome_sets, rank_mac, "REACTOME",
  title   = "Reactome Pathways — GDF15-high vs GDF15-absent Macrophages",
  out_csv = file.path(output_dir, "12_mac_fgsea_reactome.csv"),
  out_png = file.path(output_dir, "12_mac_reactome_lollipop.png"),
  n_show = 15, width = 14, height = 10)


# =============================================================================
# SECTION 13 — NEIGHBOURHOOD ZOOM PANELS
# Random GDF15-high and GDF15-low macrophage with 100 µm radius overlay.
# Each run picks new random focal cells (seeded from clock).
# Panels: cell types | GDF15 expression | TGF-β module score — shared colour scales.
# =============================================================================
message("Section 13: neighbourhood zoom panels ...")

# Seed from clock — different macrophages every run
set.seed(as.integer(Sys.time()) %% 100000L)

ZOOM_RADIUS <- 100   # µm — neighbourhood radius
ZOOM_PAD    <- 150   # µm — half-width of plot window
ZOOM_PT     <- 1.8   # background cell dot size
ZOOM_FOCAL  <- 6     # focal macrophage diamond size

# ── Pick one random macrophage from each group ────────────────────────────────
stopifnot("GDF15_status" %in% colnames(macro@meta.data))
zoom_high_bcs <- rownames(macro@meta.data)[macro@meta.data$GDF15_status == "GDF15-high"]
zoom_low_bcs  <- rownames(macro@meta.data)[macro@meta.data$GDF15_status == "GDF15-low/neg"]

focal_high_z <- sample(zoom_high_bcs, 1)
focal_low_z  <- sample(zoom_low_bcs,  1)
cat(sprintf("Zoom focal GDF15-high : %s  (x=%.1f, y=%.1f)\n",
            focal_high_z,
            xenium@meta.data[focal_high_z, "x_centroid"],
            xenium@meta.data[focal_high_z, "y_centroid"]))
cat(sprintf("Zoom focal GDF15-low  : %s  (x=%.1f, y=%.1f)\n",
            focal_low_z,
            xenium@meta.data[focal_low_z,  "x_centroid"],
            xenium@meta.data[focal_low_z,  "y_centroid"]))

# ── Attach expression + module score to metadata frame ────────────────────────
# TGFb_score was computed in Section 10c; reuse it here.
zoom_expr <- GetAssayData(xenium, assay = "Xenium", layer = "data")
zoom_meta             <- xenium@meta.data
zoom_meta$gdf15       <- as.numeric(zoom_expr["GDF15", rownames(zoom_meta)])
zoom_meta$tgfb_module <- xenium$TGFb_score[rownames(zoom_meta)]

# ── Window helper ─────────────────────────────────────────────────────────────
get_zoom_window <- function(focal_bc, meta, pad) {
  cx <- meta[focal_bc, "x_centroid"];  cy <- meta[focal_bc, "y_centroid"]
  meta[abs(meta$x_centroid - cx) <= pad & abs(meta$y_centroid - cy) <= pad, ]
}

win_hz <- get_zoom_window(focal_high_z, zoom_meta, ZOOM_PAD)
win_lz <- get_zoom_window(focal_low_z,  zoom_meta, ZOOM_PAD)

# ── Shared colour limits across both windows ──────────────────────────────────
glim <- c(0, max(c(win_hz$gdf15, win_lz$gdf15), na.rm = TRUE))
if (glim[2] < 1e-6) glim[2] <- 1
# Module score is diverging (can be negative) — use symmetric limits around 0
mod_abs_max <- max(abs(c(win_hz$tgfb_module, win_lz$tgfb_module)), na.rm = TRUE)
if (mod_abs_max < 1e-6) mod_abs_max <- 0.1
tlim <- c(-mod_abs_max, mod_abs_max)

# ── Cell-type palette (all types across both windows) ─────────────────────────
zoom_ct_levs <- sort(unique(c(as.character(win_hz$auto_celltype),
                               as.character(win_lz$auto_celltype))))
zoom_ct_cols <- setNames(scales::hue_pal()(length(zoom_ct_levs)), zoom_ct_levs)

# ── Shared zoom theme ─────────────────────────────────────────────────────────
theme_zoom <- theme_cowplot(font_size = 10) +
  theme(plot.title      = element_text(face = "bold", size = 10, hjust = 0.5),
        axis.title      = element_text(size = 8.5),
        axis.text       = element_text(size = 7.5),
        legend.title    = element_text(size = 8, face = "bold"),
        legend.text     = element_text(size = 7.5),
        legend.key.size = unit(0.45, "cm"),
        plot.background = element_rect(fill = "white", color = NA))

# ── Row-builder ───────────────────────────────────────────────────────────────
make_zoom_row <- function(focal_bc, win, focal_fill,
                           glim, tlim,
                           zoom_ct_cols, radius, pad, pt_sz, focal_sz) {
  cx <- zoom_meta[focal_bc, "x_centroid"]
  cy <- zoom_meta[focal_bc, "y_centroid"]

  theta    <- seq(0, 2 * pi, length.out = 500)
  circ_df  <- data.frame(x = cx + radius * cos(theta),
                          y = cy + radius * sin(theta))
  focal_df <- data.frame(x = cx, y = cy)
  xlims    <- c(cx - pad, cx + pad)
  ylims    <- c(cy - pad, cy + pad)

  win$ct <- factor(as.character(win$auto_celltype), levels = names(zoom_ct_cols))

  # P1 — cell types
  p1 <- ggplot(win, aes(x = x_centroid, y = y_centroid, color = ct)) +
    geom_point(size = pt_sz, alpha = 0.88, stroke = 0) +
    geom_path(data = circ_df, aes(x = x, y = y), inherit.aes = FALSE,
              color = "red", linetype = "dashed", linewidth = 1.0) +
    geom_point(data = focal_df, aes(x = x, y = y), inherit.aes = FALSE,
               shape = 23, size = focal_sz, fill = focal_fill,
               color = "black", stroke = 1.4) +
    scale_color_manual(values = zoom_ct_cols, name = "Cell type", drop = FALSE) +
    coord_fixed(xlim = xlims, ylim = ylims) +
    labs(title = "Cell types", x = "x (µm)", y = "y (µm)") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1)) +
    theme_zoom

  # P2 — GDF15 expression
  p2 <- ggplot(win[order(win$gdf15), ],
               aes(x = x_centroid, y = y_centroid, color = gdf15)) +
    geom_point(size = pt_sz, alpha = 0.90, stroke = 0) +
    geom_path(data = circ_df, aes(x = x, y = y), inherit.aes = FALSE,
              color = "red", linetype = "dashed", linewidth = 1.0) +
    geom_point(data = focal_df, aes(x = x, y = y), inherit.aes = FALSE,
               shape = 23, size = focal_sz, fill = focal_fill,
               color = "black", stroke = 1.4) +
    scale_color_gradientn(
      colors = c("grey92", "#CDEAF5", "#6BB8D4", "#1D7DB5", "#0B3060"),
      limits = glim, name = "GDF15",
      guide  = guide_colorbar(barwidth = unit(3.5, "cm"), barheight = unit(0.35, "cm"),
                               title.position = "top", title.hjust = 0.5)) +
    coord_fixed(xlim = xlims, ylim = ylims) +
    labs(title = "GDF15 expression", x = "x (µm)", y = "y (µm)") +
    theme_zoom + theme(legend.position = "bottom")

  # P3 — TGF-β module score (diverging: blue = low, white = zero, red = high)
  p3 <- ggplot(win[order(win$tgfb_module), ],
               aes(x = x_centroid, y = y_centroid, color = tgfb_module)) +
    geom_point(size = pt_sz, alpha = 0.90, stroke = 0) +
    geom_path(data = circ_df, aes(x = x, y = y), inherit.aes = FALSE,
              color = "red", linetype = "dashed", linewidth = 1.0) +
    geom_point(data = focal_df, aes(x = x, y = y), inherit.aes = FALSE,
               shape = 23, size = focal_sz, fill = focal_fill,
               color = "black", stroke = 1.4) +
    scale_color_gradient2(
      low      = "#2166AC",
      mid      = "grey95",
      high     = "#B2182B",
      midpoint = 0,
      limits   = tlim,
      name     = "TGF-β\nmodule score",
      guide    = guide_colorbar(barwidth = unit(3.5, "cm"), barheight = unit(0.35, "cm"),
                                title.position = "top", title.hjust = 0.5)) +
    coord_fixed(xlim = xlims, ylim = ylims) +
    labs(title = "TGF-β module score", x = "x (µm)", y = "y (µm)") +
    theme_zoom + theme(legend.position = "bottom")

  list(p1 = p1, p2 = p2, p3 = p3)
}

row_hz <- make_zoom_row(focal_high_z, win_hz, "#D62728",
                         glim, tlim,
                         zoom_ct_cols, ZOOM_RADIUS, ZOOM_PAD,
                         ZOOM_PT, ZOOM_FOCAL)
row_lz <- make_zoom_row(focal_low_z,  win_lz, "#4E9AC7",
                         glim, tlim,
                         zoom_ct_cols, ZOOM_RADIUS, ZOOM_PAD,
                         ZOOM_PT, ZOOM_FOCAL)

fig_zoom <- (row_hz$p1 | row_hz$p2 | row_hz$p3) /
            (row_lz$p1  | row_lz$p2  | row_lz$p3) +
  plot_annotation(
    title    = sprintf("Macrophage neighbourhood zoom  (%d µm radius)  —  ◆ = focal macrophage",
                       ZOOM_RADIUS),
    subtitle = sprintf(
      "Top row (red ◆): GDF15-high [%s]\nBottom row (blue ◆): GDF15-low/neg [%s]",
      focal_high_z, focal_low_z),
    caption  = "Colour scales shared across both rows",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13, hjust = 0),
      plot.subtitle = element_text(size = 9.5, color = "grey40", hjust = 0,
                                   margin = margin(b = 6)),
      plot.caption  = element_text(size = 7.5, color = "grey55", hjust = 1)
    )
  )

ts_zoom  <- format(Sys.time(), "%Y%m%d_%H%M%S")
zoom_out <- file.path(output_dir,
                       sprintf("13_neighborhood_zoom_%s.png", ts_zoom))
ggsave(zoom_out, fig_zoom, width = 18, height = 12, dpi = 300, bg = "white")
message(sprintf("Zoom panels saved → %s", zoom_out))
print(fig_zoom)


# =============================================================================
# SECTION 14 — SAVE OBJECTS
# =============================================================================
message("Saving Seurat objects ...")

saveRDS(macro,  file.path(output_dir, "08_macrophage_GDF15_xenium.rds"))
saveRDS(xenium, file.path(output_dir, "08_xenium_full_annotated.rds"))

message("All done! Output files saved to: ", output_dir)
message("
── Files generated ────────────────────────────────────────────────────────────
QC & clustering:
  08_qc_transcript_counts.png  |  08_umap_clusters.png
  08_umap_celltypes.png        |  08_spatial_celltypes.png

Macrophage subclustering:
  08_mac_umap_subclusters.png  |  08_mac_phenotype_heatmap.png

GDF15 by subcluster:
  08_GDF15_violin_subclusters.png  |  08_GDF15_dotplot_subclusters.png
  08_GDF15_featureplot_umap.png

Spatial (Sections 8a–8d + overlay):
  08_spatial_GDF15_macrophages.png
  08_spatial_GDF15_high_low_macrophages.png
  08_spatial_GDF15_overlay_context.png
  08_spatial_mac_subclusters.png
  08_spatial_GDF15_tissue_context.png

Neighbourhood (Section 10):
  10_neighborhood_composition.csv
  10_neighborhood_enrichment_lollipop.png
  10b_TGFB_epithelial_violin.png
  10b_TGFB_expression_stats.csv

Epithelial pseudo-bulk (Section 11):
  11_pseudobulk_DE_results.csv
  11_epi_GOBP_lollipop.png  |  11_epi_GOMF_lollipop.png
  11_epi_hallmark_lollipop.png  |  11_epi_reactome_lollipop.png

Macrophage pseudo-bulk (Section 12):
  12_mac_pseudobulk_DE.csv
  12_mac_GOBP_lollipop.png  |  12_mac_GOMF_lollipop.png
  12_mac_hallmark_lollipop.png  |  12_mac_reactome_lollipop.png

Neighbourhood zoom (Section 13):
  13_neighborhood_zoom_<timestamp>.png

Saved objects:
  08_macrophage_GDF15_xenium.rds
  08_xenium_full_annotated.rds
───────────────────────────────────────────────────────────────────────────────
")

xenium <- readRDS(file.path(output_dir, "08_xenium_full_annotated.rds"))

#SANITY CHECK 
# 1. GDF15+ percentages by cell type
read.csv(file.path(output_dir, "04b_GDF15_celltype_summary.csv"))

# 2. n per group in neighborhood analysis
cat("GDF15-high macrophages: ",
    length(WhichCells(macro, expression = GDF15_status == "GDF15-high")), "\n")
cat("GDF15-low/neg macrophages: ",
    length(WhichCells(macro, expression = GDF15_status == "GDF15-low/neg")), "\n")

# 3. Top neighborhood enrichments for GDF15-high
nc <- read.csv(file.path(output_dir, "10_neighborhood_composition.csv"))
head(nc[order(-nc$diff), ], 5)
head(nc[order( nc$diff), ], 5)

# 4. Hallmark top hits (macrophage and epithelium)
mh <- read.csv(file.path(output_dir, "12_mac_fgsea_hallmark.csv"))
mh[order(mh$padj), ][1:15, c("pathway_clean", "NES", "padj")]

eh <- read.csv(file.path(output_dir, "11_epi_fgsea_hallmark.csv"))
eh[order(eh$padj), ][1:15, c("pathway_clean", "NES", "padj")]

mac_gobp <- read.csv(file.path(output_dir, "12_mac_fgsea_GOBP.csv"))
smad <- mac_gobp[grepl("SMAD", mac_gobp$pathway_clean, ignore.case = TRUE), ]
smad[order(smad$padj), ][1:10, c("pathway_clean", "NES", "padj")]

# Also cytokine / growth-factor / signaling-receptor-regulator activities:
mac_gomf <- read.csv(file.path(output_dir, "12_mac_fgsea_GOMF.csv"))
mac_gomf[order(mac_gomf$padj), ][1:20, c("pathway_clean", "NES", "padj")]

s1_df <- stack(cell_type_markers)
names(s1_df) <- c("gene", "cell_type")
s1_df <- s1_df[, c("cell_type", "gene")]
write.csv(s1_df,
          file.path(output_dir, "S1_celltype_marker_gene_sets.csv"),
          row.names = FALSE)

macro <- readRDS("/Users/alekhyamazumdar/Documents/Xenium/Output/08_macrophage_GDF15_xenium.rds")

# Basic structure
macro
dim(macro)
table(macro$mac_subset)                  # if the subset labels exist
table(macro$GDF15_status)                # if the GDF15 stratification is stored
head(macro@meta.data)

# GDF15 distribution
gdf15 <- GetAssayData(macro, assay = "Xenium", layer = "data")["GDF15", ]
summary(gdf15)
mean(gdf15 > 0)                          # % expressing

