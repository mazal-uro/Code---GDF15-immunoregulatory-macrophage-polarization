# 04_tcell_heterogeneity.R
# Subset CD8+ T cells and characterize Trm subsets.
# Recapitulates the paper's finding that prostate Trm cells are heterogeneous
# and segregate by tissue niche.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

OUT_DIR <- "/Users/alekhyamazumdar/Documents/Xenium/Output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
xen <- readRDS(file.path(OUT_DIR, "03_xenium_niches.rds"))

# ---------- 1. Subset CD8+ T cells ----------
tc <- subset(xen, subset = cell_type == "CD8_Tcell")
cat("CD8+ T cells: ", ncol(tc), "\n")

DefaultAssay(tc) <- "Xenium"
tc <- RunPCA(tc, npcs = 20, verbose = FALSE)
tc <- RunUMAP(tc, dims = 1:20, verbose = FALSE)
tc <- FindNeighbors(tc, dims = 1:20, verbose = FALSE)
tc <- FindClusters(tc, resolution = 0.4, algorithm = 4, verbose = FALSE)

# ---------- 2. Score Trm subset signatures ----------
# Signatures derived from the paper's Trm subset definitions and standard refs.
trm_signatures <- list(
  Trm_core      = c("Cd69", "Itgae", "Runx3", "Prdm1", "Hobit", "Notch1", "Rgs1"),
  Trm_cytotoxic = c("Gzmb", "Gzma", "Gzmk", "Ifng", "Prf1", "Klrg1"),
  Trm_progenitor= c("Tcf7", "Sell", "Ccr7", "Lef1", "Bach2"),
  Trm_exhausted = c("Pdcd1", "Lag3", "Havcr2", "Tox", "Tigit"),
  Trm_tgfb_resp = c("Itgae", "Itgb7", "Smad7", "Skil"),
  Trm_il15_resp = c("Bcl2", "Mcl1", "Eomes", "Bhlhe40")
)

for (n in names(trm_signatures)) {
  feats <- intersect(trm_signatures[[n]], rownames(tc))
  if (length(feats) > 0) {
    tc <- AddModuleScore(tc, features = list(feats), name = paste0("score_", n),
                         assay = "Xenium")
  }
}

# Plot scores on UMAP
score_cols <- grep("^score_", colnames(tc@meta.data), value = TRUE)
p_scores <- FeaturePlot(tc, features = score_cols, ncol = 3, raster = FALSE) &
  scale_color_viridis_c()
ggsave(file.path(OUT_DIR, "tcell_signature_scores.png"), p_scores, width = 12, height = 8)

# ---------- 3. Subset markers ----------
tc_markers <- FindAllMarkers(tc, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.4)
write.csv(tc_markers, file.path(OUT_DIR, "tcell_subcluster_markers.csv"), row.names = FALSE)

# ---------- 4. Subset annotation (edit after inspection) ----------
tc_annotation <- c(
  "1" = "Trm_cytotoxic",
  "2" = "Trm_core",
  "3" = "Trm_progenitor",
  "4" = "Trm_exhausted",
  "5" = "Trm_proliferating"
)
tc$tcell_subset <- ifelse(
  as.character(tc$seurat_clusters) %in% names(tc_annotation),
  tc_annotation[as.character(tc$seurat_clusters)],
  paste0("tcell_cluster_", tc$seurat_clusters)
)

# Push subset labels back to full object
xen$tcell_subset <- NA_character_
xen$tcell_subset[colnames(xen) %in% colnames(tc)] <- tc$tcell_subset[match(
  colnames(xen)[colnames(xen) %in% colnames(tc)], colnames(tc)
)]

# ---------- 5. Distribution of T cell subsets across niches ----------
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
ggsave(file.path(OUT_DIR, "tcell_subset_by_niche.png"), p_dist, width = 8, height = 6)

# ---------- 6. Spatial overlay: where does each subset live? ----------
for (s in unique(xen$sample_id)) {
  p_sp <- ImageDimPlot(
    xen, fov = s, group.by = "tcell_subset",
    cells = WhichCells(xen, expression = !is.na(tcell_subset)),
    size = 0.7, dark.background = FALSE
  ) + ggtitle(paste0("CD8+ T cell subsets — ", s))
  ggsave(file.path(OUT_DIR, paste0("spatial_tcell_subset_", s, ".png")),
         p_sp, width = 8, height = 8)
}

# ---------- 7. Compare 7 dpi vs 30 dpi composition ----------
df_time <- xen@meta.data %>%
  dplyr::filter(!is.na(tcell_subset)) %>%
  dplyr::count(timepoint, tcell_subset) %>%
  dplyr::group_by(timepoint) %>%
  dplyr::mutate(frac = n / sum(n))

p_time <- ggplot(df_time, aes(x = timepoint, y = frac, fill = tcell_subset)) +
  geom_col(position = "stack") +
  labs(x = "Timepoint", y = "Fraction of CD8+ T cells", fill = "Subset") +
  theme_minimal(base_size = 11)
ggsave(file.path(OUT_DIR, "tcell_subset_by_timepoint.png"), p_time, width = 6, height = 6)

saveRDS(tc, file.path(OUT_DIR, "04_tcell_subset.rds"))
saveRDS(xen, file.path(OUT_DIR, "04_xenium_with_tcells.rds"))
message("Saved CD8 subset object and updated full object.")
