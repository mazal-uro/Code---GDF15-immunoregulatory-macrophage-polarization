# 07_macrophage_analysis.R
# Subset and characterize macrophage populations in the prostate Xenium data.
# Map M1/M2-like polarization, tissue-resident vs. infiltrating subsets, and
# their spatial relationship with CD8+ T cells and the epithelium.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(FNN)
})

OUT_DIR <- "/Users/alekhyamazumdar/Documents/Xenium/Output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
xen <- readRDS(file.path(OUT_DIR, "06_xenium_neighborhood.rds"))

# ---------- 1. Subset macrophages ----------
mac <- subset(xen, subset = cell_type == "Macrophage")
cat("Macrophages: ", ncol(mac), "\n")

DefaultAssay(mac) <- "Xenium"
mac <- RunPCA(mac, npcs = 20, verbose = FALSE)
mac <- RunUMAP(mac, dims = 1:20, verbose = FALSE)
mac <- FindNeighbors(mac, dims = 1:20, verbose = FALSE)
mac <- FindClusters(mac, resolution = 0.4, algorithm = 4, verbose = FALSE)

# ---------- 2. Macrophage subset signatures ----------
mac_sigs <- list(
  TissueResident = c("Lyve1", "Mrc1", "Folr2", "Cbr2", "Cd163", "Timd4"),
  Infiltrating   = c("Ccr2", "Ly6c2", "Plac8", "Ifitm6"),
  M1_like        = c("Nos2", "Cd86", "Tnf",  "Il1b",  "Cxcl9", "Cxcl10", "Stat1"),
  M2_like        = c("Arg1", "Mrc1", "Cd163", "Retnla", "Chil3", "Mgl2"),
  Phagocytic     = c("Mertk", "Axl", "Tyro3", "Gas6"),
  IFN_responsive = c("Isg15", "Ifit1", "Ifit3", "Irf7", "Stat1", "Mx1"),
  Lipid_assoc    = c("Trem2", "Lpl",  "Spp1", "Apoe", "Gpnmb", "Fabp5")
)

for (n in names(mac_sigs)) {
  feats <- intersect(mac_sigs[[n]], rownames(mac))
  if (length(feats) > 0) {
    mac <- AddModuleScore(mac, features = list(feats), name = paste0("mac_", n),
                          assay = "Xenium")
  }
}

p_scores <- FeaturePlot(mac, features = grep("^mac_", colnames(mac@meta.data), value = TRUE),
                        ncol = 3, raster = FALSE) & scale_color_viridis_c()
ggsave(file.path(OUT_DIR, "macrophage_signature_scores.png"), p_scores, width = 12, height = 8)

# ---------- 3. Subset markers + annotation ----------
mac_markers <- FindAllMarkers(mac, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.4)
write.csv(mac_markers, file.path(OUT_DIR, "macrophage_subcluster_markers.csv"),
          row.names = FALSE)

mac_annotation <- c(
  "1" = "TissueResident_M2",
  "2" = "Infiltrating_M1",
  "3" = "IFN_responsive",
  "4" = "Lipid_associated",
  "5" = "Phagocytic"
)
mac$mac_subset <- ifelse(
  as.character(mac$seurat_clusters) %in% names(mac_annotation),
  mac_annotation[as.character(mac$seurat_clusters)],
  paste0("mac_cluster_", mac$seurat_clusters)
)

# Push back to full object
xen$mac_subset <- NA_character_
xen$mac_subset[colnames(xen) %in% colnames(mac)] <- mac$mac_subset[match(
  colnames(xen)[colnames(xen) %in% colnames(mac)], colnames(mac)
)]

# ---------- 4. Macrophage subsets across niches and timepoints ----------
df_mac <- xen@meta.data %>%
  dplyr::filter(!is.na(mac_subset), !is.na(niche_label)) %>%
  dplyr::count(timepoint, niche_label, mac_subset) %>%
  dplyr::group_by(timepoint, niche_label) %>%
  dplyr::mutate(frac = n / sum(n))

p_mac_niche <- ggplot(df_mac,
                      aes(x = niche_label, y = frac, fill = mac_subset)) +
  geom_col(position = "stack") +
  facet_wrap(~ timepoint) +
  labs(x = "Niche", y = "Fraction of macrophages", fill = "Macrophage subset") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUT_DIR, "macrophage_by_niche_time.png"), p_mac_niche, width = 11, height = 5)

# ---------- 5. Spatial overlay ----------
for (s in unique(xen$sample_id)) {
  p_sp <- ImageDimPlot(
    xen, fov = s, group.by = "mac_subset",
    cells = WhichCells(xen, expression = !is.na(mac_subset)),
    size = 0.7, dark.background = FALSE
  ) + ggtitle(paste0("Macrophage subsets — ", s))
  ggsave(file.path(OUT_DIR, paste0("spatial_macrophage_", s, ".png")),
         p_sp, width = 8, height = 8)
}

# ---------- 6. Macrophage <-> CD8 spatial coupling ----------
# Distance from each CD8+ T cell to its nearest macrophage of each subset.
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

p_md <- ggplot(mac_dist, aes(x = mac_subset, y = distance_um, fill = tcell_subset)) +
  geom_boxplot(outlier.size = 0.3) +
  facet_wrap(~ timepoint) +
  scale_y_log10() +
  labs(x = "Macrophage subset",
       y = "CD8+ T cell distance to nearest mac (μm, log10)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUT_DIR, "tcell_macrophage_distances.png"), p_md, width = 11, height = 6)

saveRDS(mac, file.path(OUT_DIR, "07_macrophage_subset.rds"))
saveRDS(xen, file.path(OUT_DIR, "07_xenium_full.rds"))
message("Saved macrophage analysis.")
