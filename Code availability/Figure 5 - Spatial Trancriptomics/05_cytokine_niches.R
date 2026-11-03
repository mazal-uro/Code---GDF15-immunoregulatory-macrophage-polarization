# 05_cytokine_niches.R
# Map IL-15, IL-7, and TGFβ source cells in the prostate and quantify the
# spatial relationship between cytokine sources and CD8+ T cell subsets.
# This recapitulates Figures 3 and 4 of Takehara et al.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(FNN)        # nearest-neighbor search
})

OUT_DIR <- "/Users/alekhyamazumdar/Documents/Xenium/Output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
xen <- readRDS(file.path(OUT_DIR, "04_xenium_with_tcells.rds"))

# ---------- 1. Define cytokine genes on the panel ----------
cytokines <- list(
  IL15  = c("Il15", "Il15ra"),
  IL7   = c("Il7",  "Il7r"),
  TGFb  = c("Tgfb1", "Tgfb2", "Tgfb3", "Tgfbr1", "Tgfbr2")
)
# Keep only genes actually on the panel
cytokines <- lapply(cytokines, function(g) intersect(g, rownames(xen)))
print(cytokines)

# ---------- 2. Spatial expression maps per sample ----------
DefaultAssay(xen) <- "Xenium"
for (s in unique(xen$sample_id)) {
  for (cyt in names(cytokines)) {
    feats <- cytokines[[cyt]]
    if (length(feats) == 0) next
    p_sp <- ImageFeaturePlot(xen, fov = s, features = feats,
                             cols = c("lightgrey", "red"), size = 0.5) +
      plot_annotation(title = paste0(cyt, " — ", s))
    ggsave(file.path(OUT_DIR, paste0("spatial_", cyt, "_", s, ".png")),
           p_sp, width = 4 * length(feats), height = 4)
  }
}

# ---------- 3. Cytokine "source cells": top expressors per gene ----------
# Define source cells as those above the 90th percentile of expression for
# each cytokine, then characterize their cell-type composition.
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
ggsave(file.path(OUT_DIR, "cytokine_source_celltypes.png"), p_src, width = 7, height = 6)

# ---------- 4. CD8+ T cell distance to nearest cytokine source ----------
# For each CD8+ T cell, compute distance to nearest IL-15 / IL-7 / TGFβ source
# cell — operate per FOV (sample) since coords are per FOV.
distance_records <- list()

for (s in unique(xen$sample_id)) {
  sub <- subset(xen, subset = sample_id == s)
  fov <- sub@images[[s]]
  if (is.null(fov)) next
  coords <- GetTissueCoordinates(fov)
  # GetTissueCoordinates returns a data.frame with x, y, cell columns
  rownames(coords) <- coords$cell

  tc_idx  <- sub$cell_type == "CD8_Tcell"
  tc_cells <- colnames(sub)[tc_idx]
  if (length(tc_cells) == 0) next
  tc_xy <- coords[tc_cells, c("x", "y"), drop = FALSE]

  for (cyt in names(cytokines)) {
    src_col <- paste0("source_", cyt)
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
      distance_um  = nn$nn.dist[, 1]   # Xenium coords are in microns
    )
  }
}
dist_df <- bind_rows(distance_records)
write.csv(dist_df, file.path(OUT_DIR, "tcell_cytokine_distances.csv"), row.names = FALSE)

# Distance distributions per T cell subset / cytokine
p_dist <- ggplot(dist_df, aes(x = tcell_subset, y = distance_um, fill = cytokine)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ cytokine, scales = "free_y") +
  scale_y_log10() +
  labs(x = "T cell subset", y = "Distance to nearest source cell (μm, log10)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none")
ggsave(file.path(OUT_DIR, "tcell_cytokine_distances.png"), p_dist, width = 10, height = 5)

# ---------- 5. Niche-level cytokine exposure ----------
# Mean expression of each cytokine within each niche label.
niche_cyt <- xen@meta.data %>%
  dplyr::filter(!is.na(niche_label)) %>%
  dplyr::select(niche_label, starts_with("score_")) %>%
  dplyr::group_by(niche_label) %>%
  dplyr::summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")
write.csv(niche_cyt, file.path(OUT_DIR, "niche_cytokine_scores.csv"), row.names = FALSE)

saveRDS(xen, file.path(OUT_DIR, "05_xenium_cytokine.rds"))
message("Saved cytokine-niche analysis.")
