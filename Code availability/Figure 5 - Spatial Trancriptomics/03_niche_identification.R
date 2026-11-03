# 03_niche_identification.R
# Identify spatial niches (epithelial, stromal, immune-enriched, etc.) from the
# annotated Xenium object. Two complementary approaches:
#   (A) Seurat's BuildNicheAssay — clusters cells by neighborhood composition
#   (B) Banksy — spatially-aware clustering that mixes expression with
#       neighborhood expression
#
# This recapitulates Figure 5 of Takehara et al.: epithelial vs. stromal niches
# and the cell-type composition of each.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

OUT_DIR <- "/Users/alekhyamazumdar/Documents/Xenium/Output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
xen <- readRDS(file.path(OUT_DIR, "02_xenium_annotated.rds"))

# ---------- A. BuildNicheAssay (per FOV) ----------
# BuildNicheAssay computes, for each cell, the cell-type composition of its
# k nearest spatial neighbors, then clusters cells on those vectors.
# Run per sample (FOV), then concatenate.
niche_results <- list()

for (s in unique(xen$sample_id)) {
  message("Niche assay for ", s)
  sub <- subset(xen, subset = sample_id == s)
  sub <- BuildNicheAssay(
    object   = sub,
    fov      = s,
    group.by = "cell_type",
    niches.k = 5,        # number of niches to find
    neighbors.k = 30     # neighborhood size
  )
  niche_results[[s]] <- sub@meta.data[, c("niches"), drop = FALSE]
}

# Merge niche labels back onto the full object
niche_df <- do.call(rbind, lapply(niche_results, function(x) {
  data.frame(cell = rownames(x), niche = x$niches)
}))
xen$niche <- niche_df$niche[match(colnames(xen), niche_df$cell)]

# ---------- Niche composition heatmap ----------
comp <- xen@meta.data %>%
  dplyr::filter(!is.na(niche)) %>%
  dplyr::count(niche, cell_type) %>%
  dplyr::group_by(niche) %>%
  dplyr::mutate(frac = n / sum(n)) %>%
  dplyr::ungroup()

p_comp <- ggplot(comp, aes(x = factor(niche), y = cell_type, fill = frac)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(x = "Niche", y = "Cell type", fill = "Fraction\nin niche") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank())
ggsave(file.path(OUT_DIR, "niche_composition.png"), p_comp, width = 7, height = 6)

# ---------- Spatial map of niches per sample ----------
for (s in unique(xen$sample_id)) {
  p_sp <- ImageDimPlot(xen, fov = s, group.by = "niche", size = 0.5,
                       dark.background = FALSE) +
    ggtitle(s)
  ggsave(file.path(OUT_DIR, paste0("spatial_niche_", s, ".png")),
         p_sp, width = 8, height = 8)
}

# ---------- Manual interpretation -> readable niche names ----------
# Inspect niche_composition.png and rename. Typical mapping for prostate:
#   epithelial-luminal    : niche enriched for Luminal_epithelial
#   epithelial-basal      : niche enriched for Basal_epithelial
#   stromal-fibroblast    : niche enriched for Fibroblast + Smooth_muscle
#   peri-vascular         : niche enriched for Endothelial
#   immune-enriched       : niche enriched for T cells / Macrophages
niche_names <- c(
  "1" = "epithelial_luminal",
  "2" = "epithelial_basal",
  "3" = "stromal_fibroblast",
  "4" = "peri_vascular",
  "5" = "immune_enriched"
)
xen$niche_label <- ifelse(
  !is.na(xen$niche) & as.character(xen$niche) %in% names(niche_names),
  niche_names[as.character(xen$niche)],
  NA_character_
)

# ---------- B (optional). Banksy spatial-aware clustering ----------
# Uncomment if you want a second view of niches that uses both expression and
# neighborhood expression.
#
# library(Banksy)
# bnk <- xen
# bnk <- RunBanksy(bnk, lambda = 0.2, assay = "SCT", slot = "data",
#                  features = "variable", k_geom = 30, group.by = "sample_id",
#                  split.scale = TRUE)
# bnk <- RunPCA(bnk, assay = "BANKSY", reduction.name = "pca.banksy", npcs = 30)
# bnk <- FindNeighbors(bnk, reduction = "pca.banksy", dims = 1:30)
# bnk <- FindClusters(bnk, resolution = 0.5, algorithm = 4)
# bnk$banksy_cluster <- bnk$seurat_clusters

saveRDS(xen, file.path(OUT_DIR, "03_xenium_niches.rds"))
message("Saved: ", file.path(OUT_DIR, "03_xenium_niches.rds"))
