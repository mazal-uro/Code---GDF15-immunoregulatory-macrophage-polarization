# =====================================================================
# Integrate macrophages: Tabula Sapiens (healthy) + Chen (tumor)
# Using the NEWLY subsetted objects from the two robust QC scripts.
# ---------------------------------------------------------------------
# WHY THIS DIFFERS FROM THE FIRST TRIAL:
#   * First trial used SCTransform integration. SCT models raw UMI counts, but
#     the Chen object's counts are the authors' PRE-NORMALIZED (non-integer)
#     values -> SCT is invalid on it. So we integrate with LogNormalize instead
#     (both objects already carry proper log-normalized "data" layers).
#   * We intersect to COMMON genes before integrating. Chen has no MT genes, so
#     the intersection automatically removes them -> the MT batch artifact from
#     last time is handled structurally, no special-casing.
#   * RPCA anchors: more conservative than CCA, better for two different
#     technologies (atlas vs tumor).
#
# Scope: load -> harmonize genes -> integrate -> cluster -> UMAP + GDF15.
# DE / GO / heatmaps come AFTER you confirm which cluster is GDF15-high.
#
# NOTE: I could not run this (no R / data here). Seurat v5 (also works on v4).
# =====================================================================

# -----------------------------------------
# Libraries
# -----------------------------------------
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(SCpubr)
library(pheatmap)
library(clusterProfiler)   # GO enrichment
library(org.Hs.eg.db)      # human gene annotations
library(enrichplot)
library(tibble)

# Resolve dplyr verb masking by Bioconductor packages (slice/filter/etc.)
slice      <- dplyr::slice
slice_head <- dplyr::slice_head
filter     <- dplyr::filter
mutate     <- dplyr::mutate
arrange    <- dplyr::arrange
select     <- dplyr::select
rename     <- dplyr::rename
count      <- dplyr::count

# SCpubr::do_FeaturePlot errors if max.cutoff exceeds a feature's max value.
safe_cutoff <- function(obj, feature, desired, layer = "data") {
  vals <- FetchData(obj, vars = feature, layer = layer)[[1]]
  min(desired, max(vals, na.rm = TRUE))
}

# =====================================================================
# STEP 1: Load the two newly subsetted macrophage objects
# =====================================================================
sapiens <- readRDS("/Users/alekhyamazumdar/Documents/scRNA/Tabula Sapiens - Prostate/SapiensProstateMacrophages.rds")
chen    <- readRDS("/Users/alekhyamazumdar/Documents/scRNA/Chen et al/Seurat_ChenEtAlmacrophages.rds")

# Clean dataset labels (used for grouping/splitting in plots)
sapiens$dataset <- "Sapiens"
chen$dataset    <- "Chen"

# Work on the RNA assay and make sure it is a single joined layer (Seurat v5)
DefaultAssay(sapiens) <- "RNA"
DefaultAssay(chen)    <- "RNA"
sapiens <- JoinLayers(sapiens)
chen    <- JoinLayers(chen)

# =====================================================================
# STEP 2: Harmonize the gene space (common genes -> auto-removes MT)
# =====================================================================
common_genes <- intersect(rownames(sapiens), rownames(chen))
cat("Sapiens genes:", nrow(sapiens),
    "| Chen genes:", nrow(chen),
    "| common:", length(common_genes), "\n")
cat("MT genes in common set (should be 0):",
    sum(grepl("^MT-", common_genes)), "\n")

sapiens <- subset(sapiens, features = common_genes)
chen    <- subset(chen,    features = common_genes)

# =====================================================================
# STEP 3: LogNormalize + RPCA integration
# =====================================================================
# Both objects already carry log-normalized "data" layers, so we DO NOT
# re-normalize. We just pick HVGs, scale, and PCA each before finding anchors.
obj_list <- list(Sapiens = sapiens, Chen = chen)
obj_list <- lapply(obj_list, function(o) {
  DefaultAssay(o) <- "RNA"
  FindVariableFeatures(o, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
})

features <- SelectIntegrationFeatures(obj_list, nfeatures = 3000)

obj_list <- lapply(obj_list, function(o) {
  o <- ScaleData(o, features = features, verbose = FALSE)
  o <- RunPCA(o, features = features, verbose = FALSE)
  o
})

anchors <- FindIntegrationAnchors(
  object.list          = obj_list,
  anchor.features      = features,
  normalization.method = "LogNormalize",
  reduction            = "rpca",
  dims                 = 1:30
)

integrated <- IntegrateData(
  anchorset            = anchors,
  normalization.method = "LogNormalize",
  dims                 = 1:30
)

# =====================================================================
# STEP 4: Dimensionality reduction + clustering (on the integrated assay)
# =====================================================================
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
ElbowPlot(integrated, ndims = 30)
n_dims <- 20   # adjust to the elbow

integrated <- RunUMAP(integrated, dims = 1:n_dims)
integrated <- FindNeighbors(integrated, dims = 1:n_dims)
integrated <- FindClusters(integrated, resolution = 0.1)   # was 0.2 (fewer subclusters)

# -----------------------------------------
# Remove cluster 4 -- small isolated artifact (was 0 Chen / 34 Sapiens).
# Done here so it is dropped from EVERY downstream plot and analysis, and from
# the saved object. (Document this exclusion in the methods.)
# -----------------------------------------
integrated <- subset(integrated, subset = seurat_clusters != "4")
integrated$seurat_clusters <- droplevels(integrated$seurat_clusters)
Idents(integrated) <- "seurat_clusters"
cat("After removing cluster 4:", ncol(integrated), "cells | clusters:",
    paste(levels(integrated$seurat_clusters), collapse = ", "), "\n")

# Join RNA layers once so downstream expression pulls / DE work on merged data
DefaultAssay(integrated) <- "RNA"
integrated <- JoinLayers(integrated)

# =====================================================================
# STEP 5: MAIN FIGURE PANELS -- the 4 UMAP / GDF15 plots, ALL IN ONE PLACE
# ---------------------------------------------------------------------
#   1) UMAP coloured by cluster
#   2) UMAP coloured by dataset
#   3) GDF15 in Chen (PCa)
#   4) GDF15 in Tabula Sapiens (healthy)
# Each is printed on its own (no combined grid). Diagnostics follow below.
# =====================================================================
DefaultAssay(integrated) <- "RNA"; integrated <- JoinLayers(integrated)
gmax <- safe_cutoff(integrated, "GDF15", 3.0)   # shared GDF15 colour cap for 3 & 4

# One cluster palette (SCpubr-style) used for BOTH the UMAP and the bar plot,
# so the colours match. Extend the vector if you have >8 clusters.
cluster_cols <- c("#C0384F", "#8E8B2E", "#1B9E8A", "#3B2F8A", "#B24FA0",
                  "#E07B39", "#4C6FB1", "#7F7F7F")
names(cluster_cols) <- as.character(seq_along(cluster_cols) - 1)
cluster_cols <- cluster_cols[levels(factor(integrated$seurat_clusters))]

# 1) UMAP coloured by cluster (original SCpubr colours -- unchanged)
p_clusters <- SCpubr::do_DimPlot(integrated, reduction = "umap", group.by = "seurat_clusters",
                                 pt.size = 1.0, label = TRUE,
                                 plot.title = "Integrated macrophage clusters")
print(p_clusters)

# 2) UMAP coloured by dataset
p_dataset <- SCpubr::do_DimPlot(integrated, reduction = "umap", group.by = "dataset",
                                pt.size = 1.0, plot.title = "Integrated macrophages by dataset")
print(p_dataset)

# 3) GDF15 in Chen (PCa) cells only
chen_cells <- subset(integrated, subset = dataset == "Chen")
pb_chen <- SCpubr::do_FeaturePlot(chen_cells, features = "GDF15", pt.size = 1.0,
             min.cutoff = 0, max.cutoff = gmax, plot.title = "GDF15 - Chen (PCa)")
print(pb_chen)

# 4) GDF15 in Tabula Sapiens (healthy) cells only
sap_cells <- subset(integrated, subset = dataset == "Sapiens")
pb_sap <- SCpubr::do_FeaturePlot(sap_cells, features = "GDF15", pt.size = 1.0,
            min.cutoff = 0, max.cutoff = gmax, plot.title = "GDF15 - Tabula Sapiens (healthy)")
print(pb_sap)

# ---- diagnostics (not part of the 4-panel figure) ----
print(SCpubr::do_FeaturePlot(integrated, features = "GDF15", pt.size = 1.0,
        min.cutoff = 0, max.cutoff = gmax, plot.title = "GDF15 on integrated UMAP"))
print(VlnPlot(integrated, features = "GDF15", group.by = "seurat_clusters",
              split.by = "dataset", pt.size = 0))
print(table(integrated$seurat_clusters, integrated$dataset))

# ---- Stacked bar plot: cluster composition within each dataset ----
comp <- as.data.frame(table(cluster = integrated$seurat_clusters,
                            dataset = integrated$dataset))

print(
  ggplot(comp, aes(x = dataset, y = Freq, fill = cluster)) +
    geom_col(position = "fill", width = 0.8, colour = "black", linewidth = 0.5) +  # black border
    scale_fill_manual(values = cluster_cols) +          # same colours as the UMAP
    scale_y_continuous(labels = scales::percent_format()) +
    labs(title = "Cluster composition per dataset",
         x = NULL, y = "Proportion of cells", fill = "Cluster") +
    theme_minimal(base_size = 14) +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6))
)

# =====================================================================
# STEP 6: Save the integrated object
# =====================================================================
saveRDS(integrated, "/Users/alekhyamazumdar/Documents/scRNA/IntegratedMacrophages_v2.rds")

# =====================================================================
# STEP 7: Cluster of interest (the GDF15-high cluster)
# =====================================================================
# At resolution 0.1 (5 clusters), the GDF15-high population is cluster 1
# (confirmed from the GDF15 violin). Everything below reads from this variable.
cluster_of_interest <- "1"

DefaultAssay(integrated) <- "RNA"
integrated <- JoinLayers(integrated)

logfc_cutoff <- 0.5
padj_cutoff  <- 0.05

# =====================================================================
# STEP 8: DE — cluster of interest vs all other clusters (both datasets)
# =====================================================================
cluster_markers <- FindMarkers(
  integrated,
  ident.1         = cluster_of_interest,
  group.by        = "seurat_clusters",
  assay           = "RNA",
  layer           = "data",
  logfc.threshold = 0,
  min.pct         = 0.1
)
cluster_markers$gene <- rownames(cluster_markers)
head(cluster_markers)
# Dataset make-up of the cluster (is it Chen-dominant? expected yes)
print(table(integrated$dataset[integrated$seurat_clusters == cluster_of_interest]))

cluster_markers <- cluster_markers %>%
  mutate(significance = case_when(
    p_val_adj < padj_cutoff & avg_log2FC >  logfc_cutoff ~ "Up",
    p_val_adj < padj_cutoff & avg_log2FC < -logfc_cutoff ~ "Down",
    TRUE                                                 ~ "NotSig"
  ))

top_up   <- cluster_markers %>% filter(significance == "Up")   %>% arrange(p_val_adj) %>% slice_head(n = 30)
top_down <- cluster_markers %>% filter(significance == "Down") %>% arrange(p_val_adj) %>% slice_head(n = 30)
cluster_markers <- cluster_markers %>%
  mutate(label = ifelse(gene %in% c(top_up$gene, top_down$gene, "GDF15"), gene, NA))

# --- Volcano ---
print(
ggplot(cluster_markers, aes(avg_log2FC, -log10(p_val_adj), color = significance)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("Up" = "firebrick", "Down" = "navy", "NotSig" = "grey70")) +
  geom_text_repel(aes(label = label), size = 4, max.overlaps = 100, segment.color = "grey80") +
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey50") +
  theme_minimal(base_size = 14) +
  labs(title = paste0("Cluster ", cluster_of_interest, " vs All: Volcano"),
       x = paste0("Log2FC (Cluster ", cluster_of_interest, " vs Others)"),
       y = "-Log10 Adjusted P-Value", color = "Gene Status") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
)

# =====================================================================
# STEP 9: GO Biological Process (up / down genes of the cluster)
# =====================================================================
up_genes   <- cluster_markers %>% filter(significance == "Up")   %>% pull(gene)
down_genes <- cluster_markers %>% filter(significance == "Down") %>% pull(gene)
cat("Significant Up genes:", length(up_genes), "| Down genes:", length(down_genes), "\n")

safe_entrez <- function(g) if (length(g)) bitr(g, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID else character(0)
up_entrez   <- safe_entrez(up_genes)
down_entrez <- safe_entrez(down_genes)

run_go <- function(entrez) if (length(entrez))
  enrichGO(entrez, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
           pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE) else NULL
ego_up   <- run_go(up_entrez)
ego_down <- run_go(down_entrez)

if (!is.null(ego_up)   && nrow(as.data.frame(ego_up))   > 0)
  print(barplot(ego_up,   showCategory = 15, title = paste0("GO BP up  - Cluster ", cluster_of_interest)))
if (!is.null(ego_down) && nrow(as.data.frame(ego_down)) > 0)
  print(barplot(ego_down, showCategory = 15, title = paste0("GO BP down - Cluster ", cluster_of_interest)))

# =====================================================================
# STEP 10: Heatmap of top up/down genes across clusters
# =====================================================================
heatmap_genes <- unique(c(top_up$gene, top_down$gene))
expr_all      <- GetAssayData(integrated, assay = "RNA", layer = "data")
heatmap_genes <- intersect(heatmap_genes, rownames(expr_all))
expr_mat      <- as.matrix(expr_all[heatmap_genes, ])

annotation_col <- data.frame(
  Cluster = factor(integrated$seurat_clusters),
  Dataset = integrated$dataset
)
rownames(annotation_col) <- colnames(expr_mat)

pheatmap(
  expr_mat, scale = "row",
  cluster_rows = TRUE, cluster_cols = TRUE,
  annotation_col = annotation_col,
  show_rownames = TRUE, show_colnames = FALSE, fontsize_row = 6,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  main = paste0("Top DE genes in Cluster ", cluster_of_interest)
)

# =====================================================================
# STEP 11: Save DE table
# =====================================================================
write.csv(cluster_markers,
          paste0("/Users/alekhyamazumdar/Documents/scRNA/DE_cluster",
                 cluster_of_interest, "_vs_rest.csv"),
          row.names = FALSE)

# =====================================================================
# STEP 12: Reproduce the manuscript figure panels (a-f)
# ---------------------------------------------------------------------
#   a  UMAP by cluster + by dataset       -> already made in STEP 5
#   e  Volcano (cluster of interest vs all)-> already made in STEP 8
#   Below: b (GDF15 split by dataset), c (ridge), d (marker dotplot),
#          f (GO dotplot of downregulated genes)
# =====================================================================
DefaultAssay(integrated) <- "RNA"; integrated <- JoinLayers(integrated)

# ---- Panel b: the 4 UMAP/GDF15 plots now live in STEP 5 (all in one place). ----

# ---- Panel c: GDF15 ridge plot per cluster, points coloured by dataset ----
dfc <- FetchData(integrated, vars = c("GDF15", "seurat_clusters", "dataset"))
if (requireNamespace("ggridges", quietly = TRUE)) {
  print(
    ggplot(dfc, aes(x = GDF15, y = seurat_clusters, fill = seurat_clusters)) +
      ggridges::geom_density_ridges(scale = 1.1, alpha = 0.5, colour = "grey30") +
      geom_jitter(aes(colour = dataset), height = 0.12, size = 0.3, alpha = 0.5) +
      scale_colour_manual(values = c("Chen" = "#C83658FF", "Sapiens" = "#00925EFF")) +
      geom_vline(xintercept = 1.5, linetype = "dashed", colour = "grey40") +
      theme_minimal(base_size = 14) +
      labs(title = "Log-normalised GDF15 expression by cluster",
           x = "GDF15 (log-normalised)", y = "Cluster") +
      theme(legend.position = "right")
  )
} else {
  message("install.packages('ggridges') for the ridge plot; showing violin instead")
  print(VlnPlot(integrated, features = "GDF15", group.by = "seurat_clusters",
                split.by = "dataset", pt.size = 0))
}

# ---- Panel d: macrophage marker dot plot across clusters ----
# Left column of the figure: canonical pan / M1 / M2 markers (+ GDF15)
marker_genes <- c(
  "CD68","CD14","ITGAM","CSF1R","AIF1","FCGR1A","TYROBP","HLA-DRA","APOE",   # pan
  "IL1B","TNF","IL6","CD80","CD86","STAT1","CXCL9","CXCL10","NOS2","CCL2",   # M1
  "MRC1","IL10","TGFBI","CD163","ARG1","FOLR2","PPARG","IRF4","CCL18","NRP1", # M2
  "TGM2","CD200R1","IL1R2","CD274","PDCD1LG2",
  "GDF15")
marker_genes <- intersect(marker_genes, rownames(integrated))
print(
  DotPlot(integrated, features = marker_genes, group.by = "seurat_clusters") +
    coord_flip() +                                  # genes on Y, clusters on X (as in figure)
    scale_colour_gradient2(low = "navy", mid = "white", high = "firebrick", midpoint = 0) +
    labs(title = "Macrophage marker expression by cluster", x = "Genes", y = "Cluster") +
    theme(axis.text.y = element_text(size = 8))
)

# Right column of the figure: data-driven TOP markers per cluster
Idents(integrated) <- "seurat_clusters"
all_markers <- FindAllMarkers(integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
top_each <- all_markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5) %>%
  pull(gene) %>% unique()
print(
  DotPlot(integrated, features = top_each, group.by = "seurat_clusters") +
    coord_flip() +
    scale_colour_gradient2(low = "navy", mid = "white", high = "firebrick", midpoint = 0) +
    labs(title = "Top cluster markers", x = "Genes", y = "Cluster") +
    theme(axis.text.y = element_text(size = 8))
)

# ---- Panel f: GO enrichment ----
# Your figure uses a GSEA lollipop (NES, both directions), not an ORA dotplot.
# That is built in STEP 13 below (Biological Process + Molecular Function).

# =====================================================================
# STEP 13: GO GSEA lollipop plots (Biological Process + Molecular Function)
# ---------------------------------------------------------------------
# Matches the spatial-figure style: GSEA over the cluster-of-interest-vs-rest
# fold-change ranking; top 15 terms per direction; x = NES; point size =
# gene-set size; opacity = -log10(p); red = up in the GDF15 cluster, blue = down.
# (This is a GSEA, distinct from the STEP 9 over-representation test on the
#  up/down gene LISTS -- GSEA uses the whole ranking, matching your figure.)
# =====================================================================
gsea_lollipop <- function(markers, ont = c("BP", "MF"), n_top = 15, title = NULL) {
  ont <- match.arg(ont)
  ranks <- markers$avg_log2FC
  names(ranks) <- markers$gene
  ranks <- sort(ranks[!is.na(names(ranks)) & !duplicated(names(ranks))], decreasing = TRUE)
  set.seed(42)
  # pvalueCutoff = 1 -> show top-by-NES terms whether or not they pass FDR;
  # FDR<0.05 terms are then marked with a black ring + "*".
  eg <- tryCatch(gseGO(geneList = ranks, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                       ont = ont, minGSSize = 10, maxGSSize = 500,
                       pvalueCutoff = 1, verbose = FALSE),
                 error = function(e) NULL)
  df <- if (!is.null(eg)) as.data.frame(eg) else data.frame()
  if (!nrow(df)) { message("No GO:", ont, " terms"); return(invisible(NULL)) }
  up   <- head(df[order(-df$NES), ], n_top)
  down <- head(df[order(df$NES), ], n_top)
  d <- rbind(up, down); d <- d[!duplicated(d$ID), ]
  d$sig       <- !is.na(d$p.adjust) & d$p.adjust < 0.05
  d$lab       <- ifelse(d$sig, paste0(d$Description, " *"), d$Description)
  d$lab       <- factor(d$lab, levels = d$lab[order(d$NES)])
  d$direction <- ifelse(d$NES > 0, "GDF15 cluster (up)", "Rest (down)")
  d$neglogp   <- -log10(pmax(d$pvalue, 1e-300))   # avoid Inf when p == 0
  n_sig <- sum(d$sig)
  ggplot(d, aes(NES, lab)) +
    geom_segment(aes(x = 0, xend = NES, yend = lab, colour = direction), linewidth = 0.6) +
    geom_point(aes(colour = direction, size = setSize, alpha = neglogp)) +
    geom_point(data = subset(d, sig), aes(size = setSize),
               shape = 21, fill = NA, colour = "black", stroke = 1.2, show.legend = FALSE) +
    scale_colour_manual(values = c("GDF15 cluster (up)" = "firebrick",
                                   "Rest (down)" = "steelblue")) +
    scale_size(range = c(2, 8)) + scale_alpha(range = c(0.3, 1)) +
    geom_vline(xintercept = 0, colour = "grey60") +
    theme_minimal(base_size = 11) +
    labs(title = if (is.null(title)) paste0("GO:", ont) else title,
         subtitle = sprintf("Top %d/direction | size = gene-set size | opacity = -log10(p) | ring/* = FDR<0.05", n_top),
         x = "Normalised Enrichment Score (NES)", y = NULL,
         colour = NULL, size = "Gene set size", alpha = "-log10(p)",
         caption = sprintf("Black ring + '*' = FDR (BH-adjusted p) < 0.05  |  %d of %d shown terms significant",
                           n_sig, nrow(d)))
}

print(gsea_lollipop(cluster_markers, ont = "BP", title = "GO: Biological Process"))
print(gsea_lollipop(cluster_markers, ont = "MF", title = "GO: Molecular Function"))

