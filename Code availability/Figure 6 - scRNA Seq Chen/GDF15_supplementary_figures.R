# =====================================================================
# Supplementary figures: GDF15+ vs GDF15- macrophages
# ---------------------------------------------------------------------
# ALL comparisons gate on GDF15 EXPRESSION (GDF15pos vs GDF15neg), never on
# cluster identity -> cluster/integration-independent, matching the main GSEA.
# Chen (tumor) is the primary context; the enrichment panel also uses Sapiens.
#
# Figures (each printed on its own):
#   A. GDF15 expression distribution + % GDF15+ per dataset (tumor enrichment)
#   B. QC anti-artifact panel: nFeature / nCount, GDF15+ vs GDF15-
#   C. Volcano: GDF15+ vs GDF15- (within Chen)
#   D. Data-driven top DE genes (violins)
#   E. Curated metabolic / immune / SMAD dot plot by GDF15 status
#   F. Module-score box/violins (metabolic, immune, SMAD) with AUC effect sizes
#
# NOTE: curated gene panels are canonical picks -- confirm against the volcano/DE
#       before featuring specific genes. Could not run here. Seurat v5.
# =====================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(SCpubr)
library(patchwork)
slice <- dplyr::slice; filter <- dplyr::filter; mutate <- dplyr::mutate
arrange <- dplyr::arrange; slice_head <- dplyr::slice_head

# -----------------------------------------
# Load + GDF15 status on the FULL object (for the enrichment panel)
# -----------------------------------------
integrated <- readRDS("/Users/alekhyamazumdar/Documents/scRNA/IntegratedMacrophages_v2.rds")
DefaultAssay(integrated) <- "RNA"; integrated <- JoinLayers(integrated)
df_all <- FetchData(integrated, vars = c("GDF15", "dataset"))
df_all$GDF15pos <- df_all$GDF15 > 0

# =====================================================================
# FIGURE A: GDF15 distribution + % GDF15+ per dataset
# =====================================================================
# A1: histogram of GDF15 (all cells) with the >0 gate marked
print(
  ggplot(df_all, aes(GDF15)) +
    geom_histogram(bins = 50, fill = "grey50", colour = "black", linewidth = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "firebrick") +
    theme_minimal(base_size = 14) +
    labs(title = "GDF15 expression distribution (all macrophages)",
         x = "GDF15 (log-normalised)", y = "Cells")
)

# A2: fraction of GDF15+ macrophages per dataset (tumor enrichment)
prop <- df_all %>% group_by(dataset) %>%
  summarise(pct_pos = mean(GDF15pos) * 100, n = n(), .groups = "drop")
print(prop)
print(
  ggplot(prop, aes(dataset, pct_pos, fill = dataset)) +
    geom_col(width = 0.7, colour = "black", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.1f%%", pct_pos)), vjust = -0.4) +
    scale_fill_manual(values = c("Chen" = "#C83658FF", "Sapiens" = "#00925EFF")) +
    theme_minimal(base_size = 14) +
    labs(title = "GDF15+ macrophages are tumor-enriched",
         x = NULL, y = "% GDF15+ macrophages") +
    theme(legend.position = "none")
)

# =====================================================================
# Subset to Chen (tumor) and gate GDF15+ vs GDF15-
# =====================================================================
chen <- subset(integrated, subset = dataset == "Chen"); chen <- JoinLayers(chen)
g <- FetchData(chen, vars = "GDF15", layer = "data")[, 1]
pos <- g > 0; neg <- !pos
chen$GDF15_status <- factor(ifelse(pos, "GDF15pos", "GDF15neg"),
                            levels = c("GDF15neg", "GDF15pos"))
Idents(chen) <- "GDF15_status"
cat(sprintf("Within Chen: %d GDF15+  vs  %d GDF15-\n", sum(pos), sum(neg)))

status_cols <- c("GDF15neg" = "grey70", "GDF15pos" = "#08306B")

auc_of <- function(score) {
  wt <- suppressWarnings(wilcox.test(score[pos], score[neg]))
  list(auc = as.numeric(wt$statistic) / (sum(pos) * sum(neg)), p = wt$p.value)
}

# =====================================================================
# FIGURE B: QC anti-artifact panel (is GDF15+ just low-quality? -> no)
# =====================================================================
print(
  VlnPlot(chen, features = c("nFeature_RNA", "nCount_RNA"),
          group.by = "GDF15_status", cols = status_cols, pt.size = 0, ncol = 2)
)
cat("Median nFeature  GDF15+:", round(median(chen$nFeature_RNA[pos])),
    "| GDF15-:", round(median(chen$nFeature_RNA[neg])), "\n")
cat("Median nCount    GDF15+:", round(median(chen$nCount_RNA[pos])),
    "| GDF15-:", round(median(chen$nCount_RNA[neg])), "\n")

# =====================================================================
# FIGURE C: Volcano GDF15+ vs GDF15- (within Chen)
# =====================================================================
de <- FindMarkers(chen, ident.1 = "GDF15pos", ident.2 = "GDF15neg",
                  assay = "RNA", layer = "data", logfc.threshold = 0, min.pct = 0.1)
de$gene <- rownames(de)
logfc_cutoff <- 0.5; padj_cutoff <- 0.05
de <- de %>% mutate(significance = case_when(
  p_val_adj < padj_cutoff & avg_log2FC >  logfc_cutoff ~ "Up in GDF15+",
  p_val_adj < padj_cutoff & avg_log2FC < -logfc_cutoff ~ "Up in GDF15-",
  TRUE                                                 ~ "NotSig"))
lab_up   <- de %>% filter(significance == "Up in GDF15+") %>% arrange(p_val_adj) %>% slice_head(n = 20)
lab_down <- de %>% filter(significance == "Up in GDF15-") %>% arrange(p_val_adj) %>% slice_head(n = 20)
de <- de %>% mutate(label = ifelse(gene %in% c(lab_up$gene, lab_down$gene, "GDF15"), gene, NA))

print(
  ggplot(de, aes(avg_log2FC, -log10(p_val_adj), color = significance)) +
    geom_point(alpha = 0.7, size = 1.4) +
    scale_color_manual(values = c("Up in GDF15+" = "firebrick",
                                  "Up in GDF15-" = "steelblue", "NotSig" = "grey75")) +
    geom_text_repel(aes(label = label), size = 3.5, max.overlaps = 100, segment.color = "grey80") +
    geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", colour = "grey50") +
    theme_minimal(base_size = 14) +
    labs(title = "GDF15+ vs GDF15- macrophages (within Chen)",
         x = "Log2FC (GDF15+ vs GDF15-)", y = "-Log10 adjusted P", color = NULL)
)

# =====================================================================
# FIGURE D: data-driven top DE genes (violins)
# =====================================================================
top_up   <- de %>% filter(significance == "Up in GDF15+") %>% arrange(p_val_adj) %>% slice_head(n = 6) %>% pull(gene)
top_down <- de %>% filter(significance == "Up in GDF15-") %>% arrange(p_val_adj) %>% slice_head(n = 6) %>% pull(gene)
if (length(top_up))
  print(VlnPlot(chen, features = top_up, group.by = "GDF15_status",
                cols = status_cols, pt.size = 0, ncol = 3))
if (length(top_down))
  print(VlnPlot(chen, features = top_down, group.by = "GDF15_status",
                cols = status_cols, pt.size = 0, ncol = 3))

# =====================================================================
# FIGURE E: curated metabolic / immune / SMAD dot plot by GDF15 status
# =====================================================================
metabolic <- c("APOE","TREM2","SPP1","FABP5","FABP4","LPL","CD36","PPARG","PLIN2","LIPA","MARCO")
immune    <- c("HLA-DRA","HLA-DRB1","CD74","CD86","CD40","IL10","CXCL9","CXCL10","STAT1","IL1B","TNF")
smad      <- c("PMEPA1","SMAD7","SKIL","SMAD6","TGFBI","SMAD3","SMAD4")
themed <- intersect(c(metabolic, immune, smad), rownames(chen))
print(
  DotPlot(chen, features = themed, group.by = "GDF15_status") +
    coord_flip() +
    scale_colour_gradient2(low = "navy", mid = "white", high = "firebrick", midpoint = 0) +
    labs(title = "Metabolic / immune / SMAD genes by GDF15 status",
         x = "Genes", y = NULL) +
    theme(axis.text.y = element_text(size = 8))
)

# =====================================================================
# FIGURE F: module scores with AUC effect sizes
# =====================================================================
set.seed(42)
chen <- AddModuleScore(chen, features = list(intersect(metabolic, rownames(chen))), name = "Metabolic", ctrl = 50)
chen <- AddModuleScore(chen, features = list(intersect(immune,    rownames(chen))), name = "Immune",    ctrl = 50)
chen <- AddModuleScore(chen, features = list(intersect(smad,      rownames(chen))), name = "SMAD",      ctrl = 50)

score_cols <- c(Metabolic = "Metabolic1", Immune = "Immune1", SMAD = "SMAD1")
statdf <- bind_rows(lapply(names(score_cols), function(nm) {
  s <- chen[[score_cols[nm]]][, 1]
  st <- auc_of(s)
  data.frame(status = chen$GDF15_status, score = s,
             panel = sprintf("%s  (AUC=%.2f, p=%.1e)", nm, st$auc, st$p))
}))
print(
  ggplot(statdf, aes(status, score, fill = status)) +
    geom_violin(scale = "width", alpha = 0.6) +
    geom_boxplot(width = 0.15, outlier.size = 0.2) +
    scale_fill_manual(values = status_cols) +
    facet_wrap(~ panel, scales = "free_y") +
    theme_minimal(base_size = 13) +
    labs(title = "Program module scores: GDF15+ vs GDF15- (within Chen)",
         x = NULL, y = "Module score") +
    theme(legend.position = "none")
)

# =====================================================================
# Save the GDF15+/- DE table (the source for the volcano)
# =====================================================================
write.csv(de, "/Users/alekhyamazumdar/Documents/scRNA/DE_GDF15pos_vs_neg_Chen.csv",
          row.names = FALSE)
