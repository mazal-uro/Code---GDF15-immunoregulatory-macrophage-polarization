# =====================================================================
# SMAD GO-term decomposition on GDF15+ vs GDF15- macrophages (within Chen)
# ---------------------------------------------------------------------
# Your spatial figure enriched the BROAD SMAD terms (positively in GDF15-high):
#   GO:0060393 regulation of SMAD protein signal transduction   <- FIGURE TERM
#   GO:0060395 SMAD protein signal transduction
# NOT "negative regulation of SMAD". (The only NEGATIVE term enriched was of the
# TGF-beta RECEPTOR pathway, GO:0030512.)
#
# Here we score the ACTUAL GO terms (exact MSigDB membership, GDF15 removed from
# each to avoid circularity) and DECOMPOSE regulation into its positive vs
# negative arms, on the resolution-independent GDF15+ vs GDF15- design, WITHIN
# Chen (tumor) to avoid the dataset confound.
#
#   positive arm  = activating machinery (ligands/receptors/SMAD3-4)
#   negative arm  = inhibitors (SMAD6/7, SKI, PMEPA1, ...)
#   regulation    = parent (both) ; signal_transduction = the process
#
# Effect size AUC = P(GDF15+ cell scores higher than GDF15- cell).
# NOTE: could not run here (no R/data). Seurat v5.
# =====================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(SCpubr)
library(patchwork)
library(clusterProfiler)
library(enrichplot)
slice <- dplyr::slice; filter <- dplyr::filter; mutate <- dplyr::mutate; arrange <- dplyr::arrange

integrated <- readRDS("/Users/alekhyamazumdar/Documents/scRNA/IntegratedMacrophages_v2.rds")
DefaultAssay(integrated) <- "RNA"; integrated <- JoinLayers(integrated)

# =====================================================================
# The four GO terms (exact MSigDB protein-coding membership; MIRs dropped)
# =====================================================================
gene_sets <- list(
  # GO:0060395 SMAD protein signal transduction (the 84-gene process set)
  SMAD_signal = c("SMAD1","SMAD2","SMAD3","SMAD4","SMAD5","SMAD6","SMAD7","SMAD9",
    "TGFBR1","TGFBR2","TGFBR3","BMPR1A","BMPR2","ACVR1","ACVR2A","ACVRL1","ENG",
    "LRP1","LDLRAD4","TGFB1","TGFB2","TGFB3","INHBA","NODAL","AMH","GDF2","GDF5",
    "GDF6","GDF7","GDF11","BMP2","BMP4","BMP5","BMP6","BMP7","BMP10","GREM1","NOG",
    "DKK1","BMPER","TWSG1","CILP","EMILIN1","UCMA","CCN3","SKI","PMEPA1","TOB1",
    "PIN1","VEPH1","EID2","PIAS2","HIPK2","WWTR1","CITED1","RUNX2","TBX20","PPARG",
    "XBP1","FOS","JUN","JAK2","ZMIZ1","ZMIZ2","PML","DAB2","SPTBN1","ATOH8","OVOL2",
    "FAM89B","CRKL","CSNK2B","SH2B1","PARP1","HFE","GLCE","BECN1","MAP1LC3A","NUP93",
    "KIAA0319","MAGI2","RBPMS","PSG9","TTK"),

  # GO:0060393 regulation of SMAD protein signal transduction  (FIGURE TERM)
  regulation = c("ACVR1","ACVR2A","ACVRL1","AMH","ATOH8","BMP10","BMP2","BMP4","BMP5",
    "BMP6","BMP7","BMPER","BMPR1A","BMPR2","CCN3","CILP","CRKL","CSNK2B","DAB2","DKK1",
    "EID2","EMILIN1","ENG","FAM89B","GDF11","GDF2","GDF5","GDF6","GDF7","GLCE","GREM1",
    "HFE","INHBA","JAK2","KIAA0319","LDLRAD4","LRP1","NODAL","NOG","NUP93","OVOL2",
    "PARP1","PIN1","PMEPA1","PPARG","PSG9","SH2B1","SKI","SMAD3","SMAD4","SMAD6","SMAD7",
    "SPTBN1","TBX20","TGFB1","TGFB2","TGFB3","TGFBR1","TGFBR2","TGFBR3","TOB1","TTK",
    "TWSG1","UCMA","VEPH1","XBP1"),

  # GO:0060391 positive regulation (activating machinery)
  positive = c("ACVR1","ACVR2A","ACVRL1","AMH","ATOH8","BMP10","BMP2","BMP4","BMP5",
    "BMP6","BMP7","BMPER","BMPR1A","BMPR2","CSNK2B","DAB2","ENG","GDF11","GDF2","GDF5",
    "GDF6","GDF7","GLCE","HFE","INHBA","JAK2","KIAA0319","NODAL","NUP93","PARP1","PSG9",
    "SH2B1","SMAD3","SMAD4","TGFB1","TGFB2","TGFB3","TGFBR1","TGFBR2","TGFBR3","TTK","TWSG1"),

  # GO:0060392 negative regulation (inhibitors)
  negative = c("CCN3","CILP","CRKL","DKK1","EID2","EMILIN1","ENG","FAM89B","GREM1",
    "LDLRAD4","LRP1","NOG","OVOL2","PIN1","PMEPA1","PPARG","SKI","SMAD6","SMAD7","TBX20",
    "TGFBR3","UCMA","VEPH1","XBP1")
)
# Remove GDF15 from every set (safety; it is a member of several) & keep present genes
gene_sets <- lapply(gene_sets, function(g) setdiff(g, "GDF15"))

# =====================================================================
# Module scores (one column per term)
# =====================================================================
set.seed(42)
score_col <- setNames(character(0), character(0))
for (nm in names(gene_sets)) {
  present <- intersect(gene_sets[[nm]], rownames(integrated))
  integrated <- AddModuleScore(integrated, features = list(present),
                               name = paste0("S_", nm), ctrl = 50)
  score_col[nm] <- paste0("S_", nm, "1")
  cat(sprintf("%-12s: %d/%d genes present\n", nm, length(present), length(gene_sets[[nm]])))
}

# =====================================================================
# GDF15+ vs GDF15- WITHIN Chen
# =====================================================================
chen <- subset(integrated, subset = dataset == "Chen"); chen <- JoinLayers(chen)
g   <- FetchData(chen, vars = "GDF15", layer = "data")[, 1]
pos <- g > 0; neg <- !pos
cat(sprintf("\nWithin Chen: %d GDF15+  vs  %d GDF15-\n", sum(pos), sum(neg)))
chen$GDF15_status <- ifelse(pos, "GDF15pos", "GDF15neg")

auc_test <- function(score) {
  wt  <- suppressWarnings(wilcox.test(score[pos], score[neg]))
  auc <- as.numeric(wt$statistic) / (sum(pos) * sum(neg))
  data.frame(AUC = auc, p = wt$p.value,
             med_pos = median(score[pos]), med_neg = median(score[neg]))
}
res <- bind_rows(lapply(names(score_col), function(nm)
          cbind(term = nm, auc_test(chen[[score_col[nm]]][, 1]))))
cat("\n=== GDF15+ vs GDF15- (within Chen) : module-score AUC per GO term ===\n")
print(res, row.names = FALSE)

# =====================================================================
# GSEA of each GO term against the GDF15+ vs GDF15- ranking
# =====================================================================
Idents(chen) <- "GDF15_status"
de <- FindMarkers(chen, ident.1 = "GDF15pos", ident.2 = "GDF15neg",
                  assay = "RNA", layer = "data", logfc.threshold = 0, min.pct = 0.05)
de$gene <- rownames(de)
ranks <- sort(setNames(de$avg_log2FC, de$gene), decreasing = TRUE)
t2g <- bind_rows(lapply(names(gene_sets), function(nm)
          data.frame(term = nm, gene = intersect(gene_sets[[nm]], names(ranks)))))
set.seed(42)
gs <- suppressWarnings(clusterProfiler::GSEA(ranks, TERM2GENE = t2g,
         minGSSize = 5, maxGSSize = 1000, pvalueCutoff = 1, verbose = FALSE))
cat("\n=== GSEA (positive NES = up in GDF15+) ===\n")
gdf <- as.data.frame(gs)
if (nrow(gdf)) print(gdf[, c("Description","setSize","NES","pvalue","p.adjust")], row.names = FALSE)

# =====================================================================
# Plots
# =====================================================================
VlnPlot(chen, features = unname(score_col), group.by = "GDF15_status", pt.size = 0, ncol = 2)
for (i in seq_len(nrow(gdf)))
  print(enrichplot::gseaplot2(gs, geneSetID = i, title = gdf$Description[i]))

# =====================================================================
# Summary CSV (AUC merged with GSEA)
# =====================================================================
res$GSEA_NES  <- gdf$NES[match(res$term, rownames(gdf))]
res$GSEA_padj <- gdf$p.adjust[match(res$term, rownames(gdf))]
print(res, row.names = FALSE)
write.csv(res, "/Users/alekhyamazumdar/Documents/scRNA/SMAD_GOterm_decomposition_GDF15.csv",
          row.names = FALSE)
