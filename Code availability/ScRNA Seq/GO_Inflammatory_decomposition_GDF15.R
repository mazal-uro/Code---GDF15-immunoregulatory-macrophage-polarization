# =====================================================================
# GDF15+ vs GDF15- macrophages (within Chen, CLUSTER-INDEPENDENT)
# ---------------------------------------------------------------------
# ALL comparisons here gate on GDF15 EXPRESSION (GDF15pos vs GDF15neg), never on
# cluster identity -- so NOTHING here depends on the integration or clustering.
# This one script now contains, all on the same GDF15+ vs GDF15- contrast:
#   (1) Inflammatory-response GO decomposition (module score + AUC + GSEA)
#   (2) SMAD-signalling GO decomposition (module score + AUC + GSEA)
#   (3) A genome-wide GO GSEA lollipop (Biological Process + Molecular Function)
# GDF15 is a member of the SMAD terms and is STRIPPED from all scoring sets to
# avoid circularity.
# ---------------------------------------------------------------------
# Your spatial figure had INFLAMMATORY
# RESPONSE as the top GDF15-ABSENT term -> hypothesis: DOWN in GDF15+.
# The full GO:0006954 term is ~700 genes (mostly not macrophage-relevant), so
# a single module score is diluted. We therefore pull the term AND its
# positive/negative-regulation children and decompose:
#   positive reg -> pro-inflammatory drivers
#   negative reg -> anti-inflammatory / resolving regulators
# This tells us WHETHER GDF15+ cells are less inflammatory because they lack
# drivers or because they express resolvers.
#
# Gene sets pulled reproducibly with msigdbr (exact MSigDB membership) so we
# don't transcribe 700 genes by hand.
#
# Effect size AUC = P(GDF15+ cell scores higher than GDF15- cell) (<0.5 = lower).
# NOTE: could not run here. Seurat v5.
# =====================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(SCpubr)
library(patchwork)
library(clusterProfiler)
library(enrichplot)
library(msigdbr)     # install.packages("msigdbr") if needed
library(org.Hs.eg.db)   # for gseGO (genome-wide GO lollipop)
slice <- dplyr::slice; filter <- dplyr::filter; mutate <- dplyr::mutate; arrange <- dplyr::arrange

# -----------------------------------------
# Pull the GO terms from MSigDB (reproducible)
# -----------------------------------------
# msigdbr >=10 renamed args to collection/subcollection; category still works in
# many versions. If category errors, use: msigdbr(species="Homo sapiens",
# collection="C5") and filter the same way.
m <- tryCatch(msigdbr(species = "Homo sapiens", category = "C5"),
              error = function(e) msigdbr(species = "Homo sapiens", collection = "C5"))

terms <- c(
  # --- Inflammatory-response family ---
  inflammatory_response = "GOBP_INFLAMMATORY_RESPONSE",                       # GO:0006954
  positive_regulation   = "GOBP_POSITIVE_REGULATION_OF_INFLAMMATORY_RESPONSE",# GO:0050729
  negative_regulation   = "GOBP_NEGATIVE_REGULATION_OF_INFLAMMATORY_RESPONSE",# GO:0050728
  acute_inflammatory    = "GOBP_ACUTE_INFLAMMATORY_RESPONSE",                 # GO:0002526
  # --- SMAD-signalling family (GDF15 is a MEMBER of these -> stripped below) ---
  smad_signal     = "GOBP_SMAD_PROTEIN_SIGNAL_TRANSDUCTION",                     # GO:0060395
  smad_regulation = "GOBP_REGULATION_OF_SMAD_PROTEIN_SIGNAL_TRANSDUCTION",       # GO:0060390
  smad_positive   = "GOBP_POSITIVE_REGULATION_OF_SMAD_PROTEIN_SIGNAL_TRANSDUCTION", # GO:0060391
  smad_negative   = "GOBP_NEGATIVE_REGULATION_OF_SMAD_PROTEIN_SIGNAL_TRANSDUCTION"  # GO:0060392
)
gene_sets <- lapply(terms, function(t) unique(m$gene_symbol[m$gs_name == t]))
gene_sets <- lapply(gene_sets, function(g) setdiff(g, "GDF15"))   # safety (not a member, but strip)
for (nm in names(gene_sets))
  cat(sprintf("%-22s (%s): %d genes\n", nm, terms[nm], length(gene_sets[[nm]])))

# -----------------------------------------
# Load object + module-score each term
# -----------------------------------------
integrated <- readRDS("/Users/alekhyamazumdar/Documents/scRNA/IntegratedMacrophages_v2.rds")
DefaultAssay(integrated) <- "RNA"; integrated <- JoinLayers(integrated)

set.seed(42)
score_col <- setNames(character(0), character(0))
for (nm in names(gene_sets)) {
  present <- intersect(gene_sets[[nm]], rownames(integrated))
  integrated <- AddModuleScore(integrated, features = list(present),
                               name = paste0("INF_", nm), ctrl = 50)
  score_col[nm] <- paste0("INF_", nm, "1")
  cat(sprintf("  %-22s: %d/%d present\n", nm, length(present), length(gene_sets[[nm]])))
}

# -----------------------------------------
# GDF15+ vs GDF15- WITHIN Chen
# -----------------------------------------
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
cat("\n=== GDF15+ vs GDF15- (within Chen): module-score AUC (<0.5 = lower in GDF15+) ===\n")
print(res, row.names = FALSE)

# -----------------------------------------
# GO IDs + readable labels (shown in every plot so terms aren't confused)
# -----------------------------------------
go_id <- c(
  inflammatory_response = "GO:0006954",
  positive_regulation   = "GO:0050729",
  negative_regulation   = "GO:0050728",
  acute_inflammatory    = "GO:0002526",
  smad_signal     = "GO:0060395",
  smad_regulation = "GO:0060390",
  smad_positive   = "GO:0060391",
  smad_negative   = "GO:0060392"
)
pretty <- c(
  inflammatory_response = "Inflammatory response (GO:0006954)",
  positive_regulation   = "Positive reg. of inflammatory response (GO:0050729)",
  negative_regulation   = "Negative reg. of inflammatory response (GO:0050728)",
  acute_inflammatory    = "Acute inflammatory response (GO:0002526)",
  smad_signal     = "SMAD signal transduction (GO:0060395)",
  smad_regulation = "Regulation of SMAD signalling (GO:0060390)",
  smad_positive   = "Positive reg. of SMAD signalling (GO:0060391)",
  smad_negative   = "Negative reg. of SMAD signalling (GO:0060392)"
)

# -----------------------------------------
# GSEA against the GDF15+ vs GDF15- ranking
# Run EACH term on its own -> the big GO:0006954 set is guaranteed a result
# (the combined call dropped it), and each plot title carries its GO ID.
# -----------------------------------------
Idents(chen) <- "GDF15_status"
de <- FindMarkers(chen, ident.1 = "GDF15pos", ident.2 = "GDF15neg",
                  assay = "RNA", layer = "data", logfc.threshold = 0, min.pct = 0.05)
de$gene <- rownames(de)
ranks <- sort(setNames(de$avg_log2FC, de$gene), decreasing = TRUE)

gsea_list <- list()
cat("\n=== GSEA per term (negative NES = down in GDF15+) ===\n")
for (nm in names(gene_sets)) {
  t2g1 <- data.frame(term = nm, gene = intersect(gene_sets[[nm]], names(ranks)))
  set.seed(42)
  g1 <- tryCatch(suppressWarnings(clusterProfiler::GSEA(
           ranks, TERM2GENE = t2g1, minGSSize = 5, maxGSSize = 5000,
           pvalueCutoff = 1, verbose = FALSE)), error = function(e) NULL)
  gsea_list[[nm]] <- g1
  d <- if (!is.null(g1)) as.data.frame(g1) else data.frame()
  if (nrow(d))
    cat(sprintf("  %-22s %s | setSize %3d | NES %6.3f | p.adj %.3g\n",
                nm, go_id[nm], d$setSize[1], d$NES[1], d$p.adjust[1]))
  else
    cat(sprintf("  %-22s %s | no GSEA result returned\n", nm, go_id[nm]))
}

# -----------------------------------------
# fgsea directly -> guarantees stats + a plot for the big GO:0006954 term
# (clusterProfiler::GSEA drops near-null large sets; fgsea does not).
# -----------------------------------------
library(fgsea)
pathways <- lapply(gene_sets, function(g) intersect(g, names(ranks)))
set.seed(42)
fg <- fgsea::fgsea(pathways = pathways, stats = ranks, minSize = 5, maxSize = 5000)
cat("\n=== fgsea direct (all terms, including the full GO:0006954) ===\n")
print(fg[, c("pathway","size","NES","pval","padj")])

# -----------------------------------------
# Plots -- GO term named in every title.
# Nice 3-panel gseaplot2 where clusterProfiler returned the term; fall back to
# fgsea::plotEnrichment (always draws the running-ES curve) for any it dropped.
# -----------------------------------------
print(VlnPlot(chen, features = unname(score_col), group.by = "GDF15_status", pt.size = 0, ncol = 4))
for (nm in names(gene_sets)) {
  g1 <- gsea_list[[nm]]
  if (!is.null(g1) && nrow(as.data.frame(g1))) {
    print(enrichplot::gseaplot2(g1, geneSetID = 1, title = pretty[nm]))
  } else {
    st  <- fg[fg$pathway == nm, ]
    sub <- if (nrow(st)) sprintf("fgsea: NES %.2f, padj %.3g", st$NES, st$padj) else ""
    print(fgsea::plotEnrichment(pathways[[nm]], ranks) +
          ggplot2::labs(title = pretty[nm], subtitle = sub))
  }
}

# -----------------------------------------
# Summary CSV (GO IDs; GSEA stats from fgsea so the full term IS included)
# -----------------------------------------
res$GO_ID     <- go_id[res$term]
res$GSEA_NES  <- fg$NES[match(res$term, fg$pathway)]
res$GSEA_padj <- fg$padj[match(res$term, fg$pathway)]
res <- res[, c("term","GO_ID","AUC","p","med_pos","med_neg","GSEA_NES","GSEA_padj")]
print(res, row.names = FALSE)
write.csv(res, "/Users/alekhyamazumdar/Documents/scRNA/Inflammatory_SMAD_GOterm_decomposition_GDF15.csv",
          row.names = FALSE)

# =====================================================================
# Genome-wide GO GSEA lollipop (BP + MF) on the GDF15+ vs GDF15- ranking
# ---------------------------------------------------------------------
# Same lollipop style as the integration script, but here the ranking is the
# GDF15+ vs GDF15- (within Chen) fold change -- CLUSTER-INDEPENDENT. Positive
# NES = up in GDF15+; negative = up in GDF15-.
# =====================================================================
gsea_lollipop <- function(markers, ont = c("BP", "MF"), n_top = 15, title = NULL) {
  ont <- match.arg(ont)
  r <- markers$avg_log2FC; names(r) <- markers$gene
  r <- sort(r[!is.na(names(r)) & !duplicated(names(r))], decreasing = TRUE)
  set.seed(42)
  # pvalueCutoff = 1 -> return ALL tested terms so the top-by-NES terms are shown
  # whether or not they pass FDR; FDR<0.05 terms are then marked (ring + "*").
  eg <- tryCatch(gseGO(geneList = r, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                       ont = ont, minGSSize = 10, maxGSSize = 500,
                       pvalueCutoff = 1, verbose = FALSE),
                 error = function(e) NULL)
  df <- if (!is.null(eg)) as.data.frame(eg) else data.frame()
  if (!nrow(df)) { message("No GO:", ont, " terms"); return(invisible(NULL)) }
  up   <- head(df[order(-df$NES), ], n_top)
  down <- head(df[order(df$NES), ], n_top)
  d <- rbind(up, down); d <- d[!duplicated(d$ID), ]
  d$sig       <- !is.na(d$p.adjust) & d$p.adjust < 0.05        # FDR (BH-adjusted p)
  d$lab       <- ifelse(d$sig, paste0(d$Description, " *"), d$Description)
  d$lab       <- factor(d$lab, levels = d$lab[order(d$NES)])
  d$direction <- ifelse(d$NES > 0, "Up in GDF15+", "Up in GDF15-")
  d$neglogp   <- -log10(pmax(d$pvalue, 1e-300))
  n_sig <- sum(d$sig)
  ggplot(d, aes(NES, lab)) +
    geom_segment(aes(x = 0, xend = NES, yend = lab, colour = direction), linewidth = 0.6) +
    geom_point(aes(colour = direction, size = setSize, alpha = neglogp)) +
    geom_point(data = subset(d, sig), aes(size = setSize),      # black ring on FDR<0.05
               shape = 21, fill = NA, colour = "black", stroke = 1.2, show.legend = FALSE) +
    scale_colour_manual(values = c("Up in GDF15+" = "firebrick", "Up in GDF15-" = "steelblue")) +
    scale_size(range = c(2, 8)) + scale_alpha(range = c(0.3, 1)) +
    geom_vline(xintercept = 0, colour = "grey60") +
    theme_minimal(base_size = 11) +
    labs(title = if (is.null(title)) paste0("GO:", ont) else title,
         subtitle = sprintf("GDF15+ vs GDF15- (within Chen) | top %d/direction | size = set size | opacity = -log10(p) | ring/* = FDR<0.05", n_top),
         x = "Normalised Enrichment Score (NES)", y = NULL,
         colour = NULL, size = "Gene set size", alpha = "-log10(p)",
         caption = sprintf("Black ring + '*' = FDR (BH-adjusted p) < 0.05  |  %d of %d shown terms significant",
                           n_sig, nrow(d)))
}

print(gsea_lollipop(de, ont = "BP", title = "GO: Biological Process (GDF15+ vs GDF15-)"))
print(gsea_lollipop(de, ont = "MF", title = "GO: Molecular Function (GDF15+ vs GDF15-)"))

# =====================================================================
# NES summary heatmap (same style as a pan-cancer NES panel)
# ---------------------------------------------------------------------
# One column = this comparison (GDF15+ vs GDF15- within Chen); rows = the
# targeted GO terms. Cell text = NES with significance stars from the GSEA
# adjusted p-value:  *** p.adj<0.001,  ** <0.01,  * <0.05.
# Positive NES (red) = up in GDF15+ ; negative (blue) = up in GDF15-.
# (For a multi-tissue heatmap like the one you showed, run this per dataset and
#  rbind the `res` tables, adding a column for the tissue/sample.)
# =====================================================================
hm <- res
hm <- hm[!is.na(hm$GSEA_NES), ]   # drop terms with no GSEA result (e.g. full GO:0006954)
hm$stars <- as.character(cut(hm$GSEA_padj, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                             labels = c("***", "**", "*", ""), right = TRUE))
hm$stars[is.na(hm$stars)] <- ""
hm$cell  <- ifelse(is.na(hm$GSEA_NES), "NA", sprintf("%.2f%s", hm$GSEA_NES, hm$stars))
hm$label <- pretty[hm$term]
hm       <- hm[order(hm$GSEA_NES), ]                 # sort rows by NES
hm$label <- factor(hm$label, levels = hm$label)
hm$col   <- "GDF15+ vs GDF15-\n(within Chen)"

print(
  ggplot(hm, aes(x = col, y = label, fill = GSEA_NES)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = cell), size = 4) +
    scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#C0384F",
                         midpoint = 0, na.value = "grey85", name = "NES") +
    theme_minimal(base_size = 12) +
    labs(title = "GSEA NES: GDF15+ vs GDF15- (within Chen)",
         subtitle = "Positive NES = up in GDF15+  |  *** p.adj<0.001, ** <0.01, * <0.05",
         x = NULL, y = NULL) +
    theme(axis.text.x = element_text(face = "bold"), panel.grid = element_blank())
)

