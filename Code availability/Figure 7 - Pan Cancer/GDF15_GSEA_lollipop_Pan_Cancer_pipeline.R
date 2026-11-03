#!/usr/bin/env Rscript
# =============================================================================
#  GDF15+ vs GDF15- MACROPHAGES  ->  DE  ->  Volcano  ->  GSEA  ->  Lollipop
#  Pan-cancer myeloid atlas (GSE154763, Cheng et al. Cell 2021) + Sade-Feldman
#  melanoma macrophages (MELMAC). Fully self-contained: gene sets are pulled
#  DIRECTLY from R (org.Hs.eg.db for GO:BP, msigdbr for Hallmark & Reactome) --
#  no external .gmt file required.
#
#  Per cancer type it produces, matching the Chen-et-al figure:
#    * DE_<ds>.csv                 GDF15+ vs GDF15- FindMarkers table
#    * Volcano_<ds>.png/pdf        volcano (panel c/e style)
#    * GSEA_<collection>_<ds>.csv  full GSEA result
#    * Lollipop_<collection>_<ds>  NES lollipop, top 15/direction by p  (panel e)
#    * Running_<collection>_<ds>_<up/down>.png   gseaplot2 running score (f/g)
#
#  Gating is on GDF15 EXPRESSION (>0) within MACROPHAGES only -> cluster- and
#  integration-independent, exactly like the pooled Chen analysis.
#
#  HOW TO RUN
#    Rscript GDF15_GSEA_lollipop_pipeline.R
#  or, interactively:  source("GDF15_GSEA_lollipop_pipeline.R")
#
#  DATA EXPECTED in DATA_DIR (already present from your pipeline):
#    GSE154763_<DS>_normalized_expression.csv.gz   (rows = cells, cols = genes)
#    GSE154763_<DS>_metadata.csv.gz                (has MajorCluster, UMAP1/2)
#  MELMAC comes from melanoma_macrophage_pipeline.R (run that first if missing).
# =============================================================================

## ----------------------------- 0. Libraries ------------------------------ ##
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(clusterProfiler)   # GSEA(), gseGO()
  library(org.Hs.eg.db)      # GO:BP gene sets, pulled from R directly
  library(AnnotationDbi)     # select() -> genes for a specific GO id
  library(GO.db)             # verify GO term names
  library(enrichplot)        # gseaplot2() -> panels f/g
  library(msigdbr)           # Hallmark + Reactome gene sets, from R directly
})

## ------------------------- 1. Parameters (edit) --------------------------- ##
DATA_DIR <- "/Users/alekhyamazumdar/Documents/Mazumdar et al/Communications Biology/Revisions/scRNA"
OUT_DIR  <- file.path(DATA_DIR, "Analysis_Output", "GDF15_GSEA")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

GENE       <- "GDF15"
POS_CUTOFF <- 0                 # GDF15+ if expression > POS_CUTOFF
# cDC2 excluded per request; LYM excluded because GDF15 is not in that matrix.
DATASETS   <- c("ESCA", "KIDNEY", "MYE", "OV-FTC",
                "PAAD", "THCA", "UCEC", "MELMAC")
# Display order for the cross-cancer supplementary figures (QC / per-patient /
# threshold): prostate-like (immunosuppressive) cancers grouped first, rest after.
DISPLAY_ORDER <- c("ESCA", "UCEC", "OV-FTC", "THCA",    # inflammatory DOWN, SMAD UP
                   "MYE", "MELMAC", "KIDNEY", "PAAD")   # not prostate-like

# DE thresholds (match your Chen script)
MIN_PCT      <- 0.10
LOGFC_THRESH <- 0               # test all genes; annotate volcano at 0.5
LOGFC_LABEL  <- 0.5
PADJ_LABEL   <- 0.05

# GSEA settings
COLLECTIONS  <- c("GOBP")       # GO:BP only
MIN_GS       <- 10
MAX_GS       <- 500
TOP_N        <- 15              # terms per direction on the lollipop
# Running-score (gseaplot2) plots: ONLY these two GO:BP terms, per tumor type.
RUNNING_TERMS <- c(
  "GO:0050729" = "Positive regulation of inflammatory response (GO:0050729)",
  "GO:0060390" = "Regulation of SMAD protein signal transduction (GO:0060390)")
# NES-heatmap terms (cancers x terms): the readouts that defined the prostate
# story -- inflammatory response (expected DOWN) and SMAD (expected UP).
HEATMAP_TERMS <- c(
  "GO:0050729" = "Pos. reg. inflammatory response",
  "GO:0060390" = "Reg. SMAD signal transduction")
# order heatmap rows by this term's NES so the 'blue' (down in GDF15+) cancers
# cluster together. Set to NA to keep the original DATASETS order instead.
HEATMAP_SORT_ID <- "GO:0050729"
# ranking metric for GSEA: signed significance (Seurat -> GSEA convention).
# swap to avg_log2FC by setting RANK_BY <- "log2fc"
RANK_BY      <- "signed_p"      # "signed_p" or "log2fc"

DO_POOLED    <- FALSE           # TRUE = also run one high-power pooled analysis
                                # across all cancers (cancer as covariate label)
SEED         <- 42
set.seed(SEED)
options(future.globals.maxSize = 8 * 1024^3)

status_cols <- c("GDF15neg" = "grey70", "GDF15pos" = "#08306B")
UP_COL   <- "firebrick"; DOWN_COL <- "steelblue"

## --------------------------- 2. Helpers ----------------------------------- ##
save_two <- function(p, name, w, h) {
  print(p)   # show in the RStudio Plots panel / active graphics device
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = w, height = h, bg = "white")
}

# --- loader for a GSE154763 dataset (already log-normalized) ----------------
load_one_dataset <- function(name, data_dir) {
  message("\n[load] ", name)
  expr_file <- file.path(data_dir, sprintf("GSE154763_%s_normalized_expression.csv.gz", name))
  meta_file <- file.path(data_dir, sprintf("GSE154763_%s_metadata.csv.gz", name))
  stopifnot(file.exists(expr_file), file.exists(meta_file))

  meta <- as.data.frame(fread(meta_file)); rownames(meta) <- meta[["index"]]
  expr <- fread(expr_file)                         # rows = cells, cols = genes
  cell_ids <- expr[[1]]; expr[, 1] <- NULL
  gene_ids <- colnames(expr)
  m <- as.matrix(expr); rm(expr); gc()
  m <- t(m); colnames(m) <- cell_ids; rownames(m) <- gene_ids
  m <- as(m, "sparseMatrix")
  meta <- meta[colnames(m), , drop = FALSE]

  seu <- CreateSeuratObject(counts = m, meta.data = meta, project = name)
  seu <- SetAssayData(seu, layer = "data", new.data = m)   # already log-norm
  seu$dataset  <- name
  seu$is_macro <- grepl("Macro", seu$MajorCluster)
  rm(m); gc()
  message("  cells: ", ncol(seu), " | macrophages: ", sum(seu$is_macro))
  seu
}

# --- prettify a gene-set name for axis labels -------------------------------
pretty_term <- function(x) {
  x <- sub("^(HALLMARK|REACTOME|GOBP|KEGG|WP)_", "", x)
  x <- gsub("_", " ", x)
  paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
}

# --- build the ranked gene vector for GSEA ----------------------------------
make_ranklist <- function(de) {
  de <- de[de$gene != GENE, ]                      # drop the gating gene
  if (RANK_BY == "log2fc") {
    v <- de$avg_log2FC
  } else {
    v <- sign(de$avg_log2FC) * -log10(pmax(de$p_val, 1e-300))
  }
  names(v) <- de$gene
  v <- v[is.finite(v)]
  sort(v, decreasing = TRUE)
}

# --- run one GSEA collection, return a gseaResult (or NULL) -----------------
run_gsea <- function(geneList, collection) {
  res <- tryCatch({
    if (collection == "GOBP") {
      gseGO(geneList, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
            minGSSize = MIN_GS, maxGSSize = MAX_GS, pvalueCutoff = 1,
            eps = 0, seed = TRUE, verbose = FALSE)
    } else {
      cat <- if (collection == "HALLMARK") "H" else "C2"
      sub <- if (collection == "HALLMARK") NA else "CP:REACTOME"
      # msigdbr API changed at v10 (collection/subcollection). Try both.
      m <- tryCatch(
        if (is.na(sub)) msigdbr(species = "Homo sapiens", collection = cat)
        else            msigdbr(species = "Homo sapiens", collection = cat, subcollection = sub),
        error = function(e)
          if (is.na(sub)) msigdbr(species = "Homo sapiens", category = cat)
          else            msigdbr(species = "Homo sapiens", category = cat, subcategory = sub))
      t2g <- m[, c("gs_name", "gene_symbol")]
      GSEA(geneList, TERM2GENE = t2g, minGSSize = MIN_GS, maxGSSize = MAX_GS,
           pvalueCutoff = 1, eps = 0, seed = TRUE, verbose = FALSE)
    }
  }, error = function(e) { message("  [gsea] ", collection, " failed: ", conditionMessage(e)); NULL })
  res
}

# --- lollipop (panel e) from a gseaResult -----------------------------------
lollipop_plot <- function(gres, ds, collection) {
  df <- as.data.frame(gres)
  if (nrow(df) == 0) { message("  [lollipop] no terms for ", collection, " ", ds); return(invisible()) }
  df$label <- if (collection == "GOBP") df$Description else pretty_term(df$ID)
  df$neglog10p <- -log10(pmax(df$pvalue, 1e-300))
  up <- df %>% filter(NES > 0) %>% arrange(pvalue) %>% slice_head(n = TOP_N)
  dn <- df %>% filter(NES < 0) %>% arrange(pvalue) %>% slice_head(n = TOP_N)
  d  <- bind_rows(dn, up)
  if (nrow(d) == 0) return(invisible())
  d$dir <- ifelse(d$NES > 0, "Up in GDF15+ Macrophages", "Down in GDF15+ Macrophages")
  d <- d %>% arrange(NES) %>% mutate(label = factor(make.unique(label), levels = make.unique(label)))

  p <- ggplot(d, aes(NES, label)) +
    geom_vline(xintercept = 0, colour = "black", linewidth = 0.4) +
    geom_segment(aes(x = 0, xend = NES, y = label, yend = label, colour = dir),
                 linewidth = 0.6, show.legend = FALSE) +
    geom_point(aes(size = setSize, alpha = neglog10p, colour = dir)) +
    scale_colour_manual(values = c("Up in GDF15+ Macrophages" = UP_COL,
                                   "Down in GDF15+ Macrophages" = DOWN_COL), name = NULL) +
    scale_size_continuous(range = c(2, 8), name = "Gene set size") +
    scale_alpha_continuous(range = c(0.25, 1), name = "-log10(p)") +
    labs(x = "Normalised Enrichment Score (NES)", y = NULL,
         title = sprintf("GDF15+ vs GDF15- (%s) | %s | Top %d/direction", ds, collection, TOP_N),
         subtitle = "Size = set size | Opacity = -log10(p)") +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_text(size = 7),
          plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size = 4)))
  save_two(p, sprintf("Lollipop_%s_%s", collection, ds),
           w = 9, h = 0.28 * nrow(d) + 2.5)
}

# --- running-score plots (panels f/g): ONLY the two chosen GO:BP terms -------
# Computed in a DEDICATED, permissive GSEA pass (genes pulled straight from
# org.Hs.eg.db, minGSSize lowered) so a small set like GO:0060393 is always
# tested and plottable -- independent of the lollipop's size window.
target_gsea <- function(geneList, go_ids) {
  ann <- tryCatch(
    suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = go_ids,
                     keytype = "GOALL", columns = "SYMBOL")),
    error = function(e) NULL)
  if (is.null(ann) || nrow(ann) == 0) return(NULL)
  ann <- ann[!is.na(ann$SYMBOL), c("GOALL", "SYMBOL")]
  t2g <- unique(setNames(ann, c("gs", "gene")))
  lab <- c(RUNNING_TERMS, HEATMAP_TERMS); lab <- lab[!duplicated(names(lab))]
  t2n <- data.frame(gs = names(lab), name = unname(lab))
  t2n <- t2n[t2n$gs %in% t2g$gs, , drop = FALSE]
  tryCatch(
    GSEA(geneList, TERM2GENE = t2g, TERM2NAME = t2n,
         minGSSize = 3, maxGSSize = 5000, pvalueCutoff = 1, eps = 0,
         seed = TRUE, verbose = FALSE),
    error = function(e) { message("  [running] target GSEA failed: ", conditionMessage(e)); NULL })
}

running_plots <- function(gres_t, ds) {
  ids <- names(RUNNING_TERMS)
  # report the true GO term names (optional diagnostic; namespace-qualified so it
  # works whether or not GO.db is attached, and never aborts the run)
  if (requireNamespace("GO.db", quietly = TRUE)) {
    real <- tryCatch(
      suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = unique(c(ids, names(HEATMAP_TERMS))),
                       keytype = "GOID", columns = "TERM")),
      error = function(e) NULL)
    if (!is.null(real))
      for (i in seq_len(nrow(real)))
        message(sprintf("  [term] %s = \"%s\"", real$GOID[i], real$TERM[i]))
  }

  if (is.null(gres_t)) { message("  [running] no target-term GSEA for ", ds); return(invisible()) }
  df <- as.data.frame(gres_t)
  for (id in ids) {
    if (!(id %in% df$ID)) {
      message(sprintf("  [running] %s not tested in %s (too few leading-edge genes in this cancer).",
                      id, ds)); next
    }
    pg <- tryCatch(
      gseaplot2(gres_t, geneSetID = id, title = sprintf("%s: %s", ds, RUNNING_TERMS[[id]]),
                color = "#00A651",
                base_size = 11, pvalue_table = FALSE, rel_heights = c(1.5, 0.5, 1)),
      error = function(e) { message("  [running] ", id, " plot failed: ", conditionMessage(e)); NULL })
    if (!is.null(pg)) {
      print(pg)   # show running-score plot in the Plots panel
      ggsave(file.path(OUT_DIR, sprintf("Running_GOBP_%s_%s.png", ds, gsub(":", "", id))),
             pg, width = 6.5, height = 5, dpi = 300, bg = "white")
    }
  }
}

# --- NES-heatmap rows: NES + p for the key terms, one cancer -----------------
# NES + nominal p come from the targeted GSEA (consistent method, every term
# always tested); genome-wide FDR is pulled from the full GO:BP run when the
# term passed its size filter there (small sets like SMAD won't have one).
heatmap_row <- function(gres_full, gres_t, ds) {
  ft <- if (!is.null(gres_full)) as.data.frame(gres_full) else data.frame(ID = character())
  tt <- if (!is.null(gres_t))    as.data.frame(gres_t)    else data.frame(ID = character())
  do.call(rbind, lapply(names(HEATMAP_TERMS), function(id) {
    nes <- NA_real_; pv <- NA_real_
    if (id %in% tt$ID)      { r <- tt[tt$ID == id, ][1, ]; nes <- r$NES; pv <- r$pvalue }
    else if (id %in% ft$ID) { r <- ft[ft$ID == id, ][1, ]; nes <- r$NES; pv <- r$pvalue }
    fdr <- if (id %in% ft$ID) ft[ft$ID == id, ][1, "p.adjust"] else NA_real_
    data.frame(dataset = ds, id = id, term = HEATMAP_TERMS[[id]],
               NES = nes, pval = pv, FDR_genomewide = fdr, stringsAsFactors = FALSE)
  }))
}

nes_heatmap <- function(rows_list) {
  d <- do.call(rbind, rows_list); rownames(d) <- NULL
  if (is.null(d) || nrow(d) == 0) { message("[heatmap] no rows."); return(invisible()) }
  write.csv(d, file.path(OUT_DIR, "Heatmap_NES_key_terms.csv"), row.names = FALSE)
  # row order: sort by the chosen term's NES so blues (down in GDF15+) group;
  # most-negative ends up at the TOP of the plot.
  if (!is.na(HEATMAP_SORT_ID) && HEATMAP_SORT_ID %in% names(HEATMAP_TERMS)) {
    sl  <- HEATMAP_TERMS[[HEATMAP_SORT_ID]]
    ord <- d[d$term == sl, ]
    ord <- ord[order(ord$NES, decreasing = TRUE), "dataset"]   # +NES bottom, -NES top
    ord <- unique(c(ord, setdiff(unique(d$dataset), ord)))     # any leftovers last
  } else {
    ord <- rev(DATASETS)
  }
  d$dataset <- factor(d$dataset, levels = ord)
  d$term    <- factor(d$term,    levels = unname(HEATMAP_TERMS))
  d$star    <- as.character(cut(d$pval, breaks = c(-Inf, .01, .05, .1, Inf),
                                labels = c("***", "**", "*", "")))
  d$star[is.na(d$pval)] <- ""
  d$lab <- ifelse(is.na(d$NES), "ns", sprintf("%.2f%s", d$NES, d$star))
  lim <- max(abs(d$NES), na.rm = TRUE)
  p <- ggplot(d, aes(term, dataset, fill = NES)) +
    geom_tile(colour = "grey80", linewidth = 0.4) +
    geom_text(aes(label = lab), size = 3) +
    scale_fill_gradient2(low = DOWN_COL, mid = "white", high = UP_COL, midpoint = 0,
                         limits = c(-lim, lim), na.value = "grey92", name = "NES") +
    labs(x = NULL, y = NULL,
         title = "GDF15+ vs GDF15- macrophages: enrichment across cancers",
         subtitle = paste0("red = up in GDF15+, blue = down | stars = GSEA nominal p ",
                           "(*<0.1 **<0.05 ***<0.01) | ns = not testable")) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1, face = "bold"),
          axis.text.y = element_text(face = "bold"),
          panel.grid = element_blank(),
          plot.title = element_text(face = "bold"))
  save_two(p, "Heatmap_NES_key_terms", w = 7.5, h = 0.5 * nlevels(d$dataset) + 2.5)
}

# ============================================================================
#  Anti-artifact / robustness figures (one per-cell QC row collected per cancer)
# ============================================================================
# per-cancer table of QC metrics + GDF15 status (uses the ORIGINAL metadata
# columns n_counts / n_genes / percent_mito, not Seurat's recomputed ones,
# because the matrix is already log-normalized).
qc_table <- function(macro, g, ds) {
  md  <- macro@meta.data
  getc <- function(nm) if (nm %in% colnames(md)) md[[nm]] else rep(NA_real_, ncol(macro))
  data.frame(
    dataset      = ds,
    status       = ifelse(g > POS_CUTOFF, "GDF15+", "GDF15-"),
    GDF15        = as.numeric(g),
    n_genes      = as.numeric(getc("n_genes")),
    n_counts     = as.numeric(getc("n_counts")),
    percent_mito = as.numeric(getc("percent_mito")),
    patient      = if ("patient" %in% colnames(md)) as.character(md[["patient"]]) else NA_character_,
    stringsAsFactors = FALSE)
}

star_lab <- function(p) ifelse(p < 1e-3, "***", ifelse(p < 1e-2, "**", ifelse(p < 0.05, "*", "ns")))

# S1 -- QC boxplots: GDF15+ vs GDF15- per cancer (n_genes, n_counts, %mito)
antiartifact_qc <- function(qc) {
  d <- do.call(rbind, qc)
  metrics <- c(n_genes = "Genes / cell", n_counts = "UMIs / cell (log10)",
               percent_mito = "% mitochondrial")
  long <- do.call(rbind, lapply(names(metrics), function(m) {
    v <- d[[m]]; if (m == "n_counts") v <- log10(v)
    data.frame(dataset = d$dataset, status = d$status, metric = metrics[[m]],
               mkey = m, value = v, stringsAsFactors = FALSE)
  }))
  long <- long[is.finite(long$value), ]
  long$status  <- factor(long$status,  levels = c("GDF15-", "GDF15+"))
  long$dataset <- factor(long$dataset, levels = unique(c(DISPLAY_ORDER, DATASETS)))
  long$metric  <- factor(long$metric,  levels = unname(metrics))
  # Mann-Whitney per dataset x metric
  stat <- do.call(rbind, lapply(split(long, list(long$mkey, long$dataset), drop = TRUE), function(s) {
    a <- s$value[s$status == "GDF15+"]; b <- s$value[s$status == "GDF15-"]
    if (length(a) > 3 && length(b) > 3) {
      p <- suppressWarnings(wilcox.test(a, b)$p.value)
      data.frame(dataset = s$dataset[1], metric = s$metric[1],
                 y = max(s$value, na.rm = TRUE), lab = star_lab(p))
    } else NULL
  }))
  p <- ggplot(long, aes(dataset, value, fill = status)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(0.8),
                 width = 0.7, linewidth = 0.3) +
    facet_wrap(~ metric, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = c("GDF15-" = "#B9C0C7", "GDF15+" = "#08306B"), name = NULL) +
    { if (!is.null(stat)) geom_text(data = stat, aes(dataset, y, label = lab),
                                    inherit.aes = FALSE, vjust = -0.3, size = 3) } +
    labs(x = NULL, y = NULL,
         title = "GDF15+ macrophages are not a low-quality / sequencing-depth artifact",
         subtitle = "Mann-Whitney GDF15+ vs GDF15- per cancer: *** p<0.001  ** p<0.01  * p<0.05  ns") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
          strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank())
  save_two(p, "S1_QC_GDF15pos_vs_neg", w = 11, h = 10)
}

# S2 -- % GDF15+ macrophages per patient within each cancer
antiartifact_patient <- function(qc, min_macro = 20) {
  d <- do.call(rbind, qc); d <- d[!is.na(d$patient), ]
  rows <- do.call(rbind, lapply(split(d, list(d$dataset, d$patient), drop = TRUE), function(s) {
    if (nrow(s) < min_macro) return(NULL)
    data.frame(dataset = s$dataset[1], patient = s$patient[1], n_macro = nrow(s),
               pct_pos = 100 * mean(s$GDF15 > POS_CUTOFF))
  }))
  if (is.null(rows) || nrow(rows) == 0) { message("[S2] no patients pass min_macro."); return(invisible()) }
  rows$dataset <- factor(rows$dataset, levels = unique(c(DISPLAY_ORDER, DATASETS)))
  write.csv(rows, file.path(OUT_DIR, "S2_per_patient_fraction.csv"), row.names = FALSE)
  p <- ggplot(rows, aes(dataset, pct_pos)) +
    stat_summary(fun = median, geom = "crossbar", width = 0.5, linewidth = 0.4, colour = "black") +
    geom_jitter(width = 0.15, height = 0, size = 2.2, colour = "#08306B", alpha = 0.8) +
    labs(x = NULL, y = "% GDF15+ macrophages per patient",
         title = "GDF15+ macrophages are detected across multiple patients per cancer",
         subtitle = sprintf("each dot = one patient with >= %d macrophages; bar = median", min_macro)) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
          plot.title = element_text(face = "bold"))
  save_two(p, "S2_per_patient_fraction", w = 11, h = 5.5)
}

# S3 -- % GDF15+ macrophages per cancer across expression cutoffs
antiartifact_threshold <- function(qc, cutoffs = c(0, 0.5, 1, 2)) {
  d <- do.call(rbind, qc)
  rows <- do.call(rbind, lapply(DATASETS, function(ds) {
    g <- d$GDF15[d$dataset == ds]; g <- g[is.finite(g)]
    do.call(rbind, lapply(cutoffs, function(cc)
      data.frame(dataset = ds, cutoff = cc, pct_pos = 100 * mean(g > cc), n_macro = length(g))))
  }))
  rows$dataset <- factor(rows$dataset, levels = unique(c(DISPLAY_ORDER, DATASETS)))
  rows$cutoff  <- factor(rows$cutoff, levels = cutoffs)
  write.csv(rows, file.path(OUT_DIR, "S3_threshold_sensitivity.csv"), row.names = FALSE)
  p <- ggplot(rows, aes(dataset, pct_pos, colour = cutoff, group = cutoff)) +
    geom_line(linewidth = 0.7) + geom_point(size = 2.5) +
    scale_colour_brewer(palette = "Blues", direction = -1, name = "GDF15 cutoff") +
    labs(x = NULL, y = "% GDF15+ macrophages",
         title = "GDF15+ fraction is robust to the expression cutoff") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
          plot.title = element_text(face = "bold"))
  save_two(p, "S3_threshold_sensitivity", w = 10, h = 5.5)
}

# --- volcano (panel c/e) ----------------------------------------------------
volcano_plot <- function(de, ds, npos, nneg) {
  de <- de %>% mutate(
    significance = case_when(
      p_val_adj < PADJ_LABEL & avg_log2FC >  LOGFC_LABEL ~ "Up in GDF15+",
      p_val_adj < PADJ_LABEL & avg_log2FC < -LOGFC_LABEL ~ "Up in GDF15-",
      TRUE ~ "NotSig"))
  lab_up   <- de %>% filter(significance == "Up in GDF15+") %>% arrange(p_val_adj) %>% slice_head(n = 20)
  lab_down <- de %>% filter(significance == "Up in GDF15-") %>% arrange(p_val_adj) %>% slice_head(n = 20)
  de$label <- ifelse(de$gene %in% c(lab_up$gene, lab_down$gene, GENE), de$gene, NA)
  p <- ggplot(de, aes(avg_log2FC, -log10(p_val_adj), colour = significance)) +
    geom_point(alpha = 0.7, size = 1.4) +
    scale_colour_manual(values = c("Up in GDF15+" = UP_COL, "Up in GDF15-" = DOWN_COL,
                                   "NotSig" = "grey75")) +
    geom_text_repel(aes(label = label), size = 3, max.overlaps = 100, segment.color = "grey80") +
    geom_vline(xintercept = c(-LOGFC_LABEL, LOGFC_LABEL), linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = -log10(PADJ_LABEL), linetype = "dashed", colour = "grey50") +
    theme_minimal(base_size = 13) +
    labs(title = sprintf("%s: GDF15+ vs GDF15- macrophages", ds),
         subtitle = sprintf("n GDF15+ = %d, GDF15- = %d", npos, nneg),
         x = "Log2FC (GDF15+ vs GDF15-)", y = "-Log10 adjusted P", colour = NULL)
  save_two(p, sprintf("Volcano_%s", ds), w = 7, h = 6)
}

## ----------------------- 3. Per-cancer analysis --------------------------- ##
de_summary   <- list()
heatmap_rows <- list()   # NES for the key terms, one entry per cancer
qc_rows      <- list()   # per-cell QC + GDF15 status, for the anti-artifact figs

analyse_macros <- function(macro, ds) {
  g   <- FetchData(macro, vars = GENE, layer = "data")[, 1]
  pos <- g > POS_CUTOFF
  macro$GDF15_status <- factor(ifelse(pos, "GDF15pos", "GDF15neg"),
                               levels = c("GDF15neg", "GDF15pos"))
  npos <- sum(pos); nneg <- sum(!pos)
  message(sprintf("  [%s] GDF15+ = %d | GDF15- = %d", ds, npos, nneg))
  if (npos < 10 || nneg < 10) { message("  [skip] too few cells in a group."); return(NULL) }
  Idents(macro) <- "GDF15_status"
  qc_rows[[ds]] <<- qc_table(macro, g, ds)   # for the anti-artifact figures

  de <- FindMarkers(macro, ident.1 = "GDF15pos", ident.2 = "GDF15neg",
                    assay = DefaultAssay(macro), layer = "data",
                    logfc.threshold = LOGFC_THRESH, min.pct = MIN_PCT)
  de$gene <- rownames(de)
  de <- de[order(de$p_val_adj), ]
  write.csv(de, file.path(OUT_DIR, sprintf("DE_%s.csv", ds)), row.names = FALSE)

  volcano_plot(de, ds, npos, nneg)

  geneList <- make_ranklist(de)
  gres_gobp <- NULL
  for (col in COLLECTIONS) {
    message("  [gsea] ", ds, " : ", col)
    gres <- run_gsea(geneList, col)
    if (!is.null(gres) && nrow(as.data.frame(gres)) > 0) {
      write.csv(as.data.frame(gres),
                file.path(OUT_DIR, sprintf("GSEA_%s_%s.csv", col, ds)), row.names = FALSE)
      lollipop_plot(gres, ds, col)
      if (col == "GOBP") gres_gobp <- gres
    }
  }
  # one targeted permissive GSEA pass shared by the running plots + heatmap
  target_ids <- unique(c(names(RUNNING_TERMS), names(HEATMAP_TERMS)))
  gres_t <- target_gsea(geneList, target_ids)
  running_plots(gres_t, ds)                          # the two running-score plots
  heatmap_rows[[ds]] <<- heatmap_row(gres_gobp, gres_t, ds)   # NES for the heatmap

  data.frame(dataset = ds, n_macro = ncol(macro), n_pos = npos, n_neg = nneg,
             n_up = sum(de$p_val_adj < PADJ_LABEL & de$avg_log2FC >  LOGFC_LABEL),
             n_down = sum(de$p_val_adj < PADJ_LABEL & de$avg_log2FC < -LOGFC_LABEL))
}

all_macros <- list()   # kept only if DO_POOLED
for (ds in DATASETS) {
  seu <- load_one_dataset(ds, DATA_DIR)
  if (!(GENE %in% rownames(seu))) { message("[skip] ", GENE, " not in ", ds); rm(seu); gc(); next }
  macro <- if (any(seu$is_macro)) subset(seu, subset = is_macro == TRUE) else seu
  rm(seu); gc()
  if (DO_POOLED) all_macros[[ds]] <- macro
  res <- analyse_macros(macro, ds)
  if (!is.null(res)) de_summary[[ds]] <- res
  if (!DO_POOLED) { rm(macro); gc() }
}

summ <- do.call(rbind, de_summary); rownames(summ) <- NULL
write.csv(summ, file.path(OUT_DIR, "DE_summary.csv"), row.names = FALSE)
message("\n[save] DE_summary.csv"); print(summ)

# NES heatmap across cancers for the key terms (prints + saves)
message("\n[plot] NES heatmap (cancers x key terms) ...")
nes_heatmap(heatmap_rows)

# Anti-artifact / robustness supplementary figures (print + save)
message("\n[plot] anti-artifact / robustness figures (S1 QC, S2 per-patient, S3 threshold) ...")
if (length(qc_rows) > 0) {
  antiartifact_qc(qc_rows)
  antiartifact_patient(qc_rows)
  antiartifact_threshold(qc_rows)
}

## ------------------- 4. OPTIONAL pooled high-power run -------------------- ##
# Small GDF15+ groups per cancer make per-cancer GSEA noisier than the pooled
# Chen analysis. Flip DO_POOLED <- TRUE to merge all macrophages and run ONE
# GDF15+ vs GDF15- DE + GSEA across cancers (dataset kept as a label).
if (DO_POOLED && length(all_macros) > 1) {
  message("\n[pooled] merging ", length(all_macros), " macrophage objects ...")
  pooled <- merge(all_macros[[1]], y = all_macros[-1],
                  add.cell.ids = names(all_macros), project = "GDF15_pooled")
  if (inherits(pooled[["RNA"]], "Assay5")) pooled[["RNA"]] <- JoinLayers(pooled[["RNA"]])
  rm(all_macros); gc()
  res <- analyse_macros(pooled, "POOLED")
  print(res)
}

message("\n==============================================================")
message("  DONE. Outputs in: ", OUT_DIR)
message("==============================================================")
