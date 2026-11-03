# =============================================================================
# neighborhood_analysis_celltype.R
# Macrophage neighbourhood analysis using the NEW cell-type annotation
# (xenium$cell_type from annotate_celltypes_findmarkers.R), instead of the
# old module-score auto_celltype.
#
# Produces:
#   A. Neighbourhood enrichment lollipop — cell-type composition within RADIUS
#      of GDF15-high vs GDF15-low/neg macrophages (Δ proportion).
#   B. Neighbourhood zoom panels — one random GDF15-high and one GDF15-low
#      macrophage, showing cell types | GDF15 | TGFB1 (shared colour scales).
#
# PRE-REQUISITE: run annotate_celltypes_findmarkers.R first — it saves
#   Output/annot_xenium_celltype_annotated.rds  (carries xenium$cell_type).
#
# Macrophages and their GDF15-high/low split are derived FRESH from the new
# annotation here, so nothing depends on the old 08_* objects.
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)
library(Matrix)
library(scales)

# ── Paths & parameters ────────────────────────────────────────────────────────
output_dir <- "/Users/alekhyamazumdar/Documents/Xenium/Output"
annot_rds  <- file.path(output_dir, "annot_xenium_celltype_annotated.rds")

RADIUS_UM  <- 100    # µm — neighbourhood radius
PAD_UM     <- 150    # µm — half-width of zoom window (>= RADIUS_UM)
POINT_SIZE <- 1.8    # background cell dot size (zoom)
FOCAL_SIZE <- 6      # focal macrophage diamond size (zoom)
MAC_LABEL  <- "Macrophage"   # the cell_type value that denotes macrophages

theme_pub <- theme_classic(base_size = 12) +
  theme(plot.title   = element_text(face = "bold", size = 13),
        axis.text    = element_text(color = "black"),
        legend.title = element_text(face = "bold"))

# ── Load the newly annotated object ───────────────────────────────────────────
if (!file.exists(annot_rds))
  stop("Annotated object not found:\n  ", annot_rds,
       "\nRun annotate_celltypes_findmarkers.R first.")
message("Loading annotated Xenium object ...")
xenium <- readRDS(annot_rds)
stopifnot("cell_type"  %in% colnames(xenium@meta.data))
stopifnot("x_centroid" %in% colnames(xenium@meta.data))
if (!MAC_LABEL %in% levels(factor(xenium$cell_type)))
  stop(sprintf("cell_type has no '%s' label. Present: %s",
               MAC_LABEL, paste(levels(factor(xenium$cell_type)), collapse = ", ")))

# ── Expression + coordinates ──────────────────────────────────────────────────
expr_data <- GetAssayData(xenium, assay = "Xenium", layer = "data")
gdf15_all <- setNames(as.numeric(expr_data["GDF15", ]), colnames(xenium))

all_coords    <- data.frame(x = xenium$x_centroid, y = xenium$y_centroid,
                            row.names = colnames(xenium))
all_celltypes <- setNames(as.character(xenium$cell_type), colnames(xenium))

# ── Define macrophages and their GDF15-high/low split (new annotation) ────────
mac_bcs   <- colnames(xenium)[as.character(xenium$cell_type) == MAC_LABEL]
mac_gdf15 <- gdf15_all[mac_bcs]
if (sum(mac_gdf15 > 0) < 4)
  stop("Too few GDF15-expressing macrophages to define a high/low split.")
gdf15_threshold <- quantile(mac_gdf15[mac_gdf15 > 0], 0.5)   # median of expressing
mac_status <- ifelse(mac_gdf15 > gdf15_threshold, "GDF15-high", "GDF15-low/neg")
names(mac_status) <- mac_bcs

high_bcs <- names(mac_status)[mac_status == "GDF15-high"]
low_bcs  <- names(mac_status)[mac_status == "GDF15-low/neg"]
cat(sprintf("Macrophages (%s): %d\n", MAC_LABEL, length(mac_bcs)))
cat(sprintf("GDF15 threshold (median of expressing): %.3f\n", gdf15_threshold))
cat(sprintf("GDF15-high: %d  |  GDF15-low/neg: %d\n",
            length(high_bcs), length(low_bcs)))


# =============================================================================
# SECTION A — NEIGHBOURHOOD ENRICHMENT (Δ proportion by cell type)
# =============================================================================
message("[A] Neighbourhood enrichment: GDF15-high vs GDF15-low/neg ...")

# Balance sample sizes
n_compare <- min(length(high_bcs), length(low_bcs))
set.seed(42)
high_sample <- sample(high_bcs, n_compare)
low_sample  <- sample(low_bcs,  n_compare)
cat(sprintf("Comparison: %d GDF15-high vs %d GDF15-low/neg macrophages\n",
            n_compare, n_compare))

# Cell types of every non-focal cell within radius of any focal cell
get_neighborhood_types <- function(focal_barcodes, all_coords, all_celltypes, radius_um) {
  fx <- all_coords[focal_barcodes, "x"]; fy <- all_coords[focal_barcodes, "y"]
  ax <- all_coords$x; ay <- all_coords$y; r2 <- radius_um^2
  unlist(lapply(seq_along(focal_barcodes), function(i) {
    d2 <- (ax - fx[i])^2 + (ay - fy[i])^2
    all_celltypes[d2 > 0 & d2 <= r2]
  }))
}

nbr_high <- get_neighborhood_types(high_sample, all_coords, all_celltypes, RADIUS_UM)
nbr_low  <- get_neighborhood_types(low_sample,  all_coords, all_celltypes, RADIUS_UM)

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
nbr_df <- nbr_df[order(-nbr_df$diff), ]

cat("\n── Neighbourhood composition (Δ = GDF15-high − GDF15-low) ──\n")
print(nbr_df, row.names = FALSE, digits = 3)
write.csv(nbr_df, file.path(output_dir, "annot_neighborhood_composition_celltype.csv"),
          row.names = FALSE)

nbr_df$celltype_f <- factor(nbr_df$celltype,
                            levels = nbr_df$celltype[order(nbr_df$diff)])
p_nbr <- ggplot(nbr_df, aes(x = diff, y = celltype_f, color = diff > 0)) +
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = "dashed") +
  geom_segment(aes(x = 0, xend = diff, yend = celltype_f), linewidth = 1.2) +
  geom_point(size = 4) +
  scale_color_manual(
    values = c("TRUE" = "#D62728", "FALSE" = "#4E9AC7"),
    labels = c("TRUE" = "Enriched near GDF15-high",
               "FALSE" = "Enriched near GDF15-low/neg")) +
  labs(
    title = sprintf("Neighbourhood enrichment (new annotation)\nwithin %d µm, n = %d per group",
                    RADIUS_UM, n_compare),
    x = "Δ proportion  (GDF15-high − GDF15-low/neg)", y = NULL, color = NULL) +
  theme_pub
ggsave(file.path(output_dir, "annot_neighborhood_enrichment_lollipop.png"),
       p_nbr, width = 9, height = 6, dpi = 300, bg = "white")
print(p_nbr)


# =============================================================================
# SECTION B — NEIGHBOURHOOD ZOOM PANELS (cell types | GDF15 | TGFB1)
# =============================================================================
message("[B] Neighbourhood zoom panels ...")

# TGFB gene present in panel
tgfb_candidates <- intersect(c("TGFB1", "TGFB2", "TGFB3"), rownames(expr_data))
if (length(tgfb_candidates) == 0) stop("No TGFB1/2/3 found in panel.")
TGFB_GENE <- tgfb_candidates[1]
cat(sprintf("TGFB gene: %s\n", TGFB_GENE))

# Fresh random focal cells every run
set.seed(as.integer(Sys.time()) %% 100000L)
focal_high <- sample(high_bcs, 1)
focal_low  <- sample(low_bcs,  1)

all_meta <- xenium@meta.data
all_meta$gdf15 <- as.numeric(expr_data["GDF15",   rownames(all_meta)])
all_meta$tgfb  <- as.numeric(expr_data[TGFB_GENE, rownames(all_meta)])

cat(sprintf("Focal GDF15-high : %s  (x=%.1f, y=%.1f)\n", focal_high,
            all_meta[focal_high, "x_centroid"], all_meta[focal_high, "y_centroid"]))
cat(sprintf("Focal GDF15-low  : %s  (x=%.1f, y=%.1f)\n", focal_low,
            all_meta[focal_low,  "x_centroid"], all_meta[focal_low,  "y_centroid"]))

get_window <- function(focal_bc, all_meta, pad_um) {
  cx <- all_meta[focal_bc, "x_centroid"]; cy <- all_meta[focal_bc, "y_centroid"]
  all_meta[abs(all_meta$x_centroid - cx) <= pad_um &
           abs(all_meta$y_centroid - cy) <= pad_um, ]
}
win_high <- get_window(focal_high, all_meta, PAD_UM)
win_low  <- get_window(focal_low,  all_meta, PAD_UM)

# Shared colour limits across both windows
gdf15_lim <- c(0, max(c(win_high$gdf15, win_low$gdf15), na.rm = TRUE))
tgfb_lim  <- c(0, max(c(win_high$tgfb,  win_low$tgfb),  na.rm = TRUE))
if (gdf15_lim[2] < 1e-6) gdf15_lim[2] <- 1
if (tgfb_lim[2]  < 1e-6) tgfb_lim[2]  <- 1

# Cell-type palette across the union of types in both windows (NEW cell_type)
all_ct  <- sort(unique(c(as.character(win_high$cell_type),
                         as.character(win_low$cell_type))))
ct_cols <- setNames(scales::hue_pal()(length(all_ct)), all_ct)

theme_zoom <- theme_cowplot(font_size = 10) +
  theme(plot.title      = element_text(face = "bold", size = 10, hjust = 0.5),
        axis.title      = element_text(size = 8.5),
        axis.text       = element_text(size = 7.5),
        legend.title    = element_text(size = 8, face = "bold"),
        legend.text     = element_text(size = 7.5),
        legend.key.size = unit(0.45, "cm"),
        plot.background = element_rect(fill = "white", color = NA))

make_row <- function(focal_bc, group_label, focal_fill, win, ct_cols,
                     gdf15_lim, tgfb_lim, TGFB_GENE,
                     radius_um, pad_um, point_size, focal_size) {
  cx <- all_meta[focal_bc, "x_centroid"]; cy <- all_meta[focal_bc, "y_centroid"]
  theta    <- seq(0, 2 * pi, length.out = 500)
  circ     <- data.frame(x = cx + radius_um * cos(theta),
                         y = cy + radius_um * sin(theta))
  focal_df <- data.frame(x = cx, y = cy)
  xlims <- c(cx - pad_um, cx + pad_um)
  ylims <- c(cy - pad_um, cy + pad_um)

  win$ct <- factor(as.character(win$cell_type), levels = names(ct_cols))

  # P1 — cell types (new annotation)
  p1 <- ggplot(win, aes(x = x_centroid, y = y_centroid, color = ct)) +
    geom_point(size = point_size, alpha = 0.88, stroke = 0) +
    geom_path(data = circ, aes(x = x, y = y), inherit.aes = FALSE,
              color = "red", linetype = "dashed", linewidth = 1.0) +
    geom_point(data = focal_df, aes(x = x, y = y), inherit.aes = FALSE,
               shape = 23, size = focal_size, fill = focal_fill,
               color = "black", stroke = 1.4) +
    scale_color_manual(values = ct_cols, name = "Cell type", drop = FALSE) +
    coord_fixed(xlim = xlims, ylim = ylims) +
    labs(title = "Cell types", x = "x (µm)", y = "y (µm)") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1)) +
    theme_zoom

  # P2 — GDF15
  p2 <- ggplot(win[order(win$gdf15), ],
               aes(x = x_centroid, y = y_centroid, color = gdf15)) +
    geom_point(size = point_size, alpha = 0.90, stroke = 0) +
    geom_path(data = circ, aes(x = x, y = y), inherit.aes = FALSE,
              color = "red", linetype = "dashed", linewidth = 1.0) +
    geom_point(data = focal_df, aes(x = x, y = y), inherit.aes = FALSE,
               shape = 23, size = focal_size, fill = focal_fill,
               color = "black", stroke = 1.4) +
    scale_color_gradientn(
      colors = c("grey92", "#CDEAF5", "#6BB8D4", "#1D7DB5", "#0B3060"),
      limits = gdf15_lim, name = "GDF15",
      guide  = guide_colorbar(barwidth = unit(3.5, "cm"), barheight = unit(0.35, "cm"),
                              title.position = "top", title.hjust = 0.5)) +
    coord_fixed(xlim = xlims, ylim = ylims) +
    labs(title = "GDF15 expression", x = "x (µm)", y = "y (µm)") +
    theme_zoom + theme(legend.position = "bottom")

  # P3 — TGFB
  p3 <- ggplot(win[order(win$tgfb), ],
               aes(x = x_centroid, y = y_centroid, color = tgfb)) +
    geom_point(size = point_size, alpha = 0.90, stroke = 0) +
    geom_path(data = circ, aes(x = x, y = y), inherit.aes = FALSE,
              color = "red", linetype = "dashed", linewidth = 1.0) +
    geom_point(data = focal_df, aes(x = x, y = y), inherit.aes = FALSE,
               shape = 23, size = focal_size, fill = focal_fill,
               color = "black", stroke = 1.4) +
    scale_color_gradientn(
      colors = c("grey92", "#FFF3B0", "#FDAE61", "#D73027", "#7B0011"),
      limits = tgfb_lim, name = TGFB_GENE,
      guide  = guide_colorbar(barwidth = unit(3.5, "cm"), barheight = unit(0.35, "cm"),
                              title.position = "top", title.hjust = 0.5)) +
    coord_fixed(xlim = xlims, ylim = ylims) +
    labs(title = sprintf("%s expression", TGFB_GENE), x = "x (µm)", y = "y (µm)") +
    theme_zoom + theme(legend.position = "bottom")

  list(p1 = p1, p2 = p2, p3 = p3)
}

row_high <- make_row(focal_high, "GDF15-high", "#D62728", win_high, ct_cols,
                     gdf15_lim, tgfb_lim, TGFB_GENE,
                     RADIUS_UM, PAD_UM, POINT_SIZE, FOCAL_SIZE)
row_low  <- make_row(focal_low,  "GDF15-low/neg", "#4E9AC7", win_low, ct_cols,
                     gdf15_lim, tgfb_lim, TGFB_GENE,
                     RADIUS_UM, PAD_UM, POINT_SIZE, FOCAL_SIZE)

fig <- (row_high$p1 | row_high$p2 | row_high$p3) /
       (row_low$p1  | row_low$p2  | row_low$p3) +
  plot_annotation(
    title    = sprintf("Macrophage neighbourhood zoom  (%d µm radius)  —  focal macrophage = diamond",
                       RADIUS_UM),
    subtitle = sprintf("Top row (red): GDF15-high [%s]   •   Bottom row (blue): GDF15-low/neg [%s]\nNew cell-type annotation",
                       focal_high, focal_low),
    caption  = "Colour scales shared across both rows",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13, hjust = 0),
      plot.subtitle = element_text(size = 9.5, color = "grey40", hjust = 0,
                                   margin = margin(b = 6)),
      plot.caption  = element_text(size = 7.5, color = "grey55", hjust = 1)))

ts       <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_path <- file.path(output_dir, sprintf("annot_neighborhood_zoom_%s.png", ts))
ggsave(out_path, fig, width = 18, height = 12, dpi = 300, bg = "white")
message(sprintf("Zoom panels saved -> %s", out_path))
print(fig)

message("Neighbourhood analysis (new annotation) complete.")
