#!/usr/bin/env Rscript
# =============================================================================
# 09_zoom_neighborhood.R
# Zoom-in representative panels for one random GDF15-high and one random
# GDF15-low/neg macrophage.
#
# Each run picks NEW random macrophages (seeded from the system clock).
# Output is timestamped so previous runs are not overwritten.
#
# Layout  — 2 rows × 3 columns:
#   Row 1  GDF15-high macrophage  │ Cell types │ GDF15 expr │ TGFB1 expr │
#   Row 2  GDF15-low macrophage   │ Cell types │ GDF15 expr │ TGFB1 expr │
#
# Colour scales for GDF15 and TGFB1 are SHARED across both rows so the two
# neighbourhoods are directly comparable.
# The red dashed circle marks the 100 µm neighbourhood boundary.
# The focal macrophage is shown as a filled diamond (◆).
# =============================================================================

library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)

# ── User-adjustable parameters ────────────────────────────────────────────────
output_dir <- "/Users/alekhyamazumdar/Documents/Xenium/GDF15_analysis"
RADIUS_UM  <- 100   # µm — neighbourhood radius (must match main analysis)
PAD_UM     <- 150   # µm — half-width of zoom window (>= RADIUS_UM)
POINT_SIZE <- 1.8   # background cell dot size
FOCAL_SIZE <- 6     # focal macrophage diamond size

# ── Load saved objects ────────────────────────────────────────────────────────
message("Loading Seurat objects ...")
macro  <- readRDS(file.path(output_dir, "08_macrophage_GDF15_xenium.rds"))
xenium <- readRDS(file.path(output_dir, "08_xenium_full_annotated.rds"))

# ── Fresh random seed from clock — different macrophages every run ─────────────
set.seed(as.integer(Sys.time()) %% 100000L)

# ── Verify GDF15_status was saved in macro ────────────────────────────────────
if (!"GDF15_status" %in% colnames(macro@meta.data)) {
  stop("GDF15_status not found in macro metadata. Re-run the main script first.")
}

mac_meta     <- macro@meta.data
high_bcs     <- rownames(mac_meta)[mac_meta$GDF15_status == "GDF15-high"]
low_bcs      <- rownames(mac_meta)[mac_meta$GDF15_status == "GDF15-low/neg"]
cat(sprintf("Pool  —  GDF15-high: %d  |  GDF15-low/neg: %d\n",
            length(high_bcs), length(low_bcs)))

focal_high <- sample(high_bcs, 1)
focal_low  <- sample(low_bcs,  1)

all_meta <- xenium@meta.data   # has x_centroid, y_centroid, auto_celltype

cat(sprintf("Focal GDF15-high : %s  (x=%.1f, y=%.1f)\n",
            focal_high,
            all_meta[focal_high, "x_centroid"],
            all_meta[focal_high, "y_centroid"]))
cat(sprintf("Focal GDF15-low  : %s  (x=%.1f, y=%.1f)\n",
            focal_low,
            all_meta[focal_low,  "x_centroid"],
            all_meta[focal_low,  "y_centroid"]))

# ── Extract expression for all xenium cells ───────────────────────────────────
message("Extracting expression data ...")
expr_data <- GetAssayData(xenium, assay = "Xenium", layer = "data")

tgfb_candidates <- intersect(c("TGFB1", "TGFB2", "TGFB3"), rownames(expr_data))
if (length(tgfb_candidates) == 0) stop("No TGFB1/2/3 found in panel.")
TGFB_GENE <- tgfb_candidates[1]          # TGFB1 preferred
cat(sprintf("TGFB gene : %s\n", TGFB_GENE))

all_meta$gdf15 <- as.numeric(expr_data["GDF15",   rownames(all_meta)])
all_meta$tgfb  <- as.numeric(expr_data[TGFB_GENE, rownames(all_meta)])

# ── Collect cells in each window to compute SHARED colour limits ───────────────
get_window <- function(focal_bc, all_meta, pad_um) {
  cx <- all_meta[focal_bc, "x_centroid"]
  cy <- all_meta[focal_bc, "y_centroid"]
  all_meta[abs(all_meta$x_centroid - cx) <= pad_um &
           abs(all_meta$y_centroid - cy) <= pad_um, ]
}

win_high <- get_window(focal_high, all_meta, PAD_UM)
win_low  <- get_window(focal_low,  all_meta, PAD_UM)

# Shared limits across both windows
gdf15_lim <- c(0, max(c(win_high$gdf15, win_low$gdf15), na.rm = TRUE))
tgfb_lim  <- c(0, max(c(win_high$tgfb,  win_low$tgfb),  na.rm = TRUE))
if (gdf15_lim[2] < 1e-6) gdf15_lim[2] <- 1   # guard against all-zero
if (tgfb_lim[2]  < 1e-6) tgfb_lim[2]  <- 1

cat(sprintf("Shared GDF15 scale : [%.2f, %.2f]\n", gdf15_lim[1], gdf15_lim[2]))
cat(sprintf("Shared %s scale : [%.2f, %.2f]\n", TGFB_GENE, tgfb_lim[1], tgfb_lim[2]))

# ── Colour palette for cell types (shared set of types) ───────────────────────
all_ct   <- sort(unique(c(as.character(win_high$auto_celltype),
                           as.character(win_low$auto_celltype))))
ct_cols  <- setNames(scales::hue_pal()(length(all_ct)), all_ct)

# ── Shared zoom theme ─────────────────────────────────────────────────────────
theme_zoom <- theme_cowplot(font_size = 10) +
  theme(
    plot.title      = element_text(face = "bold", size = 10, hjust = 0.5),
    axis.title      = element_text(size = 8.5),
    axis.text       = element_text(size = 7.5),
    legend.title    = element_text(size = 8, face = "bold"),
    legend.text     = element_text(size = 7.5),
    legend.key.size = unit(0.45, "cm"),
    plot.background = element_rect(fill = "white", color = NA)
  )

# ── Helper: build three panels for one focal macrophage ───────────────────────
make_row <- function(focal_bc, group_label, focal_fill,
                     win, ct_cols, gdf15_lim, tgfb_lim, TGFB_GENE,
                     radius_um, pad_um, point_size, focal_size) {

  cx <- all_meta[focal_bc, "x_centroid"]
  cy <- all_meta[focal_bc, "y_centroid"]

  # 100 µm dashed circle
  theta  <- seq(0, 2 * pi, length.out = 500)
  circ   <- data.frame(x = cx + radius_um * cos(theta),
                        y = cy + radius_um * sin(theta))

  # Focal cell marker
  focal_df <- data.frame(x = cx, y = cy)

  xlims <- c(cx - pad_um, cx + pad_um)
  ylims <- c(cy - pad_um, cy + pad_um)

  # Factor cell types so legend order is consistent
  win$ct <- factor(as.character(win$auto_celltype), levels = names(ct_cols))

  # ── P1: Cell type ────────────────────────────────────────────────────────
  p1 <- ggplot(win, aes(x = x_centroid, y = y_centroid, color = ct)) +
    geom_point(size = point_size, alpha = 0.88, stroke = 0) +
    geom_path(data = circ, aes(x = x, y = y),
              inherit.aes = FALSE,
              color = "red", linetype = "dashed", linewidth = 1.0) +
    geom_point(data = focal_df, aes(x = x, y = y),
               inherit.aes = FALSE,
               shape = 23, size = focal_size,
               fill = focal_fill, color = "black", stroke = 1.4) +
    scale_color_manual(values = ct_cols, name = "Cell type",
                       drop = FALSE) +
    coord_fixed(xlim = xlims, ylim = ylims) +
    labs(title = "Cell types",
         x = "x (µm)", y = "y (µm)") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1),
                                ncol = 1)) +
    theme_zoom

  # ── P2: GDF15 expression (render high expressors on top) ─────────────────
  win_g <- win[order(win$gdf15), ]
  p2 <- ggplot(win_g, aes(x = x_centroid, y = y_centroid, color = gdf15)) +
    geom_point(size = point_size, alpha = 0.90, stroke = 0) +
    geom_path(data = circ, aes(x = x, y = y),
              inherit.aes = FALSE,
              color = "red", linetype = "dashed", linewidth = 1.0) +
    geom_point(data = focal_df, aes(x = x, y = y),
               inherit.aes = FALSE,
               shape = 23, size = focal_size,
               fill = focal_fill, color = "black", stroke = 1.4) +
    scale_color_gradientn(
      colors = c("grey92", "#CDEAF5", "#6BB8D4", "#1D7DB5", "#0B3060"),
      limits = gdf15_lim, name = "GDF15",
      guide  = guide_colorbar(barwidth = unit(3.5, "cm"), barheight = unit(0.35, "cm"),
                               title.position = "top", title.hjust = 0.5)) +
    coord_fixed(xlim = xlims, ylim = ylims) +
    labs(title = "GDF15 expression",
         x = "x (µm)", y = "y (µm)") +
    theme_zoom +
    theme(legend.position = "bottom")

  # ── P3: TGFB expression (render high expressors on top) ───────────────────
  win_t <- win[order(win$tgfb), ]
  p3 <- ggplot(win_t, aes(x = x_centroid, y = y_centroid, color = tgfb)) +
    geom_point(size = point_size, alpha = 0.90, stroke = 0) +
    geom_path(data = circ, aes(x = x, y = y),
              inherit.aes = FALSE,
              color = "red", linetype = "dashed", linewidth = 1.0) +
    geom_point(data = focal_df, aes(x = x, y = y),
               inherit.aes = FALSE,
               shape = 23, size = focal_size,
               fill = focal_fill, color = "black", stroke = 1.4) +
    scale_color_gradientn(
      colors = c("grey92", "#FFF3B0", "#FDAE61", "#D73027", "#7B0011"),
      limits = tgfb_lim, name = TGFB_GENE,
      guide  = guide_colorbar(barwidth = unit(3.5, "cm"), barheight = unit(0.35, "cm"),
                               title.position = "top", title.hjust = 0.5)) +
    coord_fixed(xlim = xlims, ylim = ylims) +
    labs(title = sprintf("%s expression", TGFB_GENE),
         x = "x (µm)", y = "y (µm)") +
    theme_zoom +
    theme(legend.position = "bottom")

  # Row tag
  tag <- sprintf("%s macrophage", group_label)
  list(p1 = p1, p2 = p2, p3 = p3, tag = tag)
}

# ── Build rows ────────────────────────────────────────────────────────────────
message("Building GDF15-high panels ...")
row_high <- make_row(focal_high, "GDF15-high",   "#D62728",
                     win_high, ct_cols, gdf15_lim, tgfb_lim, TGFB_GENE,
                     RADIUS_UM, PAD_UM, POINT_SIZE, FOCAL_SIZE)

message("Building GDF15-low panels ...")
row_low  <- make_row(focal_low,  "GDF15-low/neg", "#4E9AC7",
                     win_low, ct_cols, gdf15_lim, tgfb_lim, TGFB_GENE,
                     RADIUS_UM, PAD_UM, POINT_SIZE, FOCAL_SIZE)

# ── Assemble 2 × 3 figure ─────────────────────────────────────────────────────
fig <- (row_high$p1 | row_high$p2 | row_high$p3) /
       (row_low$p1  | row_low$p2  | row_low$p3)  +
  plot_annotation(
    title    = sprintf(
      "Macrophage neighbourhood zoom  (%d µm radius)  —  ◆ = focal macrophage",
      RADIUS_UM),
    subtitle = sprintf(
      "Top row (red ◆): GDF15-high  [%s]\nBottom row (blue ◆): GDF15-low/neg  [%s]",
      focal_high, focal_low),
    caption  = sprintf(
      "Colour scales shared across both rows  •  %s window  •  %s",
      paste0(PAD_UM * 2, " µm × ", PAD_UM * 2, " µm"),
      format(Sys.time(), "%Y-%m-%d %H:%M")),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13, hjust = 0),
      plot.subtitle = element_text(size = 9.5, color = "grey40", hjust = 0,
                                   margin = margin(b = 6)),
      plot.caption  = element_text(size = 7.5, color = "grey55", hjust = 1)
    )
  )

# ── Save with timestamp ───────────────────────────────────────────────────────
ts       <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_path <- file.path(output_dir, sprintf("09_neighborhood_zoom_%s.png", ts))
ggsave(out_path, fig, width = 18, height = 12, dpi = 300, bg = "white")
message(sprintf("\nSaved → %s", out_path))
print(fig)
