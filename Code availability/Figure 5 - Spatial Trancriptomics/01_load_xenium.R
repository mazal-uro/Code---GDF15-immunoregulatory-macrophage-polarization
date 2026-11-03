# 01_load_xenium.R
# Load the Takehara et al. prostate Xenium data deposited at GEO GSE296275.
#
# IMPORTANT: GEO deposits Xenium output in a *flattened* form — not the standard
# Xenium output bundle that Seurat::LoadXenium() expects. We therefore build
# the Seurat object manually:
#
#   * Read the cell x gene count matrix from the .mtx + barcodes + features triple.
#   * Read transcripts.parquet.gz to derive per-cell centroids (median x, y of
#     all transcripts assigned to each cell).
#   * Wrap the resulting (matrix + centroids) into a Seurat v5 object with one
#     FOV per sample.
#
# Inputs (in /Users/alekhyamazumdar/Documents/Xenium/GSE296275_RAW):
#   GSM8967906_D7_Prostate_1_cell_feature_matrix_barcodes.tsv.gz
#   GSM8967906_D7_Prostate_1_cell_feature_matrix_features.tsv.gz
#   GSM8967906_D7_Prostate_1_cell_feature_matrix_matrix.mtx.gz
#   GSM8967906_D7_Prostate_1_transcripts.parquet.gz
#   ... and the same five-file pattern for GSM8967907..GSM8967909.
#
# (The .ome.tif.gz morphology images are not needed for the analysis pipeline
#  and are skipped.)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(arrow)
  library(data.table)
  library(R.utils)        # gunzip
  library(future)
  library(dplyr)
})

plan("multisession", workers = 4)
options(future.globals.maxSize = 8 * 1024^3)

# ---------- Paths ----------
DATA_DIR <- "/Users/alekhyamazumdar/Documents/Xenium/GSE296275_RAW"
OUT_DIR  <- "/Users/alekhyamazumdar/Documents/Xenium/Output"
TMP_DIR  <- file.path(tempdir(), "xenium_parquet")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TMP_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------- Sample table ----------
# GSM ID -> human-readable sample ID + metadata. Order matters only for output;
# all samples are merged at the end.
sample_meta <- tibble::tribble(
  ~gsm,        ~sample_id,        ~timepoint, ~replicate,
  "GSM8967906", "D7_Prostate_1",  "7dpi",     "rep1",
  "GSM8967907", "D7_Prostate_2",  "7dpi",     "rep2",
  "GSM8967908", "D30_Prostate_1", "30dpi",    "rep1",
  "GSM8967909", "D30_Prostate_2", "30dpi",    "rep2"
)

# Helper: construct expected filename
.fname <- function(gsm, sample, suffix) {
  file.path(DATA_DIR, paste0(gsm, "_", sample, "_", suffix))
}

# ---------- Helper: read the cell x gene matrix ----------
read_xenium_matrix <- function(gsm, sample) {
  mtx_path <- .fname(gsm, sample, "cell_feature_matrix_matrix.mtx.gz")
  bc_path  <- .fname(gsm, sample, "cell_feature_matrix_barcodes.tsv.gz")
  ft_path  <- .fname(gsm, sample, "cell_feature_matrix_features.tsv.gz")

  message("  reading matrix: ", basename(mtx_path))
  m  <- readMM(gzfile(mtx_path))
  bc <- readLines(gzfile(bc_path))
  ft <- read.table(gzfile(ft_path), sep = "\t", header = FALSE,
                   stringsAsFactors = FALSE)
  # features.tsv has 3 cols: ensembl_id, gene_symbol, feature_type
  rownames(m) <- make.unique(ft$V2)         # gene symbol
  colnames(m) <- bc                          # cell barcode
  m <- as(m, "CsparseMatrix")

  # Drop Negative Control / Codeword / Unassigned probes that aren't real genes
  keep <- ft$V3 == "Gene Expression"
  m <- m[keep, , drop = FALSE]

  list(matrix = m, features = ft[keep, ])
}

# ---------- Helper: derive centroids from transcripts.parquet.gz ----------
# Decompresses to a temp file, reads with arrow, computes median x/y per cell.
derive_centroids <- function(gsm, sample, valid_cells) {
  pq_gz <- .fname(gsm, sample, "transcripts.parquet.gz")
  pq    <- file.path(TMP_DIR, paste0(sample, "_transcripts.parquet"))
  if (!file.exists(pq)) {
    message("  decompressing parquet: ", basename(pq_gz))
    gunzip(pq_gz, destname = pq, remove = FALSE, overwrite = TRUE)
  }
  message("  reading parquet schema")
  ds <- open_dataset(pq, format = "parquet")
  print(ds$schema)

  # Standard Xenium column names: cell_id, x_location, y_location, feature_name.
  # If the schema differs, edit the next 3 lines.
  cell_col <- "cell_id"
  x_col    <- "x_location"
  y_col    <- "y_location"

  message("  computing centroids (median x,y per cell_id)")
  # Arrow's lazy aggregation does not support median(). We project + filter in
  # arrow (fast, low-memory), then collect and compute medians in R with
  # data.table for speed. Also drop non-gene transcripts (control codewords).
  cols <- c(cell_col, x_col, y_col)
  if ("is_gene" %in% names(ds$schema)) {
    txs <- ds |>
      dplyr::filter(is_gene == TRUE,
                    .data[[cell_col]] %in% valid_cells) |>
      dplyr::select(dplyr::all_of(cols)) |>
      dplyr::collect()
  } else {
    txs <- ds |>
      dplyr::filter(.data[[cell_col]] %in% valid_cells) |>
      dplyr::select(dplyr::all_of(cols)) |>
      dplyr::collect()
  }
  setDT(txs)
  setnames(txs, c(cell_col, x_col, y_col), c("cell", "x", "y"))
  centroids <- txs[, .(x = median(x, na.rm = TRUE),
                       y = median(y, na.rm = TRUE)),
                   by = cell]
  as.data.frame(centroids)
}

# ---------- Build one Seurat object per sample ----------
build_one <- function(gsm, sample_id, timepoint, replicate) {
  message("=== ", sample_id, " (", gsm, ") ===")
  mtx_info <- read_xenium_matrix(gsm, sample_id)
  m <- mtx_info$matrix

  cents <- derive_centroids(gsm, sample_id, valid_cells = colnames(m))
  # Keep only cells that have both transcripts and matrix entries
  keep_cells <- intersect(colnames(m), cents$cell)
  m     <- m[, keep_cells, drop = FALSE]
  cents <- cents[match(keep_cells, cents$cell), ]
  stopifnot(all(cents$cell == colnames(m)))

  # Build Seurat v5 assay + object
  obj <- CreateSeuratObject(counts = m, assay = "Xenium", project = sample_id)
  obj$sample_id <- sample_id
  obj$timepoint <- timepoint
  obj$replicate <- replicate

  # Build FOV with cell centroids and attach
  cents_obj <- CreateCentroids(
    coords = data.frame(x = cents$x, y = cents$y, cell = cents$cell)
  )
  fov <- CreateFOV(
    coords     = list(centroids = cents_obj),
    type       = "centroids",
    assay      = "Xenium",
    key        = paste0(make.names(sample_id), "_")
  )
  obj[[sample_id]] <- fov
  obj
}

obj_list <- mapply(
  build_one,
  sample_meta$gsm,
  sample_meta$sample_id,
  sample_meta$timepoint,
  sample_meta$replicate,
  SIMPLIFY = FALSE
)
names(obj_list) <- sample_meta$sample_id

# ---------- Merge ----------
xen <- merge(obj_list[[1]], y = obj_list[-1], add.cell.ids = names(obj_list))
DefaultAssay(xen) <- "Xenium"

# ---------- Sanity check ----------
print(xen)
cat("Total cells:    ", ncol(xen), "\n")
cat("Genes in panel: ", nrow(xen), "\n")
print(table(xen$sample_id, xen$timepoint))

saveRDS(xen, file.path(OUT_DIR, "01_xenium_raw.rds"))
message("Saved: ", file.path(OUT_DIR, "01_xenium_raw.rds"))

# Optional cleanup of decompressed parquets to save disk:
# unlink(TMP_DIR, recursive = TRUE)
