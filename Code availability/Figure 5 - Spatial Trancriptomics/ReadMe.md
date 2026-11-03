Spatial Transcriptomics of GDF15⁺ Macrophages — Figure 5a–d

R scripts for the 10x Xenium Prime analysis of GDF15-expressing macrophages in human prostate adenocarcinoma (Figure 5a–d). Panels 5e–f (THP-1 pSMAD2 western blot and GDF15 qPCR) are wet-lab data and are not part of this repository.

Overview

Cells are clustered and annotated from the Xenium Prime 5K prostate dataset. Macrophages are subclustered and split by GDF15 expression, then compared for spatial neighborhood composition (100 µm radius) and pathway enrichment.

Panel	Content	Script
5a	Macrophage UMAP: GDF15, MRC1 (CD206), CD163, CD68	08_GDF15_macrophage_xenium.R (Sections 7c–7d)
5b	GDF15-high vs GDF15-low/neg macrophages in full tissue; 100 µm neighborhood zooms	08_GDF15_macrophage_xenium.R (Section 8d); 06_neighborhood_analysis_celltype.R (zooms)
5c	Neighborhood enrichment, Δ proportion (GDF15-high − GDF15-low/neg)	06_neighborhood_analysis_celltype.R
5d	GO Biological Process GSEA, GDF15-high vs GDF15-absent macrophages	08_GDF15_macrophage_xenium.R (Section 12a)
Requirements

R ≥ 4.3 with:

CRAN: Seurat (≥ 5.0), ggplot2, dplyr, patchwork, cowplot, viridis, ggrepel, corrplot, Matrix, scales
Bioconductor / other: limma, edgeR, fgsea, msigdbr
Data

10x Genomics Xenium Prime 5K, human prostate adenocarcinoma. Data are not included in this repository.

Scripts read the dataset as cell_feature_matrix/ (MTX) plus cells.csv.gz, produced by 00_convert_xenium_zarr.py from the original zarr bundle. That converter is referenced by the scripts but is not in this repository.

Pipeline

Run in this order.

Step	Script	Output
1	annotate_celltypes_findmarkers.R	Load, QC, cluster, FindAllMarkers, manual cluster-to-cell-type mapping → annot_xenium_celltype_annotated.rds
2	08_GDF15_macrophage_xenium.R	Cell-type annotation, macrophage subclustering, GDF15 status, pseudobulk DE and fgsea → 08_xenium_full_annotated.rds, 08_macrophage_GDF15_xenium.rds
3	06_neighborhood_analysis_celltype.R	Neighborhood enrichment lollipop and zoom panels using the Step 1 annotation
opt.	09_zoom_neighborhood.R	Timestamped zoom panels for random GDF15-high / GDF15-low macrophages (uses 08_* objects)

08_GDF15_macrophage_xenium.R has a RESUME_FROM_RDS toggle. When TRUE, it loads the saved 08_* objects and skips loading, clustering and annotation.

Key parameters
QC: nCount_Xenium ≥ 10, cell area 10–5000 µm²
Normalization: CLR
Clustering: 20 PCs, resolution 0.4
Cell types (Step 1): Epithelial Luminal / Basal / Other, Smooth Muscle, Fibroblast, Endothelial, Macrophage, T/NK, Mast, Schwann/Neural, Unknown
GDF15-high (5a–c): GDF15 above the median of GDF15-expressing macrophages; all others are GDF15-low/neg
Neighborhood (5c): 100 µm radius; equal numbers of GDF15-high and GDF15-low/neg macrophages, set.seed(42)
Pseudobulk (5d): GDF15-expressing (> 0) vs GDF15-absent (= 0) macrophages, summed per macrophage subcluster (≥ 3 cells); limma-voom with a subcluster covariate; genes ranked by signed −log10 p; fgsea on MSigDB GO:BP (human)
Notes
Scripts use hardcoded paths (/Users/alekhyamazumdar/Documents/Xenium/...). Edit output_dir and the data paths before running.
Zoom panels pick random macrophages (clock-seeded in 06 and 09), so the example cells differ between runs.
08 and 09 use the earlier module-score annotation (auto_celltype); 06 uses the marker-based annotation from Step 1 (cell_type).
