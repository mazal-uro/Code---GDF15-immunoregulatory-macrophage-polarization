Gdf15 Expression in Tumor-Infiltrating Immune Cells (TRAMP-C1 scRNA-seq) — Figure 2a–b

R scripts for the reanalysis of public single-cell RNA-seq data from TRAMP-C1 mouse prostate tumors (Cali et al., 2024). The analysis shows that Gdf15 is expressed mainly in the M2-like macrophage clusters (Figure 2a–b). Panels 2c–f (tumor model schematic, Gdf15/Cd206 immunofluorescence, quantification and colocalization) are experimental data and not part of this repository.

Overview
Panel	Content	Script
2a	UMAP of annotated clusters and Gdf15 expression	TRAMPC1_Untreated_Analysis.R
2b	Dot plot of lineage markers and Gdf15 by cluster	TRAMPC1_Untreated_Analysis.R ("DOTPLOT" section)

Only untreated tumors (orig.ident == "UNCTX") are used.

Requirements

R with: Seurat, dplyr, patchwork, tidyverse, data.table, R.utils, cowplot, SCpubr, ggplot2, scDblFinder, SingleCellExperiment

Data

Data are not included in this repository.

Dataset: GEO GSE272999, GSE272999_matrTRAMPC1.txt.gz (genes × cells count matrix, mouse)
Samples in the file: UNCTX, ENZA, ANTIIL23, CXCR2, CTXnoENZA. Only UNCTX is analysed.
Scripts
Script	Purpose
TRAMPC1_Untreated_Analysis.R	Main analysis: load, QC, normalize, cluster, subset UNCTX, markers, panels 2a–b → Seurat_TRAMPC1.rds
TRAMPC1_Doublet_Finder.R	Doublet detection with scDblFinder on the full dataset (see Notes)
Key parameters
Object creation: min.cells = 3, min.features = 200; mitochondrial pattern ^mt-
Normalization: LogNormalize, scale factor 10,000
Variable features: 4,000 (vst)
Full dataset: 20 PCs, clustering resolution 0.3
UNCTX subset: FindClusters at resolution 0.1, UMAP on 15 PCs
Dot plot (2b): clusters 13, 9, 10, 8, 16, 7, 3, 2, 1, shown as DC, NK, Treg, CD8⁺ T, CD4⁺ T, B, neutrophil and two M2 macrophage clusters. Average expression colour scale capped at 2.5.
Marker genes: DC (Zbtb46, Flt3, Xcr1), NK (Ncr1, Klrc1, Eomes), Treg (Foxp3, Il2ra, Ctla4), CD8 (Cd8a, Cd8b1, Runx3), CD4 (Cd4, Sell), B (Cd19, Ms4a1, Cd79a), neutrophil (S100a9, S100a8, Ly6g), M2 macrophage (Mrc1, Cd163, Cx3cr1), plus Gdf15
Usage
Download the GEO file and edit the hardcoded paths in both scripts.
Run TRAMPC1_Untreated_Analysis.R. Plots are printed in the R session; the Seurat object of the UNCTX subset is saved as Seurat_TRAMPC1.rds.
Notes
Cell-type labels in Figure 2a (e.g. "M2 Macrophages", "T-Regs") were assigned manually from the marker dot plots. They are not stored in the object.
Figure panels are printed in R and not saved by the scripts. Export them manually.
The scripts use hardcoded local paths (/Users/alekhyamazumdar/Documents/scRNA/TRAMPC1/).
