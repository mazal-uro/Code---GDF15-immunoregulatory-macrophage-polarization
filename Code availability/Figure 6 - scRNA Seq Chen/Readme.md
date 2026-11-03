GDF15⁺ Macrophages in Human Prostate Cancer and Healthy Prostate (scRNA-seq) — Figure 6

R scripts for the integrated single-cell analysis of macrophages from human prostate cancer (Chen et al., GSE141445) and healthy prostate (Tabula Sapiens). The analysis identifies a GDF15-high macrophage population and tests its pathway signature (Figure 6a–g).

Overview

Macrophages are isolated from each dataset, integrated, and clustered. GDF15 expression defines a macrophage cluster found mainly in tumor. GDF15⁺ vs GDF15⁻ macrophages are then compared by differential expression and GSEA, focusing on SMAD signalling and inflammatory response.

Panel	Content	Script
6a	Integrated UMAP by cluster and by dataset	Integrate_Macrophages_v2.R (Step 5)
6b	GDF15 UMAP, Chen vs Tabula Sapiens	Integrate_Macrophages_v2.R (Step 5)
6c	GDF15 ridge plot per cluster	Integrate_Macrophages_v2.R (Step 12)
6d	Pan / M1 / M2 marker dot plot and top cluster markers	Integrate_Macrophages_v2.R (Step 12)
6e	GO GSEA lollipop, GDF15⁺ vs GDF15⁻ (within Chen)	GO_Inflammatory_decomposition_GDF15.R (genome-wide GO lollipop)
6f	GSEA: positive regulation of inflammatory response (GO:0050729)	GO_Inflammatory_decomposition_GDF15.R
6g	GSEA: regulation of SMAD signalling (GO:0060390)	GO_Inflammatory_decomposition_GDF15.R

GO_SMAD_decomposition_GDF15.R (positive vs negative arms of SMAD regulation) and GDF15_supplementary_figures.R (GDF15 distribution, QC, volcano, module scores) are supporting analyses.

Requirements

R with: Seurat (v5), dplyr, ggplot2, patchwork, cowplot, ggrepel, ggridges, SCpubr, pheatmap, tibble, data.table, clusterProfiler, enrichplot, org.Hs.eg.db, msigdbr, zellkonverter, SingleCellExperiment

Data

Data are not included in this repository.

Chen et al.: GSE141445, sample GSM4203181 (GSM4203181_data.matrix.txt.gz). The matrix is already normalized and has no mitochondrial genes.
Tabula Sapiens, prostate: TS_Prostate.h5ad. Raw integer counts; mitochondrial genes present.
Pipeline

Run in order.

Step	Script	Output
1	Chen_QC_Macrophages.R	QC, clustering, macrophage selection → Seurat_ChenEtAlmacrophages.rds
2	Sapiens_QC_Macrophages.R	QC, clustering, macrophage selection → SapiensProstateMacrophages.rds
3	Integrate_Macrophages_v2.R	Gene harmonization, RPCA integration, clustering, panels 6a–d → IntegratedMacrophages_v2.rds
4	GO_Inflammatory_decomposition_GDF15.R	GDF15⁺ vs GDF15⁻ DE and GSEA, panels 6e–f
5	GO_SMAD_decomposition_GDF15.R	SMAD GO-term decomposition, panel 6g
opt.	GDF15_supplementary_figures.R	Supplementary figures A–F
Key parameters
QC (both datasets): MAD-based outlier filtering, 3 MADs (fixed thresholds available via qc_method)
Macrophage selection: clusters with high mean module score for canonical markers (CD68, CD14, CSF1R, C1QA/B/C, MRC1, CD163, etc.), mean + 1 SD across clusters; check and override manually if needed
Integration: common genes only (removes mitochondrial genes absent from Chen), LogNormalize, 3000 variable features, RPCA anchors, 30 dims
Clustering: 20 PCs, resolution 0.1 (cluster 1 = GDF15-high)
GDF15 status: GDF15⁺ = log-normalized GDF15 > 0, within Chen only (avoids dataset confound), independent of clusters
DE / GSEA: FindMarkers (GDF15pos vs GDF15neg, no logFC threshold, min.pct 0.05); genes ranked by avg_log2FC; clusterProfiler GSEA, seed 42. GDF15 is removed from the SMAD and inflammatory gene sets to avoid circularity
Notes
Scripts use hardcoded paths (/Users/alekhyamazumdar/Documents/scRNA/...). Edit them before running.
Panels 6e–g are GDF15⁺ vs GDF15⁻ within Chen. The cluster-vs-rest GSEA in Integrate_Macrophages_v2.R (Step 13) is a separate comparison.
Figure panels are printed in R and not saved by the scripts. Export them manually.
