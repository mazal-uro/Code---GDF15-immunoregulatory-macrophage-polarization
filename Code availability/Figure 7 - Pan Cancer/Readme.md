GDF15⁺ Macrophages Across Cancers (Pan-Cancer scRNA-seq) — Figure 7

R scripts for the pan-cancer analysis of GDF15-expressing macrophages (Figure 7). Within each cancer type, GDF15⁺ and GDF15⁻ macrophages are compared by differential expression and GSEA, testing whether the inflammatory-down / SMAD-up signature seen in prostate cancer is shared across tumors.

Overview
Panel	Content	Script
7a	GDF15 expression on macrophage UMAPs, per cancer	<script>
7b	Macrophage marker dot plot across cancers	<script>
7c	% GDF15⁺ macrophages per cancer	<script>
7d	GSEA running score: positive regulation of inflammatory response (GO:0050729), ESCA / THCA / UCEC / OV-FTC	GDF15_GSEA_lollipop_Pan_Cancer_pipeline.R
7e	GSEA running score: regulation of SMAD protein signal transduction (GO:0060390), same four cancers	GDF15_GSEA_lollipop_Pan_Cancer_pipeline.R
Requirements

R with: Seurat (v5), Matrix, data.table, dplyr, ggplot2, ggrepel, scales, clusterProfiler, enrichplot, org.Hs.eg.db, AnnotationDbi, GO.db, msigdbr

Data

Data are not included in this repository.

Pan-cancer myeloid atlas: GSE154763 (Cheng et al., Cell 2021). Per cancer type, two files in DATA_DIR:
GSE154763_<DS>_normalized_expression.csv.gz (cells × genes, log-normalized)
GSE154763_<DS>_metadata.csv.gz (includes MajorCluster, UMAP1, UMAP2)
Melanoma macrophages (MELMAC): Sade-Feldman et al. melanoma macrophages, prepared by melanoma_macrophage_pipeline.R. Run that first if the MELMAC object is missing.

Cancer types analysed: ESCA, KIDNEY, MYE, OV-FTC, PAAD, THCA, UCEC, MELMAC. cDC2 is excluded by design, and LYM is excluded because GDF15 is not in that matrix.

Usage
Edit DATA_DIR at the top of GDF15_GSEA_lollipop_Pan_Cancer_pipeline.R.
Run:
bash
   Rscript GDF15_GSEA_lollipop_Pan_Cancer_pipeline.R

or source() it interactively.

Outputs go to <DATA_DIR>/Analysis_Output/GDF15_GSEA/ as PNG and PDF where applicable:

Output	Content
DE_<ds>.csv	GDF15⁺ vs GDF15⁻ FindMarkers table
Volcano_<ds>	Volcano plot
GSEA_GOBP_<ds>.csv	Full GO:BP GSEA result
Lollipop_GOBP_<ds>	NES lollipop, top 15 terms per direction
Running_GOBP_<ds>_GO<id>.png	Running-score plots for GO:0050729 and GO:0060390 (panels 7d–e)
Heatmap_NES_key_terms	NES for the two key terms across cancers
S1_QC_GDF15pos_vs_neg, S2_per_patient_fraction, S3_threshold_sensitivity	Robustness checks: QC metrics, per-patient fraction, GDF15 cutoff sensitivity
Key parameters
Macrophages: cells whose MajorCluster contains "Macro"
GDF15 status: GDF15⁺ = log-normalized expression > 0, within macrophages only; independent of clustering. Cancers with fewer than 10 cells in either group are skipped
DE: FindMarkers, min.pct = 0.10, no logFC threshold
Ranking: signed −log10 p (RANK_BY = "signed_p"); GDF15 removed from the ranked list to avoid circularity
GSEA: GO:BP, gene-set size 10–500, clusterProfiler, seed 42
Optional: DO_POOLED <- TRUE also runs one pooled analysis across all cancers (default off)
Notes
DATA_DIR is a hardcoded local path. Edit it before running.
Cancers are grouped in the supplementary figures by whether they show the prostate-like pattern (inflammatory down, SMAD up: ESCA, UCEC, OV-FTC, THCA) or not (MYE, MELMAC, KIDNEY, PAAD).
