Flow Cytometry UMAPs — Figure 3i

R scripts for the UMAP analysis of CD45⁺ immune cells from tumor-bearing and control mice (Figure 3i): organ-wise UMAP clustering, tumor EV fluorescence (mEmerald), and CD64 / CD206 expression overlays.

Overview

Spectral flow cytometry data (Cytek Aurora, unmixed, FlowJo-compensated FJComp-* channels) from control and tumor-bearing mice are pooled across organs, gated on CD45, downsampled, and embedded with UMAP. Marker expression is overlaid per organ and condition.

Organs: bone marrow, liver, lungs, lymph nodes (contralateral and draining), spleen, primary tumor (tumor mice only) Conditions: Control, Tumor

Requirements

R with: flowCore, ggplot2, dplyr, Rtsne, umap, svglite (for SVG export)

Pipeline

Run in order, in the same R session (later scripts use objects created by earlier ones).

Step	Script	Purpose
1	STEP1_-_Loading_FCS_Files_and_Tsne.R	Load FCS files, label organ and condition, select fluorescence channels, run t-SNE (QC)
1 (alt.)	RScriptTsne.R	Same as Step 1 with the extended sample list (contralateral/draining lymph nodes, PC3 tumor)
2	STEP2_-_Convert_TSNE_to_UMAP.R	Gate CD45⁺ cells, downsample, compute UMAP, plot by organ × condition
3	UMAP_<marker>.R	Overlay marker expression on the UMAP and export SVG
–	ViolinPlots.R, Expression_distribution_plots.R	Expression distributions (supporting)

UMAP_script.R is identical to STEP2_-_Convert_TSNE_to_UMAP.R.

Processing parameters
Gate: CD45 (FJComp-APC-Cy7-A) > 50,000
Downsampling: up to 10,000 cells per organ × condition, set.seed(42)
Input channels (10): CD45, MHCII, CD64, CD11b, CD103, CD206, SiglecF, Ly6G, CD11c, mEmerald
t-SNE: perplexity 50, theta 0.5, seed 42
UMAP: umap package defaults
Marker–channel map
Marker	Channel	Script
CD45	FJComp-APC-Cy7-A	UMAP_CD45.R
MHCII	FJComp-Alexa Fluor 700-A	UMAP_MHCII.R
CD64	FJComp-BV421-A	UMAP_CD64.R
CD11b	FJComp-BV510-A	UMAP_CD11b.R
CD103	FJComp-BV605-A	UMAP_CD103.R
CD206	FJComp-BV650-A	UMAP_CD206.R
SiglecF	FJComp-PE-CF594-A	UMAP_SigF.R
Ly6G	FJComp-PerCP-Cy5.5-A	UMAP_Ly6G.R
CD11c	FJComp-PE-Cy7-A	CD11c_UMAP.R
mEmerald (tumor EVs)	FJComp-mEmerald-A	UMAP_mEmerald.R

Figure 3i panels: UMAP clustering (Step 2), tumor EV fluorescence (UMAP_mEmerald.R), CD64 (UMAP_CD64.R), CD206 (UMAP_CD206.R).

Usage
Edit setwd() and the FCS file names in Step 1 to match your data.
Run Step 1, then Step 2.
Run the UMAP_<marker>.R script for the marker of interest. Each saves UMAP_<marker>_Fixed.svg (4.5 × 6 in).
Notes
Marker overlays show only the top-expressing cells (percentile cutoff set per script); all other cells are shown in gray.
Colour scale limits are set manually per marker.
FCS files are not included in this repository.
