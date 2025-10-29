# -----------------------------------------
# STEP 1: Load Unzip h5ad data
# -----------------------------------------

unzip("/Users/alekhyamazumdar/Documents/scRNA/Tabula Sapiens - Prostate/TS_Prostate.h5ad.zip",
      exdir = "/Users/alekhyamazumdar/Documents/scRNA/Tabula Sapiens - Prostate/")

# -----------------------------------------
# STEP 2: Install packages and load libraries
# -----------------------------------------
install.packages("BiocManager")
BiocManager::install(c("zellkonverter", "SingleCellExperiment"))

library(zellkonverter)
library(SingleCellExperiment)
library(Seurat)
library(dplyr)
library(patchwork)
library(R.utils)
library(tidyverse)
library(data.table)
library(cowplot)
library(SCpubr)

# ---------------------------------------------------------------------
# STEP 3: Read the .h5ad file into R and change to a Seurat object
# ---------------------------------------------------------------------
adata <- readH5AD("/Users/alekhyamazumdar/Documents/scRNA/Tabula Sapiens - Prostate/TS_Prostate.h5ad")

# Extract raw counts
counts_matrix <- assays(adata)[["raw_counts"]]

# Extract metadata (optional but useful)
metadata <- as.data.frame(colData(adata))

# Create Seurat object manually
SapiensProstate <- CreateSeuratObject(counts = counts_matrix, meta.data = metadata)

# Optional: store existing dimensional reductions (e.g., scVI UMAP)
umap_coords <- reducedDims(adata)[["X_scvi_umap"]]
SapiensProstate[["umap"]] <- CreateDimReducObject(embeddings = umap_coords, key = "UMAP_")

# -----------------------------------------
# STEP 4: QC and visualization
# -----------------------------------------

# Calculate mitochondrial gene percentage
SapiensProstate[["percent.mt"]] <- PercentageFeatureSet(SapiensProstate, pattern = "^MT-")

# QC violin plots
VlnPlot(SapiensProstate, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Scatter plot: nCount vs nFeature
FeatureScatter(SapiensProstate, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# -----------------------------------------
# STEP 5: Filter low-quality cells
# -----------------------------------------
SapiensProstate <- subset(SapiensProstate, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA < 40000 & percent.mt < 20)
SapiensProstate <- subset(SapiensProstate, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA < 40000)
# -----------------------------------------
# STEP 6: Normalize and find variable features
# -----------------------------------------
SapiensProstate <- NormalizeData(SapiensProstate, normalization.method = "LogNormalize", scale.factor = 10000)
SapiensProstate <- FindVariableFeatures(SapiensProstate, selection.method = "vst", nfeatures = 4000)

# Identify top 10 variable genes and label on plot
top10 <- head(VariableFeatures(SapiensProstate), 10)
plot1 <- VariableFeaturePlot(SapiensProstate)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# -----------------------------------------
# STEP 7: Scale, PCA, and dimension reduction
# -----------------------------------------
all.genes <- rownames(SapiensProstate)
SapiensProstate <- ScaleData(SapiensProstate, features = all.genes)
SapiensProstate <- RunPCA(SapiensProstate, features = VariableFeatures(object = SapiensProstate))

# PCA output and visualizations
print(SapiensProstate[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(SapiensProstate, dims = 1:2, reduction = "pca")
DimPlot(SapiensProstate, reduction = "pca") + NoLegend()
DimHeatmap(SapiensProstate, dims = 1:15, cells = 500, balanced = TRUE)

# -----------------------------------------
# STEP 8: Determine dimensionality and clustering
# -----------------------------------------
ElbowPlot(SapiensProstate)

SapiensProstate <- FindNeighbors(SapiensProstate, dims = 1:15)
SapiensProstate <- FindClusters(SapiensProstate, resolution = 0.1)
SapiensProstate <- RunUMAP(SapiensProstate, dims = 1:15)

# UMAP with cluster labels
DimPlot(SapiensProstate, reduction = "umap", label = TRUE)

# -----------------------------------------
# STEP 9: Visualize Clusters for cell subtypes
# -----------------------------------------
VlnPlot(SapiensProstate, features = c("AR", "KRT19", "CD14", "CD68", "KIT", "GDF15"), pt.size = 0)


#Monolytic Cluster - No 3
FeaturePlot(SapiensProstate, features = c("LYZ", "FCGR3A", "CSF1R", "CD68", "CD163", "CD14", "HAVCR2", "GDF15" ))
FeaturePlot(SapiensProstate, features = c("LYZ" ))

#M2 Macrophage Cluster - No 6
FeaturePlot(SapiensProstate, features = c("MRC1", "CD163", "CD200", "CXCR1", "IL10", "CCL1", "GDF15" ))

#Basal/Intermediate Cluster - 
FeaturePlot(SapiensProstate, features = c("AR", "KRT19", "KRT18", "KRT8", "TP63", "KRT14", "KRT5", "GDF15" ))
FeaturePlot(SapiensProstate, features = c("KRT18"))

#Luminal Cluster
FeaturePlot(SapiensProstate, features = c("AR", "KRT19", "KRT18", "KRT8", "GDF15" ))


# -----------------------------------------
# STEP 10: Subset macrophage cluster (cluster 3)
# -----------------------------------------
Idents(SapiensProstate) <- "seurat_clusters"
macrophages <- subset(SapiensProstate, idents = "3")

# -----------------------------------------
# STEP 11: Re-analyze macrophages
# -----------------------------------------
macrophages <- FindVariableFeatures(macrophages)
macrophages <- ScaleData(macrophages)
macrophages <- RunPCA(macrophages)
macrophages <- FindNeighbors(macrophages, dims = 1:10)
macrophages <- FindClusters(macrophages, resolution = 0.1)
macrophages <- RunUMAP(macrophages, dims = 1:10)
DimPlot(macrophages, reduction = "umap", label = TRUE)

# -----------------------------------------
# STEP 12: Visualize GDF15 in macrophages
# -----------------------------------------

FeaturePlot(macrophages, features = c("IL1B", "TNF", "CD86", "CD163", "MRC1", "IL10", "GDF15"))
VlnPlot(macrophages, features = "GDF15", pt.size = 0.1)
vln_gdf15 <- VlnPlot(
  object = macrophages,
  features = "GDF15",
  pt.size = 0.1,
  layer = "data"
) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.text.x = element_text(face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 14),
    panel.grid.major = element_line(color = "grey85", linetype = "dotted", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  ) +
  labs(
    title = "GDF15 Expression Across Macrophage Clusters",
    x = "Cluster",
    y = "Log Normalized Expression"
  )

print(vln_gdf15)

DotPlot(macrophages, features = "GDF15") + RotatedAxis()

# -----------------------------------------
# STEP 13: Publication-ready DimPlot for macrophage subclusters
# -----------------------------------------

# UMAP by subclusters
p_macrophage_umap <- SCpubr::do_DimPlot(
  sample = macrophages,
  reduction = "umap",
  group.by = "seurat_clusters",
  pt.size = 1.5,
  label = TRUE,
  plot.title = "UMAP of Macrophage Subclusters",
  plot.subtitle = "From Cluster 8 of SapiensTabula Dataset"
)
print(p_macrophage_umap)
# Feature plot of GDF15 with custom color scale
p_macrophage_gdf15 <- SCpubr::do_FeaturePlot(
  sample = macrophages,
  features = "GDF15",
  pt.size = 1.5,
  plot.title = "GDF15 Expression on Macrophage UMAP",
  min.cutoff = 0,
  max.cutoff = 3.0
)
print(p_macrophage_gdf15)

# Combine UMAP and GDF15 expression plot side-by-side
plot_grid(p_macrophage_umap, p_macrophage_gdf15, labels = c("A", "B"), ncol = 2, align = "hv")

# -----------------------------------------
# STEP 14: Dot plot
# -----------------------------------------
# Define your gene list in grouped order
marker_genes <- c(
  # 1. Pan-macrophage / myeloid
  "LYZ", "FCGR3A", "CD68", "CD14", "HAVCR2",
  
  # 2. M1-like
  "IL1B", "TNF", "CXCL10", "CD80", "CD86",
  
  # 3. M2-like
  "MRC1", "IL10", "TGFBI", "CD163", "MSR1", "CSF1R",
  
  # 4. GDF15
  "GDF15"
)

# Extract dotplot data
dot_data <- DotPlot(macrophages, features = marker_genes)$data

# Set GDF15 to highlight
dot_data$gene_group <- ifelse(dot_data$features.plot == "GDF15", "GDF15", "Other")

# Preserve gene order
dot_data$features.plot <- factor(dot_data$features.plot, levels = marker_genes)

# Get cluster count for vertical box
cluster_count <- length(unique(dot_data$id))

# Plot
ggplot(dot_data, aes(x = features.plot, y = id, size = pct.exp, color = gene_group)) +
  geom_point(aes(alpha = avg.exp.scaled)) +
  
  # Group dividers between gene sets
  geom_vline(xintercept = c(5.5, 10.5, 16.5), linetype = "dashed", color = "grey40", linewidth = 0.75) +
  
  # Box around the gene axis
  annotate("segment", x = 0.5, xend = 17.5, y = 0.5, yend = 0.5, color = "grey40", linewidth = 0.75) +  # bottom
  annotate("segment", x = 0.5, xend = 17.5, y = cluster_count + 0.5, 
           yend = cluster_count + 0.5, color = "grey40", linewidth = 0.75) +  # top
  annotate("segment", x = 0.5, xend = 0.5, y = 0.5, yend = cluster_count + 0.5, 
           color = "grey40", linewidth = 0.75) +  # left
  annotate("segment", x = 17.5, xend = 17.5, y = 0.5, yend = cluster_count + 0.5, 
           color = "grey40", linewidth = 0.75) +  # right
  
  # Aesthetic controls
  scale_size(range = c(1, 6)) +
  scale_color_manual(values = c("GDF15" = "firebrick", "Other" = "darkblue")) +
  scale_alpha(range = c(0.2, 1), guide = "none") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Macrophage Marker Expression by Subcluster",
    x = "Genes",
    y = "Macrophage Subclusters",
    size = "% Expressing",
    color = "Marker"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.title = element_text(face = "bold"),
    
    # Transparent background grid
    panel.grid.major = element_line(color = "grey85", linetype = "dotted", linewidth = 0.3),
    panel.grid.minor = element_line(color = "grey92", linetype = "dotted", linewidth = 0.2)
  )

# -----------------------------------------
# STEP 15: Save the processed Seurat object
# -----------------------------------------
setwd("/Users/alekhyamazumdar/Documents/scRNA/Tabula Sapiens - Prostate")
saveRDS(macrophages, "SapiensProstateMacrophages.rds")
save.image(file = "SapiensProstate.RData")
