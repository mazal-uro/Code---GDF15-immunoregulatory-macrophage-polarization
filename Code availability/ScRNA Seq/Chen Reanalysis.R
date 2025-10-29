# Load necessary libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(R.utils)
library(tidyverse)
library(data.table)
library(cowplot)
library(SCpubr)

# -----------------------------------------
# STEP 1: Load the Chen dataset
# -----------------------------------------

file_path <- "~/Documents/scRNA/Chen et al/GSE141445_RAW/GSM4203181_data.matrix.txt.gz"
PCA_Chen <- read.delim(file_path, header=TRUE, sep='\t')
#PCA_Chen <- fread(file_path, header=TRUE, sep='\t')
PCA_Chen <- as.matrix(PCA_Chen)

# -----------------------------------------
# STEP 2: Create Seurat object
# -----------------------------------------

PCA_Chen <- CreateSeuratObject(counts = PCA_Chen, project = "Chen_Data", min.cells = 3, min.features = 200)
class(PCA_Chen)

# -----------------------------------------
# STEP 3: QC and visualization
# -----------------------------------------

PCA_Chen[["percent.mt"]] <- PercentageFeatureSet(PCA_Chen, pattern = "^mt-")
VlnPlot(PCA_Chen, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
FeatureScatter(PCA_Chen, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# -----------------------------------------
# STEP 4: Filter low-quality cells
# -----------------------------------------
PCA_ChenSub <- subset(PCA_Chen, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 
                      & nCount_RNA < 40000 & percent.mt < 10)

# -----------------------------------------
# STEP 5: Normalize and find variable features
# -----------------------------------------
PCA_ChenSub <- NormalizeData(PCA_ChenSub, normalization.method = "LogNormalize", scale.factor = 10000)
PCA_ChenSub <- FindVariableFeatures(PCA_ChenSub, selection.method = "vst", nfeatures = 4000)

#Identify the 10 most highly variable genes and plot them
top10 <- head(VariableFeatures(PCA_ChenSub), 10)
genes_to_label <- unique(c(top10, "GDF15"))  # add GDF15 if not already in top10

plot1 <- VariableFeaturePlot(PCA_ChenSub)
plot2 <- LabelPoints(plot = plot1, points = genes_to_label, repel = TRUE)
plot2
#Error- log-10 transformation introduced infinite values.

# -----------------------------------------
# STEP 6: Scale, PCA, and dimension reduction
# -----------------------------------------
all.genes <- rownames(PCA_ChenSub)
PCA_ChenSub <- ScaleData(PCA_ChenSub, features = all.genes)

PCA_ChenSub <- RunPCA(PCA_ChenSub, features = VariableFeatures(object = PCA_ChenSub))
print(PCA_ChenSub[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(PCA_ChenSub, dims = 1:2, reduction = "pca")
DimPlot(PCA_ChenSub, reduction = "pca") + NoLegend()
DimHeatmap(PCA_ChenSub, dims = 1:15, cells = 500, balanced = TRUE)

# -----------------------------------------
# STEP 7: Determine dimensionality and clustering
# -----------------------------------------
ElbowPlot(PCA_ChenSub)
PCA_ChenSub <- FindNeighbors(PCA_ChenSub, dims = 1:15)
PCA_ChenSub <- FindClusters(PCA_ChenSub, resolution = 0.1)
PCA_ChenSub <- RunUMAP(PCA_ChenSub, dims = 1:15)
DimPlot(PCA_ChenSub, reduction = "umap", label = TRUE)

# -----------------------------------------
# STEP 8: Identify cluster markers
# -----------------------------------------
cluster_markers <- FindAllMarkers(PCA_ChenSub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10_genes <- cluster_markers %>% group_by(cluster) %>% slice_max(n = 10, order_by = avg_log2FC)
print(top10_genes)

# -----------------------------------------
# STEP 9: Visualize marker genes
# -----------------------------------------
VlnPlot(PCA_ChenSub, features = c("MPO", "AR", "CD1A", "IL5RA", "CD68", "KIT", "MRC1", "CD14", "PTPRC", "KLRK1", "GDF15"), pt.size = 0)
VlnPlot(PCA_ChenSub, features = c("AR", "KRT19", "CD14", "CD68", "KIT", "GDF15"), pt.size = 0)
FeaturePlot(PCA_ChenSub, features = c("CD14", "GDF15"))
FeaturePlot(PCA_ChenSub, features = c("LYZ", "FCGR3A", "CSF1R", "CD68", "CD163", "CD14", "AR", "KRT19", "GDF15"))
DotPlot(PCA_ChenSub, features = c("LYZ", "FCGR3A", "CSF1R", "CD68", "CD163", "CD14", "GDF15")) + RotatedAxis()
#------------------------------------------------------------------------------------------------------------------------
#Tcell Cluster - No 2
FeaturePlot(PCA_ChenSub, features = c("CD3E", "CD2E", "CD7", "IL7R", "GDF15" ))
FeaturePlot(PCA_ChenSub, features = c("FOLH1"))

#Monolytic Cluster - No 9
FeaturePlot(PCA_ChenSub, features = c("LYZ", "FCGR3A", "CSF1R", "CD68", "CD163", "CD14", "HAVCR2", "GDF15" ))
FeaturePlot(PCA_ChenSub, features = c("LYZ" ))

#M2 Macrophage Cluster - No 7
FeaturePlot(PCA_ChenSub, features = c("MRC1", "CD163", "CD200", "CXCR1", "IL10", "CCL1", "GDF15" ))

#Basal/Intermediate Cluster - No 6
FeaturePlot(PCA_ChenSub, features = c("AR", "KRT19", "KRT18", "KRT8", "TP63", "KRT14", "KRT5", "GDF15" ))
FeaturePlot(PCA_ChenSub, features = c("KRT18"))

#Luminal Cluster
FeaturePlot(PCA_ChenSub, features = c("AR", "KRT19", "KRT18", "KRT8", "GDF15" ))

# -----------------------------------------
# STEP 10: Subset macrophage cluster (cluster 7)
# -----------------------------------------
Idents(PCA_ChenSub) <- "seurat_clusters"
macrophages <- subset(PCA_ChenSub, idents = "7")

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
DotPlot(macrophages, features = "GDF15") + RotatedAxis()

# -----------------------------------------
# STEP 13: Find markers for macrophage subclusters
# -----------------------------------------
macrophage_markers <- FindAllMarkers(macrophages, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5_macrophage <- macrophage_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(macrophages, features = top5_macrophage$gene) + NoLegend()

# -----------------------------------------
# STEP 14: Subset GDF15-positive subcluster (e.g., cluster 1)
# -----------------------------------------
GDF15_positive_macs <- subset(macrophages, idents = "1")
gdf15_cluster_markers <- macrophage_markers %>% filter(cluster == 1)
print(gdf15_cluster_markers)

# -----------------------------------------
# STEP 15: Annotate macrophage cluster in full dataset
# -----------------------------------------
PCA_ChenSub$celltype <- "Other"
PCA_ChenSub$celltype[PCA_ChenSub$seurat_clusters == "7"] <- "Macrophage"
DimPlot(PCA_ChenSub, group.by = "celltype", label = TRUE, repel = TRUE)


# -----------------------------------------
# STEP 17: Replot for better visuals
# -----------------------------------------

DimPlot(PCA_ChenSub, reduction = "umap", group.by = "seurat_clusters", pt.size = 1.5, label = FALSE) +
  labs(title = "UMAP of Chen Dataset", subtitle = "Clusters Colored by Seurat Clustering") +
  theme_minimal(base_size = 08) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )

p1 <- SCpubr::do_DimPlot(
  sample = PCA_ChenSub,
  reduction = "umap",
  group.by = "seurat_clusters",
  pt.size = 0.5,
  label = FALSE,
  plot.title = "UMAP of Chen Dataset",
  plot.subtitle = "Clusters Colored by Seurat Clustering"
)

p2 <- SCpubr::do_FeaturePlot(
  sample = PCA_ChenSub,
  features = "GDF15",
  pt.size = 0.5,
  plot.title = "GDF15 Expression on UMAP"
)


plot_grid(p1, p2, ncol = 2, align = "hv")

# -----------------------------------------
# STEP 18: Publication-ready DimPlot for macrophage subclusters
# -----------------------------------------

# UMAP by subclusters
p_macrophage_umap <- SCpubr::do_DimPlot(
  sample = macrophages,
  reduction = "umap",
  group.by = "seurat_clusters",
  pt.size = 1.5,
  label = TRUE,
  plot.title = "UMAP of Macrophage Subclusters",
  plot.subtitle = "From Cluster 8 of Chen Dataset"
)
print(p_macrophage_umap)
# Feature plot of GDF15
p_macrophage_gdf15 <- SCpubr::do_FeaturePlot(
  sample = macrophages,
  features = "GDF15",
  pt.size = 1.5,
  plot.title = "GDF15 Expression on Macrophage UMAP"
)
print(p_macrophage_gdf15)

# Combine UMAP and GDF15 expression plot side-by-side
plot_grid(p_macrophage_umap, p_macrophage_gdf15, labels = c("A", "B"), ncol = 2, align = "hv")

# -----------------------------------------
# STEP 19: Dot plot
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
# STEP 16: Save the processed Seurat object
# -----------------------------------------
setwd("/Users/alekhyamazumdar/Documents/scRNA/Chen et al")
saveRDS(macrophages, "Seurat_ChenEtAlmacrophages.rds")
save.image(file = "Seurat_ChenEtAl.RData")

load("Seurat_ChenEtAl.RData")
