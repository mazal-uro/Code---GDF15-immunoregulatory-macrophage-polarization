library(dplyr)
library(Seurat)
library(patchwork)
library(R.utils)
library(tidyverse)
library(data.table)
library(cowplot)
library(SCpubr)

# Step 1: Load the TRAMPC1 dataset into R and converting it to a dataframe with correct rows and columns
TRAMPC1 <- fread("/Users/alekhyamazumdar/Documents/scRNA/TRAMPC1/GSE272999_matrTRAMPC1.txt.gz")
class(TRAMPC1)
GENE <- TRAMPC1[["names"]]
TRAMPC1 <- as.data.frame(TRAMPC1)
rownames(TRAMPC1)<- GENE
TRAMPC1[["V1"]] <- NULL
TRAMPC1[["names"]] <- NULL

# Step 2: Create Array Matrix from Dataframe
TRAMPC1 <- as.matrix(TRAMPC1)
class(TRAMPC1)

# Step 3: Creating a Seurat object
TRAMPC1 <- CreateSeuratObject(counts = TRAMPC1, project = "TRAMPC1_Data", min.cells = 3, min.features = 200)
class(TRAMPC1)

#Step 4: PreProcessing steps - QC and selecting cells + visualization
TRAMPC1[["percent.mt"]] <- PercentageFeatureSet(TRAMPC1, pattern = "^mt-")
VlnPlot(TRAMPC1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(TRAMPC1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0)
plot1 <- FeatureScatter(TRAMPC1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1

# Step 5: Not Subsetting Cells to another variable
#TRAMPC1 <- subset(TRAMPC1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 
              #& nCount_RNA < 40000)

# Step 6: Normailize data
TRAMPC1 <- NormalizeData(TRAMPC1, normalization.method = "LogNormalize", scale.factor = 10000)

#Step 7: Identification of highly variable features (feature selection)
TRAMPC1 <- FindVariableFeatures(TRAMPC1, selection.method = "vst", nfeatures = 4000)

#Step 8: Identify the 10 most highly variable genes and plot them
top10 <- head(VariableFeatures(TRAMPC1), 10)
plot1 <- VariableFeaturePlot(TRAMPC1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
#top10 <- head(VariableFeatures(TRAMPC1), 10)
#genes_to_label <- unique(c(top10, "Gdf15"))  # ensures no duplicates
#plot1 <- VariableFeaturePlot(TRAMPC1)
#plot2 <- LabelPoints(plot = plot1, points = genes_to_label, repel = TRUE)
#plot2


#Step 9: Scaling the data so that the mean expression across cells is 0
all.genes <- rownames(TRAMPC1)
TRAMPC1 <- ScaleData(TRAMPC1, features = all.genes)

#Step 10: Perform linear dimensional reduction and vizualization
TRAMPC1 <- RunPCA(TRAMPC1, features = VariableFeatures(object = TRAMPC1))
print(TRAMPC1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(TRAMPC1, dims = 1:2, reduction = "pca")
DimPlot(TRAMPC1, reduction = "pca") + NoLegend()
DimHeatmap(TRAMPC1, dims = 1:15, cells = 500, balanced = TRUE)

#Step 11: Determine the ‘dimensionality’ of the dataset
ElbowPlot(TRAMPC1)

#Step 12: Cluster the cells and Visualize UMAP
TRAMPC1 <- FindNeighbors(TRAMPC1, dims = 1:20)
TRAMPC1 <- FindClusters(TRAMPC1, resolution = 0.3)
TRAMPC1 <- RunUMAP(TRAMPC1, dims = 1:20)
DimPlot(TRAMPC1, reduction = "umap", group.by = "orig.ident")

# Define identities
# TRAMPC1$highlight <- ifelse(TRAMPC1$orig.ident == "UNCTX", "UNCTX", "Other")

# Plot with manual colors
DimPlot(TRAMPC1, reduction = "umap", group.by = "highlight") +
 scale_color_manual(values = c("UNCTX" = "red", "Other" = "lightgray")) +
  ggtitle("Highlighting UNCTX") +
  theme(plot.title = element_text(hjust = 0.5))

FeaturePlot(TRAMPC1, features = "Gdf15", reduction = "umap", 
            cols = c("lightgray", "red")) +
  ggtitle("Gdf15 Expression") +
  theme(plot.title = element_text(hjust = 0.5))

FeaturePlot(TRAMPC1, features = "Mrc1", reduction = "umap", 
            cols = c("lightgray", "red")) +
  ggtitle("MRC1 Expression") +
  theme(plot.title = element_text(hjust = 0.5))

VlnPlot(TRAMPC1, features = "Gdf15", group.by = "orig.ident", pt.size = 0) +
  ggtitle("Gdf15 Expression by Sample") +
  theme(plot.title = element_text(hjust = 0.5))

VlnPlot(TRAMPC1, features = "Mrc1", group.by = "orig.ident", pt.size = 0) +
  ggtitle("Mrc1 Expression by Sample") +
  theme(plot.title = element_text(hjust = 0.5))


# GDF15 Expression


# Assuming your Seurat object is named 'seurat_obj'
# Specify the orig.ident values you want to keep
idents_to_keep <- c("UNCTX", "ENZA", "ANTIIL23", "CXCR2", "CTXnoENZA")
idents_to_keep <- c("UNCTX")

# Subset the Seurat object

TRAMPC1Sub<- subset(TRAMPC1, subset = orig.ident %in% idents_to_keep)
TRAMPC1Sub <- FindClusters(TRAMPC1Sub, resolution = 0.1)
TRAMPC1Sub <- RunUMAP(TRAMPC1Sub, dims = 1:15)
DimPlot(TRAMPC1Sub, reduction = "umap", group.by = "orig.ident")
DimPlot(TRAMPC1Sub, reduction = "umap", pt.size = 1)
rm(TRAMPC1)

#Find all markers
TRAMPC1Sub.markers <- FindAllMarkers(TRAMPC1Sub, only.pos = TRUE)
 TRAMPC1Sub.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# Subset the top 5 markers per cluster
top_markers <- TRAMPC1Sub.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# Generate a heatmap for the top markers
DoHeatmap(TRAMPC1Sub, features = top_markers$gene) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme_minimal() +
  labs(title = "Heatmap of Top Upregulated Markers", x = "Cells", y = "Genes")

DimPlot(TRAMPC1Sub, reduction = "umap", group.by = "seurat_clusters", pt.size = 1.5, label = FALSE, label.size = 5) +
  labs(title = "UMAP of UNCTX Cells", subtitle = "Colored by Seurat Clusters") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )

SCpubr::do_DimPlot(
  sample = TRAMPC1Sub,
  reduction = "umap",
  group.by = "seurat_clusters",
  pt.size = 1.5,
  label = FALSE,
  plot.title = "UMAP of UNCTX Cells",
  plot.subtitle = "Colored by Seurat Clusters"
)

library(cowplot)

p1 <- SCpubr::do_DimPlot(sample = TRAMPC1Sub, group.by = "seurat_clusters")
p2 <- SCpubr::do_FeaturePlot(
  sample = TRAMPC1Sub,
  features = "Xcr1",
  pt.size = 1.5,
  plot.title = "Gdf15 Expression on UMAP"
)
plot_grid(p2)
plot_grid(p1, p2, labels = c("A", "B"), ncol = 2, align = "hv")


#ALL GENES --------------------------------------------------------------------

#Dendtitic Marker-13
DotPlot(TRAMPC1Sub, features = c("Zbtb46", "Flt3", "Xcr1"), dot.scale = 5)

#NK Marker-9
DotPlot(TRAMPC1Sub, features = c("Ncr1", "Klrc1", "Eomes"), dot.scale = 5)

#TReg-10
DotPlot(TRAMPC1Sub, features = c("Foxp3", "Il2ra", "Ctla4"), dot.scale = 5)

#Cd8-8
DotPlot(TRAMPC1Sub, features = c("Cd8a", "Cd8b1", "Runx3"), dot.scale = 5)
FeaturePlot(TRAMPC1Sub, features = c("Cd8a"), pt.size = 1)

#Cd4-16/5
DotPlot(TRAMPC1Sub, features = c("Cd4", "Sell"), dot.scale = 5)
FeaturePlot(TRAMPC1Sub, features = c("Cd4"), pt.size = 1)

#Bcells-7
DotPlot(TRAMPC1Sub, features = c("Cd19", "Ms4a1", "Cd79a"), dot.scale = 5)
FeaturePlot(TRAMPC1Sub, features = c("Cd4"), pt.size = 1)

#Neutrophils-3
DotPlot(TRAMPC1Sub, features = c("S100a9", "S100a8", "Ly6g"), dot.scale = 5)
FeaturePlot(TRAMPC1Sub, features = c("Ly6g"), pt.size = 1)



#M2 Macrophages-1,2
DotPlot(TRAMPC1Sub, features = c("Mrc1", "Cd163", "Cx3cr1"), dot.scale = 5)
FeaturePlot(TRAMPC1Sub, features = c("Mrc1"), pt.size = 1)

#GDF15
DotPlot(TRAMPC1Sub, features = c("Gdf15"), dot.scale = 5)



#Monocytes-3
DotPlot(TRAMPC1Sub, features = c("Fcgr3", "Nr4a1", "Cx3cr1", "Ccr2", "Ly6c2", "Cd14", "Trem1", "Plac8", "Itgam"), dot.scale = 5)
DotPlot(TRAMPC1Sub, features = c("Fcgr3", "Cd14", "Trem1"), dot.scale = 5)
FeaturePlot(TRAMPC1Sub, features = c("Fcgr3", "Cd14", "Trem1", ), pt.size = 1)

#Neuronal Marker
FeaturePlot(TRAMPC1Sub, features = c("Nrxn1"), pt.size = 1)

#SMC Marker
FeaturePlot(TRAMPC1Sub, features = c("Myh11"))
FeaturePlot(TRAMPC1Sub, features = c("Tagln"))


#AR Cell Marker
FeaturePlot(TRAMPC1Sub, features = c("Ar"), pt.size = 1)

#Endothelial cell Marker
FeaturePlot(TRAMPC1Sub, features = c("Vwf"), pt.size = 1)

#B cell Marker
FeaturePlot(TRAMPC1Sub, features = c("Ms4a1"), pt.size = 1)
FeaturePlot(TRAMPC1Sub, features = c("Cd19"), pt.size = 1)

#Epithelial Markers
FeaturePlot(TRAMPC1Sub, features = c("Cdh1"), pt.size = 1)
FeaturePlot(TRAMPC1Sub, features = c("Pbsn"), pt.size = 1)

#T cell Marker
FeaturePlot(TRAMPC1Sub, features = c("Il2ra"), pt.size = 1)
FeaturePlot(TRAMPC1Sub, features = c("Cd3d"), pt.size = 1)
FeaturePlot(TRAMPC1Sub, features = c("Cd4"), pt.size = 1)
FeaturePlot(TRAMPC1Sub, features = c("Foxp3"), pt.size = 1)
FeaturePlot(TRAMPC1Sub, features = c("Cd8a"), pt.size = 1)

#NK cell Marker
FeaturePlot(TRAMPC1Sub, features = c("Ncr1"), pt.size = 1)
FeaturePlot(TRAMPC1Sub, features = c("Klrb1c"), pt.size = 1)

#Monocytes cell Marker
FeaturePlot(TRAMPC1Sub, features = c("Cx3cr1"), pt.size = 1)
FeaturePlot(TRAMPC1Sub, features = c("Ccr2"), pt.size = 1)

#Macrophages cell Marker
FeaturePlot(TRAMPC1Sub, features = c("Adgre1"), pt.size = 1)
FeaturePlot(TRAMPC1Sub, features = c("Mrc1"), pt.size = 1)

#M1 Marker
FeaturePlot(TRAMPC1Sub, features = c("Nos2"), pt.size = 1)

#DC Marker
FeaturePlot(TRAMPC1Sub, features = c("Xcr1"), pt.size = 1)

#Neutrophil Marker
FeaturePlot(TRAMPC1Sub, features = c("S100a8"), pt.size = 1)
FeaturePlot(TRAMPC1Sub, features = c("Ly6g"), pt.size = 1)

#Mast Cells
FeaturePlot(TRAMPC1Sub, features = c("Cpa3"), pt.size = 1)
FeaturePlot(TRAMPC1Sub, features = c("Fcer1a"), pt.size = 1)

#Macrophage DC Marker
FeaturePlot(TRAMPC1Sub, features = c("Fn1"), pt.size = 1)
FeaturePlot(TRAMPC1Sub, features = c("Slc37a2"), pt.size = 1)
FeaturePlot(TRAMPC1Sub, features = c("Des"), pt.size = 1)
FeaturePlot(TRAMPC1Sub, features = c("Ptprc"), pt.size = 1)

#MDSCs Marker
FeaturePlot(TRAMPC1Sub, features = c("Arg1"), pt.size = 1)

#Platelets Marker
FeaturePlot(TRAMPC1Sub, features = c("Hba-a1"), pt.size = 1)

#Endo Marker
FeaturePlot(TRAMPC1Sub, features = c("Vwf"), pt.size = 1)

#LEC Marker
FeaturePlot(TRAMPC1Sub, features = c("Elf3"), pt.size = 1)

#MED16 Marker
FeaturePlot(TRAMPC1Sub, features = c("Med16"), pt.size = 1)

#GDF15 Marker and M2 Macrophages
FeaturePlot(TRAMPC1Sub, features = c("Gdf15"), pt.size = 1)
FeaturePlot(TRAMPC1Sub, features = c("Mrc1"), pt.size = 1)
FeaturePlot(TRAMPC1Sub, features = c("Il10"), pt.size = 1)

FeaturePlot(TRAMPC1Sub, features = c("Mrc1", "Cd163", "Cd200", "Cxcr2", "Il10", "Ccl1", "Ar", "Gdf15" ))
FeaturePlot(TRAMPC1Sub, features = c("Mrc1", "Ar", "Cd163", "Gdf15"), split.by = "orig.ident")
VlnPlot(TRAMPC1Sub, features = c("Mrc1", "Ar", "Cd163", "Gdf15"), split.by = "orig.ident", pt.size = 0)
VlnPlot(TRAMPC1Sub, features = c("Cd163", "Gdf15"), split.by = "orig.ident", idents = c("1","3","12","18"), pt.size = 0)

VlnPlot(TRAMPC1Sub, features = c("Cd163", "Gdf15"), pt.size = 0)
FeaturePlot(TRAMPC1Sub, features = c("Cd163", "Gdf15"), pt.size = 0)
FeaturePlot(TRAMPC1Sub, features = c("Cd68", "Gdf15"), pt.size = 0)
FeaturePlot(TRAMPC1Sub, features = c("Fscn1", "Gdf15"), pt.size = 0)
VlnPlot(TRAMPC1Sub, features = c("Cd163", "Gdf15"), pt.size = 0)
VlnPlot(TRAMPC1Sub, features = c("Cd68", "Gdf15"), pt.size = 0)
VlnPlot(TRAMPC1Sub, features = c("Fscn1", "Gdf15"), pt.size = 0)



#Paper Figures Create the plot

myPalette <- c(
  "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", 
  "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
  "#CAB2D6", "#6A3D9A"
)


# 2. Create your two plots as usual
plot1 <- DotPlot(
  object  = TRAMPC1Sub,
  features= "Xcr1",
  pt.size = 0
) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 3))
print(plot1)

plot2 <- DotPlot(
  object  = TRAMPC1Sub,
  features= "Ncr1",
  idents  = c("0", "1", "3", "5", "13", "8", "9", "10", "7"),
  pt.size = 0
) +
  coord_flip() +
  scale_fill_manual(values = myPalette) +
  scale_y_continuous(limits = c(0, 3))


# 3. Print or save
print(plot2)

#DOTPLOT


library(Seurat)
library(ggplot2)

# Marker genes
all_markers <- c(
  "Zbtb46", "Flt3", "Xcr1",       # Dendritic
  "Ncr1", "Klrc1", "Eomes",       # NK
  "Foxp3", "Il2ra", "Ctla4",      # Treg
  "Cd8a", "Cd8b1", "Runx3",       # CD8
  "Cd4", "Sell",                  # CD4
  "Cd19", "Ms4a1", "Cd79a",       # B cells
  "S100a9", "S100a8", "Ly6g",     # Neutrophils
  "Mrc1", "Cd163", "Cx3cr1",      # M2 macrophages
  "Gdf15"                         # GDF15
)

# Subset to selected clusters
selected_clusters <- c("13", "9", "10", "8", "16", "3", "7", "2", "1")
TRAMPC1Sub_subset <- subset(TRAMPC1Sub, idents = selected_clusters)

# Set identity levels in desired order (bottom to top plotting order)
Idents(TRAMPC1Sub_subset) <- factor(
  Idents(TRAMPC1Sub_subset),
  levels = rev(c("13", "9", "10", "8", "16", "7", "3", "2", "1"))
)

# Create dot plot
dot <- DotPlot(TRAMPC1Sub_subset, features = all_markers, dot.scale = 5) +
  scale_color_gradientn(
    colors = c("#FFFFFF", "#FF9999", "#CC0000"),
    limits = c(0, 2.5),
    oob = scales::squish,
    name = "Average Expression"
  ) +
  scale_size(
    range = c(1, 8),
    limits = c(0, 100),
    name = "Percent Expressed"
  ) +
  RotatedAxis() +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Show plot
dot




#Saving the data
setwd("/Users/alekhyamazumdar/Documents/scRNA/TRAMPC1/")
saveRDS(TRAMPC1Sub, "Seurat_TRAMPC1.rds")
save.image(file = "TRAMPC1_Dimplot.RData")

readRDS("/Users/alekhyamazumdar/Documents/scRNA/TRAMPC1/Seurat_TRAMPC1.rds")
