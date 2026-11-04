# -----------------------------------------
# STEP 0: Load libraries
# -----------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(R.utils)
library(tidyverse)
library(data.table)
library(cowplot)
library(SCpubr)
library(SingleCellExperiment)
library(scDblFinder)

# -----------------------------------------
# STEP 1: Load and format TRAMPC1 data
# -----------------------------------------
TRAMPC1 <- fread("/Users/alekhyamazumdar/Documents/scRNA/TRAMPC1/GSE272999_matrTRAMPC1.txt.gz")
GENE <- TRAMPC1[["names"]]
TRAMPC1 <- as.data.frame(TRAMPC1)
rownames(TRAMPC1) <- GENE
TRAMPC1[["V1"]] <- NULL
TRAMPC1[["names"]] <- NULL
TRAMPC1 <- as.matrix(TRAMPC1)

# -----------------------------------------
# STEP 2: Create Seurat object
# -----------------------------------------
TRAMPC1 <- CreateSeuratObject(counts = TRAMPC1, project = "TRAMPC1_Data", min.cells = 3, min.features = 200)

# -----------------------------------------
# STEP 3: QC and visualization
# -----------------------------------------
TRAMPC1[["percent.mt"]] <- PercentageFeatureSet(TRAMPC1, pattern = "^mt-")
VlnPlot(TRAMPC1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
FeatureScatter(TRAMPC1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# -----------------------------------------
# STEP 4: Normalize and find variable features
# -----------------------------------------
TRAMPC1 <- NormalizeData(TRAMPC1, normalization.method = "LogNormalize", scale.factor = 10000)
TRAMPC1 <- FindVariableFeatures(TRAMPC1, selection.method = "vst", nfeatures = 4000)

# -----------------------------------------
# STEP 5: Scale and PCA
# -----------------------------------------
all.genes <- rownames(TRAMPC1)
TRAMPC1 <- ScaleData(TRAMPC1, features = all.genes)
TRAMPC1 <- RunPCA(TRAMPC1, features = VariableFeatures(object = TRAMPC1))

# -----------------------------------------
# STEP 6: Run scDblFinder to identify doublets
# -----------------------------------------
sce <- as.SingleCellExperiment(TRAMPC1)
sce <- scDblFinder(sce)
TRAMPC1$doublet_class <- sce$scDblFinder.class
TRAMPC1$doublet_score <- sce$scDblFinder.score

# Subset singlets
TRAMPC1 <- subset(TRAMPC1, subset = doublet_class == "singlet")

# -----------------------------------------
# STEP 7: Dimensionality reduction and clustering
# -----------------------------------------
TRAMPC1 <- FindNeighbors(TRAMPC1, dims = 1:20)
TRAMPC1 <- FindClusters(TRAMPC1, resolution = 0.3)
TRAMPC1 <- RunUMAP(TRAMPC1, dims = 1:20)
DimPlot(TRAMPC1, reduction = "umap", group.by = "seurat_clusters")