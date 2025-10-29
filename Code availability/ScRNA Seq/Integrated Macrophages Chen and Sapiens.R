# -----------------------------------------
# Load necessary libraries
# -----------------------------------------

library(dplyr)
library(Seurat)
library(patchwork)
library(R.utils)
library(tidyverse)
library(data.table)
library(cowplot)
library(SCpubr)
library(ggnewscale)
library(pheatmap)
library(ggrepel) 
library(clusterProfiler)
library(tibble)
library(SCpubr)
library(Azimuth)
library(SeuratData)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)  # Human gene annotations
library(enrichplot)    # Pretty plots for enrichment

# -----------------------------------------
# STEP 1: Load the Two Datasets
# -----------------------------------------
SapiensMacrophages <- readRDS("/Users/alekhyamazumdar/Documents/scRNA/Tabula Sapiens - Prostate/SapiensProstateMacrophages.rds")
ChenMacrophages <- readRDS("/Users/alekhyamazumdar/Documents/scRNA/Chen et al/Seurat_ChenEtAlmacrophages.rds")

# ----------------------------------------
# STEP 2. Create Seurat list and Cluster after SCTransform
# ----------------------------------------
seurat_list <- list(SapiensMacrophages, ChenMacrophages)

# ----------------------------------------
# 2. Calculate percent.mt and percent.ribo if needed
for (i in 1:length(seurat_list)) {
  seurat_list[[i]][["percent.mt"]] <- PercentageFeatureSet(seurat_list[[i]], pattern = "^MT-")
  seurat_list[[i]][["percent.ribo"]] <- PercentageFeatureSet(seurat_list[[i]], pattern = "^RPS|^RPL")
}

# ----------------------------------------
# 3. SCTransform with regression
for (i in 1:length(seurat_list)) {
  seurat_list[[i]] <- SCTransform(
    seurat_list[[i]],
    vars.to.regress = c("percent.mt", "percent.ribo"),
    verbose = FALSE
  )
}

# ----------------------------------------
# 4. Select integration features
features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)

# ----------------------------------------
# 5. Prepare for SCT integration
seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features)

# ----------------------------------------
# 6. Find anchors
anchors <- FindIntegrationAnchors(
  object.list = seurat_list,
  normalization.method = "SCT",
  anchor.features = features
)

# ----------------------------------------
# 7. Integrate data
integrated_obj <- IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT"
)

# ----------------------------------------
# 8. Dimensionality reduction
integrated_obj <- RunPCA(integrated_obj)
integrated_obj <- RunUMAP(integrated_obj, dims = 1:20)

# ----------------------------------------
# 9. Clustering and visualization
integrated_obj <- FindNeighbors(integrated_obj, dims = 1:20)
integrated_obj <- FindClusters(integrated_obj, resolution = 0.1)

# UMAP plots
DimPlot(integrated_obj, group.by = "orig.ident", pt.size = 1.5)
DimPlot(integrated_obj, group.by = "seurat_clusters", label = TRUE, pt.size = 1.5)

# -----------------------------------------
# STEP 3: UMAP Visualizations for publication
# -----------------------------------------
# -----------------------------------------
# UMAP of Integrated Clusters
# -----------------------------------------
#UMAP per clustering
p_integrated_umap <- SCpubr::do_DimPlot(
  sample = integrated_obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  pt.size = 1.5,
  label = TRUE,
  plot.title = "UMAP of Integrated Macrophage Subclusters",
  plot.subtitle = "Sapiens + Chen Integration"
)
print(p_integrated_umap)

#UMAP per Group 
p_integrated_umap <- SCpubr::do_DimPlot(
  sample = integrated_obj,
  reduction = "umap",
  group.by = "orig.ident",  # Grouping by dataset/source ID
  pt.size = 1.5,
  label = TRUE,
  plot.title = "UMAP of Integrated Macrophages by Sample",
  plot.subtitle = "Grouped by orig.ident (Sapiens vs. Chen)"
)

print(p_integrated_umap)
# -----------------------------------------
# GDF15 Feature Plot with Controlled Color Scale
# -----------------------------------------
p_integrated_gdf15 <- SCpubr::do_FeaturePlot(
  sample = integrated_obj,
  features = "GDF15",
  pt.size = 1.5,
  plot.title = "GDF15 Expression on Integrated UMAP",
  min.cutoff = 0,
  max.cutoff = 3.0
)
print(p_integrated_gdf15)
# -----------------------------------------
# Combine into One Cowplot Panel
# -----------------------------------------
cowplot::plot_grid(
  p_integrated_umap,
  p_integrated_gdf15,
  labels = c("A", "B"),
  ncol = 2,
  align = "hv"
)

# -----------------------------------------
# STEP 4: GDF15 Expression Plots (Split by Dataset)
# -----------------------------------------
# -----------------------------------------
# Set RNA assay for expression plots
# -----------------------------------------
DefaultAssay(integrated_obj) <- "RNA"
table(integrated_obj$orig.ident)
integrated_obj$orig.ident <- recode(integrated_obj$orig.ident,
                                    "SeuratProject" = "Sapiens",
                                    "Chen_Data" = "Chen")
# -----------------------------------------
# Subset the integrated object by dataset
# -----------------------------------------
sapiens_only <- subset(integrated_obj, subset = orig.ident == "Sapiens")
chen_only    <- subset(integrated_obj, subset = orig.ident == "Chen")

# -----------------------------------------
# Plot GDF15 expression for Sapiens
# -----------------------------------------
p_gdf15_sapiens <- SCpubr::do_FeaturePlot(
  sample = sapiens_only,
  features = "GDF15",
  pt.size = 1.5,
  plot.title = "GDF15 Expression (Sapiens)",
  min.cutoff = 0.0,
  max.cutoff = 3.0,
  assay = "RNA"  # explicitly use RNA values
)
print(p_gdf15_sapiens)

# -----------------------------------------
# Plot GDF15 expression for Chen
# -----------------------------------------
p_gdf15_chen <- SCpubr::do_FeaturePlot(
  sample = chen_only,
  features = "GDF15",
  pt.size = 1.5,
  plot.title = "GDF15 Expression (Chen)",
  min.cutoff = 0.0,
  max.cutoff = 3.0,
  assay = "RNA"
)
print(p_gdf15_chen)

# -----------------------------------------
# Combine plots side-by-side using cowplot
# -----------------------------------------
cowplot::plot_grid(
  p_gdf15_sapiens,
  p_gdf15_chen,
  labels = c("A", "B"),
  ncol = 2,
  align = "hv"
)

# -----------------------------------------
# Violin Plot 
# -----------------------------------------

# Set RNA assay
DefaultAssay(integrated_obj) <- "RNA"

# Base violin plot (split, no points)
base_vln <- VlnPlot(
  object = integrated_obj,
  features = "GDF15",
  group.by = "seurat_clusters",
  split.by = "orig.ident",
  split.plot = TRUE,
  pt.size = 0.5,
  cols = c("Chen" = "#C83658FF", "Sapiens" = "#00925EFF")
)

# Fetch data and prepare offset
vln_data <- FetchData(integrated_obj, vars = c("GDF15", "seurat_clusters", "orig.ident"))
vln_data$seurat_clusters <- as.factor(vln_data$seurat_clusters)

# Filter out missing orig.ident cells
vln_data <- vln_data %>% filter(!is.na(orig.ident))

# Flip the jitter offset: Chen LEFT, Sapiens RIGHT
vln_data$x_dodge <- as.numeric(vln_data$seurat_clusters) +
  ifelse(vln_data$orig.ident == "Chen", -0.2, 0.2)

# Build final plot
vln_custom <- base_vln +
  geom_point(
    data = vln_data,
    aes(
      x = x_dodge,
      y = GDF15,
      fill = orig.ident,
      color = orig.ident
    ),
    shape = 21,
    stroke = 1,
    size = 1,
    alpha = 0.7,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = c("Chen" = "#C83658FF", "Sapiens" = "#00925EFF")) +
  scale_color_manual(values = c("Chen" = "#C83658FF", "Sapiens" = "#00925EFF")) +
  
  # Horizontal line
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "grey40", linewidth = 1) +
  
  # Clean theme
  theme_minimal(base_size = 14) +
  labs(
    title = "GDF15 Expression by Cluster (RNA Assay)",
    x = "Cluster",
    y = "Log Normalized Expression"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(face = "bold", size = 12),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.major = element_line(color = "grey85", linetype = "dotted", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# Show the beautiful plot
print(vln_custom)

# -----------------------------------------
# DOT Plot (Clusters on X, Features on Y, Group Separators)
# -----------------------------------------
DefaultAssay(integrated_obj) <- "RNA"

marker_genes <- c(
  "CD68", "CD14", "ITGAM", "CSF1R", "AIF1", "FCGR1A", "TYROBP", "HLA-DRA", "APOE",
  "IL1B", "TNF", "IL6", "CD80", "CD86", "STAT1", "CXCL9", "CXCL10", "NOS2", "CCL2",
  "MRC1", "IL10", "TGFBI", "CD163", "ARG1", "FOLR2", "PPARG", "IRF4", "CCL18", "NRP1", "TGM2", "CD200R1", "IL1R2",
  "CD274", "PDCD1LG2"
)

integrated_obj$cluster_only <- paste0("Cluster ", integrated_obj$seurat_clusters)

dot_data <- DotPlot(
  integrated_obj,
  features = marker_genes,
  group.by = "cluster_only"
)$data

dot_data$features.plot <- factor(dot_data$features.plot, levels = rev(marker_genes))
dot_data$id <- factor(dot_data$id, levels = unique(dot_data$id))

num_clusters <- length(levels(dot_data$id))
num_features <- length(marker_genes)

# Separator positions after APOE (9th) and CCL2 (19th)
separator_positions <- c(9.5, 19.5)

# Plot
ggplot(dot_data, aes(x = id, y = features.plot)) +
  geom_point(aes(fill = avg.exp.scaled, size = pct.exp), shape = 21, stroke = 0.2, color = "black") +
  scale_fill_gradient2(
    low = "navy",
    mid = "white",
    high = "firebrick",
    midpoint = 0,
    limits = c(-2, 2),
    name = "Avg. Exp. (Scaled)"
  ) +
  scale_size(range = c(1, 6), name = "% Expressing") +
  
  # Vertical separators
  geom_hline(yintercept = separator_positions, color = "grey40", linetype = "dashed", linewidth = 0.75) +
  
  # Outer box
  annotate("segment", x = 0.5, xend = num_clusters + 0.5, y = 0.5, yend = 0.5, color = "black") +
  annotate("segment", x = 0.5, xend = num_clusters + 0.5, y = num_features + 0.5, yend = num_features + 0.5, color = "black") +
  annotate("segment", x = 0.5, xend = 0.5, y = 0.5, yend = num_features + 0.5, color = "black") +
  annotate("segment", x = num_clusters + 0.5, xend = num_clusters + 0.5, y = 0.5, yend = num_features + 0.5, color = "black") +
  
  theme_minimal(base_size = 14) +
  labs(
    title = "Macrophage Marker Expression by Cluster",
    x = "Cluster",
    y = "Genes"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right"
  )

# -----------------------------------------
# Volcano Plot and Pheatmap
# -----------------------------------------

# -----------------------------------------
# Step 1: Switch to RNA assay
# -----------------------------------------
DefaultAssay(integrated_obj) <- "RNA"

# -----------------------------------------
# Step 2: Calculate average expression per cluster
# -----------------------------------------
# This uses log-normalized RNA values

cluster_averages <- AverageExpression(
  integrated_obj,
  group.by = "seurat_clusters",
  assays = "RNA",
  slot = "data"  # log-normalized expression
)$RNA  # pull only RNA output

# Check the matrix
dim(cluster_averages)  # genes × clusters

# -----------------------------------------
# Step 3: Select top variable genes (optional)
# -----------------------------------------
# Calculate variance across clusters
gene_variance <- apply(cluster_averages, 1, var)

# Select top 100 most variable genes
top_genes <- names(sort(gene_variance, decreasing = TRUE))[1:100]

# Subset matrix
cluster_averages_top <- cluster_averages[top_genes, ]


# -----------------------------------------
# Step 2: Differential Expression: Cluster 1 vs all other clusters
# -----------------------------------------
DefaultAssay(integrated_obj) <- "RNA"
integrated_obj <- JoinLayers(integrated_obj)  # VERY important for Seurat v5

markers_cluster1 <- FindMarkers(
  integrated_obj,
  ident.1 = "1",    # cluster 1
  group.by = "seurat_clusters",
  assay = "RNA",
  slot = "data",    # log-normalized
  logfc.threshold = 0,   # no threshold, get all genes
  min.pct = 0.1     # only consider genes expressed in at least 10% of cells
)

# View top markers
head(markers_cluster1)
# Check distribution inside Cluster 1
table(integrated_obj$orig.ident[integrated_obj$seurat_clusters == "1"])

# -----------------------------------------
# Add significance and labeling
# -----------------------------------------
logfc_cutoff <- 0.5
padj_cutoff <- 0.05

# Add gene names and significance category
markers_cluster1 <- markers_cluster1 %>%
  mutate(
    gene = rownames(.),
    significance = case_when(
      p_val_adj < padj_cutoff & avg_log2FC > logfc_cutoff  ~ "Up",
      p_val_adj < padj_cutoff & avg_log2FC < -logfc_cutoff ~ "Down",
      TRUE                                                ~ "NotSig"
    )
  )

# Identify Top 50 Up and Down regulated genes to label
top_up <- markers_cluster1 %>%
  filter(significance == "Up") %>%
  arrange(p_val_adj) %>%
  slice_head(n = 50)

top_down <- markers_cluster1 %>%
  filter(significance == "Down") %>%
  arrange(p_val_adj) %>%
  slice_head(n = 50)

# Add labeling column
markers_cluster1 <- markers_cluster1 %>%
  mutate(
    label = ifelse(gene %in% c(top_up$gene, top_down$gene), gene, NA)
  )

# -----------------------------------------
# Volcano Plot with Labels and Highlighted GDF15
# -----------------------------------------

ggplot(markers_cluster1, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
  
  geom_point(alpha = 0.7, size = 1.5) +
  
  scale_color_manual(values = c(
    "Up" = "firebrick",
    "Down" = "navy",
    "NotSig" = "grey70"
  )) +
  
  # Label Top 30 Up/Down genes
  geom_text_repel(
    aes(label = label),
    size = 4,
    box.padding = 0.5,
    point.padding = 0.2,
    max.overlaps = 100,
    segment.color = "grey80"
  ) +
  
  theme_minimal(base_size = 14) +
  labs(
    title = "Cluster 1 vs All: Volcano Plot (Top 30 Up/Down + GDF15)",
    x = "Log2 Fold Change (Cluster 1 vs Others)",
    y = "-Log10 Adjusted P-Value",
    color = "Gene Status"
  ) +
  
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey50") +
  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  
  guides(color = guide_legend(override.aes = list(size = 4)))


# -----------------------------------------
# 1. Differential Expression: Cluster 1 vs All Clusters
# -----------------------------------------
DefaultAssay(integrated_obj) <- "RNA"
integrated_obj <- JoinLayers(integrated_obj)

markers_cluster1 <- FindMarkers(
  integrated_obj,
  ident.1 = "1",
  group.by = "seurat_clusters",
  assay = "RNA",
  slot = "data",
  logfc.threshold = 0,
  min.pct = 0.5
)
markers_cluster1$gene <- rownames(markers_cluster1)

# -----------------------------------------
# 2. Classify Upregulated and Downregulated Genes
# -----------------------------------------
logfc_cutoff <- 0.5
padj_cutoff <- 0.05

markers_cluster1 <- markers_cluster1 %>%
  mutate(
    regulation = case_when(
      avg_log2FC >= logfc_cutoff & p_val_adj < padj_cutoff ~ "Up",
      avg_log2FC <= -logfc_cutoff & p_val_adj < padj_cutoff ~ "Down",
      TRUE ~ "NotSig"
    )
  )

up_genes <- markers_cluster1 %>% filter(regulation == "Up") %>% pull(gene)
down_genes <- markers_cluster1 %>% filter(regulation == "Down") %>% pull(gene)

# -----------------------------------------
# 3. Convert Gene Symbols to Entrez IDs
# -----------------------------------------
up_genes_entrez <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
down_genes_entrez <- bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

# -----------------------------------------
# 4. GO Enrichment Analysis (Biological Process)
# -----------------------------------------
ego_up <- enrichGO(
  gene = up_genes_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

ego_down <- enrichGO(
  gene = down_genes_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# -----------------------------------------
# 5. Barplots for GO Enrichment Results
# -----------------------------------------
barplot(ego_up, showCategory = 15, title = "GO BP - Upregulated Genes")
barplot(ego_down, showCategory = 15, title = "GO BP - Downregulated Genes")

# -----------------------------------------
# 6. Dotplots with Custom Scaling
# -----------------------------------------
dotplot(ego_up, showCategory = 12, title = "GO BP - Upregulated Genes") +
  scale_color_gradient(low = "red", high = "blue", trans = "log10") +
  xlim(0, 0.1)

ego_up@result

dotplot(ego_down, showCategory = 12, title = "GO BP - Downregulated Genes") +
  scale_color_gradient(low = "red", high = "blue", trans = "log10") +
  xlim(0, 0.15)

ego_down@result


# -----------------------------------------
# STEP 1: HEATMAP of all clusters
# -----------------------------------------

# -----------------------------------------
# STEP 1: Select Top DE Genes (Up + Down)
# -----------------------------------------
top_genes <- c(top_up$gene, top_down$gene)

# -----------------------------------------
# STEP 2: Extract Log-Normalized Expression Matrix
# -----------------------------------------
# Ensure RNA assay and log-normalized data
DefaultAssay(integrated_obj) <- "RNA"

# Extract normalized data (log1p counts)
expression_matrix <- GetAssayData(integrated_obj, slot = "data")[top_genes, ]

# Optionally subset to cluster 1 or all clusters
# Example: Show all cells ordered by cluster
cell_order <- integrated_obj$seurat_clusters

# -----------------------------------------
# STEP 3: Create Annotation for Clusters
# -----------------------------------------
annotation_col <- data.frame(
  Cluster = factor(integrated_obj$seurat_clusters)
)
rownames(annotation_col) <- colnames(expression_matrix)

# -----------------------------------------
# STEP 4: Draw the Heatmap
# -----------------------------------------
pheatmap(
  expression_matrix,
  scale = "row",  # gene-wise scaling
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 6,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  main = "Top 50 Up and Down Regulated Genes in Cluster 1"
)

# -----------------------------------------
# Cluster 1 Chen vs All Sapiens Cells analysis
# -----------------------------------------
DefaultAssay(integrated_obj) <- "RNA"
integrated_obj <- JoinLayers(integrated_obj)

# -----------------------------------------
# Step 2: Subset Cells
# -----------------------------------------

# Subset Cluster 1 Chen cells
chen_cluster1_cells <- subset(integrated_obj, subset = (seurat_clusters == "1" & orig.ident == "Chen"))

# Subset all Sapiens cells
sapiens_all_cells <- subset(integrated_obj, subset = (orig.ident == "Sapiens"))

# -----------------------------------------
# Step 3: Merge Chen Cluster 1 and All Sapiens into new object
# -----------------------------------------
comparison_obj <- merge(chen_cluster1_cells, y = sapiens_all_cells)

# -----------------------------------------
# Step 4: Create a new label for comparison
# -----------------------------------------
comparison_obj$comparison_group <- ifelse(comparison_obj$orig.ident == "Chen", "Chen_Cluster1", "Sapiens_All")

# Set Idents
Idents(comparison_obj) <- comparison_obj$comparison_group

# -----------------------------------------
# Step 5: Differential Expression Analysis
# -----------------------------------------
comparison_obj <- JoinLayers(comparison_obj)
de_results <- FindMarkers(
  comparison_obj,
  ident.1 = "Chen_Cluster1",
  ident.2 = "Sapiens_All",
  assay = "RNA",
  slot = "data",
  logfc.threshold = 0,   # get all genes
  min.pct = 0.1          # only keep genes expressed in 10% of cells
)

# View top genes
head(de_results)


# -----------------------------------------
# Step 1: Select Top 100 DE Genes: Cluster 1 Chen vs All Sapiens Cells analysis
# -----------------------------------------

# Sort DE results by adjusted p-value (or log2FC if you prefer)
top_genes <- de_results %>%
  arrange(p_val_adj) %>%
  slice(1:100) %>%
  pull(gene)  # Pull gene names

# -----------------------------------------
# Step 2: Fetch expression data
# -----------------------------------------

# Make sure RNA assay is active
DefaultAssay(comparison_obj) <- "RNA"

# Fetch expression matrix
expression_matrix <- GetAssayData(comparison_obj, slot = "data")

# Subset to top genes
expression_matrix_top <- expression_matrix[top_genes, ]

# -----------------------------------------
# Step 3: Create Annotation Bar (Optional, for groups)
# -----------------------------------------
annotation_col <- data.frame(
  Group = comparison_obj$comparison_group
)
rownames(annotation_col) <- colnames(expression_matrix_top)

# -----------------------------------------
# Step 4: Plot Heatmap: Cluster 1 Chen vs All Sapiens Cells analysis
# -----------------------------------------


pheatmap(
  expression_matrix_top,
  scale = "row",  # normalize gene-wise
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,  # color bar by group (Chen_Cluster1 vs Sapiens_All)
  show_rownames = TRUE,    # show gene names (y-axis)
  show_colnames = FALSE,   # <--- hide x-axis labels
  fontsize_row = 6,        # Adjust font size
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  main = "Top 100 DE Genes: Chen Cluster 1 vs Sapiens"
)


# -----------------------------------------
# Step 1: Select Top 100 DE Genes
# -----------------------------------------
top100_genes <- de_results %>%
  arrange(p_val_adj) %>%
  slice(1:100) %>%
  pull(gene)

# -----------------------------------------
# Step 2: Add Highlight Column
# -----------------------------------------
de_results <- de_results %>%
  mutate(
    gene = rownames(.),
    label = ifelse(gene %in% top100_genes, gene, NA),  # label only Top 100
    highlight = ifelse(gene %in% top100_genes, "Top50", "Other")
  )

# -----------------------------------------
# Step 1: Define Upregulated / Downregulated / Non-significant
# -----------------------------------------

# Set thresholds
logfc_cutoff <- 0.25  # You can adjust
padj_cutoff <- 0.005   # 0.5% FDR

# Classify genes
de_results <- de_results %>%
  mutate(
    gene = rownames(.),
    significance = case_when(
      p_val_adj < padj_cutoff & avg_log2FC > logfc_cutoff ~ "Up",
      p_val_adj < padj_cutoff & avg_log2FC < -logfc_cutoff ~ "Down",
      TRUE ~ "NotSig"
    )
  )

# -----------------------------------------
# Step 2: Select Top 30 Up and Top 30 Down Genes for Labeling
# -----------------------------------------

# Top Upregulated
top30_up <- de_results %>%
  filter(significance == "Up") %>%
  arrange(p_val_adj) %>%
  slice(1:30) %>%
  pull(gene)

# Top Downregulated
top30_down <- de_results %>%
  filter(significance == "Down") %>%
  arrange(p_val_adj) %>%
  slice(1:30) %>%
  pull(gene)

# Merge top genes to label
genes_to_label <- c(top30_up, top30_down)

# Add label column
de_results <- de_results %>%
  mutate(
    label = ifelse(gene %in% genes_to_label, gene, NA)
  )

# -----------------------------------------
# Step 3: Volcano Plot
# -----------------------------------------

ggplot(de_results, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
  
  geom_point(alpha = 0.7, size = 1.5) +
  
  scale_color_manual(values = c(
    "Up" = "firebrick",
    "Down" = "navy",
    "NotSig" = "grey70"
  )) +
  
  # Label Top 30 Up/Down genes
  geom_text_repel(
    aes(label = label),
    size = 4,
    box.padding = 0.5,
    point.padding = 0.2,
    max.overlaps = 100,
    segment.color = "grey80"
  ) +
  
  # Special label for GDF15
  geom_text_repel(
    data = subset(de_results, gene == "GDF15"),
    aes(label = gene),
    color = "black",        # Black font
    size = 5,               # Match other labels
    box.padding = 0.5,
    point.padding = 0.2,
    segment.color = "grey60",  # Dash line color same as others
    segment.size = 0.5,
    max.overlaps = Inf
  )+
  
  theme_minimal(base_size = 14) +
  labs(
    title = "Chen Cluster 1 vs All Sapiens: Volcano Plot (Top 30 Up & Down Genes + GDF15)",
    x = "Log2 Fold Change (Chen Cluster 1 vs Sapiens All)",
    y = "-Log10 Adjusted P-Value",
    color = "Gene Status"
  ) +
  
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey50") +
  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  
  guides(color = guide_legend(override.aes = list(size = 4)))

# -----------------------------------------
# Step 3: GO Biological plots
# -----------------------------------------

# -----------------------------------------
# Step 1: Prepare Up and Downregulated Genes
# -----------------------------------------

# Filter upregulated genes
up_genes <- de_results %>%
  filter(significance == "Up") %>%
  pull(gene)

# Filter downregulated genes
down_genes <- de_results %>%
  filter(significance == "Down") %>%
  pull(gene)

# -----------------------------------------
# Step 2: Convert Gene Symbols to Entrez IDs
# -----------------------------------------
# clusterProfiler works best with Entrez IDs

up_genes_entrez <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
down_genes_entrez <- bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

# -----------------------------------------
# Step 3: GO Biological Process Enrichment
# -----------------------------------------

# Upregulated genes GO
ego_up <- enrichGO(
  gene = up_genes_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "BP",    # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# Downregulated genes GO
ego_down <- enrichGO(
  gene = down_genes_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# -----------------------------------------
# Step 4: Plot GO Enrichment Results
# -----------------------------------------

# Top GO terms - Upregulated
barplot(ego_up, showCategory = 15, title = "GO BP Enrichment - Upregulated Genes")

# Top GO terms - Downregulated
barplot(ego_down, showCategory = 15, title = "GO BP Enrichment - Downregulated Genes")

# -----------------------------------------
# Dotplot for GO Enrichment - Upregulated Genes
# -----------------------------------------
dotplot(ego_up, showCategory = 12, title = "GO BP Enrichment - Upregulated Genes") +
  
  xlim(0, 0.1) +  # Set x-axis limits for GeneRatio
  
  scale_color_continuous(
    low = "red",
    high = "blue",
    trans = "log10",    # log scale for p.adjust values
  ) +
  
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold")
  )

# -----------------------------------------
# Dotplot for GO Enrichment - Downregulated Genes
# -----------------------------------------
dotplot(ego_down, showCategory = 12, title = "GO BP Enrichment - Downregulated Genes") +
  
  xlim(0, 0.1) +  # Set x-axis limits for GeneRatio
  
  scale_color_continuous(
    low = "red",
    high = "blue",
    trans = "log10",    # log scale for p.adjust values
  ) +
  
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold")
  )


# Assign the tibble into ego_up_df
ego_up_df <- tibble(
  Description = as.character(ego_up@result$Description),
  GeneRatio = sapply(
    strsplit(as.character(ego_up@result$GeneRatio), "/"),
    function(x) as.numeric(x[1]) / as.numeric(x[2])
  ),
  p.adjust = as.numeric(ego_up@result$p.adjust),
  Count = as.numeric(ego_up@result$Count)
)
# Select Top 12 by direct indexing without dplyr
ego_up_df_top12 <- ego_up_df[order(ego_up_df$p.adjust), ][1:12, ]

ggplot(ego_up_df_top12, aes(x = GeneRatio, y = reorder(Description, GeneRatio))) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "navy", high = "firebrick", trans = "log10") +
  xlim(0, 0.1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "GO BP Enrichment - Top 12 Upregulated Genes",
    x = "Gene Ratio",
    y = NULL,
    color = "Adjusted P-value",
    size = "Gene Count"
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    axis.text.x = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold")
)



setwd('/Users/alekhyamazumdar/Documents/scRNA/Tabula Sapiens - Prostate')
save.image(file = "IntegratedMacrophages.RData")

