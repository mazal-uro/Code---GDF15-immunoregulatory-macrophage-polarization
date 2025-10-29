# Filter cells with CD45 expression greater than 50,000
CD45_Positive_Cells <- AllFCS %>% 
  filter(`FJComp-APC-Cy7-A` > 50000)

# Downsample to 10,000 cells per organ and condition
set.seed(42)  # For reproducibility

# Group by organ and condition, then sample 10,000 cells per group
CD45_Positive_Downsampled <- CD45_Positive_Cells %>%
  group_by(Organ, Condition) %>%
  sample_n(size = min(10000, n()), replace = FALSE) %>%  # Avoid errors if there are fewer than 10,000 cells
  ungroup()

# Verify the number of cells sampled from each group
table(CD45_Positive_Downsampled$Organ, CD45_Positive_Downsampled$Condition)

# Define fluorescence channels for t-SNE
columns_to_keep <- c(
  "FJComp-APC-Cy7-A",  # CD45
  "FJComp-Alexa Fluor 700-A",
  "FJComp-BV421-A",
  "FJComp-BV510-A",
  "FJComp-BV605-A",
  "FJComp-BV650-A",
  "FJComp-PE-CF594-A",
  "FJComp-PerCP-Cy5.5-A",
  "FJComp-PE-Cy7-A",
  "FJComp-mEmerald-A"
)

# Filter and keep only the necessary columns
TSNEData_CD45_Positive <- CD45_Positive_Downsampled[, intersect(columns_to_keep, colnames(CD45_Positive_Downsampled)), drop = FALSE]

# Retain metadata for plotting
TSNEData_CD45_Positive$Organ <- CD45_Positive_Downsampled$Organ
TSNEData_CD45_Positive$Condition <- CD45_Positive_Downsampled$Condition


# Set seed for reproducibility
set.seed(42)

# Run t-SNE on the filtered dataset
TSNE_result <- Rtsne(TSNEData_CD45_Positive[, -c(ncol(TSNEData_CD45_Positive)-1, ncol(TSNEData_CD45_Positive))], 
                     perplexity = min(50, nrow(TSNEData_CD45_Positive)/3), 
                     theta = 0.5, verbose = TRUE, check_duplicates = FALSE)

# Create a data frame for visualization
TSNE_df_CD45_Positive <- data.frame(
  X = TSNE_result$Y[,1],
  Y = TSNE_result$Y[,2],
  Organ = TSNEData_CD45_Positive$Organ,
  Condition = TSNEData_CD45_Positive$Condition
)

ggplot(TSNE_df_CD45_Positive, aes(x = X, y = Y, color = Condition)) +
  geom_point(alpha = 0.6, size = 0.7) +
  facet_grid(Organ ~ Condition) +  # Organ as rows, Condition as columns
  theme_minimal() +
  scale_color_manual(values = c("Tumor" = "#E63946", "Control" = "#457B9D")) +
  ggtitle("t-SNE: CD45+ Cells (>50,000) with 10,000 Cells per Organ and Condition") +
  labs(color = "Condition")


library(umap)
# Prepare data for UMAP (exclude metadata columns)
UMAP_data <- TSNEData_CD45_Positive[, -c(ncol(TSNEData_CD45_Positive)-1, ncol(TSNEData_CD45_Positive))]

# Run UMAP
umap_result <- umap(UMAP_data)

# Create a dataframe with UMAP results and metadata
UMAP_df_CD45_Positive <- data.frame(
  X = umap_result$layout[, 1],
  Y = umap_result$layout[, 2],
  Organ = TSNEData_CD45_Positive$Organ,
  Condition = TSNEData_CD45_Positive$Condition
)

ggplot(UMAP_df_CD45_Positive, aes(x = X, y = Y, color = Condition)) +
  geom_point(alpha = 0.2, size = 0.1) +
  facet_grid(Organ ~ Condition) +  # Organ in rows, Condition in columns
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.background = element_rect(fill = "white", color = NA),  # White panel background
    plot.background = element_rect(fill = "white", color = NA)    # White overall background
  ) +
  scale_color_manual(values = c("Tumor" = "#E63946", "Control" = "#457B9D")) +
  ggtitle("UMAP: CD45+ Cells (>50,000) with 10,000 Cells per Organ and Condition") +
  labs(color = "Condition")

#Saving Everything in R
save.image("UMAP_workspace.RData")


load("UMAP_workspace.RData")