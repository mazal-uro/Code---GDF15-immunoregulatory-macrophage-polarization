# Load required libraries
library(flowCore)
library(ggplot2)
library(dplyr)
library(Rtsne)

# Set working directory
setwd('/Users/alekhyamazumdar/Documents/Cytek Aurora/AM006 ReAnalysis with Beads/Unmixed/Clean Data/18thFEB Analysis')

# Define FCS files for tumor-bearing and control mice
fcs_files <- list(
  "Spleen_Tumor" = "Tumor Spleen_1.fcs",
  "Lungs_Tumor" = "Tumor Lungs_1.fcs",
  "Bone_Tumor" = "Tum BM_1.fcs",
  "LymphNodesCL_Tumor" = "Tum CLN_1.fcs",
  "LymphNodesDL_Tumor" = "Tumor DLN_1.fcs",
  "Liver_Tumor" = "Tum Liver_1.fcs",
  "Tumor_PC3" = "Tumor PC3_1.fcs",
  "Spleen_Control" = "Ctrl Spleen_1.fcs",
  "Lungs_Control" = "Ctrl Lungs_1.fcs",
  "Bone_Control" = "Ctrl BM_1.fcs",
  "LymphNodesCL_Control" = "Ctrl CLN_1.fcs",
  "LymphNodesDL_Control" = "Ctrl DLN_1.fcs",
  "Liver_Control" = "Ctrl Liver_1.fcs"
)

# Initialize an empty dataframe
AllFCS <- data.frame()

# Loop through files, read data, and assign metadata
for (sample in names(fcs_files)) {
  file_path <- fcs_files[[sample]]
  
  # Read FCS file
  fcs_data <- read.FCS(file_path, transformation = FALSE, truncate_max_range = FALSE)
  
  # Convert to dataframe
  df <- as.data.frame(exprs(fcs_data))
  
  # Extract Organ and Condition from sample name
  organ <- unlist(strsplit(sample, "_"))[1]  # Extract organ name
  condition <- ifelse(grepl("Tumor", sample), "Tumor", "Control")  # Label condition
  
  # Add metadata columns
  df$Organ <- organ
  df$Condition <- condition
  
  # Append to combined dataset
  AllFCS <- bind_rows(AllFCS, df)
}

# Print structure to confirm metadata
str(AllFCS)

# Print structure to confirm metadata is correct
table(AllFCS$Organ, AllFCS$Condition)

# Define fluorophores to keep
columns_to_keep <- c(
  "FJComp-APC-Cy7-A",
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

# Keep selected columns
TSNEData <- AllFCS[, intersect(columns_to_keep, colnames(AllFCS)), drop = FALSE]

# Retain metadata (Organ, Condition)
TSNEData$Organ <- AllFCS$Organ
TSNEData$Condition <- AllFCS$Condition

# Define the number of cells to keep (e.g., 10,000 instead of 40,000)
num_cells <- 10000  # Adjust as needed

# Set seed for reproducibility
set.seed(42)

# Randomly sample rows while keeping proportions
TSNE_sampled <- TSNEData %>%
  group_by(Organ, Condition) %>%  # Ensure balanced sampling from each group
  sample_n(size = num_cells / n_distinct(Organ, Condition), replace = FALSE) %>%
  ungroup()



# Run t-SNE on the sampled dataset
TSNE_result <- Rtsne(TSNE_sampled[, -c(ncol(TSNE_sampled)-1, ncol(TSNE_sampled))],
                     perplexity = 30, theta = 0.5, verbose = TRUE, check_duplicates = FALSE)

# Convert results to dataframe
TSNE_df <- data.frame(
  X = TSNE_result$Y[,1],
  Y = TSNE_result$Y[,2],
  Organ = TSNE_sampled$Organ,
  Condition = TSNE_sampled$Condition
)

table(TSNE_df$Organ, TSNE_df$Condition)

ggplot(TSNE_df, aes(x = X, y = Y, color = Condition)) +
  geom_point(alpha = 0.6, size = 0.7) +
  facet_grid(Organ ~ Condition) +  # Organ as rows, Condition as columns
  theme_minimal() +
  scale_color_manual(values = c("Tumor" = "#E63946", "Control" = "#457B9D")) +
  ggtitle("t-SNE: Organ-wise Tumor vs Control Comparison") +
  labs(color = "Condition")


colnames(AllFCS) 
str(AllFCS$`FJComp-APC-Cy7-A`)
# Visualize CD45 expression distribution
ggplot(AllFCS, aes(x = `FJComp-APC-Cy7-A`)) +  # Use backticks for special characters
  geom_density(fill = "blue", alpha = 0.4) +
  ggtitle("CD45 Expression Density Plot") +
  xlab("FJComp-APC-Cy7-A") +
  xlim(-1e6, 2e6) +  # Limit x-axis from 0 to 2,000,000
  #xlim(50000, 2e6) +  # Limit x-axis from 0 to 2,000,000
  theme_minimal()
