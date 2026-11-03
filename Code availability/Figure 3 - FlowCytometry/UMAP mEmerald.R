# Check if the column exists in your dataset
colnames(UMAP_df_CD45_Positive)

# Add CD11c expression from the filtered dataset to the UMAP dataframe
UMAP_df_CD45_Positive$mEM_Expression <- TSNEData_CD45_Positive$`FJComp-mEmerald-A`


# Define a threshold for high CD206 expression (top 20%)
mEM_threshold <- quantile(UMAP_df_CD45_Positive$mEM_Expression, 0.97)

# Highlight only high CD206-expressing cells
ggplot(UMAP_df_CD45_Positive, aes(x = X, y = Y)) +
  geom_point(color = "gray", alpha = 0.2, size = 0.1) +  # Background
  geom_point(data = subset(UMAP_df_CD45_Positive, mEM_Expression > mEM_threshold),
             aes(color = mEM_Expression), alpha = 0.5, size = 0.5) +
  facet_grid(Organ ~ Condition) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank()
  ) +
  scale_color_viridis_c(option = "inferno", limits = c(2e5, 7e6)) +
  ggtitle("UMAP Highlighting High mEMerald Expressing Cells") +
  labs(color = "mEmeraldExpression (Top 3%)")


# Export plot directly using ggsave() with explicit background settings
ggsave(
  filename = "UMAP_mEM_Fixed.svg",
  plot = last_plot(),  # Saves the most recent plot
  width = 4.5,
  height = 6,
  units = "in",
  bg = "white"  # Explicitly set background to white
)
