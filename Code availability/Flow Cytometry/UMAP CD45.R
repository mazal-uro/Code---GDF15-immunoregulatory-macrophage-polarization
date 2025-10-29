# Check if the column exists in your dataset
colnames(UMAP_df_CD45_Positive)
colnames(TSNEData_CD45_Positive)

# Add CD11c expression from the filtered dataset to the UMAP dataframe
UMAP_df_CD45_Positive$CD45_Expression <- TSNEData_CD45_Positive$`FJComp-APC-Cy7-A`


# Define a threshold for high Ly6G expression (top 10%)
CD45_threshold <- quantile(UMAP_df_CD45_Positive$CD45_Expression)

# Highlight only high Ly6G-expressing cells
ggplot(UMAP_df_CD45_Positive, aes(x = X, y = Y)) +
  geom_point(color = "gray", alpha = 0.2, size = 0.1) +  # Background
  geom_point(data = subset(UMAP_df_CD45_Positive, CD45_Expression >= CD45_threshold),
             aes(color = CD45_Expression), alpha = 0.5, size = 0.5) +
  facet_grid(Organ ~ Condition) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),  # White panel background
    plot.background = element_rect(fill = "white", color = NA)    # White overall background
  ) +
  scale_color_viridis_c(option = "inferno", limits = c(0, 1e6)) +
  ggtitle("UMAP Highlighting High CD45 Expressing Cells") +
  labs(color = "CD45Expression (Top 10%)")

# Export plot directly using ggsave() with explicit background settings
ggsave(
  filename = "UMAP_CD45_Fixed.svg",
  plot = last_plot(),  # Saves the most recent plot
  width = 4.5,
  height = 6,
  units = "in",
  bg = "white"  # Explicitly set background to white
)
