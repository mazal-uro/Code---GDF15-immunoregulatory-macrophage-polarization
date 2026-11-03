# Check if the column exists in your dataset
colnames(UMAP_df_CD45_Positive)
colnames(TSNEData_CD45_Positive)

# Add CD11c expression from the filtered dataset to the UMAP dataframe
UMAP_df_CD45_Positive$CD11b_Expression <- TSNEData_CD45_Positive$`FJComp-BV510-A`


# Define a threshold for high Ly6G expression (top 20%)
CD11b_threshold <- quantile(UMAP_df_CD45_Positive$CD11b_Expression, 0.70)

# Highlight only high Ly6G-expressing cells
ggplot(UMAP_df_CD45_Positive, aes(x = X, y = Y)) +
  geom_point(color = "gray", alpha = 0.2, size = 0.1) +  # Background points
  geom_point(data = subset(UMAP_df_CD45_Positive, CD11b_Expression > CD11b_threshold),
             aes(color = CD11b_Expression), alpha = 0.5, size = 0.5) +
  facet_grid(Organ ~ Condition) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.background = element_rect(fill = "white", color = NA),  # White panel background
    plot.background = element_rect(fill = "white", color = NA)    # White overall background
  ) +
  scale_color_viridis_c(option = "inferno", limits = c(0, 5e5)) +
  ggtitle("UMAP Highlighting High CD11b Expressing Cells") +
  labs(color = "CD11b Expression (Top 30%)")

# Export plot directly using ggsave() with explicit background settings
ggsave(
  filename = "UMAP_CD11b_Fixed.svg",
  plot = last_plot(),  # Saves the most recent plot
  width = 4.5,
  height = 6,
  units = "in",
  bg = "white"  # Explicitly set background to white
)
