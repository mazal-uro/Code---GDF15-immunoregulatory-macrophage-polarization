# Check if the column exists in your dataset
colnames(UMAP_df_CD45_Positive)
colnames(TSNEData_CD45_Positive)

# Add CD11c expression from the filtered dataset to the UMAP dataframe
UMAP_df_CD45_Positive$CD103_Expression <- TSNEData_CD45_Positive$`FJComp-BV605-A`


# Define a threshold for high Ly6G expression (top 10%)
CD103_threshold <- quantile(UMAP_df_CD45_Positive$CD103_Expression, 0.90)

# Highlight only high Ly6G-expressing cells
ggplot(UMAP_df_CD45_Positive, aes(x = X, y = Y)) +
  geom_point(color = "gray", alpha = 0.2, size = 0.1) +  # Background
  geom_point(data = subset(UMAP_df_CD45_Positive, CD103_Expression > CD103_threshold),
             aes(color = CD103_Expression), alpha = 0.5, size = 0.5) +
  facet_grid(Organ ~ Condition) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank()
  ) +
  scale_color_viridis_c(option = "inferno", limits = c(0, 1e5)) +
  ggtitle("UMAP Highlighting High CD103 Expressing Cells") +
  labs(color = "CD103Expression (Top 10%)")


# Export plot directly using ggsave() with explicit background settings
ggsave(
  filename = "UMAP_CD103_Fixed.svg",
  plot = last_plot(),  # Saves the most recent plot
  width = 4.5,
  height = 6,
  units = "in",
  bg = "white"  # Explicitly set background to white
)
