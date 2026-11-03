# Check if the column exists in your dataset
colnames(UMAP_df_CD45_Positive)

# Add CD206 expression from the filtered dataset to the UMAP dataframe
UMAP_df_CD45_Positive$CD206_Expression <- TSNEData_CD45_Positive$`FJComp-BV650-A`


# Define a threshold for high CD206 expression (top 10%)
cd206_threshold <- quantile(UMAP_df_CD45_Positive$CD206_Expression, 0.90)

# Highlight only high CD206-expressing cells
ggplot(UMAP_df_CD45_Positive, aes(x = X, y = Y)) +
  geom_point(color = "gray", alpha = 0.2, size = 0.1) +  # Background
  geom_point(data = subset(UMAP_df_CD45_Positive, CD206_Expression > cd206_threshold),
             aes(color = CD206_Expression), alpha = 0.5, size = 0.5) +
  facet_grid(Organ ~ Condition) +
  theme_minimal() +
  theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  plot.background = element_blank()
  ) +
  scale_color_viridis_c(option = "inferno", limits = c(0, 3e5)) +
  ggtitle("UMAP Highlighting High CD206 Expressing Cells") +
  labs(color = "CD206 Expression (Top 10%)")

# Export plot directly using ggsave() with explicit background settings
ggsave(
  filename = "UMAP_CD206_Fixed.svg",
  plot = last_plot(),  # Saves the most recent plot
  width = 4.5,
  height = 6,
  units = "in",
  bg = "white"  # Explicitly set background to white
)
