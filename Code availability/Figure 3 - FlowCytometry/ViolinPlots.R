# Generate a violin plot for CD206 expression
ggplot(UMAP_df_CD45_Positive, aes(x = Condition, y = `CD206_Expression`, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.7) +  # Violin plot without trimming
  geom_boxplot(width = 0.1, color = "black", alpha = 0.5) +  # Add boxplot inside violin plot
  facet_wrap(~Organ) +  # Separate plots for each organ
  theme_minimal() +
  scale_fill_manual(values = c("Tumor" = "#E63946", "Control" = "#457B9D")) + # Custom colors for conditions
  scale_y_continuous(limits = c(-1e5, 2e6), breaks = seq(0000, 5e5, by = 5e5)) +  # Y-axis from 0 to 2e6 with breaks every 500,000
  ggtitle("Violin Plot of CD206 Expression by Organ and Condition") +
  ylab("CD206 Expression (FJComp-BV650-A)") +
  xlab("Condition") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    legend.position = "none"
  )
