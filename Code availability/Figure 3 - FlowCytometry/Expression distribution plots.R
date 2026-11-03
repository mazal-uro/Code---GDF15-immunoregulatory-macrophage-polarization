colnames(AllFCS) 
str(AllFCS$`FJComp-APC-Cy7-A`)
# Visualize CD45 expression distribution
ggplot(AllFCS, aes(x = `FJComp-APC-Cy7-A`)) +  # Use backticks for special characters
  geom_density(fill = "blue", alpha = 0.4) +
  ggtitle("CD45 Expression Density Plot") +
  xlab("FJComp-APC-Cy7-A") +
  xlim(-1e5, 2e6) +  # Limit x-axis from 0 to 2,000,000
  #xlim(50000, 2e6) +  # Limit x-axis from 0 to 2,000,000
  theme_minimal()

# Visualize mEM expression distribution
ggplot(AllFCS, aes(x = `FJComp-mEmerald-A`)) +  # Use backticks for special characters
  geom_density(fill = "green", alpha = 0.4) +
  ggtitle("mEM Expression Density Plot") +
  xlab("FJComp-mEmerald-A") +
  xlim(-1e5, 1e7) +  # Limit x-axis from 0 to 2,000,000
  #xlim(50000, 2e6) +  # Limit x-axis from 0 to 2,000,000
  theme_minimal()

# Visualize CD206 expression distribution
ggplot(AllFCS, aes(x = `FJComp-BV650-A`)) +  # Use backticks for special characters
  geom_density(fill = "red", alpha = 0.4) +
  ggtitle("CD206 Expression Density Plot") +
  xlab("FJComp-BV650-A") +
  xlim(-1e5, 4e5) +  # Limit x-axis from 0 to 2,000,000
  #xlim(50000, 2e6) +  # Limit x-axis from 0 to 2,000,000
  theme_minimal()

# Visualize CD11b expression distribution
ggplot(AllFCS, aes(x = `FJComp-BV510-A`)) +  # Use backticks for special characters
  geom_density(fill = "violet", alpha = 0.4) +
  ggtitle("CD11b Expression Density Plot") +
  xlab("FJComp-BV510-A") +
  xlim(1e4, 1e6) +  # Limit x-axis from 0 to 2,000,000
  #xlim(50000, 2e6) +  # Limit x-axis from 0 to 2,000,000
  theme_minimal()

# Visualize CD11c expression distribution
ggplot(AllFCS, aes(x = `FJComp-PE-Cy7-A`)) +  # Use backticks for special characters
  geom_density(fill = "orange", alpha = 0.4) +
  ggtitle("CD11c Expression Density Plot") +
  xlab("FJComp-PE-Cy7-A") +
  xlim(-1e5, 1e7) +  # Limit x-axis from 0 to 2,000,000
  #xlim(50000, 2e6) +  # Limit x-axis from 0 to 2,000,000
  theme_minimal()
