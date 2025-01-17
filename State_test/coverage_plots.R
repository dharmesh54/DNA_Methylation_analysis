library(dplyr)
library(ggplot2)
library(tidyr)
library(pheatmap)
# List of variable names
variable_names <- c("cfdna_m1_r1", "cfdna_m2_r1", "cfdna_m3_r1",
                    "cfdna_m1_r2", "cfdna_m2_r2", "cfdna_m3_r2",
                    "islet_m1", "islet_m2", "islet_m3")

# Combine all datasets into a single data frame
combined_data <- bind_rows(
  lapply(variable_names, function(var_name) {
    data <- get(var_name)
    data$tissue <- var_name  # Add tissue name as a column
    data
  })
)

# Summarize data (optional step if needed)
summary_stats <- combined_data %>%
  group_by(tissue) %>%
  summarize(
    mean_total_count = mean(N, na.rm = TRUE),
    median_total_count = median(N, na.rm = TRUE),
    mean_methylation = mean(X, na.rm = TRUE),
    median_methylation = median(X, na.rm = TRUE)
  )

print(summary_stats)

# Plot total read count distribution
ggplot(combined_data, aes(x = N, fill = tissue)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  theme_minimal() +
  labs(title = "Total Read Count Distribution",
       x = "Total Read Count",
       y = "Frequency") +
  scale_fill_brewer(palette = "Set3")

# Plot methylation level distribution
ggplot(combined_data, aes(x = X, fill = tissue)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  theme_minimal() +
  labs(title = "Methylation Level Distribution",
       x = "Methylation Level",
       y = "Frequency") +
  scale_fill_brewer(palette = "Set3")

# Boxplot of total read counts by tissue
ggplot(combined_data, aes(x = tissue, y = N, fill = tissue)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Total Read Counts by Tissue",
       x = "Tissue",
       y = "Total Read Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Boxplot of methylation levels by tissue
ggplot(combined_data, aes(x = tissue, y = X, fill = tissue)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Methylation Levels by Tissue",
       x = "Tissue",
       y = "Methylation Level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Combine all datasets into a single data frame
combined_data2 <- bind_rows(
  lapply(variable_names, function(var_name) {
    data <- get(var_name)
    data$tissue <- var_name  # Add tissue name as a column
    data
  })
)

# Reshape data: `pos` as rows, `variable_name` as columns, and methylation levels as values
reshaped_data <- combined_data2 %>%
  select(pos, X, tissue) %>%
  pivot_wider(names_from = tissue, values_from = X)


reshaped_data <- na.omit(reshaped_data)

# Extract methylation matrix for clustering
methylation_matrix <- as.matrix(reshaped_data[,-1])  # Remove `pos` column

#Compute distance matrix and hierarchical clustering
distance_matrix <- dist(t(methylation_matrix), method = "euclidean")  # Transpose for tissue clustering
hclust_result <- hclust(distance_matrix, method = "ward.D2")

# Plot dendrogram
plot(hclust_result, main = "Hierarchical Clustering of Tissues Based on Methylation",
     xlab = "Tissues", sub = "", cex = 0.8)


