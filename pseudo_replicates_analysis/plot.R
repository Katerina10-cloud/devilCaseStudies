
library(tidyverse)

final_results = readRDS("multi_dataset_results.rds")

# =============================================================================
# VISUALIZATION
# =============================================================================

# Calculate method ordering based on mean DEGs across all cell types
method_order <- final_results %>%
  group_by(method) %>%
  summarise(mean_deg = mean(n_deg, na.rm = TRUE)) %>%
  arrange(mean_deg) %>%
  pull(method)

# Create comprehensive plot
p1 <- final_results %>% 
  mutate(
    method = factor(method, levels = method_order),
    expression_group = factor(expression_group, levels = c("low", "medium", "high")),
    seq = ifelse(grepl("Pb", method), "Pseudo-bulk", "Single-cell")
  ) %>% 
  ggplot(aes(x = method, y = n_deg, fill = expression_group)) +
  geom_boxplot() +
  facet_grid( ~ experiment_type, scales = "free") +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Differential Expression Analysis Results by Cell Type",
    subtitle = "Comparison of Bio-replicates vs Pseudo-replicates",
    x = "Method",
    y = "Number of DEGs",
    fill = "Expression Group"
  ) +
  theme(
    strip.background = element_rect(fill = "lightgray"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p1)

# Alternative plot - focus on differences between experiment types
p2 <- final_results %>% 
  select(cell_type, experiment_type, method, expression_group, n_deg) %>%
  pivot_wider(names_from = experiment_type, values_from = n_deg) %>%
  mutate(
    deg_difference = pseudo_replicates - bio_replicates,
    method = factor(method, levels = method_order)
  ) %>%
  ggplot(aes(x = method, y = deg_difference, fill = expression_group)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~cell_type, scales = "free") +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Difference in DEG Detection: Pseudo-replicates - Bio-replicates",
    subtitle = "Positive values indicate more DEGs detected with pseudo-replicates",
    x = "Method",
    y = "Difference in Number of DEGs",
    fill = "Expression Group"
  ) +
  theme(
    strip.background = element_rect(fill = "lightgray"),
    legend.position = "bottom"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)

print(p2)

# Summary statistics
print("=== SUMMARY STATISTICS ===")
summary_stats <- final_results %>%
  group_by(cell_type, experiment_type, method, expression_group) %>%
  summarise(
    n_deg = first(n_deg),
    prop_deg = first(prop_deg),
    .groups = "drop"
  )

print(summary_stats)

# Summary by cell type
print("=== SUMMARY BY CELL TYPE ===")
celltype_summary <- final_results %>%
  group_by(cell_type, experiment_type) %>%
  summarise(
    total_deg = sum(n_deg),
    mean_prop_deg = mean(prop_deg),
    .groups = "drop"
  )

print(celltype_summary)

# Save results (optional)
# write_csv(final_results, "de_analysis_celltype_results.csv")
# write_csv(summary_stats, "de_analysis_celltype_summary.csv")
