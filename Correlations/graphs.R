# =========================================
# VISUALIZACIÓN DE RESULTADOS DE CORRELACIÓN
# =========================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(viridis)
library(patchwork)

# -------------------------
# 1. CARGAR DATOS
# -------------------------

# Cargar los archivos RDS de resultados
results_sp <- readRDS("/home/db/Documentos/David_Hito2/Correlations/Rds/results_sp.rds")
results_xp <- readRDS("/home/db/Documentos/David_Hito2/Correlations/Rds/results_xp.rds")
results_vina <- readRDS("/home/db/Documentos/David_Hito2/Correlations/Rds/results_vina.rds")
results_ad4 <- readRDS("/home/db/Documentos/David_Hito2/Correlations/Rds/results_ad4.rds")

# Combinar todos los datos
all_results <- bind_rows(results_sp, results_xp, results_vina, results_ad4)

# Filtrar solo pIC50 (más interpretable que IC50)
results_pic50 <- all_results %>% 
  filter(ExpVariable == "pIC50")

cat("Total de combinaciones analizadas:", nrow(results_pic50), "\n")
cat("Métodos:", paste(unique(results_pic50$Method), collapse = ", "), "\n\n")

# -------------------------
# 2. HEATMAP: Mejor correlación por Proteína x Método
# -------------------------

# Encontrar la mejor correlación (Pearson) para cada combinación Proteína-Método
best_by_protein_method <- results_pic50 %>%
  group_by(Protein, Method) %>%
  slice_max(abs(Pearson), n = 1, with_ties = FALSE) %>%
  ungroup()

p1 <- ggplot(best_by_protein_method, aes(x = Method, y = Protein, fill = Pearson)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f", Pearson)), color = "white", fontface = "bold", size = 5) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, 
                       limits = c(-1, 1), name = "Pearson r") +
  labs(title = "Mejor Correlación por Proteína y Método",
       subtitle = "Valor más alto (abs) entre todos los frames, poses y energy terms",
       x = "Método de Docking", y = "Proteína") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

print(p1)

#ggsave("01_Heatmap_Best_Correlation_by_Protein_Method.png", p1, width = 8, height = 6, dpi = 300)
cat("✓ Guardado: 01_Heatmap_Best_Correlation_by_Protein_Method.png\n")

# -------------------------
# 3. BOXPLOT: Distribución de correlaciones por Método
# -------------------------

p2 <- ggplot(results_pic50 %>% filter(!is.na(Pearson)), 
             aes(x = Method, y = Pearson, fill = Method)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_viridis_d(option = "plasma", begin = 0.2, end = 0.8) +
  labs(title = "Distribución de Correlaciones de Pearson por Método",
       subtitle = "Todos los frames, poses y energy terms",
       x = "Método de Docking", y = "Coeficiente de Pearson") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "none"
  )

print(p2)

#ggsave("02_Boxplot_Pearson_by_Method.png", p2, width = 10, height = 6, dpi = 300)
cat("✓ Guardado: 02_Boxplot_Pearson_by_Method.png\n")

# -------------------------
# 4. VIOLIN PLOT: Pearson vs Spearman por Método
# -------------------------

correlation_long <- results_pic50 %>%
  filter(!is.na(Pearson) & !is.na(Spearman)) %>%
  select(Method, Protein, Frame, Pearson, Spearman) %>%
  pivot_longer(cols = c(Pearson, Spearman), names_to = "Correlation_Type", values_to = "Value")

p3 <- ggplot(correlation_long, aes(x = Method, y = Value, fill = Correlation_Type)) +
  geom_violin(alpha = 0.7, position = position_dodge(width = 0.9)) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = c("Pearson" = "#E64B35", "Spearman" = "#4DBBD5"),
                    name = "Tipo de\nCorrelación") +
  labs(title = "Comparación: Pearson vs Spearman por Método",
       x = "Método de Docking", y = "Coeficiente de Correlación") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

print(p3)

#ggsave("03_Violin_Pearson_vs_Spearman.png", p3, width = 12, height = 6, dpi = 300)
cat("✓ Guardado: 03_Violin_Pearson_vs_Spearman.png\n")

# -------------------------
# 5. LÍNEA: Correlación vs Número de Poses
# -------------------------

corr_by_pose <- results_pic50 %>%
  filter(!is.na(Pearson)) %>%
  group_by(Method, Pose) %>%
  summarise(
    Mean_Pearson = mean(Pearson, na.rm = TRUE),
    SD_Pearson = sd(Pearson, na.rm = TRUE),
    .groups = "drop"
  )

p4 <- ggplot(corr_by_pose, aes(x = Pose, y = Mean_Pearson, color = Method, group = Method)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = Mean_Pearson - SD_Pearson, 
                  ymax = Mean_Pearson + SD_Pearson, fill = Method), 
              alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_viridis_d(option = "plasma", begin = 0.2, end = 0.8) +
  scale_fill_viridis_d(option = "plasma", begin = 0.2, end = 0.8) +
  labs(title = "Efecto del Número de Poses en la Correlación",
       subtitle = "Media ± SD del coeficiente de Pearson",
       x = "Número de Poses", y = "Pearson (promedio)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  )

print(p4)

#ggsave("04_Line_Correlation_vs_Poses.png", p4, width = 10, height = 6, dpi = 300)
cat("✓ Guardado: 04_Line_Correlation_vs_Poses.png\n")

# -------------------------
# 6. FACET HEATMAP: Data Fusion Methods por Método de Docking
# -------------------------

fusion_comparison <- results_pic50 %>%
  filter(!is.na(Pearson)) %>%
  group_by(Method, DataFusion) %>%
  summarise(Mean_Pearson = mean(Pearson, na.rm = TRUE), .groups = "drop")

p5 <- ggplot(fusion_comparison, aes(x = DataFusion, y = Method, fill = Mean_Pearson)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f", Mean_Pearson)), color = "white", fontface = "bold") +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0,
                       limits = c(-1, 1), name = "Pearson\n(promedio)") +
  labs(title = "Rendimiento de Métodos de Data Fusion",
       subtitle = "Promedio de correlación de Pearson",
       x = "Método de Data Fusion", y = "Método de Docking") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p5)

#ggsave("05_Heatmap_DataFusion_Methods.png", p5, width = 10, height = 6, dpi = 300)
cat("✓ Guardado: 05_Heatmap_DataFusion_Methods.png\n")

# -------------------------
# 7. BAR PLOT: Top 10 Mejores Combinaciones
# -------------------------

top_combinations <- results_pic50 %>%
  filter(!is.na(Pearson)) %>%
  arrange(desc(abs(Pearson))) %>%
  head(10) %>%
  mutate(
    Combination = paste0(Protein, " | ", Method, " | ", DataFusion, 
                         "\nFrame:", Frame, " | Pose:", Pose, " | ", EnergyTerm),
    Combination = factor(Combination, levels = rev(Combination))
  )

p6 <- ggplot(top_combinations, aes(x = Pearson, y = Combination, fill = Pearson)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.3f", Pearson)), hjust = -0.1, size = 3.5) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
  labs(title = "Top 10 Mejores Combinaciones",
       subtitle = "Mayor correlación absoluta (Pearson) con pIC50",
       x = "Coeficiente de Pearson", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "none",
    axis.text.y = element_text(size = 8)
  ) +
  xlim(c(-1, max(top_combinations$Pearson) + 0.1))

print(p6)

#ggsave("06_Barplot_Top10_Combinations.png", p6, width = 14, height = 8, dpi = 300)
cat("✓ Guardado: 06_Barplot_Top10_Combinations.png\n")

# -------------------------
# 8. SCATTER: Pearson vs Spearman (para ver concordancia)
# -------------------------

p7 <- ggplot(results_pic50 %>% filter(!is.na(Pearson) & !is.na(Spearman)), 
             aes(x = Pearson, y = Spearman, color = Method)) +
  geom_point(alpha = 0.4, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray30") +
  scale_color_viridis_d(option = "plasma", begin = 0.2, end = 0.8) +
  labs(title = "Concordancia: Pearson vs Spearman",
       subtitle = "Línea diagonal = concordancia perfecta",
       x = "Coeficiente de Pearson", y = "Coeficiente de Spearman") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  ) +
  coord_fixed()

print(p7)

#ggsave("07_Scatter_Pearson_vs_Spearman.png", p7, width = 10, height = 8, dpi = 300)
cat("✓ Guardado: 07_Scatter_Pearson_vs_Spearman.png\n")

# -------------------------
# 9. RESUMEN ESTADÍSTICO
# -------------------------

summary_stats <- results_pic50 %>%
  filter(!is.na(Pearson)) %>%
  group_by(Method) %>%
  summarise(
    N = n(),
    Mean_Pearson = mean(Pearson),
    SD_Pearson = sd(Pearson),
    Median_Pearson = median(Pearson),
    Max_Pearson = max(Pearson),
    Min_Pearson = min(Pearson),
    .groups = "drop"
  )

write.csv(summary_stats, "Summary_Statistics.csv", row.names = FALSE)
cat("\n✓ Guardado: Summary_Statistics.csv\n")

print(summary_stats)

cat("\n✨ ¡Todas las visualizaciones generadas exitosamente!\n")