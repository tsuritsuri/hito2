################################################################################
# Script de Análisis de Energías de Docking (AutoDock) vs IC50 Experimentales
# Autor: Adaptado de scripts de ensemble docking
# Fecha: 2025
################################################################################

# Cargar librerías necesarias
library(dplyr)
library(ggplot2)
library(data.table)
library(openxlsx)
library(ggpubr)
library(scales)
library(psych)
library(flextable)
library(knitr)
library(kableExtra)

################################################################################
# 1. CONFIGURACIÓN DE RUTAS
################################################################################

# MODIFICAR ESTAS RUTAS SEGÚN TUS ARCHIVOS
path_energies <- "Resumen_AD4_Energies_all.csv"  # Tu CSV con energías
path_ic50 <- "datos_final.csv"            # Tu CSV con IC50

################################################################################
# 2. LECTURA DE DATOS
################################################################################

cat("Leyendo datos de energías...\n")
# Leer CSV de energías de AutoDock
energy_data <- fread(path_energies)

cat("Estructura de datos de energías:\n")
print(head(energy_data))
print(str(energy_data))

cat("\nLeyendo datos de IC50...\n")
# Leer CSV de IC50 experimentales
# Columna 1: Compound_ID (con el número después del guión bajo)
# Columna 7: Value_nM (IC50 en nM)
ic50_data <- fread(path_ic50)

# Renombrar columnas para claridad
names(ic50_data)[1] <- "Compound_ID"
names(ic50_data)[7] <- "IC50_nM"

# Seleccionar solo las columnas necesarias
ic50_data <- ic50_data %>%
  select(Compound_ID, IC50_nM)

cat("Estructura de datos de IC50:\n")
print(head(ic50_data))
print(str(ic50_data))

################################################################################
# 3. PREPROCESAMIENTO DE DATOS
################################################################################

cat("\n=== PREPROCESAMIENTO ===\n")

# Extraer el ID del ligando (número después del guión bajo)
# Ejemplo: "100_155547948.0" -> "155547948.0"
energy_data <- energy_data %>%
  mutate(
    Ligand_ID = sub(".*_", "", Ligand),  # Extrae todo después del primer "_"
    Ligand_Full = Ligand                  # Mantiene el nombre completo
  )

cat("Ejemplo de extracción de IDs:\n")
print(energy_data %>% 
        select(Ligand_Full, Ligand_ID) %>% 
        distinct() %>% 
        head(10))

# Limpiar formato de IDs para asegurar coincidencia
# Remover el .0 al final si existe en ambos datasets
energy_data$Ligand_ID <- sub("\\.0$", "", energy_data$Ligand_ID)
ic50_data$Compound_ID <- sub("\\.0$", "", as.character(ic50_data$Compound_ID))

# Convertir a character para asegurar coincidencia
energy_data$Ligand_ID <- as.character(energy_data$Ligand_ID)
ic50_data$Compound_ID <- as.character(ic50_data$Compound_ID)

# Verificar coincidencias
ligands_energy <- unique(energy_data$Ligand_ID)
ligands_ic50 <- unique(ic50_data$Compound_ID)

cat("\nNúmero de ligandos en datos de energía:", length(ligands_energy), "\n")
cat("Número de ligandos en datos de IC50:", length(ligands_ic50), "\n")

common_ligands <- intersect(ligands_energy, ligands_ic50)
cat("Número de ligandos en común:", length(common_ligands), "\n")

if(length(common_ligands) == 0) {
  cat("\n¡ADVERTENCIA! No hay ligandos en común.\n")
  cat("Primeros 10 IDs en energías:", head(ligands_energy, 10), "\n")
  cat("Primeros 10 IDs en IC50:", head(ligands_ic50, 10), "\n")
  stop("Verifica el formato de los nombres de ligandos en ambos archivos.")
}

# Mostrar algunos ejemplos de IDs que coinciden
cat("\nEjemplos de IDs que coinciden:\n")
print(head(common_ligands, 10))

################################################################################
# 4. CÁLCULO DE MÉTRICAS DE ENERGÍA POR LIGANDO
################################################################################

cat("\n=== CALCULANDO MÉTRICAS DE ENERGÍA ===\n")

# Filtrar solo ligandos que tienen IC50
energy_data_filtered <- energy_data %>%
  filter(Ligand_ID %in% common_ligands)

# Calcular estadísticas descriptivas por ligando
energy_metrics <- energy_data_filtered %>%
  group_by(Ligand_ID, Protein, Frame) %>%
  summarise(
    Min_Energy = min(Energy, na.rm = TRUE),
    Max_Energy = max(Energy, na.rm = TRUE),
    Mean_Energy = mean(Energy, na.rm = TRUE),
    Median_Energy = median(Energy, na.rm = TRUE),
    SD_Energy = sd(Energy, na.rm = TRUE),
    # Media Geométrica (adaptada para valores negativos)
    Geometric_Mean = sign(mean(Energy)) * exp(mean(log(abs(Energy)))),
    # Media Armónica (adaptada para valores negativos)
    Harmonic_Mean = sign(mean(Energy)) * (n() / sum(1 / abs(Energy))),
    N_Poses = n(),
    .groups = "drop"
  )

# Calcular métricas globales por ligando (promedio entre frames/proteínas)
energy_summary <- energy_metrics %>%
  group_by(Ligand_ID) %>%
  summarise(
    # Mejor energía de todas las poses y frames
    Best_Energy = min(Min_Energy, na.rm = TRUE),
    # Promedio de las energías mínimas por frame
    Avg_Min_Energy = mean(Min_Energy, na.rm = TRUE),
    # Promedio global de todas las energías
    Global_Mean_Energy = mean(Mean_Energy, na.rm = TRUE),
    # Media geométrica global
    Global_Geometric_Mean = mean(Geometric_Mean, na.rm = TRUE),
    # Media armónica global
    Global_Harmonic_Mean = mean(Harmonic_Mean, na.rm = TRUE),
    # Mediana global
    Global_Median = mean(Median_Energy, na.rm = TRUE),
    # Número total de frames/proteínas
    N_Frames = n(),
    # Número total de poses
    Total_Poses = sum(N_Poses),
    .groups = "drop"
  )

cat("Métricas calculadas para", nrow(energy_summary), "ligandos\n")
print(head(energy_summary))

################################################################################
# 5. MERGE CON DATOS DE IC50
################################################################################

cat("\n=== FUSIONANDO DATOS ===\n")

# Unir datos de energía con IC50
merged_data <- energy_summary %>%
  inner_join(ic50_data, by = c("Ligand_ID" = "Compound_ID"))

cat("Datos fusionados exitosamente:", nrow(merged_data), "ligandos\n")
print(head(merged_data))

# Guardar datos fusionados
write.csv(merged_data, "merged_energy_ic50.csv", row.names = FALSE)
cat("Datos guardados en: merged_energy_ic50.csv\n")

################################################################################
# 6. ANÁLISIS DE CORRELACIÓN
################################################################################

cat("\n=== ANÁLISIS DE CORRELACIÓN ===\n")

# Calcular correlaciones entre diferentes métricas de energía y IC50
correlation_results <- data.frame(
  Metric = c("Best_Energy", "Avg_Min_Energy", "Global_Mean_Energy", 
             "Global_Geometric_Mean", "Global_Harmonic_Mean", "Global_Median"),
  Pearson_r = NA,
  Pearson_p = NA,
  Spearman_rho = NA,
  Spearman_p = NA
)

metrics_to_test <- c("Best_Energy", "Avg_Min_Energy", "Global_Mean_Energy", 
                     "Global_Geometric_Mean", "Global_Harmonic_Mean", "Global_Median")

for (i in 1:length(metrics_to_test)) {
  metric <- metrics_to_test[i]
  
  # Pearson
  pearson <- cor.test(merged_data[[metric]], merged_data$IC50_nM, method = "pearson")
  correlation_results$Pearson_r[i] <- round(pearson$estimate, 4)
  correlation_results$Pearson_p[i] <- round(pearson$p.value, 4)
  
  # Spearman
  spearman <- cor.test(merged_data[[metric]], merged_data$IC50_nM, method = "spearman")
  correlation_results$Spearman_rho[i] <- round(spearman$estimate, 4)
  correlation_results$Spearman_p[i] <- round(spearman$p.value, 4)
}

cat("\nResultados de Correlación:\n")
print(correlation_results)

# Guardar tabla de correlaciones
write.xlsx(correlation_results, "correlation_results.xlsx")

# Crear tabla formateada
flex_corr <- correlation_results %>%
  flextable() %>%
  set_table_properties(layout = "autofit") %>%
  colformat_double(digits = 4)

print(flex_corr)

################################################################################
# 7. VISUALIZACIONES
################################################################################

cat("\n=== GENERANDO GRÁFICOS ===\n")

# Función para crear scatter plots con correlación
create_correlation_plot <- function(data, x_var, y_var, x_label, y_label, title) {
  # Calcular correlación
  cor_pearson <- cor.test(data[[x_var]], data[[y_var]], method = "pearson")
  cor_spearman <- cor.test(data[[x_var]], data[[y_var]], method = "spearman")
  
  # Etiqueta con resultados de correlación
  cor_label <- paste0(
    "Pearson r = ", round(cor_pearson$estimate, 3), 
    " (p = ", format.pval(cor_pearson$p.value, digits = 3), ")\n",
    "Spearman ρ = ", round(cor_spearman$estimate, 3),
    " (p = ", format.pval(cor_spearman$p.value, digits = 3), ")"
  )
  
  # Crear plot
  p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(alpha = 0.6, size = 3, color = "#2C5F8D") +
    geom_smooth(method = "lm", se = TRUE, color = "#D55E00", 
                fill = "#D55E00", alpha = 0.2) +
    labs(
      title = title,
      x = x_label,
      y = y_label,
      subtitle = cor_label
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(size = 10, color = "black"),
      panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
    )
  
  return(p)
}

# 1. Best Energy vs IC50
p1 <- create_correlation_plot(
  merged_data, 
  "Best_Energy", 
  "IC50_nM",
  "Best Docking Energy (kcal/mol)",
  "IC50 (nM)",
  "Correlation: Best Docking Energy vs IC50"
)

ggsave("plot_best_energy_vs_ic50.png", p1, width = 10, height = 8, dpi = 300)

# 2. Global Geometric Mean vs IC50
p2 <- create_correlation_plot(
  merged_data, 
  "Global_Geometric_Mean", 
  "IC50_nM",
  "Geometric Mean Energy (kcal/mol)",
  "IC50 (nM)",
  "Correlation: Geometric Mean Energy vs IC50"
)

ggsave("plot_geometric_mean_vs_ic50.png", p2, width = 10, height = 8, dpi = 300)

# 3. Avg Min Energy vs IC50
p3 <- create_correlation_plot(
  merged_data, 
  "Avg_Min_Energy", 
  "IC50_nM",
  "Average Minimum Energy (kcal/mol)",
  "IC50 (nM)",
  "Correlation: Average Minimum Energy vs IC50"
)

ggsave("plot_avg_min_energy_vs_ic50.png", p3, width = 10, height = 8, dpi = 300)

# 4. Panel con todas las métricas
library(gridExtra)

p4 <- create_correlation_plot(merged_data, "Global_Mean_Energy", "IC50_nM",
                              "Mean Energy", "IC50 (nM)", "Mean Energy")
p5 <- create_correlation_plot(merged_data, "Global_Harmonic_Mean", "IC50_nM",
                              "Harmonic Mean", "IC50 (nM)", "Harmonic Mean")
p6 <- create_correlation_plot(merged_data, "Global_Median", "IC50_nM",
                              "Median Energy", "IC50 (nM)", "Median Energy")

panel_plot <- grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 2)

ggsave("plot_all_correlations_panel.png", panel_plot, width = 16, height = 20, dpi = 300)

cat("\nGráficos guardados exitosamente\n")

################################################################################
# 8. RANKING DE LIGANDOS
################################################################################

cat("\n=== GENERANDO RANKINGS ===\n")

# Crear rankings para cada métrica
ranking_data <- merged_data %>%
  mutate(
    Rank_Best_Energy = rank(Best_Energy),
    Rank_Geometric_Mean = rank(Global_Geometric_Mean),
    Rank_Avg_Min = rank(Avg_Min_Energy),
    Rank_Median = rank(Global_Median),
    Rank_IC50 = rank(IC50_nM)
  ) %>%
  arrange(Best_Energy)

# Guardar rankings
write.xlsx(ranking_data, "ligand_rankings.xlsx")

cat("Rankings guardados en: ligand_rankings.xlsx\n")

################################################################################
# 9. RESUMEN ESTADÍSTICO
################################################################################

cat("\n=== RESUMEN ESTADÍSTICO FINAL ===\n")

summary_stats <- merged_data %>%
  summarise(
    N_Ligands = n(),
    Energy_Min = min(Best_Energy, na.rm = TRUE),
    Energy_Max = max(Best_Energy, na.rm = TRUE),
    Energy_Mean = mean(Best_Energy, na.rm = TRUE),
    Energy_SD = sd(Best_Energy, na.rm = TRUE),
    IC50_Min = min(IC50_nM, na.rm = TRUE),
    IC50_Max = max(IC50_nM, na.rm = TRUE),
    IC50_Mean = mean(IC50_nM, na.rm = TRUE),
    IC50_SD = sd(IC50_nM, na.rm = TRUE)
  )

cat("\nEstadísticas Descriptivas:\n")
print(summary_stats)

# Tabla formateada
flex_summary <- summary_stats %>%
  tidyr::pivot_longer(everything(), names_to = "Statistic", values_to = "Value") %>%
  flextable() %>%
  set_table_properties(layout = "autofit") %>%
  colformat_double(j = "Value", digits = 3)

print(flex_summary)

################################################################################
# 10. EXPORTAR TODO EN UN EXCEL CON MÚLTIPLES HOJAS
################################################################################

cat("\n=== EXPORTANDO RESULTADOS FINALES ===\n")

# Crear workbook
wb <- createWorkbook()

# Hoja 1: Datos fusionados
addWorksheet(wb, "Merged_Data")
writeData(wb, "Merged_Data", merged_data)

# Hoja 2: Correlaciones
addWorksheet(wb, "Correlations")
writeData(wb, "Correlations", correlation_results)

# Hoja 3: Rankings
addWorksheet(wb, "Rankings")
writeData(wb, "Rankings", ranking_data)

# Hoja 4: Resumen estadístico
addWorksheet(wb, "Summary_Stats")
writeData(wb, "Summary_Stats", summary_stats)

# Hoja 5: Métricas por frame
addWorksheet(wb, "Energy_Metrics_PerFrame")
writeData(wb, "Energy_Metrics_PerFrame", energy_metrics)

# Guardar workbook
saveWorkbook(wb, "complete_analysis_results.xlsx", overwrite = TRUE)

cat("\n¡ANÁLISIS COMPLETADO!\n")
cat("Archivo principal guardado en: complete_analysis_results.xlsx\n")
cat("\nArchivos generados:\n")
cat("  - merged_energy_ic50.csv\n")
cat("  - correlation_results.xlsx\n")
cat("  - ligand_rankings.xlsx\n")
cat("  - complete_analysis_results.xlsx\n")
cat("  - plot_best_energy_vs_ic50.png\n")
cat("  - plot_geometric_mean_vs_ic50.png\n")
cat("  - plot_avg_min_energy_vs_ic50.png\n")
cat("  - plot_all_correlations_panel.png\n")

################################################################################
# FIN DEL SCRIPT
################################################################################