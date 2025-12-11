library(tidyverse)

# Definir las columnas de energía y necesarias
energy_cols <- c(
  "docking.score",
  "glide.ligand.efficiency",
  "glide.gscore",
  "glide.emodel",
  "glide.energy",
  "glide.einternal"
)
needed_cols <- c("Title", energy_cols, "IC50")

# Función para parsear un archivo CSV
parse_docking_file <- function(filepath) {
  tryCatch({
    # Extraer información del nombre del archivo
    # Patrón: PROTEIN_FRAME_MODE.csv (ej: 1E9Y_6_SP.csv)
    filename <- basename(filepath)
    filename_parts <- str_match(filename, "^([^_]+)_(\\d+)_(SP|XP)\\.csv$")
    
    if (is.na(filename_parts[1])) {
      warning(paste("Nombre de archivo no sigue el patrón esperado:", filename))
      return(NULL)
    }
    
    protein_id <- filename_parts[2]
    frame_number <- as.integer(filename_parts[3])
    category <- filename_parts[4]
    
    # Leer el CSV
    df <- read.csv(filepath, stringsAsFactors = FALSE)
    
    # Verificar y seleccionar solo las columnas que existen
    available_cols <- intersect(needed_cols, colnames(df))
    
    if (length(available_cols) == 0) {
      warning(paste("No se encontraron columnas necesarias en:", filepath))
      return(NULL)
    }
    
    df_subset <- df[, available_cols, drop = FALSE]
    
    # Forzar tipos de datos correctos para evitar conflictos
    if ("Title" %in% colnames(df_subset)) {
      df_subset$Title <- as.character(df_subset$Title)
    }
    if ("IC50" %in% colnames(df_subset)) {
      df_subset$IC50 <- as.numeric(df_subset$IC50)
    }
    
    # Convertir columnas de energía a numérico
    for (col in energy_cols) {
      if (col %in% colnames(df_subset)) {
        df_subset[[col]] <- as.numeric(df_subset[[col]])
      }
    }
    
    # Agregar columna Pose ANTES de agregar metadatos
    # Ordenar por Title y docking.score (más negativo primero)
    # y asignar número de pose dentro de cada ligando
    if ("docking.score" %in% colnames(df_subset) && "Title" %in% colnames(df_subset)) {
      df_subset <- df_subset %>%
        arrange(Title, docking.score) %>%
        group_by(Title) %>%
        mutate(Pose = row_number()) %>%
        ungroup()
    } else {
      # Si no existe docking.score, igual asignar Pose = 1 a todo
      df_subset$Pose <- 1
      warning(paste("docking.score no encontrado en:", filepath, "- Asignando Pose = 1"))
    }
    
    # Agregar información de contexto
    df_subset$protein_id <- protein_id
    df_subset$frame <- frame_number
    df_subset$category <- category
    df_subset$source_file <- filename
    
    return(df_subset)
    
  }, error = function(e) {
    warning(paste("Error al procesar", filepath, ":", e$message))
    return(NULL)
  })
}

# Función principal para parsear todos los CSVs en ALL_CSV
parse_all_docking_data <- function(base_path = ".", folder_name = "ALL_CSV") {
  
  folder_path <- file.path(base_path, folder_name)
  
  # Verificar si la carpeta existe
  if (!dir.exists(folder_path)) {
    stop(paste("Carpeta no encontrada:", folder_path))
  }
  
  # Obtener todos los archivos CSV en la carpeta
  csv_files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)
  
  if (length(csv_files) == 0) {
    stop(paste("No se encontraron archivos CSV en:", folder_path))
  }
  
  cat("Encontrados", length(csv_files), "archivos CSV\n")
  
  # Lista para almacenar todos los dataframes
  all_data <- list()
  
  # Procesar cada archivo CSV con barra de progreso
  pb <- txtProgressBar(min = 0, max = length(csv_files), style = 3)
  
  for (i in seq_along(csv_files)) {
    csv_file <- csv_files[i]
    parsed_data <- parse_docking_file(csv_file)
    
    if (!is.null(parsed_data)) {
      all_data[[length(all_data) + 1]] <- parsed_data
    }
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  # Combinar todos los dataframes
  if (length(all_data) == 0) {
    stop("No se pudo parsear ningún archivo")
  }
  
  cat("\nCombinando datos...\n")
  combined_df <- bind_rows(all_data)
  
  # Remover filas con "minimized" en el Title
  rows_before <- nrow(combined_df)
  combined_df <- combined_df %>%
    filter(!grepl("minimized", Title, ignore.case = TRUE))
  rows_after <- nrow(combined_df)
  
  if (rows_before > rows_after) {
    cat("Removidas", rows_before - rows_after, "filas con 'minimized' en el Title\n")
  }
  
  # Reordenar columnas para mejor visualización
  # Incluir frame después de protein_id
  col_order <- c("protein_id", "frame", "category", "source_file", "Title", "Pose",
                 energy_cols[energy_cols %in% colnames(combined_df)], 
                 "IC50")
  col_order <- col_order[col_order %in% colnames(combined_df)]
  
  combined_df <- combined_df[, col_order]
  
  return(combined_df)
}

# Ejecutar el parser
# Ajusta base_path si tu carpeta ALL_CSV está en otro directorio
docking_data_glide <- parse_all_docking_data(base_path = "/home/db/Documentos/David_Hito2/Glide")

# Mostrar resumen
cat("\n=== RESUMEN DE DATOS ===\n")
cat("Total de filas:", nrow(docking_data_glide), "\n")
cat("\nDistribución por proteína:\n")
print(table(docking_data_glide$protein_id))
cat("\nDistribución por categoría:\n")
print(table(docking_data_glide$category))
cat("\nDistribución por proteína y categoría:\n")
print(table(docking_data_glide$protein_id, docking_data_glide$category))

# Resumen de frames
cat("\nRango de frames por proteína:\n")
frame_summary <- docking_data_glide %>%
  group_by(protein_id) %>%
  summarise(
    min_frame = min(frame),
    max_frame = max(frame),
    n_frames = n_distinct(frame),
    .groups = "drop"
  )
print(frame_summary)

# Resumen de ligandos con múltiples poses
ligands_multiple_poses <- docking_data_glide %>%
  group_by(protein_id, frame, category, Title) %>%
  summarise(n_poses = n(), .groups = "drop") %>%
  filter(n_poses > 1)

if (nrow(ligands_multiple_poses) > 0) {
  cat("\n=== LIGANDOS CON MÚLTIPLES POSES ===\n")
  cat("Total de ligandos con múltiples poses:", nrow(ligands_multiple_poses), "\n")
  cat("\nDistribución de número de poses:\n")
  print(table(ligands_multiple_poses$n_poses))
} else {
  cat("\n=== No hay ligandos con múltiples poses ===\n")
}

# Mostrar primeras filas
cat("\n=== PRIMERAS FILAS ===\n")
print(head(docking_data_glide, 10))

# Opcional: Guardar el dataframe combinado
# write.csv(docking_data, "docking_data_combined.csv", row.names = FALSE)