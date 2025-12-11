# ==============================================================================
# SCRIPT PARA PARSEAR LOGS DE GNINA (CNN AFFINITY)
# ==============================================================================

# 1. Cargar librerías necesarias
if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)

# 2. Definir la ruta base donde está la carpeta "gnina_out"
# Asegúrate de que tu sesión de R (getwd()) apunte a la carpeta que contiene "gnina_out"
base_dir <- "gnina_out"

# 3. Listar todos los archivos .log recursivamente
cat("Buscando archivos .log en", base_dir, "...\n")
files <- list.files(path = base_dir, pattern = "\\.log$", recursive = TRUE, full.names = TRUE)

# ------------------------------------------------------------------------------
# Función principal de parseo
# ------------------------------------------------------------------------------
parse_gnina_log <- function(filepath) {
  
  # --- A. Extraer metadatos del nombre del archivo ---
  # Divide la ruta por "/"
  parts <- str_split(filepath, "/")[[1]]
  
  # Asumiendo estructura: ... / gnina_out / Protein / Frame / Archivo
  filename <- tail(parts, 1)
  frame    <- parts[length(parts)-1]
  protein  <- parts[length(parts)-2]
  
  # --- B. Limpiar nombre del ligando ---
  # 1. Quitar extensión .log
  # 2. Quitar todo lo que hay antes del primer guión bajo (incluyendo el guión)
  ligand_raw  <- str_remove(filename, "\\.log$")
  ligand_name <- sub("^[^_]*_", "", ligand_raw) 
  
  # --- C. Leer contenido del log ---
  lines <- readLines(filepath)
  
  # Buscar la línea del encabezado "mode | affinity"
  header_line_idx <- grep("mode \\|\\s+affinity", lines)
  
  # Si no encuentra encabezado, retorna NULL (archivo vacío o corrupto)
  if (length(header_line_idx) == 0) return(NULL)
  
  # Saltar encabezados:
  # 1. "mode | affinity..."
  # 2. "     | (kcal/mol)..."
  # 3. "-----+------------..."
  # Los datos empiezan en la línea 4 después del hallazgo, por eso sumamos 3
  start_data <- header_line_idx + 3 
  
  # Extraer posibles líneas de datos
  raw_data_lines <- lines[start_data:length(lines)]
  
  # --- FILTRO DE SEGURIDAD ---
  # Quedarse SOLO con líneas que empiezan con un número.
  # Esto elimina separadores, líneas vacías finales o mensajes de error.
  data_lines <- raw_data_lines[grepl("^[[:space:]]*[0-9]", raw_data_lines)]
  
  if (length(data_lines) == 0) return(NULL)
  
  # --- D. Convertir texto a Data Frame ---
  df_log <- read.table(text = data_lines, header = FALSE, fill = TRUE)
  
  # SELECCIÓN DE COLUMNAS (Basado en tu log):
  # V1: Mode (Pose)
  # V2: Affinity kcal/mol (IGNORAR)
  # V3: Intramol
  # V4: CNN pose score
  # V5: CNN affinity (SELECCIONAR ESTA)
  
  df_clean <- df_log %>%
    select(V1, V5) %>%                  # Nos quedamos con Pose y CNN Affinity
    rename(Pose = V1, Energy = V5) %>%  # Renombramos
    mutate(
      Protein  = protein,
      Frame    = frame,
      Ligand   = ligand_name,
      Filepath = filepath               # Agregamos la ruta del archivo
    ) %>%
    # Reordenar columnas para el formato final deseado
    select(Protein, Frame, Ligand, Pose, Energy, Filepath)
  
  return(df_clean)
}

# ------------------------------------------------------------------------------
# 4. Ejecutar el parseo
# ------------------------------------------------------------------------------
cat("Procesando", length(files), "archivos... esto puede tomar unos minutos.\n")

# map_dfr aplica la función a cada archivo y combina todo en un solo DataFrame
final_df <- map_dfr(files, parse_gnina_log)

# ------------------------------------------------------------------------------
# 5. Ajustar tipos de datos y Guardar
# ------------------------------------------------------------------------------
if (nrow(final_df) > 0) {
  
  final_df <- final_df %>%
    mutate(
      Protein  = as.factor(Protein),
      Frame    = as.numeric(Frame),
      Pose     = as.integer(Pose),
      Energy   = as.numeric(Energy), # Esta es la CNN Affinity
      Filepath = as.character(Filepath)
    )
  
  # Guardar resultado
  saveRDS(final_df, file = "GNINA_CNN.rds")
  
  cat("¡Proceso terminado con éxito!\n")
  cat("Archivo guardado como: Gnina_Parsed_CNN.rds\n")
  print(head(final_df))
  
  # Verificar que no hay NAs
  n_nas <- sum(is.na(final_df$Energy))
  if(n_nas == 0) {
    cat("Verificación: No se encontraron valores NA en la energía.\n")
  } else {
    cat("Advertencia: Se encontraron", n_nas, "valores NA.\n")
  }
  
} else {
  cat("Error: No se pudieron extraer datos. Verifica las rutas.\n")
}