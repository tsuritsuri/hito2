# ---------------------------------------------
# LIBRERÍAS
# ---------------------------------------------
library(stringr)
library(dplyr)
library(readr)
library(purrr)
library(tibble)
library(fs)
library(furrr)
library(future)
library(progressr)

# ---------------------------------------------
# DIRECTORIO RAÍZ
# Estructura: root/<protein>/<frame>/*.log
# Ej: /home/db/.../AutoDock Vina/1E9Y/5/ligand_name.log
# ---------------------------------------------
root_dir <- "/home/db/Documentos/David_Hito2/AutoDock Vina/vina_out"

# Listar todos los archivos .log
all_log_files <- dir_ls(root_dir, recurse = TRUE, type = "file", glob = "*.log")

if (length(all_log_files) == 0) {
  stop("No se encontraron archivos .log en: ", root_dir)
}

cat("Encontrados", length(all_log_files), "archivos .log\n")

# ---------------------------------------------
# EXTRACCIÓN DE PROTEÍNA, FRAME Y LIGANDO
# ---------------------------------------------
parse_metadata <- function(log_file) {
  # Estructura: .../<protein>/<frame>/ligand.log
  log_dir     <- path_dir(log_file)           # .../<protein>/<frame>
  protein_dir <- path_dir(log_dir)            # .../<protein>
  
  protein_id <- path_file(protein_dir)
  frame_raw  <- path_file(log_dir)
  ligand_name <- path_ext_remove(path_file(log_file))
  
  # frame puede ser "1", "526", etc.
  frame_id <- suppressWarnings(as.integer(frame_raw))
  if (is.na(frame_id)) {
    frame_id <- suppressWarnings(as.integer(str_match(frame_raw, ".*?(\\d+)")[,2]))
  }
  
  list(
    protein_id = protein_id, 
    frame_id = frame_id,
    ligand_name = ligand_name
  )
}

# ---------------------------------------------
# FUNCIÓN PARA PROCESAR UN .log
# ---------------------------------------------
process_vina_log <- function(log_file) {
  metadata <- parse_metadata(log_file)
  
  # Leer el archivo
  log_text <- read_file(log_file)
  log_lines <- str_split(log_text, "\n")[[1]]
  
  # Buscar la tabla de resultados
  # Empieza después de la línea "-----+------------+----------+----------"
  separator_idx <- which(str_detect(log_lines, "^\\s*-----\\+"))
  
  if (length(separator_idx) == 0) {
    warning("No se encontró tabla de resultados en: ", log_file)
    return(tibble(
      Protein = character(0),
      Frame = integer(0),
      Ligand = character(0),
      Pose = integer(0),
      Energy = numeric(0),
      log_path = character(0)
    ))
  }
  
  # Tomar líneas después del separador
  start_idx <- separator_idx[1] + 1
  data_lines <- log_lines[start_idx:length(log_lines)]
  
  # Filtrar líneas que tienen datos (empiezan con espacios y números)
  data_lines <- data_lines[str_detect(data_lines, "^\\s*\\d+")]
  
  if (length(data_lines) == 0) {
    warning("No se encontraron datos de poses en: ", log_file)
    return(tibble(
      Protein = character(0),
      Frame = integer(0),
      Ligand = character(0),
      Pose = integer(0),
      Energy = numeric(0),
      log_path = character(0)
    ))
  }
  
  # Parsear cada línea
  parsed_data <- map_dfr(data_lines, function(line) {
    # Limpiar espacios múltiples y separar por espacios
    parts <- str_squish(line) %>% str_split("\\s+") %>% .[[1]]
    
    if (length(parts) >= 2) {
      mode <- suppressWarnings(as.integer(parts[1]))
      affinity <- suppressWarnings(as.numeric(parts[2]))
      
      if (!is.na(mode) && !is.na(affinity)) {
        return(tibble(
          Pose = mode,
          Energy = affinity
        ))
      }
    }
    return(NULL)
  })
  
  if (nrow(parsed_data) == 0) {
    return(tibble(
      Protein = character(0),
      Frame = integer(0),
      Ligand = character(0),
      Pose = integer(0),
      Energy = numeric(0),
      log_path = character(0)
    ))
  }
  
  # Agregar metadata
  parsed_data %>%
    mutate(
      Protein = metadata$protein_id,
      Frame = metadata$frame_id,
      Ligand = metadata$ligand_name,
      log_path = log_file,
      .before = 1
    ) %>%
    select(Protein, Frame, Ligand, Pose, Energy, log_path)
}

# ---------------------------------------------
# CHEQUEO RÁPIDO (muestra cómo parsea 5 rutas)
# ---------------------------------------------
cat("\nEjemplos de parseo (Protein / Frame / Ligand):\n")
print(
  tibble(path = all_log_files[seq_len(min(5, length(all_log_files)))]) |>
    mutate(
      meta = map(path, parse_metadata),
      Protein = map_chr(meta, "protein_id"),
      Frame = map_int(meta, "frame_id"),
      Ligand = map_chr(meta, "ligand_name")
    ) |>
    select(Protein, Frame, Ligand, path)
)

# ---------------------------------------------
# PARALELISMO (usa todos los cores menos 5)
# ---------------------------------------------
n_workers <- max(1, future::availableCores() - 5)
if (.Platform$OS.type == "windows") {
  plan(multisession, workers = n_workers)
} else {
  plan(multicore, workers = n_workers)
}

cat("\nProcesando ", length(all_log_files), " archivos .log con ", n_workers, " workers...\n", sep = "")

handlers(global = TRUE)
tictoc <- system.time({
  all_results <- with_progress({
    p <- progressor(along = all_log_files)
    future_map_dfr(
      all_log_files,
      ~{ res <- process_vina_log(.x); p(); res },
      .options = furrr_options(seed = TRUE)
    )
  })
})

plan(sequential)

# ---------------------------------------------
# LIMPIEZA DE NOMBRES DE LIGANDOS
# Remover el prefijo "N_" (ej: "1_12345" -> "12345")
# ---------------------------------------------
cat("\nLimpiando nombres de ligandos...\n")
ligands_before <- head(all_results$Ligand, 5)
cat("Ejemplos ANTES de limpiar:\n")
print(ligands_before)

all_results_vina <- all_results %>%
  mutate(Ligand = str_replace(Ligand, "^\\d+_", ""))

ligands_after <- head(all_results$Ligand, 5)
cat("\nEjemplos DESPUÉS de limpiar:\n")
print(ligands_after)

# ---------------------------------------------
# REPORTE Y GUARDADO
# ---------------------------------------------
cat("\n--- Completado ---\n")
cat(sprintf("Tiempo: %.2f s | Filas: %d | Ligandos únicos: %d | Proteínas: %d | Frames: %d\n",
            tictoc[["elapsed"]], 
            nrow(all_results),
            n_distinct(all_results$Ligand),
            n_distinct(all_results$Protein),
            n_distinct(all_results$Frame)))

cat("\nDistribución de poses por ligando:\n")
print(table(all_results %>% count(Protein, Frame, Ligand) %>% pull(n)))

cat("\n=== PRIMERAS FILAS ===\n")
print(head(all_results, 20))

# Opcional: Guardar el dataframe
# write.csv(all_results, "vina_results.csv", row.names = FALSE)