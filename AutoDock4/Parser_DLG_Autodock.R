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
# Estructura: root/<protein>/<frame>/dlg/*.dlg
# Ej: /.../autodock_out_pl4/1E9Y/1/dlg/XXXX.dlg
# ---------------------------------------------
root_dir <- "/home/db/Documentos/David_Hito2/AutoDock4/autodock_out2"

# listar .dlg SOLO si están bajo una carpeta llamada 'dlg'
all_files <- dir_ls(root_dir, recurse = TRUE, type = "file", glob = "*.dlg")
dlg_files <- all_files[path_file(path_dir(all_files)) == "dlg"]

if (length(dlg_files) == 0) stop("No se encontraron .dlg bajo '.../<protein>/<frame>/dlg/' en: ", root_dir)

# ---------------------------------------------
# EXTRACCIÓN ROBUSTA DE PROTEÍNA Y FRAME
# (siempre: filename -> dlg -> <frame> -> <protein> )
# ---------------------------------------------
parse_protein_frame <- function(dlg_file) {
  dlg_dir    <- path_dir(dlg_file)          # .../<frame>/dlg
  frame_dir  <- path_dir(dlg_dir)           # .../<frame>
  protein_dir<- path_dir(frame_dir)         # .../<protein>
  
  protein_id <- path_file(protein_dir)
  frame_raw  <- path_file(frame_dir)
  
  # frame puede ser "1", "526" o "frame526" -> sacar dígitos si hace falta
  frame_id <- suppressWarnings(as.integer(frame_raw))
  if (is.na(frame_id)) {
    frame_id <- suppressWarnings(as.integer(str_match(frame_raw, ".*?(\\d+)")[,2]))
  }
  
  list(protein_id = protein_id, frame_id = frame_id)
}

# ---------------------------------------------
# FUNCIÓN PARA PROCESAR UN .dlg
# ---------------------------------------------
process_dlg <- function(dlg_file) {
  ids <- parse_protein_frame(dlg_file)
  protein_id <- ids$protein_id
  frame_id   <- ids$frame_id
  
  # leer
  dlg_text  <- read_file(dlg_file)
  dlg_lines <- str_split(dlg_text, "\n")[[1]]
  
  # nombre del ligando (fallback: nombre de archivo)
  ligand_line <- str_extract(dlg_text, "INPUT-LIGAND-PDBQT: REMARK\\s+Name\\s*=.*")
  ligand_name <- if (!is.na(ligand_line)) {
    str_match(ligand_line, "Name\\s*=\\s*(.*)$")[,2]
  } else {
    tools::file_path_sans_ext(path_file(dlg_file))
  }
  
  # líneas 'RANKING'
  ranking_lines <- dlg_lines[str_detect(dlg_lines, "RANKING\\s*$")]
  if (length(ranking_lines) == 0) {
    return(tibble(
      Protein = character(0),
      Frame   = integer(0),
      Ligand  = character(0),
      Pose    = integer(0),
      Energy  = numeric(0),
      DLG_Path = character(0)
    ))
  }
  
  split_lines <- str_split(str_squish(ranking_lines), "\\s+")
  energies <- suppressWarnings(sapply(split_lines, function(x) as.numeric(x[4])))
  
  tibble(
    Protein = protein_id,
    Frame   = frame_id,
    Ligand  = ligand_name,
    Pose    = seq_along(energies),
    Energy  = energies,
    DLG_Path = dlg_file
  )
}

# ---------------------------------------------
# CHEQUEO RÁPIDO (muestra cómo parsea 5 rutas)
# ---------------------------------------------
cat("Ejemplos de parseo (Protein / Frame):\n")
print(
  tibble(path = dlg_files[seq_len(min(5, length(dlg_files)))]) |>
    mutate(
      Protein = map_chr(path, ~parse_protein_frame(.x)$protein_id),
      Frame   = map_int(path,  ~parse_protein_frame(.x)$frame_id)
    ) |>
    select(Protein, Frame, path)
)

# ---------------------------------------------
# PARALELISMO (usa todos los cores menos 1)
# ---------------------------------------------
n_workers <- max(1, future::availableCores() - 5)
if (.Platform$OS.type == "windows") plan(multisession, workers = n_workers) else plan(multicore, workers = n_workers)

cat("\nProcesando ", length(dlg_files), " .dlg con ", n_workers, " workers...\n", sep = "")

handlers(global = TRUE)
tictoc <- system.time({
  all_results <- with_progress({
    p <- progressor(along = dlg_files)
    future_map_dfr(
      dlg_files,
      ~{ res <- process_dlg(.x); p(); res },
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

all_results_ad4 <- all_results %>%
  mutate(Ligand = str_replace(Ligand, "^\\d+_", ""))

ligands_after <- head(all_results$Ligand, 5)
cat("\nEjemplos DESPUÉS de limpiar:\n")
print(ligands_after)

# ---------------------------------------------
# REPORTE Y GUARDADO
# ---------------------------------------------
cat("\n--- Completado ---\n")
cat(sprintf("Tiempo: %.2f s | Filas: %d | Ligandos únicos: %d | Proteínas: %d | Frames: %d\n",
            tictoc[["elapsed"]], nrow(all_results),
            n_distinct(all_results$Ligand),
            n_distinct(all_results$Protein),
            n_distinct(all_results$Frame)))

print(head(all_results, 10))

