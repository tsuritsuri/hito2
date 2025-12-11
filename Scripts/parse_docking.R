#!/usr/bin/env Rscript
# ============================================================
# Parseo de salida AutoDock (resumen_resultados.txt) a CSV
# ============================================================

# ---- 0) Paquetes ----
suppressPackageStartupMessages({
  library(tidyverse)
})

# ---- Función principal de parseo ----
parse_docking_file <- function(txt_file, out_csv) {
  tryCatch({
    # Leer texto completo
    lines <- readLines(txt_file, warn = FALSE, encoding = "UTF-8")
    
    # Detectar bloques por ligando
    ligando_indices <- grep("^Procesando ligando:", lines)
    if (length(ligando_indices) == 0) {
      cat("No se encontraron ligandos en:", txt_file, "\n")
      return(FALSE)
    }
    
    ligando_names <- gsub("^Procesando ligando:\\s*", "", lines[ligando_indices])
    block_ends <- c(ligando_indices[-1] - 1, length(lines))
    
    # Parseo
    rows <- list()
    
    for (i in seq_along(ligando_indices)) {
      ligando <- ligando_names[i]
      start <- ligando_indices[i]
      end <- block_ends[i]
      block <- lines[start:end]
      
      # Encontrar líneas con datos numéricos
      data_lines <- block[grepl("^\\s*\\d+\\s+\\d+\\s+\\d+\\s+-?\\d+\\.\\d+\\s+\\d+\\.\\d+\\s+\\d+\\.\\d+", block)]
      
      if (length(data_lines) > 0) {
        parsed_data <- str_match_all(
          data_lines,
          "^\\s*(\\d+)\\s+(\\d+)\\s+(\\d+)\\s+(-?\\d+\\.\\d+)\\s+(\\d+\\.\\d+)\\s+(\\d+\\.\\d+).*$"
        )
        
        df_blocks <- lapply(parsed_data, function(match) {
          if (nrow(match) > 0) {
            tibble(
              Ligando = ligando,
              Rank = as.integer(match[1, 2]),
              SubRank = as.integer(match[1, 3]),
              Run = as.integer(match[1, 4]),
              Binding_Energy = as.numeric(match[1, 5]),
              Cluster_RMSD = as.numeric(match[1, 6]),
              Reference_RMSD = as.numeric(match[1, 7])
            )
          }
        })
        
        df_block <- bind_rows(df_blocks)
        rows[[length(rows) + 1]] <- df_block
      }
    }
    
    if (length(rows) > 0) {
      df <- bind_rows(rows)
      write_csv(df, out_csv)
      cat("✓ Procesado:", txt_file, "->", out_csv, "(", nrow(df), "filas )\n")
      return(TRUE)
    } else {
      cat("✗ No se pudieron extraer datos de:", txt_file, "\n")
      return(FALSE)
    }
    
  }, error = function(e) {
    cat("✗ Error procesando", txt_file, ":", e$message, "\n")
    return(FALSE)
  })
}

# ---- Ejecución desde línea de comandos ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Uso: Rscript parse_docking.R <archivo_entrada.txt> <archivo_salida.csv>")
}

input_file <- args[1]
output_file <- args[2]

parse_docking_file(input_file, output_file)