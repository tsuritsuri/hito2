#!/bin/bash
# unify_by_protein.sh - Unificar resultados desde las carpetas de frames

BASE_DIR="/home/db/Documentos/David_Hito2/autodock_out"
OUTPUT_DIR="docking_by_protein"

echo "Unificando resultados por proteína desde carpetas de frames..."
mkdir -p "$OUTPUT_DIR"

# Proteínas a procesar
PROTEINAS=("1E9Y" "1E9Z" "6QSU" "6ZJA")

for proteina in "${PROTEINAS[@]}"; do
    echo "Procesando $proteina..."

    OUTPUT_FILE="$OUTPUT_DIR/docking_results_${proteina}.csv"

    # Cabecera con columna adicional para Frame
    echo "Ligando,Rank,SubRank,Run,Binding_Energy,Cluster_RMSD,Reference_RMSD,Frame,Proteina" > "$OUTPUT_FILE"

    # Buscar todos los CSV en los frames de esta proteína
    find "$BASE_DIR/$proteina" -name "docking_results_parsed.csv" | while read csv_file; do
        # Extraer nombre del frame desde la ruta: autodock_out/1E9Y/frame100/dlg/...
        frame=$(basename $(dirname $(dirname "$csv_file")))

        # Procesar el CSV y añadir columnas de Frame y Proteína
        awk -F',' -v frame="$frame" -v proteina="$proteina" '
        NR == 1 { next }  # Saltar cabecera del CSV individual
        {
            # Reconstruir línea con las nuevas columnas
            print $1 "," $2 "," $3 "," $4 "," $5 "," $6 "," $7 "," frame "," proteina
        }
        ' "$csv_file"
    done >> "$OUTPUT_FILE"

    # Estadísticas
    total_lines=$(tail -n +2 "$OUTPUT_FILE" | wc -l)
    unique_ligands=$(tail -n +2 "$OUTPUT_FILE" | cut -d',' -f1 | sort -u | wc -l)
    frames_used=$(tail -n +2 "$OUTPUT_FILE" | cut -d',' -f8 | sort -u | wc -l)

    echo "  → $total_lines resultados totales"
    echo "  → $unique_ligands ligandos únicos"
    echo "  → $frames_used frames utilizados"
done

echo "=== Resumen final ==="
for proteina in "${PROTEINAS[@]}"; do
    file="$OUTPUT_DIR/docking_results_${proteina}.csv"
    if [[ -f "$file" ]]; then
        count=$(tail -n +2 "$file" | wc -l)
        echo "$proteina: $count resultados"
    else
        echo "$proteina: NO ENCONTRADO"
    fi
done

echo "Archivos guardados en: $OUTPUT_DIR/"
