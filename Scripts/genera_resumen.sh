#!/bin/bash

# Directorio madre
BASE_DIR="/home/db/Documentos/David_Hito2"

echo "Iniciando procesamiento de carpetas dlg..."
echo "==========================================="
echo ""

# Contador de carpetas procesadas
count=0

# Buscar todas las carpetas llamadas "dlg" recursivamente
find "$BASE_DIR" -type d -name "dlg" | while read dlg_dir; do
    # Contar archivos .dlg en esta carpeta
    num_dlg=$(find "$dlg_dir" -maxdepth 1 -name "*.dlg" -type f | wc -l)
    
    if [ "$num_dlg" -eq 0 ]; then
        echo "⊘ Saltando $dlg_dir (sin archivos .dlg)"
        continue
    fi
    
    ((count++))
    
    # Obtener ruta relativa para info
    rel_path=$(echo "$dlg_dir" | sed "s|$BASE_DIR/||")
    
    echo "📁 Procesando: $rel_path"
    echo "   Archivos .dlg encontrados: $num_dlg"
    
    # Archivos de salida en la misma carpeta dlg
    SUMMARY_FILE="${dlg_dir}/resumen_resultados.txt"
    CSV_FILE="${dlg_dir}/resumen_resultados.csv"
    
    > "$SUMMARY_FILE"  # Vaciar o crear el archivo resumen TXT
    
    # Crear encabezado del CSV
    echo "Ligando,Rank,Binding_Energy,RMSD_Reference,RMSD_LowerBound,Cluster_Size" > "$CSV_FILE"
    
    # Procesar todos los archivos .dlg en esta carpeta dlg
    for dlg_file in "$dlg_dir"/*.dlg; do
        if [ -f "$dlg_file" ]; then
            # Obtener el nombre del ligando sin la extensión
            ligand_name=$(basename "$dlg_file" .dlg)
            
            # Agregar al TXT
            echo "Procesando ligando: $ligand_name" >> "$SUMMARY_FILE"
            echo "-------------------------------------" >> "$SUMMARY_FILE"
            
            # Extraer la tabla RMSD y las 22 líneas siguientes
            awk '/RMSD TABLE/ {print; c=22; next} c-- > 0 {print}' "$dlg_file" >> "$SUMMARY_FILE"
            
            # Extraer datos para el CSV
            awk -v ligand="$ligand_name" '/RMSD TABLE/,/^$/ {
                if ($1 ~ /^[0-9]+$/ && NF >= 5) {
                    print ligand "," $1 "," $4 "," $5 "," $6 "," $2
                }
            }' "$dlg_file" >> "$CSV_FILE"
            
            echo -e "\n" >> "$SUMMARY_FILE"
        fi
    done
    
    echo "   ✓ Generados: resumen_resultados.txt y resumen_resultados.csv"
    echo ""
done

echo "==========================================="
echo "Proceso completado."
echo "Se procesaron las carpetas dlg encontradas."
echo "==========================================="
