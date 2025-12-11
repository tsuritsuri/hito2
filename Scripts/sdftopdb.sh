#!/bin/bash

# Definir carpetas
INPUT_DIR="lig_sdf"
OUTPUT_DIR="lig_pdb"

# Crear carpeta de salida si no existe
mkdir -p "$OUTPUT_DIR"

# Recorrer y convertir todos los .sdf dentro de lig_sdf
for f in "$INPUT_DIR"/*.sdf; do
    fname=$(basename "$f" .sdf)
    echo "Convirtiendo $fname.sdf..."
    obabel "$f" -O "$OUTPUT_DIR/$fname.pdb"
done

echo "✅ Conversión completada. Archivos .pdb guardados en $OUTPUT_DIR"

