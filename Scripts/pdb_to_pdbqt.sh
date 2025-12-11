#!/bin/bash
# pdb_to_pdbqt_fixed.sh

# RUTAS ABSOLUTAS
PYTHONSH="/home/db/Documentos/MGLTools-1.5.7/bin/pythonsh"
PREPARE_SCRIPT="/home/db/Documentos/David_Hito2/prepare_ligand4.py"
INPUT_DIR="/home/db/Documentos/David_Hito2/lig_pdb"
OUTPUT_DIR="/home/db/Documentos/David_Hito2/pdbqt_files"

cd "$INPUT_DIR"
mkdir -p "$OUTPUT_DIR"

count=0
for file in *.pdb; do
    ((count++))
    echo "[$count] Converting: $file"
    "$PYTHONSH" "$PREPARE_SCRIPT" -l "$file" -o "$OUTPUT_DIR/${file%.pdb}.pdbqt"
done

echo "Done: $count files"
