#!/bin/bash
# process_good_ligands.sh

cd "/home/db/Documentos/David_Hito2/lig_pdb"
mkdir -p "../pdbqt_files"
mkdir -p "../bad_pdbs"

PYTHONSH="/home/db/Documentos/MGLTools-1.5.7/bin/pythonsh"
PREPARE_SCRIPT="/home/db/Documentos/David_Hito2/prepare_ligand4.py"

echo "=== FILTRANDO LIGANDOS PARA AUTODOCK4 ==="
echo "Buscando archivos problemáticos..."

# Contadores
total=$(ls *.pdb | wc -l)
bad_coords=0
bad_se=0
good=0

for file in *.pdb; do
    # Verificar coordenadas cero
    if grep -q "0.000   0.000   0.000" "$file"; then
        echo "❌ $file - COORDENADAS CERO"
        mv "$file" "../bad_pdbs/zero_coords_$file"
        ((bad_coords++))
    
    # Verificar selenio
    elif grep -q " SE " "$file" || grep -q "SE   " "$file"; then
        echo "❌ $file - CONTIENE SELENIO"
        mv "$file" "../bad_pdbs/se_$file" 
        ((bad_se++))
    
    # Archivo bueno
    else
        echo "✅ $file - PROCESANDO"
        if "$PYTHONSH" "$PREPARE_SCRIPT" -l "$file" -o "../pdbqt_files/${file%.pdb}.pdbqt" 2>/dev/null; then
            ((good++))
        else
            echo "   ⚠️  Error en conversión, moviendo a bad_pdbs"
            mv "$file" "../bad_pdbs/error_$file"
        fi
    fi
done

echo "=========================================="
echo "RESUMEN FINAL:"
echo " - Total archivos: $total"
echo " - ✅ Buenos y convertidos: $good"
echo " - ❌ Coordenadas cero: $bad_coords" 
echo " - ❌ Con selenio: $bad_se"
echo " - Archivos buenos en: ../pdbqt_files/"
echo " - Archivos malos en: ../bad_pdbs/"
