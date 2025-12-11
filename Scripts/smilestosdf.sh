#!/bin/bash
set -euo pipefail

CSV_FILE="datos.csv"
OUTPUT_DIR="molecules_output"
FORMAT="sdf"   # "sdf" o "mol2"
mkdir -p "$OUTPUT_DIR"

clean_filename(){ echo "$1" | sed 's/[^a-zA-Z0-9._-]/_/g; s/__*/_/g'; }
unquote(){ local s="${1:-}"; s="${s#\"}"; s="${s%\"}"; printf '%s' "$s"; }

have_coords(){
  # devuelve 0 si hay coords (no todo 0.0000), 1 si todo es 0.0
  local f="$1"
  local n_atoms
  n_atoms=$(awk 'NR==4{print $1}' "$f")
  # cuenta líneas con "0.0000    0.0000    0.0000" en los bloques de átomos (sólo aprox)
  local zeros
  zeros=$(grep -cE '^[[:space:]]*0\.0000[[:space:]]+0\.0000[[:space:]]+0\.0000[[:space:]]+[A-Z]' "$f" || true)
  # si el número de “cero-cero-cero” >= n_atoms → no hay coords
  [ -n "$n_atoms" ] && [ "$zeros" -ge "$n_atoms" ] && return 1 || return 0
}

echo "Iniciando conversión SMILES → $FORMAT"
LOG_FILE="conversion_log.txt"
: > "$LOG_FILE"

line_number=0 success=0 fail=0 fallback2d=0

# AVISO: parser CSV simple (sin comas internas). Si tu CSV tiene comas entre comillas, dímelo y te paso versión en Python.
tail -n +2 "$CSV_FILE" | \
while IFS=',' read -r compound_id smiles source_db organism assay_type activity_type value_nm units reference notes inchikey; do
  line_number=$((line_number+1))

  compound_id=$(unquote "$compound_id"); smiles=$(unquote "$smiles")
  source_db=$(unquote "$source_db"); organism=$(unquote "$organism")
  assay_type=$(unquote "$assay_type"); activity_type=$(unquote "$activity_type")
  value_nm=$(unquote "$value_nm"); units=$(unquote "$units")
  reference=$(unquote "$reference"); notes=$(unquote "$notes")
  inchikey=$(unquote "$inchikey")

  if [ -z "$smiles" ] || [ "$smiles" = "NULL" ] || [ "$smiles" = "N/A" ]; then
    echo "L$line_number: SMILES vacío, salto" | tee -a "$LOG_FILE"
    fail=$((fail+1)); continue
  fi

  if [ -n "$compound_id" ] && [ "$compound_id" != "NULL" ]; then
    out="${OUTPUT_DIR}/${line_number}_$(clean_filename "$compound_id").${FORMAT}"
  else
    out="${OUTPUT_DIR}/compound_${line_number}.${FORMAT}"
  fi

  echo "L$line_number: $compound_id"

  # 1) Intento 3D
  err3d=$(mktemp)
  if obabel -:"$smiles" -o"$FORMAT" --gen3D \
      --title "$compound_id" \
      --property "Compound_ID" "$compound_id" \
      --property "Source_DB" "$source_db" \
      --property "Organism" "$organism" \
      --property "Assay_Type" "$assay_type" \
      --property "Activity_Type" "$activity_type" \
      --property "Value_nM" "$value_nm" \
      --property "Units" "$units" \
      --property "Reference" "$reference" \
      --property "Notes" "$notes" \
      --property "InChIKey" "$inchikey" \
      --property "Line_Number" "$line_number" \
      -O "$out" 2> "$err3d"; then

    if grep -q "3D coordinate generation failed" "$err3d" || ! have_coords "$out"; then
      # 2) Fallback a 2D
      err2d=$(mktemp)
      if obabel -:"$smiles" -o"$FORMAT" --gen2D -O "$out" 2> "$err2d"; then
        if [ "$FORMAT" = "sdf" ]; then
          tmp=$(mktemp); head -n -1 "$out" > "$tmp"
          { echo; echo "> <IC50>"; echo "${value_nm:-}"; echo; printf "\$\$\$\$\n"; } >> "$tmp"
          mv "$tmp" "$out"
        fi
        echo "  ↳ 3D FALLÓ, guardado en 2D → $out" | tee -a "$LOG_FILE"
        echo "  (stderr 3D):"; sed 's/^/    /' "$err3d" >> "$LOG_FILE"
        echo "  (stderr 2D):"; sed 's/^/    /' "$err2d" >> "$LOG_FILE"
        fallback2d=$((fallback2d+1))
      else
        echo "  ✗ Error: ni 3D ni 2D generaron coords → $compound_id" | tee -a "$LOG_FILE"
        sed 's/^/    /' "$err3d" >> "$LOG_FILE"
        sed 's/^/    /' "$err2d" >> "$LOG_FILE"
        rm -f "$out"
        fail=$((fail+1))
      fi
      rm -f "$err2d"
    else
      # 3D OK → añade IC50 si SDF
      if [ "$FORMAT" = "sdf" ]; then
        tmp=$(mktemp); head -n -1 "$out" > "$tmp"
        { echo; echo "> <IC50>"; echo "${value_nm:-}"; echo; printf "\$\$\$\$\n"; } >> "$tmp"
        mv "$tmp" "$out"
      fi
      echo "  ✓ 3D OK → $out" | tee -a "$LOG_FILE"
      success=$((success+1))
    fi
  else
    echo "  ✗ Error de obabel al convertir SMILES" | tee -a "$LOG_FILE"
    sed 's/^/    /' "$err3d" >> "$LOG_FILE"
    fail=$((fail+1))
  fi
  rm -f "$err3d"
done

echo "----------------------------------------"
echo "Hecho. 3D OK: $success | 2D fallback: $fallback2d | Fallos: $fail"
echo "Log: $LOG_FILE"


