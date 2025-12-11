#!/usr/bin/env bash
# Vina batch docking (AD Vina, robust)

set -euo pipefail

# --------------- CONFIG ----------------
BASE="/home/db/Documentos/David_Hito2"
PROTEIN_ROOT="$BASE/frames_pdbqt"      # subdirs: 1E9Y, 1E9Z, 6QSU, 6ZJA
LIGAND_DIR="$BASE/ligandos/lig_pdbqt"           # ligandos .pdbqt
OUT_ROOT="$BASE/vina_out"
LIGAND_EXT="pdbqt"

VINA_BIN="vina"                        # o ruta absoluta si no está en PATH
MAX_JOBS=40
SEED=42
EXHAUSTIVENESS=16
NUM_MODES=10
NUM_CPUS=1

# Centros por proteína (mismo para todos los frames)
declare -A CX CY CZ
CX["1E9Y"]=-30.429; CY["1E9Y"]=8.276;   CZ["1E9Y"]=-15.299
CX["1E9Z"]=-30.429; CY["1E9Z"]=8.276;   CZ["1E9Z"]=-15.299
CX["6QSU"]=4.233;   CY["6QSU"]=-7.758;  CZ["6QSU"]=33.454
CX["6ZJA"]=-12.675; CY["6ZJA"]=2.876;   CZ["6ZJA"]=-32.740

# Tamaños de caja (constantes aquí; cámbialos si necesitas por proteína)
SIZE_X=15
SIZE_Y=15
SIZE_Z=15
# ---------------------------------------

# ---- helpers ----
trap 'echo "ERROR on line $LINENO"; exit 1' ERR

command -v "$VINA_BIN" >/dev/null 2>&1 || {
  echo "ERROR: no encuentro '$VINA_BIN' en PATH"; exit 1; }

mkdir -p "$OUT_ROOT"

# Ligandos (lista estable)
mapfile -d '' -t ligands < <(find "$LIGAND_DIR" -maxdepth 1 -type f -name "*.${LIGAND_EXT}" -print0)
(( ${#ligands[@]} > 0 )) || { echo "No hay ligandos *.$LIGAND_EXT"; exit 1; }
TOTAL_LIGANDS=${#ligands[@]}

# Cuenta segura de archivos con patrón (sin usar ls)
count_matches() {
  local dir="$1" ext="$2"
  local -a arr=()
  shopt -s nullglob
  arr=( "$dir"/*."$ext" )
  shopt -u nullglob
  echo "${#arr[@]}"
}

wait_for_slot() {
  while (( $(jobs -rp | wc -l) >= MAX_JOBS )); do
    # Espera al menos a que termine uno
    wait -n 2>/dev/null || true
  done
}

echo "===> Vina batch"
echo "Proteínas: $PROTEIN_ROOT"
echo "Ligandos:  $LIGAND_DIR  (TOTAL=$TOTAL_LIGANDS)"
echo "Salida:    $OUT_ROOT"
echo

# Recorre estructuras (subdirectorios de PROTEIN_ROOT)
while IFS= read -r -d '' struct_dir; do
  struct_name="$(basename "$struct_dir")"
  echo "---- Proteína: $struct_name ----"

  # centro definido?
  if [[ -z "${CX[$struct_name]+_}" ]]; then
    echo "  ERROR: no definiste centro para $struct_name"; continue
  fi

  # frames .pdbqt del receptor
  shopt -s nullglob
  frames=( "$struct_dir"/*.pdbqt )
  shopt -u nullglob
  (( ${#frames[@]} > 0 )) || { echo "  (sin frames .pdbqt)"; continue; }

  STRUCT_OUT="$OUT_ROOT/$struct_name"; mkdir -p "$STRUCT_OUT"

  for REC_FILE in "${frames[@]}"; do
    frame_base="$(basename "$REC_FILE" .pdbqt)"
    FRAME_OUT="$STRUCT_OUT/$frame_base"
    mkdir -p "$FRAME_OUT"

    have_count="$(count_matches "$FRAME_OUT" "$LIGAND_EXT")"
    if (( have_count >= TOTAL_LIGANDS )); then
      echo "  > Frame $frame_base [skip] ya completo ($have_count/$TOTAL_LIGANDS)"
      continue
    fi

    echo "  > Frame: $frame_base"
    for LIG in "${ligands[@]}"; do
      LIG_NAME="$(basename "$LIG" ."$LIGAND_EXT")"
      OUT_POSE="$FRAME_OUT/${LIG_NAME}.pdbqt"
      OUT_LOG="$FRAME_OUT/${LIG_NAME}.log"
      [[ -f "$OUT_POSE" ]] && { echo "    [skip] $LIG_NAME"; continue; }

      (
        "$VINA_BIN" \
          --receptor "$REC_FILE" \
          --ligand   "$LIG" \
          --center_x "${CX[$struct_name]}" \
          --center_y "${CY[$struct_name]}" \
          --center_z "${CZ[$struct_name]}" \
          --size_x   "$SIZE_X" \
          --size_y   "$SIZE_Y" \
          --size_z   "$SIZE_Z" \
          --exhaustiveness "$EXHAUSTIVENESS" \
          --num_modes "$NUM_MODES" \
          --seed "$SEED" \
          --cpu "$NUM_CPUS" \
          --out "$OUT_POSE" \
          --log "$OUT_LOG"
      ) &

      wait_for_slot
    done

    echo "    Esperando jobs del frame $frame_base..."
    while (( $(jobs -rp | wc -l) > 0 )); do
      wait -n 2>/dev/null || true
    done

    have_count="$(count_matches "$FRAME_OUT" "$LIGAND_EXT")"
    echo "  > Frame $frame_base OK: $have_count / $TOTAL_LIGANDS"
  done

  echo "---- Fin: $struct_name ----"
  echo
done < <(find "$PROTEIN_ROOT" -mindepth 1 -maxdepth 1 -type d -print0)

echo "===> Resultados: $OUT_ROOT"
