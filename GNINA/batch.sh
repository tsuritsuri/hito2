#!/usr/bin/env bash
# GNINA batch docking (Dual-GPU Optimized)

set -euo pipefail

# --------------- CONFIG ----------------
BASE="/home/db/Documentos/David_Hito2"
PROTEIN_ROOT="$BASE/frames_pdbqt"
LIGAND_DIR="$BASE/ligandos/lig_pdbqt"
OUT_ROOT="$BASE/gnina_out"
LIGAND_EXT="pdbqt"
OUT_EXT="sdf"

VINA_BIN="gnina"
MAX_JOBS=8               # <--- LIMIT THIS to prevent GPU crash (approx 8 per GPU)
NUM_CPUS=5                # <--- INCREASE THIS: 16 jobs * 3 cores = 48 cores total usage
SEED=42
EXHAUSTIVENESS=16         # Balanced for speed/accuracy

# Coordinates (Same as before)
declare -A CX CY CZ
CX["1E9Y"]=-30.429; CY["1E9Y"]=8.276;   CZ["1E9Y"]=-15.299
CX["1E9Z"]=-30.429; CY["1E9Z"]=8.276;   CZ["1E9Z"]=-15.299
CX["6QSU"]=4.233;   CY["6QSU"]=-7.758;  CZ["6QSU"]=33.454
CX["6ZJA"]=-12.675; CY["6ZJA"]=2.876;   CZ["6ZJA"]=-32.740

SIZE_X=15; SIZE_Y=15; SIZE_Z=15
# ---------------------------------------

command -v "$VINA_BIN" >/dev/null 2>&1 || { echo "ERROR: gnina not found"; exit 1; }
mkdir -p "$OUT_ROOT"

# Load Ligands
mapfile -d '' -t ligands < <(find "$LIGAND_DIR" -maxdepth 1 -type f -name "*.${LIGAND_EXT}" -print0)
TOTAL_LIGANDS=${#ligands[@]}
(( TOTAL_LIGANDS > 0 )) || { echo "No ligands found."; exit 1; }

wait_for_slot() {
  while (( $(jobs -rp | wc -l) >= MAX_JOBS )); do
    wait -n 2>/dev/null || true
  done
}

job_counter=0

echo "===> GNINA Dual-GPU Batch"
echo "     CPUs: $NUM_CPUS per job | Max Jobs: $MAX_JOBS"
echo

# Recorre estructuras
while IFS= read -r -d '' struct_dir; do
  struct_name="$(basename "$struct_dir")"

  if [[ -z "${CX[$struct_name]+_}" ]]; then continue; fi

  shopt -s nullglob
  frames=( "$struct_dir"/*.pdbqt )
  shopt -u nullglob

  STRUCT_OUT="$OUT_ROOT/$struct_name"; mkdir -p "$STRUCT_OUT"

  for REC_FILE in "${frames[@]}"; do
    frame_base="$(basename "$REC_FILE" .pdbqt)"
    FRAME_OUT="$STRUCT_OUT/$frame_base"
    mkdir -p "$FRAME_OUT"

    # Skip if done
    if (( $(find "$FRAME_OUT" -name "*.$OUT_EXT" | wc -l) >= TOTAL_LIGANDS )); then
      echo "  > Frame $frame_base [skip] complete."
      continue
    fi

    echo "  > Processing Frame: $frame_base"

    for LIG in "${ligands[@]}"; do
      LIG_NAME="$(basename "$LIG" ."$LIGAND_EXT")"
      OUT_POSE="$FRAME_OUT/${LIG_NAME}.${OUT_EXT}"
      OUT_LOG="$FRAME_OUT/${LIG_NAME}.log"

      [[ -f "$OUT_POSE" ]] && continue

      # --- GPU SELECTION FIX ---
      # Use safe math to prevent 'set -e' crash
      GPU_ID=$(( job_counter % 2 ))

      # Safe increment (avoiding the ((...)) shorthand which can return error on 0)
      job_counter=$((job_counter + 1))

      (
        # Export variable ONLY for this sub-process
        export CUDA_VISIBLE_DEVICES=$GPU_ID

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
          --num_modes 10 \
          --seed "$SEED" \
          --cpu "$NUM_CPUS" \
          --out "$OUT_POSE" \
          --log "$OUT_LOG" \
          --cnn_scoring "rescore" \
          --cnn "crossdock_default2018" >/dev/null 2>&1
      ) &

      wait_for_slot
    done

    # Wait for current frame to finish before moving to next
    wait
    echo "  > Frame $frame_base Done."
  done
done < <(find "$PROTEIN_ROOT" -mindepth 1 -maxdepth 1 -type d -print0)
