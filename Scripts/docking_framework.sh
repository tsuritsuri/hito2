#!/bin/bash
# AD4 batch docking (AutoGrid4 + AutoDock4)

set -o pipefail

# ------------------------- CONFIG -------------------------
BASE="/home/db/Documentos/David_Hito2"
PROTEIN_ROOT="$BASE/frames_pdbqt"           # 1E9Y, 1E9Z, 6QSU, 6ZJA
LIGAND_DIR="$BASE/pdbqt_files"       # ligandos .pdbqt
DPF_TEMPLATE="$BASE/dpf_template.dpf"
OUT_ROOT="$BASE/autodock_out2"
LIGAND_EXT="pdbqt"

MODE="medium"                               # fast | medium | heavy
MAX_JOBS=35
SCALE_BY_TORSDOF=1                          # 1=on, 0=off

LIGAND_TYPES=(A Br C Cl F HD N NA OA P S SA)
LIGAND_TYPES_STR="$(printf "%s " "${LIGAND_TYPES[@]}")"; LIGAND_TYPES_STR="${LIGAND_TYPES_STR% }"

PARAM_NI="$BASE/AD4.1_bound_Ni.dat"

AUTOGRID_BIN="autogrid4"
AUTODOCK_BIN="autodock4"
# --------------------- FIN CONFIG -------------------------

# -------- MASTER LOG --------
MASTER_LOG="$BASE/docking_master.log"
{
  echo "==== AD4 batch run $(date '+%F %T') ===="
  echo "HOST: $(hostname) | MODE=$MODE | MAX_JOBS=$MAX_JOBS | SCALE_BY_TORSDOF=$SCALE_BY_TORSDOF"
  echo "PROTEIN_ROOT=$PROTEIN_ROOT"
  echo "LIGAND_DIR=$LIGAND_DIR"
  echo "OUT_ROOT=$OUT_ROOT"
  echo
} >> "$MASTER_LOG"
exec > >(tee -a "$MASTER_LOG") 2>&1
# ----------------------------

# ---------- checks ----------
for d in "$PROTEIN_ROOT" "$LIGAND_DIR"; do
  [[ -d "$d" ]] || { echo "ERROR: No existe directorio: $d"; exit 1; }
done
[[ -f "$DPF_TEMPLATE" ]] || { echo "ERROR: Falta DPF_TEMPLATE: $DPF_TEMPLATE"; exit 1; }

if [[ ! -f "$PARAM_NI" ]]; then
  echo "ADVERTENCIA: buscando AD4.1_bound_Ni.dat en $BASE ..."
  cand=$(find "$BASE" -maxdepth 3 -type f -name 'AD4.1_bound_Ni.dat' 2>/dev/null | head -n1)
  [[ -n "$cand" ]] || { echo "ERROR: no se encontró AD4.1_bound_Ni.dat"; exit 1; }
  PARAM_NI="$cand"
  echo "Usando PARAM_NI: $PARAM_NI"
fi

for exe in "$AUTOGRID_BIN" "$AUTODOCK_BIN"; do
  command -v "$exe" >/dev/null 2>&1 || echo "ADVERTENCIA: '$exe' no está en PATH."
done

# Ligandos
mapfile -d '' -t ligands < <(find "$LIGAND_DIR" -maxdepth 1 -type f -name "*.${LIGAND_EXT}" -print0)
[[ ${#ligands[@]} -gt 0 ]] || { echo "ERROR: No hay ligandos *.$LIGAND_EXT en $LIGAND_DIR"; exit 1; }
TOTAL_LIGANDS=${#ligands[@]}
echo "Ligandos detectados: $TOTAL_LIGANDS"

mkdir -p "$OUT_ROOT"

# ---------- helpers ----------
compute_ga() {
  local tors=${1:-6}
  local base_evals base_runs pop swits
  case "$MODE" in
    fast)   base_evals=250000;  base_runs=10; pop=100; swits=150;;
    medium) base_evals=1000000; base_runs=10; pop=150; swits=300;;
    heavy)  base_evals=2000000; base_runs=10; pop=150; swits=300;;
    *)      base_evals=1000000; base_runs=10; pop=150; swits=300;;
  esac
  local mult=1
  if (( SCALE_BY_TORSDOF )); then
    if   (( tors <= 10 )); then mult=1
    else                         mult=2
    fi
  fi
  GA_EVALS=$(( base_evals * mult ))
  GA_RUNS=$base_runs
  GA_POP=$pop
  SW_ITS=$swits
}

wait_for_slot() {
  while (( $(jobs -rp | wc -l) >= MAX_JOBS )); do
    wait -n 2>/dev/null || true
  done
}

# ---------- MAPS con gridfld ABSOLUTO ----------
# crea/valida: <frame>.{A,Br,C,Cl,F,HD,N,NA,OA,P,S,SA}.map + <frame>.e.map + <frame>.d.map + <frame>.maps.fld
ensure_maps() {
  local maps_dir="$1" gpf_file="$2" frame_base="$3" protein_file="$4" param_ni="$5"
  mkdir -p "$maps_dir"
  local fld_target="${frame_base}.maps.fld"
  local fld_abs="${maps_dir}/${fld_target}"
  local glg_log="${frame_base}.glg"

  local missing=0
  for t in "${LIGAND_TYPES[@]}"; do
    [[ -f "$maps_dir/${frame_base}.${t}.map" ]] || missing=1
  done
  [[ -f "$maps_dir/${frame_base}.e.map" ]] || missing=1
  [[ -f "$maps_dir/${frame_base}.d.map" ]] || missing=1
  [[ -f "$fld_abs" ]] || missing=1

  if (( missing )); then
    echo "    (Re)generando MAPS en $maps_dir ..."
    pushd "$maps_dir" >/dev/null
    local GPF_TMP="$maps_dir/${frame_base}.tmp.gpf"
    sed -E \
      -e "s|^parameter_file\b.*|parameter_file $param_ni|g" \
      -e "s|^receptor\b.*|receptor $protein_file|g" \
      -e "s|^gridfld\b.*|gridfld ${fld_abs}|g" \
      -e "s|^ligand_types\b.*|ligand_types $LIGAND_TYPES_STR|g" \
      -e "/^map\b.*/d" \
      -e "/^elecmap\b.*/d" \
      -e "/^(desolvmap|desolvmap)\b.*/d" \
      "$gpf_file" > "$GPF_TMP"

    if ! grep -qE '^ligand_types\b' "$GPF_TMP"; then
      sed -i "1iligand_types $LIGAND_TYPES_STR" "$GPF_TMP"
    fi
    {
      for t in "${LIGAND_TYPES[@]}"; do
        echo "map ${frame_base}.${t}.map"
      done
      echo "elecmap ${frame_base}.e.map"
      echo "desolvmap ${frame_base}.d.map"
    } >> "$GPF_TMP"

    "$AUTOGRID_BIN" -p "$GPF_TMP" -l "$glg_log"
    local ag_status=$?
    rm -f "$GPF_TMP"
    popd >/dev/null

    if [[ $ag_status -ne 0 || ! -f "$fld_abs" ]]; then
      echo "    ERROR: AutoGrid falló o falta $fld_abs (ver $maps_dir/$glg_log)"
      return 1
    fi
    for t in "${LIGAND_TYPES[@]}"; do
      [[ -f "$maps_dir/${frame_base}.${t}.map" ]] || { echo "    ERROR: falta ${frame_base}.${t}.map"; return 1; }
    done
    [[ -f "$maps_dir/${frame_base}.e.map" && -f "$maps_dir/${frame_base}.d.map" ]] || { echo "    ERROR: faltan e/d maps"; return 1; }
  fi

  # sanity: si por alguna razón quedaron relativos, normalizamos headers de .map
  sed -i -E "s|^gridfld\b.*|gridfld ${fld_abs}|" "$maps_dir"/*.map 2>/dev/null || true

  FLD_PATH="$fld_abs"
  return 0
}

echo "===> Iniciando batch docking"
echo "Proteínas: $PROTEIN_ROOT"
echo "Ligandos:  $LIGAND_DIR  (TOTAL=$TOTAL_LIGANDS)"
echo "Salida:    $OUT_ROOT"
echo "Modo:      $MODE  | MAX_JOBS=$MAX_JOBS | scale_by_torsdof=$SCALE_BY_TORSDOF"
echo

# ---------- main loop ----------
while IFS= read -r -d '' struct_dir; do
  struct_name="$(basename "$struct_dir")"
  echo "---- Estructura: $struct_name ----"

  shopt -s nullglob; frames=( "$struct_dir"/*.pdbqt ); shopt -u nullglob
  [[ ${#frames[@]} -gt 0 ]] || { echo "  (sin frames .pdbqt)"; continue; }

  STRUCT_OUT="$OUT_ROOT/$struct_name"; mkdir -p "$STRUCT_OUT"
  FALLBACK_GPF="$struct_dir/grid.gpf"

  for PROTEIN_FILE in "${frames[@]}"; do
    frame_base="$(basename "$PROTEIN_FILE" .pdbqt)"
    echo "  > Frame: $frame_base"

    FRAME_OUT="$STRUCT_OUT/${frame_base}/dlg"
    mkdir -p "$FRAME_OUT"
    have_count=$(ls "$FRAME_OUT"/*.dlg 2>/dev/null | wc -l)
    if [[ "$have_count" -ge "$TOTAL_LIGANDS" ]]; then
      echo "    [skip] frame $frame_base ya completo ($have_count/$TOTAL_LIGANDS)"
      continue
    fi

    FRAME_GPF="$struct_dir/$frame_base.gpf"
    if   [[ -f "$FRAME_GPF"    ]]; then GPF_FILE="$FRAME_GPF"
    elif [[ -f "$FALLBACK_GPF" ]]; then GPF_FILE="$FALLBACK_GPF"
    else echo "    ERROR: no hay .gpf (ni $FRAME_GPF ni $FALLBACK_GPF)"; continue
    fi

    MAPS_DIR="$STRUCT_OUT/${frame_base}/maps"
    mkdir -p "$MAPS_DIR"

    if ensure_maps "$MAPS_DIR" "$GPF_FILE" "$frame_base" "$PROTEIN_FILE" "$PARAM_NI"; then
      echo "    MAPS OK: $MAPS_DIR"
    else
      echo "    ERROR: MAPS no listos para $frame_base"; continue
    fi

    # ---- loop ligandos con paralelismo ----
    for LIGAND_FILE in "${ligands[@]}"; do
      LIGAND_NAME="$(basename "$LIGAND_FILE" ."$LIGAND_EXT")"
      DPF_FILE="$FRAME_OUT/${LIGAND_NAME}.dpf"
      DLG_FILE="$FRAME_OUT/${LIGAND_NAME}.dlg"

      [[ -f "$DLG_FILE" ]] && { echo "    [skip] $LIGAND_NAME (ya existe)"; continue; }

      cp "$DPF_TEMPLATE" "$DPF_FILE"
      sed -i 's/\r$//' "$DPF_FILE"

      # center desde .fld
      if [[ -f "$MAPS_DIR/${frame_base}.maps.fld" ]]; then
        grid_center=$(grep "^#CENTER" "$MAPS_DIR/${frame_base}.maps.fld" | head -1 | awk '{print $2, $3, $4}')
        if [[ -n "$grid_center" ]]; then
          if grep -q "^about " "$DPF_FILE"; then
            sed -i -E "s|^about .*|about $grid_center|g" "$DPF_FILE"
          else
            sed -i -E "1i\about $grid_center" "$DPF_FILE"
          fi

          sed -i -E '/^tran0 /d; /^quaternion0 /d; /^quat0 /d' "$DPF_FILE"
          {
            echo "tran0 $grid_center"
            echo "quat0 0.0 0.0 0.0 1.0"
          } >> "$DPF_FILE"

          echo "    ✓ Grid center: $grid_center"
        else
          echo "    ✗ ERROR: No se pudo extraer #CENTER del .fld"; continue
        fi
      else
        echo "    ✗ ERROR: No existe $MAPS_DIR/${frame_base}.maps.fld"; continue
      fi

      # ligand_types
      if grep -qE '^ligand_types\b' "$DPF_FILE"; then
        sed -i -E "s|^ligand_types\b.*|ligand_types $LIGAND_TYPES_STR|g" "$DPF_FILE"
      else
        sed -i -E "1i\ligand_types $LIGAND_TYPES_STR" "$DPF_FILE"
      fi

      # fld y maps (maps insertados tras fld)
      if grep -qE '^fld\b' "$DPF_FILE"; then
        sed -i -E "s|^fld\b.*|fld $MAPS_DIR/${frame_base}.maps.fld|g" "$DPF_FILE"
      else
        sed -i -E "/^ligand_types\b.*/a fld $MAPS_DIR/${frame_base}.maps.fld" "$DPF_FILE"
      fi
      sed -i -E "/^(map|elecmap|desolvmap|desolvmap)\b.*/d" "$DPF_FILE"
      awk -v md="$MAPS_DIR" -v fb="$frame_base" -v types="$LIGAND_TYPES_STR" '
        BEGIN{n=split(types,T," ")}
        {print; if($1=="fld"){for(i=1;i<=n;i++)print "map " md "/" fb "." T[i] ".map"; print "elecmap " md "/" fb ".e.map"; print "desolvmap " md "/" fb ".d.map"}}' \
        "$DPF_FILE" > "$DPF_FILE.tmp" && mv "$DPF_FILE.tmp" "$DPF_FILE"

      # ligand / protein
      sed -i -E "s|^move\b.*|move $LIGAND_FILE|g" "$DPF_FILE"
      if grep -qE '^protein\s*=' "$DPF_FILE"; then
        sed -i -E "s|^protein\s*=.*|protein = $PROTEIN_FILE|g" "$DPF_FILE"
      fi

      # presupuesto GA
      tors=$(awk '/^TORSDOF/{print $2; exit}' "$LIGAND_FILE"); [[ -z "$tors" ]] && tors=6
      compute_ga "$tors"
      if grep -qE '^ga_num_evals\b' "$DPF_FILE"; then sed -i -E "s|^ga_num_evals\b.*|ga_num_evals $GA_EVALS|g" "$DPF_FILE"; else echo "ga_num_evals $GA_EVALS" >> "$DPF_FILE"; fi
      if grep -qE '^ga_pop_size\b' "$DPF_FILE"; then sed -i -E "s|^ga_pop_size\b.*|ga_pop_size $GA_POP|g" "$DPF_FILE"; else echo "ga_pop_size $GA_POP" >> "$DPF_FILE"; fi
      if grep -qE '^sw_max_its\b' "$DPF_FILE"; then sed -i -E "s|^sw_max_its\b.*|sw_max_its $SW_ITS|g" "$DPF_FILE"; else echo "sw_max_its $SW_ITS" >> "$DPF_FILE"; fi
      if grep -qE '^ga_run\b' "$DPF_FILE"; then sed -i -E "s|^ga_run\b.*|ga_run $GA_RUNS|g" "$DPF_FILE"; else echo "ga_run $GA_RUNS" >> "$DPF_FILE"; fi
      sed -i -E "/^torsdof\b.*/d" "$DPF_FILE"

      echo "    Dockeando: $LIGAND_NAME -> $frame_base (tors=$tors, evals=$GA_EVALS, runs=$GA_RUNS)"

      (
        if "$AUTODOCK_BIN" -p "$DPF_FILE" -l "$DLG_FILE" 2>/dev/null; then
          rm -f "$DPF_FILE"
        else
          mv "$DPF_FILE" "$DPF_FILE.ERROR" 2>/dev/null
          echo "    ERROR: AutoDock falló o no generó $DLG_FILE" >&2
        fi
      ) &
      wait_for_slot
    done

    echo "    Esperando jobs del frame $frame_base..."
    while (( $(jobs -rp | wc -l) > 0 )); do
      wait -n 2>/dev/null || true
    done

    have_count=$(ls "$FRAME_OUT"/*.dlg 2>/dev/null | wc -l)
    echo "  > Frame $frame_base terminado. Progreso: $have_count / $TOTAL_LIGANDS DLGs"
  done

  echo "---- Fin estructura: $struct_name ----"
  echo
done < <(find "$PROTEIN_ROOT" -mindepth 1 -maxdepth 1 -type d -print0)

echo "==== FINISHED run $(date '+%F %T') ===="
echo "===> Resultados en: $OUT_ROOT"
