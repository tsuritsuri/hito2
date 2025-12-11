#!/bin/bash

# Configuración
BASE_GRID_PATH="/home/db/Documentos/David_Hito2/Glide/grids_frames_glide"
LIGANDFILE="/home/db/Documentos/David_Hito2/Glide/ligprep_1/ligprep_1-out.maegz"
POSES_PER_LIG=10
POSTDOCK_NPOSE=10
NJOBS=20  # Número de trabajos paralelos
NCPUS=40  # Número total de CPUs disponibles

# Proteínas a procesar
PROTEINS=("1E9Y" "1E9Z")

# Directorio de trabajo base
BASE_WORK_DIR="$(pwd)/docking_results"
mkdir -p "$BASE_WORK_DIR"

# Log file principal
MAIN_LOG="$BASE_WORK_DIR/docking_complete.log"
echo "=== Inicio del docking batch completo: $(date) ===" | tee "$MAIN_LOG"
echo "" | tee -a "$MAIN_LOG"

# Función para correr dockings con una precisión específica
run_docking_precision() {
    local PRECISION=$1

    echo "╔════════════════════════════════════════╗" | tee -a "$MAIN_LOG"
    echo "║  INICIANDO DOCKINGS CON PRECISION: $PRECISION  ║" | tee -a "$MAIN_LOG"
    echo "╚════════════════════════════════════════╝" | tee -a "$MAIN_LOG"
    echo "" | tee -a "$MAIN_LOG"

    # Directorio de trabajo para esta precisión
    WORK_DIR="$BASE_WORK_DIR/${PRECISION}"
    mkdir -p "$WORK_DIR"

    # Log file para esta precisión
    LOG_FILE="$WORK_DIR/docking_progress_${PRECISION}.log"
    echo "=== Inicio docking $PRECISION: $(date) ===" | tee "$LOG_FILE"
    echo "Precision: $PRECISION" | tee -a "$LOG_FILE"
    echo "Poses por ligando: $POSES_PER_LIG" | tee -a "$LOG_FILE"
    echo "" | tee -a "$LOG_FILE"

    # Contador de trabajos
    TOTAL_JOBS=0
    COMPLETED_JOBS=0
    SKIPPED_JOBS=0

    # Iterar sobre cada proteína
    for PROTEIN in "${PROTEINS[@]}"; do
        GRID_FOLDER="$BASE_GRID_PATH/grids_${PROTEIN}"

        # Verificar que la carpeta de grillas existe
        if [ ! -d "$GRID_FOLDER" ]; then
            echo "⚠️  Carpeta no encontrada: $GRID_FOLDER" | tee -a "$LOG_FILE"
            continue
        fi

        echo "📁 Procesando proteína: $PROTEIN" | tee -a "$LOG_FILE"

        # Encontrar todos los archivos .zip en la carpeta
        GRID_FILES=($(find "$GRID_FOLDER" -name "*.zip" | sort))

        if [ ${#GRID_FILES[@]} -eq 0 ]; then
            echo "   ⚠️  No se encontraron archivos .zip en $GRID_FOLDER" | tee -a "$LOG_FILE"
            continue
        fi

        echo "   Encontrados ${#GRID_FILES[@]} archivos de grilla" | tee -a "$LOG_FILE"

        # Procesar cada archivo de grilla
        for GRID_FILE in "${GRID_FILES[@]}"; do
            # Extraer el nombre del archivo sin extensión
            GRID_BASENAME=$(basename "$GRID_FILE" .zip)

            # Crear nombre del directorio de salida
            OUTPUT_DIR="$WORK_DIR/${PROTEIN}-${GRID_BASENAME}"

            # Verificar si ya fue procesado
            if [ -f "$OUTPUT_DIR/docking_pv.maegz" ]; then
                echo "   ✓ Ya procesado: ${PROTEIN}-${GRID_BASENAME}" | tee -a "$LOG_FILE"
                ((SKIPPED_JOBS++))
                continue
            fi

            # Crear directorio de salida
            mkdir -p "$OUTPUT_DIR"

            # Crear archivo .in temporal
            INPUT_FILE="$OUTPUT_DIR/docking.in"
            cat > "$INPUT_FILE" << EOF
GRIDFILE $GRID_FILE
LIGANDFILE $LIGANDFILE
PRECISION $PRECISION
POSES_PER_LIG $POSES_PER_LIG
POSTDOCK_NPOSE $POSTDOCK_NPOSE
WRITEREPT YES
EOF

            echo "   🚀 Iniciando docking: ${PROTEIN}-${GRID_BASENAME}" | tee -a "$LOG_FILE"
            echo "      Grid: $GRID_FILE" | tee -a "$LOG_FILE"
            echo "      Output: $OUTPUT_DIR" | tee -a "$LOG_FILE"
            echo "      Inicio: $(date)" | tee -a "$LOG_FILE"

            ((TOTAL_JOBS++))

            # Cambiar al directorio de salida y correr Glide
            cd "$OUTPUT_DIR"

            # Correr Glide
            "${SCHRODINGER}/glide" docking.in -OVERWRITE -NJOBS $NJOBS -HOST localhost:$NCPUS -TMPLAUNCHDIR -WAIT

            # Verificar si terminó exitosamente
            if [ -f "docking_pv.maegz" ]; then
                echo "      ✓ Completado exitosamente: ${PROTEIN}-${GRID_BASENAME}" | tee -a "$LOG_FILE"
                echo "      Fin: $(date)" | tee -a "$LOG_FILE"
                ((COMPLETED_JOBS++))
            else
                echo "      ❌ ERROR: No se generó el archivo de salida" | tee -a "$LOG_FILE"
            fi

            echo "" | tee -a "$LOG_FILE"

            # Volver al directorio de trabajo
            cd "$WORK_DIR"
        done

        echo "" | tee -a "$LOG_FILE"
    done

    # Resumen para esta precisión
    echo "=== RESUMEN $PRECISION ===" | tee -a "$LOG_FILE"
    echo "Total de dockings programados: $TOTAL_JOBS" | tee -a "$LOG_FILE"
    echo "Completados exitosamente: $COMPLETED_JOBS" | tee -a "$LOG_FILE"
    echo "Ya procesados (saltados): $SKIPPED_JOBS" | tee -a "$LOG_FILE"
    echo "Fin $PRECISION: $(date)" | tee -a "$LOG_FILE"

    # Agregar al log principal
    echo "=== RESUMEN $PRECISION ===" | tee -a "$MAIN_LOG"
    echo "Completados: $COMPLETED_JOBS / $TOTAL_JOBS" | tee -a "$MAIN_LOG"
    echo "Fin: $(date)" | tee -a "$MAIN_LOG"
    echo "" | tee -a "$MAIN_LOG"
}

# ============================================
# EJECUTAR DOCKINGS: PRIMERO XP, LUEGO SP
# ============================================

echo "Iniciando workflow completo: XP → SP" | tee -a "$MAIN_LOG"
echo "" | tee -a "$MAIN_LOG"

# Primero XP
run_docking_precision "XP"

echo "" | tee -a "$MAIN_LOG"
echo "╔════════════════════════════════════════╗" | tee -a "$MAIN_LOG"
echo "║   XP COMPLETADO. INICIANDO SP...       ║" | tee -a "$MAIN_LOG"
echo "╚════════════════════════════════════════╝" | tee -a "$MAIN_LOG"
echo "" | tee -a "$MAIN_LOG"

# Luego SP
run_docking_precision "SP"

# Resumen final completo
echo "" | tee -a "$MAIN_LOG"
echo "╔════════════════════════════════════════╗" | tee -a "$MAIN_LOG"
echo "║  TODOS LOS DOCKINGS COMPLETADOS! 🎉   ║" | tee -a "$MAIN_LOG"
echo "╚════════════════════════════════════════╝" | tee -a "$MAIN_LOG"
echo "" | tee -a "$MAIN_LOG"
echo "Resultados guardados en:" | tee -a "$MAIN_LOG"
echo "  - XP: $BASE_WORK_DIR/XP/" | tee -a "$MAIN_LOG"
echo "  - SP: $BASE_WORK_DIR/SP/" | tee -a "$MAIN_LOG"
echo "" | tee -a "$MAIN_LOG"
echo "Fin total: $(date)" | tee -a "$MAIN_LOG"

# Listar estructura de resultados
echo "" | tee -a "$MAIN_LOG"
echo "📂 Estructura de resultados:" | tee -a "$MAIN_LOG"
tree -L 2 "$BASE_WORK_DIR" 2>/dev/null || find "$BASE_WORK_DIR" -maxdepth 2 -type d
