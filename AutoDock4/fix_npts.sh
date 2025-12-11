#!/bin/bash

# Definir la carpeta base
BASE_DIR="gpfs_2"

# Verificar si el directorio existe
if [ ! -d "$BASE_DIR" ]; then
    echo "Error: La carpeta '$BASE_DIR' no existe en este directorio."
    exit 1
fi

echo "Iniciando modificación de archivos .gpf en '$BASE_DIR'..."

# Explicación del comando:
# 1. find "$BASE_DIR" ... : Busca dentro de gpfs_2
# 2. -name "*.gpf" : Solo archivos que terminen en .gpf
# 3. -exec sed -i ... : Ejecuta sed en modo "in-place" (edición directa) sobre cada archivo encontrado
# 4. 's/^npts .*/.../' : Busca cualquier línea que empiece (^) con "npts" seguido de cualquier cosa (.*)
#                        y la sustituye por "npts 40 40 40" conservando el comentario típico.

find "$BASE_DIR" -type f -name "*.gpf" -exec sed -i 's/^npts .*/npts 40 40 40                        # num.grid points in xyz/' {} +

echo "¡Proceso terminado! Todos los archivos .gpf ahora tienen npts 40 40 40."
