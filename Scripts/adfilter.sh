#!/bin/bash
set -euo pipefail
IN="lig_sdf"
GOOD="molecules_ok"
BAD="molecules_unsupported"
mkdir -p "$GOOD" "$BAD"

# elementos “seguros” para AD4/Vina (ajusta si usas halógenos/metales comunes)
ALLOWED='^(H|C|N|O|S|P|F|Cl|Br|I)$'

for f in "$IN"/*.sdf; do
  # extrae símbolos atómicos del bloque de átomos
  syms=$(awk 'NR==4{n=$1} NR>4 && NR<=4+n {print $4}' "$f" | sort -u)
  ok=1
  while read s; do
    [[ "$s" =~ $ALLOWED ]] || { ok=0; break; }
  done <<< "$syms"
  if (( ok )); then
    mv "$f" "$GOOD/"
  else
    echo "Excluyendo $(basename "$f") por elemento no soportado: $syms"
    mv "$f" "$BAD/"
  fi
done

echo "OK en:   $GOOD"
echo "Excluidos: $BAD"

