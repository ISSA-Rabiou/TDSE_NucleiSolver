#!/bin/bash
# Compile et calcule le spectre.
#
# Sans argument :
#   ./spectre.sh
# utilise ./auto_cor.txt et crée ./file_spectrum.txt.
#
# Avec arguments :
#   ./spectre.sh fichier_auto_cor fichier_spectre dE Emin Emax n_filter
# Exemple :
#   ./spectre.sh /donnees/auto_cor.txt /resultats/file_spectrum.txt 1.0e-5 0.0 0.5 1
#
# Si dE est omis ou vaut -1, le programme le calcule automatiquement.

set -e

# Dossier contenant spectre.f90 et le programme compilé.
PROGRAM_DIR=$(cd "$(dirname "$0")" && pwd)

# Les fichiers de données sont cherchés dans le dossier courant de l'utilisateur.
INPUT_NAME=${1:-auto_cor.txt}
OUTPUT_NAME=${2:-file_spectrum.txt}
DE=${3:--1.0}
EMIN=${4:-0.0}
EMAX=${5:--1.0}
NFILTER=${6:-1}

# Résolution portable du chemin d'entrée : compatible Linux et macOS.
INPUT_DIR=$(cd "$(dirname "$INPUT_NAME")" && pwd)
INPUT_FILE="$INPUT_DIR/$(basename "$INPUT_NAME")"

if [ ! -f "$INPUT_FILE" ]; then
    echo "Erreur : fichier introuvable : $INPUT_FILE"
    echo "Placez auto_cor.txt dans le dossier courant ou donnez son chemin."
    exit 1
fi

# Résolution portable du chemin de sortie.
OUTPUT_DIR=$(dirname "$OUTPUT_NAME")
mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR=$(cd "$OUTPUT_DIR" && pwd)
OUTPUT_FILE="$OUTPUT_DIR/$(basename "$OUTPUT_NAME")"
PARAMETER_FILE="$OUTPUT_DIR/spectrum_parameters.txt"

if ! command -v gfortran >/dev/null 2>&1; then
    echo "Erreur : gfortran n'est pas installe."
    exit 1
fi

mkdir -p "$PROGRAM_DIR/bin"

# Compilation automatique.
gfortran -O2 -std=f2008 -Wall "$PROGRAM_DIR/spectre.f90" -o "$PROGRAM_DIR/bin/spectre"

# Parametres transmis au programme Fortran.
PARAM_FILE="${TMPDIR:-/tmp}/spectre_$$.nml"
cat > "$PARAM_FILE" <<EOF
&param
  file_auto       = '$INPUT_FILE'
  file_spectrum   = '$OUTPUT_FILE'
  file_parameters = '$PARAMETER_FILE'
  dE              = $DE
  Emin            = $EMIN
  Emax            = $EMAX
  conv_E          = 1.0
  conv_t          = 1.0
  n_filter        = $NFILTER
/
EOF

"$PROGRAM_DIR/bin/spectre" < "$PARAM_FILE"
rm -f "$PARAM_FILE"

# Informations utilisées par le notebook.
python3 - "$INPUT_FILE" "$OUTPUT_FILE" "$PARAMETER_FILE" <<'PY'
import json
import sys
from pathlib import Path

manifest = Path(sys.argv[2]).parent / "spectrum_run.json"
manifest.write_text(json.dumps({
    "autocorrelation_file": sys.argv[1],
    "spectrum_file": sys.argv[2],
    "parameter_file": sys.argv[3],
}, indent=2, ensure_ascii=False) + "\n")
print("Manifeste cree :", manifest)
PY
