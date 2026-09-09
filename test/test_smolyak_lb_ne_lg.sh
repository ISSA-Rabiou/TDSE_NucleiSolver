#!/usr/bin/env bash
#
# Test Smolyak avec LB != LG pour H2 et H2O.
#
# Pré-requis :
#   - exécutable ../TD_SCHROD.x présent depuis le répertoire du projet ;
#   - inputs0/input_H2*_developed_*.inp présents ;
#   - les sources acceptent uniquement les paramètres Smolyak dédiés :
#       A_b_smol, B_b_smol, A_g_smol, B_g_smol.
#
set -Eeuo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
EXE="$ROOT/TD_SCHROD.x"
INPUT_DIR="$ROOT/test/inputs0"
OUT="$ROOT/test/smolyak_tests_lb_ne_lg"

if [[ ! -x "$EXE" ]]; then
  echo "ERREUR : exécutable absent ou non exécutable : $EXE" >&2
  exit 1
fi

mkdir -p "$OUT"
rm -rf "$OUT"/*

make_case() {
  local molecule="$1"
  local method="$2"
  local src="$3"
  local tf="$4"
  local dt="$5"
  local lb="$6"
  local lg="$7"
  local name="${molecule}_${method}"
  local raw_input="$INPUT_DIR/$src"
  local input="$OUT/$name.inp"
  local run_dir="$OUT/$name"
  local stdout="$OUT/$name.stdout"

  if [[ ! -f "$raw_input" ]]; then
    echo "ERREUR : input absent : $raw_input" >&2
    return 1
  fi

  # On extrait uniquement le bloc de données situé entre << ** et **.
  # Puis on adapte les options historiques aux interfaces compilées.
  awk 'BEGIN { on=0 }
       /<< \*\*/ { on=1; next }
       /^\*\*$/  { if (on) exit }
       on         { print }' "$raw_input" \
    | sed \
        -e "s/tf = 50\./tf = $tf/" \
        -e "s/tf=50/tf=$tf/" \
        -e "s/delta_t = 0.05/delta_t = $dt/" \
        -e "s/delta_t=0.05/delta_t=$dt/" \
        -e "s/delta_t = 0.10/delta_t = $dt/" \
        -e "s/delta_t=0.10/delta_t=$dt/" \
        -e "s/name='dp' nb_basis=/name='smolyak' nb_basis=/" \
        -e "s/nb_basis=2/nb_basis=2 LB=$lb LG=$lg A_b_smol=1 B_b_smol=1 A_g_smol=1 B_g_smol=1/" \
        -e "s/nb_basis=4/nb_basis=4 LB=$lb LG=$lg A_b_smol=1 B_b_smol=1 A_g_smol=1 B_g_smol=1/" \
        -e "/Cart_Type[[:space:]]*=/d" \
        -e "s/option=3/option=1/" \
    > "$input"

  mkdir -p "$run_dir"
  echo "RUN  $name  LB=$lb  LG=$lg  tf=$tf  dt=$dt"
  (
    cd "$run_dir"
    "$EXE" < "$input" > "$stdout" 2>&1
  )

  if grep -Eiq 'STOP ERROR|Fortran runtime error|Segmentation fault|Backtrace for this error|(^|[[:space:]])ERROR([[:space:]:]|$)' "$stdout"; then
    echo "FAIL $name : erreur détectée dans $stdout" >&2
    return 1
  fi
  if ! grep -q 'validation max_norm_drift' "$stdout"; then
    echo "FAIL $name : validation de norme absente" >&2
    return 1
  fi
  if ! grep -q 'validation max_energy_drift' "$stdout"; then
    echo "FAIL $name : validation d’énergie absente" >&2
    return 1
  fi

  echo "PASS $name"
}

make_case H2  standard input_H2_developed_std.inp    0.15 0.05  4  6
make_case H2  hagedorn input_H2_developed_hagedorn.inp 0.15 0.05  4  6
make_case H2O standard input_H2O_developed_std.inp   0.30 0.10 10 12
make_case H2O hagedorn input_H2O_developed_hagedorn.inp 0.30 0.10 12 14

printf '\n=== validation LB != LG réussie ===\n'
for stdout in "$OUT"/*.stdout; do
  echo "--- $stdout"
  grep -E 'order_LB|order_LG|nb_SG|nq_SG|<psi\|H\|psi>|validation max_|Norm diff|Taylor condition is fulfilled' "$stdout" | tail -80 || true
done
