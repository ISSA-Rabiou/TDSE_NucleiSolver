#!/usr/bin/env bash
set -u
ROOT=$(cd "$(dirname "$0")/.." && pwd)
OUT="$ROOT/test/smolyak_tests"
rm -rf "$OUT"
mkdir -p "$OUT"
make_case() {
  local molecule=$1 method=$2 src=$3 tf=$4 dt=$5 lb=$6 lg=$7
  local name=${molecule}_${method}
  awk 'BEGIN{on=0} /<< \*\*/{on=1;next} /^\*\*$/{if(on) exit} on{print}' "$ROOT/test/inputs0/$src" \
   | sed -e "s/tf = 50\./tf = $tf/" -e "s/tf=50/tf=$tf/" \
         -e "s/tf = 50/tf = $tf/" -e "s/tf=50/tf=$tf/" \
         -e "s/delta_t = 0.05/delta_t = $dt/" -e "s/delta_t=0.05/delta_t=$dt/" \
         -e "s/delta_t = 0.10/delta_t = $dt/" -e "s/delta_t=0.10/delta_t=$dt/" \
         -e "s/name='dp' nb_basis=/name='smolyak' nb_basis=/" \
         -e "s/nb_basis=2/nb_basis=2 LB=$lb LG=$lg A_smol=1 B_smol=1/" \
         -e "s/nb_basis=4/nb_basis=4 LB=$lb LG=$lg A_smol=1 B_smol=1/" \
         -e "s/propa_name='non_hagedorn'/propa_name='non_hagedorn'/" \
         -e "s/propa_name='hagedorn'/propa_name='hagedorn'/" > "$OUT/$name.inp"
  mkdir -p "$OUT/$name"
  (cd "$OUT/$name" && "$ROOT/TD_SCHROD.x" < "$OUT/$name.inp" > "$OUT/$name.stdout" 2>&1) || true
}
make_case H2 standard input_H2_developed_std.inp 0.15 0.05 4 4
make_case H2 hagedorn input_H2_developed_hagedorn.inp 0.15 0.05 4 4
make_case H2O standard input_H2O_developed_std.inp 0.30 0.10 10 12
make_case H2O hagedorn input_H2O_developed_hagedorn.inp 0.30 0.10 10 12
printf '%s\n' '=== résumé ==='
for f in "$OUT"/*.stdout; do
  echo "--- $f"
  grep -E 'nuclear_sparse_basis_size|nuclear_sparse_grid_size|Diagonalisation|<psi\|H\|psi>|validation max_|Norm diff|Taylor condition is fulfilled|ERROR|error termination|Fortran runtime error' "$f" | tail -80 || true
done
