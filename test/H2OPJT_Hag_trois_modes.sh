#!/usr/bin/env bash
set -euo pipefail

# Valeurs optimisées de reference :
Q1=1.8102076015368764
Q2=1.8102076015368764
Q3=1.8238629061636986

# Deplacements du centre initial pour combiner les trois modes.
delta_sym=0.08
delta_antisym=0.08
delta_flexion=0.1

Qeq1=$(echo "$Q1 + $delta_sym + $delta_antisym" | bc -l)
Qeq2=$(echo "$Q2 + $delta_sym - $delta_antisym" | bc -l)
Qeq3=$(echo "$Q3 + $delta_flexion" | bc -l)

echo "Qeq1=$Qeq1"
echo "Qeq2=$Qeq2"
echo "Qeq3=$Qeq3"

mode="trois_modes"
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd -- "$SCRIPT_DIR/.." && pwd)"
if [[ $# -ne 2 ]]; then
  echo "Usage: $0 NB DT" >&2
  exit 2
fi
nb=$1
dt=$2
nq=$((nb + 5))
result_dir="$SCRIPT_DIR/results_H2OPJT_hagedorn_${mode}_${nb}_${dt}"
mkdir -p "$result_dir"
cd "$result_dir"

"$ROOT_DIR/TD_SCHROD.x" > "resultat_H2OPJT_Hagedorn_${mode}_${nb}_${dt}.log" 2>&1 <<EOF

&coordQMLpot
coord_type = 'CURV'
construction_type = 'developed'
sym_type   = 'SYM'
model_type = 'QML'
pot0       = 0.0
/

&system
  printlev = 2
/

&variables
  Old_Qtransfo = F
  nb_Qtransfo  = 2
  WriteT = F
  Cart_Type = "BF"
  JJ=0
  Tana=t
  FortranForm=f
/

&Coord_transfo name_transfo='zmat' nat=3 /
O
H 1
H 1 2

&Coord_transfo name_transfo='active' /
1  1  1

&minimum /
R1 1.
R2 1.
a2 1.83

&potential
  pot_name  = 'H2O_PJT'
  option=1
  adiabatic = f
  nsurf = 1
  ndim = 3
/

&basis_nd name='dp' nb_basis=4 /
&basis_nd name='HO' nb=$nb nq=$nq use_initial_gwp=T debug=F /
&basis_nd name='HO' nb=$nb nq=$nq use_initial_gwp=T debug=F /
&basis_nd name='HO' nb=$nb nq=$nq use_initial_gwp=T debug=F /
&basis_nd name='el' nb=1 /

&defGWP ndim=3 Elecindex=1 Coef=(1.,0.) /
&defWP0 sigma=0.25556361271643913 Beta=0.0 Qeq=$Qeq1 imp_k=0.0 gamma=0.0 /
&defWP0 sigma=0.25556361271643907 Beta=0.0 Qeq=$Qeq2 imp_k=0.0 gamma=0.0 /
&defWP0 sigma=0.30696837793008830 Beta=0.0 Qeq=$Qeq3 imp_k=0.0 gamma=0.0 /

&prop t0=0. tf=100000.0 delta_t=$dt eps=1.e-20
  max_iter=500 propa_name='hagedorn_taylor' Beta=T P=T renorm=T write_E=T write_auto_cor=T debug=F /
EOF
