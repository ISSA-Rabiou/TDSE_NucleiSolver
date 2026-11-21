#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd -- "$SCRIPT_DIR/.." && pwd)"
if [[ $# -ne 2 ]]; then
  echo "Usage: $0 NB DT" >&2
  exit 2
fi
nb=$1
dt=$2
nq=$((nb + 5))
result_dir="$SCRIPT_DIR/results_H2OPJT_Hag_${nb}_${dt}"
mkdir -p "$result_dir"
cd "$result_dir"
rm -f E.txt Norm.txt Qt.txt SQt.txt Imp_k.txt pop.txt alpha.txt auto_cor.txt E_13.txt Norm_13.txt validation.txt

"$ROOT_DIR/TD_SCHROD.x" > "resultat_H2OPJT_Hag_${nb}_${dt}.log" 2>&1 <<EOF

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
&defWP0 sigma=0.25556361271643913 Beta=0.0 Qeq=1.9 imp_k=0.0 gamma=0.0 /
&defWP0 sigma=0.25556361271643907 Beta=0.0 Qeq=1.9 imp_k=0.0 gamma=0.0 /
&defWP0 sigma=0.30696837793008830 Beta=0.0 Qeq=1.9 imp_k=0.0 gamma=0.0 /

&prop t0=0. tf=10000.0 delta_t=$dt eps=1.e-20
  max_iter=500 propa_name='hagedorn_taylor' Beta=T P=T renorm=f write_E=T write_auto_cor=T debug=F /
EOF

