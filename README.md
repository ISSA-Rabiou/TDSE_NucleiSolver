# TD_SCHROD — Quantum molecular dynamics

TD_SCHROD* is a free software under the MIT Licence.

## Authors, contributors, and scientific context

```
    Copyright (c) 2026 R. ISSA [1]

      with contributions of:
        LAUVERGNAT David [2]
        SODOGA Komi [1]
```

[1]: Université de Lomé, Togo, FDS, Department of Physics
[2]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France


`TD_SCHROD` is a Fortran code for quantum molecular dynamics. It computes the time evolution of a wavefunction by solving the time-dependent Schrödinger equation (**TDSE**) and stationary states by solving the time-independent Schrödinger equation (**TISE**).

The code uses the QDUtil, Tnum-Tana, and QuantumModelLib, etc. libraries and, depending on the configuration, LAPACK/BLAS and ARPACK-NG. Potentials can be provided by the built-in models or by QuantumModelLib (QML).

## 1. General calculation organization

The program `TD_SCHROD.x` reads input from standard input. An input file can therefore be run with shell redirection:

```bash
./TD_SCHROD.x < mon_calcul.inp > mon_calcul.log 2>&1
```

The program follows these steps :

1. read the potential model and molecular coordinates ;
2. read and construct the primitive basis ;
3. read the initial Gaussian wavepacket GWP ;
4. read the namelist `&prop` ;
5. branch to the calculation `TISE` or `TDSE` according to `dynamics` ;

The choice `dynamics='TISE'` is handled before propagation. There is therefore no implicit diagonalization during propagation TDSE. Conversely, `dynamics='TDSE'` launches only time propagation.

## 2. Compilation

### 2.1 Linux compilation

From the project root :

```bash
make -j1 FC=gfortran OMP=1 LAPACK=1 ARPACK=1 QML=1
```

Main options :

| Option | values | Function |
|---|---:|---|
| `FC` | `gfortran`, `ifort`, `ifx`, ... | Fortran compiler |
| `OMP` | `0` or `1` | disables or enables OpenMP |
| `LAPACK` | `0` or `1` | disables or enables LAPACK/BLAS |
| `ARPACK` | `0` or `1` | disables or enables ARPACK |
| `QML` | according to the project scripts | uses the QuantumModelLib libraries |
| `OPT` | `0` or `1` | project optimization level |
| `INT` | `4` or `8` | default integer |
| `RKIND` | `real32`, `real64`, `real128` | real precision |

LAPACK is required for dense diagonalization and small projected problems. ARPACK is required only for the TISE `arpack` method.

When `ARPACK=1`, the Makefile automatically downloads ARPACK-NG 3.9.1,
compiles it with `FC`, and installs it locally in `EXT_Lib/arpack-ng_loc`.
The same command therefore works on Linux, macOS and other Unix-like systems;
it does not use Homebrew or another system package manager:

```bash
make -j1 FC=gfortran OMP=1 LAPACK=1 ARPACK=1
```

The first build requires a Fortran compiler, CMake, `curl` or `wget`, `tar`,
and BLAS/LAPACK. With `ARPACK=0`, ARPACK is neither downloaded nor linked.
On macOS, the Makefile automatically uses `-force_load` for the single
monolithic `libTnum-TanaFull.a` archive. This is necessary because Apple `ld`
may not rescan all archive members and can otherwise report unresolved QDUtil
symbols. QDUtil must not be linked a second time: it is already included in
that archive, and linking it separately would create duplicate symbols. On
other Unix-like systems, normal archive extraction is used. If necessary,
`FORCE_LOAD=0` or `FORCE_LOAD=1` can still be supplied explicitly.

The main target is :

```text
TD_SCHROD.x
```

To clean local objects :

```bash
make clean
```

To also remove the executable:

```bash
make distclean
```

To rebuild external libraries when necessary :

```bash
make cleanall RECCLEAN=1
make FC=gfortran OMP=1 LAPACK=1 ARPACK=1 QML=1
```

## 2. Minimal input-file structure

An input file generally contains the following blocks, in this order:

```text
&coordQMLpot
  ...
/

&potential
  ...
/

&basis_nd ... /
&basis_nd ... /

&defGWP ... /
&defWP0 ... /
&defWP0 ... /

&prop
  ...
/
```

For molecular models using Tnum-Tana, the blocks `&system`, `&variables`, `&Coord_transfo`, and optionally `&minimum` are placed before `&potential`, as shown by the H2O and ClH2+ examples.

## 3. Potential and coordinates

### 3.1 Coordinate representations and QML models

The `&coordQMLpot` namelist selects the coordinate representation used by the operator and by the QML potential interface. `coord_type='CC'` denotes Cartesian or normal coordinates in the constant/standard kinetic representation implemented by this code. `coord_type='CURV'` denotes curvilinear coordinates: Tnum/Tana supplies coordinate transformations, metric information, Jacobian/volume factors, and the corresponding kinetic-energy data.

The `sym_type` value describes the coordinate representation expected by the model. It must be interpreted together with `pot_itQtransfo`; `sym_type` alone does not perform a transformation. For a symmetrized coordinate representation, use the relevant Tnum transformation stage, commonly `sym_type='SYM'` with `pot_itQtransfo=1`. To keep the active-coordinate representation, use `pot_itQtransfo=-1`. More generally, non-negative values select a Tnum transformation stage. `pot_itQtransfo` therefore controls where the QML potential is evaluated along the coordinate-transformation chain: it must match the coordinates expected by the QML model.

A CURV example using a symmetrized representation is:

```fortran
&coordQMLpot
  coord_type       = 'CURV'
  sym_type         = 'SYM'
  pot_itQtransfo   = 1
  construction_type= 'developed'
  model_type       = 'QML'
  pot0             = 0.0
/
```

A CURV example that evaluates the potential directly in active coordinates is:

```fortran
&coordQMLpot
  coord_type       = 'CURV'
  sym_type         = 'SYM'
  pot_itQtransfo   = -1
  construction_type= 'developed'
  model_type       = 'QML'
  pot0             = 0.0
/
```

The important fields are `coord_type`, `pot_itQtransfo`, `construction_type`, `model_type='QML'`, `pot_name`, `option`, `adiabatic`, `nsurf`, and `ndim`. The exact `pot_name` and `option` must be supported by QuantumModelLib.

### 3.2 CURV kinetic-energy action formulations

For `coord_type='CURV'`, `src/Op_m.f90` implements three alternative numerical formulations for applying the same physical kinetic-energy operator to a grid wavefunction. They are alternative discretization/action routes, not three different Hamiltonians:

- `construction_type='developed'` is the expanded F2/F1/Vep form. It stores and uses the coefficient arrays F2, F1, `rho`, and the effective potential, including diagonal and mixed derivative contributions. This is the established/reference route and the default.
- `construction_type='divergence'` is the direct divergence/flux form. It uses the metric tensor G(Q) and `rho`, applies derivatives to complete flux vectors, and avoids hand-expanded mixed-term bookkeeping.
- `construction_type='similarity'` is the similarity-transformed/direct flux form. It evaluates metric/Jacobian data on demand and applies the corresponding `rho`/Jacobian factors in the transformed action.

Select the formulation in `&coordQMLpot`, for example `construction_type='divergence'`. Use `developed` for the established production/reference route; use `divergence` or `similarity` as independent cross-checks and for curvilinear-operator validation. Compare basis/grid convergence, energies, norms, and observable histories rather than comparing a single run.

The `split` propagator is restricted by the code to `coord_type='CC'` with a fixed periodic Fourier grid and a constant Cartesian kinetic operator. It cannot be combined with a moving Hagedorn basis.

### 3.3 Molecular coordinates, units, and isotopes

For a molecular model, a Z-matrix transformation can be defined as follows:

```fortran
&Coord_transfo name_transfo='zmat' nat=3 /
35Cl
1H 1
1H 1 2

&Coord_transfo name_transfo='active' /
1 1 1
```

The notation `35Cl` or `1H` explicitly specifies the isotope. In atomic-unit molecular calculations, typical coordinate conventions are bohr for lengths and radians for angles; the Tnum input should state isotopes explicitly when a particular isotopologue is intended. In the Botschwina ClH2+ model, the active coordinates are `Q(1)=R1` in bohr, `Q(2)=R2` in bohr, and `Q(3)=theta` in radians.

## 4. Basis construction

The basis is read with the namelist `&basis_nd`. Construction is hierarchical : a composite basis is declared with `nb_basis`, then its sub-bases are read from the following blocks.

### 4.1 Direct-product basis

For `ndim=3` and one electronic surface :

```fortran
&basis_nd name='dp' nb_basis=4 /
&basis_nd name='HO' nb=20 nq=30 Q0=0.0 scaleQ=1.0 Imp_k=0.0 Alpha=(1.0,0.0) /
&basis_nd name='HO' nb=20 nq=30 Q0=0.0 scaleQ=1.0 Imp_k=0.0 Alpha=(1.0,0.0) /
&basis_nd name='HO' nb=20 nq=30 Q0=0.0 scaleQ=1.0 Imp_k=0.0 Alpha=(1.0,0.0) /
&basis_nd name='el' nb=1 /
```

`nb_basis=4` means here three nuclear bases and one electronic basis. The electronic basis must be placed last in this composite structure.

The dimension of the product basis is essentially the product of the nuclear `nb` values multiplied by the number of electronic surfaces. A larger basis rapidly increases the cost of propagation and especially the cost of dense diagonalization.

### 4.2 Available primitive bases

The code currently constructs the following primitive types :

| `name` | Description | Main parameters |
|---|---|---|
| `HO` or `hag` | harmonic/Hagedorn basis | `nb`, `nq`, `Q0`, `scaleQ`, `Imp_k`, `Alpha` |
| `fourier` | periodic Fourier basis | `nb`, `nq` |
| `boxab` | sine basis in `[A,B]` | `nb`, `nq`, `A`, `B` |
| `el` | basis electronic | `nb` |

Names are converted to lowercase by the code.

### 4.3 Primitive-basis parameters

- `nb` : number of basis functions ;
- `nq` : number of grid/quadrature points ;
- `Q0` : basis center or origin ;
- `scaleQ` : coordinate scale factor ;
- `Imp_k` : initial momentum ;
- `Alpha=(a,b)` : complex parameter of the Hagedorn basis, with real part `a` and imaginary part `b` ;
- `A`, `B` : bounds for `boxab` ;
- `use_initial_gwp=T` : requests that the geometric parameters of the HO/Hagedorn basis be determined from the initial GWP;
- `debug=F` : disables detailed diagnostics for this basis. Use `debug=T` only to inspect its construction.

### 4.4 QML, initial GWP, and HO/Hagedorn basis coupling

In the current version, the HO/Hagedorn basis is not necessarily defined independently of the potential. For a QML potential, the program follows this sequence:

1. `InitQModelBasisParam` reads the initial-geometry parameters associated with the QML model;
2. `Read_tab_GWP` reads the initial wavepacket and its parameters `Qeq`, `sigma`, `Beta`, `imp_k` and `gamma` ;
3. if the GWP comes from the QML model, QML centers and widths are injected into the GWP ;
4. `Apply_Initial_GWP_To_Basis` transfers GWP parameters to primitive bases that have `use_initial_gwp=T` ;
5. `construct_primitive_basis` finally constructs the functions HO/Hagedorn with these updated parameters.

Thus, `use_initial_gwp=T` must be used when a basis must follow the initial wavefunction linked to the QML potential. The recommended value for usual calculations is:

```fortran
&basis_nd name='HO' nb=12 nq=12 use_initial_gwp=T debug=F /
```

In this case, the values `Q0`, `scaleQ`, `Imp_k`, and `Alpha` of the basis must not be considered independent of the GWP. The explicit values written in `&defWP0` (`Qeq`, `sigma`, etc.) and the parameters returned by QML become the physical reference for the basis. This convention is used in the H2O, HONO, and ClH2+ tests.

Example for three nuclear coordinates and one electronic surface:

```fortran
&basis_nd name='dp' nb_basis=4 /
&basis_nd name='HO' nb=12 nq=12 use_initial_gwp=T debug=F /
&basis_nd name='HO' nb=12 nq=12 use_initial_gwp=T debug=F /
&basis_nd name='HO' nb=12 nq=12 use_initial_gwp=T debug=F /
&basis_nd name='el' nb=1 /
```

The keyword `debug=F` is used to print more or less information during basis reading and construction. It does not disable the choice `debug` of the namelist `&prop`, which must also be set separately.

### 4.5 Smolyak sparse basis

A sparse composite basis can be requested with `name='smolyak'` or `name='smoly'` :

```fortran
&basis_nd name='smolyak' nb_basis=4 LB=5 LG=5 \
  A_b_smol=1 B_b_smol=1 A_g_smol=1 B_g_smol=1 /
&basis_nd name='HO' nb=20 nq=25 Q0=0.0 scaleQ=1.0 /
&basis_nd name='HO' nb=20 nq=25 Q0=0.0 scaleQ=1.0 /
&basis_nd name='HO' nb=20 nq=25 Q0=0.0 scaleQ=1.0 /
&basis_nd name='el' nb=1 /
```

The parameters `LB` and `LG` control the basis and grid levels, respectively. The parameters `A_b_smol`, `B_b_smol`, `A_g_smol`, and `B_g_smol` control the selection rules. The `retinal_smoly*.inp` examples and `test_smolyak_*.sh` scripts show complete configurations.

The sparse basis reduces the effective dimension, but requires checking potential, grid, and operator compatibility. The program calls a dedicated validation when the Smolyak representation is active.

## 5. Initial Gaussian wavepacket

### 5.1 Global definition

```fortran
&defGWP ndim=3 Elecindex=1 Coef=(1.0,0.0) /
```

- `ndim` : number of nuclear coordinates;
- `Elecindex` : initial electronic surface;
- `Coef` : initial complex coefficient.

### 5.2 Per-coordinate definition

A `&defWP0` block must be provided for each nuclear coordinate:

```fortran
&defWP0 sigma=0.20 Beta=0.0 Qeq=2.46 imp_k=0.0 gamma=0.0 /
&defWP0 sigma=0.20 Beta=0.0 Qeq=2.46 imp_k=0.0 gamma=0.0 /
&defWP0 sigma=0.20 Beta=0.0 Qeq=1.64 imp_k=0.0 gamma=0.0 /
```

- `sigma` : Gaussian width ;
- `Beta` : imaginary part related to the width ;
- `Qeq` : initial position ;
- `imp_k` : initial momentum ;
- `gamma` : initial phase.

The program always projects the initial GWP onto the selected basis. The basis can therefore differ from the initial wavepacket, especially for Hagedorn calculations.

## 6. Choosing the calculation type: TDSE or TISE

The choice is made in `&prop` with `dynamics` :

```fortran
dynamics='TDSE'
```

or :

```fortran
dynamics='TISE'
```

Values are case-insensitive. The default value is `TDSE`.

## 6. TDSE calculation: time propagation

### 6.1 Time parameters

```fortran
&prop
  t0=0.0
  tf=60.0
  delta_t=0.1
  eps=1.e-20
  max_iter=500
  dynamics='TDSE'
  propa_name='hagedorn_taylor'
  Beta=T
  P=T
  renorm=F
  debug=F
/
```

- `t0`, `tf` : initial and final times ;
- `delta_t` : external time step ;
- `eps` : internal propagator tolerance ;
- `max_iter` : maximum number of internal iterations ;
- `Beta` : evolution of the imaginary part of the Hagedorn parameter;
- `P` : evolution of the Hagedorn-basis momentum;
- `renorm` : wavefunction renormalization.

`tf=t0` is accepted to perform initialization and diagnostics without performing a propagation step.

### 6.2 Propagator selection

The recommended syntax is a compound name :

```fortran
propa_name='standard_taylor'
```

or :

```fortran
propa_name='hagedorn_taylor'
```

The part before `_` selects the propagated basis type :

- `standard` : fixed/non-Hagedorn basis ;
- `hagedorn` : mobile Hagedorn basis.

The part after `_` selects the integrator. The integrators accepted by the code are :

| Integrator | value of `propa_name` |
|---|---|
| Runge-Kutta order 1 | `rk1` |
| Runge-Kutta order 2 | `rk2` |
| Runge-Kutta order 3 | `rk3` |
| Runge-Kutta order 4 | `rk4` |
| Taylor | `taylor` |
| SIL | `sil` |
| Crank-Nicolson | `cn` |
| Split operator | `split` |
| Chebyshev | `cheby` |
| VP | `vp` |
| ITP | `itp` |
| GSSIL | `gssil` |

Examples:

```fortran
propa_name='standard_taylor'
propa_name='hagedorn_rk4'
propa_name='standard_cn'
propa_name='standard_split'
```

The old two-word syntax remains accepted:

```fortran
propa_name='hagedorn'
propa_name2='taylor'
```

For the `split` propagator, use a fixed periodic Fourier basis. `split_order` can be `1`, `2`, `4`, or `6` and is meaningful only for `split`:

```fortran
propa_name='standard_split'
split_order=4
```

The `split` propagator must not be combined with a moving Hagedorn basis.

### 7.3 TDSE output files

The following flags control output :

```fortran
write_Qt=T
write_SQt=T
write_Pt=T
write_E=T
write_Norm=T
write_pop=T
write_alpha=T
write_auto_cor=T
write_validation=T
```

Exact filenames are handled by the analysis and propagation routines. `debug=T` adds detailed diagnostics to standard output.

Restart options are :

```fortran
restart=F
checkpoint=T
checkpoint_every=10
checkpoint_file='TD_SCHROD.restart'
```

## 8. Propagation controls, diagnostics, and restart

The `propa_t` type is populated from the `&prop` namelist in `src/sub_propa_m.f90`. The time controls are `t0`, `tf`, and `delta_t`; `eps` is the convergence tolerance used by iterative propagators, while `max_iter` limits internal iterations. `E_min` and `E_max` bound the Chebyshev energy interval. `split_order` selects order 1, 2, 4, or 6 for `split` and must remain at its default value for other propagators.

`propa_name` is a combined selector such as `standard_taylor` or `hagedorn_taylor`: the prefix selects a fixed (`standard`) or moving Hagedorn basis, and the suffix selects the integrator. Supported suffixes include `rk1`, `rk2`, `rk3`, `rk4`, `taylor`, `sil`, `cn`, `split`, `cheby`, `vp`, `itp`, and `gssil`. The legacy `propa_name2` field remains accepted only for compatibility with the older two-field syntax; new inputs should use the combined `propa_name`. `Beta`, `P`, and `renorm` control Hagedorn-parameter evolution, Hagedorn momentum evolution, and optional renormalization.

The calculation selector is `dynamics='TDSE'` or `dynamics='TISE'`. `tise_method` selects `dense`/`lapack`, `lanczos`, `davidson`, or `arpack`; `tise_nstates`, `tise_max_iter`, and `tise_tol` control the stationary solver. `energy_unit` accepts `au` or `cm-1` (also `cm1` as an input alias). `tise_energy_file` and `tise_transition_file` set the stationary-state and transition output files. Set `tise_relax=T` to precede the TISE solve with imaginary-time relaxation; its three controls are `tise_relax_delta_t`, `tise_relax_max_iter`, and `tise_relax_tol`.

`action_H=T` requests Hamiltonian-action diagnostics and `writeV=T` requests potential output where supported. `debug` controls propagation diagnostics. The observable switches are `write_Qt`, `write_SQt`, `write_Pt`, `write_E`, `write_Norm`, `write_pop`, `write_alpha`, `write_auto_cor`, `write_validation`, `write_E_13`, and `write_Norm_13`. They write coordinate, width, momentum, energy, norm, population, Hagedorn-parameter, autocorrelation, validation, and auxiliary diagnostic histories, respectively.

### 8.1 Checkpoint and restart

With `checkpoint=T`, the real-time propagation writes `checkpoint_file` periodically and again at the final step. `checkpoint_every` is the interval in propagation steps. With `restart=T`, the program loads that file and continues from the saved time and step; observable histories are reconciled and truncated before continuation so that the resumed run does not duplicate later records. The checkpoint stores propagation state and reference diagnostics, and compatibility is checked using the input signature, propagation method/integrator, `t0`/`delta_t`/`eps`, numerical controls, basis dimensions, and associated metadata.

A restart must use the same pre-`&prop` input: the same potential, coordinates and Tnum transformation, basis, and initial wavepacket. Propagation settings must remain compatible. `tf` may be extended to continue beyond the saved endpoint, but it must not be earlier than the checkpoint time and must remain aligned with `delta_t`. ITP and GSSIL are standalone modes and cannot be restarted as a normal real-time propagation. This mechanism makes long calculations recoverable after interruption and extendable without restarting from `t0`.

A realistic restartable TDSE namelist is:

```fortran
&prop
  t0=0.0, tf=200.0, delta_t=0.05, eps=1.e-20
  E_min=-0.5, E_max=2.0, max_iter=500
  propa_name='hagedorn_taylor', Beta=T, P=T, renorm=F
  dynamics='TDSE', debug=F
  write_E=T, write_Norm=T, write_validation=T
  action_H=F, writeV=F
  checkpoint=T, checkpoint_every=20
  checkpoint_file='TD_SCHROD.restart', restart=T
/
```

To start a new calculation, use `restart=F` with the same checkpoint filename or choose a new filename. To extend a completed run, keep the compatible input and increase only `tf` as needed.

## 9. TISE calculation: stationary states

### 9.1 Minimal example

```fortran
&prop
  t0=0.0
  tf=0.0
  delta_t=1.0
  eps=1.e-10
  max_iter=500
  propa_name='standard_taylor'
  Beta=F
  P=F
  renorm=F
  dynamics='TISE'
  tise_method='arpack'
  tise_nstates=3
  tise_max_iter=500
  tise_tol=1.e-10
  energy_unit='cm-1'
  tise_energy_file='TISE_energies.dat'
  tise_transition_file='TISE_transitions.txt'
  debug=F
/
```

In TISE mode, `propa_name` is not used to propagate the wavefunction. It is retained in the input to preserve the common structure of the `&prop` namelist.

### 9.2 TISE methods

The module `TISE_m` provides the following methods :

| Method | Keyword | Usage |
|---|---|---|
| Dense LAPACK diagonalization | `dense` or `lapack` | small and medium bases; constructs `H%CMat` |
| Lanczos | `lanczos` | large bases; Hamiltonian applied in matrix-free form |
| Block-preconditioned Davidson | `davidson` | multiple states; diagonal preconditioner |
| ARPACK | `arpack` | large bases ; reverse communication and Hamiltonian matrix-free |

Main parameters :

- `tise_nstates` : number of requested states;
- `tise_max_iter` : maximum number of iterations of the iterative method;
- `tise_tol` : residual tolerance ;
- `energy_unit='au'` : energies in atomic units;
- `energy_unit='cm-1'` : energies in cm⁻¹ ; `cm1` is also accepted ;
- `tise_energy_file` : energy file ;
- `tise_transition_file` : transition file.

For a basis of dimension `N`, the method dense stores a matrix of size `N x N`. For a large basis, prefer `lanczos`, `davidson` or `arpack`.

### 9.3 Imaginary-time relaxation

Ground-state search by imaginary-time relaxation is enabled with :

```fortran
tise_relax=T
tise_relax_delta_t=1.e-3
tise_relax_max_iter=500
tise_relax_tol=1.e-10
```

At each iteration, the code applies the update:

```text
psi <- psi - delta_tau H psi
```

then renormalizes `psi`. After relaxation, the TISE solver is called with the obtained state as its initial state.

A value of `tise_relax_delta_t` that is too large can make relaxation unstable. Start with `1.e-3` or `1.e-4` for a new model.

### 9.4 Files produced by TISE

The solver writes :

- the file specified by `tise_energy_file`, containing the state index and energy;
- the file specified by `tise_transition_file`, containing the difference `E_i-E_0` and transition weight;
- `TISE_eigenvectors.dat`, containing the computed eigenvectors.

The weight written to the transition file is based on :

```text
| <psi_initial | psi_i> |^2
```

## 10. Complete calculation examples

In the examples below, `TDSI` is understood as `TISE`: this is the name used by the code for the time-independent Schrödinger equation.

### 10.1 TDSE example with a basis linked to the initial GWP

The file `test/default.inp` illustrates TDSE propagation of the Hénon-Heiles model. The basis is a product of two HO bases and one electronic basis:

```fortran
&basis_nd name='dp' nb_basis=3 /
&basis_nd name='HO' nb=20 nq=25 use_initial_gwp=T debug=F /
&basis_nd name='HO' nb=20 nq=25 use_initial_gwp=T debug=F /
&basis_nd name='el' nb=1 /

&defGWP ndim=2 Elecindex=1 Coef=(1.0,0.0) /
&defWP0 sigma=1.29 Beta=0.0 Qeq=2.0 imp_k=0.0 gamma=0.0 /
&defWP0 sigma=1.41 Beta=0.0 Qeq=0.0 imp_k=0.0 gamma=0.0 /

&prop
  t0=0.0 tf=60.0 delta_t=0.1
  eps=1.e-20 max_iter=500
  dynamics='TDSE'
  propa_name='hagedorn_taylor'
  Beta=T P=T renorm=F
  debug=F
/
```

The corresponding command is :

```bash
make FC=gfortran OMP=1 LAPACK=1 ARPACK=0 QML=1
bash test/default.inp 20
```

The script computes `nq=nb+5`, creates a results directory, and launches `TD_SCHROD.x`. The initial wavefunction defined in `&defGWP`/`&defWP0` is projected onto the HO/Hagedorn basis before propagation begins.

For propagation on a fixed basis, replace only the propagator choice with :

```fortran
propa_name='standard_taylor'
Beta=F P=F renorm=F debug=F
```

### 10.2 TISE example (stationary solution) with H2O

The ready-to-use file is :

```bash
bash test/TISE_H2O.inp
```

The essential part of its input is :

```fortran
&basis_nd name='dp' nb_basis=4 /
&basis_nd name='HO' nb=10 nq=12 use_initial_gwp=T debug=F /
&basis_nd name='HO' nb=10 nq=12 use_initial_gwp=T debug=F /
&basis_nd name='HO' nb=10 nq=12 use_initial_gwp=T debug=F /
&basis_nd name='el' nb=1 debug=F /

&defGWP ndim=3 Elecindex=1 Coef=(1.0,0.0) /
&defWP0 sigma=0.2645 Beta=0.0 Qeq=1.8108 imp_k=1.0 gamma=0.0 /
&defWP0 sigma=0.2645 Beta=0.0 Qeq=1.8108 imp_k=1.0 gamma=0.0 /
&defWP0 sigma=0.3390 Beta=0.0 Qeq=1.8231 imp_k=1.0 gamma=0.0 /

&prop
  t0=0.0 tf=0.0 delta_t=1.0
  eps=1.e-10 max_iter=500
  dynamics='TISE'
  tise_method='dense'
  tise_nstates=3
  tise_max_iter=20
  tise_tol=1.e-8
  energy_unit='cm-1'
  tise_energy_file='H2O_TISE_energies.dat'
  tise_transition_file='H2O_TISE_transitions.txt'
  tise_relax=F
  debug=F
/
```

In this example, the three HO bases have `use_initial_gwp=T`. The H2O QML parameters and the initial GWP parameters are therefore used to construct the basis before the stationary solution. `tise_method='dense'` is suitable for this small basis. The energy file, transition file, and `TISE_eigenvectors.dat` are written to the current directory.

### 10.3 TISE example with ARPACK and imaginary-time relaxation

For a larger basis, replace the block `&prop` by :

```fortran
&prop
  t0=0.0 tf=0.0 delta_t=1.0
  eps=1.e-10 max_iter=500
  dynamics='TISE'
  tise_method='arpack'
  tise_nstates=12
  tise_max_iter=500
  tise_tol=1.e-10
  energy_unit='cm-1'
  tise_relax=T
  tise_relax_delta_t=1.e-3
  tise_relax_max_iter=500
  tise_relax_tol=1.e-10
  debug=F
/
```

Compile and run this calculation :

```bash
make FC=gfortran OMP=1 LAPACK=1 ARPACK=1 QML=1
./TD_SCHROD.x < mon_TISE.inp > mon_TISE.log 2>&1
```

The relaxation first searches for an approximation to the ground state from the initial GWP. ARPACK then computes the requested `tise_nstates` eigenstates.

### 10.4 HONO TISE with ARPACK

```bash
bash test/TISE_HONO.inp
```

This example uses `nb=3`, `nq=4` for the bases HO, `use_initial_gwp=T`, `debug=F` and ARPACK to compute several states.

### 10.5 Botschwina ClH2+ TISE

The script accepts :

```bash
bash test/TISE_ClH2p_Botschwina.inp [option] [method] [nb] [nq]
```

Reference calculation example :

```bash
bash test/TISE_ClH2p_Botschwina.inp 2 arpack 12 12
```

- `option=1` : CEPA-1 ;
- `option=2` : corrected CEPA-1;
- `method` : `dense`, `lanczos`, `davidson` or `arpack` ;
- `nb` : number of HO functions per coordinate ;
- `nq` : number of quadrature points per coordinate.

For a quick test :

```bash
bash test/TISE_ClH2p_Botschwina.inp 2 arpack 2 4
```

The isotope variant 37Cl is available in :

```bash
test/TISE_ClH2p_Botschwina_37Cl.inp
```

## 11. Parameter-selection advice

- Start with a small basis (`nb` and `nq` modest) to check the syntax and potential.
- Increase `nb` before `nq` to study convergence of the vibrational states; the strategy nevertheless depends on the model.
- Use `dense` to validate a small basis and compare iterative methods.
- Use `arpack` or `lanczos` for large dimensions and a small number of states.
- Use `davidson` when several low-lying states are requested and a diagonal preconditioner is useful.
- In TDSE propagation, reduce `delta_t` if the Taylor propagator does not converge.
- Enable `debug=T` and `debug=T` in the `&basis_nd` blocks to diagnose basis construction.
- For an isotopologue, write the masses explicitly in the Z-matrix (`35Cl`, `37Cl`, `1H`, etc.).

## 12. Common-error troubleshooting

### `unsupported propagation method`

Check the spelling of `propa_name`. Recommended names include `standard_taylor`, `hagedorn_taylor`, `standard_rk4`, `standard_cn` and `standard_split`.

### `dynamics must be TDSE or TISE`

Use exactly :

```fortran
dynamics='TDSE'
```

or :

```fortran
dynamics='TISE'
```

### ARPACK not found during linking

The project installs ARPACK locally. Re-run the standard command:

```bash
make -j1 FC=gfortran OMP=1 LAPACK=1 ARPACK=1
```

If this fails before compilation, install the generic prerequisites `cmake`,
`curl` (or `wget`), `tar`, BLAS/LAPACK and a Fortran compiler. No
Homebrew-specific command is needed.

### `tise_max_iter` insufficient

Increase `tise_max_iter`. It must be at least equal to the number of requested states:

```fortran
tise_nstates=12
tise_max_iter=500
```

### Different results for two isotopologues

This is expected : the masses change the kinetic term. Compare only calculations using the same isotopes, the same potential option and the same basis.

## Validation and troubleshooting

For a new model, first run a small basis and short propagation. Check the initial norm, energy and grid-domain coverage, then increase `nb` and `nq` systematically. For CURV calculations, compare `construction_type='developed'` against `construction_type='divergence'` or `construction_type='similarity'` and inspect energy/norm convergence. For TISE, compare `dense` with an iterative solver on a small basis before scaling up. Use `write_validation=T` and `debug=T` when diagnosing norm, energy, basis-parameter, and autocorrelation behavior.

If the QML potential is evaluated at unexpected coordinates, verify `coord_type`, `sym_type`, and especially `pot_itQtransfo`; a non-negative stage transforms active coordinates before the QML call, whereas `pot_itQtransfo=-1` keeps them active. If a restart is refused, restore the original pre-`&prop` input and compatible integrator/numerical controls, use the exact `checkpoint_file`, and only then extend `tf`.

## 13. Useful directory tree

```text
app/TD_SCHROD.f90       Main program
src/Basis_m.f90         Basis reading and construction
src/Op_m.f90            Hamiltonian operator and application H|psi>
src/Propa_m.f90         Propagation loop
src/sub_propa_m.f90     `&prop` namelist and propagators
src/TISE_m.f90          Stationary solvers and imaginary-time relaxation
src/psi_m.f90           Wavefunctions and operations on psi
Makefile                Compilation and linking
EXT_Lib/                External libraries
 test/                  Examples and validation scripts
```

## 14. Quick summary

For a TDSE calculation :

```bash
make FC=gfortran OMP=1 LAPACK=1 ARPACK=0 QML=1
./TD_SCHROD.x < test/default.inp > tdse.log 2>&1
```

For a TISE calculation with ARPACK :

```bash
make FC=gfortran OMP=1 LAPACK=1 ARPACK=1 QML=1
bash test/TISE_ClH2p_Botschwina.inp 2 arpack 12 12
```

The basis is selected in the blocks `&basis_nd`, the dynamics is selected with `dynamics`, and the solver or propagator is selected in `&prop`.

## Contributing :
Any contribution is welcome, just contact me at issagoudouya@gmail.com, I'm open to collaboration and feedback.
