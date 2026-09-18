# Quantum Molecular Dynamics Simulator  
**Hagedorn Wavepacket Propagation in Fortran**  

A high-performance tool for simulating quantum molecular dynamics using generalized Hagedorn wavepackets with support for mixed primitive basis sets (e.g., Fourier, Hagedorn).  

---

## Table of Contents  
- [Key Features](#key-features)  
- [Theoretical Background](#theoretical-background)  
- [Installation](#installation)  
- [Default.comp file configuration](#configuration)
  - [Main program complilation](#Compilation)  
  - [Potential Setup](#potential-setup)  
  - [Basis Parameters](#basis-parameters)  
  - [Initial Wavepacket](#initial-wavepacket)  
  - [Propagation Settings](#propagation-settings)  
- [Testing](#Testing)  
- [License](#license)  

---

## Key Features  

### Advanced Wavepacket Propagation  
- Generalized Hagedorn basis (multi-dimensional)  
- Hybrid primitive bases (Hagedorn/Fourier)  
- Non-variational 3-step algorithm:  
  1. Standard propagation (asuming basis is time independent)  
  2. Basis parameter evolution  
  3. Projection onto updated basis  


### Customization Options  
- Tunable propagators (`Hagedorn`/`Taylor`)  
- Basis control flags (`B`, `P`, `renorm`)  
- Adiabatic/diabatic representations  

---

### Extended Hagedorn Formulation  
- Multi-dimensional systems  
- Mixed basis construction.

---

## configuration
The main program compilation need two directories one for object(obj) et second for external librairies.
They are created automatically.
The external libraries (stored in EXT_lib directory) are compilled automatically.

## Potential Setup (`&potential` Namelist)

### Basic Syntax
```fortran
&potential
    pot_name     ! Potential model name
    option       ! Model variant
    adiabatic    ! Adiabatic representation
    nsurf        ! Electronic surfaces
    ndim         ! Degrees of freedom
/
```

## Basis Parameters (`&basis_nd` Namelist)

### Global Basis Definition
```fortran
&basis_nd name='dp' nb_basis=3 /  ! Direct product of 3 primitive bases
```
###  Example of Basis Configuration
```fortran
! Fourier basis example (1D)
&basis_nd name='fourier' nb=64 nq=128 A=-5.0 B=5.0 scaleQ=1.0 /

! Hagedorn basis example (1D)
&basis_nd name='HO' nb=32 nq=64 Q0=0.0 Imp_k=1.0 Alpha=(1.0,0.1) /

! Electronic states basis
&basis_nd name='el' nb=1 /

```
---

### Basis Parameters Reference Table
| Parameter  | Type      | Description                                  | Required For          | Default Value |
|------------|-----------|----------------------------------------------|-----------------------|---------------|
| `name`     | string    | Basis type (`dp`,`fourier`,`hagedorn`,`el`)  | All                   | -             |
| `nb_basis` | integer   | Total number of primitive bases              | Global definition only| -             |
| `nb`       | integer   | Number of basis functions                    | Fourier/Hagedorn      | -             |
| `nq`       | integer   | Number of grid points                        | Fourier/Hagedorn      | -             |
| `A`, `B`   | real      | Grid boundaries [a.u.]                       | Fourier basis only    | -             |
| `Q0`       | real      | Initial center position [a.u.]               | Hagedorn basis        | 0.0           |
| `Imp_k`    | real      | Initial momentum [a.u.]                      | Optional              | 0.0           |
| `Alpha`    | complex   | Gaussian width (α = a + ib)                  | Hagedorn basis        | (1.0,0.0)     |
| `scaleQ`   | real      | Coordinate scaling factor                    | Optional              | 1.0           |
---

#### Legend:
- `dp` = Direct product basis
- `a.u.` = atomic units
- Complex numbers format: `(real_part,imaginary_part)`

#### Usage Notes:
1. **Global definition must come first**:
   ```fortran
   &basis_nd name='dp' nb_basis=3 /  ! Before any primitive basis

###  Complete 3D Example with 2 electronique surfaces.

```fortran
! Global definition
&basis_nd name='dp' nb_basis=4 /

! X-coordinate (Fourier)
&basis_nd name='fourier' nb=64 nq=128 A=-10.0 B=10.0 scaleQ=1.2 /

! Y-coordinate (Hagedorn)
&basis_nd name='hagedorn' nb=48 nq=96 Q0=1.5 Imp_k=0.5 Alpha=(0.8,0.2) /

! Z-coordinate (Hagedorn)
&basis_nd name='hagedorn' nb=32 nq=64 Q0=-2.0 Alpha=(1.2,0.0) /

! Electronic states (2 surfaces)
&basis_nd name='el' nb=2 /
```

## Initial Wavepacket Configuration

### 1. Global Wavepacket Parameters (`&defGWP`)

```fortran
&defGWP ndim=3 Elecindex=1 Coef=(1.0,0.0) /
```

| Parameter  | Type     | Description                          | Default    |
|------------|----------|--------------------------------------|------------|
| `ndim`     | integer  | Number of physical dimensions        | Required   |
| `Elecindex`| integer  | Initial electronic state (1=ground)  | 1          |
| `Coef`     | complex  | Wavepacket coefficient (real,imag)   | (1.0,0.0)  |

---

### 2. Per-Dimension Parameters (`&defWP0`)

#### Standard Syntax:
```fortran
&defWP0 sigma=0.5 Beta=0.1 Qeq=0.0 imp_k=0.0 gamma=0.0 /
```

| Parameter | Units  | Physical Meaning                  | Default |
|-----------|--------|-----------------------------------|---------|
| `sigma`   | a.u.   | Gaussian width (σ)                | 1.0     |
| `Beta`    | a.u.   | Imaginary part of width parameter | 0.0     |
| `Qeq`     | a.u.   | Initial position                  | 0.0     |
| `imp_k`   | a.u.   | Initial momentum                  | 0.0     |
| `gamma`   | rad    | Initial phase angle               | 0.0     |

---

## Propagation Settings (`&prop` Namelist)

### Basic Syntax
```fortran
&prop
  t0 tf delta_t
  propa_name propa_name2
  Beta P renorm
  eps max_iter
/
```
### Parameter Reference Table

| Parameter      | Type    | Description                          | Valid Options                 | Default     |
|----------------|---------|--------------------------------------|-------------------------------|-------------|
| `t0`           | real    | Initial time (a.u.)                  | ≥ 0.0                         | 0.0         |
| `tf`           | real    | Final time (a.u.)                    | > `t0`                        | Required    |
| `delta_t`      | real    | Time step (a.u.)                     | > 0.0                         | Required    |
| `propa_name`   | string  | Propagation type                     | `'hagedorn'`, `'no_hagedorn'` | Required    |
| `propa_name2`  | string  | Integrator method                    | `'Taylor'`, `'SIL'`, `'RK4'`  | `'Taylor'`  |
| `Beta`         | logical | Include imag(α) in basis evolution   | `.true.`, `.false.`           | `.false.`   |
| `P`            | logical | Include momentum in basis            | `.true.`, `.false.`           | `.false.`   |
| `renorm`       | logical | Renormalize wavefunction             | `.true.`, `.false.`           | `.false.`   |
| `eps`          | real    | Convergence threshold                | > 0.0                         | 1.0E-20     |
| `max_iter`     | integer | Max iterations per step              | ≥ 1                           | 500         |
---


## Compilation and Test

Compilation and propagation are handled together using the `default.inp` script.

```bash
# Build with a specific basis size (e.g., 10)
# Go to the test directory ./test and run
./default.inp 10
```

## License

See the `LICENSE` file for more information.