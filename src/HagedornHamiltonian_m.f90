module HagedornHamiltonian_m
  !!==========================================================================
  !!  HagedornHamiltonian_m
  !!
  !!  Analytic Hamiltonian action in the Hagedorn basis.
  !!
  !!  Provided capabilities
  !!    * kinetic operator from exact second derivatives,
  !!    * exact action of quadratic scalar potentials,
  !!    * Galerkin projection of a general grid potential back to a chosen
  !!      Hagedorn basis, including multi-surface matrix-valued potentials.
  !!
  !!  This design mirrors the role of Op_m.f90 in TDSE_NucleiSolver but replaces
  !!  finite-difference style thinking by exact Hagedorn operator algebra.
  !!==========================================================================
  USE QDUtil_m
  use HagedornLinalg_m
  use HagedornBasis_m
  implicit none
  private

  public :: apply_kinetic_operator
  public :: apply_quadratic_potential
  public :: apply_hamiltonian_quadratic
  public :: build_hamiltonian_matrix_quadratic
  public :: apply_hamiltonian_grid_projection
  public :: apply_hamiltonian_grid_projection_multi

contains

  subroutine apply_kinetic_operator(basis, coeff, mass_inv, ext2_basis, coeff_t)
    !!---------------------------------------------------------------------
    !! T ψ = -1/2 ε⁴ Σ_{i,j} M^{-1}_{ij} ∂_i ∂_j ψ
    !!---------------------------------------------------------------------
    type(HagedornBasis_t), intent(in) :: basis
    complex(KIND=Rkind), intent(in) :: coeff(:)
    real(KIND=Rkind), intent(in) :: mass_inv(:,:)
    type(HagedornBasis_t), intent(inout) :: ext2_basis
    complex(KIND=Rkind), allocatable, intent(out) :: coeff_t(:)

    complex(KIND=Rkind), allocatable :: coeff_hess(:,:,:)
    integer :: i, j

    if (size(mass_inv,1) /= basis%ndim .or. size(mass_inv,2) /= basis%ndim) then
      error stop 'apply_kinetic_operator: mass_inv has wrong size.'
    end if

    call apply_second_derivative_operator(basis, coeff, ext2_basis, coeff_hess)
    allocate(coeff_t(ext2_basis%basis_size))
    coeff_t = CZERO

    do j = 1, basis%ndim
      do i = 1, basis%ndim
        coeff_t = coeff_t - HALF * basis%eps**4 * mass_inv(i,j) * coeff_hess(:,i,j)
      end do
    end do
  end subroutine apply_kinetic_operator

  subroutine apply_quadratic_potential(basis, coeff, v0, gradV, hessV, ext2_basis, coeff_v)
    !!---------------------------------------------------------------------
    !! Apply the scalar quadratic potential
    !!    V(x) = v0 + gradV^T x + 1/2 x^T hessV x
    !! exactly in coefficient space.
    !!---------------------------------------------------------------------
    type(HagedornBasis_t), intent(in) :: basis
    complex(KIND=Rkind), intent(in) :: coeff(:)
    complex(KIND=Rkind), intent(in) :: v0
    complex(KIND=Rkind), intent(in) :: gradV(:)
    complex(KIND=Rkind), intent(in) :: hessV(:,:)
    type(HagedornBasis_t), intent(inout) :: ext2_basis
    complex(KIND=Rkind), allocatable, intent(out) :: coeff_v(:)

    type(HagedornBasis_t) :: ext1_basis, tmp_basis
    complex(KIND=Rkind), allocatable :: coeff_embed(:), coeff_x(:,:), coeff_x_embed(:)
    complex(KIND=Rkind), allocatable :: coeff_xx(:,:)
    integer :: i, j

    if (size(gradV) /= basis%ndim) error stop 'apply_quadratic_potential: gradV has wrong size.'
    if (size(hessV,1) /= basis%ndim .or. size(hessV,2) /= basis%ndim) then
      error stop 'apply_quadratic_potential: hessV has wrong size.'
    end if

    call apply_position_operator(basis, coeff, ext1_basis, coeff_x)
    call init_extended_hagedorn_basis(ext2_basis, ext1_basis, 1)
    allocate(coeff_v(ext2_basis%basis_size))
    coeff_v = CZERO

    call embed_coefficients(basis, ext2_basis, coeff, coeff_embed)
    coeff_v = coeff_v + v0 * coeff_embed

    do i = 1, basis%ndim
      call embed_coefficients(ext1_basis, ext2_basis, coeff_x(:,i), coeff_x_embed)
      coeff_v = coeff_v + gradV(i) * coeff_x_embed
      deallocate(coeff_x_embed)
    end do

    do j = 1, basis%ndim
      call apply_position_operator(ext1_basis, coeff_x(:,j), tmp_basis, coeff_xx)
      do i = 1, basis%ndim
        coeff_v = coeff_v + HALF * hessV(i,j) * coeff_xx(:,i)
      end do
      call destroy_hagedorn_basis(tmp_basis)
      if (allocated(coeff_xx)) deallocate(coeff_xx)
    end do
  end subroutine apply_quadratic_potential

  subroutine apply_hamiltonian_quadratic(basis, coeff, mass_inv, v0, gradV, hessV, target_basis, coeff_h)
    !!---------------------------------------------------------------------
    !! Full scalar Hamiltonian H = T + V for a quadratic potential, with the
    !! result projected onto target_basis.
    !!---------------------------------------------------------------------
    type(HagedornBasis_t), intent(in) :: basis, target_basis
    complex(KIND=Rkind), intent(in) :: coeff(:)
    real(KIND=Rkind), intent(in) :: mass_inv(:,:)
    complex(KIND=Rkind), intent(in) :: v0
    complex(KIND=Rkind), intent(in) :: gradV(:)
    complex(KIND=Rkind), intent(in) :: hessV(:,:)
    complex(KIND=Rkind), allocatable, intent(out) :: coeff_h(:)

    type(HagedornBasis_t) :: ext2_t, ext2_v
    complex(KIND=Rkind), allocatable :: coeff_t(:), coeff_v(:), coeff_t_proj(:), coeff_v_proj(:)

    call apply_kinetic_operator(basis, coeff, mass_inv, ext2_t, coeff_t)
    call apply_quadratic_potential(basis, coeff, v0, gradV, hessV, ext2_v, coeff_v)

    call project_coefficients(ext2_t, target_basis, coeff_t, coeff_t_proj)
    call project_coefficients(ext2_v, target_basis, coeff_v, coeff_v_proj)

    allocate(coeff_h(target_basis%basis_size))
    coeff_h = coeff_t_proj + coeff_v_proj
  end subroutine apply_hamiltonian_quadratic

  subroutine build_hamiltonian_matrix_quadratic(basis, mass_inv, v0, gradV, hessV, hmat)
    type(HagedornBasis_t), intent(in) :: basis
    real(KIND=Rkind), intent(in) :: mass_inv(:,:)
    complex(KIND=Rkind), intent(in) :: v0
    complex(KIND=Rkind), intent(in) :: gradV(:)
    complex(KIND=Rkind), intent(in) :: hessV(:,:)
    complex(KIND=Rkind), allocatable, intent(out) :: hmat(:,:)

    complex(KIND=Rkind), allocatable :: ej(:), hj(:)
    integer :: j

    allocate(hmat(basis%basis_size, basis%basis_size))
    hmat = CZERO
    allocate(ej(basis%basis_size))

    do j = 1, basis%basis_size
      ej = CZERO
      ej(j) = CONE
      call apply_hamiltonian_quadratic(basis, ej, mass_inv, v0, gradV, hessV, basis, hj)
      hmat(:,j) = hj
      if (allocated(hj)) deallocate(hj)
    end do
  end subroutine build_hamiltonian_matrix_quadratic

  subroutine apply_hamiltonian_grid_projection(basis, coeff, mass_inv, points, weights, potential_values, target_basis, coeff_h)
    !!---------------------------------------------------------------------
    !! Scalar Hamiltonian with arbitrary scalar potential sampled on a grid.
    !!
    !! Kinetic part: analytic in coefficient space.
    !! Potential part: evaluate ψ on grid, multiply by V(x), then project by
    !! quadrature onto target_basis.
    !!---------------------------------------------------------------------
    type(HagedornBasis_t), intent(in) :: basis, target_basis
    complex(KIND=Rkind), intent(in) :: coeff(:)
    real(KIND=Rkind), intent(in) :: mass_inv(:,:)
    real(KIND=Rkind), intent(in) :: points(:,:), weights(:)
    complex(KIND=Rkind), intent(in) :: potential_values(:)
    complex(KIND=Rkind), allocatable, intent(out) :: coeff_h(:)

    type(HagedornBasis_t) :: ext2_basis
    complex(KIND=Rkind), allocatable :: coeff_t(:), coeff_t_proj(:)
    complex(KIND=Rkind), allocatable :: phi_target(:,:), psi_vals(:), hpsi_vals(:)
    integer :: npts, ib

    npts = size(points,2)
    if (size(weights) /= npts) error stop 'apply_hamiltonian_grid_projection: weights size mismatch.'
    if (size(potential_values) /= npts) error stop 'apply_hamiltonian_grid_projection: potential size mismatch.'

    call apply_kinetic_operator(basis, coeff, mass_inv, ext2_basis, coeff_t)
    call project_coefficients(ext2_basis, target_basis, coeff_t, coeff_t_proj)

    call evaluate_wavepacket_on_grid(basis, coeff, points, psi_vals)
    hpsi_vals = potential_values * psi_vals

    call evaluate_basis_on_grid(target_basis, points, phi_target)

    allocate(coeff_h(target_basis%basis_size))
    coeff_h = coeff_t_proj
    do ib = 1, target_basis%basis_size
      coeff_h(ib) = coeff_h(ib) + sum(conjg(phi_target(ib,:)) * weights * hpsi_vals)
    end do
  end subroutine apply_hamiltonian_grid_projection

  subroutine apply_hamiltonian_grid_projection_multi(basis, coeff, mass_inv, points, weights, potential_values, target_basis, coeff_h)
    !!---------------------------------------------------------------------
    !! Multi-surface analogue of apply_hamiltonian_grid_projection.
    !!
    !! Input:
    !!   coeff(nb, nsurf)
    !!   potential_values(npts, nsurf, nsurf)
    !! Output:
    !!   coeff_h(ntarget, nsurf)
    !!---------------------------------------------------------------------
    type(HagedornBasis_t), intent(in) :: basis, target_basis
    complex(KIND=Rkind), intent(in) :: coeff(:,:)
    real(KIND=Rkind), intent(in) :: mass_inv(:,:)
    real(KIND=Rkind), intent(in) :: points(:,:), weights(:)
    complex(KIND=Rkind), intent(in) :: potential_values(:,:,:)
    complex(KIND=Rkind), allocatable, intent(out) :: coeff_h(:,:)

    type(HagedornBasis_t) :: ext2_basis
    complex(KIND=Rkind), allocatable :: coeff_t(:), coeff_t_proj(:)
    complex(KIND=Rkind), allocatable :: phi_target(:,:), psi_vals(:,:), vpsi(:,:), psi_tmp(:)
    integer :: npts, nsurf, ib, is, js

    if (size(coeff,1) /= basis%basis_size) error stop 'apply_hamiltonian_grid_projection_multi: bad coefficient size.'
    npts = size(points,2)
    nsurf = size(coeff,2)
    if (size(weights) /= npts) error stop 'apply_hamiltonian_grid_projection_multi: weights size mismatch.'
    if (size(potential_values,1) /= npts .or. size(potential_values,2) /= nsurf .or. size(potential_values,3) /= nsurf) then
      error stop 'apply_hamiltonian_grid_projection_multi: potential tensor has wrong size.'
    end if

    allocate(coeff_h(target_basis%basis_size, nsurf))
    coeff_h = CZERO

    ! Kinetic term acts independently on each electronic surface.
    do is = 1, nsurf
      call apply_kinetic_operator(basis, coeff(:,is), mass_inv, ext2_basis, coeff_t)
      call project_coefficients(ext2_basis, target_basis, coeff_t, coeff_t_proj)
      coeff_h(:,is) = coeff_h(:,is) + coeff_t_proj
      if (allocated(coeff_t)) deallocate(coeff_t)
      if (allocated(coeff_t_proj)) deallocate(coeff_t_proj)
    end do

    allocate(psi_vals(npts, nsurf), vpsi(npts, nsurf))
    psi_vals = CZERO
    vpsi = CZERO

    do is = 1, nsurf
      call evaluate_wavepacket_on_grid(basis, coeff(:,is), points, psi_tmp)
      psi_vals(:,is) = psi_tmp
      if (allocated(psi_tmp)) deallocate(psi_tmp)
    end do

    do is = 1, nsurf
      do js = 1, nsurf
        vpsi(:,is) = vpsi(:,is) + potential_values(:,is,js) * psi_vals(:,js)
      end do
    end do

    call evaluate_basis_on_grid(target_basis, points, phi_target)
    do is = 1, nsurf
      do ib = 1, target_basis%basis_size
        coeff_h(ib,is) = coeff_h(ib,is) + sum(conjg(phi_target(ib,:)) * weights * vpsi(:,is))
      end do
    end do
  end subroutine apply_hamiltonian_grid_projection_multi

end module HagedornHamiltonian_m
