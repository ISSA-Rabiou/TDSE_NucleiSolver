module HagedornDiracFrenkel_m
  !!=========================================================================
  !!  HagedornDiracFrenkel_m
  !!
  !!  Variational propagation of a single homogeneous Hagedorn wavepacket
  !!  using the Dirac-Frenkel / variational thawed Gaussian equations in the
  !!  Hagedorn parametrization
  !!
  !!      dq/dt = M^{-1} p
  !!      dp/dt = - <grad V>
  !!      dQ/dt = M^{-1} P
  !!      dP/dt = - <Hess V> Q
  !!      dS/dt = T(p) - V0,
  !!
  !!  with
  !!      V0 = <V> - 1/2 Tr( <Hess V> Sigma )
  !!      Sigma = (eps^2 / 2) Q Q^H.
  !!
  !!  The module implements a symmetric TVT splitting, which is time-reversible
  !!  and follows the geometric-integration recommendations of the literature.
  !!
  !!  Scope of the present module:
  !!    * mathematically consistent propagation of one moving Hagedorn packet
  !!      (ground-state coefficient in the moving basis),
  !!    * expectation values and diagnostics,
  !!    * direct reuse of the existing Hagedorn basis implementation for
  !!      evaluation and exact kinetic moments.
  !!
  !!  References
  !!    * G. A. Hagedorn, Semiclassical quantum mechanics I--IV.
  !!    * E. Faou, V. Gradinaru, C. Lubich, Computing semiclassical quantum
  !!      dynamics with Hagedorn wavepackets, SIAM J. Sci. Comput. 31 (2009).
  !!    * J. J. L. Vanicek et al., Family of Gaussian wavepacket dynamics
  !!      methods from the perspective of a nonlinear Schrodinger equation,
  !!      J. Chem. Phys. 159, 014114 (2023).
  !!=========================================================================
  use QDUtil_m
  use HagedornLinalg_m
  use HagedornQuadrature_m, only: gauss_hermite_rule, tensor_product_grid
  use NDindex_hagedorn_m
  use HagedornBasis_m
  implicit none
  private

  public :: HagedornDFState_t
  public :: init_hagedorn_df_state
  public :: destroy_hagedorn_df_state
  public :: refresh_hagedorn_df_basis
  public :: propagate_hagedorn_df_tvt_step
  public :: compute_hagedorn_df_observables
  public :: compute_hagedorn_df_potential_moments
  public :: evaluate_hagedorn_df_on_grid
  public :: export_hagedorn_df_density_1d

  abstract interface
    subroutine scalar_potential_with_derivatives(x, v, gradV, hessV)
      import :: Rkind
      real(kind=Rkind), intent(in)  :: x(:)
      real(kind=Rkind), intent(out) :: v
      real(kind=Rkind), intent(out) :: gradV(:)
      real(kind=Rkind), intent(out) :: hessV(:,:)
    end subroutine scalar_potential_with_derivatives
  end interface

  interface
    subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
      import :: Rkind
      character(len=1), intent(in) :: jobz, uplo
      integer, intent(in) :: n, lda, lwork
      real(kind=Rkind), intent(inout) :: a(lda, *)
      real(kind=Rkind), intent(out) :: w(*)
      real(kind=Rkind), intent(inout) :: work(*)
      integer, intent(out) :: info
    end subroutine dsyev
  end interface

  type, public :: HagedornDFState_t
    integer :: ndim = 0
    real(kind=Rkind) :: eps = ONE
    real(kind=Rkind) :: time = ZERO
    type(HagedornBasisShape_t) :: shape
    type(HagedornBasis_t) :: basis
    complex(kind=Rkind), allocatable :: coeff(:)
    complex(kind=Rkind), allocatable :: q(:)
    complex(kind=Rkind), allocatable :: p(:)
    complex(kind=Rkind), allocatable :: qmat(:,:)
    complex(kind=Rkind), allocatable :: pmat(:,:)
    complex(kind=Rkind) :: actionS = CZERO
  end type HagedornDFState_t

contains

  subroutine init_hagedorn_df_state(state, q, p, qmat, pmat, actionS, eps)
    type(HagedornDFState_t), intent(inout) :: state
    complex(kind=Rkind), intent(in) :: q(:), p(:), qmat(:,:), pmat(:,:)
    complex(kind=Rkind), intent(in), optional :: actionS
    real(kind=Rkind), intent(in), optional :: eps

    call destroy_hagedorn_df_state(state)

    state%ndim = size(q)
    if (size(p) /= state%ndim) error stop 'init_hagedorn_df_state: incompatible q/p sizes.'
    if (size(qmat,1) /= state%ndim .or. size(qmat,2) /= state%ndim) error stop 'init_hagedorn_df_state: Q has wrong size.'
    if (size(pmat,1) /= state%ndim .or. size(pmat,2) /= state%ndim) error stop 'init_hagedorn_df_state: P has wrong size.'

    allocate(state%q(state%ndim), state%p(state%ndim), state%qmat(state%ndim,state%ndim), state%pmat(state%ndim,state%ndim))
    state%q = q
    state%p = p
    state%qmat = qmat
    state%pmat = pmat
    if (present(actionS)) state%actionS = actionS
    if (present(eps)) state%eps = eps

    call init_simplex_shape(state%shape, state%ndim, 0)
    allocate(state%coeff(1))
    state%coeff(1) = CONE

    call refresh_hagedorn_df_basis(state)
  end subroutine init_hagedorn_df_state


  subroutine destroy_hagedorn_df_state(state)
    type(HagedornDFState_t), intent(inout) :: state

    call destroy_hagedorn_basis(state%basis)
    call destroy_hagedorn_shape(state%shape)
    if (allocated(state%coeff)) deallocate(state%coeff)
    if (allocated(state%q)) deallocate(state%q)
    if (allocated(state%p)) deallocate(state%p)
    if (allocated(state%qmat)) deallocate(state%qmat)
    if (allocated(state%pmat)) deallocate(state%pmat)

    state%ndim = 0
    state%eps = ONE
    state%time = ZERO
    state%actionS = CZERO
  end subroutine destroy_hagedorn_df_state


  subroutine refresh_hagedorn_df_basis(state)
    type(HagedornDFState_t), intent(inout) :: state

    if (.not. allocated(state%q)) error stop 'refresh_hagedorn_df_basis: state not initialized.'
    call init_hagedorn_basis(state%basis, state%shape, state%q, state%p, state%qmat, state%pmat, state%actionS, state%eps)
    if (.not. state%basis%admissible) then
      error stop 'refresh_hagedorn_df_basis: Hagedorn admissibility lost.'
    end if
  end subroutine refresh_hagedorn_df_basis


  subroutine propagate_hagedorn_df_tvt_step(state, mass_inv, potential_cb, dt, quadrature_order)
    type(HagedornDFState_t), intent(inout) :: state
    real(kind=Rkind), intent(in) :: mass_inv(:,:)
    procedure(scalar_potential_with_derivatives) :: potential_cb
    real(kind=Rkind), intent(in) :: dt
    integer, intent(in), optional :: quadrature_order

    integer :: nq

    nq = max(8, 2*state%ndim + 6)
    if (present(quadrature_order)) nq = quadrature_order

    call kinetic_subflow(state, mass_inv, HALF*dt)
    call potential_subflow(state, potential_cb, nq, dt)
    call kinetic_subflow(state, mass_inv, HALF*dt)
    state%time = state%time + dt
    call refresh_hagedorn_df_basis(state)
  end subroutine propagate_hagedorn_df_tvt_step


  subroutine kinetic_subflow(state, mass_inv, dt)
    type(HagedornDFState_t), intent(inout) :: state
    real(kind=Rkind), intent(in) :: mass_inv(:,:)
    real(kind=Rkind), intent(in) :: dt
    complex(kind=Rkind) :: tkin

    if (size(mass_inv,1) /= state%ndim .or. size(mass_inv,2) /= state%ndim) then
      error stop 'kinetic_subflow: mass_inv has wrong size.'
    end if

    state%q    = state%q    + dt * matmul(mass_inv, state%p)
    state%qmat = state%qmat + dt * matmul(mass_inv, state%pmat)
    tkin = HALF * dot_product(state%p, matmul(mass_inv, state%p))
    state%actionS = state%actionS + dt * tkin
  end subroutine kinetic_subflow


  subroutine potential_subflow(state, potential_cb, quadrature_order, dt)
    type(HagedornDFState_t), intent(inout) :: state
    procedure(scalar_potential_with_derivatives) :: potential_cb
    integer, intent(in) :: quadrature_order
    real(kind=Rkind), intent(in) :: dt

    real(kind=Rkind) :: vavg
    real(kind=Rkind), allocatable :: gradavg(:), hessavg(:,:), sigma(:,:)
    complex(kind=Rkind) :: v0eff

    call compute_hagedorn_df_potential_moments(state, potential_cb, quadrature_order, vavg, gradavg, hessavg)
    sigma = covariance_matrix(state)

    state%p    = state%p    - dt * cmplx(gradavg, ZERO, kind=Rkind)
    state%pmat = state%pmat - dt * matmul(cmplx(hessavg, ZERO, kind=Rkind), state%qmat)

    v0eff = cmplx(vavg - HALF * trace_r(matmul(hessavg, sigma)), ZERO, kind=Rkind)
    state%actionS = state%actionS - dt * v0eff
  end subroutine potential_subflow


  subroutine compute_hagedorn_df_potential_moments(state, potential_cb, quadrature_order, vavg, gradavg, hessavg)
    type(HagedornDFState_t), intent(in) :: state
    procedure(scalar_potential_with_derivatives) :: potential_cb
    integer, intent(in) :: quadrature_order
    real(kind=Rkind), intent(out) :: vavg
    real(kind=Rkind), allocatable, intent(out) :: gradavg(:), hessavg(:,:)

    real(kind=Rkind), allocatable :: points(:,:), rho_w(:)
    real(kind=Rkind), allocatable :: gradV(:), hessV(:,:), x(:)
    real(kind=Rkind) :: vloc
    integer :: ipt, npts

    call build_density_quadrature(state, quadrature_order, points, rho_w)
    npts = size(points,2)

    allocate(gradavg(state%ndim), hessavg(state%ndim,state%ndim), gradV(state%ndim), hessV(state%ndim,state%ndim), x(state%ndim))
    gradavg = ZERO
    hessavg = ZERO
    vavg = ZERO

    do ipt = 1, npts
      x = points(:,ipt)
      call potential_cb(x, vloc, gradV, hessV)
      vavg = vavg + rho_w(ipt) * vloc
      gradavg = gradavg + rho_w(ipt) * gradV
      hessavg = hessavg + rho_w(ipt) * hessV
    end do
  end subroutine compute_hagedorn_df_potential_moments


  subroutine compute_hagedorn_df_observables(state, mass_inv, potential_cb, quadrature_order, norm, xavg, pavg, energy, kinetic, potential)
    type(HagedornDFState_t), intent(in) :: state
    real(kind=Rkind), intent(in) :: mass_inv(:,:)
    procedure(scalar_potential_with_derivatives) :: potential_cb
    integer, intent(in) :: quadrature_order
    real(kind=Rkind), intent(out) :: norm
    real(kind=Rkind), allocatable, intent(out) :: xavg(:), pavg(:)
    real(kind=Rkind), intent(out) :: energy, kinetic, potential

    type(HagedornBasis_t) :: ext_basis
    complex(kind=Rkind), allocatable :: coeff_embed(:), coeff_x(:,:), coeff_p(:,:)
    real(kind=Rkind) :: vavg
    real(kind=Rkind), allocatable :: gradavg(:), hessavg(:,:)
    complex(kind=Rkind) :: ztmp
    integer :: i, j

    if (size(mass_inv,1) /= state%ndim .or. size(mass_inv,2) /= state%ndim) then
      error stop 'compute_hagedorn_df_observables: mass_inv has wrong size.'
    end if

    norm = real(dot_product(state%coeff, state%coeff), kind=Rkind)

    call apply_position_operator(state%basis, state%coeff, ext_basis, coeff_x)
    call embed_coefficients(state%basis, ext_basis, state%coeff, coeff_embed)
    allocate(xavg(state%ndim))
    do i = 1, state%ndim
      xavg(i) = real(dot_product(coeff_embed, coeff_x(:,i)), kind=Rkind)
    end do
    call destroy_hagedorn_basis(ext_basis)
    if (allocated(coeff_x)) deallocate(coeff_x)
    if (allocated(coeff_embed)) deallocate(coeff_embed)

    call apply_momentum_operator(state%basis, state%coeff, ext_basis, coeff_p)
    call embed_coefficients(state%basis, ext_basis, state%coeff, coeff_embed)
    allocate(pavg(state%ndim))
    do i = 1, state%ndim
      pavg(i) = real(dot_product(coeff_embed, coeff_p(:,i)), kind=Rkind)
    end do

    ztmp = CZERO
    do j = 1, state%ndim
      do i = 1, state%ndim
        ztmp = ztmp + HALF * mass_inv(i,j) * dot_product(coeff_p(:,i), coeff_p(:,j))
      end do
    end do
    kinetic = real(ztmp, kind=Rkind)

    call compute_hagedorn_df_potential_moments(state, potential_cb, quadrature_order, vavg, gradavg, hessavg)
    potential = vavg
    energy = kinetic + potential
  end subroutine compute_hagedorn_df_observables


  subroutine evaluate_hagedorn_df_on_grid(state, points, psi)
    type(HagedornDFState_t), intent(in) :: state
    real(kind=Rkind), intent(in) :: points(:,:)
    complex(kind=Rkind), allocatable, intent(out) :: psi(:)

    call evaluate_wavepacket_on_grid(state%basis, state%coeff, points, psi)
  end subroutine evaluate_hagedorn_df_on_grid


  subroutine export_hagedorn_df_density_1d(state, xgrid, filename)
    type(HagedornDFState_t), intent(in) :: state
    real(kind=Rkind), intent(in) :: xgrid(:)
    character(len=*), intent(in) :: filename

    real(kind=Rkind), allocatable :: points(:,:)
    complex(kind=Rkind), allocatable :: psi(:)
    integer :: i, unitno

    if (state%ndim /= 1) error stop 'export_hagedorn_df_density_1d: only ndim=1 supported.'

    allocate(points(1,size(xgrid)))
    points(1,:) = xgrid
    call evaluate_hagedorn_df_on_grid(state, points, psi)

    open(newunit=unitno, file=filename, status='replace', action='write')
    write(unitno,'(A)') '# x   Re(psi)   Im(psi)   |psi|^2'
    do i = 1, size(xgrid)
      write(unitno,'(ES24.16,1X,ES24.16,1X,ES24.16,1X,ES24.16)') xgrid(i), real(psi(i),kind=Rkind), aimag(psi(i)), abs(psi(i))**2
    end do
    close(unitno)
  end subroutine export_hagedorn_df_density_1d


  subroutine build_density_quadrature(state, quadrature_order, points, rho_weights)
    type(HagedornDFState_t), intent(in) :: state
    integer, intent(in) :: quadrature_order
    real(kind=Rkind), allocatable, intent(out) :: points(:,:), rho_weights(:)

    real(kind=Rkind), allocatable :: nodes_1d(:), w_1d(:), xi(:,:), wprod(:), rootG(:,:), qreal(:)
    integer :: npts, ipt

    call gauss_hermite_rule(quadrature_order, nodes_1d, w_1d)
    call tensor_product_grid(nodes_1d, w_1d, state%ndim, xi, wprod)

    rootG = matrix_square_root_spd(real(matmul(state%qmat, conjg(transpose(state%qmat))), kind=Rkind))
    npts = size(xi,2)
    allocate(points(state%ndim,npts), rho_weights(npts), qreal(state%ndim))
    qreal = real(state%q, kind=Rkind)

    points = spread(qreal, 2, npts) + state%eps * matmul(rootG, xi)
    rho_weights = wprod / (PI**(HALF * real(state%ndim, kind=Rkind)))

    do ipt = 1, npts
      if (rho_weights(ipt) < ZERO) error stop 'build_density_quadrature: negative quadrature weight.'
    end do
  end subroutine build_density_quadrature


  function covariance_matrix(state) result(sigma)
    type(HagedornDFState_t), intent(in) :: state
    real(kind=Rkind), allocatable :: sigma(:,:)
    real(kind=Rkind), allocatable :: gmat(:,:)

    gmat = real(matmul(state%qmat, conjg(transpose(state%qmat))), kind=Rkind)
    allocate(sigma(state%ndim,state%ndim))
    sigma = HALF * state%eps**2 * symmetrize_r(gmat)
  end function covariance_matrix


  function matrix_square_root_spd(a) result(root)
    real(kind=Rkind), intent(in) :: a(:,:)
    real(kind=Rkind), allocatable :: root(:,:)
    real(kind=Rkind), allocatable :: evec(:,:), umat(:,:), eval(:), work(:)
    real(kind=Rkind) :: query(1)
    integer :: n, lwork, info, i

    n = size(a,1)
    if (size(a,2) /= n) error stop 'matrix_square_root_spd: matrix must be square.'

    allocate(evec(n,n), umat(n,n), eval(n))
    evec = symmetrize_r(a)

    lwork = -1
    call dsyev('V', 'U', n, evec, n, eval, query, lwork, info)
    if (info /= 0) error stop 'matrix_square_root_spd: DSYEV query failed.'
    lwork = max(1, int(query(1)))
    allocate(work(lwork))
    call dsyev('V', 'U', n, evec, n, eval, work, lwork, info)
    if (info /= 0) error stop 'matrix_square_root_spd: DSYEV failed.'
    if (minval(eval) <= ZERO) error stop 'matrix_square_root_spd: matrix is not positive definite.'

    umat = evec
    do i = 1, n
      evec(:,i) = evec(:,i) * sqrt(eval(i))
    end do

    allocate(root(n,n))
    root = matmul(evec, transpose(umat))
    root = symmetrize_r(root)
  end function matrix_square_root_spd


  pure function symmetrize_r(a) result(asym)
    real(kind=Rkind), intent(in) :: a(:,:)
    real(kind=Rkind) :: asym(size(a,1), size(a,2))
    asym = HALF * (a + transpose(a))
  end function symmetrize_r


  pure function trace_r(a) result(val)
    real(kind=Rkind), intent(in) :: a(:,:)
    real(kind=Rkind) :: val
    integer :: i, n

    n = min(size(a,1), size(a,2))
    val = ZERO
    do i = 1, n
      val = val + a(i,i)
    end do
  end function trace_r

end module HagedornDiracFrenkel_m
