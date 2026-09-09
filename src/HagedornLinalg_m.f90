module HagedornLinalg_m
  !!==========================================================================
  !!  HagedornLinalg_m
  !!
  !!  Low-level numerical utilities used by the Hagedorn basis library.
  !!  The module is intentionally self-contained and depends only on BLAS/LAPACK
  !!  interfaces that are available on standard scientific Fortran platforms.
  !!
  !!  Main services
  !!    * numerical kind and constants
  !!    * complex matrix inversion and determinant through LU factorization
  !!    * symmetric eigensolver wrapper for Gauss-Hermite quadrature
  !!    * Frobenius norms and admissibility diagnostics for (Q,P)
  !!    * tensor-product Gauss-Hermite quadrature grids
  !!
  !!  Mathematical conventions
  !!    The 1D Gauss-Hermite rule produced here integrates
  !!       ∫_{-∞}^{∞} exp(-x^2) f(x) dx
  !!    with the physicists' Hermite normalization.
  !!
  !!  References
  !!    * Golub & Welsch, Calculation of Gauss quadrature rules, Math. Comp. 1969.
  !!    * G. A. Hagedorn, Semiclassical quantum mechanics I--IV.
  !!==========================================================================
  USE QDUtil_m
  implicit none
  private

  public :: eye_c, eye_r
  public :: frob_norm_c, frob_norm_r
  public :: matrix_inverse_c, determinant_c
  public :: check_hagedorn_admissibility
  public :: gauss_hermite_rule
  public :: tensor_product_grid
  public :: max_abs_c

  interface
    subroutine zgetrf(m, n, a, lda, ipiv, info)
      import :: Rkind
      integer, intent(in) :: m, n, lda
      complex(Rkind), intent(inout) :: a(lda, *)
      integer, intent(out) :: ipiv(*)
      integer, intent(out) :: info
    end subroutine zgetrf

    subroutine zgetri(n, a, lda, ipiv, work, lwork, info)
      import :: Rkind
      integer, intent(in) :: n, lda, lwork
      complex(Rkind), intent(inout) :: a(lda, *)
      integer, intent(in) :: ipiv(*)
      complex(Rkind), intent(inout) :: work(*)
      integer, intent(out) :: info
    end subroutine zgetri

    subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
      import :: Rkind
      character(len=1), intent(in) :: jobz, uplo
      integer, intent(in) :: n, lda, lwork
      real(KIND=Rkind), intent(inout) :: a(lda, *)
      real(KIND=Rkind), intent(out) :: w(*)
      real(KIND=Rkind), intent(inout) :: work(*)
      integer, intent(out) :: info
    end subroutine dsyev
  end interface

contains

  pure function eye_c(n) result(iout)
    integer, intent(in) :: n
    complex(KIND=Rkind) :: iout(n,n)
    integer :: k
    iout =  CZERO
    do k = 1, n
      iout(k,k) = CONE
    end do
  end function eye_c

  pure function eye_r(n) result(iout)
    integer, intent(in) :: n
    real(KIND=Rkind) :: iout(n,n)
    integer :: k
    iout = ZERO
    do k = 1, n
      iout(k,k) = ONE
    end do
  end function eye_r

  pure function frob_norm_c(a) result(val)
    complex(KIND=Rkind), intent(in) :: a(:,:)
    real(KIND=Rkind) :: val
    val = sqrt(sum(abs(a)**2))
  end function frob_norm_c

  pure function frob_norm_r(a) result(val)
    real(KIND=Rkind), intent(in) :: a(:,:)
    real(KIND=Rkind) :: val
    val = sqrt(sum(a**2))
  end function frob_norm_r

  pure function max_abs_c(a) result(val)
    complex(KIND=Rkind), intent(in) :: a(:,:)
    real(KIND=Rkind) :: val
    if (size(a) == 0) then
      val = ZERO
    else
      val = maxval(abs(a))
    end if
  end function max_abs_c

  function matrix_inverse_c(a) result(ainv)
    complex(KIND=Rkind), intent(in) :: a(:,:)
    complex(KIND=Rkind), allocatable :: ainv(:,:)
    integer :: n, info, lwork
    integer, allocatable :: ipiv(:)
    complex(KIND=Rkind), allocatable :: work(:)
    complex(KIND=Rkind) :: query(1)

    n = size(a,1)
    if (size(a,2) /= n) error stop 'matrix_inverse_c: matrix must be square.'

    allocate(ainv(n,n), ipiv(n))
    ainv = a

    call zgetrf(n, n, ainv, n, ipiv, info)
    if (info /= 0) error stop 'matrix_inverse_c: LU factorization failed.'

    lwork = -1
    call zgetri(n, ainv, n, ipiv, query, lwork, info)
    lwork = max(1, int(real(query(1), KIND=Rkind)))
    allocate(work(lwork))

    call zgetri(n, ainv, n, ipiv, work, lwork, info)
    if (info /= 0) error stop 'matrix_inverse_c: matrix inversion failed.'
  end function matrix_inverse_c

  function determinant_c(a) result(det_a)
    complex(KIND=Rkind), intent(in) :: a(:,:)
    complex(KIND=Rkind) :: det_a
    complex(KIND=Rkind), allocatable :: lu(:,:)
    integer, allocatable :: ipiv(:)
    integer :: n, info, k, sign

    n = size(a,1)
    if (size(a,2) /= n) error stop 'determinant_c: matrix must be square.'

    allocate(lu(n,n), ipiv(n))
    lu = a
    call zgetrf(n, n, lu, n, ipiv, info)
    if (info /= 0) error stop 'determinant_c: LU factorization failed.'

    sign = 1
    do k = 1, n
      if (ipiv(k) /= k) sign = -sign
    end do

    det_a = cmplx(real(sign,KIND=Rkind), ZERO, KIND=Rkind)
    do k = 1, n
      det_a = det_a * lu(k,k)
    end do
  end function determinant_c

  subroutine check_hagedorn_admissibility(qmat, pmat, ok, residual_symp, residual_norm, tol)
    !!---------------------------------------------------------------------
    !! Verify the Hagedorn admissibility conditions
    !!    Q^T P - P^T Q = 0,
    !!    Q^H P - P^H Q = 2 i I.
    !!---------------------------------------------------------------------
    complex(KIND=Rkind), intent(in) :: qmat(:,:), pmat(:,:)
    logical, intent(out) :: ok
    real(KIND=Rkind), intent(out) :: residual_symp, residual_norm
    real(KIND=Rkind), intent(in), optional :: tol
    complex(KIND=Rkind), allocatable :: r1(:,:), r2(:,:)
    real(KIND=Rkind) :: threshold
    integer :: d

    if (size(qmat,1) /= size(qmat,2)) error stop 'check_hagedorn_admissibility: Q must be square.'
    if (size(pmat,1) /= size(pmat,2)) error stop 'check_hagedorn_admissibility: P must be square.'
    if (size(qmat,1) /= size(pmat,1)) error stop 'check_hagedorn_admissibility: Q and P dimensions mismatch.'

    d = size(qmat,1)
    threshold = ONETENTH**10
    if (present(tol)) threshold = tol

    allocate(r1(d,d), r2(d,d))
    r1 = matmul(transpose(qmat), pmat) - matmul(transpose(pmat), qmat)
    r2 = matmul(conjg(transpose(qmat)), pmat) - matmul(conjg(transpose(pmat)), qmat) - TWO*EYE*eye_c(d)

    residual_symp = frob_norm_c(r1)
    residual_norm = frob_norm_c(r2)
    ok = (residual_symp <= threshold) .and. (residual_norm <= threshold)
  end subroutine check_hagedorn_admissibility

  subroutine gauss_hermite_rule(n, x, w)
    !!---------------------------------------------------------------------
    !! Compute the n-point Gauss-Hermite quadrature rule for the weight
    !! exp(-x^2) by diagonalizing the Jacobi matrix of Hermite polynomials.
    !!---------------------------------------------------------------------
    integer, intent(in) :: n
    real(KIND=Rkind), allocatable, intent(out) :: x(:), w(:)
    real(KIND=Rkind), allocatable :: jmat(:,:), evals(:), work(:)
    real(KIND=Rkind) :: query(1)
    integer :: info, lwork, k

    if (n <= 0) error stop 'gauss_hermite_rule: n must be positive.'

    allocate(jmat(n,n), evals(n), x(n), w(n))
    jmat = ZERO
    do k = 1, n-1
      jmat(k, k+1) = sqrt(real(k,Rkind) / 2.0_Rkind)
      jmat(k+1, k) = jmat(k, k+1)
    end do

    lwork = -1
    call dsyev('V', 'U', n, jmat, n, evals, query, lwork, info)
    lwork = max(1, int(query(1)))
    allocate(work(lwork))
    call dsyev('V', 'U', n, jmat, n, evals, work, lwork, info)
    if (info /= 0) error stop 'gauss_hermite_rule: DSYEV failed.'

    x = evals
    do k = 1, n
      w(k) = sqrt(PI) * jmat(1,k)**2
    end do
  end subroutine gauss_hermite_rule

  subroutine tensor_product_grid(nodes_1d, weights_1d, ndim, points, weights)
    !!---------------------------------------------------------------------
    !! Build the tensor-product grid from a common 1D rule.
    !!
    !! Input:
    !!   nodes_1d, weights_1d : 1D quadrature nodes and weights.
    !!   ndim                 : physical dimension.
    !!
    !! Output:
    !!   points(ndim, npts)   : tensor-product grid.
    !!   weights(npts)        : product weights.
    !!---------------------------------------------------------------------
    real(KIND=Rkind), intent(in) :: nodes_1d(:), weights_1d(:)
    integer, intent(in) :: ndim
    real(KIND=Rkind), allocatable, intent(out) :: points(:,:), weights(:)

    integer :: n1d, npts, ipt
    integer, allocatable :: idx(:)
    logical :: done

    if (size(nodes_1d) /= size(weights_1d)) error stop 'tensor_product_grid: inconsistent 1D rule.'
    if (ndim <= 0) error stop 'tensor_product_grid: ndim must be positive.'

    n1d = size(nodes_1d)
    npts = n1d**ndim
    allocate(points(ndim, npts), weights(npts), idx(ndim))

    idx = 1
    ipt = 0
    done = .false.
    do while (.not. done)
      ipt = ipt + 1
      weights(ipt) = ONE
      points(:,ipt) = ZERO
      points(:,ipt) = nodes_1d(idx)
      weights(ipt) = product(weights_1d(idx))

      call increment_multi_index(idx, n1d, done)
    end do
  contains
    subroutine increment_multi_index(ind, upper, finished)
      integer, intent(inout) :: ind(:)
      integer, intent(in) :: upper
      logical, intent(out) :: finished
      integer :: d

      finished = .false.
      ind(1) = ind(1) + 1
      do d = 1, size(ind)-1
        if (ind(d) > upper) then
          ind(d) = 1
          ind(d+1) = ind(d+1) + 1
        end if
      end do
      if (ind(size(ind)) > upper) finished = .true.
    end subroutine increment_multi_index
  end subroutine tensor_product_grid

end module HagedornLinalg_m
