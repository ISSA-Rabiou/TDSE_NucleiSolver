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
  !!    * Hermitian positive-definite matrix square root
  !!    * Frobenius norms and admissibility diagnostics for (Q,P)
  !!    * generic matrix norms and admissibility diagnostics for (Q,P)
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
  public :: matrix_inverse_c, determinant_c, matrix_sqrt_hermitian
  public :: check_hagedorn_admissibility
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

    subroutine zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
      import :: Rkind
      character(len=1), intent(in) :: jobz, uplo
      integer, intent(in) :: n, lda, lwork
      complex(KIND=Rkind), intent(inout) :: a(lda,*)
      real(KIND=Rkind), intent(out) :: w(*)
      complex(KIND=Rkind), intent(inout) :: work(*)
      real(KIND=Rkind), intent(inout) :: rwork(*)
      integer, intent(out) :: info
    end subroutine zheev
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

  function matrix_sqrt_hermitian(A) result(S)
    complex(KIND=Rkind), intent(in) :: A(:,:)
    complex(KIND=Rkind), allocatable :: S(:,:)

    integer :: n, info, lwork, i
    complex(KIND=Rkind), allocatable :: U(:,:), work(:)
    real(KIND=Rkind), allocatable :: eig(:), rwork(:)
    complex(KIND=Rkind) :: work_query(1)

    n = size(A,1)
    if (size(A,2) /= n) error stop 'matrix_sqrt_hermitian: matrix must be square.'

    allocate(U(n,n), eig(n), rwork(max(1,3*n-2)))
    U = A

    lwork = -1
    call zheev('V','U',n,U,n,eig,work_query,lwork,rwork,info)
    if (info /= 0) error stop 'matrix_sqrt_hermitian: workspace query failed.'

    lwork = max(1, int(real(work_query(1), KIND=Rkind)))
    allocate(work(lwork))

    call zheev('V','U',n,U,n,eig,work,lwork,rwork,info)
    if (info /= 0) error stop 'matrix_sqrt_hermitian: diagonalization failed.'

    do i = 1, n
      if (eig(i) <= ZERO) error stop 'matrix_sqrt_hermitian: matrix not positive definite.'
    end do

    allocate(S(n,n))
    S = CZERO
    do i = 1, n
      S = S + sqrt(eig(i)) * outer_product(U(:,i), conjg(U(:,i)))
    end do

  contains
    function outer_product(x, y) result(mat)
      complex(KIND=Rkind), intent(in) :: x(:), y(:)
      complex(KIND=Rkind) :: mat(size(x), size(y))
      integer :: i, j
      do i = 1, size(x)
        do j = 1, size(y)
          mat(i,j) = x(i) * y(j)
        end do
      end do
    end function outer_product
  end function matrix_sqrt_hermitian


end module HagedornLinalg_m
