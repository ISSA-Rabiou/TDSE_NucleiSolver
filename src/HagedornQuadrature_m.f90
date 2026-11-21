module HagedornQuadrature_m
  !! Gauss-Hermite rules and tensor-product grids for Hagedorn quadrature.
  use QDUtil_m
  implicit none
  private
  public :: gauss_hermite_rule, tensor_product_grid

  interface
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


end module HagedornQuadrature_m
