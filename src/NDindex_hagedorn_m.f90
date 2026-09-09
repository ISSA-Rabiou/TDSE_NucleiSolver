module NDindex_hagedorn_m
   USE QDUtil_m
  implicit none
  private

  integer, parameter, public :: SHAPE_HYPERCUBIC = 1
  integer, parameter, public :: SHAPE_SIMPLEX    = 2
  integer, parameter, public :: SHAPE_HYPERBOLIC = 3

  type, public :: HagedornBasisShape_t
    integer :: ndim = 0
    integer :: shape_kind = 0
    integer, allocatable :: limits(:)  ! Hypercubic: k_d < limits(d)
    integer :: nmax = -1               ! Simplex: sum_d k_d <= nmax
    integer :: k_hyp = -1              ! Hyperbolic: prod_d (1+k_d) <= k_hyp
    integer(KIND=ILkind) :: basis_size = 0_ILkind
    integer, allocatable :: nodes(:,:) ! (ndim, basis_size)
  end type HagedornBasisShape_t

  public :: init_hypercubic_shape
  public :: init_simplex_shape
  public :: init_hyperbolic_shape
  public :: destroy_hagedorn_shape
  public :: get_basis_size
  public :: get_multiindex
  public :: find_multiindex
  public :: is_in_basis
  public :: build_hypercubic_indices
  public :: build_simplex_indices
  public :: build_hyperbolic_indices

contains

  subroutine init_hypercubic_shape(shape, limits)
    type(HagedornBasisShape_t), intent(inout) :: shape
    integer, intent(in) :: limits(:)

    if (size(limits) <= 0) error stop 'init_hypercubic_shape: empty limits array.'
    if (any(limits <= 0)) error stop 'init_hypercubic_shape: all limits must be positive.'

    call destroy_hagedorn_shape(shape)

    shape%ndim       = size(limits)
    shape%shape_kind = SHAPE_HYPERCUBIC
    allocate(shape%limits(shape%ndim))
    shape%limits = limits
    shape%nmax  = -1
    shape%k_hyp = -1

    call build_hypercubic_indices(limits, shape%nodes)
    shape%basis_size = int(size(shape%nodes, 2), ILkind)
  end subroutine init_hypercubic_shape


  subroutine init_simplex_shape(shape, ndim, nmax)
    type(HagedornBasisShape_t), intent(inout) :: shape
    integer, intent(in) :: ndim, nmax

    if (ndim <= 0) error stop 'init_simplex_shape: ndim must be > 0.'
    if (nmax < 0)  error stop 'init_simplex_shape: nmax must be >= 0.'

    call destroy_hagedorn_shape(shape)

    shape%ndim       = ndim
    shape%shape_kind = SHAPE_SIMPLEX
    allocate(shape%limits(ndim))
    shape%limits = nmax + 1
    shape%nmax  = nmax
    shape%k_hyp = -1

    call build_simplex_indices(ndim, nmax, shape%nodes)
    shape%basis_size = int(size(shape%nodes, 2), ILkind)
  end subroutine init_simplex_shape


  subroutine init_hyperbolic_shape(shape, ndim, k_hyp)
    type(HagedornBasisShape_t), intent(inout) :: shape
    integer, intent(in) :: ndim, k_hyp

    if (ndim <= 0) error stop 'init_hyperbolic_shape: ndim must be > 0.'
    if (k_hyp <= 0) error stop 'init_hyperbolic_shape: k_hyp must be > 0.'

    call destroy_hagedorn_shape(shape)

    shape%ndim       = ndim
    shape%shape_kind = SHAPE_HYPERBOLIC
    allocate(shape%limits(ndim))
    shape%limits = k_hyp
    shape%nmax  = -1
    shape%k_hyp = k_hyp

    call build_hyperbolic_indices(ndim, k_hyp, shape%nodes)
    shape%basis_size = int(size(shape%nodes, 2), ILkind)
  end subroutine init_hyperbolic_shape


  subroutine destroy_hagedorn_shape(shape)
    type(HagedornBasisShape_t), intent(inout) :: shape

    if (allocated(shape%limits)) deallocate(shape%limits)
    if (allocated(shape%nodes))  deallocate(shape%nodes)
    shape%ndim       = 0
    shape%shape_kind = 0
    shape%nmax       = -1
    shape%k_hyp      = -1
    shape%basis_size = 0_ILkind
  end subroutine destroy_hagedorn_shape


  integer(KIND=ILkind) function get_basis_size(shape) result(n)
    type(HagedornBasisShape_t), intent(in) :: shape
    n = shape%basis_size
  end function get_basis_size


  subroutine get_multiindex(shape, i, k)
    type(HagedornBasisShape_t), intent(in) :: shape
    integer, intent(in) :: i
    integer, intent(out) :: k(:)

    if (.not. allocated(shape%nodes)) error stop 'get_multiindex: shape was not initialized.'
    if (size(k) /= shape%ndim) error stop 'get_multiindex: wrong output size.'
    if (i < 1 .or. i > size(shape%nodes, 2)) error stop 'get_multiindex: invalid basis index.'

    k = shape%nodes(:, i)
  end subroutine get_multiindex


  integer function find_multiindex(shape, k) result(idx)
    type(HagedornBasisShape_t), intent(in) :: shape
    integer, intent(in) :: k(:)
    integer :: j

    idx = 0
    if (.not. allocated(shape%nodes)) return
    if (size(k) /= shape%ndim) return

    do j = 1, size(shape%nodes, 2)
      if (all(shape%nodes(:, j) == k)) then
        idx = j
        return
      end if
    end do
  end function find_multiindex


  logical function is_in_basis(shape, k) result(ok)
    type(HagedornBasisShape_t), intent(in) :: shape
    integer, intent(in) :: k(:)

    ok = (find_multiindex(shape, k) > 0)
  end function is_in_basis


  subroutine build_hypercubic_indices(limits, nodes)
    integer, intent(in) :: limits(:)
    integer, allocatable, intent(out) :: nodes(:,:)

    integer :: ndim
    integer(KIND=ILkind) :: n_nodes, pos
    integer, allocatable :: current(:)

    ndim = size(limits)
    if (ndim <= 0) error stop 'build_hypercubic_indices: ndim must be > 0.'
    if (any(limits <= 0)) error stop 'build_hypercubic_indices: all limits must be positive.'

    n_nodes = count_hypercubic_indices(limits)
    allocate(nodes(ndim, int(n_nodes, ILkind)))
    allocate(current(ndim))
    current = 0
    pos = 0_ILkind

    call rec_fill_hypercubic(ndim, limits, current, nodes, pos)
  end subroutine build_hypercubic_indices


  subroutine build_simplex_indices(ndim, nmax, nodes)
    integer, intent(in) :: ndim, nmax
    integer, allocatable, intent(out) :: nodes(:,:)

    integer(KIND=ILkind) :: n_nodes, pos
    integer, allocatable :: current(:)

    if (ndim <= 0) error stop 'build_simplex_indices: ndim must be > 0.'
    if (nmax < 0)  error stop 'build_simplex_indices: nmax must be >= 0.'

    n_nodes = count_simplex_indices(ndim, nmax)
    allocate(nodes(ndim, int(n_nodes)))
    allocate(current(ndim))
    current = 0
    pos = 0_ILkind

    call rec_fill_simplex(ndim, nmax, current, nodes, pos)
  end subroutine build_simplex_indices


  subroutine build_hyperbolic_indices(ndim, k_hyp, nodes)
    integer, intent(in) :: ndim, k_hyp
    integer, allocatable, intent(out) :: nodes(:,:)

    integer(KIND=ILkind) :: n_nodes, pos
    integer, allocatable :: current(:)

    if (ndim <= 0) error stop 'build_hyperbolic_indices: ndim must be > 0.'
    if (k_hyp <= 0) error stop 'build_hyperbolic_indices: k_hyp must be > 0.'

    n_nodes = count_hyperbolic_indices(ndim, k_hyp)
    allocate(nodes(ndim, int(n_nodes)))
    allocate(current(ndim))
    current = 0
    pos = 0_ILkind

    call rec_fill_hyperbolic(ndim, k_hyp, current, nodes, pos)
  end subroutine build_hyperbolic_indices


  recursive integer(KIND=ILkind) function count_hypercubic_indices(limits) result(n)
    integer, intent(in) :: limits(:)
    integer :: d

    n = 1_ILkind
    do d = 1, size(limits)
      n = n * int(limits(d), KIND=ILkind)
    end do
  end function count_hypercubic_indices


  recursive integer(KIND=ILkind) function count_simplex_indices(ndim, nmax) result(n)
    integer, intent(in) :: ndim, nmax
    integer :: kd

    if (ndim == 1) then
      n = int(nmax + 1, KIND=ILkind)
    else
      n = 0_ILkind
      do kd = 0, nmax
        n = n + count_simplex_indices(ndim - 1, nmax - kd)
      end do
    end if
  end function count_simplex_indices


  recursive integer(KIND=ILkind) function count_hyperbolic_indices(ndim, k_hyp) result(n)
    integer, intent(in) :: ndim, k_hyp
    integer :: kd, next_k

    if (ndim == 1) then
      n = int(k_hyp, KIND=ILkind)
    else
      n = 0_ILkind
      do kd = 0, k_hyp - 1
        next_k = k_hyp / (kd + 1)
        if (next_k >= 1) n = n + count_hyperbolic_indices(ndim - 1, next_k)
      end do
    end if
  end function count_hyperbolic_indices


  recursive subroutine rec_fill_hypercubic(dim, limits, current, nodes, pos)
    integer, intent(in) :: dim
    integer, intent(in) :: limits(:)
    integer, intent(inout) :: current(:)
    integer, intent(inout) :: nodes(:,:)
    integer(KIND=ILkind), intent(inout) :: pos
    integer :: kd

    if (dim == 1) then
      do kd = 0, limits(1) - 1
        current(1) = kd
        pos = pos + 1_ILkind
        nodes(:, int(pos)) = current
      end do
    else
      do kd = 0, limits(dim) - 1
        current(dim) = kd
        call rec_fill_hypercubic(dim - 1, limits, current, nodes, pos)
      end do
    end if
  end subroutine rec_fill_hypercubic


  recursive subroutine rec_fill_simplex(dim, remaining_sum, current, nodes, pos)
    integer, intent(in) :: dim, remaining_sum
    integer, intent(inout) :: current(:)
    integer, intent(inout) :: nodes(:,:)
    integer(KIND=ILkind), intent(inout) :: pos
    integer :: kd

    if (dim == 1) then
      do kd = 0, remaining_sum
        current(1) = kd
        pos = pos + 1_ILkind
        nodes(:, int(pos)) = current
      end do
    else
      do kd = 0, remaining_sum
        current(dim) = kd
        call rec_fill_simplex(dim - 1, remaining_sum - kd, current, nodes, pos)
      end do
    end if
  end subroutine rec_fill_simplex


  recursive subroutine rec_fill_hyperbolic(dim, remaining_k, current, nodes, pos)
    integer, intent(in) :: dim, remaining_k
    integer, intent(inout) :: current(:)
    integer, intent(inout) :: nodes(:,:)
    integer(KIND=ILkind), intent(inout) :: pos
    integer :: kd, next_k

    if (dim == 1) then
      do kd = 0, remaining_k - 1
        current(1) = kd
        pos = pos + 1_ILkind
        nodes(:, int(pos)) = current
      end do
    else
      do kd = 0, remaining_k - 1
        current(dim) = kd
        next_k = remaining_k / (kd + 1)
        if (next_k >= 1) call rec_fill_hyperbolic(dim - 1, next_k, current, nodes, pos)
      end do
    end if
  end subroutine rec_fill_hyperbolic

end module NDindex_hagedorn_m
