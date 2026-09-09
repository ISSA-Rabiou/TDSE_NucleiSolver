module HagedornOverlap_m
  USE QDUtil_m
  USE HagedornBasis_m
  USE NDindex_hagedorn_m
  USE polyortho_m
  implicit none
  private


  type, public :: HagedornPacket_t
    integer :: ndim = 0
    real(KIND=Rkind) :: hbar = ONE
    complex(KIND=Rkind), allocatable :: qvec(:)
    complex(KIND=Rkind), allocatable :: pvec(:)
    complex(KIND=Rkind), allocatable :: qmat(:,:)
    complex(KIND=Rkind), allocatable :: pmat(:,:)
    complex(KIND=Rkind)              :: actionS = CZERO
    complex(KIND=Rkind), allocatable :: inv_qmat(:,:)
    complex(KIND=Rkind), allocatable :: bigX(:,:)
    complex(KIND=Rkind), allocatable :: smallx(:)
  end type HagedornPacket_t

  type, public :: DetPhaseTracker_t
    real(KIND=Rkind) :: phase = ZERO
  end type DetPhaseTracker_t

  type :: OverlapAux_t
    complex(KIND=Rkind), allocatable :: bigU(:,:), bigV(:,:), bigW(:,:), bigW2(:,:)
    complex(KIND=Rkind), allocatable :: bigM1(:,:), bigN1(:,:), bigM2(:,:), bigN2(:,:)
    complex(KIND=Rkind), allocatable :: v1(:), v2(:), w1(:), w2(:), z1(:), z2(:)
  end type OverlapAux_t

  public :: init_hagedorn_packet
  public :: destroy_hagedorn_packet
  public :: overlap00
  public :: build_overlap_matrix
  public :: contract_overlap_matrix, build_legacy_overlap_matrix_1d
  public :: Cmat_to_QP, hermitian,maxabs_mat,Calc_S_Bourquin522
  public :: Calc_S_algebraic_genfunc,Calc_S_numeric_realgrid
  interface
    subroutine zgetrf(m, n, a, lda, ipiv, info)
      import :: Rkind
      integer, intent(in) :: m, n, lda
      complex(KIND=Rkind), intent(inout) :: a(lda, *)
      integer, intent(out) :: ipiv(*)
      integer, intent(out) :: info
    end subroutine zgetrf

    subroutine zgetri(n, a, lda, ipiv, work, lwork, info)
      import :: Rkind
      integer, intent(in) :: n, lda, lwork
      complex(KIND=Rkind), intent(inout) :: a(lda, *)
      integer, intent(in) :: ipiv(*)
      complex(KIND=Rkind), intent(inout) :: work(*)
      integer, intent(out) :: info
    end subroutine zgetri

    subroutine zgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
      import :: Rkind
      integer, intent(in) :: n, nrhs, lda, ldb
      complex(KIND=Rkind), intent(inout) :: a(lda, *)
      integer, intent(out) :: ipiv(*)
      complex(KIND=Rkind), intent(inout) :: b(ldb, *)
      integer, intent(out) :: info
    end subroutine zgesv

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

  subroutine init_hagedorn_packet(packet, qvec_in, pvec_in, qmat_in, pmat_in, actionS_in, hbar)
    type(HagedornPacket_t), intent(inout) :: packet
    complex(KIND=Rkind), intent(in) :: qvec_in(:), pvec_in(:), qmat_in(:,:), pmat_in(:,:)
    complex(KIND=Rkind), intent(in) :: actionS_in
    real(KIND=Rkind), intent(in), optional :: hbar
    integer :: ndim

    ndim = size(qvec_in)
    if (size(pvec_in) /= ndim) error stop 'init_hagedorn_packet: q and p have different sizes.'
    if (size(qmat_in,1) /= ndim .or. size(qmat_in,2) /= ndim) error stop 'init_hagedorn_packet: Q has wrong size.'
    if (size(pmat_in,1) /= ndim .or. size(pmat_in,2) /= ndim) error stop 'init_hagedorn_packet: P has wrong size.'

    call destroy_hagedorn_packet(packet)

    packet%ndim = ndim
    if (present(hbar)) packet%hbar = hbar

    allocate(packet%qvec(ndim), packet%pvec(ndim), packet%qmat(ndim,ndim), packet%pmat(ndim,ndim))
    packet%qvec = qvec_in
    packet%pvec = pvec_in
    packet%qmat = qmat_in
    packet%pmat = pmat_in
    packet%actionS = actionS_in

    call update_packet_cache(packet)
  end subroutine init_hagedorn_packet

  subroutine destroy_hagedorn_packet(packet)
    type(HagedornPacket_t), intent(inout) :: packet

    if (allocated(packet%qvec))      deallocate(packet%qvec)
    if (allocated(packet%pvec))      deallocate(packet%pvec)
    if (allocated(packet%qmat))      deallocate(packet%qmat)
    if (allocated(packet%pmat))      deallocate(packet%pmat)
    if (allocated(packet%inv_qmat))  deallocate(packet%inv_qmat)
    if (allocated(packet%bigX))      deallocate(packet%bigX)
    if (allocated(packet%smallx))    deallocate(packet%smallx)

    packet%ndim = 0
    packet%hbar = ONE
    packet%actionS = CZERO
  end subroutine destroy_hagedorn_packet

  subroutine update_packet_cache(packet)
    type(HagedornPacket_t), intent(inout) :: packet
    integer :: ndim

    ndim = packet%ndim
    packet%inv_qmat = matrix_inverse(packet%qmat)

    if (allocated(packet%bigX)) deallocate(packet%bigX)
    allocate(packet%bigX(2*ndim, ndim))
    packet%bigX(1:ndim, :)        = packet%qmat
    packet%bigX(ndim+1:2*ndim, :) = packet%pmat

    if (allocated(packet%smallx)) deallocate(packet%smallx)
    allocate(packet%smallx(2*ndim))
    packet%smallx(1:ndim)        = packet%qvec
    packet%smallx(ndim+1:2*ndim) = packet%pvec
  end subroutine update_packet_cache

  function overlap00(packet1, packet2, tracker) result(val)
    type(HagedornPacket_t), intent(in) :: packet1, packet2
    type(DetPhaseTracker_t), intent(inout), optional :: tracker
    complex(KIND=Rkind) :: val

    complex(KIND=Rkind), allocatable :: A1(:,:), A2(:,:), dA(:,:)
    complex(KIND=Rkind), allocatable :: xi1(:), xi2(:), dXi(:), rhs(:,:)
    complex(KIND=Rkind) :: eta1, eta2, dEta, exponent_term, pf
    complex(KIND=Rkind), allocatable :: matdet(:,:)
    real(KIND=Rkind) :: phase_local, det_abs
    type(DetPhaseTracker_t) :: local_tracker
    integer :: ndim

    ndim = packet1%ndim
    if (packet2%ndim /= ndim) error stop 'overlap00: incompatible dimensions.'

    A1 = matmul(packet1%pmat, packet1%inv_qmat)
    A2 = matmul(packet2%pmat, packet2%inv_qmat)
    dA = A2 - conjg(A1)

    xi1 = packet1%pvec - matmul(A1, packet1%qvec)
    eta1 = packet1%actionS - HALF * dot_noconj(xi1 + packet1%pvec, packet1%qvec)

    xi2 = packet2%pvec - matmul(A2, packet2%qvec)
    eta2 = packet2%actionS - HALF * dot_noconj(xi2 + packet2%pvec, packet2%qvec)

    dXi = xi2 - conjg(xi1)
    dEta = eta2 - conjg(eta1)

    matdet = (-EYE/TWO) * (matmul(hermitian(packet1%qmat), packet2%pmat) - &
                             matmul(hermitian(packet1%pmat), packet2%qmat))

    if (present(tracker)) then
      det_abs = detphase_abs(matdet, tracker)
      phase_local = tracker%phase
    else
      local_tracker%phase = ZERO
      det_abs = detphase_abs(matdet, local_tracker)
      phase_local = local_tracker%phase
    end if

    pf = sqrt(cmplx(det_abs, ZERO, KIND=Rkind))

    rhs = reshape(dXi, [ndim, 1])
    rhs = solve_linear_system(dA, rhs)
    exponent_term = -HALF * dot_noconj(dXi, rhs(:,1)) + dEta

    val = exp(-EYE * phase_local / TWO) / pf * exp(EYE / packet1%hbar * exponent_term)
  end function overlap00


  subroutine build_overlap_matrix(packet_bra, packet_ket, basis_bra, basis_ket, Smat, tracker)
    type(HagedornPacket_t), intent(in) :: packet_bra, packet_ket
    integer, intent(in) :: basis_bra(:,:), basis_ket(:,:)
    complex(KIND=Rkind), allocatable, intent(out) :: Smat(:,:)
    type(DetPhaseTracker_t), intent(inout), optional :: tracker

    type(OverlapAux_t) :: aux
    integer :: ndim, nbra, nket
    integer :: ib, ik, level, max_level
    integer, allocatable :: deg_bra(:), deg_ket(:)
    integer :: i0_bra, i0_ket

    ndim = packet_bra%ndim
    if (packet_ket%ndim /= ndim) error stop 'build_overlap_matrix: incompatible packet dimensions.'
    if (size(basis_bra,1) /= ndim) error stop 'build_overlap_matrix: wrong bra basis dimension.'
    if (size(basis_ket,1) /= ndim) error stop 'build_overlap_matrix: wrong ket basis dimension.'

    nbra = size(basis_bra, 2)
    nket = size(basis_ket, 2)

    allocate(Smat(nbra, nket))
    Smat = CZERO

    allocate(deg_bra(nbra), deg_ket(nket))
    do ib = 1, nbra
      deg_bra(ib) = sum(basis_bra(:, ib))
    end do
    do ik = 1, nket
      deg_ket(ik) = sum(basis_ket(:, ik))
    end do

    i0_bra = find_zero_column(basis_bra)
    i0_ket = find_zero_column(basis_ket)
    if (i0_bra <= 0 .or. i0_ket <= 0) error stop 'build_overlap_matrix: both bases must contain the zero multi-index.'

    call assign_overlap_aux(packet_bra, packet_ket, aux)
    Smat(i0_bra, i0_ket) = overlap00(packet_bra, packet_ket, tracker)

    max_level = maxval(deg_bra) + maxval(deg_ket)

    do level = 1, max_level
      do ib = 1, nbra
        do ik = 1, nket
          if (deg_bra(ib) + deg_ket(ik) /= level) cycle
          Smat(ib, ik) = overlap_jk_from_matrix(basis_bra(:,ib), basis_ket(:,ik), &
                                                basis_bra, basis_ket, Smat, aux)
        end do
      end do
    end do
  end subroutine build_overlap_matrix

  function contract_overlap_matrix(coeff_bra, coeff_ket, Smat) result(val)
    complex(KIND=Rkind), intent(in) :: coeff_bra(:), coeff_ket(:)
    complex(KIND=Rkind), intent(in) :: Smat(:,:)
    complex(KIND=Rkind) :: val
    integer :: ib, ik

    if (size(coeff_bra) /= size(Smat,1)) error stop 'contract_overlap_matrix: coeff_bra has wrong size.'
    if (size(coeff_ket) /= size(Smat,2)) error stop 'contract_overlap_matrix: coeff_ket has wrong size.'

    val = CZERO
    do ib = 1, size(Smat,1)
      do ik = 1, size(Smat,2)
        val = val + conjg(coeff_bra(ib)) * coeff_ket(ik) * Smat(ib,ik)
      end do
    end do
  end function contract_overlap_matrix

  subroutine assign_overlap_aux(packet1, packet2, aux)
    type(HagedornPacket_t), intent(in) :: packet1, packet2
    type(OverlapAux_t), intent(out) :: aux

    integer :: ndim
    complex(KIND=Rkind), allocatable :: J(:,:), delta_x(:), UdagU(:,:), UbarUT(:,:)
    real(KIND=Rkind) :: scalar_w

    ndim = packet1%ndim
    J = build_symplectic_matrix(ndim)

    aux%bigU = (EYE/TWO) * symp_form_mm(packet2%bigX, conjg(packet1%bigX), J)
    aux%bigV = (EYE/TWO) * symp_form_mm(packet2%bigX, packet1%bigX, J)

    delta_x = packet1%smallx - packet2%smallx
    aux%v1 = (EYE / sqrt(TWO * packet1%hbar)) * symp_form_mv(packet2%bigX, delta_x, J)
    aux%v2 = (EYE / sqrt(TWO * packet1%hbar)) * symp_form_mv(packet1%bigX, -delta_x, J)

    aux%w1 = -matmul(transpose(aux%bigV), conjg(aux%v1)) + aux%v2
    aux%w2 =  matmul(conjg(aux%bigV), aux%v2) + conjg(aux%v1)

    if (ndim == 1) then
      allocate(aux%bigW(1,1), aux%bigW2(1,1), aux%bigM1(1,1), aux%bigN1(1,1), aux%bigM2(1,1), aux%bigN2(1,1))
      scalar_w = ONE / abs(aux%bigU(1,1))**2
      aux%bigW(1,1)  = cmplx(scalar_w, ZERO, KIND=Rkind)
      aux%bigW2(1,1) = aux%bigW(1,1)
      aux%bigM1(1,1) = aux%bigW(1,1) * conjg(aux%bigU(1,1))
      aux%bigN1(1,1) = aux%bigW(1,1) * aux%bigV(1,1) * conjg(aux%bigU(1,1))
      ! Direct translation of the validated Python scalar branch.
      aux%bigM2(1,1) = aux%bigW(1,1) * aux%bigU(1,1)
      aux%bigN2(1,1) = aux%bigW(1,1) * conjg(aux%bigV(1,1)) * conjg(aux%bigU(1,1))
      aux%z1 = aux%bigW(1,1) * aux%w1
      aux%z2 = aux%bigW(1,1) * aux%w2
    else
      UdagU = matmul(hermitian(aux%bigU), aux%bigU)
      UbarUT = matmul(conjg(aux%bigU), transpose(aux%bigU))

      aux%bigW  = matrix_inverse(UdagU)
      aux%bigW2 = matrix_inverse(UbarUT)

      aux%bigM1 = matmul(aux%bigW, hermitian(aux%bigU))
      aux%bigN1 = matmul(aux%bigW, matmul(transpose(aux%bigV), conjg(aux%bigU)))
      aux%bigM2 = matmul(aux%bigW2, conjg(aux%bigU))
      aux%bigN2 = matmul(aux%bigW2, matmul(conjg(aux%bigV), hermitian(aux%bigU)))

      aux%z1 = matmul(aux%bigW,  aux%w1)
      aux%z2 = matmul(aux%bigW2, aux%w2)
    end if
  end subroutine assign_overlap_aux

  function overlap_jk_from_matrix(j, k, basis_bra, basis_ket, Smat, aux) result(val)
    integer, intent(in) :: j(:), k(:)
    integer, intent(in) :: basis_bra(:,:), basis_ket(:,:)
    complex(KIND=Rkind), intent(in) :: Smat(:,:)
    type(OverlapAux_t), intent(in) :: aux
    complex(KIND=Rkind) :: val

    integer :: ndim, d, i, idx_jred, idx_kred, idx_dep, idx_j, idx_k
    integer, allocatable :: jred(:), kred(:), tmp(:)
    complex(KIND=Rkind), allocatable :: overlap_k_m(:), overlap_j_m(:)
    complex(KIND=Rkind), allocatable :: row_k(:), row_j(:), zvec(:)
    integer :: sumj, sumk
    real(KIND=Rkind) :: denom

    ndim = size(j)
    sumj = sum(j)
    sumk = sum(k)

    idx_j = find_column_vector(basis_bra, j)
    idx_k = find_column_vector(basis_ket, k)

    if (sumj == 0 .and. sumk == 0) then
      val = Smat(idx_j, idx_k)
      return
    end if

    allocate(jred(ndim), kred(ndim), tmp(ndim), overlap_k_m(ndim), overlap_j_m(ndim))

    if (sumj > 0 .and. (sumk == 0 .or. sumj <= sumk)) then
      d = smallest_positive_component_temp(j)
      jred = j
      jred(d) = jred(d) - 1
      idx_jred = find_column_vector(basis_bra, jred)
      if (idx_jred <= 0) error stop 'overlap_jk_from_matrix: missing reduced bra index.'

      do i = 1, ndim
        if (k(i) > 0) then
          tmp = k
          tmp(i) = tmp(i) - 1
          idx_dep = find_column_vector(basis_ket, tmp)
          if (idx_dep > 0) then
            ! Python recurrence uses the already reduced bra index j_ here.
            overlap_k_m(i) = sqrt(real(k(i), Rkind)) * Smat(idx_jred, idx_dep)
          else
            overlap_k_m(i) = CZERO
          end if
        else
          overlap_k_m(i) = CZERO
        end if
      end do

      do i = 1, ndim
        if (jred(i) > 0) then
          tmp = jred
          tmp(i) = tmp(i) - 1
          idx_dep = find_column_vector(basis_bra, tmp)
          if (idx_dep > 0) then
            overlap_j_m(i) = sqrt(real(jred(i), Rkind)) * Smat(idx_dep, idx_k)
          else
            overlap_j_m(i) = CZERO
          end if
        else
          overlap_j_m(i) = CZERO
        end if
      end do

      row_k = aux%bigM1(d, :)
      row_j = -aux%bigN1(d, :)
      zvec  = aux%z1
      denom = sqrt(real(j(d), Rkind))

      val = dot_noconj(row_k, overlap_k_m) + dot_noconj(row_j, overlap_j_m) + &
            zvec(d) * Smat(idx_jred, idx_k)
      val = val / denom

    else
      d = smallest_positive_component_temp(k)
      kred = k
      kred(d) = kred(d) - 1
      idx_kred = find_column_vector(basis_ket, kred)
      if (idx_kred <= 0) error stop 'overlap_jk_from_matrix: missing reduced ket index.'

      do i = 1, ndim
        if (kred(i) > 0) then
          tmp = kred
          tmp(i) = tmp(i) - 1
          idx_dep = find_column_vector(basis_ket, tmp)
          if (idx_dep > 0) then
            overlap_k_m(i) = sqrt(real(kred(i), Rkind)) * Smat(idx_j, idx_dep)
          else
            overlap_k_m(i) = CZERO
          end if
        else
          overlap_k_m(i) = CZERO
        end if
      end do

      do i = 1, ndim
        if (j(i) > 0) then
          tmp = j
          tmp(i) = tmp(i) - 1
          idx_dep = find_column_vector(basis_bra, tmp)
          if (idx_dep > 0) then
            ! Python recurrence uses the already reduced ket index k_ here.
            overlap_j_m(i) = sqrt(real(j(i), Rkind)) * Smat(idx_dep, idx_kred)
          else
            overlap_j_m(i) = CZERO
          end if
        else
          overlap_j_m(i) = CZERO
        end if
      end do

      row_k = aux%bigN2(d, :)
      row_j = aux%bigM2(d, :)
      zvec  = aux%z2
      denom = sqrt(real(k(d), Rkind))

      val = dot_noconj(row_k, overlap_k_m) + dot_noconj(row_j, overlap_j_m) + &
            zvec(d) * Smat(idx_j, idx_kred)
      val = val / denom
    end if
  end function overlap_jk_from_matrix

  function build_symplectic_matrix(ndim) result(J)
    integer, intent(in) :: ndim
    complex(KIND=Rkind), allocatable :: J(:,:)

    allocate(J(2*ndim, 2*ndim))
    J = CZERO
    J(1:ndim, ndim+1:2*ndim) = identity_matrix(ndim)
    J(ndim+1:2*ndim, 1:ndim) = -identity_matrix(ndim)
  end function build_symplectic_matrix

  function identity_matrix(n) result(I)
    integer, intent(in) :: n
    complex(KIND=Rkind) :: I(n,n)
    integer :: d

    I = CZERO
    do d = 1, n
      I(d,d) = CONE
    end do
  end function identity_matrix

  function hermitian(A) result(Ah)
    complex(KIND=Rkind), intent(in) :: A(:,:)
    complex(KIND=Rkind), allocatable :: Ah(:,:)
    Ah = transpose(conjg(A))
  end function hermitian

  function symp_form_mm(a, b, J) result(c)
    complex(KIND=Rkind), intent(in) :: a(:,:), b(:,:), J(:,:)
    complex(KIND=Rkind), allocatable :: c(:,:)
    c = matmul(transpose(a), matmul(J, b))
  end function symp_form_mm

  function symp_form_mv(a, b, J) result(c)
    complex(KIND=Rkind), intent(in) :: a(:,:), b(:), J(:,:)
    complex(KIND=Rkind), allocatable :: c(:)
    c = matmul(transpose(a), matmul(J, b))
  end function symp_form_mv

  function dot_noconj(x, y) result(val)
    complex(KIND=Rkind), intent(in) :: x(:), y(:)
    complex(KIND=Rkind) :: val
    integer :: i

    val = CZERO
    do i = 1, size(x)
      val = val + x(i) * y(i)
    end do
  end function dot_noconj

  function matrix_inverse(A) result(Ainv)
    complex(KIND=Rkind), intent(in) :: A(:,:)
    complex(KIND=Rkind), allocatable :: Ainv(:,:)
    complex(KIND=Rkind), allocatable :: work(:)
    integer, allocatable :: ipiv(:)
    integer :: n, info, lwork
    complex(KIND=Rkind) :: work_query(1)

    n = size(A,1)
    if (size(A,2) /= n) error stop 'matrix_inverse: matrix must be square.'

    allocate(Ainv(n,n), ipiv(n))
    Ainv = A

    call zgetrf(n, n, Ainv, n, ipiv, info)
    if (info /= 0) error stop 'matrix_inverse: LU factorization failed.'

    lwork = -1
    call zgetri(n, Ainv, n, ipiv, work_query, lwork, info)
    lwork = max(1, int(real(work_query(1), KIND=Rkind)))
    allocate(work(lwork))
    call zgetri(n, Ainv, n, ipiv, work, lwork, info)
    if (info /= 0) error stop 'matrix_inverse: matrix inversion failed.'
  end function matrix_inverse

  function solve_linear_system(A, B) result(X)
    complex(KIND=Rkind), intent(in) :: A(:,:), B(:,:)
    complex(KIND=Rkind), allocatable :: X(:,:)
    complex(KIND=Rkind), allocatable :: Acopy(:,:)
    integer, allocatable :: ipiv(:)
    integer :: n, nrhs, info

    n = size(A,1)
    if (size(A,2) /= n) error stop 'solve_linear_system: A must be square.'
    if (size(B,1) /= n) error stop 'solve_linear_system: incompatible right-hand side.'

    nrhs = size(B,2)
    allocate(Acopy(n,n), X(n,nrhs), ipiv(n))
    Acopy = A
    X = B

    call zgesv(n, nrhs, Acopy, n, ipiv, X, n, info)
    if (info /= 0) error stop 'solve_linear_system: linear solve failed.'
  end function solve_linear_system

  function complex_determinant(A) result(detA)
    complex(KIND=Rkind), intent(in) :: A(:,:)
    complex(KIND=Rkind) :: detA
    complex(KIND=Rkind), allocatable :: Acopy(:,:)
    integer, allocatable :: ipiv(:)
    integer :: n, info, i

    n = size(A,1)
    if (size(A,2) /= n) error stop 'complex_determinant: matrix must be square.'

    allocate(Acopy(n,n), ipiv(n))
    Acopy = A
    call zgetrf(n, n, Acopy, n, ipiv, info)
    if (info /= 0) error stop 'complex_determinant: LU factorization failed.'

    detA = CONE
    do i = 1, n
      detA = detA * Acopy(i,i)
      if (ipiv(i) /= i) detA = -detA
    end do
  end function complex_determinant

  real(KIND=Rkind) function detphase_abs(A, tracker) result(absdet)
    complex(KIND=Rkind), intent(in) :: A(:,:)
    type(DetPhaseTracker_t), intent(inout) :: tracker
    complex(KIND=Rkind) :: detA
    real(KIND=Rkind) :: phi, phase_diff, turns

    detA = complex_determinant(A)
    phi = atan2(aimag(detA), real(detA, KIND=Rkind))
    turns = floor(tracker%phase / (TWO * PI))

    if (phi < ZERO) phi = phi + TWO * PI
    phi = phi + turns * TWO * PI

    phase_diff = phi - tracker%phase
    tracker%phase = phi

    if (abs(phase_diff) > PI) then
      if (phase_diff > ZERO) then
        tracker%phase = tracker%phase - TWO * PI
      else
        tracker%phase = tracker%phase + TWO * PI
      end if
    end if

    absdet = abs(detA)
  end function detphase_abs

  integer function find_zero_column(basis) result(idx)
    integer, intent(in) :: basis(:,:)
    integer :: j

    idx = 0
    do j = 1, size(basis, 2)
      if (all(basis(:,j) == 0)) then
        idx = j
        return
      end if
    end do
  end function find_zero_column

  integer function find_column_vector(basis, target) result(idx)
    integer, intent(in) :: basis(:,:)
    integer, intent(in) :: target(:)
    integer :: j

    idx = 0
    do j = 1, size(basis, 2)
      if (all(basis(:,j) == target)) then
        idx = j
        return
      end if
    end do
  end function find_column_vector

  integer function smallest_positive_component_temp(k) result(pos)
    integer, intent(in) :: k(:)
    integer :: i, best

    pos = 0
    best = huge(1)
    do i = 1, size(k)
      if (k(i) > 0 .and. k(i) < best) then
        best = k(i)
        pos = i
      end if
    end do
    if (pos == 0) error stop 'smallest_positive_component_temp: no positive component found.'
  end function smallest_positive_component_temp

  subroutine Cmat_to_QP(Cmat, Qmat, Pmat)
    complex(KIND=Rkind), intent(in) :: Cmat(:,:)
    complex(KIND=Rkind), allocatable, intent(out) :: Qmat(:,:), Pmat(:,:)

    integer :: ndim
    complex(KIND=Rkind), allocatable :: ImC(:,:), sqrtImC(:,:)

    ndim = size(Cmat,1)
    if (size(Cmat,2) /= ndim) error stop 'Cmat_to_QP: matrix must be square.'

    allocate(ImC(ndim,ndim))
    ImC = (Cmat - conjg(transpose(Cmat))) / (TWO*EYE)

    sqrtImC = matrix_sqrt_hermitian(ImC)
    Qmat = matrix_inverse(sqrtImC)

    allocate(Pmat(ndim,ndim))
    Pmat = matmul(Cmat,Qmat)
  end subroutine Cmat_to_QP

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

    pure subroutine overlap_core_coeffs(q0, alpha0, p0, qt, alphat, pt, Acoef, Bcoef, Ccoef, s0, st)

    real(kind=Rkind), intent(in) :: q0, p0, qt, pt

    complex(kind=Rkind), intent(in) :: alpha0, alphat

    complex(kind=Rkind), intent(out) :: Acoef, Bcoef, Ccoef

    real(kind=Rkind), intent(out) :: s0, st

    real(kind=Rkind) :: a0, b0, at, bt

    a0 = real(alpha0, kind=Rkind)

    b0 = aimag(alpha0)

    at = real(alphat, kind=Rkind)

    bt = aimag(alphat)

    s0 = sqrt(a0)

    st = sqrt(at)

    Acoef = HALF * ( cmplx(a0, b0, kind=Rkind) + cmplx(at, -bt, kind=Rkind) )

    Bcoef = cmplx(a0, b0, kind=Rkind) * q0 + cmplx(at, -bt, kind=Rkind) * qt + EYE * (p0 - pt)

    Ccoef = -HALF*cmplx(a0, b0, kind=Rkind) * q0*q0 - HALF*cmplx(at, -bt, kind=Rkind) * qt*qt &

            - EYE*p0*q0 + EYE*pt*qt

  end subroutine overlap_core_coeffs


    pure real(kind=Rkind) function fact_int(n) result(f)

    integer, intent(in) :: n

    integer :: k

    f = ONE

    if (n <= 1) return

    do k = 2, n

      f = f * real(k, kind=Rkind)

    end do

  end function fact_int



  pure real(kind=Rkind) function maxabs_mat(A) result(v)

    complex(kind=Rkind), intent(in) :: A(:,:)

    v = maxval(abs(A))

  end function maxabs_mat



    subroutine Calc_S_numeric_realgrid(S, nb, nq, Q0, SQ0, A0, P0, Qt, SQt, At, Pt)

    integer, intent(in) :: nb, nq

    real(kind=Rkind), intent(in) :: Q0, SQ0, P0, Qt, SQt, Pt

    complex(kind=Rkind), intent(in) :: A0, At

    complex(kind=Rkind), intent(out) :: S(nb,nb)

    integer :: ib, iq

    real(kind=Rkind) :: SQeq, Qeq

    real(kind=Rkind), allocatable :: Q(:), W(:)

    complex(kind=Rkind), allocatable :: d0gb0(:,:), d0gb1(:,:), d0bgw(:,:)



    SQeq = sqrt(SQ0*SQ0 + SQt*SQt) / sqrt(TWO)

    Qeq  = (SQ0*SQ0*Q0 + SQt*SQt*Qt) / (SQ0*SQ0 + SQt*SQt)



    allocate(Q(nq), W(nq), d0gb0(nq,nb), d0gb1(nq,nb), d0bgw(nb,nq))

    call hermgauss_plain(nq, Q, W)

    W = W / SQeq

    Q = Qeq + Q / SQeq



    do iq = 1, nq

      do ib = 1, nb

        call basis_value_legacy(Q(iq), Q0, A0, P0, ib-1, d0gb0(iq,ib))

        call basis_value_legacy(Q(iq), Qt, At, Pt, ib-1, d0gb1(iq,ib))

      end do

    end do



    d0bgw = transpose(d0gb1)

    do ib = 1, nb

      d0bgw(ib,:) = d0bgw(ib,:) * W(:)

    end do



    S = matmul(conjg(d0bgw), d0gb0)

    deallocate(Q, W, d0gb0, d0gb1, d0bgw)

  end subroutine Calc_S_numeric_realgrid




    subroutine Calc_S_algebraic_genfunc(S, nb, q0, sq0_in, alpha0, p0, qt, sqt_in, alphat, pt)

    integer, intent(in) :: nb

    real(kind=Rkind), intent(in) :: q0, sq0_in, p0, qt, sqt_in, pt

    complex(kind=Rkind), intent(in) :: alpha0, alphat

    complex(kind=Rkind), intent(out) :: S(nb,nb)



    integer :: n, j, k, l, m

    real(kind=Rkind) :: s0, st

    complex(kind=Rkind) :: Acoef, Bcoef, Ccoef, Kpref

    complex(kind=Rkind) :: au, av, cuv, du, dv, tmp

    complex(kind=Rkind), allocatable :: cu(:), cv(:)



    call overlap_core_coeffs(q0, alpha0, p0, qt, alphat, pt, Acoef, Bcoef, Ccoef, s0, st)

    if (abs(s0 - sq0_in) > 1.0e-12_Rkind .or. abs(st - sqt_in) > 1.0e-12_Rkind) then

      write(out_unit,'(A)') 'WARNING: SQ0/SQt are not consistent with sqrt(real(alpha)); exact methods follow alpha.'

    end if



    Kpref = sqrt( cmplx(s0*st, ZERO, kind=Rkind) / Acoef ) * exp(Ccoef + Bcoef*Bcoef/(FOUR*Acoef))

    au    = cmplx(st*st, ZERO, kind=Rkind) / (TWO*Acoef) - HALF

    av    = cmplx(s0*s0, ZERO, kind=Rkind) / (TWO*Acoef) - HALF

    cuv   = cmplx(s0*st, ZERO, kind=Rkind) / Acoef

    du    = cmplx(st/sqrt(TWO), ZERO, kind=Rkind) * (Bcoef/Acoef - TWO*qt)

    dv    = cmplx(s0/sqrt(TWO), ZERO, kind=Rkind) * (Bcoef/Acoef - TWO*q0)



    allocate(cu(0:nb-1), cv(0:nb-1))

    cu = CZERO ; cv = CZERO



    do n = 0, nb-1

      do m = 0, n/2

        cu(n) = cu(n) + au**m / fact_int(m) * du**(n-2*m) / fact_int(n-2*m)

        cv(n) = cv(n) + av**m / fact_int(m) * dv**(n-2*m) / fact_int(n-2*m)

      end do

    end do



    S = CZERO

    do j = 0, nb-1

      do k = 0, nb-1

        tmp = CZERO

        do l = 0, min(j,k)

          tmp = tmp + cuv**l / fact_int(l) * cu(j-l) * cv(k-l)

        end do

        S(j+1,k+1) = Kpref * tmp * sqrt( fact_int(j) * fact_int(k) )

      end do

    end do



    deallocate(cu, cv)

  end subroutine Calc_S_algebraic_genfunc


    subroutine Calc_S_Bourquin522(S, nb, q0, sq0_in, alpha0, p0, qt, sqt_in, alphat, pt, nq_exact)

    integer, intent(in) :: nb

    real(kind=Rkind), intent(in) :: q0, sq0_in, p0, qt, sqt_in, pt

    complex(kind=Rkind), intent(in) :: alpha0, alphat

    integer, intent(in), optional :: nq_exact

    complex(kind=Rkind), intent(out) :: S(nb,nb)



    integer :: nq, i, j, k

    real(kind=Rkind), allocatable :: z(:), w(:)

    real(kind=Rkind) :: s0, st

    complex(kind=Rkind) :: Acoef, Bcoef, Ccoef, mu, sqrtA, pref, x

    complex(kind=Rkind), allocatable :: pbra(:), pket(:)



    call overlap_core_coeffs(q0, alpha0, p0, qt, alphat, pt, Acoef, Bcoef, Ccoef, s0, st)

    if (abs(s0 - sq0_in) > 1.0e-12_Rkind .or. abs(st - sqt_in) > 1.0e-12_Rkind) then

      write(out_unit,'(A)') 'WARNING: SQ0/SQt are not consistent with sqrt(real(alpha)); exact methods follow alpha.'

    end if



    nq = nb

    if (present(nq_exact)) nq = max(nq_exact, nb)



    allocate(z(nq), w(nq), pbra(nb), pket(nb))

    call hermgauss_standard(nq, z, w)



    mu    = Bcoef / (TWO*Acoef)

    sqrtA = sqrt(Acoef)

    pref  = sqrt( cmplx(s0*st, ZERO, kind=Rkind) / Acoef ) * exp(Ccoef + Bcoef*Bcoef/(FOUR*Acoef))



    S = CZERO

    do i = 1, nq

      x = mu + cmplx(z(i), ZERO, kind=Rkind) / sqrtA

      do j = 1, nb

        pbra(j) = poly_hermite_c( cmplx(st,ZERO,kind=Rkind) * (x - qt), j-1 )

        pket(j) = poly_hermite_c( cmplx(s0,ZERO,kind=Rkind) * (x - q0), j-1 )

      end do

      do j = 1, nb

        do k = 1, nb

          S(j,k) = S(j,k) + pref * w(i) * pbra(j) * pket(k)

        end do

      end do

    end do



    deallocate(z, w, pbra, pket)

  end subroutine Calc_S_Bourquin522





  subroutine build_legacy_overlap_matrix_1d(Smat, nb, q0, sq0, alpha0, p0, qt, sqt, alphat, pt)

    !!-------------------------------------------------------------------

    !! Exact 1D overlap matrix between two legacy HO/Hagedorn direct factors.

    !!

    !! Legacy convention in Basis_m:

    !!   phi_n(x) = a^(1/4) H_n(sqrt(a)(x-q)) exp[-a (x-q)^2/2]

    !!              * exp[-i b (x-q)^2/2 + i p (x-q)]

    !! with alpha = a - i b, scaleQ = sqrt(a).

    !!

    !! Mapping to the homogeneous Hagedorn packet is

    !!   q = q0,

    !!   p = p0,

    !!   Q = 1 / scaleQ,

    !!   P = (i * conjg(alpha)) * Q,

    !!   S = 0.

    !!-------------------------------------------------------------------

    complex(KIND=Rkind), allocatable, intent(out) :: Smat(:,:)

    integer, intent(in) :: nb

    real(KIND=Rkind), intent(in) :: q0, sq0, p0, qt, sqt, pt

    complex(KIND=Rkind), intent(in) :: alpha0, alphat

    type(HagedornPacket_t) :: bra_packet, ket_packet

    type(HagedornBasisShape_t) :: shape

    complex(KIND=Rkind) :: qvec0(1), pvec0(1), qmat0(1,1), pmat0(1,1)

    complex(KIND=Rkind) :: qvect(1), pvect(1), qmatt(1,1), pmatt(1,1)

    if (nb <= 0) error stop 'build_legacy_overlap_matrix_1d: nb must be positive.'

    if (sq0 <= ZERO .or. sqt <= ZERO) error stop 'build_legacy_overlap_matrix_1d: non-positive scaleQ.'

    qvec0(1) = cmplx(q0, ZERO, KIND=Rkind)

    pvec0(1) = cmplx(p0, ZERO, KIND=Rkind)

    qmat0(1,1) = cmplx(ONE / sq0, ZERO, KIND=Rkind)

    pmat0(1,1) = EYE * conjg(alpha0) * qmat0(1,1)

    qvect(1) = cmplx(qt, ZERO, KIND=Rkind)

    pvect(1) = cmplx(pt, ZERO, KIND=Rkind)

    qmatt(1,1) = cmplx(ONE / sqt, ZERO, KIND=Rkind)

    pmatt(1,1) = EYE * conjg(alphat) * qmatt(1,1)

    call init_hypercubic_shape(shape, [nb])

    call init_hagedorn_packet(bra_packet, qvec0, pvec0, qmat0, pmat0, CZERO, hbar=ONE)

    call init_hagedorn_packet(ket_packet, qvect, pvect, qmatt, pmatt, CZERO, hbar=ONE)

    call build_overlap_matrix(bra_packet, ket_packet, shape%nodes, shape%nodes, Smat)

    call destroy_hagedorn_packet(bra_packet)

    call destroy_hagedorn_packet(ket_packet)

    call destroy_hagedorn_shape(shape)

  end subroutine build_legacy_overlap_matrix_1d


end module HagedornOverlap_m
