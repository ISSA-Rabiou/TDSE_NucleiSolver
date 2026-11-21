module Hagedorn_m
   USE QDUtil_m
   Use Basis_m
   Use Psi_m
   Use Op_m
   USE Ana_psi_m
   USE polyortho_m
   USE sub_propa_m
   USE HagedornOverlap_m
   USE HagedornBasis_m, ONLY: Cmat_to_QP
   USE HagedornLinalg_m, ONLY: max_abs_c
   implicit none
   PRIVATE
     !> Public interfaces for Hagedorn basis construction and projection operations
   PUBLIC :: Calc_Basis_parameters ! Compute basis parameters from a psi structure
   PUBLIC :: march_Global                     ! Global propagation step
   
     !> Hagedorn propagation and testing routines
   PUBLIC :: Hagedorn                         ! Hagedorn routine
   PUBLIC :: Test_Calc_S       ! Test routine for computing matrix S
   PUBLIC :: Test_projection_hagedorn         ! Test of projection method using Hagedorn basis
   
     !> Utility routines for basis construction and manipulation
   PUBLIC :: build_identity_matrix            ! Build identity matrix for trivial bases
   PUBLIC :: build_overlap_matrices         ! Compute overlap matrices
   
     !> Hagedorn basis construction routines
   PUBLIC :: Construct_Hagedorn_none_Variational_Basis       ! Construct non-variational Hagedorn basis
   
   PUBLIC :: Test_Compare_Overlap_ThreeMethods        ! Test routine for computing matrix S with three methods
   
CONTAINS


subroutine build_identity_matrix(Id, n)
  implicit none
  integer, intent(in)                            :: n
  complex(kind=Rkind), allocatable, intent(out)  :: Id(:,:)
  integer                                        :: ib

  ! Allocate the identity matrix
  allocate(Id(n,n))

  ! Initialize all elements to zero
  Id(:,:) = CZERO

  ! Set diagonal elements to 1
  do ib = 1, n
     Id(ib,ib) = CONE
  end do

end subroutine build_identity_matrix

SUBROUTINE Construct_Hagedorn_none_Variational_Basis(Basis, Qt, SQt, At, Pt)
  USE QDUtil_m
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE

  ! === Interface ===
  TYPE(Basis_t),            INTENT(INOUT)  :: Basis         ! Basis to be modified
  REAL(kind=Rkind),         INTENT(IN)     :: Qt(:)         ! <Q>
  REAL(kind=Rkind),         INTENT(IN)     :: SQt(:)        ! sqrt(<Q²> - <Q>²)
  COMPLEX(kind=Rkind),      INTENT(IN)     :: At(:)         ! Complex width A
  REAL(kind=Rkind),         INTENT(IN)     :: Pt(:)         ! <P>

  ! === Local variables ===
  REAL(kind=Rkind), ALLOCATABLE :: Q0(:), SQ0(:), P0(:)
  COMPLEX(kind=Rkind), ALLOCATABLE :: A0(:)
  INTEGER :: ndim, ib, nb, nq

  ! === Determine dimensionality and validate parameters before rebuilding
  IF (.NOT. ALLOCATED(Basis%tab_basis)) ERROR STOP 'Hagedorn basis has no primitive bases.'
  ndim = SIZE(Basis%tab_basis) - 1
  IF (ndim < 1) ERROR STOP 'Hagedorn basis must contain at least one nuclear dimension.'
  IF (SIZE(Qt) /= ndim .OR. SIZE(SQt) /= ndim .OR. SIZE(At) /= ndim .OR. SIZE(Pt) /= ndim) &
     ERROR STOP 'Hagedorn basis parameter arrays do not match the basis dimension.'
  DO ib = 1, ndim
    IF (.NOT. IEEE_IS_FINITE(Qt(ib)) .OR. .NOT. IEEE_IS_FINITE(SQt(ib)) .OR. &
        .NOT. IEEE_IS_FINITE(Pt(ib)) .OR. .NOT. IEEE_IS_FINITE(REAL(At(ib),KIND=Rkind)) .OR. &
        .NOT. IEEE_IS_FINITE(AIMAG(At(ib)))) ERROR STOP 'Non-finite Hagedorn basis parameter.'
    IF (SQt(ib) <= SQRT(TINY(ONE)) .OR. REAL(At(ib),KIND=Rkind) <= ZERO) &
       ERROR STOP 'Hagedorn basis width must be strictly positive.'
  END DO

  ! === Allocate arrays to store current basis parameters
  ALLOCATE(Q0(ndim), SQ0(ndim), P0(ndim), A0(ndim))
  !<--- Initialize arrays to zero
  Q0(:) = ZERO
  SQ0(:) = ZERO
  P0(:) = ZERO
  A0(:) = CZERO

  ! === Get original basis parameters
  CALL Get_Basis_Parameters(Basis, Q0, SQ0, A0, P0)

  ! === Change the parameters of the basis to new values
  CALL Change_Basis_Parameters(Basis, Qt, SQt, At, Pt)

  ! === Rebuild the primitive basis
  CALL construct_primitive_basis_temp(Basis)

  ! === Loop over all dimensions to compute S matrix
  DO ib = 1, ndim
    nb = Basis%tab_basis(ib)%nb
    nq = Basis%tab_basis(ib)%nq

    IF (Basis%tab_basis(ib)%Basis_name == 'herm' .OR. &
        Basis%tab_basis(ib)%Basis_name == 'ho') THEN
      !Using algebric routine
      IF (.NOT. ALLOCATED(Basis%tab_basis(ib)%S)) ALLOCATE(Basis%tab_basis(ib)%S(nb,nb))
      IF (SIZE(Basis%tab_basis(ib)%S,1) /= nb .OR. SIZE(Basis%tab_basis(ib)%S,2) /= nb) THEN
        DEALLOCATE(Basis%tab_basis(ib)%S)
        ALLOCATE(Basis%tab_basis(ib)%S(nb,nb))
      END IF
             call Calc_S_algebraic_genfunc(Basis%tab_basis(ib)%S, nb, Q0(ib), SQ0(ib), A0(ib), P0(ib), &
                                           Qt(ib), SQt(ib), At(ib), Pt(ib))
  
      !CALL Calc_S(Basis%tab_basis(ib)%S, nb, nq, &
      !            Q0(ib), SQ0(ib), A0(ib), P0(ib), &
      !            Qt(ib), SQt(ib), At(ib), Pt(ib))
    ELSE
      CALL build_identity_matrix(Basis%tab_basis(ib)%S, nb)
    END IF
  END DO

  ! === Cleanup: deallocate local arrays
  DEALLOCATE(Q0, SQ0, P0, A0)

END SUBROUTINE Construct_Hagedorn_none_Variational_Basis

!=======================================================================
!  SUBROUTINE: Calc_Basis_parameters
!
!  PURPOSE:
!    This routine computes key observables of a quantum wavefunction ψ
!    expressed in a given basis. The following quantities are evaluated:
!      - Qt(:)   : Expectation value(s) of position <Q>
!      - SQt(:)  : Standard deviation(s) of position
!      - Pt(:)   : Expectation value(s) of momentum <P>
!      - At(:)   : Complex width parameter "A" (related to squeezing/β)
!
!    The behavior depends on whether propagation settings (propa) are given:
!
!      1) Without "propa" (diagnostic mode):
!         → Compute all quantities (Qt, SQt, Pt, At).
!
!      2) With "propa" (propagation mode):
!         → Initialize Pt=0, At=SQt^2 by default.
!         → If propa%P    = .TRUE., compute <P>.
!         → If propa%Beta = .TRUE., compute "A".
!
!    This design allows a unified interface:
!       - In analysis mode (no "propa"): always compute everything.
!       - In propagation mode: compute only what is needed, saving CPU time.
!
!  INPUT:
!    psi   : Wavefunction object of type psi_t (contains basis + coefficients).
!    propa : (optional) Propagation settings of type propa_t.
!
!  OUTPUT (modified in place):
!    Qt(:), SQt(:), Pt(:), At(:)
!
!  NOTES:
!    - Uses an integration index table Tab_Iq, allocated internally.
!    - Ensures clean memory handling (Tab_Iq is always deallocated).
!=======================================================================
SUBROUTINE Calc_Basis_parameters(psi, Qt, SQt, At, Pt, propa, need_Pt, need_At)
   USE psi_m
   IMPLICIT NONE

   ! === Interfaces ===
   TYPE(psi_t),          INTENT(IN)     :: psi        ! Wavefunction input
   REAL(kind=Rkind),     INTENT(INOUT)  :: Qt(:)      ! Position expectation value <Q>
   REAL(kind=Rkind),     INTENT(INOUT)  :: SQt(:)     ! Std deviation of position
   COMPLEX(kind=Rkind),  INTENT(INOUT)  :: At(:)      ! "A" parameter (complex width)
   REAL(kind=Rkind),     INTENT(INOUT)  :: Pt(:)      ! Momentum expectation value <P>
   TYPE(propa_t), OPTIONAL, INTENT(IN)  :: propa      ! Propagation settings (optional)
   LOGICAL, OPTIONAL, INTENT(IN) :: need_Pt, need_At

   ! === Local variables ===
   INTEGER, ALLOCATABLE                :: Tab_Iq(:, :)
   LOGICAL :: compute_Pt, compute_At

   ! ---------------------------------------------------------------
   ! Build integration index table (used for <ψ|Q|ψ>, <ψ|Q²|ψ>, etc.)
   ! ---------------------------------------------------------------
   CALL Calc_tab_Iq0(Tab_Iq, psi%Basis)

   ! ---------------------------------------------------------------
   ! Compute expectation values of Q and its standard deviation
   ! ---------------------------------------------------------------
   CALL Calc_AVQ_SQ_nD(psi, Qt, SQt, Tab_Iq)

   compute_Pt = .TRUE.
   compute_At = .TRUE.
   IF (PRESENT(propa)) THEN
      compute_Pt = propa%P
      compute_At = propa%Beta
   END IF
   IF (PRESENT(need_Pt)) compute_Pt = need_Pt
   IF (PRESENT(need_At)) compute_At = need_At

   Pt(:) = ZERO
   At(:) = SQt(:) * SQt(:)
   IF (compute_Pt) CALL Calc_Av_imp_k_nD(psi, Pt)
   IF (compute_At) CALL Calc_Avg_A_nD(psi, At)

   ! ---------------------------------------------------------------
   ! Cleanup
   ! ---------------------------------------------------------------
   DEALLOCATE(Tab_Iq)

END SUBROUTINE Calc_Basis_parameters

SUBROUTINE Calc_S(S,nb,nq,Q0,SQ0,A0,P0,Qt,SQt,At,Pt)

  IMPLICIT NONE
  ! Input/Output
  complex(kind=Rkind), intent(inout) :: S(:,:)

  ! Input parameters
  complex(kind=Rkind), intent(in)    :: A0, At
  real(kind=Rkind), intent(in)       :: Q0, SQ0, P0
  real(kind=Rkind), intent(in)       :: Qt, SQt, Pt
  integer, intent(in)                :: nb, nq

  ! Local variables
  integer :: ib, iq
  real(kind=Rkind) :: SQeq, Qeq
  real(kind=Rkind), allocatable :: Q(:), W(:)
  real(kind=Rkind) :: B0, Bt
  complex(kind=Rkind), allocatable :: d0gb0(:,:),d0gb1(:,:), d0bgw(:,:)

  !--------------------------------------------------------------
  ! Compute effective scaling factor SQeq and center Qeq
  !--------------------------------------------------------------
  SQeq = sqrt(SQ0*SQ0 + SQt*SQt) / sqrt(TWO)
  Qeq  = (SQ0*SQ0*Q0 + SQt*SQt*Qt) / (SQ0*SQ0 + SQt*SQt)

  !--------------------------------------------------------------
  ! Extract imaginary parts of At and A0 (b_t and b_0)
  !--------------------------------------------------------------
  Bt = AIMAG(At)
  B0 = AIMAG(A0)

  !--------------------------------------------------------------
  ! Allocate grid points and weights for Hermite quadrature
  !--------------------------------------------------------------
  allocate(Q(nq), W(nq))
  call hercom(nq, Q(:), W(:))

  ! Scale quadrature weights and shift grid points
  W(:) = W(:) / SQeq
  Q(:) = Qeq + Q(:) / SQeq

  !--------------------------------------------------------------
  ! Allocate basis function arrays
  !--------------------------------------------------------------
  allocate(d0gb0(nq, nb))
  allocate(d0gb1(nq, nb))
  allocate(d0bgw(nb, nq))

  !--------------------------------------------------------------
  ! Debug prints for checking parameters
  !--------------------------------------------------------------
  !print*, 'Q0, A0, P0 =', Q0, A0, P0
  !print*, 'Qt, At, Pt =', Qt, At, Pt
  !print*, 'SQeq, Qeq =', SQeq, Qeq

  !--------------------------------------------------------------
  ! Compute Hagedorn basis functions at grid points Q(:)
  ! for both initial (Q0, A0, P0) and target (Qt, At, Pt)
  !--------------------------------------------------------------
  DO iq = 1, nq
      DO ib = 1, nb
          call d0poly_Hermite_exp_cplx(Q(iq), Q0, A0, P0, ib-1, d0gb0(iq, ib))
          call d0poly_Hermite_exp_cplx(Q(iq), Qt, At, Pt, ib-1, d0gb1(iq, ib))
      END DO
  END DO

  d0bgw = transpose(d0gb1)
   DO ib = 1, nb
      d0bgw(ib, :) = d0bgw(ib, :)*W(:)
   END DO

  !--------------------------------------------------------------
  ! Compute overlap matrix S = <d0bgw|d0gb>
  !--------------------------------------------------------------
  S = MATMUL(conjg(d0bgw), d0gb0)

  !--------------------------------------------------------------
  ! Print overlap matrix S for verification
  !--------------------------------------------------------------
  !call Write_VecMat(S, out_unit, 5, info='S')

  !--------------------------------------------------------------
  ! Deallocate arrays
  !--------------------------------------------------------------
  deallocate(d0gb0, d0gb1, d0bgw, W, Q)

End SUBROUTINE

SUBROUTINE Test_Calc_S()
  USE QDUtil_m
  IMPLICIT NONE

  integer, parameter :: nb = 3
  complex(kind=Rkind) :: S(nb, nb)
  real(kind=Rkind) :: Q0, SQ0, P0, Qt, SQt, Pt
  complex(kind=Rkind) :: A0, At

  !--------------------------------------------------------------
  ! Initialize Parameters for the Test
  !--------------------------------------------------------------

  A0  = CMPLX(1.0_Rkind, 0.1_Rkind, kind=Rkind)
  Q0  = 0.5_Rkind
  P0  = 0.25_Rkind
  SQ0 = 1.0_Rkind

  At  = CMPLX(1.0_Rkind, 0.2_Rkind, kind=Rkind)
  Qt  = 0.49_Rkind
  SQt = 1.0_Rkind
  Pt  = 0.24_Rkind
  print*, 'Q0, SQ0, A0, P0 =', Q0, SQ0, A0, P0
  print*, 'Qt, SQt, At, Pt =', Qt, SQt, At, Pt

  !--------------------------------------------------------------
  ! Call Calc_S to compute the overlap matrix S
  !--------------------------------------------------------------
  !Using algebraic formula
  call Calc_S_algebraic_genfunc(S, nb, Q0, SQ0, A0, P0, Qt, SQt, At, Pt)
  !call Calc_S(S, nb, nq, Q0, SQ0, A0, P0, Qt, SQt, At, Pt)
    call  Write_VecMat(S, out_unit, nb,  info='S')
  !--------------------------------------------------------------
  ! End of test subroutine
  !--------------------------------------------------------------
  print*, 'Test_Calc_S finished.'

END SUBROUTINE Test_Calc_S


  subroutine Test_Compare_Overlap_ThreeMethods()

    integer, parameter :: nb = 10

    integer, parameter :: n_scenarios = 6

    integer, parameter :: nq_values(4) = [10, 25, 50, 100]

    real(kind=Rkind), parameter :: tol_exact = 1.0e-12_Rkind

    real(kind=Rkind), parameter :: tol_num   = 1.0e-10_Rkind



    type :: Scenario_t

      real(kind=Rkind) :: Q0, SQ0, P0, Qt, SQt, Pt

      complex(kind=Rkind) :: A0, At

    end type Scenario_t



    type(Scenario_t) :: scenarios(n_scenarios)

    complex(kind=Rkind) :: Snum(nb,nb), Salg(nb,nb), Sbourq(nb,nb)

    real(kind=Rkind) :: e_num_alg, e_num_bq, e_alg_bq

    integer :: i, iq

    logical :: all_ok



   ! Scénario 1 (déjà cohérent)

  scenarios(1) = Scenario_t(0.5_Rkind, 1.0_Rkind, 0.25_Rkind, &

                            0.49_Rkind, 1.0_Rkind, 0.24_Rkind, &

                            cmplx(1.0_Rkind, 0.1_Rkind, kind=Rkind), &

                            cmplx(1.0_Rkind, 0.2_Rkind, kind=Rkind))



! Scénario 2 : 

  scenarios(2) = Scenario_t(0.3_Rkind, 1.0_Rkind, 0.1_Rkind, &

                            0.35_Rkind, 1.0_Rkind, 0.12_Rkind, &

                            cmplx(1.0_Rkind, 0.0_Rkind, kind=Rkind), &

                            cmplx(1.0_Rkind, 0.0_Rkind, kind=Rkind))



! Scénario 3 : 

  scenarios(3) = Scenario_t(0.2_Rkind, sqrt(0.5_Rkind), 0.3_Rkind, &

                            0.25_Rkind, sqrt(0.6_Rkind), 0.28_Rkind, &

                            cmplx(0.5_Rkind, 2.0_Rkind, kind=Rkind), &

                            cmplx(0.6_Rkind, 1.5_Rkind, kind=Rkind))



! Scénario 4 (déjà cohérent)

  scenarios(4) = Scenario_t(0.0_Rkind, 1.0_Rkind, 0.0_Rkind, &

                            0.1_Rkind, 1.0_Rkind, 0.1_Rkind, &

                            cmplx(1.0_Rkind, 0.0_Rkind, kind=Rkind), &

                            cmplx(1.0_Rkind, 0.0_Rkind, kind=Rkind))



! Scénario 5 (déjà cohérent)

  scenarios(5) = Scenario_t(0.0_Rkind, 1.0_Rkind, 0.0_Rkind, &

                            2.0_Rkind, 1.0_Rkind, 0.0_Rkind, &

                            cmplx(1.0_Rkind, 0.0_Rkind, kind=Rkind), &

                            cmplx(1.0_Rkind, 0.0_Rkind, kind=Rkind))



! Scénario 6 : 

  scenarios(6) = Scenario_t(0.7_Rkind, sqrt(0.8_Rkind), 0.4_Rkind, &

                            0.6_Rkind, sqrt(0.9_Rkind), 0.35_Rkind, &

                            cmplx(0.8_Rkind, 0.3_Rkind, kind=Rkind), &

                            cmplx(0.9_Rkind, 0.4_Rkind, kind=Rkind))



    write(out_unit,'(A)') ''

    write(out_unit,'(A)') '=== Test de comparaison sur trois méthodes ==='

    write(out_unit,'(A)') '1) Calc_S_numeric_realgrid  : quadrature reelle historique'

    write(out_unit,'(A)') '2) Calc_S_algebraic_genfunc : formule algebrique exacte (fonction generatrice)'

    write(out_unit,'(A)') '3) Calc_S_Bourquin522       : quadrature complexe exacte de type Bourquin 5.2.2'

    write(out_unit,'(A,I0)') 'nb = ', nb

    write(out_unit,'(A,4(I0,1X))') 'nq testes pour la methode numerique : ', nq_values

    write(out_unit,'(A)') ''



    all_ok = .true.

    do i = 1, n_scenarios

      call Calc_S_algebraic_genfunc(Salg, nb, scenarios(i)%Q0, scenarios(i)%SQ0, scenarios(i)%A0, scenarios(i)%P0, &

                                         scenarios(i)%Qt, scenarios(i)%SQt, scenarios(i)%At, scenarios(i)%Pt)

      call Calc_S_Bourquin522(Sbourq, nb, scenarios(i)%Q0, scenarios(i)%SQ0, scenarios(i)%A0, scenarios(i)%P0, &

                                   scenarios(i)%Qt, scenarios(i)%SQt, scenarios(i)%At, scenarios(i)%Pt, nq_exact=nb)



      e_alg_bq = max_abs_c(Salg - Sbourq)

      if (e_alg_bq > tol_exact) all_ok = .false.



      write(out_unit,'(A,I0,A)') 'Scenario ', i, ' :'

      write(out_unit,'(A,2(" (",F10.6,",",F10.6,")"))') '  A0 / At =', real(scenarios(i)%A0), aimag(scenarios(i)%A0), &

                                                            real(scenarios(i)%At), aimag(scenarios(i)%At)

      write(out_unit,'(A,3F12.6)') '  Q0, SQ0, P0 =', scenarios(i)%Q0, scenarios(i)%SQ0, scenarios(i)%P0

      write(out_unit,'(A,3F12.6)') '  Qt, SQt, Pt =', scenarios(i)%Qt, scenarios(i)%SQt, scenarios(i)%Pt

      write(out_unit,'(A,ES12.4)') '  ||S_alg - S_bourquin||_max = ', e_alg_bq

      write(out_unit,'(A)') '  Convergence de la methode numerique reelle :'



      do iq = 1, size(nq_values)

        call Calc_S_numeric_realgrid(Snum, nb, nq_values(iq), scenarios(i)%Q0, scenarios(i)%SQ0, scenarios(i)%A0, scenarios(i)%P0, &

                                               scenarios(i)%Qt, scenarios(i)%SQt, scenarios(i)%At, scenarios(i)%Pt)

        e_num_alg = max_abs_c(Snum - Salg)

        e_num_bq  = max_abs_c(Snum - Sbourq)

        write(out_unit,'(A,I4,A,ES12.4,A,ES12.4)') '    nq =', nq_values(iq), '  |num-alg| =', e_num_alg, '  |num-bq| =', e_num_bq

      end do

      ! Make the highest-resolution numerical comparison an assertion for
      ! every scenario (the previous test only checked the final scenario).
      call Calc_S_numeric_realgrid(Snum, nb, 100, scenarios(i)%Q0, scenarios(i)%SQ0, scenarios(i)%A0, scenarios(i)%P0, &
                                   scenarios(i)%Qt, scenarios(i)%SQt, scenarios(i)%At, scenarios(i)%Pt)
      e_num_alg = max_abs_c(Snum - Salg)
      e_num_bq  = max_abs_c(Snum - Sbourq)
      if (e_num_alg > tol_num .or. e_num_bq > tol_num) all_ok = .false.

      write(out_unit,'(A)') ''

    end do



    if (all_ok) then

      write(out_unit,'(A)') 'RESULTAT GLOBAL : OK - les deux methodes exactes coincident au bruit machine,'

      write(out_unit,'(A)') 'et la methode numerique reelle converge vers elles quand nq augmente.'

    else

      write(out_unit,'(A)') 'RESULTAT GLOBAL : ECHEC - verifier les conventions de phase ou les parametrisations.'

      error stop 'Test_Compare_Overlap_ThreeMethods failed.'

    end if

  end subroutine Test_Compare_Overlap_ThreeMethods


subroutine build_overlap_matrices(B1, B2)
  implicit none

  ! Input/output
  TYPE(Basis_t), intent(inout)         :: B1, B2

  ! Local variables
  integer                              :: ib, nb, nq, ndim
  character(len=10)                    :: bname
  real(kind=Rkind)                     :: Q0, SQ0, P0, Qt, SQt, Pt
  complex(kind=Rkind)                  :: A0, At

  ! Determine number of dimensions (suppose même taille pour B1 et B2)
  ndim = size(B1%tab_basis)

  ! Loop over each coordinate
  do ib = 1, ndim
     bname = B1%tab_basis(ib)%Basis_name

     ! Extract parameters of basis B1
     nb   = B1%tab_basis(ib)%nb
     nq   = B1%tab_basis(ib)%nq
     Q0   = B1%tab_basis(ib)%Q0
     SQ0  = B1%tab_basis(ib)%scaleQ
     A0   = B1%tab_basis(ib)%alpha
     P0   = B1%tab_basis(ib)%Imp_k

     ! Extract parameters of basis B2
     Qt   = B2%tab_basis(ib)%Q0
     SQt  = B2%tab_basis(ib)%scaleQ
     At   = B2%tab_basis(ib)%alpha
     Pt   = B2%tab_basis(ib)%Imp_k
      
     print*, 'Q0,Qt', Q0,Qt
     print*, 'SQ0,SQt', SQ0,SQt
     print*, 'A0,At', A0,At
     print*, 'P0,Pt', P0,Pt
     print*, 'nb,nq', nb,nq 
     ! Build S according to the basis type
     if (trim(bname) == 'hag' .or. trim(bname) == 'ho') then
          !Using algebric formula
        call Calc_S_algebraic_genfunc(B1%tab_basis(ib)%S, nb, Q0, SQ0, A0, P0, Qt, SQt, At, Pt)
        !call Calc_S(B1%tab_basis(ib)%S, nb, nq, Q0, SQ0, A0, P0, Qt, SQt, At, Pt)
     else
        call build_identity_matrix(B1%tab_basis(ib)%S, nb)
     end if
     call Write_VecMat(B1%tab_basis(ib)%S, out_unit, nb, info='S')
  end do

end subroutine build_overlap_matrices

subroutine Test_projection_hagedorn(Basis,psi)
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  implicit none

  TYPE(Basis_t), intent(in)                 :: Basis
  TYPE(psi_t), intent(in)                   :: psi

  ! Local variables
  TYPE(psi_t)                               :: psi2, psi1
  TYPE(Basis_t)                             :: B1, B2
  type(Op_t)                                :: H1,H2
  integer, allocatable                      :: Tab_iq1(:,:)
  integer, allocatable                      :: Tab_iq2(:,:)
  
  real(kind=Rkind)                          :: norm1  , norm2, E1, E2
  integer                                   :: ib
  complex(kind=Rkind)                       :: alpha(2)
  real(kind=Rkind)                          :: q(2), p(2), sq(2)

  ! Initialization of parameters
  alpha = [1.0_Rkind + 0.1_Rkind*EYE, 1.0_Rkind + 0.1_Rkind*EYE]
  q     = [1.9_Rkind, 0.1_Rkind]
  p     = [0.1_Rkind, 0.1_Rkind]
  sq    = [1.0_Rkind, 1.0_Rkind]



  ! Copy Basis into B1 and B2
  call init_Basis1_TO_Basis2(B1, Basis)
  call init_Basis1_TO_Basis2(B2, Basis)

  ! Modify the parameters of B2
  B2%tab_basis(1)%alpha    = alpha(1)
  B2%tab_basis(2)%alpha    = alpha(2)
  B2%tab_basis(1)%q0       = q(1)
  B2%tab_basis(2)%q0       = q(2)
  B2%tab_basis(1)%Imp_k    = p(1)
  B2%tab_basis(2)%Imp_k    = p(2)
  B2%tab_basis(1)%scaleQ   = sq(1)
  B2%tab_basis(2)%scaleQ   = sq(2)

  ! Construct both primitive bases
  call construct_primitive_basis(B1)
  call construct_primitive_basis(B2)

  call Calc_tab_Iq0(Tab_Iq1,B1)
  call Calc_tab_Iq0(Tab_Iq2,B2)

  call Set_Op(H1, B1,Tab_Iq1,nio=in_unit)
  call Set_Op(H2, B2,Tab_Iq2,nio=in_unit)

  ! Initialize psi2 with B2 and psi1 with B1
  call init_psi(psi2, B2, cplx=.true., grid=.false.)
  call init_psi(psi1, B1, cplx=.true., grid=.false.)

  call build_overlap_matrices(B1,B2)

  ! Fill psi1 with CZERO except for the first coefficient
  psi1%CVec(:) = psi%CVec(:)


  ! Initialize psi2 to zero
  psi2%CVec(:) = CZERO

  ! Print the parameters of both bases for verification
  print *, 'B1'
  print *, "Basis is allocated", Basis_IS_allocated(B1)
  print *, "B1%tab_basis(1)%alpha(:)",    B1%tab_basis(1)%alpha
  print *, "B1%tab_basis(2)%alpha(:)",    B1%tab_basis(2)%alpha
  print *, "B1%tab_basis(1)%q0(:)",       B1%tab_basis(1)%q0
  print *, "B1%tab_basis(2)%q0(:)",       B1%tab_basis(2)%q0
  print *, "B1%tab_basis(1)%imp_k(:)",    B1%tab_basis(1)%Imp_k
  print *, "B1%tab_basis(2)%imp_k(:)",    B1%tab_basis(2)%Imp_k
  print *, "B1%tab_basis(1)%scaleQ(:)",  B1%tab_basis(1)%scaleQ
  print *, "B1%tab_basis(2)%scaleQ(:)",  B1%tab_basis(2)%scaleQ

  print *, 'B2'
  print *, "Basis is allocated", Basis_IS_allocated(B2)
  print *, "B2%tab_basis(1)%alpha(:)",    B2%tab_basis(1)%alpha
  print *, "B2%tab_basis(2)%alpha(:)",    B2%tab_basis(2)%alpha
  print *, "B2%tab_basis(1)%q0(:)",       B2%tab_basis(1)%q0
  print *, "B2%tab_basis(2)%q0(:)",       B2%tab_basis(2)%q0
  print *, "B2%tab_basis(1)%imp_k(:)",    B2%tab_basis(1)%Imp_k
  print *, "B2%tab_basis(2)%imp_k(:)",    B2%tab_basis(2)%Imp_k
  print *, "B2%tab_basis(1)%scaleQ(:)",  B2%tab_basis(1)%scaleQ
  print *, "B2%tab_basis(2)%scaleQ(:)",  B2%tab_basis(2)%scaleQ


   ! Compute the norm and energy of psi1
   call Calc_Av_E(E1, psi1, H1)
   call Calc_Norm_OF_psi(psi1, norm1)
   print *, 'Norm, E1', norm1, E1
  ! Display psi1 coefficients
  print *, '--------------------------------------------------------'
  do ib = 1, size(psi1%CVec)
     write(out_unit, *) ib, abs(psi1%CVec(ib))**2
  end do


  ! Project psi1 (constructed with B1) onto B2; result is stored in psi2
  call Projection(psi2, psi1)

  ! Display the coefficients of psi2
  do ib = 1, size(psi2%CVec)
     write(out_unit, *) ib, abs(psi2%CVec(ib))**2
  end do

  ! Compute the norm and energy of psi2
  call Calc_Av_E(E2, psi2, H2)
  call Calc_Norm_OF_psi(psi2, norm2)
  print *, 'Norm, E2', norm2, E2

  IF (.NOT. IEEE_IS_FINITE(norm1) .OR. .NOT. IEEE_IS_FINITE(E1) .OR. &
      .NOT. IEEE_IS_FINITE(norm2) .OR. .NOT. IEEE_IS_FINITE(E2)) &
     ERROR STOP 'Test_projection_hagedorn: non-finite norm or energy after projection.'
  IF (norm2 <= SQRT(TINY(ONE))) &
     ERROR STOP 'Test_projection_hagedorn: projected state has zero/vanishing norm.'
  IF (.NOT. ALL(IEEE_IS_FINITE(REAL(psi2%CVec,KIND=Rkind))) .OR. &
      .NOT. ALL(IEEE_IS_FINITE(AIMAG(psi2%CVec)))) &
     ERROR STOP 'Test_projection_hagedorn: projected coefficients contain non-finite values.'

  ! End of projection test
end subroutine test_projection_hagedorn

 
 SUBROUTINE Hagedorn(psi_main, psi_temp, propa, use_psi0_params)
   !============================================================================
   ! PURPOSE:
   !   Update a wavefunction in a Hagedorn basis by constructing a new non-variational
   !   Hagedorn basis and projecting the temporary wavefunction onto it.
   !
   !   - If use_psi0_params = .false., the basis parameters are computed from psi_temp
   !     (forward propagation, like Hagedorn).
   !   - If use_psi0_params = .true., the basis parameters are taken from psi_main
   !     (inverse projection, like Hagedorn_inv).
   !
   ! INPUT/OUTPUT:
   !   psi_main          : main wavefunction to update (INOUT)
   !   psi_temp          : temporary wavefunction from propagation (INOUT)
   !   propa             : propagation parameters (IN), contains renorm and calculation options
   !   use_psi0_params   : logical, TRUE => take parameters from psi_main (inverse)
   !
   !============================================================================
   USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
   IMPLICIT NONE
 
   TYPE(psi_t),     INTENT(INOUT) :: psi_main
   TYPE(psi_t),     INTENT(INOUT) :: psi_temp
   TYPE(propa_t),   INTENT(IN)    :: propa
   LOGICAL,         INTENT(IN)    :: use_psi0_params
 
   ! Local Gaussian/Hagedorn parameters
   COMPLEX(KIND=Rkind), ALLOCATABLE :: At(:)
   REAL(KIND=Rkind), ALLOCATABLE    :: Qt(:), SQt(:), Pt(:)
   REAL(KIND=Rkind)                  :: Norm0, norm
   INTEGER                           :: ndim, idim
 
   !-------------------------------- Initialization --------------------------------
   IF (.NOT. ASSOCIATED(psi_main%Basis)) ERROR STOP 'Hagedorn: psi_main has no associated basis.'
   IF (.NOT. ASSOCIATED(psi_temp%Basis)) ERROR STOP 'Hagedorn: psi_temp has no associated basis.'
   IF (.NOT. ALLOCATED(psi_main%CVec) .OR. .NOT. ALLOCATED(psi_temp%CVec)) &
      ERROR STOP 'Hagedorn: coefficient vectors must be allocated.'
   IF (.NOT. ALLOCATED(psi_main%Basis%tab_basis) .OR. .NOT. ALLOCATED(psi_temp%Basis%tab_basis)) &
      ERROR STOP 'Hagedorn: primitive basis metadata must be allocated.'
   IF (SIZE(psi_main%CVec) /= SIZE(psi_temp%CVec)) &
      ERROR STOP 'Hagedorn: source and destination coefficient dimensions differ.'
   IF (SIZE(psi_main%Basis%tab_basis) /= SIZE(psi_temp%Basis%tab_basis)) &
      ERROR STOP 'Hagedorn: source and destination basis dimensions differ.'

   ndim = SIZE(psi_main%Basis%tab_basis) - 1
   IF (ndim < 1) ERROR STOP 'Hagedorn: at least one nuclear degree of freedom is required.'
   ALLOCATE(Qt(ndim), SQt(ndim), Pt(ndim), At(ndim))
 
   !-------------------------------- Compute Parameters --------------------------------
   IF (.NOT. use_psi0_params) THEN
      ! Forward propagation: compute from psi_temp
      CALL Calc_Basis_parameters(psi_temp, Qt, SQt, At, Pt, propa)
   ELSE
      ! Inverse projection: copy from psi_main
      CALL Get_Basis_Parameters(psi_main%Basis, Qt, SQt, At, Pt)
   END IF

   DO idim = 1, ndim
      IF (.NOT. IEEE_IS_FINITE(Qt(idim)) .OR. .NOT. IEEE_IS_FINITE(SQt(idim)) .OR. &
          .NOT. IEEE_IS_FINITE(Pt(idim)) .OR. .NOT. IEEE_IS_FINITE(REAL(At(idim),KIND=Rkind)) .OR. &
          .NOT. IEEE_IS_FINITE(AIMAG(At(idim)))) ERROR STOP 'Hagedorn: non-finite packet parameter.'
      IF (SQt(idim) <= SQRT(TINY(ONE)) .OR. REAL(At(idim),KIND=Rkind) <= ZERO) &
         ERROR STOP 'Hagedorn: packet width must be strictly positive.'
   END DO
 
   !------------------------ Construct new Hagedorn basis ---------------------------
   CALL Construct_Hagedorn_none_Variational_Basis(psi_temp%Basis, Qt, SQt, At, Pt)

   ! Projection coefficients are expressed in the new basis. Do not rely on
   ! psi_main and psi_temp having happened to share a basis pointer.
   psi_main%Basis => psi_temp%Basis
 
   !--------------------------- Save norm before projection --------------------------
   IF (propa%debug) CALL Calc_Norm_OF_psi(psi_temp, Norm0)
 
   !-------------------------- Project temporary wavefunction ------------------------
   CALL Projection(psi_main, psi_temp)
 
   !----------------------------- Optional renormalization ---------------------------
   IF (propa%renorm) THEN
      CALL Calc_Norm_OF_psi(psi_main, norm)
      psi_main%CVec(:) = psi_main%CVec(:) / norm
   END IF
 
   !----------------------------- Final norm check -----------------------------------
   IF (propa%debug) THEN
      CALL Calc_Norm_OF_psi(psi_main, norm)
      WRITE(out_unit,*) 'Norm diff =', ABS(Norm0 - norm), 'Norm0 =', Norm0, 'Norm =', norm
   END IF
 
   !----------------------------- Clean-up -------------------------------------------
   DEALLOCATE(Qt, SQt, At, Pt)
 
 END SUBROUTINE Hagedorn

 SUBROUTINE march_Global(psi, psi_dt, t, propa, H, u_Norm_13, u_E_13)
  IMPLICIT NONE

  ! Input/output arguments
  TYPE(propa_t), INTENT(IN)     :: propa     ! Propagation parameters
  TYPE(psi_t),    INTENT(INOUT) :: psi       ! Current wavefunction (updated if needed)
  TYPE(psi_t),    INTENT(INOUT) :: psi_dt    ! Temporary wavefunction (used for propagation)
  REAL(kind=Rkind), INTENT(IN)  :: t         ! Current time
  TYPE(Op_t),     INTENT(INOUT) :: H         ! Hamiltonian operator (includes potential, kinetic, etc.)
  INTEGER,        INTENT(IN)    :: u_Norm_13, u_E_13 ! Automatically allocated output units

  ! Local variables
  COMPLEX(kind=Rkind), ALLOCATABLE :: At(:)      ! Complex width parameters for each degree of freedom
  REAL(kind=Rkind),    ALLOCATABLE :: Qt(:)      ! Gaussian centers in position
  REAL(kind=Rkind),    ALLOCATABLE :: SQt(:)     ! Widths (real part of square root of At)
  REAL(kind=Rkind),    ALLOCATABLE :: Pt(:)      ! Gaussian centers in momentum
  INTEGER, ALLOCATABLE             :: Tab_iq(:, :) ! Index mapping for operator construction

  REAL(kind=Rkind) :: Norm0, norm      ! Norms before and after propagation (for conservation check)
  REAL(kind=Rkind) :: E0, E            ! Energies before and after propagation
  INTEGER          :: ndim             ! Number of degrees of freedom (basis size - 1)
  LOGICAL :: debug
  CHARACTER(len=Name_len)                  :: name

  ! Validate required associations before dereferencing metadata, regardless of
  ! the debug setting. This turns invalid state into a useful diagnostic rather
  ! than a segmentation fault in a later basis/projection routine.
  IF (.NOT. ASSOCIATED(psi%Basis)) ERROR STOP 'march_Global: psi has no associated basis.'
  IF (.NOT. ASSOCIATED(psi_dt%Basis)) ERROR STOP 'march_Global: psi_dt has no associated basis.'
  IF (.NOT. ALLOCATED(psi%CVec) .OR. .NOT. ALLOCATED(psi_dt%CVec)) &
     ERROR STOP 'march_Global: coefficient vectors must be allocated.'
  IF (.NOT. ALLOCATED(psi%Basis%tab_basis) .OR. .NOT. ALLOCATED(psi_dt%Basis%tab_basis)) &
     ERROR STOP 'march_Global: basis metadata must be allocated.'
  IF (SIZE(psi%CVec) /= SIZE(psi_dt%CVec)) &
     ERROR STOP 'march_Global: psi and psi_dt coefficient dimensions differ.'
  IF (SIZE(psi%Basis%tab_basis) < 2 .OR. SIZE(psi_dt%Basis%tab_basis) < 2) &
     ERROR STOP 'march_Global: basis must contain nuclear and electronic components.'

  ! Determine number of degrees of freedom from the basis
  debug = propa%debug
  ndim = SIZE(psi%Basis%tab_basis) - 1

  IF (debug) THEN
     WRITE(out_unit,*) 'BEGIN march_Global'
     WRITE(out_unit,*) 'ndim = ', ndim
     IF(.NOT. ASSOCIATED(psi%Basis)) THEN
        WRITE(out_unit,*) 'psi%Basis is not associated in march_Global routine'
        STOP 'psi%Basis is not associated in march_Global routine'
     END IF
     IF(.NOT. ASSOCIATED(psi_dt%Basis)) THEN
        WRITE(out_unit,*) 'psi_dt%Basis is not associated in march_Global routine'
        STOP 'psi_dt%Basis is not associated in march_Global routine'
     END IF
     IF(.NOT. ALLOCATED(psi%CVec)) THEN
        WRITE(out_unit,*) 'psi%CVec is not allocated in march_Global routine'
        STOP 'psi%CVec is not allocated in march_Global routine'
     ELSE
        WRITE(out_unit,*) 'psi%CVec is allocated in march_Global routine'
        WRITE(out_unit,*) 'size(psi%CVec) = ', size(psi%CVec)
        FLUSH(out_unit)
     END IF
     IF(.NOT. ALLOCATED(psi_dt%CVec)) THEN
        WRITE(out_unit,*) 'psi_dt%CVec is not allocated in march_Global routine'
        STOP 'psi_dt%CVec is not allocated in march_Global routine'
     ELSE
        WRITE(out_unit,*) 'psi_dt%CVec is allocated in march_Global routine'
        WRITE(out_unit,*) 'size(psi_dt%CVec) = ', size(psi_dt%CVec)
        FLUSH(out_unit)
     END IF
     FLUSH(out_unit)
  END IF

  ! Allocate and initialize Gaussian parameters
  ALLOCATE(Qt(ndim), SQt(ndim), Pt(ndim), At(ndim))
  Qt(:)  = ZERO         ! Initial position centers
  SQt(:) = ONE          ! Initial widths (sqrt of Re(At))
  Pt(:)  = ZERO         ! Initial momenta
  At(:)  = ONE          ! Initial complex widths (simplified to 1 for default)

  ! If the propagator is Hagedorn-based, perform full Hagedorn propagation
  IF (propa%propa_name == 'hagedorn') THEN
     ! Three-phase Hagedorn update:
     ! 1) Propagate psi -> psi_dt with H and the old basis.
     CALL march(psi, psi_dt, t, propa, H)

     ! Compute diagnostics only when their corresponding output is requested.
     IF (propa%write_E_13) CALL Calc_Av_E(E0, psi_dt, H)
     IF (propa%write_Norm_13) CALL Calc_Norm_OF_psi(psi_dt, Norm0)

     ! 2) Hagedorn derives the packet parameters from psi_dt, rebuilds the
     !    moving basis, and constructs the old-to-new overlap matrices.
     ! 3) It projects psi_dt into that new basis; then rebuild H for the next
     !    time step. The overlap S is <new_basis|old_basis>, applied to old
     !    coefficient vectors from the left along every nuclear dimension.
     CALL Hagedorn(psi, psi_dt, propa,use_psi0_params=.FALSE.)

     ! Reconstruct the operator with the newly installed basis.
     IF (ALLOCATED(H%Scal_pot)) DEALLOCATE(H%Scal_pot)
     CALL Calc_tab_Iq(Tab_iq, psi%Basis)
     CALL Set_Op(H, psi%Basis, Tab_iq,nio=in_unit)

     ! Recalculate and write only the requested diagnostics.
     IF (propa%write_E_13) THEN
        CALL Calc_Av_E(E, psi, H)
        WRITE(u_E_13,*) t, ABS(E0 - E), E0, E
        FLUSH(u_E_13)
     END IF
     IF (propa%write_Norm_13) THEN
        CALL Calc_Norm_OF_psi(psi, norm)
        WRITE(u_Norm_13,*) t, ABS(Norm0 - norm), Norm0, norm
        FLUSH(u_Norm_13)
     END IF
  ELSE
     ! Default propagation: simple time step without Hagedorn basis
     CALL march(psi, psi_dt, t, propa, H)

     ! Directly update psi with temporary psi_dt
     !stop 'Hagedorn basis not implemented'
     IF(debug) THEN
     WRITE(out_unit,*) 'size(psi_dt%CVec), size(psi%CVec) =:', size(psi_dt%CVec), size(psi%CVec)
     WRITE(out_unit,*) 'size(psi_dt%CVec) == size(psi%CVec)=:', size(psi_dt%CVec)==size(psi%CVec)
     IF(.NOT. size(psi_dt%CVec)==size(psi%CVec)) WRITE(out_unit,*) 'psi_dt%CVec and psi%CVec  In march_Global have different sizes'
     END IF
      
     name = propa%propa_name2 
     CALL string_uppercase_TO_lowercase(name)
     IF (name=='vp')THEN
       CALL Set_Op(H, psi_dt%Basis, Tab_Iq,alloc=.FALSE.,nio=in_unit)
      CALL init_psi(psi,psi_dt%Basis,cplx=.TRUE.,grid=.FALSE.)
     END IF 
     psi%CVec(:) = psi_dt%CVec(:)
     FLUSH(out_unit)
  END IF

END SUBROUTINE march_Global

 SUBROUTINE d0poly_Hermite_exp(Q,ib,d0)
    USE QDUtil_m
    IMPLICIT NONE
    complex(kind=Rkind),intent(inout) :: d0
    real(kind=Rkind),intent(in)       :: Q
    integer,intent(in)                :: ib
   
      d0 = poly_Hermite(Q,ib)*exp(-HALF*Q*Q)
    RETURN
 END SUBROUTINE


 SUBROUTINE d0poly_Hermite_exp_cplx(Q, Qt, At, Pt, ib, d0)
  USE QDUtil_m
  IMPLICIT NONE

  ! Input/output arguments
  complex(kind=Rkind), INTENT(INOUT) :: d0           ! Output: result of Hermite polynomial * exponential (complex)
  real(kind=Rkind),    INTENT(IN)    :: Q            ! Position variable
  real(kind=Rkind),    INTENT(IN)    :: Qt           ! Center of the Gaussian
  real(kind=Rkind),    INTENT(IN)    :: Pt           ! Momentum
  complex(kind=Rkind), INTENT(IN)    :: At           ! Complex width parameter 
  integer,            INTENT(IN)     :: ib           ! Hermite polynomial order

  ! Local variables
  complex(kind=Rkind) :: w                          ! Complex phase factor
  real(kind=Rkind)    :: SQ                         ! Square root of the real part of At
  real(kind=Rkind)    :: DQ                         ! Scaled displacement from Qt
  real(kind=Rkind)    :: Bt                         ! Imaginary part of At

  ! Compute square root of the real part of At (assumed to be positive)
  SQ = sqrt(real(At, kind=Rkind))

  ! Compute scaled coordinate (DQ = sqrt(Re(At)) * (Q - Qt))
  DQ = SQ * (Q - Qt)

  ! Extract the imaginary part of At
  Bt = aimag(At)

  ! Compute the complex exponential factor:
  !   exp[-(i/2) * Bt * (Q - Qt)^2 + i * Pt * (Q - Qt)]
  w = exp(-HALF * EYE * Bt * (Q - Qt)**2 + EYE * Pt * (Q - Qt))

  ! Call routine to compute the Hermite polynomial times Gaussian part
  CALL d0poly_Hermite_exp(DQ, ib, d0)

  ! Multiply the result by sqrt(SQ) and the complex exponential
  d0 = sqrt(SQ) * d0 * w

  ! Optional: write debug output
  ! WRITE(out_unit,*) 'l x,p,beta, d0 :', ib, Q, Pt, Bt, d0

  RETURN
END SUBROUTINE d0poly_Hermite_exp_cplx


end module Hagedorn_m
