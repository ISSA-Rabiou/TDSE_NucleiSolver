MODULE Sub_Vp_m
   USE QDUtil_m
   USE NDindex_m
   USE op_m
   USE Basis_m
   USE psi_m
   IMPLICIT NONE

   PRIVATE
   PUBLIC :: Psi_Lambda_Transform, Test_Psi_Lambda_Transform
   PUBLIC :: Calc_S_nD, Calc_M_nD,Calc_F_nD,all_deriv_nD
   PUBLIC :: Test_All_Derivatives,Calc_S_MClambda_MLambdaLambda
   PUBLIC :: Test_Matrices,Update_Basis_From_Lambda,Update_Psi_From_DeltaLambda


CONTAINS

SUBROUTINE Calc_S_nD(nb, S)
  !-------------------------------------------------------------
  ! Purpose:
  !   Construct the nD overlap matrix S for the Hagedorn basis.
  !   For an orthonormal basis, this is simply an identity matrix
  !   of size equal to the number of coefficients in psi%CVec.
  !-------------------------------------------------------------
  USE QDUtil_m
  IMPLICIT NONE

  ! Inputs
  INTEGER, INTENT(IN) :: nb

  ! Output
  COMPLEX(KIND=Rkind), ALLOCATABLE, INTENT(OUT) :: S(:,:)

  ! Local variables
  INTEGER :: i

  !-------------------------------------------------------------
  ! Allocate S as identity matrix
  !-------------------------------------------------------------
  ALLOCATE(S(nb,nb))
  S = CZERO
  DO i = 1, nb
     S(i,i) = CONE
  END DO

END SUBROUTINE Calc_S_nD

SUBROUTINE Test_All_Derivatives(psi0)
  USE QDUtil_m
  USE Basis_m
  IMPLICIT NONE

  TYPE(psi_t), INTENT(IN), TARGET :: psi0

  TYPE(psi_t), ALLOCATABLE :: dAllPsi(:)
  INTEGER :: ndim, nDeriv, inb, nio, idx
  INTEGER :: i, expected_size
  CHARACTER(LEN=32) :: label

  !----------------------------------------------------------------
  ! Initialiser le test : dimension, nombre de dérivées
  !----------------------------------------------------------------
  ndim   = SIZE(psi0%Basis%tab_basis) - 1
  nDeriv = 4 * ndim

  !----------------------------------------------------------------
  ! Calcul de toutes les dérivées
  !----------------------------------------------------------------
  CALL all_deriv_nD(psi0, dAllPsi)

  expected_size = nDeriv
  IF (SIZE(dAllPsi) /= expected_size) THEN
    WRITE(*,*) "Erreur Test_All_Derivatives : taille dAllPsi = ", SIZE(dAllPsi), &
               " mais attendu = ", expected_size
    STOP "ERROR in Test_All_Derivatives"
  END IF

  !----------------------------------------------------------------
  ! Affichage
  !----------------------------------------------------------------
  DO i = 1, nDeriv
    ! Calculer inb et nio correspondant à cet indice i
    WRITE(out_unit, *) "dAllPsi(", i, ")%CVec"
    CALL Write_VecMat(dAllPsi(i)%CVec, out_unit,10,info="dAllPsi%CVec")
  END DO

  !----------------------------------------------------------------
  ! Nettoyage
  !----------------------------------------------------------------
  DO i = 1, nDeriv
    CALL dealloc_psi(dAllPsi(i))
  END DO

  DEALLOCATE(dAllPsi)

END SUBROUTINE Test_All_Derivatives



SUBROUTINE deriv_k_nD(psi0, dppsi, ib, nio)
  USE QDUtil_m
  IMPLICIT NONE
  TYPE(psi_t), INTENT(IN),    TARGET :: psi0
  TYPE(psi_t), INTENT(INOUT), TARGET :: dppsi
  INTEGER, INTENT(IN) :: ib, nio
  LOGICAL, PARAMETER :: debug = .FALSE.

  COMPLEX(KIND=Rkind), POINTER :: temp(:,:,:), dBB(:,:,:)
  COMPLEX(KIND=Rkind), POINTER :: d1bb(:,:)
  INTEGER :: Ib1, Ib2, Ib3
  INTEGER :: Inb, Ndim, i1, i3
  INTEGER :: total_size

  Ndim = SIZE(psi0%Basis%tab_basis) - 1
  total_size = SIZE(psi0%CVec)

     ! --- calcul des dimensions locales ---
     IF (ib == 1) THEN
        Ib1 = 1
        Ib2 = psi0%Basis%tab_basis(1)%nb
        Ib3 = PRODUCT(psi0%Basis%tab_basis(2:Ndim)%nb) * psi0%Basis%tab_basis(Ndim+1)%nb

     ELSE IF (ib == Ndim) THEN
        Ib1 = PRODUCT(psi0%Basis%tab_basis(1:Ndim-1)%nb)
        Ib2 = psi0%Basis%tab_basis(Ndim)%nb
        Ib3 = psi0%Basis%tab_basis(Ndim+1)%nb

     ELSE
        Ib1 = PRODUCT(psi0%Basis%tab_basis(1:ib-1)%nb)
        Ib2 = psi0%Basis%tab_basis(ib)%nb
        Ib3 = PRODUCT(psi0%Basis%tab_basis(ib+1:Ndim)%nb) * psi0%Basis%tab_basis(Ndim+1)%nb
     END IF

     IF (Ib1*Ib2*Ib3 /= total_size) THEN
        WRITE(out_unit, *) 'ERROR: shape mismatch in deriv_k_nD: ', Ib1*Ib2*Ib3, ' /= ', total_size
        STOP 'ERROR deriv_k_nD shape mismatch'
     END IF

     ! --- associer la mémoire linéaire CVec à des vues 3D (bounds/rank remapping) ---
     ! NOTE: requires psi0%CVec and dppsi%CVec to have TARGET attribute
     temp(1:Ib1,1:Ib2,1:Ib3) => psi0%CVec
     dBB (1:Ib1,1:Ib2,1:Ib3) => dppsi%CVec

     ! --- choisir la matrice de dérivation correspondant à 'nio' ---
     SELECT CASE (nio)
     CASE (2)
        d1bb(1:Ib2,1:Ib2) => psi0%Basis%tab_basis(ib)%d1bb(:,:,2)
     CASE (3)
        d1bb(1:Ib2,1:Ib2) => psi0%Basis%tab_basis(ib)%d1bb(:,:,3)
     CASE (4)
        d1bb(1:Ib2,1:Ib2) => psi0%Basis%tab_basis(ib)%d1bb(:,:,4)
     CASE (5)
        d1bb(1:Ib2,1:Ib2) => psi0%Basis%tab_basis(ib)%d1bb(:,:,5)
     CASE DEFAULT
        WRITE(out_unit, *) 'ERROR: Invalid nio value'
        STOP 'ERROR in deriv_k_nD'
     END SELECT

     ! --- application de la dérivation (produit matrice x vecteur sur chaque "slice") ---
     DO i3 = 1, Ib3
        DO i1 = 1, Ib1
           dBB(i1,:,i3) = dBB(i1,:,i3) + MATMUL(d1bb, temp(i1,:,i3))
        END DO
     END DO

  IF (debug) THEN
     CALL Write_VecMat(dppsi%CVec, out_unit, 5, info='dppsi%CVec')
     FLUSH(out_unit)
  END IF

END SUBROUTINE deriv_k_nD

SUBROUTINE all_deriv_nD(psi0, dAllPsi)
  USE QDUtil_m
  USE Basis_m
  TYPE(psi_t), INTENT(IN), TARGET       :: psi0
  TYPE(psi_t), ALLOCATABLE, INTENT(OUT) :: dAllPsi(:)
  INTEGER :: ndim, inb, nio, nDeriv, idx

  ndim   = SIZE(psi0%Basis%tab_basis) - 1
  nDeriv = 4 * ndim

  ALLOCATE(dAllPsi(nDeriv))

  idx = 0
  DO inb = 1, ndim
     DO nio = 2, 5
        idx = idx + 1
        CALL init_psi(dAllPsi(idx), psi0%Basis, cplx=.TRUE., grid=.FALSE.)
        CALL deriv_k_nD(psi0, dAllPsi(idx), inb, nio)
     END DO
  END DO

  IF (idx /= nDeriv) THEN
     WRITE(out_unit, *) 'ERROR: unexpected idx count in all_deriv_nD: ', idx, ' expected ', nDeriv
     STOP 'ERROR all_deriv_nD'
  END IF
END SUBROUTINE all_deriv_nD

SUBROUTINE Calc_S_MClambda_MLambdaLambda(psi, S, M_Clambda, M_lambdaLambda)
  !===============================================================
  ! SUBROUTINE Calc_S_MClambda_MLambdaLambda
  !
  ! DESCRIPTION:
  !   Build simultaneously:
  !     - S:      the overlap (recouvrement) matrix between basis functions
  !     - M_Cλ:   coefficient-parameter coupling matrix
  !     - M_λλ:   parameter-parameter Gram matrix of derivative vectors
  !
  ! DIMENSIONS:
  !   Let nC = size of psi%CVec (total number of basis components)
  !   Let nDeriv = 4 * ndim  (ndim = number of DOF)
  !
  ! INPUT:
  !   psi   : wavefunction, includes psi%CVec, psi%Basis, etc.
  !
  ! OUTPUT:
  !   S            : complex(nC, nC) overlap matrix
  !   M_Clambda    : complex(nC, nDeriv) coupling matrix
  !   M_lambdaLambda : complex(nDeriv, nDeriv) Gram matrix of derivatives
  !
  ! REQUIREMENTS:
  !   all_deriv_nD and Calc_S_nD are available
  !   psi%CVec is allocated
  !
  !===============================================================
  USE QDUtil_m
  IMPLICIT NONE

  TYPE(psi_t), INTENT(IN) :: psi
  COMPLEX(KIND=Rkind), ALLOCATABLE, INTENT(OUT) :: S(:,:), M_Clambda(:,:), M_lambdaLambda(:,:)

  ! Local variables
  TYPE(psi_t), ALLOCATABLE :: dAllPsi(:)
  COMPLEX(KIND=Rkind), ALLOCATABLE :: V(:,:)  ! to hold derivative vectors
  INTEGER :: ndim, nDeriv, nC
  INTEGER :: i, j
  COMPLEX(KIND=Rkind) :: conj_vi_dot_S_vj  ! temporary for M_lambdaLambda
  LOGICAL, PARAMETER :: debug = .TRUE.
  !===============================================================
  ! Input checks and setup
  !===============================================================
  ! Check psi%CVec
  IF (.NOT. ALLOCATED(psi%CVec)) THEN
    WRITE(out_unit,*) 'ERROR: Calc_S_MClambda_MLambdaLambda: psi%CVec not allocated'
    STOP 'Calc_S_MClambda_MLambdaLambda'
  END IF

  ndim = SIZE(psi%Basis%tab_basis) - 1
  IF (ndim < 1) THEN
    WRITE(out_unit,*) 'ERROR: Calc_S_MClambda_MLambdaLambda: ndim < 1'
    STOP 'Calc_S_MClambda_MLambdaLambda'
  END IF

  nDeriv = 4 * ndim
  nC = SIZE(psi%CVec)

  IF (debug) THEN
    WRITE(out_unit,*) 'Calc_S_MClambda_MLambdaLambda: ndim=', ndim, ' nDeriv=', nDeriv, ' nC=', nC
  END IF

  !===============================================================
  ! Compute all derivative vectors v_alpha
  !===============================================================
  CALL all_deriv_nD(psi, dAllPsi)
  IF (SIZE(dAllPsi) /= nDeriv) THEN
    WRITE(out_unit,*) 'ERROR: all_deriv_nD returned ', SIZE(dAllPsi), ' items, expected ', nDeriv
    STOP 'Calc_S_MClambda_MLambdaLambda'
  END IF

  !===============================================================
  ! Allocate V(nC, nDeriv), fill its columns by v_alpha
  !===============================================================
  ALLOCATE(V(nC, nDeriv))
  V = CZERO
  DO i = 1, nDeriv
    IF (.NOT. ALLOCATED(dAllPsi(i)%CVec)) THEN
      WRITE(out_unit,*) 'ERROR: derivative vector ', i, ' not allocated'
      STOP 'Calc_S_MClambda_MLambdaLambda'
    END IF
    IF (SIZE(dAllPsi(i)%CVec) /= nC) THEN
      WRITE(out_unit,*) 'ERROR: derivative vector ', i, ' has wrong size: ', SIZE(dAllPsi(i)%CVec), ' vs ', nC
      STOP 'Calc_S_MClambda_MLambdaLambda'
    END IF
    V(:,i) = dAllPsi(i)%CVec
  END DO

  !===============================================================
  ! Build S
  !===============================================================
  CALL Calc_S_nD(nC, S)
  IF (.NOT. ALLOCATED(S)) THEN
    WRITE(out_unit,*) 'ERROR: Calc_S_MClambda_MLambdaLambda: S not allocated by Calc_S_nD'
    STOP 'Calc_S_MClambda_MLambdaLambda'
  END IF
  IF (SIZE(S,1) /= nC .OR. SIZE(S,2) /= nC) THEN
    WRITE(out_unit,*) 'ERROR: S size mismatch'
    STOP 'Calc_S_MClambda_MLambdaLambda'
  END IF

  !===============================================================
  ! Allocate output matrices
  !===============================================================
  IF (ALLOCATED(M_Clambda)) DEALLOCATE(M_Clambda)
  ALLOCATE(M_Clambda(nC, nDeriv))

  IF (ALLOCATED(M_lambdaLambda)) DEALLOCATE(M_lambdaLambda)
  ALLOCATE(M_lambdaLambda(nDeriv, nDeriv))

  !===============================================================
  ! Compute M_Cλ = S * V
  !===============================================================
  M_Clambda = MATMUL(S, V)

  !===============================================================
  ! Compute M_λλ (Gram matrix of derivative vectors, using S)
  !===============================================================
  DO i = 1, nDeriv
    DO j = 1, nDeriv
      ! (M_lambdaLambda)(i,j) = conj( V(:,i) )^T * ( S * V(:,j) )
      conj_vi_dot_S_vj = DOT_PRODUCT( CONJG( V(:,i) ), MATMUL( S, V(:,j) ) )
      !M_lambdaLambda(i,j) = DOT_PRODUCT(CONJG(V(:,i)), V(:,j))!conj_vi_dot_S_vj
      M_lambdaLambda(i,j) = conj_vi_dot_S_vj
    END DO
  END DO

  !===============================================================
  ! Debug printing
  !===============================================================
  IF (debug) THEN
    out_unit = 6
    WRITE(out_unit,*) '==== S matrix ===='
    CALL Write_VecMat(S, out_unit, 5, info='S')
    WRITE(out_unit,*) '==== M_Cλ matrix ===='
    CALL Write_VecMat(M_Clambda, out_unit, 5, info='M_Clambda')
    WRITE(out_unit,*) '==== M_λλ matrix ===='
    CALL Write_VecMat(M_lambdaLambda, out_unit, 5, info='M_lambdaLambda')
    FLUSH(out_unit)
  END IF

  !===============================================================
  ! Cleanup
  !===============================================================
  DO i = 1, nDeriv
    CALL dealloc_psi(dAllPsi(i))
  END DO
  DEALLOCATE(dAllPsi)
  IF (ALLOCATED(V)) DEALLOCATE(V)

END SUBROUTINE Calc_S_MClambda_MLambdaLambda


SUBROUTINE Calc_F_nD(psi, H, F_nD,FC_nD,F_L)
  USE QDUtil_m
  IMPLICIT NONE

  TYPE(psi_t), INTENT(IN) :: psi
  TYPE(Op_t),   INTENT(IN) :: H
  COMPLEX(KIND=Rkind), ALLOCATABLE,OPTIONAL, INTENT(OUT) :: F_nD(:)
  COMPLEX(KIND=Rkind), ALLOCATABLE,OPTIONAL, INTENT(OUT) :: FC_nD(:)
  COMPLEX(KIND=Rkind), ALLOCATABLE,OPTIONAL, INTENT(OUT) :: F_L(:)
  TYPE(psi_t), ALLOCATABLE :: dAllPsi(:)
  TYPE(psi_t) :: Hpsi, tmp
  INTEGER :: ndim, nDeriv, nC, i, Ntot
  LOGICAL, PARAMETER :: debug = .TRUE.

  !--------------------------------------------------------------------
  ! Pré-conditions
  !--------------------------------------------------------------------
  IF (.NOT. ALLOCATED(psi%CVec)) THEN
    WRITE(out_unit,*) 'ERROR Calc_F_nD: psi%CVec not allocated'
    STOP 'ERROR Calc_F_nD'
  END IF

  IF (.NOT. ASSOCIATED(psi%Basis)) THEN
    WRITE(out_unit,*) 'ERROR Calc_F_nD: psi%Basis not associated'
    STOP 'ERROR Calc_F_nD'
  END IF

  ndim = SIZE(psi%Basis%tab_basis) - 1
  IF (ndim < 1) THEN
    WRITE(out_unit,*) 'ERROR Calc_F_nD: ndim < 1'
    STOP 'ERROR Calc_F_nD'
  END IF

  nDeriv = 4 * ndim
  nC = SIZE(psi%CVec)
  Ntot = nC + nDeriv

  !--------------------------------------------------------------------
  ! Compute Hpsi = H |Psi>
  !--------------------------------------------------------------------
  CALL init_psi(Hpsi, psi%Basis, cplx = .TRUE., grid = .FALSE.)
  CALL init_psi(tmp,  psi%Basis, cplx = .TRUE., grid = .FALSE.)

  ! Vérifier que init_psi a bien alloué CVec et qu'elle a la bonne taille
  IF (.NOT. ALLOCATED(Hpsi%CVec)) THEN
    WRITE(out_unit,*) 'ERROR Calc_F_nD: init_psi did not allocate Hpsi%CVec'
    STOP 'ERROR Calc_F_nD'
  END IF
  IF (SIZE(Hpsi%CVec) /= nC) THEN
    WRITE(out_unit,*) 'ERROR Calc_F_nD: Hpsi%CVec size mismatch (', SIZE(Hpsi%CVec), ' vs ', nC, ')'
    STOP 'ERROR Calc_F_nD'
  END IF

  CALL calc_Oppsi(H, psi, Hpsi)

  !--------------------------------------------------------------------
  ! Compute derivative vectors v_alpha
  !--------------------------------------------------------------------
  CALL all_deriv_nD(psi, dAllPsi)
  IF (.NOT. ALLOCATED(dAllPsi)) THEN
    WRITE(out_unit,*) 'ERROR Calc_F_nD: all_deriv_nD did not allocate dAllPsi'
    STOP 'ERROR Calc_F_nD'
  END IF

  IF (SIZE(dAllPsi) /= nDeriv) THEN
    WRITE(out_unit,*) 'ERROR Calc_F_nD: all_deriv_nD returned ', SIZE(dAllPsi), ' items; expected ', nDeriv
    STOP 'ERROR Calc_F_nD'
  END IF

  DO i = 1, nDeriv
    IF (.NOT. ALLOCATED(dAllPsi(i)%CVec)) THEN
      WRITE(out_unit,*) 'ERROR Calc_F_nD: derivative ', i, ' not allocated'
      STOP 'ERROR Calc_F_nD'
    END IF
    IF ( SIZE(dAllPsi(i)%CVec) /= nC ) THEN
      WRITE(out_unit,*) 'ERROR Calc_F_nD: derivative ', i, ' size mismatch vs psi%CVec'
      STOP 'ERROR Calc_F_nD'
    END IF
  END DO

  !--------------------------------------------------------------------
  ! Allocate et initialise F_nD OR FC_nD OR F_lambda
  !--------------------------------------------------------------------
  
  IF (PRESENT(F_nD)) THEN
    IF (ALLOCATED(F_nD)) THEN
      DEALLOCATE(F_nD)
    END IF
    ALLOCATE(F_nD(Ntot))
    F_nD = CZERO
  END IF

  IF (PRESENT(FC_nD)) THEN
    IF (ALLOCATED(FC_nD)) THEN
      DEALLOCATE(FC_nD)
    END IF
    ALLOCATE(FC_nD(nC))
    FC_nD = CZERO
  END IF

  IF (PRESENT(F_L)) THEN
    IF (ALLOCATED(F_L)) THEN
      DEALLOCATE(F_L)
    END IF
    ALLOCATE(F_L(nDeriv))
    F_L = CZERO
  END IF

  !--------------------------------------------------------------------
  ! Compute F_C = < phi_m | Hpsi >
  !--------------------------------------------------------------------
  IF (.NOT. ALLOCATED(tmp%CVec)) THEN
    WRITE(out_unit,*) 'ERROR Calc_F_nD: tmp%CVec not allocated by init_psi'
    STOP 'ERROR Calc_F_nD'
  END IF
  IF (SIZE(tmp%CVec) /= nC) THEN
    WRITE(out_unit,*) 'ERROR Calc_F_nD: tmp%CVec size mismatch (', SIZE(tmp%CVec), ' vs ', nC, ')'
    STOP 'ERROR Calc_F_nD'
  END IF

  DO i = 1, nC
    tmp%CVec = CZERO         ! mise à zéro entière
    tmp%CVec(i) = CONE
    IF (PRESENT(F_nD)) THEN
      F_nD(i) = Hpsi%CVec(i)!DOT_PRODUCT(CONJG(tmp%CVec), Hpsi%CVec)
    END IF
    IF (PRESENT(FC_nD)) THEN
      FC_nD(i) = Hpsi%CVec(i)!DOT_PRODUCT(CONJG(tmp%CVec), Hpsi%CVec)
    END IF
  END DO

  !--------------------------------------------------------------------
  ! Compute F_lambda = < v_a | Hpsi >
  !--------------------------------------------------------------------
  DO i = 1, nDeriv
    IF (PRESENT(F_nD)) THEN
      F_nD(nC + i) = DOT_PRODUCT( CONJG(dAllPsi(i)%CVec), Hpsi%CVec )
    END IF
    IF (PRESENT(F_L)) THEN
      F_L(i) = DOT_PRODUCT( CONJG(dAllPsi(i)%CVec), Hpsi%CVec )
    END IF
  END DO

  !--------------------------------------------------------------------
  ! Debug printing
  !--------------------------------------------------------------------
  IF (debug) THEN
    IF (PRESENT(F_nD)) THEN
      WRITE(out_unit,*) 'F_nD vector (size ', SIZE(F_nD), '):'
      CALL Write_VecMat( F_nD, out_unit, 5, info = 'F_nD' )
      FLUSH(out_unit)
    END IF
    IF (PRESENT(FC_nD)) THEN
      WRITE(out_unit,*) 'FC_nD vector (size ', SIZE(FC_nD), '):'
      CALL Write_VecMat( FC_nD, out_unit, 5, info = 'FC_nD' )
      FLUSH(out_unit)
    END IF
    IF (PRESENT(F_L)) THEN
      WRITE(out_unit,*) 'F_lambda vector (size ', SIZE(F_L), '):'
      CALL Write_VecMat( F_L, out_unit, 5, info = 'F_L' )
      FLUSH(out_unit)
    END IF
  END IF

  !--------------------------------------------------------------------
  ! Normal exit: cleanup et return
  !--------------------------------------------------------------------
  ! Cleanup sûr : on déalloue seulement si alloué
  IF (ALLOCATED(Hpsi%CVec)) CALL dealloc_psi(Hpsi)
  IF (ALLOCATED(tmp%CVec))  CALL dealloc_psi(tmp)

  IF (ALLOCATED(dAllPsi)) THEN
    DO i = 1, SIZE(dAllPsi)
      IF (ALLOCATED(dAllPsi(i)%CVec)) CALL dealloc_psi(dAllPsi(i))
    END DO
    DEALLOCATE(dAllPsi)
  END IF

END SUBROUTINE Calc_F_nD


 !=======================================================================
! SUBROUTINE: Update_Basis_From_Lambda
!
! PURPOSE:
!   Update the parameters of a multidimensional basis according to
!   the increments DeltaLambda_L over a timestep dt.
!
! INPUT:
!   Basis         : the multidimensional basis to update
!   DeltaLambda_L : vector of increments (q0, p, Re(alpha), Im(alpha), ...)
!   dt            : timestep
!
! OUTPUT:
!   Basis         : updated basis with new parameters
!
! NOTES:
!   - Assumes DeltaLambda_L is ordered as:
!       (q01,p1,Re(a1),Im(a1), ..., q0k,pk,Re(ak),Im(ak), ...)
!=======================================================================
SUBROUTINE Update_Basis_From_Lambda(Basis, DeltaLambda_L, dt)
  USE QDUtil_m
  IMPLICIT NONE

  TYPE(Basis_t), INTENT(INOUT)       :: Basis
  COMPLEX(KIND=Rkind), INTENT(IN)    :: DeltaLambda_L(:)
  REAL(KIND=Rkind), INTENT(IN)       :: dt

  INTEGER :: ndim, k, idx
  REAL(KIND=Rkind) :: dq, dp
  COMPLEX(KIND=Rkind) :: dAlpha
  LOGICAL ,PARAMETER   :: debug=.TRUE.

  ! Number of dimensions
  ndim = SIZE(Basis%tab_basis)-1

  IF(debug) THEN 
    WRITE(out_unit,*) 'DeltaLambda_L :'
    CALL Write_VecMat(DeltaLambda_L, out_unit, 6, info= 'DeltaLambda_L')
  END IF

  IF (SIZE(DeltaLambda_L) /= 4*ndim) THEN
     WRITE(*,*) "ERROR: DeltaLambda_L size does not match number of dimensions"
     STOP
  END IF

  ! Loop over dimensions
  DO k = 1, ndim
     idx = 4*(k-1) + 1

     ! Increment for Q0 and P (real)
     dq = REAL(DeltaLambda_L(idx))*(-REAL(dt))    ! factor -i*dt mapped to -dt here
     dp = REAL(DeltaLambda_L(idx+1))*(-REAL(dt))

     Basis%tab_basis(k)%Q0 = Basis%tab_basis(k)%Q0 + dq
     Basis%tab_basis(k)%imp_k = Basis%tab_basis(k)%imp_k + dp

     ! Increment for alpha (complex)
     dAlpha = DeltaLambda_L(idx+2) + EYE*DeltaLambda_L(idx+3)
     Basis%tab_basis(k)%Alpha = Basis%tab_basis(k)%Alpha + (-EYE*dt)*dAlpha
  END DO

  ! Reconstruct the primitive basis (grid points, weights, derivatives)
  CALL construct_primitive_basis_temp(Basis)

END SUBROUTINE Update_Basis_From_Lambda


!=======================================================================
! SUBROUTINE: Update_Psi_From_DeltaLambda
!
! PURPOSE:
!   Update the propagated wavefunction psi_dt from psi, applying both:
!     - coefficient increments ΔΛ_C
!     - variational parameter increments ΔΛ_λ to the basis
!
! INPUT:
!   psi          : reference wavefunction (type psi_t)
!   DeltaLambda_C: increments for the coefficients C (size nC)
!   DeltaLambda_L: increments for the λ-parameters (size nDeriv)
!   dt           : timestep
!
! OUTPUT:
!   psi_dt       : updated wavefunction with new coefficients and basis
!
! REQUIREMENTS:
!   - all_deriv_nD must generate the derivative vectors v_α
!   - Update_Basis_From_Lambda must update the basis according to ΔΛ_λ
!=======================================================================
SUBROUTINE Update_Psi_From_DeltaLambda(psi, psi_dt, D_C, D_L, dt)
  USE QDUtil_m
  IMPLICIT NONE

  TYPE(psi_t), INTENT(IN)             :: psi
  TYPE(psi_t), INTENT(INOUT)          :: psi_dt
  COMPLEX(KIND=Rkind), INTENT(IN)     :: D_C(:)
  COMPLEX(KIND=Rkind), INTENT(IN)     :: D_L(:)
  REAL(KIND=Rkind), INTENT(IN)        :: dt

  TYPE(psi_t), ALLOCATABLE            :: dAllPsi(:)
  INTEGER                             :: i, nC, nDeriv

  !------------------------------------------------------------
  ! Step 0. Initialize psi_dt
  !------------------------------------------------------------
  ! init_psi(psi_dt, psi%Basis, cplx=.TRUE., grid=.FALSE.)
  psi_dt%CVec = psi%CVec

  nC     = SIZE(psi%CVec)
  nDeriv = SIZE(D_L)

  ! Safety check
  IF (SIZE(D_C) /= nC) THEN
     WRITE(*,*) "ERROR Update_Psi_From_DeltaLambda: wrong size for DeltaLambda_C"
     STOP "ERROR Update_Psi_From_DeltaLambda: wrong size for DeltaLambda_C"
  END IF

  !------------------------------------------------------------
  ! Step 1. Update coefficients C
  !------------------------------------------------------------
  psi_dt%CVec = psi_dt%CVec +(-EYE * dt) * D_C

  !------------------------------------------------------------
  ! Step 3. Apply the contribution of ΔΛ_λ to ψ via derivative vectors
  !         Δψ_λ = Σ ΔΛ_λ(α) v_α
  !------------------------------------------------------------
  CALL all_deriv_nD(psi, dAllPsi)
  IF (SIZE(dAllPsi) /= nDeriv) THEN
     WRITE(*,*) "ERROR Update_Psi_From_DeltaLambda: mismatch nDeriv"
     STOP
  END IF

  DO i = 1, nDeriv
     IF (.NOT. ALLOCATED(dAllPsi(i)%CVec)) THEN
        WRITE(*,*) "ERROR: derivative ", i, " not allocated"
        STOP
     END IF
     psi_dt%CVec = psi_dt%CVec + (-EYE * dt) * D_L(i) * dAllPsi(i)%CVec
  END DO


  DO i =1,nC
    WRITE(out_unit,*) i,psi_dt%CVec(i) 
  END DO

  stop 'ici'
!------------------------------------------------------------
! Step 2. Update the basis using ΔΛ_λ
!         This modifies centers, widths, etc.
!------------------------------------------------------------
CALL Update_Basis_From_Lambda(psi_dt%Basis, D_L, dt)






  !------------------------------------------------------------
  ! Step 4. Cleanup derivative vectors
  !------------------------------------------------------------
  DO i = 1, nDeriv
     CALL dealloc_psi(dAllPsi(i))
  END DO
  DEALLOCATE(dAllPsi)

END SUBROUTINE Update_Psi_From_DeltaLambda


SUBROUTINE Test_Matrices(psi, H)
  !===============================================================
  ! SUBROUTINE Test_Matrices
  !
  ! DESCRIPTION:
  !   Cette sous-routine construit les matrices S, M_Cλ et M_λλ
  !   à partir de l'objet psi et les affiche à l'écran.
  !
  ! INPUT:
  !   psi   : objet de type psi_t contenant les informations nécessaires.
  !
  !===============================================================
  USE QDUtil_m
  IMPLICIT NONE

  TYPE(psi_t), INTENT(IN) :: psi
  COMPLEX(KIND=Rkind), ALLOCATABLE :: S(:,:), M_Clambda(:,:), M_lambdaLambda(:,:), M_nD(:,:), F_nD(:)
  TYPE(Op_t) , INTENT(IN) :: H
  INTEGER :: nC, nDeriv

  !===============================================================
  ! Déterminer la taille des matrices
  !===============================================================
  nC = SIZE(psi%CVec)
  nDeriv = 4 * (SIZE(psi%Basis%tab_basis) - 1)

  !===============================================================
  ! Construire les matrices
  !===============================================================
  CALL Calc_S_MClambda_MLambdaLambda(psi, S, M_Clambda, M_lambdaLambda)
  CALL Calc_M_nD(psi, M_nD)
  CALL Calc_F_nD(psi, H, F_nD)

  !===============================================================
  ! Afficher les matrices
  !===============================================================
  WRITE(out_unit,*) 'Matrice de recouvrement S :'
  CALL Write_VecMat(S, out_unit, 6, info= 'S')

  WRITE(out_unit,*) 'Matrice M_Cλ (coefficient-paramètre) :'
  CALL Write_VecMat(M_Clambda, out_unit, 6, info='M_Cλ')

  WRITE(out_unit,*) 'Matrice M_λλ (paramètre-paramètre) :'
  CALL Write_VecMat(M_lambdaLambda, out_unit, 6, info='M_λλ')

  WRITE(out_unit,*) 'Matrice M_nD (variational matrix) :'
  CALL Write_VecMat(M_nD, out_unit, 6, info='M_nD')

  WRITE(out_unit,*) 'Matrice F_nD (force-like vector) :'
  CALL Write_VecMat(F_nD, out_unit, 6, info='F_nD')

END SUBROUTINE Test_Matrices

SUBROUTINE Calc_M_nD(psi, M_nD)
  !====================================================================
  ! SUBROUTINE Calc_M_nD
  !
  ! DESCRIPTION:
  !   Constructs the full variational matrix M for psi, using
  !   Calc_S_MClambda_MLambdaLambda to build the blocks:
  !
  !     M = | S         M_Cλ     |
  !         | M_Cλ†     M_λλ     |
  !
  ! INPUT:
  !   psi   : wavefunction object (with psi%CVec, psi%Basis etc.)
  !
  ! OUTPUT:
  !   M_nD  : the full complex matrix of size (nC + nDeriv) × (nC + nDeriv)
  !
  !====================================================================

  USE QDUtil_m
  IMPLICIT NONE

  TYPE(psi_t), INTENT(IN) :: psi
  COMPLEX(KIND=Rkind), ALLOCATABLE, INTENT(OUT) :: M_nD(:,:)

  ! Blocks
  COMPLEX(KIND=Rkind), ALLOCATABLE :: S(:,:), M_Cl(:,:), M_Ll(:,:)
  ! Local sizes
  INTEGER :: ndim, nDeriv, nC, Ntot, Nmat
  LOGICAL, PARAMETER :: debug = .TRUE.

  !====================================================================
  ! Input and dimension checks
  !====================================================================
  IF (.NOT. ALLOCATED(psi%CVec)) THEN
    WRITE(*,*) 'ERROR Calc_M_nD: psi%CVec not allocated'
    STOP 'Calc_M_nD'
  END IF

  IF (.NOT. ASSOCIATED(psi%Basis)) THEN
    WRITE(*,*) 'ERROR Calc_M_nD: psi%Basis not associated'
    STOP 'Calc_M_nD'
  END IF

  ndim = SIZE(psi%Basis%tab_basis) - 1
  IF (ndim < 1) THEN
    WRITE(*,*) 'ERROR Calc_M_nD: ndim < 1'
    STOP 'Calc_M_nD'
  END IF

  nDeriv = 4 * ndim
  nC = SIZE(psi%CVec)

  IF (debug) THEN
    WRITE(out_unit,*) 'Calc_M_nD: ndim=', ndim, ' nDeriv=', nDeriv, ' nC=', nC
    FLUSH(out_unit)
  END IF

  !====================================================================
  ! Compute S, M_Cλ, M_λλ via composite routine
  !====================================================================
  CALL Calc_S_MClambda_MLambdaLambda(psi, S, M_Cl, M_ll)

  ! Check allocations and sizes
  IF (.NOT. ALLOCATED(S)) THEN
    WRITE(*,*) 'ERROR Calc_M_nD: S not allocated by Calc_S_MClambda_MLambdaLambda'
    STOP 'Calc_M_nD'
  END IF
  IF (SIZE(S,1) /= nC .OR. SIZE(S,2) /= nC) THEN
    WRITE(*,*) 'ERROR Calc_M_nD: S size mismatch: got (', SIZE(S,1),',',SIZE(S,2),') exp (',nC,',',nC,')'
    STOP 'Calc_M_nD'
  END IF

  IF (.NOT. ALLOCATED(M_Cl)) THEN
    WRITE(*,*) 'ERROR Calc_M_nD: M_Cl not allocated'
    STOP 'Calc_M_nD'
  END If
  IF (SIZE(M_Cl,1) /= nC) THEN
    WRITE(*,*) 'ERROR Calc_M_nD: M_Cl rows mismatch: got ', SIZE(M_Cl,1), ' expected ', nC
    STOP 'Calc_M_nD'
  End If

  ! Determine total λ-block size
  Ntot = SIZE(M_Cl,2)  ! number of columns = nDeriv

  IF (SIZE(M_ll,1) /= Ntot .OR. SIZE(M_ll,2) /= Ntot) THEN
    WRITE(*,*) 'ERROR Calc_M_nD: M_ll size mismatch: got (', SIZE(M_ll,1),',',SIZE(M_ll,2),') expected (', Ntot,',', Ntot,')'
    STOP 'Calc_M_nD'
  End If

  !====================================================================
  ! Assemble full M_nD
  !====================================================================
  Nmat = nC + Ntot

  IF (ALLOCATED(M_nD)) DEALLOCATE(M_nD)
  ALLOCATE(M_nD(Nmat, Nmat))
  M_nD = CZERO

  ! Upper-left: S
  M_nD(1:nC, 1:nC) = S

  ! Upper-right: M_Cλ
  M_nD(1:nC, nC+1:Nmat) = M_Cl

  ! Lower-left: conjugate transpose
  M_nD(nC+1:Nmat, 1:nC) = CONJG( TRANSPOSE( M_Cl ) )

  ! Lower-right: M_λλ
  M_nD(nC+1:Nmat, nC+1:Nmat) = M_ll

  !====================================================================
  ! Debug print
  !====================================================================
  IF (debug) THEN
    WRITE(out_unit,*) 'Calc_M_nD: Full M_nD size =', Nmat, '×', Nmat
    CALL Write_VecMat(M_nD, out_unit, 5, info='M_nD')
    FLUSH(out_unit)
  END IF

  !====================================================================
  ! Clean up
  !====================================================================
  IF (ALLOCATED(S)) DEALLOCATE(S)
  IF (ALLOCATED(M_Cl)) DEALLOCATE(M_Cl)
  IF (ALLOCATED(M_ll)) DEALLOCATE(M_ll)

END SUBROUTINE Calc_M_nD

SUBROUTINE Psi_Lambda_Transform(psi, Lambda, toLambda)
  USE QDUtil_m
  IMPLICIT NONE
  TYPE(psi_t), INTENT(INOUT)  :: psi
  COMPLEX(KIND=Rkind), ALLOCATABLE, INTENT(INOUT) :: Lambda(:)
  LOGICAL, INTENT(IN)         :: toLambda
  LOGICAL, PARAMETER :: debug = .FALSE.

  INTEGER :: nb, ndof, i, pos

  ! Nombre de coefficients et DOFs
  nb = SIZE(psi%CVec)
  ndof = SIZE(psi%Basis%tab_basis) - 1  


  IF (toLambda) THEN
     ! =====================================================
     ! psi -> Lambda
     ! =====================================================
     IF (.NOT. ALLOCATED(Lambda)) ALLOCATE(Lambda(nb + 4*ndof))
     IF(debug) THEN
        WRITE(out_unit,*) 'psi->Lambda'
        WRITE(out_unit,*) 'size(psi%CVec) = ', SIZE(psi%CVec)
        WRITE(out_unit,*) 'size(Lambda) = ', SIZE(Lambda)
     END IF
     

     ! 1. Copier les coefficients C
     Lambda(1:nb) = psi%CVec(:)

     ! 2. Ajouter les paramètres (q0, p, a, b) pour chaque DOF
     pos = nb
     DO i = 1, ndof
        pos = pos + 1
        Lambda(pos) = CMPLX(psi%Basis%tab_basis(i)%Q0, ZERO, KIND=Rkind)     ! q0
        pos = pos + 1
        Lambda(pos) = CMPLX(psi%Basis%tab_basis(i)%Imp_k, ZERO, KIND=Rkind)  ! p
        pos = pos + 1
        Lambda(pos) = CMPLX(REAL(psi%Basis%tab_basis(i)%Alpha), ZERO, KIND=Rkind)  ! a
        pos = pos + 1
        Lambda(pos) = CMPLX(AIMAG(psi%Basis%tab_basis(i)%Alpha), ZERO, KIND=Rkind) ! b
     END DO

  ELSE
     ! =====================================================
     ! Lambda -> psi
     ! =====================================================
     IF (.NOT. ALLOCATED(psi%CVec)) THEN
        WRITE(out_unit,*) 'ERROR Psi_Lambda_Transform: psi%CVec is not allocated'
        STOP 'Psi_Lambda_Transform : psi%CVec is not allocated'
     END IF

     ! 1. Extraire les coefficients
     IF(debug) THEN
        WRITE(out_unit,*) 'Lambda->psi'
        WRITE(out_unit,*) 'size(psi%CVec) = ', SIZE(psi%CVec)
        WRITE(out_unit,*) 'size(Lambda) = ', SIZE(Lambda)
     END IF

     IF(debug) THEN
     DO i = 1, nb
        WRITE(out_unit, *) i,Lambda(i)
     END DO
     END IF
     
     psi%CVec(1:nb) = Lambda(1:nb)
     DO i = 1, nb
        WRITE(out_unit, *) i,psi%CVec(i)
     END DO
        
     ! 2. Extraire les paramètres (q0, p, a, b)
     pos = nb
     DO i = 1, ndof
        pos = pos + 1
        psi%Basis%tab_basis(i)%Q0    = REAL(Lambda(pos))
        pos = pos + 1
        psi%Basis%tab_basis(i)%Imp_k = REAL(Lambda(pos))
        pos = pos + 1
        psi%Basis%tab_basis(i)%Alpha = CMPLX(REAL(Lambda(pos)), AIMAG(psi%Basis%tab_basis(i)%Alpha), KIND=Rkind)
        pos = pos + 1
        psi%Basis%tab_basis(i)%Alpha = CMPLX(REAL(psi%Basis%tab_basis(i)%Alpha), REAL(Lambda(pos)), KIND=Rkind)
     END DO

     CALL construct_primitive_basis_temp(psi%Basis)

  END IF

END SUBROUTINE Psi_Lambda_Transform


SUBROUTINE Test_Psi_Lambda_Transform(psi)
  USE QDUtil_m
  IMPLICIT NONE
  TYPE(psi_t), INTENT(INOUT) :: psi
  TYPE(psi_t)             :: psi2
  COMPLEX(KIND=Rkind), ALLOCATABLE :: Lambda(:)
  INTEGER :: i

  WRITE(out_unit, *) ">>> Test_Psi_Lambda_Transform"

  ! ------------------------------
  ! psi -> Lambda
  ! ------------------------------
  CALL Psi_Lambda_Transform(psi, Lambda, .TRUE.)

  WRITE(out_unit, *) "Lambda vector (from psi):"
  DO i = 1, SIZE(Lambda)
     WRITE(out_unit, '(I3,2X,F12.6,2X,F12.6)') i, REAL(Lambda(i)), AIMAG(Lambda(i))
  END DO
  ! ------------------------------
  ! Lambda -> psi2
  ! ------------------------------
   CALL init_psi(psi2, psi%Basis, .TRUE., .FALSE.)
  CALL Psi_Lambda_Transform(psi2, Lambda, .FALSE.)

  WRITE(out_unit, *) "psi%CVec (original):"
  DO i = 1, SIZE(psi%CVec)
     WRITE(out_unit, '(I3,2X,F12.6,2X,F12.6)') i, REAL(psi%CVec(i)), AIMAG(psi%CVec(i))
  END DO

  WRITE(out_unit, *) "psi2%CVec (reconstructed):"
  DO i = 1, SIZE(psi2%CVec)
     WRITE(out_unit, '(I3,2X,F12.6,2X,F12.6)') i, REAL(psi2%CVec(i)), AIMAG(psi2%CVec(i))
  END DO

  WRITE(out_unit, *) "<<< End Test_Psi_Lambda_Transform"

  !Clear memory
  CALL dealloc_psi(psi2)

END SUBROUTINE Test_Psi_Lambda_Transform

END MODULE 
