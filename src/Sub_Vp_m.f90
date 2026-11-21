MODULE Sub_Vp_m
   USE QDUtil_m
   USE NDindex_m
   USE Op_m
   USE SmolyakSparseBasisAlgebra_m, ONLY: DotSmolyakBasis
   USE Basis_m
   USE psi_m
   USE polyortho_m, ONLY : hercom, poly_Hermite
   IMPLICIT NONE

   PRIVATE
   PUBLIC :: Psi_Lambda_Transform, Test_Psi_Lambda_Transform
   PUBLIC :: Calc_S_nD, Calc_M_nD, Calc_F_nD, all_deriv_nD
   PUBLIC :: Test_All_Derivatives, Calc_S_MClambda_MLambdaLambda
   PUBLIC :: Test_Matrices, Update_Basis_From_Lambda, Update_Psi_From_DeltaLambda

CONTAINS

!=======================================================================
! Build the coefficient block of the TDVP metric.
! For the current code path, the basis is assumed orthonormal at fixed time,
! therefore this block is the identity matrix.
!=======================================================================
SUBROUTINE Calc_S_nD(nb, S)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nb
  COMPLEX(KIND=Rkind), ALLOCATABLE, INTENT(OUT) :: S(:,:)
  INTEGER :: i

  ALLOCATE(S(nb,nb))
  S = CZERO
  DO i = 1, nb
     S(i,i) = CONE
  END DO

END SUBROUTINE Calc_S_nD

!=======================================================================
! Diagnostic routine: print all tangent vectors d|psi>/d lambda_alpha.
!=======================================================================
SUBROUTINE Test_All_Derivatives(psi0)
  IMPLICIT NONE

  TYPE(psi_t), INTENT(IN), TARGET :: psi0
  TYPE(psi_t), ALLOCATABLE :: dAllPsi(:)
  INTEGER :: ndim, nDeriv, i

  ndim   = SIZE(psi0%Basis%tab_basis) - 1
  nDeriv = 4 * ndim

  CALL all_deriv_nD(psi0, dAllPsi)

  IF (SIZE(dAllPsi) /= nDeriv) THEN
     WRITE(out_unit,*) 'ERROR in Test_All_Derivatives: invalid number of derivatives.'
     STOP 'ERROR in Test_All_Derivatives'
  END IF

  DO i = 1, nDeriv
     WRITE(out_unit, *) 'dAllPsi(', i, ')%CVec'
     CALL Write_VecMat(dAllPsi(i)%CVec, out_unit, 10, info='dAllPsi%CVec')
  END DO

  DO i = 1, nDeriv
     CALL dealloc_psi(dAllPsi(i))
  END DO
  DEALLOCATE(dAllPsi)

END SUBROUTINE Test_All_Derivatives

!=======================================================================
! Apply one parameter-derivative matrix d1bb(:,:,nio) along a selected
! multidimensional basis direction ib.
!
! The output vector is dppsi = (d / d lambda_{ib,nio}) |psi0> expressed in
! the same coefficient basis as psi0.
!=======================================================================
SUBROUTINE deriv_k_nD(psi0, dppsi, ib, nio)
  IMPLICIT NONE

  TYPE(psi_t), INTENT(IN),    TARGET :: psi0
  TYPE(psi_t), INTENT(INOUT), TARGET :: dppsi
  INTEGER, INTENT(IN) :: ib, nio

  COMPLEX(KIND=Rkind), POINTER :: d1bb(:,:)
  COMPLEX(KIND=Rkind), ALLOCATABLE :: vin(:), vout(:)
  INTEGER, ALLOCATABLE :: Ib1(:), Ib2(:), Ib3(:)
  INTEGER :: n1, n2, n3
  INTEGER :: i1, i2, i3, idx

  IF (.NOT. ASSOCIATED(psi0%Basis)) THEN
     WRITE(out_unit,*) 'ERROR deriv_k_nD: psi0%Basis is not associated.'
     STOP 'ERROR deriv_k_nD'
  END IF
  IF (.NOT. ALLOCATED(psi0%CVec)) THEN
     WRITE(out_unit,*) 'ERROR deriv_k_nD: psi0%CVec is not allocated.'
     STOP 'ERROR deriv_k_nD'
  END IF
  IF (.NOT. ALLOCATED(dppsi%CVec)) THEN
     WRITE(out_unit,*) 'ERROR deriv_k_nD: dppsi%CVec is not allocated.'
     STOP 'ERROR deriv_k_nD'
  END IF

  CALL Calc_index(Ib1=Ib1, Ib2=Ib2, Ib3=Ib3, Basis=psi0%Basis)
  n1 = Ib1(ib)
  n2 = Ib2(ib)
  n3 = Ib3(ib)

  SELECT CASE (nio)
  CASE (2)
     d1bb => psi0%Basis%tab_basis(ib)%d1bb(:,:,2)
  CASE (3)
     d1bb => psi0%Basis%tab_basis(ib)%d1bb(:,:,3)
  CASE (4)
     d1bb => psi0%Basis%tab_basis(ib)%d1bb(:,:,4)
  CASE (5)
     d1bb => psi0%Basis%tab_basis(ib)%d1bb(:,:,5)
  CASE DEFAULT
     WRITE(out_unit,*) 'ERROR deriv_k_nD: invalid nio = ', nio
     STOP 'ERROR deriv_k_nD'
  END SELECT

  ALLOCATE(vin(n2), vout(n2))
  dppsi%CVec = CZERO

  DO i3 = 1, n3
     DO i1 = 1, n1
        DO i2 = 1, n2
           idx = i1 + (i2-1)*n1 + (i3-1)*n1*n2
           vin(i2) = psi0%CVec(idx)
        END DO

        vout = MATMUL(d1bb, vin)

        DO i2 = 1, n2
           idx = i1 + (i2-1)*n1 + (i3-1)*n1*n2
           dppsi%CVec(idx) = vout(i2)
        END DO
     END DO
  END DO

  DEALLOCATE(vin, vout, Ib1, Ib2, Ib3)

END SUBROUTINE deriv_k_nD

!=======================================================================
! Compute the full tangent set ordered as:
!   (dq0_1, dp_1, da_1, db_1, ..., dq0_ndim, dp_ndim, da_ndim, db_ndim)
!=======================================================================
SUBROUTINE all_deriv_nD(psi0, dAllPsi)
  IMPLICIT NONE

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

END SUBROUTINE all_deriv_nD

!=======================================================================
! Build the TDVP blocks:
!   S              = <phi_i | phi_j>
!   M_Clambda      = <phi_i | dpsi/dlambda_a>
!   M_lambdaLambda = <dpsi/dlambda_a | dpsi/dlambda_b>
!
! With the present orthonormal basis conventions, S is the identity matrix.
!=======================================================================
SUBROUTINE Calc_S_MClambda_MLambdaLambda(psi, S, M_Clambda, M_lambdaLambda)
  IMPLICIT NONE

  TYPE(psi_t), INTENT(IN) :: psi
  COMPLEX(KIND=Rkind), ALLOCATABLE, INTENT(OUT) :: S(:,:), M_Clambda(:,:), M_lambdaLambda(:,:)

  TYPE(psi_t), ALLOCATABLE :: dAllPsi(:)
  COMPLEX(KIND=Rkind), ALLOCATABLE :: V(:,:)
  INTEGER :: ndim, nDeriv, nC, i

  IF (.NOT. ASSOCIATED(psi%Basis)) THEN
     WRITE(out_unit,*) 'ERROR Calc_S_MClambda_MLambdaLambda: psi%Basis is not associated.'
     STOP 'Calc_S_MClambda_MLambdaLambda'
  END IF
  IF (.NOT. ALLOCATED(psi%CVec)) THEN
     WRITE(out_unit,*) 'ERROR Calc_S_MClambda_MLambdaLambda: psi%CVec is not allocated.'
     STOP 'Calc_S_MClambda_MLambdaLambda'
  END IF

  ndim   = SIZE(psi%Basis%tab_basis) - 1
  nDeriv = 4 * ndim
  nC     = SIZE(psi%CVec)

  CALL all_deriv_nD(psi, dAllPsi)

  ALLOCATE(V(nC, nDeriv))
  DO i = 1, nDeriv
     V(:,i) = dAllPsi(i)%CVec
  END DO

  CALL Calc_S_nD(nC, S)

  ALLOCATE(M_Clambda(nC, nDeriv))
  ALLOCATE(M_lambdaLambda(nDeriv, nDeriv))

  M_Clambda      = MATMUL(S, V)
  M_lambdaLambda = MATMUL(CONJG(TRANSPOSE(V)), M_Clambda)

  DO i = 1, nDeriv
     CALL dealloc_psi(dAllPsi(i))
  END DO
  DEALLOCATE(dAllPsi, V)

END SUBROUTINE Calc_S_MClambda_MLambdaLambda

!=======================================================================
! Build the right-hand side of the TDVP linear system.
!
! F_C(i) = <phi_i | H | psi>
! F_L(a) = <dpsi/dlambda_a | H | psi>
!
! If F_nD is requested, the routine returns the concatenated vector
! [F_C, F_L].
!=======================================================================
SUBROUTINE Calc_F_nD(psi, H, F_nD, FC_nD, F_L)
  IMPLICIT NONE

  TYPE(psi_t), INTENT(IN) :: psi
  TYPE(Op_t),  INTENT(IN) :: H
  COMPLEX(KIND=Rkind), ALLOCATABLE, OPTIONAL, INTENT(OUT) :: F_nD(:)
  COMPLEX(KIND=Rkind), ALLOCATABLE, OPTIONAL, INTENT(OUT) :: FC_nD(:)
  COMPLEX(KIND=Rkind), ALLOCATABLE, OPTIONAL, INTENT(OUT) :: F_L(:)

  TYPE(psi_t), ALLOCATABLE :: dAllPsi(:)
  TYPE(psi_t) :: Hpsi
  INTEGER :: ndim, nDeriv, nC, Ntot, i

  IF (.NOT. ASSOCIATED(psi%Basis)) THEN
     WRITE(out_unit,*) 'ERROR Calc_F_nD: psi%Basis is not associated.'
     STOP 'ERROR Calc_F_nD'
  END IF
  IF (.NOT. ALLOCATED(psi%CVec)) THEN
     WRITE(out_unit,*) 'ERROR Calc_F_nD: psi%CVec is not allocated.'
     STOP 'ERROR Calc_F_nD'
  END IF

  ndim   = SIZE(psi%Basis%tab_basis) - 1
  nDeriv = 4 * ndim
  nC     = SIZE(psi%CVec)
  Ntot   = nC + nDeriv

  CALL init_psi(Hpsi, psi%Basis, cplx=.TRUE., grid=.FALSE.)
  CALL calc_OpPsi(H, psi, Hpsi)
  CALL all_deriv_nD(psi, dAllPsi)

  IF (PRESENT(F_nD)) THEN
     ALLOCATE(F_nD(Ntot))
     F_nD = CZERO
     F_nD(1:nC) = Hpsi%CVec
     DO i = 1, nDeriv
        F_nD(nC+i) = DOT_PRODUCT(CONJG(dAllPsi(i)%CVec), Hpsi%CVec)
     END DO
  END IF

  IF (PRESENT(FC_nD)) THEN
     ALLOCATE(FC_nD(nC))
     FC_nD = Hpsi%CVec
  END IF

  IF (PRESENT(F_L)) THEN
     ALLOCATE(F_L(nDeriv))
     DO i = 1, nDeriv
        F_L(i) = DOT_PRODUCT(CONJG(dAllPsi(i)%CVec), Hpsi%CVec)
     END DO
  END IF

  CALL dealloc_psi(Hpsi)
  DO i = 1, nDeriv
     CALL dealloc_psi(dAllPsi(i))
  END DO
  DEALLOCATE(dAllPsi)

END SUBROUTINE Calc_F_nD

!=======================================================================
! Update the basis parameters using the TDVP derivative vector D_L.
!
! D_L contains time-derivatives of the real coordinates
! (q0, p0, Re(alpha), Im(alpha)).
! The actual Euler increment is
!   delta = -i * dt * D_L
! and each real coordinate is updated with the real part of its increment.
!=======================================================================
SUBROUTINE Update_Basis_From_Lambda(Basis, DeltaLambda_L, dt)
  IMPLICIT NONE

  TYPE(Basis_t), INTENT(INOUT)    :: Basis
  COMPLEX(KIND=Rkind), INTENT(IN) :: DeltaLambda_L(:)
  REAL(KIND=Rkind),    INTENT(IN) :: dt

  INTEGER :: ndim, k, idx
  COMPLEX(KIND=Rkind) :: dq, dp, da, db

  ndim = SIZE(Basis%tab_basis) - 1

  IF (SIZE(DeltaLambda_L) /= 4*ndim) THEN
     WRITE(out_unit,*) 'ERROR Update_Basis_From_Lambda: inconsistent lambda size.'
     STOP 'ERROR Update_Basis_From_Lambda'
  END IF

  DO k = 1, ndim
     idx = 4*(k-1)

     dq = (-EYE * dt) * DeltaLambda_L(idx + 1)
     dp = (-EYE * dt) * DeltaLambda_L(idx + 2)
     da = (-EYE * dt) * DeltaLambda_L(idx + 3)
     db = (-EYE * dt) * DeltaLambda_L(idx + 4)

     Basis%tab_basis(k)%Q0    = Basis%tab_basis(k)%Q0    + REAL(dq, KIND=Rkind)
     Basis%tab_basis(k)%Imp_k = Basis%tab_basis(k)%Imp_k + REAL(dp, KIND=Rkind)
     Basis%tab_basis(k)%Alpha = Basis%tab_basis(k)%Alpha + &
                                CMPLX(REAL(da, KIND=Rkind), REAL(db, KIND=Rkind), KIND=Rkind)

     IF (REAL(Basis%tab_basis(k)%Alpha, KIND=Rkind) <= ZERO) THEN
        WRITE(out_unit,*) 'ERROR Update_Basis_From_Lambda: Re(alpha) must remain positive.'
        STOP 'ERROR Update_Basis_From_Lambda'
     END IF
     Basis%tab_basis(k)%scaleQ = SQRT(REAL(Basis%tab_basis(k)%Alpha, KIND=Rkind))
  END DO

  CALL construct_primitive_basis_temp(Basis)

END SUBROUTINE Update_Basis_From_Lambda

!=======================================================================
! Two-step TDVP update compatible with the existing project structure.
!
! Step 1:
!   Build a first-order tangent update in the OLD basis
!      psi_tmp = psi + (-i dt) [ D_C + sum_a D_L(a) dpsi/dlambda_a ]
!
! Step 2:
!   Build the NEW basis from D_L, compute the overlap matrices between the
!   old and new Hagedorn bases, and project psi_tmp onto the new basis.
!
! This keeps the spirit of the existing three-step code while using the
! tangent-space information produced by the TDVP linear system.
!=======================================================================
SUBROUTINE Update_Psi_From_DeltaLambda(psi, psi_dt, D_C, D_L, dt)
  IMPLICIT NONE

  TYPE(psi_t), INTENT(IN)         :: psi
  TYPE(psi_t), INTENT(INOUT)      :: psi_dt
  COMPLEX(KIND=Rkind), INTENT(IN) :: D_C(:)
  COMPLEX(KIND=Rkind), INTENT(IN) :: D_L(:)
  REAL(KIND=Rkind),    INTENT(IN) :: dt

  TYPE(psi_t), ALLOCATABLE :: dAllPsi(:)
  TYPE(psi_t)              :: psi_tmp, psi_proj
  TYPE(Basis_t)            :: Basis_new
  INTEGER :: i, nC, nDeriv

  IF (.NOT. ASSOCIATED(psi%Basis)) THEN
     WRITE(out_unit,*) 'ERROR Update_Psi_From_DeltaLambda: psi%Basis is not associated.'
     STOP 'ERROR Update_Psi_From_DeltaLambda'
  END IF
  IF (.NOT. ASSOCIATED(psi_dt%Basis)) THEN
     WRITE(out_unit,*) 'ERROR Update_Psi_From_DeltaLambda: psi_dt%Basis is not associated.'
     STOP 'ERROR Update_Psi_From_DeltaLambda'
  END IF
  IF (.NOT. ALLOCATED(psi%CVec)) THEN
     WRITE(out_unit,*) 'ERROR Update_Psi_From_DeltaLambda: psi%CVec is not allocated.'
     STOP 'ERROR Update_Psi_From_DeltaLambda'
  END IF
  IF (.NOT. ALLOCATED(psi_dt%CVec)) THEN
     WRITE(out_unit,*) 'ERROR Update_Psi_From_DeltaLambda: psi_dt%CVec must be allocated by caller.'
     STOP 'ERROR Update_Psi_From_DeltaLambda'
  END IF

  nC     = SIZE(psi%CVec)
  nDeriv = SIZE(D_L)

  IF (SIZE(D_C) /= nC) THEN
     WRITE(out_unit,*) 'ERROR Update_Psi_From_DeltaLambda: wrong size for D_C.'
     STOP 'ERROR Update_Psi_From_DeltaLambda'
  END IF

  CALL init_psi(psi_tmp, psi%Basis, cplx=.TRUE., grid=.FALSE.)
  psi_tmp%CVec = psi%CVec + (-EYE * dt) * D_C

  CALL all_deriv_nD(psi, dAllPsi)
  IF (SIZE(dAllPsi) /= nDeriv) THEN
     WRITE(out_unit,*) 'ERROR Update_Psi_From_DeltaLambda: wrong number of tangent vectors.'
     STOP 'ERROR Update_Psi_From_DeltaLambda'
  END IF

  DO i = 1, nDeriv
     psi_tmp%CVec = psi_tmp%CVec + (-EYE * dt) * D_L(i) * dAllPsi(i)%CVec
  END DO

  CALL init_Basis1_TO_Basis2(Basis_new, psi%Basis)
  CALL Update_Basis_From_Lambda(Basis_new, D_L, dt)
  CALL Build_Overlap_Matrices_Local(psi%Basis, Basis_new)

  CALL init_psi(psi_proj, Basis_new, cplx=.TRUE., grid=.FALSE.)
  psi_proj%CVec = CZERO
  CALL Projection(psi_proj, psi_tmp)

  CALL Update_Basis_From_Lambda(psi_dt%Basis, D_L, dt)
  psi_dt%CVec = psi_proj%CVec

  DO i = 1, nDeriv
     CALL dealloc_psi(dAllPsi(i))
  END DO
  DEALLOCATE(dAllPsi)

  CALL dealloc_psi(psi_tmp)
  CALL dealloc_psi(psi_proj)
  CALL Deallocate_Basis(Basis_new)

END SUBROUTINE Update_Psi_From_DeltaLambda

!=======================================================================
! Diagnostic helper.
!=======================================================================
SUBROUTINE Test_Matrices(psi, H)
  IMPLICIT NONE

  TYPE(psi_t), INTENT(IN) :: psi
  TYPE(Op_t),  INTENT(IN) :: H
  COMPLEX(KIND=Rkind), ALLOCATABLE :: S(:,:), M_Clambda(:,:), M_lambdaLambda(:,:), M_nD(:,:), F_nD(:)

  CALL Calc_S_MClambda_MLambdaLambda(psi, S, M_Clambda, M_lambdaLambda)
  CALL Calc_M_nD(psi, M_nD)
  CALL Calc_F_nD(psi, H, F_nD)

  WRITE(out_unit,*) 'Overlap block S:'
  CALL Write_VecMat(S, out_unit, 6, info='S')

  WRITE(out_unit,*) 'Coefficient/parameter block M_Clambda:'
  CALL Write_VecMat(M_Clambda, out_unit, 6, info='M_Clambda')

  WRITE(out_unit,*) 'Parameter block M_lambdaLambda:'
  CALL Write_VecMat(M_lambdaLambda, out_unit, 6, info='M_lambdaLambda')

  WRITE(out_unit,*) 'Full TDVP matrix M_nD:'
  CALL Write_VecMat(M_nD, out_unit, 6, info='M_nD')

  WRITE(out_unit,*) 'Right-hand-side F_nD:'
  CALL Write_VecMat(F_nD, out_unit, 6, info='F_nD')

  DEALLOCATE(S, M_Clambda, M_lambdaLambda, M_nD, F_nD)

END SUBROUTINE Test_Matrices

!=======================================================================
! Assemble the full TDVP metric.
!=======================================================================
SUBROUTINE Calc_M_nD(psi, M_nD)
  IMPLICIT NONE

  TYPE(psi_t), INTENT(IN) :: psi
  COMPLEX(KIND=Rkind), ALLOCATABLE, INTENT(OUT) :: M_nD(:,:)

  COMPLEX(KIND=Rkind), ALLOCATABLE :: S(:,:), M_Cl(:,:), M_Ll(:,:)
  INTEGER :: nC, Ntot, Nmat

  IF (.NOT. ASSOCIATED(psi%Basis)) THEN
     WRITE(out_unit,*) 'ERROR Calc_M_nD: psi%Basis is not associated.'
     STOP 'Calc_M_nD'
  END IF
  IF (.NOT. ALLOCATED(psi%CVec)) THEN
     WRITE(out_unit,*) 'ERROR Calc_M_nD: psi%CVec is not allocated.'
     STOP 'Calc_M_nD'
  END IF

  CALL Calc_S_MClambda_MLambdaLambda(psi, S, M_Cl, M_Ll)

  nC   = SIZE(S,1)
  Ntot = SIZE(M_Cl,2)
  Nmat = nC + Ntot

  ALLOCATE(M_nD(Nmat, Nmat))
  M_nD = CZERO

  M_nD(1:nC, 1:nC)         = S
  M_nD(1:nC, nC+1:Nmat)    = M_Cl
  M_nD(nC+1:Nmat, 1:nC)    = CONJG(TRANSPOSE(M_Cl))
  M_nD(nC+1:Nmat, nC+1:Nmat) = M_Ll

  DEALLOCATE(S, M_Cl, M_Ll)

END SUBROUTINE Calc_M_nD

!=======================================================================
! Transform between the explicit split representation
!   Lambda = [ C ; q0 ; p0 ; Re(alpha) ; Im(alpha) ]
! and the pair (psi coefficients, basis parameters).
!=======================================================================
SUBROUTINE Psi_Lambda_Transform(psi, Lambda, toLambda)
  IMPLICIT NONE

  TYPE(psi_t), INTENT(INOUT) :: psi
  COMPLEX(KIND=Rkind), ALLOCATABLE, INTENT(INOUT) :: Lambda(:)
  LOGICAL, INTENT(IN) :: toLambda

  INTEGER :: nb, ndof, i, pos, expected_size

  IF (.NOT. ASSOCIATED(psi%Basis)) THEN
     WRITE(out_unit,*) 'ERROR Psi_Lambda_Transform: psi%Basis is not associated.'
     STOP 'Psi_Lambda_Transform'
  END IF
  IF (.NOT. ALLOCATED(psi%CVec)) THEN
     WRITE(out_unit,*) 'ERROR Psi_Lambda_Transform: psi%CVec is not allocated.'
     STOP 'Psi_Lambda_Transform'
  END IF

  nb   = SIZE(psi%CVec)
  ndof = SIZE(psi%Basis%tab_basis) - 1
  expected_size = nb + 4*ndof

  IF (toLambda) THEN
     IF (ALLOCATED(Lambda)) THEN
        IF (SIZE(Lambda) /= expected_size) DEALLOCATE(Lambda)
     END IF
     IF (.NOT. ALLOCATED(Lambda)) ALLOCATE(Lambda(expected_size))

     Lambda(1:nb) = psi%CVec

     pos = nb
     DO i = 1, ndof
        pos = pos + 1; Lambda(pos) = CMPLX(psi%Basis%tab_basis(i)%Q0,    ZERO, KIND=Rkind)
        pos = pos + 1; Lambda(pos) = CMPLX(psi%Basis%tab_basis(i)%Imp_k, ZERO, KIND=Rkind)
        pos = pos + 1; Lambda(pos) = CMPLX(REAL(psi%Basis%tab_basis(i)%Alpha,  KIND=Rkind), ZERO, KIND=Rkind)
        pos = pos + 1; Lambda(pos) = CMPLX(AIMAG(psi%Basis%tab_basis(i)%Alpha), ZERO, KIND=Rkind)
     END DO

  ELSE
     IF (.NOT. ALLOCATED(Lambda)) THEN
        WRITE(out_unit,*) 'ERROR Psi_Lambda_Transform: Lambda is not allocated.'
        STOP 'Psi_Lambda_Transform'
     END IF
     IF (SIZE(Lambda) /= expected_size) THEN
        WRITE(out_unit,*) 'ERROR Psi_Lambda_Transform: invalid Lambda size.'
        STOP 'Psi_Lambda_Transform'
     END IF

     psi%CVec = Lambda(1:nb)

     pos = nb
     DO i = 1, ndof
        pos = pos + 1; psi%Basis%tab_basis(i)%Q0    = REAL(Lambda(pos), KIND=Rkind)
        pos = pos + 1; psi%Basis%tab_basis(i)%Imp_k = REAL(Lambda(pos), KIND=Rkind)
        pos = pos + 1
        psi%Basis%tab_basis(i)%Alpha = CMPLX( &
             REAL(Lambda(pos), KIND=Rkind), &
             AIMAG(psi%Basis%tab_basis(i)%Alpha), KIND=Rkind)
        pos = pos + 1
        psi%Basis%tab_basis(i)%Alpha = CMPLX( &
             REAL(psi%Basis%tab_basis(i)%Alpha, KIND=Rkind), &
             REAL(Lambda(pos), KIND=Rkind), KIND=Rkind)

        IF (REAL(psi%Basis%tab_basis(i)%Alpha, KIND=Rkind) <= ZERO) THEN
           WRITE(out_unit,*) 'ERROR Psi_Lambda_Transform: Re(alpha) must remain positive.'
           STOP 'Psi_Lambda_Transform'
        END IF
        psi%Basis%tab_basis(i)%scaleQ = SQRT(REAL(psi%Basis%tab_basis(i)%Alpha, KIND=Rkind))
     END DO

     CALL construct_primitive_basis_temp(psi%Basis)
  END IF

END SUBROUTINE Psi_Lambda_Transform

!=======================================================================
! Round-trip test of Psi_Lambda_Transform.
!=======================================================================
SUBROUTINE Test_Psi_Lambda_Transform(psi)
  IMPLICIT NONE

  TYPE(psi_t), INTENT(INOUT) :: psi
  TYPE(psi_t)                :: psi2
  COMPLEX(KIND=Rkind), ALLOCATABLE :: Lambda(:)
  INTEGER :: i

  WRITE(out_unit, *) '>>> Test_Psi_Lambda_Transform'

  CALL Psi_Lambda_Transform(psi, Lambda, .TRUE.)

  WRITE(out_unit, *) 'Lambda vector (from psi):'
  DO i = 1, SIZE(Lambda)
     WRITE(out_unit, '(I3,2X,F12.6,2X,F12.6)') i, REAL(Lambda(i)), AIMAG(Lambda(i))
  END DO

  CALL init_psi(psi2, psi%Basis, .TRUE., .FALSE.)
  CALL Psi_Lambda_Transform(psi2, Lambda, .FALSE.)

  WRITE(out_unit, *) 'psi%CVec (original):'
  DO i = 1, SIZE(psi%CVec)
     WRITE(out_unit, '(I3,2X,F12.6,2X,F12.6)') i, REAL(psi%CVec(i)), AIMAG(psi%CVec(i))
  END DO

  WRITE(out_unit, *) 'psi2%CVec (reconstructed):'
  DO i = 1, SIZE(psi2%CVec)
     WRITE(out_unit, '(I3,2X,F12.6,2X,F12.6)') i, REAL(psi2%CVec(i)), AIMAG(psi2%CVec(i))
  END DO

  WRITE(out_unit, *) '<<< End Test_Psi_Lambda_Transform'

  CALL dealloc_psi(psi2)
  IF (ALLOCATED(Lambda)) DEALLOCATE(Lambda)

END SUBROUTINE Test_Psi_Lambda_Transform

!=======================================================================
! Build overlap matrices between two direct-product bases.
! The overlap for each moving Hagedorn/HO factor is stored in B_old%tab_basis(i)%S.
! For static factors (Fourier, box, electronic), the identity is used.
!=======================================================================
SUBROUTINE Build_Overlap_Matrices_Local(B_old, B_new)
  IMPLICIT NONE

  TYPE(Basis_t), INTENT(INOUT) :: B_old
  TYPE(Basis_t), INTENT(IN)    :: B_new

  INTEGER :: ib, ndim, nb, nq
  REAL(KIND=Rkind) :: Q0, SQ0, P0, Qt, SQt, Pt
  COMPLEX(KIND=Rkind) :: A0, At
  CHARACTER(LEN=:), ALLOCATABLE :: bname

  ndim = SIZE(B_old%tab_basis)

  DO ib = 1, ndim
     bname = TRIM(B_old%tab_basis(ib)%Basis_name)
     nb = B_old%tab_basis(ib)%nb

     IF (bname == 'el') THEN
        CALL Build_Identity_Local(B_old%tab_basis(ib)%S, nb)
        CYCLE
     END IF

     IF (bname == 'hag' .OR. bname == 'ho') THEN
        nq  = B_old%tab_basis(ib)%nq
        Q0  = B_old%tab_basis(ib)%Q0
        SQ0 = B_old%tab_basis(ib)%scaleQ
        A0  = B_old%tab_basis(ib)%Alpha
        P0  = B_old%tab_basis(ib)%Imp_k

        Qt  = B_new%tab_basis(ib)%Q0
        SQt = B_new%tab_basis(ib)%scaleQ
        At  = B_new%tab_basis(ib)%Alpha
        Pt  = B_new%tab_basis(ib)%Imp_k

        CALL Calc_Cross_Overlap_Hagedorn(B_old%tab_basis(ib)%S, nb, nq, Q0, SQ0, A0, P0, Qt, SQt, At, Pt)
     ELSE
        CALL Build_Identity_Local(B_old%tab_basis(ib)%S, nb)
     END IF
  END DO

END SUBROUTINE Build_Overlap_Matrices_Local

!=======================================================================
! Identity matrix helper.
!=======================================================================
SUBROUTINE Build_Identity_Local(S, n)
  IMPLICIT NONE

  COMPLEX(KIND=Rkind), ALLOCATABLE, INTENT(INOUT) :: S(:,:)
  INTEGER, INTENT(IN) :: n
  INTEGER :: i

  IF (ALLOCATED(S)) DEALLOCATE(S)
  ALLOCATE(S(n,n))
  S = CZERO
  DO i = 1, n
     S(i,i) = CONE
  END DO

END SUBROUTINE Build_Identity_Local

!=======================================================================
! Numerical overlap between two 1D Hagedorn/HO bases.
! This is the same quantity needed by Projection(old -> new).
!=======================================================================
SUBROUTINE Calc_Cross_Overlap_Hagedorn(S, nb, nq, Q0, SQ0, A0, P0, Qt, SQt, At, Pt)
  IMPLICIT NONE

  COMPLEX(KIND=Rkind), ALLOCATABLE, INTENT(INOUT) :: S(:,:)
  INTEGER,             INTENT(IN) :: nb, nq
  REAL(KIND=Rkind),    INTENT(IN) :: Q0, SQ0, P0, Qt, SQt, Pt
  COMPLEX(KIND=Rkind), INTENT(IN) :: A0, At

  INTEGER :: ib, iq
  REAL(KIND=Rkind) :: SQeq, Qeq
  REAL(KIND=Rkind), ALLOCATABLE :: Q(:), W(:)
  COMPLEX(KIND=Rkind), ALLOCATABLE :: d0gb0(:,:), d0gb1(:,:), d0bgw(:,:)

  SQeq = SQRT(SQ0*SQ0 + SQt*SQt) / SQRT(TWO)
  Qeq  = (SQ0*SQ0*Q0 + SQt*SQt*Qt) / (SQ0*SQ0 + SQt*SQt)

  ALLOCATE(Q(nq), W(nq))
  CALL hercom(nq, Q, W)

  W = W / SQeq
  Q = Qeq + Q / SQeq

  ALLOCATE(d0gb0(nq, nb), d0gb1(nq, nb), d0bgw(nb, nq))

  DO iq = 1, nq
     DO ib = 1, nb
        CALL Eval_Hagedorn_Function(Q(iq), Q0, A0, P0, ib-1, d0gb0(iq, ib))
        CALL Eval_Hagedorn_Function(Q(iq), Qt, At, Pt, ib-1, d0gb1(iq, ib))
     END DO
  END DO

  d0bgw = TRANSPOSE(d0gb1)
  DO ib = 1, nb
     d0bgw(ib,:) = d0bgw(ib,:) * W(:)
  END DO

  IF (ALLOCATED(S)) DEALLOCATE(S)
  ALLOCATE(S(nb, nb))
  S = MATMUL(CONJG(d0bgw), d0gb0)

  DEALLOCATE(d0gb0, d0gb1, d0bgw, Q, W)

END SUBROUTINE Calc_Cross_Overlap_Hagedorn

!=======================================================================
! Evaluate one Hagedorn/HO basis function.
!=======================================================================
SUBROUTINE Eval_Hagedorn_Function(Q, Qt, At, Pt, ib, d0)
  IMPLICIT NONE

  REAL(KIND=Rkind),    INTENT(IN)  :: Q, Qt, Pt
  COMPLEX(KIND=Rkind), INTENT(IN)  :: At
  INTEGER,             INTENT(IN)  :: ib
  COMPLEX(KIND=Rkind), INTENT(OUT) :: d0

  COMPLEX(KIND=Rkind) :: phase
  REAL(KIND=Rkind) :: SQ, DQ

  SQ = SQRT(REAL(At, KIND=Rkind))
  DQ = SQ * (Q - Qt)

  phase = EXP(-HALF * EYE * AIMAG(At) * (Q - Qt)**2 + EYE * Pt * (Q - Qt))
  d0 = SQRT(SQ) * poly_Hermite(DQ, ib) * EXP(-HALF * DQ*DQ) * phase

END SUBROUTINE Eval_Hagedorn_Function

END MODULE Sub_Vp_m
