SUBROUTINE Calc_dngg_grid(Basis, analytic_d1bb)
   USE QDUtil_m
   IMPLICIT NONE
   TYPE(Basis_t), INTENT(INOUT) :: Basis
   LOGICAL, OPTIONAL, INTENT(IN) :: analytic_d1bb
   INTEGER :: nb, nq, ib, m, n
   LOGICAL :: use_analytic
   LOGICAL, PARAMETER :: debug = .TRUE.
   REAL(kind=Rkind) :: a, b, p, q0
   CHARACTER(len=Name_len) :: name
 
   IF (debug) THEN
      WRITE(out_unit,*) 'BEGINNING Calc_dngg_grid'
      WRITE(out_unit,*) 'Basis name: ', TRIM(Basis%Basis_name)
      WRITE(out_unit,*) 'nb: ', Basis%nb
      WRITE(out_unit,*) 'nq: ', Basis%nq
      FLUSH(out_unit)
   END IF
 

    name = Basis%Basis_name 
   call string_uppercase_TO_lowercase(name)
   !-------------------------------------------------------------------------
   ! Determine calculation method for d1bb
   IF (PRESENT(analytic_d1bb)) THEN
      use_analytic = analytic_d1bb
      IF (debug) THEN
         WRITE(out_unit,*) 'd1bb calculation method: ', &
              MERGE('ANALYTIC', 'NUMERIC ', use_analytic)
      END IF
   ELSE
      use_analytic = .FALSE. ! Default: numerical calculation
      IF (debug) THEN
         WRITE(out_unit,*) 'd1bb calculation method: NUMERIC (default)'
      END IF
   END IF
 
   !-------------------------------------------------------------------------
   ! Retrieve dimensions
   nb = Basis%nb
   nq = Basis%nq
 
   !-------------------------------------------------------------------------
   ! Deallocate previous arrays if allocated
   IF (ALLOCATED(Basis%d1gg)) DEALLOCATE(Basis%d1gg)
   IF (ALLOCATED(Basis%d2gg)) DEALLOCATE(Basis%d2gg)
   IF (ALLOCATED(Basis%d1bb)) DEALLOCATE(Basis%d1bb)
 
   !-------------------------------------------------------------------------
   ! Allocate arrays
   ALLOCATE(Basis%d1gg(nq, nq, 1))
   ALLOCATE(Basis%d2gg(nq, nq, 1, 1))
 
   !-------------------------------------------------------------------------
   ! Skip computation for 'el' basis type
   IF (name == 'el') RETURN
 
   !-------------------------------------------------------------------------
   ! Grid-grid matrices
   Basis%d1gg(:,:,1)    = MATMUL(Basis%d1gb(:,:,1), CONJG(Basis%d0bgw))
   Basis%d2gg(:,:,1,1)  = MATMUL(Basis%d2gb(:,:,1,1), CONJG(Basis%d0bgw))
 
   IF(name == 'ho') THEN
      !-------------------------------------------------------------------------
      ! Basis-basis matrices (weighted overlaps)
      ! d1bb(:,:,k) = < b_m | ∂_k b_n >
      ALLOCATE(Basis%d1bb(nb, nb, 5))
      
      IF (use_analytic) THEN
         ! Analytical calculation of d1bb  
         Basis%d1bb(:,:,:) = CZERO
         
         ! Retrieve basis parameters
         a = REAL(Basis%alpha,kind=Rkind)
         b = IMG(Basis%alpha)
         p = Basis%imp_k
         q0 = Basis%Q0
         
         DO n = 1, nb
            DO m = 1, nb
               ! d1bb(m,n,1) = < b_m | ∂_x b_n > (not provided in analytical form)
               ! Will be calculated numerically if needed
               
               ! d1bb(m,n,2) = < b_m | ∂_{q0} b_n >
               IF (m == n) THEN
                  Basis%d1bb(m,n,2) = -EYE * p
               ELSE IF (m == n+1) THEN
                  Basis%d1bb(m,n,2) = (a + EYE*b) / SQRT(TWO * a * n)
               ELSE IF (m == n-1) THEN
                  Basis%d1bb(m,n,2) = -(a - EYE*b) / SQRT(TWO * a * (n-1))
               END IF
               
               ! d1bb(m,n,3) = < b_m | ∂_p b_n >
               IF (m == n-1) THEN
                  Basis%d1bb(m,n,3) = (EYE / SQRT(a)) * SQRT((n-ONE)/TWO)
               ELSE IF (m == n+1) THEN
                  Basis%d1bb(m,n,3) = (EYE / SQRT(a)) * SQRT(n/TWO)
               END IF
               
               ! d1bb(m,n,4) = < b_m | ∂_a b_n >
               IF (m == n-2) THEN
                  Basis%d1bb(m,n,4) = SQRT((n-ONE)*(n-TWO)) / (FOUR * a)
               ELSE IF (m == n+2) THEN
                  Basis%d1bb(m,n,4) = -SQRT(n*(n+ONE)) / (FOUR * a)
               END IF
               
               ! d1bb(m,n,5) = < b_m | ∂_b b_n >
               IF (m == n-2) THEN
                  Basis%d1bb(m,n,5) = -EYE * SQRT((n-ONE)*(n-TWO)) / (FOUR * a)
               ELSE IF (m == n) THEN
                  Basis%d1bb(m,n,5) = -EYE *(TWO*n - ONE) / (FOUR * a)
               ELSE IF (m == n+2) THEN
                  Basis%d1bb(m,n,5) = -EYE * SQRT(n*(n+ONE)) / (FOUR * a)
               END IF
            END DO
         END DO
         
      ELSE
         ! Numerical calculation (original method)
         Basis%d1bb(:,:,1) = MATMUL(Basis%d0bgw, Basis%d1gb(:,:,1)) ! d/dx
         Basis%d1bb(:,:,2) = MATMUL(Basis%d0bgw, Basis%d1gb(:,:,2)) ! ∂/∂q0
         Basis%d1bb(:,:,3) = MATMUL(Basis%d0bgw, Basis%d1gb(:,:,3)) ! ∂/∂p
         Basis%d1bb(:,:,4) = MATMUL(Basis%d0bgw, Basis%d1gb(:,:,4)) ! ∂/∂a
         Basis%d1bb(:,:,5) = MATMUL(Basis%d0bgw, Basis%d1gb(:,:,5)) ! ∂/∂b
      END IF
   END IF
 
   !-------------------------------------------------------------------------
   ! Debug output
   IF (debug) THEN
      CALL Write_VecMat(Basis%d1gg(:,:,1), out_unit, 5, info='d1gg')
      WRITE(out_unit,*)
      CALL Write_VecMat(Basis%d2gg(:,:,1,1), out_unit, 5, info='d2gg')
      IF(Basis%Basis_name == 'ho') THEN
         DO ib = 2, 5
            WRITE(out_unit,*) 'd1bb(:,:,', ib, ') - ', &
                 MERGE('ANALYTIC', 'NUMERIC ', use_analytic)
            CALL Write_VecMat(Basis%d1bb(:,:,ib), out_unit, 5, info='d1bb')
         END DO
      END IF
      WRITE(out_unit,*) 'END Calc_dngg_grid'
      FLUSH(out_unit)
   END IF
 
 END SUBROUTINE Calc_dngg_grid




 SUBROUTINE march_VP(psi, psi_dt, t, propa, H)
   USE QDUtil_m
   IMPLICIT NONE
   TYPE(psi_t),    INTENT(IN)    :: psi
   TYPE(psi_t),    INTENT(INOUT) :: psi_dt
   TYPE(propa_t),  INTENT(IN)    :: propa
   TYPE(Op_t),     INTENT(IN)    :: H
   REAL(KIND=Rkind), INTENT(IN)  :: t
 
   ! local variables
   COMPLEX(KIND=Rkind), ALLOCATABLE   :: Lambda(:), Lambda_new(:)
   TYPE(psi_t)                        :: tmp
   COMPLEX(KIND=Rkind), ALLOCATABLE   :: M(:,:), F(:), D_L(:)
   REAL(KIND=Rkind)                   :: Norm0, Norm
   TYPE(psi_t), ALLOCATABLE           :: dAllPsi(:)
   INTEGER :: n, nC, nDeriv, ib, ndim
   LOGICAL :: error_flag = .FALSE.
 
   LOGICAL, PARAMETER                 :: debug = .TRUE.
 
   CALL init_psi(tmp, psi%Basis, cplx=.TRUE., grid=.TRUE.)
   tmp%CVec = psi%CVec

   !------------------------------------------------------------
   ! Step 1. Convert psi -> Lambda
   !------------------------------------------------------------
   CALL Psi_Lambda_Transform(tmp, Lambda, .TRUE.)
   n = SIZE(Lambda)
   ndim = SIZE(psi%Basis%tab_basis) - 1
   nDeriv = 4 * ndim
   nC = SIZE(psi%CVec)
   IF (ALLOCATED(dAllPsi)) DEALLOCATE(dAllPsi)
   ALLOCATE(dAllPsi(nDeriv))
   CALL Calc_Norm_OF_Psi(psi, Norm0)

   !------------------------------------------------------------
   ! Step 2. Build M and F at time t
   !------------------------------------------------------------
   CALL Calc_M_nD(psi=psi, M_nD=M)     ! generalized metric matrix
   CALL Calc_F_nD(psi=psi, H=H, F_nD=F)! force-like vector
   CALL all_deriv_nD(psi, dAllPsi)

   !------------------------------------------------------------
   ! Step 3. Solve ΔΛ = -i dt M^(-1) F
   !------------------------------------------------------------
   ALLOCATE(Lambda_new(n))
   Lambda_new = Lambda

   M = inv_OF_Mat(M, inv_type=0, epsi=1.d-20)   ! M = (M)^(-1)

   ! Update Lambda
   Lambda_new = Lambda - EYE * propa%delta_t * MATMUL(M,F)
   ALLOCATE(D_L(nDeriv))
   D_L(:) = Lambda_new(nC+1:n)

   DO ib = 1, nDeriv
      IF (.NOT. ALLOCATED(dAllPsi(ib)%CVec)) THEN
         WRITE(*,*) "ERROR in march_VP: derivative ", ib, " not allocated"
         error_flag = .TRUE.
         EXIT
      END IF
      Lambda_new(1:nC) = Lambda_new(1:nC) + (-EYE * propa%delta_t) * D_L(ib) * dAllPsi(ib)%CVec(:)
   END DO

   !------------------------------------------------------------
   ! Step 4. Transform back: Lambda -> psi_dt (if no error)
   !------------------------------------------------------------
   IF (.NOT. error_flag) THEN
      CALL Psi_Lambda_Transform(psi_dt, Lambda_new, .FALSE.)
   END IF

   !------------------------------------------------------------
   ! Debug info
   !------------------------------------------------------------
   IF (debug .AND. .NOT. error_flag) THEN
      WRITE(out_unit,*) "BEGIN march_VP at t=", t, " dt=", propa%delta_t
      CALL Calc_Norm_OF_Psi(psi_dt, Norm)
      WRITE(out_unit, *) '<psi_dt|psi_dt> = ', Norm, ' |ΔNorm| = ', abs(Norm0 - Norm)
      WRITE(out_unit,*) "END march_VP"
      FLUSH(out_unit)
   END IF

   !------------------------------------------------------------
   ! Clean up
   !------------------------------------------------------------
   CALL dealloc_psi(tmp)
   DEALLOCATE(Lambda, Lambda_new, M, F, D_L)
   DO ib = 1, nDeriv
      IF (ALLOCATED(dAllPsi(ib)%CVec)) CALL dealloc_psi(dAllPsi(ib))
   END DO
   IF (ALLOCATED(dAllPsi)) DEALLOCATE(dAllPsi)

   RETURN
END SUBROUTINE march_VP 
