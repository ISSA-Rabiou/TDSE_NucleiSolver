!----------------------------------------------------------------------------------
!> MODULE: lanczos_m
!> Description: This module implements the Short Iterative Lanczos (SIL) method
!>              for the time propagation of the Schrödinger equation.
!>              It uses Lanczos iterations to build a tridiagonal Hamiltonian
!>              matrix in a Krylov basis, then diagonalizes this small matrix
!>              to compute the propagation operator.
!----------------------------------------------------------------------------------
MODULE lanczos_m
  USE QDUtil_m           
  USE Op_m
  USE psi_m
  USE Basis_m

  IMPLICIT NONE
  PRIVATE
  
  ! Public interface declarations
  PUBLIC :: march_mSIL, find_ground_state_SIL

CONTAINS


    !----------------------------------------------------------------------------------
    !> SUBROUTINE: gram_schmidt
    !> Description: Orthogonalizes and normalizes a new vector `q_new` with respect to a given
    !>              set of orthonormal vectors `Q`. The vector `q_new` is modified
    !>              in place. This is a robust method to maintain orthogonality.
    !> Arguments:
    !>   - Q      [in]      (complex(kind=Rkind), allocatable) : Orthonormal basis vectors.
    !>   - q_new  [inout]   (complex(kind=Rkind), allocatable) : Vector to be orthogonalized and normalized.
    !>   - k      [in]      (integer) : Number of vectors in the basis `Q`.
    !----------------------------------------------------------------------------------
    SUBROUTINE gram_schmidt(Q, q_new, k)
      IMPLICIT NONE
      COMPLEX(KIND=Rkind), INTENT(IN)     :: Q(:,:)
      COMPLEX(KIND=Rkind), INTENT(INOUT)  :: q_new(:)
      INTEGER, INTENT(IN)                 :: k
      REAL(KIND=Rkind)                    :: COEFF1,COEFF2  
    
      ! Local variables
      INTEGER :: j
    
      ! Orthogonalization loop: q_new ⟂ {Q(:,1)...Q(:,k)}
      DO j = 1, k
        CALL Dotproduct(Q(:,j), q_new(:),COEFF1)
        !CALL Dotproduct(Q(:,j), Q(:,j),COEFF2)
        q_new(:) = q_new(:) - COEFF1 * Q(:,j)
      END DO
    
    END SUBROUTINE gram_schmidt


!----------------------------------------------------------------
! Build tridiagonal matrix T from alpha/beta and diagonalize it.
! - alpha(1:nb_k) : diagonal elements (REAL)
! - beta(2:nb_k)  : off-diagonal elements (REAL), beta(1) unused
! Outputs:
! - eigvals(1:nb_k)  : REAL eigenvalues
! - eigvecs(1:nb_k,1:nb_k) : REAL eigenvectors (columns)
! - info : integer status (0 = ok)
!----------------------------------------------------------------
SUBROUTINE build_tridiag(alpha, beta, nb_k, eigvals, eigvecs, info)
  IMPLICIT NONE
  INTEGER, INTENT(IN)                   :: nb_k
  REAL(KIND=Rkind), INTENT(IN)          :: alpha(:)
  REAL(KIND=Rkind), INTENT(IN)          :: beta(:)
  REAL(KIND=Rkind), INTENT(OUT)         :: eigvals(:)
  REAL(KIND=Rkind), INTENT(OUT)         :: eigvecs(:,:)
  INTEGER, INTENT(OUT)                  :: info

  REAL(KIND=Rkind), ALLOCATABLE         :: T(:,:)
  INTEGER                               :: i
  LOGICAL, PARAMETER                    :: debug=.FALSE.

  info = 0

  IF (SIZE(alpha) < nb_k .OR. SIZE(eigvals) < nb_k .OR. &
      SIZE(eigvecs,1) < nb_k .OR. SIZE(eigvecs,2) < nb_k) THEN
    WRITE(out_unit, *) 'Error: build_tridiag - array sizes too small.'
    info = -1
    RETURN
  END IF

  ALLOCATE(T(nb_k, nb_k))
  T = ZERO

  DO i = 1, nb_k
    T(i,i) = alpha(i)
    IF (i < nb_k) THEN
      T(i,i+1) = beta(i+1)
      T(i+1,i) = beta(i+1)
    END IF
  END DO

  ! Optional: print T for debug
  IF (debug) THEN
    CALL Write_VecMat(T, out_unit, nb_k/2, info='T (tridiagonal)')
    WRITE(out_unit, *) " "
  END IF

  ! Diagonalisation : utilise ta routine existante 'diagonalization'
  ! signature attendue : diagonalization(T, eigvals, eigvecs, nb_k)
  CALL diagonalization(T, eigvals, eigvecs, nb_k)

  ! NOTE: diagonalization doit remplir eigvals (REAL) et eigvecs (REAL, colonnes).
  ! Si elle renvoie un code d'erreur, remonte-le via info (ici on suppose success)
  info = 0

  DEALLOCATE(T)
END SUBROUTINE build_tridiag


!----------------------------------------------------------------
! Build propagator U = V * diag(exp(-i*eigvals*dt)) * V^T
! - eigvals(1:nb_k) : REAL eigenvalues
! - eigvecs(1:nb_k,1:nb_k) : REAL eigenvectors (columns)
! Outputs:
! - U_prop(1:nb_k,1:nb_k) : COMPLEX propagator in Krylov subspace
!----------------------------------------------------------------
SUBROUTINE build_propagator(eigvals, eigvecs, nb_k, dt, U_prop)
  IMPLICIT NONE
  INTEGER, INTENT(IN)                    :: nb_k
  REAL(KIND=Rkind), INTENT(IN)           :: eigvals(:)
  REAL(KIND=Rkind), INTENT(IN)           :: eigvecs(:,:)
  REAL(KIND=Rkind), INTENT(IN)           :: dt
  COMPLEX(KIND=Rkind), INTENT(OUT)       :: U_prop(:,:)
 
  COMPLEX(KIND=Rkind), ALLOCATABLE       :: tmpc(:,:)
  INTEGER                                :: j
  LOGICAL, PARAMETER                     :: debug=.FALSE.


  IF (SIZE(eigvals) < nb_k .OR. SIZE(eigvecs,1) < nb_k .OR. &
      SIZE(eigvecs,2) < nb_k .OR. SIZE(U_prop,1) < nb_k .OR. SIZE(U_prop,2) < nb_k) THEN
    WRITE(out_unit, *) 'Error: build_propagator - array sizes too small.'
    RETURN
  END IF

  ALLOCATE(tmpc(nb_k, nb_k))
  tmpc = CZERO
  U_prop = CZERO

  ! Convert columns of eigvecs to complex and scale by the exponential factor:
  ! tmpc(:,j) = V(:,j) * exp(-i * eigvals(j) * dt)
  DO j = 1, nb_k
    tmpc(:, j) = CMPLX(eigvecs(:, j), ZERO, KIND=Rkind) * EXP(-EYE * eigvals(j) * dt)
  END DO

  ! U = tmpc * V^T  where V^T = CONJG(TRANSPOSE(CMPLX(eigvecs,0)))
  U_prop = MATMUL(tmpc, CONJG(TRANSPOSE(CMPLX(eigvecs, ZERO, KIND=Rkind))))

  ! Debug printing
  IF (debug) THEN
    CALL Write_VecMat(U_prop, out_unit, nb_k/2, info='U_prop (Krylov propagator)')
    WRITE(out_unit, *) " "
  END IF

  DEALLOCATE(tmpc)
END SUBROUTINE build_propagator

 SUBROUTINE check_unitarity(U, n)
   IMPLICIT NONE
   INTEGER, INTENT(IN)                   :: n
   COMPLEX(KIND=Rkind), INTENT(IN)       :: U(:,:)
 
   COMPLEX(KIND=Rkind), ALLOCATABLE      :: M(:,:), Iden(:,:)
   REAL(KIND=Rkind)                      :: err_max
   INTEGER                               :: i
 
   ! Allocation
   ALLOCATE(M(n,n), Iden(n,n))
 
   ! M = U * U^T
   M = MATMUL(U, CONJG(TRANSPOSE(U)))
 
   ! Identité
   Iden = CZERO
   DO i = 1, n
     Iden(i,i) = CONE
   END DO
 
   ! Erreur max (norme infinie, ici max |M - I|)
   err_max = MAXVAL(ABS(M - Iden))
   print *,"err_max :" , err_max
 
   DEALLOCATE(M, Iden)
 END SUBROUTINE check_unitarity


!----------------------------------------------------------------
! Complex dot product <a|b> = sum(conjg(a)*b)  (returns complex)
!----------------------------------------------------------------
 SUBROUTINE Dotproduct(W1, W2,result)
  IMPLICIT NONE
  COMPLEX(KIND=Rkind), INTENT(IN)       :: W1(:) ,W2(:)
  REAL(KIND=Rkind), INTENT(INOUT)       :: result

  INTEGER                               :: i,nb

  nb =SIZE(W1)
  result  = ZERO 
  result = SUM(CONJG(W1) * W2)
 !print*, "res :", result
END SUBROUTINE Dotproduct



!--------------------------------------------------------------------------------
!> SUBROUTINE: Lanczos_step (with external Gram-Schmidt Orthogonalization)
!> Description: Executes a single step of the Lanczos algorithm. It computes
!>              the new vector candidate and orthogonalizes it against the
!>              existing Krylov basis.
!> Arguments:
!>   - H      [in]      (Op_t) : The Hamiltonian operator.
!>   - Q      [in]      (complex(kind=Rkind), allocatable) : The current orthonormal Krylov basis.
!>   - w_new  [inout]   (complex(kind=Rkind), allocatable) : The new, un-normalized Lanczos vector candidate.
!>   - alpha  [out]     (real(kind=Rkind)) : Diagonal element of the tridiagonal matrix.
!>   - k      [in]      (integer) : Current size of the Krylov basis Q.
!>   - Basis  [in]      (Basis_t) : Basis of the wavefunction.
!--------------------------------------------------------------------------------
    SUBROUTINE Lanczos_step(H, Q, w_new, alpha, beta, k, Basis)
      IMPLICIT NONE
      TYPE(Op_t),     INTENT(IN)         :: H
      COMPLEX(KIND=Rkind), INTENT(IN)    :: Q(:,:)
      COMPLEX(KIND=Rkind), INTENT(OUT)   :: w_new(:)
      REAL(KIND=Rkind), INTENT(INOUT)    :: alpha(:), beta(:)
      INTEGER, INTENT(IN)                :: k
      TYPE(Basis_t), INTENT(IN)          :: Basis
    
      ! Local variables
      TYPE(psi_t)       :: Hpsi, q_k_psi
      REAL(KIND=Rkind)  :: norm_w
    
      ! Initialisation des vecteurs de travail
      CALL init_psi(Hpsi,   Basis, cplx=.TRUE., grid=.FALSE.)
      CALL init_psi(q_k_psi,Basis, cplx=.TRUE., grid=.FALSE.)
    
      ! q_k = colonne k de Q
      q_k_psi%CVec(:) = Q(:,k)
    
      ! H|q_k>
      CALL calc_OpPsi(H, q_k_psi, Hpsi)
    
      ! α_k = <q_k|H|q_k> (réel si H hermitien)
      CALL Dotproduct(Q(:,k),Hpsi%CVec(:),alpha(k))
    
      ! w = H|q_k> - α_k |q_k>
      w_new(:) = Hpsi%CVec(:) - alpha(k) * Q(:,k)
      If (k > 1) THEN
        w_new(:) = w_new(:) - beta(k-1)*Q(:,k-1)
      END IF
    
      ! Orthogonalisation par Gram–Schmidt (contre q_1...q_k)
      CALL gram_schmidt(Q, w_new, k)   
      ! β_{k+1} = ||w||
       CALL Dotproduct(w_new,w_new,norm_w)
      norm_w =  SQRT(norm_w)
      beta(k+1)  = norm_w

      !WRITE(out_unit, *) 'alpha and beta : ', alpha, beta
    
      ! Normalisation : q_{k+1} = w / β
      IF (norm_w > ZERO) THEN
        w_new(:) = w_new(:) / norm_w
      END IF


    
      CALL dealloc_psi(Hpsi)
      CALL dealloc_psi(q_k_psi)
    
    END SUBROUTINE Lanczos_step  
    
    

  !----------------------------------------------------------------
  !> @brief Estimate the norm of the wavefunction difference
  !> @details Estimate the norm of the wavefunction difference
  !>          Δψ using the asymptotic error formula:
  !>
  !>          ||Δψ|| ≈ (dt)^(n+1) * (β1 * β2 * ... * βn) / (n+1)!
  !>
  !> @param dt       - REAL, time step
  !> @param betas    - REAL array of β coefficients (β1, β2, ..., βn)
  !> @param n        - INTEGER, order of the method (length of betas)
  !> @param debug_in - LOGICAL, optional: if .TRUE., prints debug info
  !>
  !> @return norm_estimate - REAL, estimated norm of Δψ
  !
  !  Inputs   :
  !     dt       - REAL, time step
  !     betas    - REAL array of β coefficients (β1, β2, ..., βn)
  !     n        - INTEGER, order of the method (length of betas)
  !     debug_in - LOGICAL, optional: if .TRUE., prints debug info
  !
  !  Output   :
  !     norm_estimate - REAL, estimated norm of Δψ
  !
  !----------------------------------------------------------------
FUNCTION EstimateDeltaPsiNorm(dt, betas, n) RESULT(norm_estimate)
  IMPLICIT NONE
  ! Arguments
  REAL(KIND=Rkind), INTENT(IN)           :: dt
  REAL(KIND=Rkind), INTENT(IN)           :: betas(:)
  INTEGER, INTENT(IN)                    :: n
  REAL(KIND=Rkind)                       :: norm_estimate

  ! Local variables
  REAL(KIND=Rkind) :: prod_betas
  LOGICAL         :: debug
  INTEGER :: i

  ! Initialize debug flag
  debug = .TRUE.


  !-----------------------------------------------
  ! Compute product of betas
  !-----------------------------------------------
  prod_betas = ONE
  DO i = 1, n
     prod_betas = prod_betas * betas(i)
  END DO


  !-----------------------------------------------
  ! Final estimate
  !-----------------------------------------------
  norm_estimate = (dt**(n)) * prod_betas / Factorial(n)

  !-----------------------------------------------
  ! Debug output
  !-----------------------------------------------
  IF (debug) THEN
     WRITE(out_unit,*) "-----------------------------------"
     WRITE(out_unit,*) " Debug info: EstimateDeltaPsiNorm"
     WRITE(out_unit,*) " dt         =", dt
     WRITE(out_unit,*) " n          =", n
     WRITE(out_unit,*) " betas      =", betas
     WRITE(out_unit,*) " prod_betas =", prod_betas
     WRITE(out_unit,*) " result     =", norm_estimate
     WRITE(out_unit,*) "-----------------------------------"
  END IF

END FUNCTION EstimateDeltaPsiNorm

    
    !--------------------------------------------------------------------------------
    !> SUBROUTINE: march_mSIL (Updated for Gram-Schmidt)
    !> Description:  This is the core routine for time propagation using the
    !>              Short Iterative Lanczos (SIL) method. It orchestrates the
    !>              steps of building the Krylov basis, diagonalizing the projected
    !>              tridiagonal matrix, and calculating the evolution operator to
    !>              propagate the wavefunction for a single time step.
    !> Arguments:
    !>   psi         [inout] (psi_t) : Wavefunction to propagate.
    !>   H           [in]    (Op_t) : Hamiltonian operator.
    !>   tol   [in]    (real(kind=Rkind)) : Precision required for convergence.
    !>   deltat_t    [in]    (real(kind=Rkind)) : Time step for the propagation.
    !>   max_krylov  [in]    (integer) : Maximum dimension of the Krylov space.
    !--------------------------------------------------------------------------------
    SUBROUTINE march_mSIL(H, psi, Basis, dt, max_krylov, tol)
      IMPLICIT NONE
      TYPE(Op_t),     INTENT(IN)           :: H
      TYPE(psi_t),    INTENT(INOUT)        :: psi
      TYPE(Basis_t),  INTENT(IN)           :: Basis
      REAL(KIND=Rkind), INTENT(IN)         :: dt
      INTEGER, INTENT(IN)                  :: max_krylov
      REAL(KIND=Rkind), INTENT(IN)         :: tol
    
      ! --- Variables locales ---
      COMPLEX(KIND=Rkind), ALLOCATABLE     :: K_basis(:,:), w_new(:)
      REAL(KIND=Rkind),    ALLOCATABLE     :: alpha(:), beta(:), eigvals(:)
      REAL(KIND=Rkind),    ALLOCATABLE     :: eigvecs(:,:)
      COMPLEX(KIND=Rkind), ALLOCATABLE     :: U_prop(:,:), cvec_tmp(:)
      INTEGER                              :: n, iter, nb_k, info
      REAL(KIND=Rkind)                     :: norm_psi,estimate_error
      LOGICAL, PARAMETER                   :: debug=.TRUE.
    
      n = SIZE(psi%CVec)
    
      ! Allocation
      ALLOCATE(K_basis(n,max_krylov), w_new(n))
      ALLOCATE(alpha(max_krylov), beta(max_krylov))
    
      ! --- Normalisation initiale ---
      call Calc_Norm_OF_Psi(psi, norm_psi)
      IF (norm_psi > ZERO) THEN
        K_basis(:,1) = psi%CVec / norm_psi
      ELSE
        WRITE(out_unit,*) "Error: initial psi is zero vector"
        STOP "march_mSIL error"
      END IF
    
      ! --- Lanczos iterations ---
      DO iter = 1, max_krylov
        CALL Lanczos_step(H, K_basis, w_new, alpha, beta, iter, Basis)
    
        ! Critère d’arrêt : si la norme de w_new est trop petite
         CALL Dotproduct(w_new(:),w_new(:),estimate_error)
         estimate_error = SQRT(estimate_error)*(dt**iter)/Factorial(iter)
        IF (iter < max_krylov .AND. estimate_error > tol) THEN
          K_basis(:, iter+1) = w_new
        ELSE 
            nb_k = iter
            WRITE(out_unit, *) 'Converged after', nb_k, 'iterations.'
            EXIT
        END IF

        IF (iter == max_krylov) THEN
          nb_k = max_krylov   ! si pas de convergence avant la limite
          WRITE(out_unit,*) "Warning: max_krylov :",max_krylov," reached"
          EXIT
        END IF
      END DO
    
      ALLOCATE(eigvals(nb_k), eigvecs(nb_k,nb_k))
      ALLOCATE(U_prop(nb_k,nb_k), cvec_tmp(nb_k))
    
      ! --- Construction et diagonalisation de T ---
      CALL build_tridiag(alpha, beta, nb_k, eigvals, eigvecs, info)
      IF (info /= 0) THEN
        WRITE(out_unit,*) "Error in build_tridiag, info=", info
        STOP "march_mSIL error"
      END IF
    
      ! --- Construction du propagateur réduit ---
      CALL build_propagator(eigvals, eigvecs, nb_k, dt, U_prop)

      IF(debug) CALL check_unitarity(U_prop, nb_k)
    
      ! --- Projection de psi dans la base de Krylov ---
      cvec_tmp(1:nb_k) = MATMUL(CONJG(TRANSPOSE(K_basis(:,1:nb_k))), psi%CVec(:))
    
      ! --- Propagation ---
      cvec_tmp(1:nb_k) = MATMUL(U_prop(1:nb_k,1:nb_k), cvec_tmp(1:nb_k))
    
      ! --- Reconstruction ---
      psi%CVec(:) = MATMUL(K_basis(:,1:nb_k), cvec_tmp(1:nb_k))
    
      ! --- Nettoyage ---
      DEALLOCATE(K_basis, w_new, alpha, beta, eigvals, eigvecs, U_prop, cvec_tmp)
    
    END SUBROUTINE march_mSIL


   !--------------------------------------------------------------------------------
   !> SUBROUTINE: find_ground_state_SIL
   !> Description: This routine computes the ground state energy and corresponding
   !>              wave function of a quantum system using the Short Iterative
   !>              Lanczos (SIL) method. It builds a Krylov subspace from an
   !>              initial guess vector, represents the Hamiltonian in this
   !>              subspace as a small tridiagonal matrix, and diagonalizes it to
   !>              find the lowest eigenvalue and eigenvector. The ground state
   !>              wave function is then reconstructed in the original basis.
   !> Arguments:
   !>   H                      [in]    (Op_t) : Hamiltonian operator.
   !>   psi                    [inout] (psi_t) : Initial guess wave function.
   !>   Basis                  [in]    (Basis_t) : Basis information.
   !>   max_krylov             [in]    (integer) : Maximum dimension of the Krylov space.
   !>   tol                    [in]    (real(kind=Rkind)) : Tolerance for convergence.
   !>   Gs_Energy              [out]   (real(kind=Rkind)) : Ground state energy.
   !>   Gs_Psi                 [out]   (psi_t) : Ground state wave function.
   !--------------------------------------------------------------------------------
   SUBROUTINE find_ground_state_SIL(Gs_Energy, Gs_Psi, H, psi, Basis, max_krylov, tol)
     IMPLICIT NONE
     TYPE(Op_t),     INTENT(IN)              :: H
     TYPE(psi_t),    INTENT(IN)              :: psi
     TYPE(Basis_t),  INTENT(IN)              :: Basis
     INTEGER, INTENT(IN)                     :: max_krylov
     REAL(KIND=Rkind), INTENT(INOUT)         :: Gs_Energy
     TYPE(psi_t),    INTENT(INOUT)           :: Gs_Psi
     REAL(KIND=Rkind), INTENT(IN)            :: tol
   
     ! --- Variables locales -----------------------------------------
     COMPLEX(KIND=Rkind), ALLOCATABLE        :: K_basis(:,:), w_new(:)
     REAL(KIND=Rkind),    ALLOCATABLE        :: alpha(:), beta(:), eigvals(:)
     REAL(KIND=Rkind),    ALLOCATABLE        :: eigvecs(:,:)
     COMPLEX(KIND=Rkind), ALLOCATABLE        :: phi_0(:)
     INTEGER                                 :: n, iter, nb_k, info
     REAL(KIND=Rkind)                        :: norm_psi
     REAL(KIND=Rkind)                        :: min_eigval
     INTEGER                                 :: min_idx
     REAL(KIND=Rkind)                        :: norm_w_new
     LOGICAL, PARAMETER                      :: debug=.TRUE.
   
     n = SIZE(psi%CVec)
   
     IF (debug) THEN
       WRITE(out_unit,*) "--- Starting find_ground_state_SIL ---"
       WRITE(out_unit,*) "max_krylov: ", max_krylov
       WRITE(out_unit,*) "tolerance: ", tol
     END IF
   
     ! Allocation ------------------------------------------------
     ALLOCATE(K_basis(n,max_krylov), w_new(n))
     ALLOCATE(alpha(max_krylov), beta(max_krylov))
   
     ! --- Normalisation initiale ---
     CALL Calc_Norm_OF_Psi(psi, norm_psi)
     IF (norm_psi > ZERO) THEN
       K_basis(:,1) = psi%CVec / norm_psi
       IF (debug) WRITE(out_unit,*) "Initial psi normalized. Norm = ", norm_psi
     ELSE
       WRITE(out_unit,*) "Error: initial psi is zero vector"
       STOP "find_ground_state_SIL error"
     END IF
   
     ! --- Lanczos iterations ---
     DO iter = 1, max_krylov
       CALL Lanczos_step(H, K_basis, w_new, alpha, beta, iter, Basis)
   
       CALL Dotproduct(w_new, w_new, norm_w_new)
       norm_w_new = SQRT(norm_w_new)/Factorial(iter)
       
       IF (debug) WRITE(out_unit,*) "Iteration ", iter, ": norm(w_new) = ", norm_w_new
   
       ! Stop condition for convergence
       IF (norm_w_new < tol) THEN
         nb_k = iter
         IF (debug) WRITE(out_unit, *) 'Converged after', nb_k, 'iterations.'
         EXIT
       END IF
   
       IF (iter < max_krylov) THEN
         K_basis(:, iter+1) = w_new
       END IF
   
       IF (iter == max_krylov) THEN
         nb_k = max_krylov   ! if no convergence before the limit
         IF (debug) WRITE(out_unit,*) "Warning: max_krylov :",max_krylov," reached, stopping."
         EXIT
       END IF
     END DO
   
     ALLOCATE(eigvals(nb_k), eigvecs(nb_k,nb_k))
   
     ! --- Construction and diagonalisation of T ---
     IF (debug) WRITE(out_unit,*) "Building and diagonalizing tridiagonal matrix of size ", nb_k
     CALL build_tridiag(alpha, beta, nb_k, eigvals, eigvecs, info)
     IF (info /= 0) THEN
       WRITE(out_unit,*) "Error in build_tridiag, info=", info
       STOP "find_ground_state_SIL error"
     END IF
   
     ! --- Extraction of the ground state ---
     min_eigval = eigvals(1)
     min_idx = 1
     DO iter = 2, nb_k
       IF (eigvals(iter) < min_eigval) THEN
         min_eigval = eigvals(iter)
         min_idx = iter
       END IF
     END DO
   
     Gs_Energy = min_eigval
     IF (debug) WRITE(out_unit, '(A,F20.15)') "Ground state energy found: ", Gs_Energy
   
     ALLOCATE(phi_0(nb_k))
     phi_0(:) = eigvecs(:, min_idx)
   
     ! --- Reconstruction of the ground state vector ---
     Gs_Psi%CVec(:) = MATMUL(K_basis(:,1:nb_k), phi_0(:))
   
     IF (debug) WRITE(out_unit,*) "--- Exiting find_ground_state_SIL ---"
   
     ! --- Nettoyage ---
     DEALLOCATE(K_basis, w_new, alpha, beta, eigvals, eigvecs, phi_0)
   
   END SUBROUTINE find_ground_state_SIL

end module lanczos_m
