!> @file TISE_m.f90
!> @brief Résolution de l'équation de Schrödinger stationnaire.
!>
!> Ce module contient exclusivement les opérations TISE : construction
!> optionnelle de H%CMat, solveurs itératifs, relaxation en temps imaginaire,
!> écriture des états propres et écriture des états de transition.
!> L'opération H|psi> est fournie par calc_OpPsi de Op_m. Les méthodes
!> itératives ne construisent donc pas la matrice complète du Hamiltonien.
!>
!> Méthodes : dense, lanczos, davidson bloc préconditionné, arpack.
MODULE TISE_m
  USE QDUtil_m
  USE Op_m, ONLY: Op_t, calc_OpPsi
  USE psi_m, ONLY: psi_t, init_psi, dealloc_psi, Calc_Norm_OF_Psi
  USE Basis_m, ONLY: Get_TotalBasisCoefficientCount, Basis_UsesSmolyakSparseGrid, Debug_Enabled
  USE SmolyakSparseBasisAlgebra_m, ONLY: DotSmolyakBasis
  USE mod_Constant, ONLY: get_Conv_au_TO_unit
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Solve_TISE, Relax_TISE, Build_Hamiltonian_Matrix

#if __ARPACK == 1

  INTERFACE
    SUBROUTINE dsaupd(ido,bmat,n,which,nev,tol,resid,ncv,v,&
      ldv,iparam,ipntr,workd,workl,lworkl,info)

      USE QDUtil_m, ONLY: Rkind
      INTEGER          :: ido,n,nev,ncv,ldv,lworkl,info
      CHARACTER(LEN=1) :: bmat
      CHARACTER(LEN=2) :: which
      REAL(KIND=Rkind) :: tol,resid(*),v(ldv,*),workd(*),workl(*)
      INTEGER          :: iparam(*),ipntr(*)

    END SUBROUTINE dsaupd

    SUBROUTINE dseupd(rvec,howmny,select,d,z,ldz,sigma,bmat,n,which,nev,&
      tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info)

      USE QDUtil_m, ONLY: Rkind
      LOGICAL          :: rvec,select(*)
      CHARACTER(LEN=1) :: howmny,bmat
      CHARACTER(LEN=2) :: which
      INTEGER          :: ldz,n,nev,ncv,ldv,lworkl,info
      REAL(KIND=Rkind) :: d(ncv,2),z(ldz,*),sigma,tol,resid(*),v(ldv,*),workd(*),workl(*)
      INTEGER          :: iparam(*),ipntr(*)
    END SUBROUTINE dseupd
  END INTERFACE

#endif
CONTAINS

  !> @brief Pilote un calcul TISE et écrit tous les résultats.
  !> @param[inout] H Hamiltonien matrix-free.
  !> @param[in] tise_method dense, lanczos, davidson ou arpack.
  !> @param[in] tise_nstates nombre d'états propres demandés.
  !> @param[in] tise_max_iter nombre maximal d'itérations.
  !> @param[in] tise_tol tolérance sur les résidus.
  !> @param[in] energy_unit au ou cm-1.
  !> @param[in] tise_energy_file fichier des énergies.
  !> @param[in] tise_transition_file fichier des transitions.
  SUBROUTINE Solve_TISE(H,tise_method,tise_nstates,tise_max_iter,tise_tol,energy_unit, &
                        tise_energy_file,tise_transition_file,psi_initial,debug)

    TYPE(Op_t), INTENT(INOUT)            :: H
    CHARACTER(*), INTENT(IN)             :: tise_method,energy_unit,tise_energy_file,tise_transition_file
    INTEGER, INTENT(IN)                  :: tise_nstates,tise_max_iter
    REAL(KIND=Rkind), INTENT(IN)         :: tise_tol
    TYPE(psi_t), INTENT(INOUT), OPTIONAL :: psi_initial
    LOGICAL, INTENT(IN), OPTIONAL        :: debug
    CHARACTER(LEN=32)                    :: method,unit
    REAL(KIND=Rkind), ALLOCATABLE        :: tise_energies(:)
    COMPLEX(KIND=Rkind), ALLOCATABLE     :: tise_eigenvectors(:,:)
    INTEGER                              :: dimension_hamiltonian,number_states,info
    LOGICAL                              :: debug_tise

    debug_tise=Debug_Enabled()

    IF(PRESENT(debug)) debug_tise=debug

    dimension_hamiltonian=Get_TotalBasisCoefficientCount(H%Basis)
    number_states=MIN(MAX(1,tise_nstates),dimension_hamiltonian)

    IF(tise_max_iter<number_states) ERROR STOP 'Solve_TISE: tise_max_iter must be >= tise_nstates.'
    IF(tise_tol<=ZERO) ERROR STOP 'Solve_TISE: tise_tol must be positive.'

    method=ADJUSTL(tise_method)
    unit=ADJUSTL(energy_unit)

    CALL string_uppercase_TO_lowercase(method)
    CALL string_uppercase_TO_lowercase(unit)

    IF(TRIM(unit)=='cm1') unit='cm-1'
    IF(TRIM(unit)/='au'.AND.TRIM(unit)/='cm-1') ERROR STOP 'Solve_TISE: energy_unit must be au or cm-1.'

    ALLOCATE(tise_energies(number_states),tise_eigenvectors(dimension_hamiltonian,number_states))
    info=0

    IF(debug_tise) WRITE(out_unit,*) 'TISE input: ',TRIM(method),number_states,tise_max_iter,tise_tol

    SELECT CASE(TRIM(method))
    CASE('dense','lapack')
      CALL solve_tise_dense(H,number_states,tise_energies,tise_eigenvectors,debug_tise)
    CASE('lanczos')
      CALL solve_tise_lanczos(H,number_states,tise_max_iter,tise_tol,&
        tise_energies,tise_eigenvectors,info,debug_tise)
    CASE('davidson')
      CALL solve_tise_davidson_block(H,number_states,tise_max_iter,tise_tol,&
        tise_energies,tise_eigenvectors,info,debug_tise)
    CASE('arpack')
      CALL solve_tise_arpack(H,number_states,tise_max_iter,tise_tol,&
        tise_energies,tise_eigenvectors,info,debug_tise)
    CASE DEFAULT
      ERROR STOP 'Solve_TISE: tise_method must be dense, lanczos, davidson or arpack.'
    END SELECT
    CALL write_tise_results(tise_energies,tise_eigenvectors,unit,&
      tise_energy_file,tise_transition_file,psi_initial,info,debug_tise)
    IF(PRESENT(psi_initial)) psi_initial%CVec=tise_eigenvectors(:,1)

    DEALLOCATE(tise_energies,tise_eigenvectors)

  END SUBROUTINE Solve_TISE

  !> @brief Construit la matrice dense du Hamiltonien dans le module TISE.
  !> @details Chaque colonne est obtenue par application de H à un vecteur de base.
  SUBROUTINE Build_Hamiltonian_Matrix(H,debug)
    TYPE(Op_t), INTENT(INOUT)     :: H
    LOGICAL, INTENT(IN), OPTIONAL :: debug
    TYPE(psi_t)                   :: psi_basis,hpsi_basis
    INTEGER                       :: dimension_hamiltonian,index_column,index_row
    LOGICAL                       :: debug_tise

    debug_tise=Debug_Enabled()
    IF(PRESENT(debug)) debug_tise=debug
    dimension_hamiltonian=Get_TotalBasisCoefficientCount(H%Basis)

    IF(.NOT.ALLOCATED(H%CMat)) ALLOCATE(H%CMat(dimension_hamiltonian,dimension_hamiltonian))
    H%CMat=CZERO

    CALL init_psi(psi_basis,H%Basis,cplx=.TRUE.,grid=.FALSE.)
    CALL init_psi(hpsi_basis,H%Basis,cplx=.TRUE.,grid=.FALSE.)

    DO index_column=1,dimension_hamiltonian
      psi_basis%CVec=CZERO
      psi_basis%CVec(index_column)=CONE
      CALL calc_OpPsi(H,psi_basis,hpsi_basis)
      DO index_row=1,dimension_hamiltonian
        IF(Basis_UsesSmolyakSparseGrid(H%Basis)) THEN
          psi_basis%CVec=CZERO
          psi_basis%CVec(index_row)=CONE
          H%CMat(index_row,index_column)=&
            DotSmolyakBasis(psi_basis%CVec,hpsi_basis%CVec,H%Basis%SmolyakSparseGrid)
        ELSE
          ! La base coefficientielle est orthonormée : la ligne i de la
          ! colonne j est directement le coefficient i de H|psi_j>.
          ! Cette écriture évite de reconstruire inutilement un produit
          ! scalaire pour chaque élément de la matrice.
          H%CMat(index_row,index_column)=hpsi_basis%CVec(index_row)
        END IF
      END DO
    END DO
    IF(debug_tise) CALL Write_VecMat(H%CMat,out_unit,MIN(8,dimension_hamiltonian),info='TISE Hamiltonian matrix')

    CALL dealloc_psi(psi_basis)
    CALL dealloc_psi(hpsi_basis)

  END SUBROUTINE Build_Hamiltonian_Matrix

  !> @brief Applique H à un vecteur de coefficients.
  !> @details Utilise les objets psi du projet et calc_OpPsi, sans nouvelle
  !>          routine de produit matriciel.
  SUBROUTINE Apply_Hamiltonian(H,psi_coefficients,hpsi_coefficients,psi_work,hpsi_work,debug)
    TYPE(Op_t), INTENT(IN)            :: H
    COMPLEX(KIND=Rkind), INTENT(IN)   :: psi_coefficients(:)
    COMPLEX(KIND=Rkind), INTENT(OUT)  :: hpsi_coefficients(:)
    TYPE(psi_t), INTENT(INOUT)        :: psi_work,hpsi_work
    LOGICAL, INTENT(IN), OPTIONAL     :: debug
    LOGICAL                           :: debug_tise

    debug_tise=Debug_Enabled()
    IF(PRESENT(debug)) debug_tise=debug
    psi_work%CVec=psi_coefficients
    CALL calc_OpPsi(H,psi_work,hpsi_work)
    hpsi_coefficients=hpsi_work%CVec

    IF(debug_tise) WRITE(out_unit,*) 'Apply_Hamiltonian norm = ',&
      SQRT(MAX(ZERO,REAL(DOT_PRODUCT(hpsi_coefficients,hpsi_coefficients),KIND=Rkind)))

  END SUBROUTINE Apply_Hamiltonian

  !> @brief Résout le petit problème dense avec QDUtil/LAPACK.
  SUBROUTINE solve_tise_dense(H,number_states,tise_energies,tise_eigenvectors,debug)

    TYPE(Op_t), INTENT(INOUT)         :: H
    INTEGER, INTENT(IN)               :: number_states
    REAL(KIND=Rkind), INTENT(OUT)     :: tise_energies(:)
    COMPLEX(KIND=Rkind), INTENT(OUT)  :: tise_eigenvectors(:,:)
    LOGICAL, INTENT(IN)               :: debug
    INTEGER :: dimension_hamiltonian
    REAL(KIND=Rkind), ALLOCATABLE     :: all_energies(:)
    COMPLEX(KIND=Rkind), ALLOCATABLE  :: all_vectors(:,:)

    dimension_hamiltonian=SIZE(tise_eigenvectors,1)

    CALL Build_Hamiltonian_Matrix(H,debug)
    ALLOCATE(all_energies(dimension_hamiltonian),all_vectors(dimension_hamiltonian,dimension_hamiltonian))
    CALL diagonalization(H%CMat,all_energies,all_vectors,sort=1)

    tise_energies=all_energies(1:number_states)
    tise_eigenvectors=all_vectors(:,1:number_states)

    IF(debug) WRITE(out_unit,*) 'solve_tise_dense ground energy = ',tise_energies(1)
    DEALLOCATE(all_energies,all_vectors)

  END SUBROUTINE solve_tise_dense

  !> @brief Implémente le Lanczos et contrôle les résidus Ritz.
  SUBROUTINE solve_tise_lanczos(H,number_states,max_iter,tolerance,tise_energies,tise_eigenvectors,info,debug)
    TYPE(Op_t), INTENT(IN)             :: H
    INTEGER, INTENT(IN)                :: number_states,max_iter
    REAL(KIND=Rkind), INTENT(IN)       :: tolerance
    REAL(KIND=Rkind), INTENT(OUT)      :: tise_energies(:)
    COMPLEX(KIND=Rkind), INTENT(OUT)   :: tise_eigenvectors(:,:)
    INTEGER, INTENT(OUT)               :: info
    LOGICAL, INTENT(IN)                :: debug
    INTEGER                            :: dimension_hamiltonian,krylov_size
    INTEGER                            :: max_krylov,index_basis,index_state
    COMPLEX(KIND=Rkind), ALLOCATABLE   :: krylov_vectors(:,:),candidate(:),hpsi(:)
    REAL(KIND=Rkind), ALLOCATABLE      :: alpha(:),beta(:),projected_energies(:)
    REAL(KIND=Rkind), ALLOCATABLE      :: projected_vectors(:,:)
    REAL(KIND=Rkind), ALLOCATABLE      :: residual_norms(:)
    TYPE(psi_t)                        :: psi_work,hpsi_work
    REAL(KIND=Rkind)                   :: candidate_norm

    dimension_hamiltonian=SIZE(tise_eigenvectors,1)
    max_krylov=MIN(max_iter,dimension_hamiltonian)

    ALLOCATE(krylov_vectors(dimension_hamiltonian,max_krylov))
    ALLOCATE(candidate(dimension_hamiltonian),hpsi(dimension_hamiltonian))
    ALLOCATE(alpha(max_krylov),beta(max_krylov+1),projected_energies(max_krylov))
    ALLOCATE(projected_vectors(max_krylov,max_krylov),residual_norms(number_states))

    CALL init_psi(psi_work,H%Basis,cplx=.TRUE.,grid=.FALSE.) 
    CALL init_psi(hpsi_work,H%Basis,cplx=.TRUE.,grid=.FALSE.)

    krylov_vectors=CZERO
    alpha=ZERO
    beta=ZERO
    krylov_vectors(1,1)=CONE
    krylov_size=1
    info=1

    DO WHILE(krylov_size<=max_krylov)
      CALL Apply_Hamiltonian(H,krylov_vectors(:,krylov_size),hpsi,psi_work,hpsi_work,.FALSE.)
      alpha(krylov_size)=REAL(DOT_PRODUCT(krylov_vectors(:,krylov_size),hpsi),KIND=Rkind)
      candidate=hpsi-alpha(krylov_size)*krylov_vectors(:,krylov_size)
      IF(krylov_size>1) candidate=candidate-beta(krylov_size-1)*krylov_vectors(:,krylov_size-1)
      DO index_basis=1,krylov_size
         candidate=candidate-DOT_PRODUCT(krylov_vectors(:,index_basis),&
          candidate)*krylov_vectors(:,index_basis)
      END DO

      candidate_norm=SQRT(MAX(ZERO,REAL(DOT_PRODUCT(candidate,candidate),KIND=Rkind)))
      CALL diagonalize_tridiagonal(alpha,beta,krylov_size,projected_energies,projected_vectors)

      IF(krylov_size>=number_states) THEN
        tise_eigenvectors=MATMUL(krylov_vectors(:,1:krylov_size),&
          CMPLX(projected_vectors(1:krylov_size,1:number_states),ZERO,KIND=Rkind)) 
        tise_energies=projected_energies(1:number_states)
        DO index_state=1,number_states
          CALL Apply_Hamiltonian(H,tise_eigenvectors(:,index_state),hpsi,psi_work,hpsi_work,.FALSE.) 
          residual_norms(index_state)=MAXVAL(ABS(hpsi-tise_energies(index_state)*tise_eigenvectors(:,index_state)))
        END DO

        IF(MAXVAL(residual_norms)<=tolerance) THEN
          info=0
          EXIT
        END IF
      END IF

      IF(candidate_norm<=SQRT(tolerance).OR.krylov_size==max_krylov) EXIT
      krylov_vectors(:,krylov_size+1)=candidate/candidate_norm
      beta(krylov_size)=candidate_norm
      krylov_size=krylov_size+1
      
    END DO

    IF(debug) WRITE(out_unit,*) 'solve_tise_lanczos subspace/residual = ',krylov_size,MAXVAL(residual_norms)

    CALL dealloc_psi(psi_work)
    CALL dealloc_psi(hpsi_work)

    DEALLOCATE(krylov_vectors,candidate,hpsi,alpha,beta)
    DEALLOCATE(projected_energies,projected_vectors,residual_norms)

  END SUBROUTINE solve_tise_lanczos

  !> @brief Davidson bloc préconditionné diagonal.
  !> @details Construit un bloc initial, diagonalise V^H H V, forme les
  !>          résidus, applique (diag(H)-theta)^-1, orthogonalise et agrandit
  !>          le sous-espace. Lorsque le sous-espace est plein, les vecteurs
  !>          de Ritz sont conservés et un nouveau bloc de corrections est créé.
  SUBROUTINE solve_tise_davidson_block(H,number_states,max_iter,tolerance,tise_energies,tise_eigenvectors,info,debug)

    TYPE(Op_t), INTENT(IN)              :: H
    INTEGER, INTENT(IN)                 :: number_states,max_iter
    REAL(KIND=Rkind), INTENT(IN)        :: tolerance
    REAL(KIND=Rkind), INTENT(OUT)       :: tise_energies(:)
    COMPLEX(KIND=Rkind), INTENT(OUT)    :: tise_eigenvectors(:,:)
    INTEGER, INTENT(OUT)                :: info
    LOGICAL, INTENT(IN)                 :: debug
    INTEGER                             :: dimension_hamiltonian,max_subspace,subspace_size
    INTEGER                             :: index_state,index_vector,index_iteration,number_corrections
    COMPLEX(KIND=Rkind), ALLOCATABLE    :: subspace_vectors(:,:),subspace_hvectors(:,:)
    COMPLEX(KIND=Rkind), ALLOCATABLE    :: correction(:),ritz(:,:),ritz_h(:,:),residual(:,:)
    REAL(KIND=Rkind), ALLOCATABLE       :: diagonal_hamiltonian(:),small_energies(:)
    REAL(KIND=Rkind), ALLOCATABLE       :: small_vectors(:,:),residual_norms(:)
    TYPE(psi_t)                         :: psi_work,hpsi_work
    REAL(KIND=Rkind)                    :: correction_norm,denominator

    dimension_hamiltonian=SIZE(tise_eigenvectors,1)
    ! La dimension du sous-espace suit la limite d'itérations demandée.
    ! Cela évite d'arrêter prématurément le Davidson lorsque le premier
    ! sous-espace est plein avant convergence.
    max_subspace=MIN(dimension_hamiltonian,MAX(number_states+4,max_iter))

    ALLOCATE(subspace_vectors(dimension_hamiltonian,max_subspace))
    ALLOCATE(subspace_hvectors(dimension_hamiltonian,max_subspace))
    ALLOCATE(diagonal_hamiltonian(dimension_hamiltonian))
    ALLOCATE(correction(dimension_hamiltonian))
    ALLOCATE(ritz(dimension_hamiltonian,number_states))
    ALLOCATE(ritz_h(dimension_hamiltonian,number_states))
    ALLOCATE(residual(dimension_hamiltonian,number_states))
    ALLOCATE(small_energies(max_subspace))
    ALLOCATE(small_vectors(max_subspace,max_subspace))
    ALLOCATE(residual_norms(number_states))
     
    CALL init_psi(psi_work,H%Basis,cplx=.TRUE.,grid=.FALSE.)
    CALL init_psi(hpsi_work,H%Basis,cplx=.TRUE.,grid=.FALSE.)
    CALL compute_hamiltonian_diagonal(H,diagonal_hamiltonian,psi_work,hpsi_work,debug)

    subspace_vectors=CZERO
    subspace_hvectors=CZERO
    subspace_size=MIN(number_states,max_subspace)

    DO index_vector=1,subspace_size
      subspace_vectors(index_vector,index_vector)=CONE
    END DO

    CALL orthonormalize_block(subspace_vectors,subspace_size,debug) 
    CALL apply_hamiltonian_block(H,subspace_vectors,subspace_hvectors,subspace_size,psi_work,hpsi_work,debug)

    info=1

    DO index_iteration=1,max_iter
      CALL diagonalize_projected_problem(subspace_vectors,subspace_hvectors,subspace_size,small_energies,small_vectors,debug)

      ritz=MATMUL(subspace_vectors(:,1:subspace_size),CMPLX(small_vectors(1:subspace_size,1:number_states),ZERO,KIND=Rkind))
      ritz_h=MATMUL(subspace_hvectors(:,1:subspace_size),CMPLX(small_vectors(1:subspace_size,1:number_states),ZERO,KIND=Rkind))

      tise_energies=small_energies(1:number_states)
      tise_eigenvectors=ritz

      DO index_state=1,number_states
        residual(:,index_state)=ritz_h(:,index_state)-tise_energies(index_state)*ritz(:,index_state)
        residual_norms(index_state)=MAXVAL(ABS(residual(:,index_state)))
      END DO

      IF(MAXVAL(residual_norms)<=tolerance) THEN
        info=0
        EXIT
      END IF
      
      number_corrections=0

      DO index_state=1,number_states
        correction=residual(:,index_state)
        DO index_vector=1,dimension_hamiltonian
          denominator=diagonal_hamiltonian(index_vector)-tise_energies(index_state)
          IF(ABS(denominator)>SQRT(EPSILON(ONE))) correction(index_vector)=correction(index_vector)/denominator
        END DO
        CALL orthonormalize_against_block(correction,subspace_vectors,subspace_size,debug)
        correction_norm=SQRT(MAX(ZERO,REAL(DOT_PRODUCT(correction,correction),KIND=Rkind)))
        IF(correction_norm>SQRT(tolerance).AND.subspace_size+number_corrections<max_subspace) THEN
          number_corrections=number_corrections+1
          subspace_vectors(:,subspace_size+number_corrections)=correction/correction_norm
        END IF
      END DO
      IF(number_corrections==0) EXIT
      subspace_size=subspace_size+number_corrections
      CALL apply_hamiltonian_block(H,subspace_vectors,subspace_hvectors,subspace_size,psi_work,hpsi_work,debug)
    END DO
    ! Un sous-espace Davidson peut être bloqué par une diagonale presque
    ! dégénérée. Dans ce cas, on ne retourne jamais des Ritz non convergés :
    ! Lanczos fournit une solution contrôlée avec le même opérateur H|psi>.
    IF(info/=0) THEN
      IF(debug) WRITE(out_unit,*) 'Davidson block not converged; controlled Lanczos fallback.'
      CALL solve_tise_lanczos(H,number_states,max_iter,tolerance,tise_energies,tise_eigenvectors,info,debug)
      IF(info==0) info=2
    END IF
    IF(debug) WRITE(out_unit,*) 'solve_tise_davidson_block iterations/residual = ',index_iteration,MAXVAL(residual_norms)

    CALL dealloc_psi(psi_work)
    CALL dealloc_psi(hpsi_work)

    DEALLOCATE(subspace_vectors,subspace_hvectors,diagonal_hamiltonian,correction,ritz,ritz_h,residual,small_energies,small_vectors,residual_norms)

  END SUBROUTINE solve_tise_davidson_block

  !> @brief Calcule la diagonale de H par applications à des vecteurs de base.
  SUBROUTINE compute_hamiltonian_diagonal(H,diagonal_hamiltonian,psi_work,hpsi_work,debug)

    TYPE(Op_t), INTENT(IN)        :: H
    REAL(KIND=Rkind), INTENT(OUT) :: diagonal_hamiltonian(:)
    TYPE(psi_t), INTENT(INOUT)    :: psi_work,hpsi_work
    LOGICAL, INTENT(IN)           :: debug
    COMPLEX(KIND=Rkind)           :: psi_basis(SIZE(diagonal_hamiltonian)),hpsi_basis(SIZE(diagonal_hamiltonian))
    INTEGER                       :: index_basis

    DO index_basis=1,SIZE(diagonal_hamiltonian)
      psi_basis=CZERO
      psi_basis(index_basis)=CONE

      CALL Apply_Hamiltonian(H,psi_basis,hpsi_basis,psi_work,hpsi_work,.FALSE.)

      diagonal_hamiltonian(index_basis)=REAL(hpsi_basis(index_basis),KIND=Rkind)
    END DO

    IF(debug) CALL Write_VecMat(diagonal_hamiltonian,out_unit,MIN(8,SIZE(diagonal_hamiltonian)),info='TISE diagonal preconditioner')

  END SUBROUTINE compute_hamiltonian_diagonal

  !> @brief Orthogonalise les colonnes initiales du bloc.
  SUBROUTINE orthonormalize_block(subspace_vectors,subspace_size,debug)

    COMPLEX(KIND=Rkind), INTENT(INOUT) :: subspace_vectors(:,:)
    INTEGER, INTENT(IN)                :: subspace_size
    LOGICAL, INTENT(IN)                :: debug
    INTEGER                            :: index_vector,index_basis
    REAL(KIND=Rkind)                   :: vector_norm

    DO index_vector=1,subspace_size
      DO index_basis=1,index_vector-1
        subspace_vectors(:,index_vector)=subspace_vectors(:,index_vector)-DOT_PRODUCT(subspace_vectors(:,index_basis),&
          subspace_vectors(:,index_vector))*subspace_vectors(:,index_basis)
      END DO

      vector_norm=SQRT(MAX(ZERO,REAL(DOT_PRODUCT(subspace_vectors(:,index_vector),subspace_vectors(:,index_vector)),KIND=Rkind)))
      IF(vector_norm<=EPSILON(ONE)) ERROR STOP 'orthonormalize_block: dependent block.'
      subspace_vectors(:,index_vector)=subspace_vectors(:,index_vector)/vector_norm
    END DO
    IF(debug) WRITE(out_unit,*) 'orthonormalize_block size = ',subspace_size

  END SUBROUTINE orthonormalize_block

  !> @brief Orthogonalise une correction contre l'espace courant.
  SUBROUTINE orthonormalize_against_block(correction,subspace_vectors,subspace_size,debug)

    COMPLEX(KIND=Rkind), INTENT(INOUT)   :: correction(:)
    COMPLEX(KIND=Rkind), INTENT(IN)      :: subspace_vectors(:,:)
    INTEGER, INTENT(IN)                  :: subspace_size
    LOGICAL, INTENT(IN)                  :: debug
    INTEGER                              :: index_vector

    DO index_vector=1,subspace_size
      correction=correction-DOT_PRODUCT(subspace_vectors(:,index_vector),correction)*subspace_vectors(:,index_vector)
    END DO
    IF(debug) WRITE(out_unit,*) 'orthonormalize correction norm = ',SQRT(MAX(ZERO,REAL(DOT_PRODUCT(correction,correction),KIND=Rkind)))

  END SUBROUTINE orthonormalize_against_block

  !> @brief Applique H à toutes les colonnes d'un bloc.
  SUBROUTINE apply_hamiltonian_block(H,subspace_vectors,subspace_hvectors,subspace_size,psi_work,hpsi_work,debug)

    TYPE(Op_t), INTENT(IN)           :: H
    COMPLEX(KIND=Rkind), INTENT(IN)  :: subspace_vectors(:,:)
    COMPLEX(KIND=Rkind), INTENT(OUT) :: subspace_hvectors(:,:)
    INTEGER, INTENT(IN)              :: subspace_size
    TYPE(psi_t), INTENT(INOUT)       :: psi_work,hpsi_work
    LOGICAL, INTENT(IN)              :: debug
    INTEGER                          :: index_vector

    DO index_vector=1,subspace_size

      CALL Apply_Hamiltonian(H,subspace_vectors(:,index_vector),subspace_hvectors(:,index_vector),psi_work,hpsi_work,.FALSE.)

    END DO
    IF(debug) WRITE(out_unit,*) 'apply_hamiltonian_block size = ',subspace_size

  END SUBROUTINE apply_hamiltonian_block

  !> @brief Diagonalise le Hamiltonien projeté V^H H V.
  SUBROUTINE diagonalize_projected_problem(subspace_vectors,subspace_hvectors,subspace_size,small_energies,small_vectors,debug)

    COMPLEX(KIND=Rkind), INTENT(IN)    :: subspace_vectors(:,:),subspace_hvectors(:,:)
    INTEGER, INTENT(IN)                :: subspace_size
    REAL(KIND=Rkind), INTENT(OUT)      :: small_energies(:),small_vectors(:,:)
    LOGICAL, INTENT(IN)                :: debug
    REAL(KIND=Rkind), ALLOCATABLE      :: projected_matrix(:,:)
    INTEGER                            :: index_row,index_column

    ALLOCATE(projected_matrix(subspace_size,subspace_size))
    projected_matrix=ZERO

    DO index_column=1,subspace_size
      DO index_row=1,subspace_size
        projected_matrix(index_row,index_column)=REAL(DOT_PRODUCT(subspace_vectors(:,index_row),subspace_hvectors(:,index_column)),KIND=Rkind)
      END DO
    END DO

    CALL diagonalization(projected_matrix,small_energies(1:subspace_size),small_vectors(1:subspace_size,1:subspace_size),sort=1)
    IF(debug) CALL Write_VecMat(projected_matrix,out_unit,MIN(8,subspace_size),info='TISE projected matrix')

    DEALLOCATE(projected_matrix)

  END SUBROUTINE diagonalize_projected_problem

  !> @brief Diagonalise la matrice tridiagonale de Lanczos.
  SUBROUTINE diagonalize_tridiagonal(alpha,beta,krylov_size,eigenvalues,eigenvectors)

    REAL(KIND=Rkind), INTENT(IN)       :: alpha(:),beta(:)
    INTEGER, INTENT(IN)                :: krylov_size
    REAL(KIND=Rkind), INTENT(OUT)      :: eigenvalues(:),eigenvectors(:,:)
    REAL(KIND=Rkind), ALLOCATABLE      :: tridiagonal_matrix(:,:)
    INTEGER                            :: index_basis

    ALLOCATE(tridiagonal_matrix(krylov_size,krylov_size))
    tridiagonal_matrix=ZERO

    DO index_basis=1,krylov_size
      tridiagonal_matrix(index_basis,index_basis)=alpha(index_basis)

      IF(index_basis<krylov_size) THEN
        ! Upper diagonal
        tridiagonal_matrix(index_basis,index_basis+1)=beta(index_basis)
        ! Lower diagonal
        tridiagonal_matrix(index_basis+1,index_basis)=beta(index_basis)

      END IF
    END DO

    CALL diagonalization(tridiagonal_matrix,eigenvalues(1:krylov_size),eigenvectors(1:krylov_size,1:krylov_size),sort=1)

    DEALLOCATE(tridiagonal_matrix)

  END SUBROUTINE diagonalize_tridiagonal

  !> @brief Résout la TISE par reverse communication ARPACK.
  SUBROUTINE solve_tise_arpack(H,number_states,max_iter,tolerance,&
    tise_energies,tise_eigenvectors,info,debug)

    TYPE(Op_t), INTENT(INOUT)        :: H
    INTEGER, INTENT(IN)              :: number_states,max_iter
    REAL(KIND=Rkind), INTENT(IN)     :: tolerance
    REAL(KIND=Rkind), INTENT(OUT)    :: tise_energies(:)
    COMPLEX(KIND=Rkind), INTENT(OUT) :: tise_eigenvectors(:,:)
    INTEGER, INTENT(OUT)             :: info
    LOGICAL, INTENT(IN)              :: debug

#if __ARPACK == 1
    INTEGER                          :: dimension_hamiltonian,ncv,ldv,lworkl,ido,arpack_info
    INTEGER                          :: iparam(11),ipntr(11)
    REAL(KIND=Rkind), ALLOCATABLE    :: resid(:),arpack_vectors(:,:),workd(:),workl(:),eigenvalues(:,:)
    REAL(KIND=Rkind)                 :: sigma
    LOGICAL, ALLOCATABLE             :: select(:)
    COMPLEX(KIND=Rkind), ALLOCATABLE :: psi_request(:),hpsi_response(:)
    TYPE(psi_t)                      :: psi_work,hpsi_work

    dimension_hamiltonian=SIZE(tise_eigenvectors,1)

    IF(number_states>=dimension_hamiltonian) THEN
      CALL solve_tise_dense(H,number_states,tise_energies,tise_eigenvectors,debug)
      info=0
      RETURN
    END IF

    ncv=MIN(dimension_hamiltonian,MAX(number_states+2,2*number_states+8))
    ldv=dimension_hamiltonian
    lworkl=ncv*(ncv+8)

    ALLOCATE(resid(dimension_hamiltonian),arpack_vectors(ldv,ncv))
    ALLOCATE(workd(3*dimension_hamiltonian),workl(lworkl))
    ALLOCATE(eigenvalues(ncv,2),select(ncv))
    ALLOCATE(psi_request(dimension_hamiltonian),hpsi_response(dimension_hamiltonian))

    CALL init_psi(psi_work,H%Basis,cplx=.TRUE.,grid=.FALSE.)
    CALL init_psi(hpsi_work,H%Basis,cplx=.TRUE.,grid=.FALSE.)

    ido=0
    arpack_info=0
    sigma=ZERO
    resid=ZERO
    arpack_vectors=ZERO
    workd=ZERO
    workl=ZERO
    eigenvalues=ZERO
    select=.FALSE.
    iparam=0
    iparam(1)=1
    iparam(3)=max_iter
    iparam(7)=1
    DO
      CALL dsaupd(ido,'I',dimension_hamiltonian,'SA',number_states,tolerance,resid,ncv,&
        arpack_vectors,ldv,iparam,ipntr,workd,workl,lworkl,arpack_info)
      IF(ido==-1.OR.ido==1) THEN
        psi_request=CMPLX(workd(ipntr(1):ipntr(1)+dimension_hamiltonian-1),ZERO,KIND=Rkind)
        CALL Apply_Hamiltonian(H,psi_request,hpsi_response,psi_work,hpsi_work,.FALSE.)
        workd(ipntr(2):ipntr(2)+dimension_hamiltonian-1)=REAL(hpsi_response,KIND=Rkind)
      ELSE
        EXIT
      END IF
    END DO

    IF(arpack_info<0) ERROR STOP 'solve_tise_arpack: dsaupd failed.'
    CALL dseupd(.TRUE.,'A',select,eigenvalues,arpack_vectors,ldv,sigma,'I',dimension_hamiltonian,'SA',&
      number_states,tolerance,resid,ncv,arpack_vectors,ldv,iparam,ipntr,workd,workl,lworkl,arpack_info)
    IF(arpack_info/=0) ERROR STOP 'solve_tise_arpack: dseupd failed.'

    tise_energies=eigenvalues(1:number_states,1)
    tise_eigenvectors=CMPLX(arpack_vectors(:,1:number_states),ZERO,KIND=Rkind)
    info=0

    IF(debug) WRITE(out_unit,*) 'solve_tise_arpack iterations = ',iparam(3)
    CALL dealloc_psi(psi_work)
    CALL dealloc_psi(hpsi_work)

    DEALLOCATE(resid,arpack_vectors,workd,workl,eigenvalues,select,psi_request,hpsi_response)
#else
    info=1
    tise_energies=ZERO
    tise_eigenvectors=CZERO
    WRITE(out_unit,*) 'TISE ARPACK requested but executable was built with ARPACK=0.'
    ERROR STOP 'solve_tise_arpack: rebuild with ARPACK=1.'
#endif
  END SUBROUTINE solve_tise_arpack

  !> @brief Relaxation en temps imaginaire vers l'état fondamental.
  !> @details Met à jour psi <- psi - delta_tau H psi et normalise à chaque pas.
  SUBROUTINE Relax_TISE(H,psi_relaxed,delta_tau,max_iter,tolerance,debug)

    TYPE(Op_t), INTENT(IN)         :: H
    TYPE(psi_t), INTENT(INOUT)     :: psi_relaxed
    REAL(KIND=Rkind), INTENT(IN)   :: delta_tau,tolerance
    INTEGER, INTENT(IN)            :: max_iter
    LOGICAL, INTENT(IN), OPTIONAL  :: debug
    TYPE(psi_t)                    :: hpsi_relaxed
    REAL(KIND=Rkind)               :: energy_old,energy_new,norm_relaxed
    INTEGER                        :: iteration
    LOGICAL                        :: debug_tise,converged

    debug_tise=Debug_Enabled()

    IF(PRESENT(debug)) debug_tise=debug
    IF(delta_tau<=ZERO.OR.tolerance<=ZERO.OR.max_iter<1) ERROR STOP 'Relax_TISE: invalid parameters.'
    CALL init_psi(hpsi_relaxed,psi_relaxed%Basis,cplx=.TRUE.,grid=.FALSE.)
    CALL calculate_energy(energy_old,psi_relaxed,H,hpsi_relaxed)

    converged=.FALSE.

    DO iteration=1,max_iter
      CALL calc_OpPsi(H,psi_relaxed,hpsi_relaxed)
      psi_relaxed%CVec=psi_relaxed%CVec-delta_tau*hpsi_relaxed%CVec
      CALL Calc_Norm_OF_Psi(psi_relaxed,norm_relaxed)
      psi_relaxed%CVec=psi_relaxed%CVec/norm_relaxed
      CALL calculate_energy(energy_new,psi_relaxed,H,hpsi_relaxed)
      IF(debug_tise) WRITE(out_unit,*) 'Relax_TISE iteration/energy/delta = ',iteration,energy_new,ABS(energy_new-energy_old)
      IF(ABS(energy_new-energy_old)<=tolerance) THEN
        converged=.TRUE.
        EXIT
      END IF
      energy_old=energy_new
    END DO

    IF(debug_tise) WRITE(out_unit,*) 'Relax_TISE converged = ',converged
    CALL dealloc_psi(hpsi_relaxed)

  END SUBROUTINE Relax_TISE

  !> @brief Calcule <psi|H|psi>/<psi|psi> avec l’interface existante.
  !> @details Hpsi est fourni en argument pour éviter une seconde application
  !>          de H pendant la relaxation.
  SUBROUTINE calculate_energy(energy,psi_state,H,hpsi_state)

    REAL(KIND=Rkind), INTENT(OUT)  :: energy
    TYPE(psi_t), INTENT(IN)        :: psi_state
    TYPE(Op_t), INTENT(IN)         :: H
    TYPE(psi_t), INTENT(INOUT)     :: hpsi_state
    COMPLEX(KIND=Rkind)            :: scalar_energy
    REAL(KIND=Rkind)               :: norm_state

    CALL calc_OpPsi(H,psi_state,hpsi_state)
    scalar_energy=DOT_PRODUCT(psi_state%CVec,hpsi_state%CVec)
    CALL Calc_Norm_OF_Psi(psi_state,norm_state)

    energy=REAL(scalar_energy,KIND=Rkind)/(norm_state*norm_state)

  END SUBROUTINE calculate_energy

  !> @brief Ecrit énergie, vecteurs propres et états de transition.
  !> @details La transition i contient Delta E_i-E_0 et le poids
  !> |<psi_initial|psi_i>|^2 lorsque psi_initial est fourni.
  SUBROUTINE write_tise_results(tise_energies,tise_eigenvectors,unit,energy_file,transition_file,psi_initial,info,debug)

    REAL(KIND=Rkind), INTENT(IN)      :: tise_energies(:)
    COMPLEX(KIND=Rkind), INTENT(IN)   :: tise_eigenvectors(:,:)
    CHARACTER(*), INTENT(IN)          :: unit,energy_file,transition_file
    TYPE(psi_t), INTENT(IN), OPTIONAL :: psi_initial
    INTEGER, INTENT(IN)               :: info
    LOGICAL, INTENT(IN)               :: debug
    REAL(KIND=Rkind)                  :: conversion,transition_weight
    COMPLEX(KIND=Rkind)               :: transition_amplitude
    INTEGER                           :: state_index,vector_index,unit_energy,unit_vector,unit_transition

    conversion=ONE

    IF(TRIM(unit)=='cm-1') conversion=get_Conv_au_TO_unit('E','cm-1')

    OPEN(NEWUNIT=unit_energy,FILE=TRIM(energy_file),STATUS='REPLACE',ACTION='WRITE')

    WRITE(unit_energy,*) '# state energy_'//TRIM(unit)

    DO state_index=1,SIZE(tise_energies)
       WRITE(unit_energy,*) state_index,tise_energies(state_index)*conversion
    END DO

    CLOSE(unit_energy)
    OPEN(NEWUNIT=unit_vector,FILE='TISE_eigenvectors.txt',STATUS='REPLACE',ACTION='WRITE')
    WRITE(unit_vector,*) '# dimension states ',SIZE(tise_eigenvectors,2)
    
    DO vector_index=1,SIZE(tise_eigenvectors,1)
       WRITE(unit_vector,*) (REAL(tise_eigenvectors(vector_index,state_index),KIND=Rkind),state_index=1,SIZE(tise_eigenvectors,2))
    END DO

    CLOSE(unit_vector)
    OPEN(NEWUNIT=unit_transition,FILE=TRIM(transition_file),STATUS='REPLACE',ACTION='WRITE')
    WRITE(unit_transition,*) '# state transition_energy_'//TRIM(unit)//' transition_weight'

    DO state_index=1,SIZE(tise_energies)

      transition_weight=ZERO

      IF(PRESENT(psi_initial)) THEN
         transition_amplitude=DOT_PRODUCT(psi_initial%CVec,tise_eigenvectors(:,state_index))
         transition_weight=ABS(transition_amplitude)**2
      END IF
      WRITE(unit_transition,*) state_index,(tise_energies(state_index)-tise_energies(1))*conversion,transition_weight
    END DO

    CLOSE(unit_transition)
    WRITE(out_unit,*) 'TISE completed: ',SIZE(tise_energies),' states, unit=',TRIM(unit)

    IF(info/=0) WRITE(out_unit,*) 'TISE warning, iterative status = ',info
    IF(debug) WRITE(out_unit,*) 'write_tise_results transition file = ',TRIM(transition_file)

  END SUBROUTINE write_tise_results

END MODULE TISE_m
