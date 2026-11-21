!----------------------------------------------------------------------------------
!> MODULE: sub_propa_m
!> Description: This module manages the time propagation of wavefunctions
!>              in the context of molecular dynamics. It provides an interface
!>              for various time integration methods, encapsulates simulation
!>              parameters in a data type, and offers utility routines for
!>              analysis and data handling.
!----------------------------------------------------------------------------------
MODULE sub_propa_m
   USE QDUtil_m           
   USE psi_m
   USE Op_m
   USE Ana_psi_m
   USE lanczos_m
   USE Sub_Vp_m
   USE Basis_m
   USE SmolyakSparseBasisAlgebra_m, ONLY: DotSmolyakBasis
   USE TISE_m, ONLY: Build_Hamiltonian_Matrix
   
   implicit none

   !----------------------------------------------------------------------------------
   !> TYPE: propa_t
   !> Description: Contains the parameters for the time propagation.
   !>   t0            : Starting time of the simulation.
   !>   tf            : Final time of the simulation.
   !>   delta_t       : Time step for the integration.
   !>   eps           : Tolerance for convergence criteria.
   !>   E_min         : Minimum energy in Chebyshev expansion.
   !>   E_max         : Maximum energy in Chebyshev expansion.
   !>   max_iter      : Maximum number of iterations.
   !>   propa_name    : Name of the propagation method (standard/Hagedorn).
   !>   propa_name2   : Name of the propagation method, used in the 'march' routine.
   !>   Beta          : Logical flag for beta calculations.
   !>   P             : Logical flag for propagation.
   !>   renorm        : Logical flag for renormalization.
   !----------------------------------------------------------------------------------
   TYPE propa_t
      REAL(kind=Rkind)               :: t0
      REAL(kind=Rkind)               :: tf
      REAL(kind=Rkind)               :: delta_t
      REAL(kind=Rkind)               :: eps
      REAL(kind=Rkind)               :: E_min
      REAL(kind=Rkind)               :: E_max
      INTEGER                        :: max_iter
      INTEGER                        :: split_order = 2
      CHARACTER(len=:), ALLOCATABLE  :: propa_name
      CHARACTER(len=:), ALLOCATABLE  :: propa_name2
      LOGICAL                        :: Beta 
      LOGICAL                        :: P
      LOGICAL                        :: renorm
      LOGICAL                        :: action_H = .FALSE.
      LOGICAL                        :: writeV = .FALSE.
      CHARACTER(LEN=16)              :: dynamics = 'TDSE'
      CHARACTER(LEN=16)              :: tise_method = 'dense'
      CHARACTER(LEN=16)              :: energy_unit = 'au'
      CHARACTER(LEN=512)             :: tise_energy_file = 'TISE_energies.dat'
      CHARACTER(LEN=512)             :: tise_transition_file = 'TISE_transitions.txt'
      INTEGER                        :: tise_nstates = 1
      INTEGER                        :: tise_max_iter = 100
      REAL(KIND=Rkind)               :: tise_tol = ONETENTH**8
      LOGICAL                        :: tise_relax = .FALSE.
      REAL(KIND=Rkind)               :: tise_relax_delta_t = ONETENTH**3
      INTEGER                        :: tise_relax_max_iter = 500
      REAL(KIND=Rkind)               :: tise_relax_tol = ONETENTH**10
      LOGICAL                        :: debug = .FALSE.
      LOGICAL                        :: write_Qt = .FALSE.
      LOGICAL                        :: write_SQt = .FALSE.
      LOGICAL                        :: write_Pt = .FALSE.
      LOGICAL                        :: write_E = .FALSE.
      LOGICAL                        :: write_Norm = .FALSE.
      LOGICAL                        :: write_pop = .FALSE.
      LOGICAL                        :: write_alpha = .FALSE.
      LOGICAL                        :: write_auto_cor = .FALSE.
      LOGICAL                        :: write_validation = .FALSE.
      LOGICAL                        :: write_E_13 = .FALSE.
      LOGICAL                        :: write_Norm_13 = .FALSE.
      LOGICAL                        :: restart = .FALSE.
      LOGICAL                        :: checkpoint = .TRUE.
      INTEGER(KIND=Ikind)            :: checkpoint_every = 1_Ikind
      CHARACTER(LEN=512)             :: checkpoint_file = 'TD_SCHROD.restart'
   END TYPE propa_t

   ! Define public interfaces
   PUBLIC :: march_taylor, march
   PUBLIC :: write_propa, Analyse, open_propa_file, read_propa
   PUBLIC :: CopyPsi

CONTAINS


   !----------------------------------------------------------------------------------
   !> SUBROUTINE: march
   !> Description: This is a dispatcher subroutine that selects the appropriate
   !>              time propagation routine based on the method name specified in
   !>              the 'propa_t' type. It provides a centralized interface for
   !>              starting the propagation step.
   !> Arguments:
   !>   - psi       : (IN)    Initial wavefunction.
   !>   - psi_dt    : (INOUT) Wavefunction after propagation.
   !>   - t         : (IN)    Current simulation time.
   !>   - propa     : (IN)    Structure containing propagation parameters.
   !>   - H         : (IN)    Hamiltonian operator.
   !----------------------------------------------------------------------------------
   SUBROUTINE march(psi, psi_dt, t, propa, H)
      TYPE(propa_t), INTENT(IN)                :: propa
      TYPE(psi_t), INTENT(IN)                  :: psi
      TYPE(psi_t), INTENT(INOUT)               :: psi_dt
      TYPE(Op_t), INTENT(IN)                   :: H
      REAL(KIND=Rkind), INTENT(IN)             :: t
      CHARACTER(len=Name_len)                  :: name

      name = propa%propa_name2 
      CALL string_uppercase_TO_lowercase(name)
      SELECT CASE (name)
      CASE ('rk1') 
        CALL march_RK1(psi, psi_dt, t, propa, H)
      CASE ('rk2') 
        CALL march_RK2(psi, psi_dt, t, propa, H)
      CASE ('rk3') 
        CALL march_RK3(psi, psi_dt, t, propa, H)
      CASE ('rk4') 
        CALL march_RK4(psi, psi_dt, t, propa, H)
      CASE ('taylor')
         CALL march_taylor(psi, psi_dt, t, propa, H)
      CASE ('sil') 
        CALL march_SIL(psi, psi_dt, t, propa,H)
      CASE ('split')
        CALL march_SPLIT(psi, psi_dt, t, propa, H)
      CASE ('cn')
        CALL march_CN(psi, psi_dt, t, propa, H)
      CASE ('cheby')
        CALL march_Chebychev(psi, psi_dt, t, propa, H)
      CASE ('vp')
         CALL march_VP(psi, psi_dt, t, propa, H)
         !CALL march_VP_RK4(psi, psi_dt, t, propa, H)
         !CALL march_VPTMP(psi, psi_dt, t, propa, H)
      CASE DEFAULT
         ERROR STOP 'march: unsupported propagation method (read_propa should validate it).'
      END SELECT      
   END SUBROUTINE


   !----------------------------------------------------------------------------------
   !> SUBROUTINE: CopyPsi
   !> Description: Copies the data of one wavefunction (`psi1`) to another (`psi2`).
   !>              The user-provided implementation contained non-functional code
   !>              that was likely left over from previous development. It has been
   !>              cleaned up and simplified to perform a straightforward copy
   !>              of the coefficient vectors, as this is the standard and most
   !>              efficient way to copy wavefunction data.
   !> Arguments:
   !>   - psi2 : (INOUT) Destination wavefunction.
   !>   - psi1 : (IN)    Source wavefunction.
   !----------------------------------------------------------------------------------
   SUBROUTINE CopyPsi(psi2,psi1)
      IMPLICIT NONE
      TYPE(psi_t), INTENT(IN)                    :: psi1
      TYPE(psi_t), INTENT(INOUT)                 :: psi2
      ! The original routine was non-functional; for a simple copy,
      ! the coefficient vector is copied directly.
      psi2%CVec(:) = psi1%CVec(:)
   END SUBROUTINE


   !----------------------------------------------------------------------------------
   !> SUBROUTINE: Analyse
   !> Description: Performs an analysis of the wavefunction at a given time 't'.
   !>              It calculates the population, average energy, and norm, and
   !>              writes these results to specific file units.
   !> Arguments:
   !>   - psi : (IN) Wavefunction to be analyzed.
   !>   - t   : (IN) Current time.
   !----------------------------------------------------------------------------------
   SUBROUTINE Analyse(psi, t)
      IMPLICIT NONE
      TYPE(psi_t), INTENT(IN)              :: psi
      REAL(KIND=Rkind), INTENT(IN)         :: t

      !local variables --------------------------------------------------------
      REAL(KIND=Rkind), ALLOCATABLE        :: pop(:), Qm(:), Qp(:)
      REAL(KIND=Rkind)                     :: Norm, E
      INTEGER                              :: Ndim
      INTEGER, SAVE                        :: u_pop_an = -1, u_Qm_an = -1, u_Qp_an = -1
      LOGICAL, SAVE                        :: analysis_files_open = .FALSE.

      !allocate variables ----------------------------------------------------
      Ndim = size(psi%Basis%tab_basis)
      ALLOCATE (Pop(psi%Basis%tab_basis(Ndim)%nb))
      ALLOCATE (Qp(psi%Basis%tab_basis(Ndim)%nb))
      ALLOCATE (Qm(Ndim - 1))

      IF (.NOT. analysis_files_open) THEN
         CALL open_propa_file(u_pop_an, 'population')
         CALL open_propa_file(u_Qm_an, 'Qm')
         CALL open_propa_file(u_Qp_an, 'Qp')
         analysis_files_open = .TRUE.
      END IF

      ! Initialize variables --------------------------------------------------
      pop(:) = ZERO
      Qm(:) = ZERO
      Qp(:) = ZERO
      E = ZERO; Norm = ZERO
      
      !--- Start analysis ---
      call Population(psi, pop)
      call Calc_average_energy(psi, E)
      call Calc_Norm_OF_psi(psi, Norm)
      call Qpop(psi, Qp)
      
      WRITE (u_pop_an, *) t, E, Norm, pop
      WRITE (u_Qm_an, *) t, Qm
      WRITE (u_Qp_an, *) t, Qp
      
      !deallocate variables --------------------------------------------------
      DEALLOCATE (pop, Qm, Qp)
      !--- End analysis ---

   END SUBROUTINE Analyse

   !=====================================================================
   !> SUBROUTINE: march_taylor
   !! Description: Performs time propagation of the wavefunction `psi`
   !! using the Taylor series expansion method. The procedure is iterative,
   !! applying the Hamiltonian and adding scaled terms to the wavefunction
   !! until a convergence threshold is met or the maximum number of
   !! iterations is reached.
   !!
   !! Arguments:
   !!   - psi     : (IN)     Initial wavefunction
   !!   - psi_dt  : (INOUT)  Wavefunction after propagation
   !!   - t       : (IN)     Current time
   !!   - propa   : (IN)     Propagation parameters
   !!   - H       : (IN)     Hamiltonian operator
   !=====================================================================
   SUBROUTINE march_taylor(psi, psi_dt, t, propa, H)
      TYPE(psi_t), INTENT(INOUT)       :: psi_dt
      TYPE(psi_t), INTENT(IN)          :: psi
      TYPE(propa_t), INTENT(IN)        :: propa
      TYPE(Op_t), INTENT(IN)           :: H
      REAL(KIND=Rkind), INTENT(IN)     :: t

      TYPE(psi_t)                      :: psi0
      TYPE(psi_t)                      :: Hpsi
      REAL(KIND=Rkind)                 :: alpha, Norm, Norm0, dt_sub
      INTEGER                          :: kk, isub, nsub, retry
      INTEGER, PARAMETER               :: max_retry = 12
      LOGICAL                          :: debug, attempt_ok, converged

      CALL init_psi(Hpsi, psi%basis, cplx=.true., grid=.FALSE.)
      CALL init_psi(psi0, psi%basis, cplx=.true., grid=.FALSE.)
      alpha = TEN**15
      debug = propa%debug
      nsub = 1

      IF (debug) THEN
         WRITE(out_unit, *) 'BEGINNING march_taylor ', t, ' delta_t=', propa%delta_t
         FLUSH(out_unit)
      END IF

      ! Try the requested step first. If the Taylor terms become too large or
      ! the truncation tolerance is not reached, retry the whole step with
      ! twice as many internal substeps. The caller still observes one result
      ! at the requested outer time grid.
      DO retry = 0, max_retry
         dt_sub = propa%delta_t / REAL(nsub, KIND=Rkind)
         psi_dt%CVec = psi%CVec
         attempt_ok = .TRUE.

         DO isub = 1, nsub
            psi0%CVec = psi_dt%CVec
            converged = .FALSE.

            DO kk = 1, propa%max_iter
               CALL mEyeHPsi(psi0, Hpsi, H)
               psi0%CVec(:) = Hpsi%CVec(:) * (dt_sub / REAL(kk, KIND=Rkind))
               psi_dt%CVec(:) = psi_dt%CVec(:) + psi0%CVec(:)
               CALL Calc_Norm_OF_Psi(psi0, Norm)
               IF (debug) WRITE(out_unit, *) 'Taylor substep/order/norm =', isub, kk, Norm

               IF (Norm >= alpha) THEN
                  attempt_ok = .FALSE.
                  EXIT
               ELSE IF (Norm <= propa%eps) THEN
                  converged = .TRUE.
                  EXIT
               END IF
            END DO

            IF (.NOT. attempt_ok .OR. .NOT. converged) THEN
               attempt_ok = .FALSE.
               EXIT
            END IF
         END DO

         IF (attempt_ok) EXIT
         IF (debug) WRITE(out_unit, *) 'Retrying Taylor step with', 2*nsub, 'substeps'
         nsub = 2 * nsub
      END DO

      IF (.NOT. attempt_ok) THEN
         ERROR STOP 'march_taylor: no stable/converged substep; reduce delta_t or increase max_iter'
      END IF

      IF (debug) THEN
         CALL Calc_Norm_OF_Psi(psi, Norm0)
         CALL Calc_Norm_OF_Psi(psi_dt, Norm)
         WRITE(out_unit, *) 'Taylor internal substeps =', nsub
         WRITE(out_unit, *) '<psi_dt|psi_dt> = ', Norm, ' |delta Norm| =', ABS(Norm0 - Norm)
         WRITE(out_unit, *) 'END march_taylor'
         FLUSH(out_unit)
      END IF

      CALL dealloc_psi(psi0)
      CALL dealloc_psi(Hpsi)
   END SUBROUTINE march_taylor
  

   !----------------------------------------------------------------------------------
   !> SUBROUTINE: march_RK1 (Forward Euler)
   !! Description: First-order explicit integrator.
   !!              Mathematical form: psi_dt = psi + Δt * f(psi)
   !!              with f(psi) = -iHpsi.
   !!
   !! Algorithmic steps:
   !!   1. Evaluate K1 = f(psi) (one application of H).
   !!   2. Update psi_dt = psi + dt * K1.
   !!
   !! Arguments:
   !!   - psi     : (IN)    Initial wavefunction.
   !!   - psi_dt  : (INOUT) Wavefunction after propagation.
   !!   - t       : (IN)    Current time.
   !!   - propa   : (IN)    Propagation parameters.
   !!   - H       : (IN)    Hamiltonian operator.
   !----------------------------------------------------------------------------------
   
   SUBROUTINE march_RK1(psi, psi_dt, t, propa, H)
      TYPE(psi_t), INTENT(INOUT)       :: psi_dt
      TYPE(psi_t), INTENT(IN)          :: psi
      TYPE(propa_t), INTENT(IN)        :: propa
      TYPE(Op_t), INTENT(IN)           :: H
      REAL(KIND=Rkind), INTENT(IN)     :: t
   
      !local variables-------------------------------------------------------------------------------
      TYPE(psi_t) :: K1
      REAL(KIND=Rkind) :: Norm0, Norm
      LOGICAL :: debug

      debug = propa%debug
      
      IF (debug) THEN
      WRITE (out_unit, *) 'BEGINNING march_RK1  ', 't=', t, 'Δt=', propa%delta_t
      WRITE (out_unit, *) ''
      END IF
   
      CALL init_psi(K1, psi%basis, cplx=.true., grid=.false.)
   
      ! Compute derivative K1 = -i* H * psi
      CALL mEyeHPsi(psi, K1, H)
   
      ! Euler update: psi_dt = psi + dt * K1
      psi_dt%CVec(:) = psi%CVec(:) + propa%delta_t * K1%CVec(:)
   
      ! Optional: compute norms for diagnostics
      IF (debug) THEN
         CALL Calc_Norm_OF_Psi(psi, Norm0)
         CALL Calc_Norm_OF_Psi(psi_dt, Norm)
         WRITE (out_unit, *) '<psi_dt|psi_dt> = ', Norm, ' |ΔNorm| =', abs(Norm0 - Norm)
         WRITE (out_unit, *) 'END march_RK1'
         WRITE (out_unit, *) ''
         FLUSH(out_unit)
      END IF
       
      CALL dealloc_psi(K1)
   
   END SUBROUTINE march_RK1   

 !--------------------------------------------------------------
 ! RK2 (Midpoint) — second-order explicit integrator
 ! Description:
 !   - Mathematical form: two-stage midpoint RK. Let f(psi) = -i H psi.
 !     k1 = f(psi), psi_mid = psi + (Δt/2)*k1, k2 = f(psi_mid),
 !     psi_dt = psi + Δt * k2.
 !   - Algorithmic steps:
 !       1) Compute K1 = f(psi).
 !       2) Build psi_mid = psi + (Δt/2)*K1.
 !       3) Compute K2 = f(psi_mid).
 !       4) Update psi_dt = psi + Δt * K2.
 !   
 !--------------------------------------------------------------
SUBROUTINE march_RK2(psi, psi_dt, t, propa, H)
   TYPE(psi_t), INTENT(INOUT)       :: psi_dt
   TYPE(psi_t), INTENT(IN)          :: psi
   TYPE(propa_t), INTENT(IN)        :: propa
   TYPE(Op_t), INTENT(IN)           :: H
   REAL(KIND=Rkind), INTENT(IN)     :: t
   
   !local variables-------------------------------------------------------------------------------
   TYPE(psi_t) :: K1, K2, psi_mid
   REAL(KIND=Rkind) :: Norm0, Norm
   LOGICAL :: debug

   debug = propa%debug
   
   IF (debug) THEN
   WRITE (out_unit, *) 'BEGINNING march_RK2  ', 't=', t, 'Δt=', propa%delta_t
   WRITE (out_unit, *) ''
   END IF

   CALL init_psi(K1, psi%basis, cplx=.TRUE., grid=.FALSE.)
   CALL init_psi(K2, psi%basis, cplx=.TRUE., grid=.FALSE.)
   CALL init_psi(psi_mid, psi%basis, cplx=.TRUE., grid=.FALSE.)

   ! k1 = f(psi)
   CALL mEyeHPsi(psi, K1, H)

   ! mid point: psi_mid = psi + (Δt/2) * k1
   psi_mid%CVec(:) = psi%CVec(:) + HALF * propa%delta_t * K1%CVec(:)

   ! k2 = f(psi_mid)
   CALL mEyeHPsi(psi_mid, K2, H)

   ! update: psi_dt = psi + Δt * k2
   psi_dt%CVec(:) = psi%CVec(:) + propa%delta_t * K2%CVec(:)
   
   IF (debug) THEN
   CALL Calc_Norm_OF_Psi(psi, Norm0)
   CALL Calc_Norm_OF_Psi(psi_dt, Norm)
   WRITE (out_unit, *) '<psi_dt|psi_dt> = ', Norm, ' |ΔNorm| =', abs(Norm0 - Norm)
   WRITE (out_unit, *) 'END march_RK2'
   WRITE (out_unit, *) ''
   FLUSH(out_unit)
   END IF

   CALL dealloc_psi(K1)
   CALL dealloc_psi(K2)
   CALL dealloc_psi(psi_mid)

END SUBROUTINE march_RK2



 !--------------------------------------------------------------
 ! RK3 (Kutta third-order) — third-order explicit integrator
 ! Description:
 !   - Mathematical form: classical 3-stage Kutta rule implemented here.
 !     It uses a combination of k1,k2,k3 to obtain third-order accuracy.
 !   - Algorithmic steps (as implemented):
 !       1) K1 = f(psi).
 !       2) K2 = f(psi + Δt/2 * K1).
 !       3) K3 = f(psi - Δt * K1 + 2 Δt * K2).
 !       4) psi_dt = psi + (Δt/6) * (K1 + 4 K2 + K3).
 !--------------------------------------------------------------
SUBROUTINE march_RK3(psi, psi_dt, t, propa, H)
   TYPE(psi_t), INTENT(INOUT)       :: psi_dt
   TYPE(psi_t), INTENT(IN)          :: psi
   TYPE(propa_t), INTENT(IN)        :: propa
   TYPE(Op_t), INTENT(IN)           :: H
   REAL(KIND=Rkind), INTENT(IN)     :: t
   
   !local variables-------------------------------------------------------------------------------
   TYPE(psi_t) :: K1, K2, K3, psi_tmp
   REAL(KIND=Rkind) :: Norm0, Norm
   LOGICAL :: debug

   debug = propa%debug

   IF (debug) THEN
   WRITE (out_unit, *) 'BEGINNING march_RK3  ', 't=', t, 'Δt=', propa%delta_t
   WRITE (out_unit, *) ''
   END IF

   CALL init_psi(K1, psi%basis, cplx=.TRUE., grid=.FALSE.)
   CALL init_psi(K2, psi%basis, cplx=.TRUE., grid=.FALSE.)
   CALL init_psi(K3, psi%basis, cplx=.TRUE., grid=.FALSE.)
   CALL init_psi(psi_tmp, psi%basis, cplx=.TRUE., grid=.FALSE.)

   ! k1 = f(psi)
   CALL mEyeHPsi(psi, K1, H)

   ! k2 = f(psi + Δt/2 * k1)
   psi_tmp%CVec(:) = psi%CVec(:) + HALF * propa%delta_t * K1%CVec(:)
   CALL mEyeHPsi(psi_tmp, K2, H)

   ! k3 = f(psi - Δt * k1 + 2 Δt * k2)
   psi_tmp%CVec(:) = psi%CVec(:) - propa%delta_t * K1%CVec(:) + TWO * propa%delta_t * K2%CVec(:)
   CALL mEyeHPsi(psi_tmp, K3, H)

   ! combine
   psi_dt%CVec(:) = psi%CVec(:) + (propa%delta_t*SIXTH) * (K1%CVec(:) + FOUR*K2%CVec(:) + K3%CVec(:))
    
   IF (debug) THEN
   CALL Calc_Norm_OF_Psi(psi, Norm0)
   CALL Calc_Norm_OF_Psi(psi_dt, Norm)
   WRITE (out_unit, *) '<psi_dt|psi_dt> = ', Norm, ' |ΔNorm| =', abs(Norm0 - Norm)
   WRITE (out_unit, *) 'END march_RK3'
   WRITE (out_unit, *) ''
   FLUSH(out_unit)
   END IF

   ! deallocation
   CALL dealloc_psi(K1)
   CALL dealloc_psi(K2)
   CALL dealloc_psi(K3)
   CALL dealloc_psi(psi_tmp)

END SUBROUTINE march_RK3

!--------------------------------------------------------------
! RK4 (Classical fourth-order Runge-Kutta)
! Description:
!   - Mathematical form: four-stage RK with weights (1,2,2,1)/6.
!   - Algorithmic steps:
!       1) K1 = f(psi)
!       2) K2 = f(psi + Δt/2 * K1)
!       3) K3 = f(psi + Δt/2 * K2)
!       4) K4 = f(psi + Δt * K3)
!       5) psi_dt = psi + Δt/6 * (K1 + 2 K2 + 2 K3 + K4)
!--------------------------------------------------------------
SUBROUTINE march_RK4(psi, psi_dt, t, propa, H)
   TYPE(psi_t), INTENT(INOUT)       :: psi_dt
   TYPE(psi_t), INTENT(IN)          :: psi
   TYPE(propa_t), INTENT(IN)        :: propa
   TYPE(Op_t), INTENT(IN)           :: H
   REAL(KIND=Rkind), INTENT(IN)     :: t

   !local variables-------------------------------------------------------------------------------
   LOGICAL :: debug
   TYPE(psi_t) :: K1, K2, K3, K4, psi_tmp
   REAL(KIND=Rkind) :: Norm0, Norm

   CALL init_psi(K1, psi%basis, cplx=.TRUE., grid=.FALSE.)
   CALL init_psi(K2, psi%basis, cplx=.TRUE., grid=.FALSE.)
   CALL init_psi(K3, psi%basis, cplx=.true., grid=.FALSE.)
   CALL init_psi(K4, psi%basis, cplx=.true., grid=.FALSE.)
   CALL init_psi(psi_tmp, psi%basis, cplx=.true., grid=.FALSE.)

   debug = propa%debug

   IF (debug) THEN
   write (out_unit, *) 'BEGINNING march_RK4  ', 't=', t, 'Δt=', propa%delta_t
   WRITE (out_unit, *) ''
   END IF

   ! k1 = f(psi)
   CALL mEyeHPsi(psi, K1, H)

   ! k2 = f(psi + Δt/2 * k1)
   psi_tmp%CVec(:) = psi%CVec(:) + HALF * propa%delta_t * K1%CVec(:)
   CALL mEyeHPsi(psi_tmp, K2, H)

   ! k3 = f(psi + Δt/2 * k2)
   psi_tmp%CVec(:) = psi%CVec(:) + HALF * propa%delta_t * K2%CVec(:)
   CALL mEyeHPsi(psi_tmp, K3, H)

   ! k4 = f(psi + Δt * k3)
   psi_tmp%CVec(:) = psi%CVec(:) + propa%delta_t * K3%CVec(:)
   CALL mEyeHPsi(psi_tmp, K4, H)

   ! combine
   psi_dt%CVec(:) = psi%CVec(:) + (propa%delta_t*SIXTH) * (K1%CVec(:) + TWO*K2%CVec(:) + TWO*K3%CVec(:) + K4%CVec(:))

   IF (debug) THEN
   CALL Calc_Norm_OF_Psi(psi, Norm0)
   CALL Calc_Norm_OF_Psi(psi_dt, Norm)
   WRITE (out_unit, *) '<psi_dt|psi_dt> = ', Norm, ' |ΔNorm| =', abs(Norm0 - Norm)
   WRITE (out_unit, *) 'END march_RK4'
   WRITE (out_unit, *) ''
   FLUSH(out_unit)
   END IF

   ! deallocation
   CALL dealloc_psi(K1)
   CALL dealloc_psi(K2)
   CALL dealloc_psi(K3)
   CALL dealloc_psi(K4)
   CALL dealloc_psi(psi_tmp)

END SUBROUTINE march_RK4

SUBROUTINE march_SPLIT(psi, psi_dt, t, propa, H)
   ! Lie-Trotter, Strang, or symmetric higher-order composition on a
   ! periodic Fourier grid; the selected order is stored in propa%split_order.
   ! This is intentionally restricted to uniform periodic Fourier bases and
   ! constant Cartesian kinetic metrics; a Gaussian quadrature grid is not a
   ! Fourier grid and must not be silently treated as one.
   IMPLICIT NONE
   TYPE(psi_t), INTENT(IN) :: psi
   TYPE(psi_t), INTENT(INOUT) :: psi_dt
   REAL(KIND=Rkind), INTENT(IN) :: t
   TYPE(propa_t), INTENT(IN) :: propa
   TYPE(Op_t), INTENT(IN) :: H
   TYPE(psi_t) :: psi_grid
   INTEGER, ALLOCATABLE :: dims(:)
   REAL(KIND=Rkind), ALLOCATABLE :: GGdef(:,:), dx(:)
   INTEGER :: ndim, nq, nsurf, d, iq, stride
   CHARACTER(LEN=Name_len) :: basis_name
   REAL(KIND=Rkind) :: tol, norm0, norm1
   LOGICAL :: debug

   debug = propa%debug
   basis_name = H%Tdef%coord_type
   CALL string_uppercase_TO_lowercase(basis_name)
   IF (TRIM(basis_name) /= 'cc') &
      ERROR STOP 'march_SPLIT requires coord_type=CC (constant Cartesian kinetic operator).'
   IF (.NOT. ALLOCATED(psi%Basis%tab_basis)) ERROR STOP 'march_SPLIT: basis is not initialized.'
   IF (.NOT. ALLOCATED(H%Scal_pot)) ERROR STOP 'march_SPLIT: scalar potential is not initialized.'

   ndim = SIZE(psi%Basis%tab_basis)-1
   IF (ndim < 1) ERROR STOP 'march_SPLIT: at least one nuclear dimension is required.'
   nq = psi%Basis%nq
   nsurf = psi%Basis%tab_basis(ndim+1)%nb
   IF (nq < 1 .OR. nsurf < 1) ERROR STOP 'march_SPLIT: invalid grid/surface dimensions.'
   IF (SIZE(H%Scal_pot,1) /= nq .OR. SIZE(H%Scal_pot,2) /= nsurf .OR. &
       SIZE(H%Scal_pot,3) /= nsurf) ERROR STOP 'march_SPLIT: potential dimensions do not match the wavefunction.'
   IF (.NOT. ALLOCATED(psi%CVec) .OR. .NOT. ALLOCATED(psi_dt%CVec)) &
      ERROR STOP 'march_SPLIT: wavefunction coefficients are not allocated.'
   IF (SIZE(psi%CVec) /= nq*nsurf .OR. SIZE(psi_dt%CVec) /= nq*nsurf) &
      ERROR STOP 'march_SPLIT: wavefunction/grid dimensions do not match.'

   ALLOCATE(dims(ndim), GGdef(ndim,ndim), dx(ndim))
   stride = 1
   DO d = 1, ndim
      basis_name = psi%Basis%tab_basis(d)%Basis_name
      CALL string_uppercase_TO_lowercase(basis_name)
      IF (TRIM(basis_name) /= 'fourier') &
         ERROR STOP 'march_SPLIT requires name=fourier for every nuclear degree of freedom.'
      IF (.NOT. ALLOCATED(psi%Basis%tab_basis(d)%x) .OR. .NOT. ALLOCATED(psi%Basis%tab_basis(d)%w)) &
         ERROR STOP 'march_SPLIT: Fourier grid or quadrature weights are not allocated.'
      dims(d) = psi%Basis%tab_basis(d)%nq
      IF (dims(d) < 2 .OR. SIZE(psi%Basis%tab_basis(d)%x) /= dims(d) .OR. &
          SIZE(psi%Basis%tab_basis(d)%w) /= dims(d)) ERROR STOP 'march_SPLIT: invalid Fourier grid size.'
      dx(d) = psi%Basis%tab_basis(d)%x(2)-psi%Basis%tab_basis(d)%x(1)
      IF (dx(d) <= ZERO) ERROR STOP 'march_SPLIT: Fourier grid must be increasing.'
      tol = 100*EPSILON(ONE)*MAX(ONE,ABS(dx(d)))
      IF (ABS(psi%Basis%tab_basis(d)%w(1)-dx(d)) > tol) &
         ERROR STOP 'march_SPLIT: Fourier grid weights must equal the grid spacing.'
      DO iq = 2, dims(d)
         IF (ABS((psi%Basis%tab_basis(d)%x(iq)-psi%Basis%tab_basis(d)%x(iq-1))-dx(d)) > tol) &
            ERROR STOP 'march_SPLIT: nonuniform position grid is not supported.'
         IF (ABS(psi%Basis%tab_basis(d)%w(iq)-dx(d)) > tol) &
            ERROR STOP 'march_SPLIT: Fourier grid weights must equal the grid spacing.'
      END DO
      IF (dims(d)>HUGE(stride)/stride) &
         ERROR STOP 'march_SPLIT: tensor-grid size exceeds the supported integer indexing range.'
      stride = stride*dims(d)
   END DO
   IF (stride /= nq) ERROR STOP 'march_SPLIT: only a complete tensor-product Fourier grid is supported.'
   CALL get_Qmodel_GGdef(GGdef)

   IF (debug) THEN
      WRITE(out_unit,*) 'BEGINNING march_SPLIT (unitary tensor FFT)', &
                        ' order=',propa%split_order,' t=',t,' dt=',propa%delta_t
      WRITE(out_unit,*) 'dimensions=',dims,' surfaces=',nsurf
      FLUSH(out_unit)
   END IF

   CALL init_psi(psi_grid, psi%Basis, cplx=.TRUE., grid=.TRUE.)
   IF (psi%Grid) THEN
      psi_grid%CVec = psi%CVec
   ELSE
      CALL BasisTOGrid_nD_cplx(psi_grid%CVec, psi%CVec, psi%Basis)
   END IF
   CALL ApplySplitStepGrid(psi_grid%CVec,H%Scal_pot,dims,dx,GGdef, &
                           propa%delta_t,nsurf,propa%split_order)

   IF (psi_dt%Grid) THEN
      psi_dt%CVec = psi_grid%CVec
   ELSE
      CALL GridTOBasis_nD_cplx(psi_dt%CVec, psi_grid%CVec, psi%Basis)
   END IF

   IF (debug) THEN
      CALL Calc_Norm_OF_Psi(psi,norm0)
      CALL Calc_Norm_OF_Psi(psi_dt,norm1)
      WRITE(out_unit,*) 'norm before/after=',norm0,norm1,' drift=',ABS(norm0-norm1)
      WRITE(out_unit,*) 'END march_SPLIT'
      FLUSH(out_unit)
   END IF
   CALL dealloc_psi(psi_grid)
   DEALLOCATE(dims,GGdef,dx)
END SUBROUTINE march_SPLIT

SUBROUTINE ApplySplitStepGrid(data, V, dims, dx, GGdef, dt, nsurf, order)
   ! Fourier split propagation with selectable temporal order 1, 2, 4 or 6.
   ! The default remains second-order Strang for backward compatibility.
   IMPLICIT NONE
   COMPLEX(KIND=Rkind), INTENT(INOUT) :: data(:)
   REAL(KIND=Rkind), INTENT(IN) :: V(:,:,:), dx(:), GGdef(:,:), dt
   INTEGER, INTENT(IN) :: dims(:), nsurf
   INTEGER, OPTIONAL, INTENT(IN) :: order
   REAL(KIND=Rkind), ALLOCATABLE :: weights(:)
   INTEGER :: ndim, nq, split_order, stage

   ndim=SIZE(dims)
   IF (ndim<1 .OR. SIZE(dx)/=ndim .OR. &
       SIZE(GGdef,1)/=ndim .OR. SIZE(GGdef,2)/=ndim) &
      ERROR STOP 'ApplySplitStepGrid: invalid spatial dimensions.'
   IF (ANY(dims<2) .OR. nsurf<1) ERROR STOP 'ApplySplitStepGrid: invalid grid or surface dimensions.'
   nq=SplitGridSize(dims)
   IF (nq>HUGE(nq)/nsurf) ERROR STOP 'ApplySplitStepGrid: grid/surface size exceeds integer indexing range.'
   IF (SIZE(data)/=nq*nsurf .OR. SIZE(V,1)/=nq .OR. &
       SIZE(V,2)/=nsurf .OR. SIZE(V,3)/=nsurf) &
      ERROR STOP 'ApplySplitStepGrid: inconsistent grid/potential sizes.'
   split_order=2
   IF (PRESENT(order)) split_order=order

   SELECT CASE(split_order)
   CASE(1)
      CALL ApplySplitFirstOrderStep(data,V,dims,dx,GGdef,dt,nsurf)
   CASE(2)
      CALL ApplySplitStrangStep(data,V,dims,dx,GGdef,dt,nsurf)
   CASE(4,6)
      CALL GetSplitStrangWeights(split_order,weights)
      DO stage=1,SIZE(weights)
         CALL ApplySplitStrangStep(data,V,dims,dx,GGdef,dt*weights(stage),nsurf)
      END DO
      DEALLOCATE(weights)
   CASE DEFAULT
      ERROR STOP 'ApplySplitStepGrid: order must be 1, 2, 4 or 6.'
   END SELECT
END SUBROUTINE ApplySplitStepGrid

SUBROUTINE ApplySplitFirstOrderStep(data,V,dims,dx,GGdef,dt,nsurf)
   IMPLICIT NONE
   COMPLEX(KIND=Rkind), INTENT(INOUT) :: data(:)
   REAL(KIND=Rkind), INTENT(IN) :: V(:,:,:),dx(:),GGdef(:,:),dt
   INTEGER, INTENT(IN) :: dims(:),nsurf
   COMPLEX(KIND=Rkind), ALLOCATABLE :: Ufull(:,:,:)

   ! Lie-Trotter ordering: apply V for dt, then T for dt.
   CALL BuildSplitPotentialHalf(V,TWO*dt,Ufull)
   CALL ApplySplitPotentialMatrix(data,Ufull)
   CALL ApplySplitKineticStep(data,dims,dx,GGdef,dt,nsurf)
   DEALLOCATE(Ufull)
END SUBROUTINE ApplySplitFirstOrderStep

SUBROUTINE ApplySplitStrangStep(data,V,dims,dx,GGdef,dt,nsurf)
   IMPLICIT NONE
   COMPLEX(KIND=Rkind), INTENT(INOUT) :: data(:)
   REAL(KIND=Rkind), INTENT(IN) :: V(:,:,:),dx(:),GGdef(:,:),dt
   INTEGER, INTENT(IN) :: dims(:),nsurf
   COMPLEX(KIND=Rkind), ALLOCATABLE :: Uhalf(:,:,:)

   CALL BuildSplitPotentialHalf(V,dt,Uhalf)
   CALL ApplySplitPotentialMatrix(data,Uhalf)
   CALL ApplySplitKineticStep(data,dims,dx,GGdef,dt,nsurf)
   CALL ApplySplitPotentialMatrix(data,Uhalf)
   DEALLOCATE(Uhalf)
END SUBROUTINE ApplySplitStrangStep

SUBROUTINE ApplySplitKineticStep(data,dims,dx,GGdef,dt,nsurf)
   IMPLICIT NONE
   COMPLEX(KIND=Rkind), INTENT(INOUT) :: data(:)
   REAL(KIND=Rkind), INTENT(IN) :: dx(:),GGdef(:,:),dt
   INTEGER, INTENT(IN) :: dims(:),nsurf
   REAL(KIND=Rkind), ALLOCATABLE :: pvec(:)
   REAL(KIND=Rkind) :: tkin
   INTEGER :: ndim,nq,iq,d,stride,linear,digit,isurf

   ndim=SIZE(dims); nq=SplitGridSize(dims)
   ALLOCATE(pvec(ndim))
   CALL SplitFourierTransform(data,dims,nsurf,.FALSE.)
   DO iq=1,nq
      linear=iq-1; stride=1
      DO d=1,ndim
         digit=MOD(linear/stride,dims(d))
         IF (digit>(dims(d)-1)/2) digit=digit-dims(d)
         pvec(d)=TWO*pi*REAL(digit,Rkind)/(REAL(dims(d),Rkind)*dx(d))
         stride=stride*dims(d)
      END DO
      tkin=HALF*DOT_PRODUCT(pvec,MATMUL(GGdef,pvec))
      DO isurf=1,nsurf
         data(iq+(isurf-1)*nq)=data(iq+(isurf-1)*nq)*EXP(-EYE*tkin*dt)
      END DO
   END DO
   CALL SplitFourierTransform(data,dims,nsurf,.TRUE.)
   DEALLOCATE(pvec)
END SUBROUTINE ApplySplitKineticStep

SUBROUTINE GetSplitStrangWeights(order,weights)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: order
   REAL(KIND=Rkind), ALLOCATABLE, INTENT(OUT) :: weights(:)
   REAL(KIND=Rkind), ALLOCATABLE :: current(:),next(:)
   REAL(KIND=Rkind) :: z1,z0,p
   INTEGER :: current_order,n

   IF (order/=4 .AND. order/=6) ERROR STOP 'GetSplitStrangWeights: supported orders are 4 and 6.'
   ALLOCATE(current(1)); current(1)=ONE; current_order=2
   DO WHILE(current_order<order)
      p=REAL(current_order+1,Rkind)
      z1=ONE/(TWO-TWO**(ONE/p))
      z0=ONE-TWO*z1
      n=SIZE(current); ALLOCATE(next(3*n))
      next(1:n)=z1*current
      next(n+1:2*n)=z0*current
      next(2*n+1:3*n)=z1*current
      CALL MOVE_ALLOC(next,current)
      current_order=current_order+2
   END DO
   CALL MOVE_ALLOC(current,weights)
END SUBROUTINE GetSplitStrangWeights

SUBROUTINE SplitFourierTransform(data,dims,nsurf,inverse)
   ! Separable unitary multidimensional FFT; Bluestein handles arbitrary sizes.
   IMPLICIT NONE
   COMPLEX(KIND=Rkind), TARGET, INTENT(INOUT) :: data(:)
   INTEGER, INTENT(IN) :: dims(:),nsurf
   LOGICAL, INTENT(IN) :: inverse
   COMPLEX(KIND=Rkind), POINTER :: grid(:,:)
   COMPLEX(KIND=Rkind), ALLOCATABLE :: work(:,:),next(:,:),line(:)
   INTEGER :: ndim,nq,axis,n,stride,block,nblocks
   INTEGER :: iblock,inner,j,idx,isurf

   ndim=SIZE(dims)
   IF (ndim<1) ERROR STOP 'SplitFourierTransform: at least one dimension is required.'
   IF (ANY(dims<2) .OR. nsurf<1) ERROR STOP 'SplitFourierTransform: invalid dimensions.'
   nq=SplitGridSize(dims)
   IF (nq>HUGE(nq)/nsurf) ERROR STOP 'SplitFourierTransform: grid/surface size exceeds integer indexing range.'
   IF (SIZE(data)/=nq*nsurf) ERROR STOP 'SplitFourierTransform: array size mismatch.'
   grid(1:nq,1:nsurf)=>data
   ALLOCATE(work(nq,nsurf),next(nq,nsurf)); work=grid
   stride=1
   DO axis=1,ndim
      n=dims(axis); block=stride*n; nblocks=nq/block
      ALLOCATE(line(n)); next=work
      DO iblock=0,nblocks-1
         DO inner=0,stride-1
            DO isurf=1,nsurf
               DO j=0,n-1
                  idx=iblock*block+inner+j*stride+1
                  line(j+1)=work(idx,isurf)
               END DO
               CALL SplitFFT1D(line,inverse)
               DO j=0,n-1
                  idx=iblock*block+inner+j*stride+1
                  next(idx,isurf)=line(j+1)
               END DO
            END DO
         END DO
      END DO
      work=next; stride=block; DEALLOCATE(line)
   END DO
   grid=work; NULLIFY(grid); DEALLOCATE(work,next)
END SUBROUTINE SplitFourierTransform

INTEGER FUNCTION SplitGridSize(dims) RESULT(nq)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: dims(:)
   INTEGER :: d

   IF (SIZE(dims)<1 .OR. ANY(dims<1)) ERROR STOP 'SplitGridSize: invalid tensor dimensions.'
   nq=1
   DO d=1,SIZE(dims)
      IF (dims(d)>HUGE(nq)/nq) &
         ERROR STOP 'SplitGridSize: tensor-grid size exceeds the supported integer indexing range.'
      nq=nq*dims(d)
   END DO
END FUNCTION SplitGridSize

SUBROUTINE SplitFFT1D(data,inverse)
   IMPLICIT NONE
   COMPLEX(KIND=Rkind), INTENT(INOUT) :: data(:)
   LOGICAL, INTENT(IN) :: inverse
   COMPLEX(KIND=Rkind), ALLOCATABLE :: a(:),b(:)
   COMPLEX(KIND=Rkind) :: phase
   REAL(KIND=Rkind) :: direction,angle,norm
   INTEGER :: n,m,j

   n=SIZE(data)
   IF (n<1) ERROR STOP 'SplitFFT1D: empty transform.'
   direction=-ONE
   IF (inverse) direction=ONE
   norm=ONE/SQRT(REAL(n,Rkind))
   IF (IAND(n,n-1)==0) THEN
      CALL SplitFFTRadix2(data,inverse)
      data=data*norm
      RETURN
   END IF

   ! Bluestein chirp-z convolution gives an O(M log M) transform for every n;
   ! M is the next power of two >= 2*n-1.
   IF (n>(HUGE(n)-1)/2) ERROR STOP 'SplitFFT1D: transform too large.'
   m=1
   DO WHILE(m<2*n-1)
      IF (m>HUGE(m)/2) ERROR STOP 'SplitFFT1D: convolution workspace overflow.'
      m=2*m
   END DO
   ALLOCATE(a(m),b(m)); a=CZERO; b=CZERO
   DO j=0,n-1
      angle=pi*REAL(j,Rkind)*REAL(j,Rkind)/REAL(n,Rkind)
      phase=CMPLX(COS(direction*angle),SIN(direction*angle),KIND=Rkind)
      a(j+1)=data(j+1)*phase
      phase=CONJG(phase)
      b(j+1)=phase
      IF (j>0) b(m-j+1)=phase
   END DO
   CALL SplitFFTRadix2(a,.FALSE.)
   CALL SplitFFTRadix2(b,.FALSE.)
   a=a*b
   CALL SplitFFTRadix2(a,.TRUE.)
   DO j=0,n-1
      angle=pi*REAL(j,Rkind)*REAL(j,Rkind)/REAL(n,Rkind)
      phase=CMPLX(COS(direction*angle),SIN(direction*angle),KIND=Rkind)
      data(j+1)=a(j+1)*phase*(norm/REAL(m,Rkind))
   END DO
   DEALLOCATE(a,b)
END SUBROUTINE SplitFFT1D

SUBROUTINE SplitFFTRadix2(data,inverse)
   IMPLICIT NONE
   COMPLEX(KIND=Rkind), INTENT(INOUT) :: data(:)
   LOGICAL, INTENT(IN) :: inverse
   COMPLEX(KIND=Rkind) :: temp,w,wm,u,v
   REAL(KIND=Rkind) :: direction,angle
   INTEGER :: n,i,j,bit,len,half,start,k

   n=SIZE(data)
   IF (n<1 .OR. IAND(n,n-1)/=0) ERROR STOP 'SplitFFTRadix2: length must be a power of two.'
   j=1
   DO i=1,n-1
      IF (i<j) THEN
         temp=data(i); data(i)=data(j); data(j)=temp
      END IF
      bit=n/2
      DO WHILE(bit>=1 .AND. j>bit)
         j=j-bit; bit=bit/2
      END DO
      j=j+bit
   END DO
   direction=-ONE
   IF (inverse) direction=ONE
   len=2
   DO WHILE(len<=n)
      half=len/2; angle=direction*TWO*pi/REAL(len,Rkind)
      wm=CMPLX(COS(angle),SIN(angle),KIND=Rkind)
      DO start=1,n,len
         w=CMPLX(ONE,ZERO,KIND=Rkind)
         DO k=0,half-1
            u=data(start+k); v=w*data(start+k+half)
            data(start+k)=u+v; data(start+k+half)=u-v
            w=w*wm
         END DO
      END DO
      IF (len==n) EXIT
      len=2*len
   END DO
END SUBROUTINE SplitFFTRadix2

SUBROUTINE BuildSplitPotentialHalf(V, dt, Uhalf)
   ! Precompute exp(-i V(q) dt/2) once per grid point; the same matrices
   ! are then reused for both halves of this Strang step.
   IMPLICIT NONE
   REAL(KIND=Rkind), INTENT(IN) :: V(:,:,:), dt
   COMPLEX(KIND=Rkind), ALLOCATABLE, INTENT(OUT) :: Uhalf(:,:,:)
   REAL(KIND=Rkind), ALLOCATABLE :: vmat(:,:), eval(:), evec(:,:)
   INTEGER :: nq, nsurf, iq, i, j

   nq=SIZE(V,1); nsurf=SIZE(V,2)
   IF (nsurf<1 .OR. SIZE(V,3)/=nsurf) &
      ERROR STOP 'BuildSplitPotentialHalf: invalid potential dimensions.'
   ALLOCATE(Uhalf(nq,nsurf,nsurf))
   IF (nsurf==1) THEN
      DO iq=1,nq
         Uhalf(iq,1,1)=EXP(-EYE*HALF*dt*V(iq,1,1))
      END DO
      RETURN
   END IF

   ALLOCATE(vmat(nsurf,nsurf),eval(nsurf),evec(nsurf,nsurf))
   DO iq=1,nq
      vmat(:,:)=V(iq,:,:)
      IF (MAXVAL(ABS(vmat-TRANSPOSE(vmat))) > &
          100*EPSILON(ONE)*MAX(ONE,MAXVAL(ABS(vmat)))) &
         ERROR STOP 'BuildSplitPotentialHalf: the electronic potential must be real symmetric.'
      CALL diagonalization(vmat,eval,evec,sort=1)
      DO i=1,nsurf
         DO j=1,nsurf
            Uhalf(iq,i,j)=SUM(evec(i,:)*EXP(-EYE*HALF*dt*eval(:))*evec(j,:))
         END DO
      END DO
   END DO
   DEALLOCATE(vmat,eval,evec)
END SUBROUTINE BuildSplitPotentialHalf

SUBROUTINE ApplySplitPotentialMatrix(data,Uhalf)
   IMPLICIT NONE
   COMPLEX(KIND=Rkind), TARGET, INTENT(INOUT) :: data(:)
   COMPLEX(KIND=Rkind), INTENT(IN) :: Uhalf(:,:,:)
   COMPLEX(KIND=Rkind), POINTER :: grid(:,:)
   COMPLEX(KIND=Rkind), ALLOCATABLE :: ctmp(:)
   INTEGER :: nq,nsurf,iq

   nq=SIZE(Uhalf,1); nsurf=SIZE(Uhalf,2)
   IF (SIZE(Uhalf,3)/=nsurf .OR. SIZE(data)/=nq*nsurf) &
      ERROR STOP 'ApplySplitPotentialMatrix: inconsistent dimensions.'
   grid(1:nq,1:nsurf)=>data
   ALLOCATE(ctmp(nsurf))
   DO iq=1,nq
      ctmp=MATMUL(Uhalf(iq,:,:),grid(iq,:))
      grid(iq,:)=ctmp
   END DO
   DEALLOCATE(ctmp)
   NULLIFY(grid)
END SUBROUTINE ApplySplitPotentialMatrix


 !--------------------------------------------------------------
! march_VP
! Description:
!   - Time propagation of the variational parameters Λ
!     using the implicit finite-difference scheme:
!
!       i M Λ̇ = F   →   Λ(t+dt) = Λ(t) - i dt (M^(-1)) F
!
! Inputs:
!   psi   : Current wavefunction structure
!   propa : Propagator parameters (contains delta_t, etc.)
!   H     : Hamiltonian operator
!   t     : Current time
!
! Outputs:
!   psi_dt : Updated wavefunction after one step
!--------------------------------------------------------------
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
   COMPLEX(KIND=Rkind), ALLOCATABLE   :: M(:,:), F(:)
   REAL(KIND=Rkind)                   :: Norm0, Norm
   INTEGER :: info, n
 
   LOGICAL :: debug
 

   CALL init_psi(tmp, psi%Basis, cplx=.TRUE., grid=.FALSE.)
   tmp%CVec = psi%CVec
   !------------------------------------------------------------
   ! Step 1. Convert psi -> Lambda
   !------------------------------------------------------------
   CALL Psi_Lambda_Transform(tmp, Lambda, .TRUE.)
   n = SIZE(Lambda)
   CALL Calc_Norm_OF_Psi(psi, Norm0)
   !------------------------------------------------------------
   ! Step 2. Build M and F at time t
   !------------------------------------------------------------
   CALL Calc_M_nD(psi=psi, M_nD=M)   ! generalized metric matrix
   CALL Calc_F_nD(psi=psi, H=H, F_nD=F)! force-like vector
 
   !------------------------------------------------------------
   ! Step 3. Solve ΔΛ = -i dt M^(-1) F
   !------------------------------------------------------------
   ALLOCATE(Lambda_new(n))
   Lambda_new = Lambda
 
   M = inv_OF_Mat(M,inv_type=0,epsi=1.d-20)   ! M = (M)^(-1)
 
   ! Update Lambda
   Lambda_new = Lambda -  EYE* propa%delta_t * MATMUL(M,F)
 
   !------------------------------------------------------------
   ! Step 4. Transform back: Lambda -> psi_dt
   !------------------------------------------------------------
   CALL Psi_Lambda_Transform(psi_dt, Lambda_new, .FALSE.)
 
   !------------------------------------------------------------
   ! Debug info
   !------------------------------------------------------------
   debug = propa%debug
   IF (debug) THEN
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
   DEALLOCATE(Lambda, Lambda_new, M, F)
   CALL dealloc_psi(tmp)
 
 END SUBROUTINE march_VP


 !--------------------------------------------------------------
! march_VP_RK4
! Description:
!   - Time propagation of the variational parameters Λ
!     using the Runge-Kutta 4th order method:
!
!       i M Λ̇ = F   →   Λ̇ = -i M^(-1) F
!
! Inputs:
!   psi   : Current wavefunction structure
!   propa : Propagator parameters (contains delta_t, etc.)
!   H     : Hamiltonian operator
!   t     : Current time
!
! Outputs:
!   psi_dt : Updated wavefunction after one step
!--------------------------------------------------------------
SUBROUTINE march_VP_RK4(psi, psi_dt, t, propa, H)
   USE QDUtil_m
   IMPLICIT NONE
   TYPE(psi_t),    INTENT(IN)    :: psi
   TYPE(psi_t),    INTENT(INOUT) :: psi_dt
   TYPE(propa_t),  INTENT(IN)    :: propa
   TYPE(Op_t),     INTENT(IN)    :: H
   REAL(KIND=Rkind), INTENT(IN)  :: t
 
   ! local variables
   COMPLEX(KIND=Rkind), ALLOCATABLE :: Lambda(:), Lambda_new(:)
   COMPLEX(KIND=Rkind), ALLOCATABLE :: k1(:), k2(:), k3(:), k4(:)
   COMPLEX(KIND=Rkind), ALLOCATABLE :: Lambda_temp(:)
   TYPE(psi_t)                      :: tmp_psi
   TYPE(Op_t)                       :: tmp_op
   COMPLEX(KIND=Rkind), ALLOCATABLE :: M(:,:), F(:)
   REAL(KIND=Rkind)                 :: Norm0, Norm, half_dt
   INTEGER, ALLOCATABLE             :: Tab_Iq(:,:)
   INTEGER :: nC,i
   REAL(KIND=Rkind)                 :: scale_factor
 
   LOGICAL :: debug
 
   half_dt = propa%delta_t * HALF
 
   !------------------------------------------------------------
   ! Step 1. Convert psi -> Lambda
   !------------------------------------------------------------
   CALL init_psi(tmp_psi, psi%Basis, cplx=.TRUE., grid=.FALSE.)
   CALL Set_Op(tmp_op, psi%Basis, Tab_Iq,alloc=.FALSE.,nio=in_unit)
   tmp_psi%CVec = psi%CVec
   CALL Psi_Lambda_Transform(tmp_psi, Lambda, .TRUE.)
   nC = SIZE(Lambda)
   debug = propa%debug
   IF(debug) THEN
      CALL Calc_Norm_OF_Psi(psi, Norm0)
      DO i = 1, nC
         WRITE(out_unit, *) i, Lambda(i)
      END DO
   END IF
   ! Allocate arrays for RK4
   ALLOCATE(Lambda_new(nC), k1(nC), k2(nC), k3(nC), k4(nC))
   ALLOCATE(Lambda_temp(nC))
 
   !------------------------------------------------------------
   ! Step 2. Calculate k1 = f(t, Lambda)
   !------------------------------------------------------------
   CALL Calc_M_nD(psi=psi, M_nD=M)
   CALL Calc_F_nD(psi=psi, H=tmp_op, F_nD=F)
   scale_factor = ONE / MAX(MAXVAL(ABS(M)), MAXVAL(ABS(F)), 1.d-20)
   
   ! Solve M * k1 = F → k1 = M^(-1) F
   k1 = -EYE * scale_factor * MATMUL(inv_OF_Mat(M,inv_type=0,epsi=1.d-20) , F)  ! f(t, Lambda) = -i M^(-1) F
 
   !------------------------------------------------------------
   ! Step 3. Calculate k2 = f(t + dt/2, Lambda + dt/2 * k1)
   !------------------------------------------------------------
   Lambda_temp = Lambda + half_dt * k1
   
   ! Create temporary wavefunction for intermediate step
   CALL init_psi(tmp_psi, psi%Basis, cplx=.TRUE., grid=.FALSE.)
   CALL Psi_Lambda_Transform(tmp_psi, Lambda_temp, .FALSE.)
   
   CALL Set_Op(tmp_op, tmp_psi%Basis, Tab_Iq,alloc=.FALSE.,nio=in_unit)
   CALL Calc_M_nD(psi=tmp_psi, M_nD=M)
   CALL Calc_F_nD(psi=tmp_psi, H=tmp_op, F_nD=F)

   ! Solve M * k2 = F → k2 = M^(-1) F
   k2 = -EYE * MATMUL(inv_OF_Mat(M,inv_type=0,epsi=1.d-20) , F)  ! f(t, Lambda) = -i M^(-1) F
 
   !------------------------------------------------------------
   ! Step 4. Calculate k3 = f(t + dt/2, Lambda + dt/2 * k2)
   !------------------------------------------------------------
   Lambda_temp = Lambda + half_dt * k2
   
   CALL Psi_Lambda_Transform(tmp_psi, Lambda_temp, .FALSE.)
   
   CALL Set_Op(tmp_op, tmp_psi%Basis, Tab_Iq,alloc=.FALSE.,nio=in_unit)
   CALL Calc_M_nD(psi=tmp_psi, M_nD=M)
   CALL Calc_F_nD(psi=tmp_psi, H=tmp_op, F_nD=F)
   
   ! Solve M * k3 = F → k3 = M^(-1) F
   k3 = -EYE * MATMUL(inv_OF_Mat(M,inv_type=0,epsi=1.d-20) , F) ! f(t, Lambda) = -i M^(-1) F
 
   !------------------------------------------------------------
   ! Step 5. Calculate k4 = f(t + dt, Lambda + dt * k3)
   !------------------------------------------------------------
   Lambda_temp = Lambda + propa%delta_t * k3
   
   CALL Psi_Lambda_Transform(tmp_psi, Lambda_temp, .FALSE.)
   
   CALL Set_Op(tmp_op, tmp_psi%Basis, Tab_Iq,alloc=.FALSE.,nio=in_unit)
   CALL Calc_M_nD(psi=tmp_psi, M_nD=M)
   CALL Calc_F_nD(psi=tmp_psi, H=tmp_op, F_nD=F)
   
   ! Solve M * k4 = F → k4 = M^(-1) F
   k4 = -EYE *MATMUL(inv_OF_Mat(M,inv_type=0,epsi=1.d-20) , F)  ! f(t, Lambda) = -i M^(-1) F
 
   !------------------------------------------------------------
   ! Step 6. Combine all slopes for final update
   !------------------------------------------------------------
   Lambda_new = Lambda + propa%delta_t * SIXTH*(k1 + TWO*k2 + TWO*k3 + k4)
 
   !------------------------------------------------------------
   ! Step 7. Transform back: Lambda_new -> psi_dt
   !------------------------------------------------------------
   CALL Psi_Lambda_Transform(psi_dt, Lambda_new, .FALSE.)
 
   !------------------------------------------------------------
   ! Debug info
   !------------------------------------------------------------
   IF (debug) THEN

      CALL Write_Vec(psi_dt%CVec, out_unit, 1,info='psi_dt')
      WRITE(out_unit,*) "BEGIN march_VP_RK4 at t=", t, " dt=", propa%delta_t
      CALL Calc_Norm_OF_Psi(psi_dt, Norm)
      WRITE(out_unit, *) '<psi_dt|psi_dt> = ', Norm, ' |ΔNorm| = ', abs(Norm0 - Norm)
      WRITE(out_unit,*) "END march_VP_RK4"
      FLUSH(out_unit)
   END IF
   !------------------------------------------------------------
   ! Clean up
   !------------------------------------------------------------
   CALL dealloc_psi(tmp_psi)
   CALL dealloc_op(tmp_op) 
   DEALLOCATE(Lambda, Lambda_new, k1, k2, k3, k4, Lambda_temp, M, F)
 
 END SUBROUTINE march_VP_RK4



!=======================================================================
! SUBROUTINE: march_VPTMP
!
! PURPOSE:
!   Time propagation of the wavefunction ψ(t) using the implicit
!   finite-difference scheme in the Parametric Variational Dynamics
!   with Trajectories (PVDT).
!
!   Governing equation:
!       i M dΛ/dt = F,
!   where the block structure of M allows for efficient solution
!   via the Schur complement:
!
!       S = M_CC - M_Cλ M_λλ^{-1} M_λC
!
!   Algorithm:
!     1. Build matrices (S, M_Cλ, M_λλ) and forces (F_C, F_λ)
!     2. Invert M_λλ with inv_OF_Mat
!     3. Solve reduced system S ΔΛ_C = F_C - M_Cλ M_λλ^{-1} F_λ
!     4. Back-substitute to recover ΔΛ_λ
!     5. Update ψ(t) → ψ(t+Δt) using (ΔΛ_C, ΔΛ_λ)
!
! INPUT:
!   psi     : Current wavefunction object
!   psi_dt  : Updated wavefunction object (output)
!   t       : Current time
!   propa   : Propagation parameters (contains time step Δt, etc.)
!   H       : Hamiltonian operator
!
! OUTPUT:
!   psi_dt  : Updated wavefunction after one propagation step
!
! DEPENDENCIES:
!   - Calc_S_MClambda_MLambdaLambda: builds S, M_Cλ, M_λλ
!   - Calc_F_nD: builds F_C and F_λ
!   - inv_OF_Mat: inverts a square matrix with stabilization
!
!=======================================================================

 SUBROUTINE march_VPTMP(psi, psi_dt, t, propa, H)
   USE QDUtil_m
   IMPLICIT NONE

   !--------------------------------------------------------------------
   ! Arguments
   !--------------------------------------------------------------------
   TYPE(psi_t),    INTENT(IN)    :: psi
   TYPE(psi_t),    INTENT(INOUT) :: psi_dt
   TYPE(propa_t),  INTENT(IN)    :: propa
   TYPE(Op_t),     INTENT(IN)    :: H
   REAL(KIND=Rkind), INTENT(IN)  :: t

   !--------------------------------------------------------------------
   ! Local variables
   !--------------------------------------------------------------------
   COMPLEX(KIND=Rkind), ALLOCATABLE :: FC(:), FL(:)
   COMPLEX(KIND=Rkind), ALLOCATABLE :: M_CC(:,:), M_CL(:,:), M_LL(:,:)
   COMPLEX(KIND=Rkind), ALLOCATABLE :: D_L(:), D_C(:)
   COMPLEX(KIND=Rkind), ALLOCATABLE :: Schur_M(:,:), temp_M(:,:), temp_V(:)
   REAL(KIND=Rkind)                 :: Norm0, Norm
   INTEGER                          :: nC, nL, i, j, alloc_stat
   LOGICAL :: debug
   CHARACTER(LEN=:), ALLOCATABLE    :: error_msg

   !--------------------------------------------------------------------
   ! 0. Initial norm of ψ(t)
   !--------------------------------------------------------------------
   debug = propa%debug
   IF (debug) THEN
      CALL Calc_Norm_OF_Psi(psi, Norm0)
      WRITE(out_unit,*) 'Norm of initial psi in March_VPTMP: ', Norm0
   END IF

   !--------------------------------------------------------------------
   ! 1. Build matrices and force vectors
   !--------------------------------------------------------------------
   CALL Calc_S_MClambda_MLambdaLambda(psi, M_CC, M_CL, M_LL)
   CALL Calc_F_nD(psi=psi, H=H, FC_nD=FC, F_L=FL)

   nC = SIZE(FC)
   nL = SIZE(FL)

   ! Vérification des dimensions
   IF (SIZE(M_CC,1) /= nC .OR. SIZE(M_CC,2) /= nC .OR. &
       SIZE(M_CL,1) /= nC .OR. SIZE(M_CL,2) /= nL .OR. &
       SIZE(M_LL,1) /= nL .OR. SIZE(M_LL,2) /= nL) THEN
      WRITE(out_unit,*) 'ERROR: Dimension mismatch in march_VPTMP'
      STOP
   END IF

   !--------------------------------------------------------------------
   ! 2. Invert M_λλ with error handling
   !--------------------------------------------------------------------
   ALLOCATE(temp_M(nL,nL), stat=alloc_stat)
   IF (alloc_stat /= 0) THEN
      error_msg = 'Allocation failed for temp_M in march_VPTMP'
      WRITE(out_unit,*) error_msg
      STOP
   END IF
   temp_M = M_LL
   M_LL=inv_OF_Mat(temp_M, inv_type=0, epsi=1.d-20)
   DEALLOCATE(temp_M, stat=alloc_stat)
   IF (alloc_stat /= 0) THEN
      WRITE(out_unit,*) 'Warning: Deallocation failed for temp_M'
   END IF

   !--------------------------------------------------------------------
   ! 3. Compute Schur Matrix with safe operations
   !--------------------------------------------------------------------
   ALLOCATE(Schur_M(nC,nC), stat=alloc_stat)
   IF (alloc_stat /= 0) THEN
      error_msg = 'Allocation failed for Schur_M in march_VPTMP'
      WRITE(out_unit,*) error_msg
      STOP
   END IF
   
   ALLOCATE(temp_M(nC,nL), stat=alloc_stat)
   IF (alloc_stat /= 0) THEN
      error_msg = 'Allocation failed for temp_M in march_VPTMP'
      WRITE(out_unit,*) error_msg
      STOP
   END IF
   
   ! M_CL * M_LL
   temp_M = MATMUL(M_CL, M_LL)
   
   ! M_CL * M_LL * TRANSPOSE(M_CL)
   Schur_M = MATMUL(temp_M, CONJG(TRANSPOSE(M_CL)))
   
   ! Schur complement: M_CC - M_CL * M_LL^{-1} * M_CL^T
   Schur_M = M_CC - Schur_M
   
   DEALLOCATE(temp_M, stat=alloc_stat)
   IF (alloc_stat /= 0) THEN
      WRITE(out_unit,*) 'Warning: Deallocation failed for temp_M'
   END IF

   !--------------------------------------------------------------------
   ! 4. Solve reduced system with Schur complement
   !--------------------------------------------------------------------
   ALLOCATE(D_C(nC), stat=alloc_stat)
   IF (alloc_stat /= 0) THEN
      error_msg = 'Allocation failed for D_C in march_VPTMP'
      WRITE(out_unit,*) error_msg
      STOP
   END IF
   
   ALLOCATE(temp_V(nL), stat=alloc_stat)
   IF (alloc_stat /= 0) THEN
      error_msg = 'Allocation failed for temp_V in march_VPTMP'
      WRITE(out_unit,*) error_msg
      STOP
   END IF
   
   ! temp_V = M_LL * FL
   temp_V = MATMUL(M_LL, FL)
   
   ! D_C = FC - M_CL * temp_V
   D_C = FC - MATMUL(M_CL, temp_V)
   
   ! Invert Schur matrix
   ALLOCATE(temp_M(nC,nC), stat=alloc_stat)
   IF (alloc_stat /= 0) THEN
      error_msg = 'Allocation failed for temp_M in Schur inversion'
      WRITE(out_unit,*) error_msg
      STOP
   END IF
   temp_M = Schur_M
   Schur_M = inv_OF_Mat(temp_M, inv_type=0, epsi=1.d-20)

   DEALLOCATE(temp_M, stat=alloc_stat)
   IF (alloc_stat /= 0) THEN
      WRITE(out_unit,*) 'Warning: Deallocation failed for temp_M'
   END IF
   
   ! D_C = Schur_M^{-1} * D_C
   D_C = MATMUL(Schur_M, D_C)

   !--------------------------------------------------------------------
   ! 5. Compute D_L
   !--------------------------------------------------------------------
   ALLOCATE(D_L(nL), stat=alloc_stat)
   IF (alloc_stat /= 0) THEN
      error_msg = 'Allocation failed for D_L in march_VPTMP'
      WRITE(out_unit,*) error_msg
      STOP
   END IF
   
   ! D_L = FL - M_CL^T * D_C
   D_L = FL - MATMUL(CONJG(TRANSPOSE(M_CL)), D_C)
   
   ! D_L = M_LL * D_L
   D_L = MATMUL(M_LL, D_L)

   DEALLOCATE(temp_V, stat=alloc_stat)
   IF (alloc_stat /= 0) THEN
      WRITE(out_unit,*) 'Warning: Deallocation failed for temp_V'
   END IF

   !--------------------------------------------------------------------
   ! 6. Update wavefunction
   !--------------------------------------------------------------------
   CALL Update_Psi_From_DeltaLambda(psi, psi_dt, D_C, D_L, propa%delta_t)

   !--------------------------------------------------------------------
   ! 7. Debug output
   !--------------------------------------------------------------------
   IF (debug) THEN
      WRITE(out_unit,*) "============================================="
      WRITE(out_unit,*) "BEGIN march_VP at t=", t, " dt=", propa%delta_t
      CALL Calc_Norm_OF_Psi(psi_dt, Norm)
      WRITE(out_unit, *) "<psi_dt|psi_dt> = ", Norm, " |ΔNorm| = ", ABS(Norm0 - Norm)
      WRITE(out_unit,*) "END march_VP"
      WRITE(out_unit,*) "============================================="
      FLUSH(out_unit)
   END IF

   !--------------------------------------------------------------------
   ! 8. Cleanup avec vérifications complètes
   !--------------------------------------------------------------------
   CALL safe_deallocate_C1D(D_C, 'D_C')
   CALL safe_deallocate_C1D(D_L, 'D_L')
   CALL safe_deallocate_C1D(FC, 'FC')
   CALL safe_deallocate_C1D(FL, 'FL')
   CALL safe_deallocate_C2D(Schur_M, 'Schur_M')
   CALL safe_deallocate_C2D(M_CC, 'M_CC')
   CALL safe_deallocate_C2D(M_CL, 'M_CL')
   CALL safe_deallocate_C2D(M_LL, 'M_LL')

CONTAINS

   ! Sous-routines internes pour la désallocation sécurisée
   SUBROUTINE safe_deallocate_C1D(array, array_name)
      COMPLEX(KIND=Rkind), ALLOCATABLE, INTENT(INOUT) :: array(:)
      CHARACTER(LEN=*), INTENT(IN) :: array_name
      INTEGER :: stat
      
      IF (ALLOCATED(array)) THEN
         DEALLOCATE(array, stat=stat)
         IF (stat /= 0) THEN
            WRITE(out_unit,*) 'Warning: Failed to deallocate ', TRIM(array_name)
         END IF
      END IF
   END SUBROUTINE safe_deallocate_C1D

   SUBROUTINE safe_deallocate_C2D(array, array_name)
      COMPLEX(KIND=Rkind), ALLOCATABLE, INTENT(INOUT) :: array(:,:)
      CHARACTER(LEN=*), INTENT(IN) :: array_name
      INTEGER :: stat
      
      IF (ALLOCATED(array)) THEN
         DEALLOCATE(array, stat=stat)
         IF (stat /= 0) THEN
            WRITE(out_unit,*) 'Warning: Failed to deallocate ', TRIM(array_name)
         END IF
      END IF
   END SUBROUTINE safe_deallocate_C2D

END SUBROUTINE march_VPTMP

!----------------------------------------------------------------------------------
!> SUBROUTINE: march_CN
!> Description: Wrapper for time propagation using the Crank–Nicolson (CN) method.
!>
!>              This routine builds the Crank–Nicolson propagator U and applies it
!>              to the wavefunction. 
!>
!> Arguments:
!>   - psi       : (IN)     Wavefunction at time t (unchanged).
!>   - psi_dt    : (INOUT)  Wavefunction at time t + delta_t (modified).
!>   - t         : (IN)     Current time.
!>   - propa     : (IN)     Propagation parameters (delta_t, ...).
!>   - H         : (IN)     Hamiltonian operator.
!----------------------------------------------------------------------------------
SUBROUTINE march_CN(psi, psi_dt, t, propa, H)
   IMPLICIT NONE

   !-------------- Arguments ------------------------------------------------
   TYPE(psi_t),   INTENT(IN)     :: psi
   TYPE(psi_t),   INTENT(INOUT)  :: psi_dt
   REAL(KIND=Rkind), INTENT(IN)  :: t
   TYPE(propa_t), INTENT(IN)     :: propa
   TYPE(Op_t),    INTENT(IN)     :: H

   !-------------- Local variables --------------------------------------
   REAL(KIND=Rkind)                   :: Norm, Norm0
   COMPLEX(KIND=Rkind), ALLOCATABLE   :: U(:,:)
   LOGICAL :: debug
   INTEGER, ALLOCATABLE               :: Tab_Iq(:, :) 
   TYPE(Op_t)                         :: Op
   ! --- Begin propagation ---------------------------------------------
   debug = propa%debug
   IF (debug) THEN
      WRITE(out_unit,*) 'BEGINNING march_CN ', t,'Δt=', propa%delta_t
      FLUSH(out_unit)
   END IF

   !initialize  Op
   CALL set_op(Op, psi%basis,Tab_Iq,alloc=.TRUE., nio=in_unit)
   CALL Build_Hamiltonian_Matrix(Op,debug)

   ! Build Crank–Nicolson propagator U
   CALL Build_CN_Propagator(Op%CMat, propa%delta_t, U)

   ! Apply propagation: psi_dt = U * psi
   psi_dt%CVec(:) = MATMUL(U, psi%CVec)

   ! --- Debug: check norm preservation --------------------------------
   IF (debug) THEN
      CALL Calc_Norm_OF_Psi(psi, Norm0)
      CALL Calc_Norm_OF_Psi(psi_dt, Norm)
      WRITE(out_unit,*) '<psi_dt|psi_dt> = ', Norm,'  |ΔNorm|= ', ABS(Norm0 - Norm)
      WRITE(out_unit,*) 'END march_CN'
      FLUSH(out_unit)
   END IF

END SUBROUTINE march_CN

!----------------------------------------------------------------------------------
!> SUBROUTINE: march_CN_With_GMRES
!>
!> Purpose:
!>   Wrapper for time propagation using the Crank–Nicolson (CN) scheme combined
!>   with a GMRES iterative solver. This routine constructs the CN propagator
!>   implicitly and applies it to the wavefunction without forming the full matrix.
!>
!> Description:
!>   Given an input wavefunction |psi(t)>, this subroutine computes the propagated
!>   wavefunction |psi(t+Δt)> using the Crank–Nicolson method. The linear system
!>   associated with CN is solved using GMRES, making the method suitable for
!>   large-scale Hamiltonians where explicit matrices are not available.
!>
!> Arguments:
!>   - psi     : (IN)     Wavefunction at time t (not modified).
!>   - psi_dt  : (INOUT)  Wavefunction at time t + Δt (overwritten on exit).
!>   - t       : (IN)     Current simulation time.
!>   - propa   : (IN)     Propagation parameters (delta_t, max_iter, eps, etc.).
!>   - H       : (IN)     Hamiltonian operator.
!>
!> Notes:
!>   - The routine checks the norm preservation of the propagated wavefunction.
!>   - Debugging mtiseages are printed if the local variable 'debug' is set to .TRUE.
!----------------------------------------------------------------------------------

!============================================================
! SUBROUTINE: solve_tdse_cn_gmres
! Purpose: propagate TDSE one time step using CN+GMRES
!============================================================
SUBROUTINE solve_tdse_cn_gmres(psi_n, H_op, psi_np1, dt, max_gmres, gmres_tol, debug_enabled)
   IMPLICIT NONE
   TYPE(psi_t),  INTENT(IN)  :: psi_n
   TYPE(op_t),   INTENT(IN)  :: H_op
   TYPE(psi_t),  INTENT(INOUT) :: psi_np1
   REAL(KIND=Rkind),INTENT(IN)  :: dt
   INTEGER,      INTENT(IN)  :: max_gmres
   REAL(KIND=Rkind),INTENT(IN)  :: gmres_tol
   LOGICAL, INTENT(IN) :: debug_enabled
 
   TYPE(psi_t) :: b
   COMPLEX(KIND=Rkind) :: cn_fac
 
   ! Compute CN factor: i*dt/2
   cn_fac = CMPLX(0.0_Rkind, 0.5_Rkind*dt, kind=Rkind)
 
   ! Initialize RHS b = (I - cn_fac*H)*psi_n
   CALL init_psi(b, psi_n%Basis, cplx=.TRUE., grid=.FALSE.)
   CALL matvec_cn_op(H_op, psi_n, b, -cn_fac)  ! b = psi_n - cn_fac*H*psi_n
 
   ! Initialize gutise psi_np1 = 0
   CALL init_psi(psi_np1, psi_n%Basis, cplx=.TRUE., grid=.FALSE.)
   psi_np1%CVec = CZERO
 
   ! Solve (I + cn_fac*H) psi_np1 = b using GMRES
   CALL gmres_solver(H_op, b, psi_np1, max_gmres, gmres_tol, cn_fac, debug_enabled)
 
   ! Clean up
   CALL dealloc_psi(b)
 END SUBROUTINE solve_tdse_cn_gmres
 
 !============================================================
 ! SUBROUTINE: march_CN_With_GMRES
 ! Wrapper: propagate one CN step
 !============================================================
 SUBROUTINE march_CN_With_GMRES(psi, psi_dt, t, propa, H)
   IMPLICIT NONE
   TYPE(psi_t),   INTENT(IN)     :: psi
   TYPE(psi_t),   INTENT(INOUT)  :: psi_dt
   REAL(KIND=Rkind), INTENT(IN)  :: t
   TYPE(propa_t), INTENT(IN)     :: propa
   TYPE(Op_t),    INTENT(IN)     :: H
 
   REAL(KIND=Rkind) :: Norm0, Norm
   LOGICAL :: debug

   debug = propa%debug
 
   IF(debug) THEN
      WRITE(out_unit,*) '>>> BEGIN march_CN_With_GMRES at time=', t
      IF(.NOT. ASSOCIATED(psi%Basis)) STOP 'psi%Basis not associated'
      IF(.NOT. ASSOCIATED(psi_dt%Basis)) STOP 'psi_dt%Basis not associated'
   END IF
 
   CALL solve_tdse_cn_gmres(psi, H, psi_dt, propa%delta_t, propa%max_iter, propa%eps, propa%debug)
 
   ! Norm check
   IF(debug) THEN
      CALL Calc_Norm_OF_Psi(psi, Norm0)
      CALL Calc_Norm_OF_Psi(psi_dt, Norm)
      WRITE(out_unit,*) 'Norm psi_n = ', Norm0, '  Norm psi_np1 = ', Norm
      WRITE(out_unit,*) '>>> END march_CN_With_GMRES'
      FLUSH(out_unit)
   END IF
 END SUBROUTINE march_CN_With_GMRES


!>-------------------------------------------------------------------------------------
! @brief 
! @details 
!   This module provides tools to propagate the Time-Dependent Schrödinger
!   Equation (TDSE) using the Crank–Nicolson scheme combined with the GMRES
!   (Generalized Minimal Residual) iterative solver.
!
!   - Crank–Nicolson (CN) discretizes the TDSE into a linear system:
!       (I + i*dt/(2*hbar)*H) ψ(t+dt) = (I - i*dt/(2*hbar)*H) ψ(t)
!
!   - GMRES is then used to solve the linear system iteratively without
!     explicitly forming large matrices. This is efficient for high-dimensional
!     quantum problems where H is large and sparse.
 
 !--------------------------------------------------------------------------------------
 ! SUBROUTINE: matvec_cn_op
 !
 ! PURPOSE:
 !   Apply the Crank–Nicolson operator to a wavefunction:
 !       out_psi = (I + fac*H) psi
 !
 !   This routine is the "matrix-vector product" needed by GMRES.
 !   It avoids forming the full matrix explicitly: instead it applies H to ψ
 !   using the external operator routine `calc_OpPsi`.
 !
 ! ARGUMENTS:
 !   H_op     (IN)  : Hamiltonian operator object
 !   psi      (IN)  : Input wavefunction
 !   out_psi (OUT)  : Result of applying CN operator
 !   fac      (IN)  : Complex scaling factor (± i*dt/(2*hbar))
 !--------------------------------------------------------------------------------------
 SUBROUTINE matvec_cn_op(H_op, psi, out_psi, fac)
   TYPE(op_t),    INTENT(IN)  :: H_op
   TYPE(psi_t),   INTENT(IN)  :: psi
   TYPE(psi_t),   INTENT(INOUT) :: out_psi
   COMPLEX(KIND=Rkind), INTENT(IN) :: fac
   TYPE(psi_t) :: h_psi_tmp

   
   IF (Debug_Enabled()) THEN
      WRITE(out_unit,*) 'BEGIN matvec_cn_op'
      IF(.NOT. ASSOCIATED(psi%Basis)) THEN
         WRITE(out_unit,*) 'psi%Basis is not associated in matvec_cn_op routine'
         STOP 'psi%Basis is not associated in matvec_cn_op routine'
      ELSE
         WRITE(out_unit,*) 'psi%Basis is associated in matvec_cn_op routine'
         FLUSH(out_unit)
      END IF
      FLUSH(out_unit)
   END IF

   IF (Debug_Enabled()) THEN
      WRITE(out_unit,*) 'BEGIN matvec_cn_op'
      IF(.NOT. ASSOCIATED(out_psi%Basis)) THEN
         WRITE(out_unit,*) 'out_psi%Basis is not associated in matvec_cn_op routine'
         STOP 'out_psi%Basis is not associated in matvec_cn_op routine'
      ELSE
         WRITE(out_unit,*) 'out_psi%Basis is associated in matvec_cn_op routine'
         FLUSH(out_unit)
      END IF
      FLUSH(out_unit)
   END IF

   
   CALL init_psi(h_psi_tmp, psi%Basis, cplx=.TRUE., grid=.FALSE.)
   CALL calc_OpPsi(H_op, psi, h_psi_tmp)
   out_psi%CVec = CZERO
   out_psi%CVec = psi%CVec + fac * h_psi_tmp%CVec

   ! Clean up
   CALL dealloc_psi(h_psi_tmp)

   IF (Debug_Enabled()) THEN
      WRITE(out_unit,*) 'END matvec_cn_op'
      FLUSH(out_unit)
   END IF
 END SUBROUTINE matvec_cn_op


 
 
 
 !------------------------------------------------------------------------------
 ! SUBROUTINE: gmres_solver
 !
 ! PURPOSE:
 !   Solve the linear system A x = b using the GMRES method, where the matrix A
 !   is defined implicitly by the CN operator:
 !       A = I + fac*H
 !
 !   - Builds a Krylov subspace using the Arnoldi proctise
 !   - Uses Modified Gram–Schmidt to orthogonalize basis vectors
 !   - Applies Givens rotations to maintain a small least-squares problem
 !   - Iteratively improves solution until residual < tol or max_iter reached
 !
 ! ARGUMENTS:
 !   H_op      (IN)    : Hamiltonian operator
 !   psi_b     (IN)    : Right-hand side vector b
 !   psi_x    (INOUT)  : On input initial gutise (often zero), on output solution
 !   max_iter  (IN)    : Maximum number of iterations
 !   tol       (IN)    : Residual tolerance
 !   cn_fac    (IN)    : Complex CN factor (i*dt/2)
 !------------------------------------------------------------------------------
 !===========================================================
! SUBROUTINE: gmres_solver_stable
! PURPOSE:
!   Solve (I + cn_fac*H) x = b using GMRES with Arnoldi + Givens rotations
!   This version avoids NaNs and instabilities.
!===========================================================


 !=============================================
! GMRES solver with Crank-Nicolson propagation
!=============================================
SUBROUTINE gmres_solver(H_op, psi_b, psi_x, max_iter, tol, cn_fac, debug_enabled)
   IMPLICIT NONE
   TYPE(op_t),   INTENT(IN)    :: H_op
   TYPE(psi_t),  INTENT(IN)    :: psi_b
   TYPE(psi_t),  INTENT(INOUT) :: psi_x
   INTEGER,      INTENT(IN)    :: max_iter
   REAL(KIND=Rkind),INTENT(IN) :: tol
   COMPLEX(KIND=Rkind), INTENT(IN):: cn_fac
 
   TYPE(psi_t), ALLOCATABLE :: Q(:)
   COMPLEX(KIND=Rkind), ALLOCATABLE :: H(:,:)
   TYPE(psi_t) :: r, v, tmp
   REAL(KIND=Rkind) :: beta, norm
   REAL(KIND=Rkind), ALLOCATABLE :: cs(:), sn(:)
   COMPLEX(KIND=Rkind), ALLOCATABLE :: g(:)
   INTEGER :: j, i
   LOGICAL, INTENT(IN) :: debug_enabled
 
   ! Check basis
   IF(.NOT. ASSOCIATED(psi_b%Basis)) STOP 'psi_b%Basis not associated'
 
   ! Allocate Krylov subspace
   ALLOCATE(Q(max_iter+1))
   DO i=1, max_iter+1
      CALL init_psi(Q(i), psi_b%Basis, cplx=.TRUE., grid=.FALSE.)
   END DO
 
   CALL init_psi(r, psi_b%Basis, cplx=.TRUE., grid=.FALSE.)
   CALL init_psi(v, psi_b%Basis, cplx=.TRUE., grid=.FALSE.)
   CALL init_psi(tmp, psi_b%Basis, cplx=.TRUE., grid=.FALSE.)
 
   ALLOCATE(H(max_iter+1,max_iter), cs(max_iter), sn(max_iter), g(max_iter+1))
   H = CZERO
   g = 0.0_Rkind
   cs = 0.0_Rkind
   sn = 0.0_Rkind
 
   ! Initial residual r = b - A*psi_x
   CALL matvec_cn_op(H_op, psi_x, tmp, cn_fac)
   r%CVec = psi_b%CVec - tmp%CVec
   CALL Calc_Norm_OF_Psi(r, beta)
   IF (beta < tol) THEN
      ! Already converged
      DO i=1,max_iter+1
         CALL dealloc_psi(Q(i))
      END DO
      CALL dealloc_psi(r); CALL dealloc_psi(v); CALL dealloc_psi(tmp)
      DEALLOCATE(H, cs, sn, g)
      RETURN
   END IF
 
   ! Normalize first Krylov vector
   Q(1)%CVec = r%CVec / beta
   g = 0.0_Rkind
   g(1) = beta
 
   ! Main Arnoldi + GMRES loop
   DO j = 1, max_iter
      ! v = A * q_j
      CALL matvec_cn_op(H_op, Q(j), v, cn_fac)
 
      ! Modified Gram-Schmidt
      DO i = 1, j
         H(i,j) = DOT_PRODUCT(CONJG(Q(i)%CVec), v%CVec)
         v%CVec = v%CVec - H(i,j) * Q(i)%CVec
      END DO
 
      ! Re-orthogonalize once
      DO i = 1, j
         H(i,j) = H(i,j) + DOT_PRODUCT(CONJG(Q(i)%CVec), v%CVec)
         v%CVec = v%CVec - DOT_PRODUCT(CONJG(Q(i)%CVec), v%CVec) * Q(i)%CVec
      END DO
 
      CALL Calc_Norm_OF_Psi(v, norm)
      H(j+1,j) = norm
      IF (ABS(H(j+1,j)) < 1.0E-14_Rkind) EXIT
      Q(j+1)%CVec = v%CVec / H(j+1,j)
 
      ! Apply Givens rotations to Htiseenberg
      DO i=1,j-1
         CALL apply_givens(cs(i), sn(i), H(i,j), H(i+1,j))
      END DO
 
      ! Compute new Givens rotation
      CALL givens_rotation_params(H(j,j), H(j+1,j), cs(j), sn(j))
      CALL apply_givens(cs(j), sn(j), H(j,j), H(j+1,j))
      CALL apply_givens(cs(j), sn(j), g(j), g(j+1))
 
      ! Update solution psi_x
      CALL solve_leastsq_givens(j, H, g, Q, psi_x)
 
      ! Check convergence
      IF (ABS(g(j+1)) < tol) EXIT
   END DO
 
   ! Cleanup
   DO i=1,max_iter+1
      CALL dealloc_psi(Q(i))
   END DO
   CALL dealloc_psi(r); CALL dealloc_psi(v); CALL dealloc_psi(tmp)
   DEALLOCATE(H, cs, sn, g)
 END SUBROUTINE gmres_solver
 
 !-------------------------------------------
 SUBROUTINE solve_leastsq_givens(j, H, g, Q, psi_x)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: j
   COMPLEX(KIND=Rkind), INTENT(IN) :: H(:,:)
   COMPLEX(KIND=Rkind), INTENT(INOUT) :: g(:)
   TYPE(psi_t), INTENT(IN) :: Q(:)
   TYPE(psi_t), INTENT(INOUT) :: psi_x
   INTEGER :: i
   COMPLEX(KIND=Rkind), ALLOCATABLE :: y(:)
   ! Solve upper-triangular system H*y = g(1:j)
   ALLOCATE(y(j))
   y = g(1:j)
   DO i=j,1,-1
      IF (ABS(H(i,i)) > 1.0E-16_Rkind) THEN
         y(i) = (y(i) - SUM(H(i,i+1:j) * y(i+1:j))) / H(i,i)
      ELSE
         y(i) = 0.0_Rkind
      END IF
   END DO
   ! Update psi_x
   psi_x%CVec = psi_x%CVec + MATMUL(SPREAD(y,2,SIZE(Q(1)%CVec)), Q(j)%CVec)
   DEALLOCATE(y)
 END SUBROUTINE solve_leastsq_givens
 
 !-------------------------------------------
 SUBROUTINE givens_rotation_params(a, b, c, s)
   IMPLICIT NONE
   COMPLEX(KIND=Rkind), INTENT(IN) :: a, b
   REAL(KIND=Rkind), INTENT(OUT) :: c, s
   REAL(KIND=Rkind) :: scale
   IF (ABS(b) < 1.0E-16_Rkind) THEN
      c = 1.0_Rkind; s = 0.0_Rkind
   ELSE IF (ABS(b) > ABS(a)) THEN
      scale = ABS(a/b)
      s = 1.0_Rkind / SQRT(1.0_Rkind + scale*scale)
      c = scale * s
   ELSE
      scale = ABS(b/a)
      c = 1.0_Rkind / SQRT(1.0_Rkind + scale*scale)
      s = scale * c
   END IF
 END SUBROUTINE givens_rotation_params
 
 !-------------------------------------------
 SUBROUTINE apply_givens(c, s, a, b)
   IMPLICIT NONE
   REAL(KIND=Rkind), INTENT(IN) :: c, s
   COMPLEX(KIND=Rkind), INTENT(INOUT) :: a, b
   COMPLEX(KIND=Rkind) :: temp
   temp = c*a + s*b
   b = -s*a + c*b
   a = temp
 END SUBROUTINE apply_givens
!-------------------------------------------------------
 ! solve_complex_linear_system : simple Gaussian elimination with partial pivoting
 ! Solves A x = b, where A is complex (n x n), b complex (n)
 ! Overwrites x with solution.
 ! -------------------------------------------------------
 SUBROUTINE solve_complex_linear_system(A, b, x, n)
   INTEGER, INTENT(IN) :: n
   COMPLEX(KIND=Rkind), INTENT(INOUT) :: A(:,:)
   COMPLEX(KIND=Rkind), INTENT(INOUT) :: b(:)
   COMPLEX(KIND=Rkind), INTENT(OUT) :: x(:)
 
   INTEGER :: i, j, k, piv
   REAL(KIND=Rkind) :: maxval, aval
   COMPLEX(KIND=Rkind) :: maxabs, tmp, factor
   REAL(KIND=Rkind) :: absval
 
   ! We'll perform in-place elimination on A and b
   DO k = 1, n-1
      ! partial pivoting: find row with max |A(i,k)|
      maxval = ZERO
      piv = k
      DO i = k, n
         absval = ABS(A(i,k))
         IF (absval > maxval) THEN
            maxval = absval
            piv = i
         END IF
      END DO
      IF (maxval == ZERO) THEN
         WRITE(*,*) 'solve_complex_linear_system: singular matrix (pivot zero)'
         STOP
      END IF
      ! swap rows k and piv if needed
      IF (piv /= k) THEN
         A([k,piv],:) = A([piv,k],:)   ! array-swap trick not standard; do manual
         tmp = b(k); b(k) = b(piv); b(piv) = tmp
         ! Manual swap to be safe:
         DO j = 1, n
            tmp = A(k,j); A(k,j) = A(piv,j); A(piv,j) = tmp
         END DO
      END IF
 
      ! elimination
      DO i = k+1, n
         factor = A(i,k) / A(k,k)
         A(i,k) = CZERO
         DO j = k+1, n
            A(i,j) = A(i,j) - factor * A(k,j)
         END DO
         b(i) = b(i) - factor * b(k)
      END DO
   END DO
 
   ! back substitution: A is upper triangular
   x = CZERO
   DO i = n, 1, -1
      tmp = b(i)
      IF (i < n) THEN
         DO j = i+1, n
            tmp = tmp - A(i,j) * x(j)
         END DO
      END IF
      x(i) = tmp / A(i,i)
   END DO
 END SUBROUTINE solve_complex_linear_system
 
 !==============================================================================
! SUBROUTINE: solve_tdse_cn_gmres
!
! PURPOSE:
!   Propagate the TDSE one time step using Crank–Nicolson with GMRES solver.
!   Equation:
!       (I + i*dt/(2*hbar)*H) ψ(t+dt) = (I - i*dt/(2*hbar)*H) ψ(t)
!
!   - Builds the right-hand side vector b = (I - fac*H) ψ_n
!   - Initializes gutise ψ_{n+1}=0
!   - Calls GMRES to solve for ψ_{n+1}
!
! ARGUMENTS:
!   psi_n    (IN)  : Wavefunction at time t
!   H_op     (IN)  : Hamiltonian operator
!   psi_np1 (OUT)  : Wavefunction at time t+dt
!   dt,hbar  (IN)  : Time step and Planck constant (ℏ)
!   max_gmres(IN)  : Maximum GMRES iterations
!   gmres_tol(IN)  : Convergence tolerance
!==============================================================================
!>-------------------------------------------------------------------------------------
!
!>@ compute_Cheby_coeffs
!
!  PURPOSE:
!    Calculates the complex-valued coefficients `a_k` of the Chebyshev
!    expansion for the time-evolution operator in the TDSE. The calculation
!    is based on the formula: a_k = (2 - delta_{k0})(-i)^k J_k(arg), where
!    J_k is the Btiseel function of the first kind and `arg` is scaled
!    time step `delta_t` and spectral radius `delta_E`. The routine
!    stops when the magnitude of the Btiseel function drops below a
!    pre-defined tolerance `limit`.
!
!  ARGUMENTS:
!
!  INPUTS:
!    delta_E    : The energy difference, related to the spectral
!                 radius of the Hamiltonian.
!    delta_t    : The time step for the propagation.
!    limit      : The tolerance threshold for stopping the
!                 calculation.
!
!  OUTPUTS:
!    num_coeffs : The number of coefficients calculated.
!    a          : An array containing the calculated
!                           Chebyshev coefficients.
!
!-----------------------------------------------------------------------

SUBROUTINE Compute_Cheby_Coeffs(delta_E, delta_t, limit, num_coeffs, a)
   USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
   IMPLICIT NONE
   REAL(KIND=Rkind), INTENT(IN) :: delta_E, delta_t, limit
   INTEGER, INTENT(OUT) :: num_coeffs
   COMPLEX(KIND=Rkind), ALLOCATABLE, INTENT(OUT) :: a(:)
   INTEGER, PARAMETER :: max_order = 100000
   INTEGER :: k
   REAL(KIND=Rkind) :: arg_bessel, abs_arg, prev_bessel, curr_bessel
   REAL(KIND=Rkind) :: ratio, tail_bound
   LOGICAL :: converged

   IF (.NOT. IEEE_IS_FINITE(delta_E) .OR. .NOT. IEEE_IS_FINITE(delta_t) .OR. &
       .NOT. IEEE_IS_FINITE(limit)) ERROR STOP 'Compute_Cheby_Coeffs: non-finite input.'
   IF (delta_E <= ZERO) ERROR STOP 'Compute_Cheby_Coeffs: E_max must be greater than E_min.'
   IF (limit <= ZERO .OR. limit >= ONE) ERROR STOP 'Compute_Cheby_Coeffs: eps must be in (0,1).'

   arg_bessel = delta_E * delta_t
   abs_arg = ABS(arg_bessel)
   IF (.NOT. IEEE_IS_FINITE(abs_arg) .OR. abs_arg > REAL(max_order-1000,Rkind)) &
      ERROR STOP 'Compute_Cheby_Coeffs: scaled time is too large for a safe expansion.'

   IF (arg_bessel == ZERO) THEN
      num_coeffs = 1
      ALLOCATE(a(1))
      a(1) = CONE
      RETURN
   END IF

   ! Omit the tail only after k > |z|, where the Btiseel coefficients are
   ! in their decreasing regime. The geometric estimate bounds 2*sum |J_k(z)|.
   num_coeffs = 1
   converged = .FALSE.
   prev_bessel = bessel_jn(0,arg_bessel)
   DO k = 1, max_order
      curr_bessel = bessel_jn(k,arg_bessel)
      IF (REAL(k,Rkind) > abs_arg .AND. REAL(k,Rkind) >= abs_arg + TWO) THEN
         IF (ABS(prev_bessel) > TINY(ONE)) THEN
            ratio = ABS(curr_bessel / prev_bessel)
         ELSE IF (ABS(curr_bessel) <= TINY(ONE)) THEN
            ratio = ZERO
         ELSE
            ratio = ONE
         END IF
         IF (ratio < ONE) THEN
            tail_bound = TWO * ABS(curr_bessel) / MAX(ONE-ratio, SQRT(EPSILON(ONE)))
            IF (tail_bound <= limit) THEN
               num_coeffs = k       ! retain orders 0,...,k-1
               converged = .TRUE.
               EXIT
            END IF
         END IF
      END IF
      prev_bessel = curr_bessel
   END DO
   IF (.NOT. converged) ERROR STOP 'Compute_Cheby_Coeffs: Btiseel tail did not converge.'

   ALLOCATE(a(num_coeffs))
   a(1) = CMPLX(bessel_jn(0,arg_bessel), ZERO, KIND=Rkind)
   DO k = 1, num_coeffs-1
      SELECT CASE (MOD(k,4))
      CASE (0)
         a(k+1) = CMPLX(TWO*bessel_jn(k,arg_bessel), ZERO, KIND=Rkind)
      CASE (1)
         a(k+1) = CMPLX(ZERO, -TWO*bessel_jn(k,arg_bessel), KIND=Rkind)
      CASE (2)
         a(k+1) = CMPLX(-TWO*bessel_jn(k,arg_bessel), ZERO, KIND=Rkind)
      CASE (3)
         a(k+1) = CMPLX(ZERO, TWO*bessel_jn(k,arg_bessel), KIND=Rkind)
      END SELECT
   END DO
END SUBROUTINE Compute_Cheby_Coeffs
 
!-------------------------------------------------------------------------
!
!  SUBROUTINE: Apply_Hprime_to_psi
!
!  PURPOSE:
!    Calculates the action of the scaled and shifted Hamiltonian H' on a
!    wavefunction psi, defined as H'|psi> = ((H - alpha)/delta_E)|psi>.
!    This routine relies on an external subroutine 'calc_Oppsi' to compute
!    the action of the un-modified Hamiltonian H on psi.
!
!  ARGUMENTS:
!
!  INPUTS:
!    H          : The Hamiltonian operator.
!    psi        : The input wavefunction.
!    alpha      : The energy shift value.
!    delta_E    : The energy scaling factor (spectral radius).
!
!  OUTPUTS:
!    Hprime_psi : The resulting wavefunction after the application of H'.
!
!-------------------------------------------------------------------------
SUBROUTINE Apply_Hprime_To_Psi(H, psi,Hprime_psi, alpha, delta_E)
   IMPLICIT NONE
   
   ! Inputs
   TYPE(Op_t), INTENT(IN) :: H
   TYPE(psi_t), INTENT(IN) :: psi
   REAL(KIND=Rkind), INTENT(IN) :: alpha, delta_E

   
   ! Outputs
   TYPE(psi_t), INTENT(INOUT) :: Hprime_psi
   
   ! Local variables
   TYPE(psi_t)                 :: Hpsi
 
   ! Initialize temporary vector Hpsi
   CALL init_psi(Hpsi, psi%Basis, cplx=.TRUE., grid=.FALSE.)
   Hpsi%CVec(:) = CZERO
   Hprime_psi%CVec(:) = CZERO
   
   IF (Debug_Enabled()) THEN
      WRITE(out_unit,*) 'BEGINNING Apply_Hprime_To_Psi alpha = ', alpha, ' ΔE = ', delta_E
      FLUSH(out_unit)
   END IF
   ! Étape 1 : Calculer H|psi> en utilisant la routine existante
   ! La routine calc_Oppsi doit être acctiseible
   call calc_Oppsi(H, psi, Hpsi)
   
   ! Étape 2 : Appliquer le décalage (H|psi> - alpha|psi>)
   ! et la mise à l'échelle (division par delta_E)
   
     Hprime_psi%CVec(:) = (Hpsi%CVec(:) - alpha * psi%CVec(:)) / delta_E
 
     ! Cleanup
     CALL dealloc_psi(Hpsi)
   
 END SUBROUTINE Apply_Hprime_To_Psi



 !>--------------------------------------------------------------------------
!>
!>  SUBROUTINE: cheby_step
!>
!>  PURPOSE:
!>    Performs one full step of the Chebyshev propagation recurrence.
!>    It calculates the action of the scaled and shifted Hamiltonian H' on
!>    the previous state phi_k, and then applies the three-term recurrence
!>    relation to compute the next state vector phi_kp1 based on phi_k and
!>    phi_km1.
!
!>  ARGUMENTS:
!
!>  INPUTS:
!>    H          : Representation of the Hamiltonian operator.
!>    phi_k      : The state vector at step k (corresponds to T_k(H')*chi).
!>    phi_km1    : The state vector at step k-1 (corresponds to T_{k-1}(H')*chi).
!>    alpha      : The energy shift value for H'.
!>    delta_E    : The energy scaling factor for H'.
!
!>  OUTPUTS:
!>    phi_kp1    : The resulting state vector at step k+1.
!
!>--------------------------------------------------------------------------
SUBROUTINE Cheby_Step(H, phi_k, phi_km1, phi_kp1, alpha, delta_E)
  
   IMPLICIT NONE
   
   ! Inputs
   TYPE(Op_t), INTENT(IN)  :: H
   TYPE(psi_t), INTENT(IN) :: phi_k
   TYPE(psi_t), INTENT(IN) :: phi_km1
   REAL(KIND=Rkind), INTENT(IN) :: alpha, delta_E
   
   ! Outputs
   TYPE(psi_t), INTENT(INOUT) :: phi_kp1
     ! Local variables
   TYPE(psi_t) :: Hprime_phi_k
   
 
   ! Initialize output wavefunction
   CALL init_psi(Hprime_phi_k, phi_k%Basis, cplx=.TRUE., grid=.FALSE.)
   Hprime_phi_k%CVec(:) = CZERO
   
   !------------------------------------------
   ! Step 1: Calculate the action of the scaled Hamiltonian H' on phi_k
   ! H' * phi_k = ((H - alpha)/delta_E) * phi_k
   !------------------------------------------
   CALL Apply_Hprime_to_psi(H, phi_k, Hprime_phi_k, alpha, delta_E)
   
   !------------------------------------------
   ! Step 2: Apply the Chebyshev recurrence relation
   ! phi_{k+1} = 2 * (H' * phi_k) - phi_{k-1}
   !------------------------------------------
     phi_kp1%CVec(:) = TWO * Hprime_phi_k%CVec(:) - phi_km1%CVec(:)
 
     ! Clean up
   CALL dealloc_psi(Hprime_phi_k)
   
 END SUBROUTINE Cheby_Step


!----------------------------------------------------------------------------------
!> SUBROUTINE: march_Chebychev
!> Description: Wrapper for time propagation using the Chebyshev method.
!
!>              This routine computes the Chebyshev expansion of the time-evolution
!>              operator and applies it to the wavefunction to propagate it over a
!>              time step delta_t.
!
! Arguments:
!   - psi       : (IN)     Wavefunction at time t (unchanged).
!   - psi_dt    : (INOUT)  Wavefunction at time t + delta_t (modified).
!   - t         : (IN)     Current time.
!   - propa     : (IN)     Propagation parameters (delta_t, alpha, delta_E).
!   - H         : (IN)     Hamiltonian operator.
!----------------------------------------------------------------------------------
SUBROUTINE march_Chebychev(psi, psi_dt, t, propa, H)
   IMPLICIT NONE

   !-------------- INPUTS ------------------------------------------------
   TYPE(psi_t),   INTENT(IN)     :: psi
   REAL(KIND=Rkind), INTENT(IN)  :: t
   TYPE(propa_t), INTENT(IN)     :: propa
   TYPE(Op_t),    INTENT(IN)     :: H

   !-------------- OUTPUTS ----------------------------------------------
   TYPE(psi_t),   INTENT(INOUT)  :: psi_dt

   !-------------- Local variables --------------------------------------
   TYPE(psi_t)                        :: phi_0, phi_1, phi_k_minus_1, phi_k, phi_k_plus_1
   REAL(KIND=Rkind)                   :: Norm, Norm0, alpha, delta_E
   COMPLEX(KIND=Rkind), ALLOCATABLE   :: cheby_coeffs(:)
   INTEGER                            :: num_coeffs, k
   COMPLEX(KIND=Rkind)                :: exp_alpha_term
   LOGICAL :: debug

   ! --- Begin propagation ---------------------------------------------
   debug = propa%debug
   IF (debug) THEN
      WRITE(out_unit,*) 'BEGINNING march_Chebychev ', t, ' Δt = ', propa%delta_t
      FLUSH(out_unit)
   END IF

   ! 1. Extract propagation parameters and calculate exponential term
   alpha = HALF*(propa%E_min+propa%E_max)
   delta_E = HALF*(propa%E_max-propa%E_min)
   exp_alpha_term = EXP(-EYE*alpha * propa%delta_t)

   ! 2. Calculate Btiseel coefficients a_k
   CALL compute_Cheby_coeffs(delta_E, propa%delta_t, propa%eps, &
                             num_coeffs, cheby_coeffs)

   ! 3. Initialize state vectors
   CALL init_psi(phi_0, psi%Basis, cplx=.TRUE., grid=.FALSE.)
   CALL init_psi(phi_1, psi%Basis, cplx=.TRUE., grid=.FALSE.)
   CALL init_psi(phi_k_minus_1, psi%Basis, cplx=.TRUE., grid=.FALSE.)
   CALL init_psi(phi_k, psi%Basis, cplx=.TRUE., grid=.FALSE.)
   CALL init_psi(phi_k_plus_1, psi%Basis, cplx=.TRUE., grid=.FALSE.)

   phi_0%CVec(:) = CZERO
   phi_1%CVec(:) = CZERO
   phi_k_minus_1%CVec(:) = CZERO
   phi_k%CVec(:) = CZERO
   phi_k_plus_1%CVec(:) = CZERO



   ! phi_0 = psi (initial state)
    phi_0%CVec(:) = psi%CVec(:)

   ! 4. Initialise l'accumulation de la solution propagée
   psi_dt%CVec(:) = cheby_coeffs(1) * phi_0%CVec(:)  ! Terme a_0 * phi_0

   ! The expansion may contain only its k=0 term for a very small time step.
   IF (num_coeffs > 1) THEN
      CALL Apply_Hprime_To_Psi(H, phi_0, phi_1, alpha, delta_E)
      psi_dt%CVec(:) = psi_dt%CVec(:) + cheby_coeffs(2) * phi_1%CVec(:)
      phi_k_minus_1%CVec = phi_0%CVec
      phi_k%CVec = phi_1%CVec

      IF (num_coeffs > 2) THEN
         DO k = 2, num_coeffs-1
            CALL Cheby_Step(H, phi_k, phi_k_minus_1, phi_k_plus_1, alpha, delta_E)
            psi_dt%CVec(:) = psi_dt%CVec(:) + cheby_coeffs(k+1) * phi_k_plus_1%CVec(:)
            phi_k_minus_1%CVec = phi_k%CVec
            phi_k%CVec = phi_k_plus_1%CVec
         END DO
      END IF
   END IF

   ! 7. Apply the global phase factor exp(-i * alpha * delta_t)
   psi_dt%CVec(:) = exp_alpha_term * psi_dt%CVec(:)

   ! 8. Clean up local state vectors
   CALL dealloc_psi(phi_0)
   CALL dealloc_psi(phi_1)
   CALL dealloc_psi(phi_k_minus_1)
   CALL dealloc_psi(phi_k)
   CALL dealloc_psi(phi_k_plus_1)
   DEALLOCATE(cheby_coeffs)

   ! --- Debug: check norm preservation --------------------------------
   IF (debug) THEN
      CALL Calc_Norm_OF_Psi(psi, Norm0)
      CALL Calc_Norm_OF_Psi(psi_dt, Norm)
      WRITE(out_unit,*) '<psi_dt|psi_dt> = ', Norm,'  |ΔNorm| = ', ABS(Norm0 - Norm)
      WRITE(out_unit,*) 'END march_Chebychev'
      FLUSH(out_unit)
   END IF

END SUBROUTINE march_Chebychev

 
 
 
 


!=============================================================
! Subroutine : relax_propagation
! Purpose    : Imaginary-time propagation with selectable propagator
!=============================================================

SUBROUTINE relax_propagation(psi, H, propa)
   IMPLICIT NONE

   TYPE(psi_t),   INTENT(INOUT) :: psi         ! Input/output: wavefunction
   TYPE(Op_t),    INTENT(IN)    :: H            ! Input: Hamiltonian operator
   TYPE(propa_t), INTENT(IN)    :: propa        ! Input: propagation parameters

   ! Local variables
   REAL(KIND=Rkind) :: E_old, E_new, diffE, Norm
   LOGICAL :: debug
   INTEGER :: iter, iprint
   LOGICAL :: converged
   REAL(KIND=Rkind) :: dtau

   ! Temporary vector for Hamiltonian application
   TYPE(psi_t) :: Hpsi

   ! Initialize temporary vector Hpsi
   CALL init_psi(Hpsi, psi%basis, cplx=.TRUE., grid=.FALSE.)

   ! Set propagation parameters
   dtau = propa%delta_t
   converged = .FALSE.
   iprint = 1

   ! Compute initial energy
   CALL Calc_Av_E(E_old, psi, H)
   debug = propa%debug
   IF (debug) WRITE(out_unit,*) '>>> BEGIN imaginary-time propagation'
   IF (debug) WRITE(out_unit,*) 'Initial energy = ', E_old

   ! Iterative propagation loop
   DO iter = 1, propa%max_iter

      ! Apply Hamiltonian to psi: Hpsi = -H|psi>
      CALL mEyeHPsi(psi, Hpsi, H, relax=.TRUE.)

      ! Imaginary-time propagation: psi <- psi + dtau * Hpsi
      psi%CVec(:) = psi%CVec(:) + dtau * Hpsi%CVec(:)

      ! Normalize wavefunction
      CALL Calc_Norm_OF_Psi(psi, Norm)
      psi%CVec(:) = psi%CVec(:) /Norm

      ! Compute new energy
      CALL Calc_Av_E(E_new, psi, H)
      diffE = ABS(E_new - E_old)

      ! Output progrtise
      IF (debug .AND. MOD(iter, iprint) == 0) THEN
         WRITE(out_unit,*) ' it', iter, '  E = ', E_new, '  ΔE = ', diffE
      END IF

      ! Check for convergence
      IF (diffE < propa%eps) THEN
         converged = .TRUE.
         EXIT
      END IF

      E_old = E_new
   END DO

   ! Deallocate temporary vector
   CALL dealloc_psi(Hpsi)

   ! Final result
   IF (converged) THEN
      WRITE(out_unit,*) '>>> Converged after ', iter, ' iterations.'
      WRITE(out_unit,*) '>>> Ground-state energy = ', E_new
   ELSE
      WRITE(out_unit,*) '>>> WARNING: max_iter reached. Final E=', E_new
   END IF

   IF (debug) THEN
      WRITE(out_unit,*) '>>> END imaginary-time propagation'
      FLUSH(out_unit)
   END IF

END SUBROUTINE relax_propagation




!----------------------------------------------------------------------
! SUBROUTINE: read_propa
!
! DESCRIPTION:
!   Reads propagation parameters (time interval, time step, tolerance,
!   method names, and various options) from the "prop" namelist and
!   stores them in the "propa" structure.
!   Method names are converted to lowercase.
!----------------------------------------------------------------------

SUBROUTINE read_propa(propa)
   IMPLICIT NONE

   TYPE(propa_t), INTENT(INOUT) :: propa
   REAL(KIND=Rkind)             :: t0, tf, delta_t, eps, E_min, E_max
   CHARACTER(LEN=80)            :: propa_name, propa_name2, selector
   CHARACTER(LEN=16)            :: dynamics, tise_method, energy_unit
   CHARACTER(LEN=512)           :: tise_energy_file, tise_transition_file
   INTEGER                      :: max_iter, separator, split_order
   INTEGER                      :: tise_nstates, tise_max_iter, tise_relax_max_iter
   REAL(KIND=Rkind)             :: tise_tol, tise_relax_delta_t, tise_relax_tol
   LOGICAL                      :: Beta, P, renorm, action_H, writeV, tise_relax
   LOGICAL                      :: debug, write_Qt, write_SQt, write_Pt, write_E
   LOGICAL                      :: write_Norm, write_pop, write_alpha, write_auto_cor
   LOGICAL                      :: write_validation, write_E_13, write_Norm_13
   LOGICAL                      :: restart, checkpoint
   INTEGER(KIND=Ikind)          :: checkpoint_every
   CHARACTER(LEN=512)           :: checkpoint_file

   ! propa_name accepts a single combined value, e.g. 'standard_taylor' or
   ! 'hagedorn_taylor'. propa_name2 remains accepted only as a legacy alias.
   NAMELIST /prop/ t0, tf, delta_t, max_iter, split_order, eps, E_min, E_max, propa_name, propa_name2, &
                   Beta, P, renorm, action_H, writeV, dynamics, tise_method, tise_nstates, &
                   tise_max_iter, tise_tol, energy_unit, tise_energy_file, tise_transition_file, &
                   tise_relax, tise_relax_delta_t, tise_relax_max_iter, tise_relax_tol, debug, &
                   write_Qt, write_SQt, write_Pt, write_E, write_Norm, write_pop, &
                   write_alpha, write_auto_cor, write_validation, write_E_13, write_Norm_13, &
                   restart, checkpoint, checkpoint_every, checkpoint_file

   t0 = ZERO
   tf = TEN
   delta_t = ONETENTH**3
   eps = ONETENTH**10
   E_min = ZERO
   E_max = TEN
   max_iter = 500
   split_order = 2
   propa_name = 'standard_taylor'
   propa_name2 = ''
   Beta = .TRUE.
   P = .TRUE.
   renorm = .TRUE.
   action_H = .FALSE.
   writeV = .FALSE.
   dynamics = 'TDSE'
   tise_method = 'dense'
   tise_nstates = 1
   tise_max_iter = 100
   tise_tol = ONETENTH**8
   energy_unit = 'au'
   tise_energy_file = 'TISE_energies.dat'
   tise_transition_file = 'TISE_transitions.txt'
   tise_relax = .FALSE.
   tise_relax_delta_t = ONETENTH**3
   tise_relax_max_iter = 500
   tise_relax_tol = ONETENTH**10
   debug = .FALSE.
   write_Qt = .FALSE.
   write_SQt = .FALSE.
   write_Pt = .FALSE.
   write_E = .FALSE.
   write_Norm = .FALSE.
   write_pop = .FALSE.
   write_alpha = .FALSE.
   write_auto_cor = .FALSE.
   write_validation = .FALSE.
   write_E_13 = .FALSE.
   write_Norm_13 = .FALSE.
   restart = .FALSE.
   checkpoint = .TRUE.
   checkpoint_every = 1_Ikind
   checkpoint_file = 'TD_SCHROD.restart'

   READ (*, NML=prop)

   IF (delta_t <= ZERO) THEN
      WRITE(out_unit,*) 'Invalid propagation delta_t: ', delta_t, ' (must be > 0).'
      ERROR STOP 'read_propa: invalid delta_t'
   END IF
   IF (tf < t0) THEN
      WRITE(out_unit,*) 'Invalid propagation interval: t0=', t0, ' tf=', tf
      ERROR STOP 'read_propa: tf must be >= t0'
   END IF
   IF (max_iter < 1) THEN
      WRITE(out_unit,*) 'Invalid propagation max_iter: ', max_iter
      ERROR STOP 'read_propa: max_iter must be >= 1'
   END IF
   IF (tise_nstates < 1 .OR. tise_max_iter < 1 .OR. tise_tol <= ZERO) &
      ERROR STOP 'read_propa: invalid TISE parameters.'
   IF (tise_relax_delta_t <= ZERO .OR. tise_relax_max_iter < 1 .OR. tise_relax_tol <= ZERO) &
      ERROR STOP 'read_propa: invalid TISE relaxation parameters.'
   IF (checkpoint_every < 1_Ikind) ERROR STOP 'read_propa: checkpoint_every must be >= 1.'
   IF ((checkpoint .OR. restart) .AND. LEN_TRIM(checkpoint_file) == 0) &
      ERROR STOP 'read_propa: checkpoint_file cannot be empty when checkpointing/restarting.'

   selector = TRIM(ADJUSTL(propa_name))
   CALL string_uppercase_TO_lowercase(selector)
   separator = INDEX(TRIM(selector), '_')

   IF (TRIM(selector) == 'non_hagedorn' .OR. TRIM(selector) == 'no_hagedorn' .OR. &
       TRIM(selector) == 'nonhagedorn') THEN
      propa%propa_name = 'standard'
      IF (LEN_TRIM(propa_name2) > 0) THEN
         propa%propa_name2 = TRIM(ADJUSTL(propa_name2))
         CALL string_uppercase_TO_lowercase(propa%propa_name2)
      ELSE
         propa%propa_name2 = 'taylor'
      END IF
   ELSE IF (separator > 0) THEN
      IF (selector(:separator-1) == 'non' .AND. &
          INDEX(selector(separator+1:LEN_TRIM(selector)), 'hagedorn_') == 1) THEN
         ! Combined legacy alias non_hagedorn_<integrator>.
         propa%propa_name = 'standard'
         propa%propa_name2 = selector(separator+1+LEN('hagedorn_'):LEN_TRIM(selector))
      ELSE IF (selector(:separator-1) == 'no' .AND. &
               INDEX(selector(separator+1:LEN_TRIM(selector)), 'hagedorn_') == 1) THEN
         propa%propa_name = 'standard'
         propa%propa_name2 = selector(separator+1+LEN('hagedorn_'):LEN_TRIM(selector))
      ELSE
         propa%propa_name = selector(:separator-1)
         propa%propa_name2 = selector(separator+1:LEN_TRIM(selector))
      END IF
   ELSE
      ! Legacy two-key format: propa_name='hagedorn' and propa_name2='taylor'.
      propa%propa_name = selector
      IF (LEN_TRIM(propa_name2) > 0) THEN
         propa%propa_name2 = TRIM(ADJUSTL(propa_name2))
         CALL string_uppercase_TO_lowercase(propa%propa_name2)
      ELSE
         propa%propa_name2 = 'taylor'
      END IF
   END IF

   SELECT CASE (TRIM(propa%propa_name))
   CASE ('standard', 'non_hagedorn', 'no_hagedorn', 'nonhagedorn')
      propa%propa_name = 'standard'
   CASE ('hagedorn')
      propa%propa_name = 'hagedorn'
   CASE DEFAULT
      ERROR STOP "read_propa: unsupported propagation method in propa_name"
   END SELECT
   CALL string_uppercase_TO_lowercase(propa%propa_name2)

   IF (TRIM(propa%propa_name) == 'hagedorn' .AND. &
       (TRIM(propa%propa_name2) == 'itp' .OR. TRIM(propa%propa_name2) == 'gssil')) &
      ERROR STOP 'read_propa: ITP/GSSIL are standalone modes and cannot be combined with Hagedorn.'
   IF (TRIM(propa%propa_name) == 'hagedorn' .AND. TRIM(propa%propa_name2) == 'split') &
      ERROR STOP 'read_propa: Split requires a fixed periodic Fourier grid and cannot use a moving Hagedorn basis.'
   IF (TRIM(propa%propa_name2) == 'split') THEN
      IF (split_order /= 1 .AND. split_order /= 2 .AND. split_order /= 4 .AND. split_order /= 6) &
         ERROR STOP 'read_propa: split_order must be 1, 2, 4 or 6.'
   ELSE IF (split_order /= 2) THEN
      ERROR STOP 'read_propa: split_order is only meaningful with the split propagator.'
   END IF

   propa%t0 = t0
   propa%tf = tf
   propa%delta_t = delta_t
   propa%eps = eps
   propa%E_min = E_min
   propa%E_max = E_max
   propa%max_iter = max_iter
   propa%split_order = split_order
   propa%Beta = Beta
   propa%P = P
   propa%renorm = renorm
   propa%action_H = action_H
   propa%writeV = writeV
   propa%dynamics = TRIM(ADJUSTL(dynamics))
   propa%tise_method = TRIM(ADJUSTL(tise_method))
   propa%tise_nstates = tise_nstates
   propa%tise_max_iter = tise_max_iter
   propa%tise_tol = tise_tol
   propa%tise_relax = tise_relax
   propa%tise_relax_delta_t = tise_relax_delta_t
   propa%tise_relax_max_iter = tise_relax_max_iter
   propa%tise_relax_tol = tise_relax_tol
   propa%energy_unit = TRIM(ADJUSTL(energy_unit))
   propa%tise_energy_file = TRIM(tise_energy_file)
   propa%tise_transition_file = TRIM(tise_transition_file)
   CALL string_uppercase_TO_lowercase(propa%dynamics)
   CALL string_uppercase_TO_lowercase(propa%tise_method)
   CALL string_uppercase_TO_lowercase(propa%energy_unit)
   IF (TRIM(propa%dynamics) /= 'tdse' .AND. TRIM(propa%dynamics) /= 'tise') &
      ERROR STOP 'read_propa: dynamics must be TDSE or TISE.'
   IF (TRIM(propa%energy_unit) == 'cm1') propa%energy_unit = 'cm-1'
   IF (TRIM(propa%energy_unit) /= 'au' .AND. TRIM(propa%energy_unit) /= 'cm-1') &
      ERROR STOP 'read_propa: energy_unit must be au or cm-1.'
   propa%debug = debug
   CALL Set_Propa_Debug(debug)
   propa%write_Qt = write_Qt
   propa%write_SQt = write_SQt
   propa%write_Pt = write_Pt
   propa%write_E = write_E
   propa%write_Norm = write_Norm
   propa%write_pop = write_pop
   propa%write_alpha = write_alpha
   propa%write_auto_cor = write_auto_cor
   propa%write_validation = write_validation
   propa%write_E_13 = write_E_13
   propa%write_Norm_13 = write_Norm_13
   propa%restart = restart
   propa%checkpoint = checkpoint
   propa%checkpoint_every = checkpoint_every
   propa%checkpoint_file = TRIM(checkpoint_file)

   SELECT CASE (TRIM(propa%propa_name2))
   CASE ('rk1', 'rk2', 'rk3', 'rk4', 'taylor', 'sil', 'split', 'cn', 'cheby', 'vp', 'itp', 'gssil')
      CONTINUE
   CASE DEFAULT
      ERROR STOP "read_propa: unsupported propagation method in propa_name"
   END SELECT
END SUBROUTINE read_propa

!----------------------------------------------------------------------
! SUBROUTINE: mEyeHPsi
!
! DESCRIPTION:
!   Computes the action of the Hamiltonian operator H on a wavefunction
!   psi, multiplies the result by -i if relax is false, and -1 if relax is true, and stores it in Hpsi.
!
!   Mathematically:
!      Hpsi = -i * ( H * psi ) if relax is false
!      Hpsi = -1 * ( H * psi ) if relax is true
!
! ARGUMENTS:
!   psi   : Input wavefunction.
!   Hpsi  : Output wavefunction after applying -i*H.
!   H     : Hamiltonian operator.
!----------------------------------------------------------------------

SUBROUTINE mEyeHPsi(psi, Hpsi, H, relax)
   TYPE(psi_t), INTENT(IN)    :: psi
   TYPE(psi_t), INTENT(INOUT) :: Hpsi
   TYPE(Op_t),  INTENT(IN)    :: H
   LOGICAL,     INTENT(IN), OPTIONAL :: relax

   ! Apply the operator H to psi
   CALL calc_OpPsi(H, psi, Hpsi)

   ! Check if imaginary-time propagation is requested
   IF (PRESENT(relax)) THEN
      IF (relax) THEN
         ! Relaxation (imaginary time): multiply by -1
         Hpsi%CVec(:) = -Hpsi%CVec(:)
      ELSE
         ! Real time: multiply by -i
         Hpsi%CVec(:) = -EYE * Hpsi%CVec(:)
      END IF
   ELSE
      ! Default: real time propagation
      Hpsi%CVec(:) = -EYE * Hpsi%CVec(:)
   END IF

END SUBROUTINE mEyeHPsi

!----------------------------------------------------------------------
! SUBROUTINE: write_propa
!
! DESCRIPTION:
!   Writes the current propagation parameters stored in the "propa"
!   structure to the output unit (out_unit) in a clean, aligned format.
!   This is useful for logging or verifying the propagation settings
!   before or during the simulation.
!
! ARGUMENTS:
!   propa : Propagation parameters structure (read-only here).
!----------------------------------------------------------------------

SUBROUTINE write_propa(propa)
   USE psi_m
   IMPLICIT NONE

   TYPE(propa_t), INTENT(INOUT) :: propa

   WRITE (out_unit, '(A15, " : ", ES15.8)') 't0',         propa%t0
   WRITE (out_unit, '(A15, " : ", ES15.8)') 'tf',         propa%tf
   WRITE (out_unit, '(A15, " : ", ES15.8)') 'Δt',         propa%delta_t
   WRITE (out_unit, '(A15, " : ", ES15.8)') 'eps',        propa%eps
   WRITE (out_unit, '(A15, " : ", ES15.8)') 'E_min',      propa%E_min
   WRITE (out_unit, '(A15, " : ", ES15.8)') 'E_max',      propa%E_max
   WRITE (out_unit, '(A15, " : ", I0)')     'max_iter',   propa%max_iter
   WRITE (out_unit, '(A15, " : ", I0)')     'split_order', propa%split_order
   WRITE (out_unit, '(A15, " : ", A, "_", A)') 'propagator', TRIM(propa%propa_name), TRIM(propa%propa_name2)
   WRITE (out_unit, '(A15, " : ", L1)')     'Beta',       propa%Beta
   WRITE (out_unit, '(A15, " : ", L1)')     'P',          propa%P
   WRITE (out_unit, '(A15, " : ", L1)')     'renorm',     propa%renorm
   WRITE (out_unit, '(A15, " : ", L1)')     'debug',      propa%debug
   WRITE (out_unit, '(A15, " : ", L1)')     'write_Qt', propa%write_Qt
   WRITE (out_unit, '(A15, " : ", L1)')     'write_SQt', propa%write_SQt
   WRITE (out_unit, '(A15, " : ", L1)')     'write_Pt', propa%write_Pt
   WRITE (out_unit, '(A15, " : ", L1)')     'write_E', propa%write_E
   WRITE (out_unit, '(A15, " : ", L1)')     'write_Norm', propa%write_Norm
   WRITE (out_unit, '(A15, " : ", L1)')     'write_pop', propa%write_pop
   WRITE (out_unit, '(A15, " : ", L1)')     'write_alpha', propa%write_alpha
   WRITE (out_unit, '(A15, " : ", L1)')     'write_auto_cor', propa%write_auto_cor
   WRITE (out_unit, '(A15, " : ", L1)')     'write_validation', propa%write_validation
   WRITE (out_unit, '(A15, " : ", L1)')     'write_E_13', propa%write_E_13
   WRITE (out_unit, '(A15, " : ", L1)')     'write_Norm_13', propa%write_Norm_13
END SUBROUTINE write_propa

!----------------------------------------------------------------------
! SUBROUTINE: Calc_average_energy
!
! DESCRIPTION:
!   Calculates the average energy of a wavefunction Psi using the
!   Hamiltonian operator H. The energy is computed as <Psi | H | Psi>
!   and normalized by <Psi | Psi> to ensure proper scaling.
!
! ARGUMENTS:
!   Psi   : Input wavefunction.
!   E     : Output average energy.
!----------------------------------------------------------------------

   SUBROUTINE Calc_average_energy(Psi, E)
      !>-------------------------------------------------------
      !>     E = <Psi | H | Psi>
      !>--------------------------------------------------------
      TYPE(psi_t), intent(in)                        :: psi
      REAL(kind=Rkind), intent(inout)                :: E
      TYPE(psi_t)                                    :: Hpsi, psi_b
      TYPE(Op_t)                                     :: H
      REAL(KIND=Rkind)                               :: Norm
      if (Psi%Grid) then

         IF (Debug_Enabled()) PRINT*,"psi  is on Grid"
         CALL init_psi(psi_b, psi%Basis, cplx=.TRUE., grid=.false.)
         call GridTOBasis_nD_cplx(psi_b%CVec, psi%CVec, psi%Basis)
         CALL init_psi(Hpsi, psi%Basis, cplx=.TRUE., grid=.false.)
         call calc_OpPsi(H, psi_b, Hpsi)
         !stop 'verifie le calcul d E'
         IF (Debug_Enabled()) write(out_unit,*) 'rappele toi , sur la calcule de l energie'
         E = real(dot_product(Hpsi%CVec, psi_b%CVec), kind=Rkind)

      else
         IF (Debug_Enabled()) PRINT*,"psi is on basis"
         CALL init_psi(Hpsi, psi%Basis, cplx=.TRUE., grid=.false.)
         call calc_OpPsi(H, psi, Hpsi)
         !stop 'verifie le calcul d E'
         IF (Debug_Enabled()) write(out_unit,*) 'rappele toi , sur la calcule de l energie'
         E = real(dot_product(Hpsi%CVec, psi%CVec), kind=Rkind)

      end if
      call Calc_Norm_OF_Psi(psi, Norm)
      E = E/Norm**2
      !print *, "<psi|H|psi> = ", E, "<psi|psi> =", Norm

      CALL dealloc_psi(Hpsi)
      CALL dealloc_psi(psi_b)
   End SUBROUTINE Calc_average_energy

   !----------------------------------------------------------------------
   ! SUBROUTINE: Calc_Av_E
   !
   ! DESCRIPTION:
   !   Computes the expectation value E = <Psi | H | Psi>
   !   Normalizes by <Psi|Psi> to ensure proper scaling.
   !
   ! ARGUMENTS:
   !   E     : Output average energy.
   !   psi   : Input wavefunction.
   !   H     : Hamiltonian operator.
   !----------------------------------------------------------------------

   SUBROUTINE Calc_Av_E(E, psi, H)
      !> -----------------------------------------------------------
      !> Computes the expectation value E = <Psi | H | Psi>
      !> Normalizes by <Psi|Psi> to ensure proper scaling.
      !> Includes a local debug flag to optionally print E.
      !> -----------------------------------------------------------
      IMPLICIT NONE
   
      ! ===== Arguments =====
      REAL(KIND=Rkind), INTENT(INOUT) :: E
      TYPE(psi_t), INTENT(IN)         :: psi
      TYPE(Op_t),  INTENT(IN)         :: H
   
      ! ===== Local variables =====
      TYPE(psi_t)      :: Hpsi, psi_b
      REAL(KIND=Rkind) :: Norm
      LOGICAL          :: debug  ! Debug flag
   
      ! ===== Debug control =====
      debug = Debug_Enabled()
   
      ! ===== Expectation value calculation =====
      IF (psi%Grid) THEN
         ! Convert wavefunction from grid to basis representation
         CALL init_psi(psi_b, psi%Basis, cplx=.TRUE., grid=.FALSE.)
         CALL GridTOBasis_nD_cplx(psi_b%CVec, psi%CVec, psi%Basis)
   
         ! Apply operator H to wavefunction in basis
         CALL init_psi(Hpsi, psi%Basis, cplx=.TRUE., grid=.FALSE.)
         CALL calc_OpPsi(H, psi_b, Hpsi)
   
         ! Compute <psi_b | H | psi_b> with the physical metric.
         ! A packed Smolyak vector contains overlapping term blocks, so a
         ! plain DOT_PRODUCT would count coefficients with the wrong metric.
         IF (Basis_UsesSmolyakSparseGrid(psi_b%Basis)) THEN
            E = REAL(DotSmolyakBasis(psi_b%CVec, Hpsi%CVec, psi_b%Basis%SmolyakSparseGrid), KIND=Rkind)
         ELSE
            E = REAL(DOT_PRODUCT(Hpsi%CVec, psi_b%CVec), KIND=Rkind)
         END IF
      ELSE
         ! Directly apply H to wavefunction
         CALL init_psi(Hpsi, psi%Basis, cplx=.TRUE., grid=.FALSE.)
         CALL calc_OpPsi(H, psi, Hpsi)
   
         ! Compute <psi | H | psi> with the physical metric.
         IF (Basis_UsesSmolyakSparseGrid(psi%Basis)) THEN
            E = REAL(DotSmolyakBasis(psi%CVec, Hpsi%CVec, psi%Basis%SmolyakSparseGrid), KIND=Rkind)
         ELSE
            E = REAL(DOT_PRODUCT(Hpsi%CVec, psi%CVec), KIND=Rkind)
         END IF
      END IF
   
      ! ===== Normalize by the norm of psi =====
      CALL Calc_Norm_OF_Psi(psi, Norm)
      E = E / Norm**2
   
      ! ===== Debug output =====
      IF (debug) THEN
         WRITE(*,'(A,ES16.8)') 'DEBUG: <psi|H|psi> normalized = ', E
      END IF
   
      ! ===== Cleanup =====
      CALL dealloc_psi(Hpsi)
      CALL dealloc_psi(psi_b)
   END SUBROUTINE Calc_Av_E

!===============================================================
!> @brief   Crée et ouvre un fichier texte pour l'enregistrement
!>          des résultats d'une propagation.
!>
!> @param[out] nio    Unité logique attribuée automatiquement par NEWUNIT.
!> @param[in] name    Nom de base à utiliser pour créer le fichier.
!>
!> @details
!> Le nom final est simplement TRIM(name)//'.txt'. Les noms des
!> propagateurs ne sont volontairement pas ajoutés au nom du fichier.
!===============================================================
SUBROUTINE open_propa_file(nio, name, append)
   !==== Arguments ====
   INTEGER,          INTENT(OUT) :: nio      ! Unité attribuée par NEWUNIT
   CHARACTER(*),     INTENT(IN)  :: name     ! Nom de base du fichier
   LOGICAL, OPTIONAL, INTENT(IN) :: append   ! Ajoute à l'historique si vrai
   !==== Variables locales ====
   CHARACTER(512) :: name_tot
   INTEGER :: ios
   LOGICAL :: do_append
   name_tot = TRIM(name) // '.txt'
   do_append = .FALSE.
   IF (PRESENT(append)) do_append = append
   !==== Ouverture du fichier ====
   IF (do_append) THEN
      OPEN(NEWUNIT=nio, FILE=TRIM(name_tot), STATUS='UNKNOWN', POSITION='APPEND', ACTION='WRITE', &
           FORM='FORMATTED', IOSTAT=ios)
   ELSE
      OPEN(NEWUNIT=nio, FILE=TRIM(name_tot), STATUS='REPLACE', ACTION='WRITE', &
           FORM='FORMATTED', IOSTAT=ios)
   END IF
   IF (ios /= 0) ERROR STOP 'open_propa_file: cannot open output file'

END SUBROUTINE open_propa_file

     !==============================================================
     !> @brief
     !>   Write a complex autocorrelation value at a given time 
     !>   to a file unit with different output options.
     !>
     !> @details
     !>   Controlled by the character string 'op':
     !>     - "t_z"     : write (t, z)   where z is complex
     !>     - "t_norm"  : write (t, |z|)
     !>     - "t_re_im" : write (t, Re(z), Im(z))
     !>     - "t_re"    : write (t, Re(z))
     !>     - "t_im"    : write (t, Im(z))
     !>     - "z"       : write z only
     !>     - "norm"    : write |z| only
     !>     - "re_im"   : write Re(z), Im(z)
     !>
     !> @param[in] z    Complex autocorrelation value
     !> @param[in] nio  File unit number (must be open before calling)
     !> @param[in] op   Output format option (see list above)
     !> @param[in] t    Time value (required if op begins with "t_")
     !==============================================================
   SUBROUTINE Write_Complex(z, nio, op, t)
      IMPLICIT NONE
      COMPLEX(KIND=Rkind), INTENT(in) :: z(:)
      INTEGER, INTENT(in)              :: nio
      CHARACTER(*), INTENT(in)         :: op
      REAL(KIND=Rkind), INTENT(in)     :: t

      ! Test if z is a scalar (size=1) or vector (size>1)
      IF (SIZE(z) == 1) THEN
         CALL Write_Complex_scalar(z(1), nio, op, t)
      ELSE
         CALL Write_Complex_vector(z, nio, op, t)
      END IF

   END SUBROUTINE Write_Complex

   !=============================
   ! Cas scalaire
   !=============================
   SUBROUTINE Write_Complex_scalar(z, nio, op, t)
      IMPLICIT NONE
      COMPLEX(KIND=Rkind), INTENT(in) :: z
      INTEGER,             INTENT(in) :: nio
      CHARACTER(*),        INTENT(in) :: op
      REAL(KIND=Rkind),    INTENT(in) :: t

      REAL(KIND=Rkind) :: zr, zi, znorm
      zr    = REAL(z, KIND=Rkind)
      zi    = AIMAG(z)
      znorm = ABS(z)

      SELECT CASE (TRIM(ADJUSTL(op)))
      CASE ("t_z","t_re_im")
         WRITE(nio,*) t, zr, zi
      CASE ("t_norm")
         WRITE(nio,*) t, znorm
      CASE ("t_re")
         WRITE(nio,*) t, zr
      CASE ("t_im")
         WRITE(nio,*) t, zi
      CASE ("z","re_im")
         WRITE(nio,*) zr, zi
      CASE ("norm")
         WRITE(nio,*) znorm
      CASE DEFAULT
         WRITE(*,*) "Error: Unknown option in Write_Complex_scalar -> ", op
      END SELECT
   END SUBROUTINE Write_Complex_scalar

   !=============================
   ! Cas vecteur
   !=============================
   SUBROUTINE Write_Complex_vector(z, nio, op, t)
      IMPLICIT NONE
      COMPLEX(KIND=Rkind), DIMENSION(:), INTENT(in) :: z
      INTEGER,             INTENT(in) :: nio
      CHARACTER(*),        INTENT(in) :: op
      REAL(KIND=Rkind),    INTENT(in) :: t

      INTEGER ::n
      n = SIZE(z)

      SELECT CASE (TRIM(ADJUSTL(op)))
      CASE ("t_z","t_re_im")
         WRITE(nio,*) t, REAL(z(:),KIND=Rkind), AIMAG(z(:))
      CASE ("t_norm")
         WRITE(nio,*) t, ABS(z(:))
      CASE ("t_re")
         WRITE(nio,*) t, REAL(z(:),KIND=Rkind)
      CASE ("t_im")
         WRITE(nio,*) t, AIMAG(z(:))
      CASE ("z","re_im")
         WRITE(nio,*) REAL(z(:),KIND=Rkind), AIMAG(z(:))
      CASE ("norm")
         WRITE(nio,*) ABS(z(:))
      CASE DEFAULT
         WRITE(*,*) "Error: Unknown option in Write_Complex_vector -> ", op
      END SELECT
   END SUBROUTINE Write_Complex_vector


   !----------------------------------------------------------------------------------
!> SUBROUTINE: march_SIL
!> Description: This routine acts as a wrapper for the time propagation
!>              using the Short Iterative Lanczos (SIL) method. It handles the
!>              initial copying of the wavefunction and debug diagnostics (if enabled).
!>
!>              Pedagogically, this routine illustrates the step-by-step approach
!>              for simulation. It isolates the complexity of the propagation
!>              algorithm (`march_SIL` from the `march_SIL_m` module) from the
!>              main control flow. It also ensures that the original wavefunction
!>              (`psi`) remains unchanged, as the propagation is performed on a
!>              working copy (`psi_dt`).
!>
!> Arguments:
!>   - psi       : (IN)     The wavefunction at time t. It is not modified.
!>   - psi_dt    : (INOUT)  The wavefunction at time t + delta_t. Modified by the routine.
!>   - t         : (IN)     The current time.
!>   - propa     : (IN)     The propagation parameters (delta_t, eps, max_iter).
!>   - H         : (IN)     The Hamiltonian operator.
!----------------------------------------------------------------------------------
SUBROUTINE march_SIL(psi, psi_dt, t, propa, H)

   TYPE(psi_t), INTENT(INOUT)      :: psi_dt
   TYPE(psi_t), INTENT(IN)         :: psi
   TYPE(propa_t), INTENT(IN)       :: propa
   TYPE(Op_t), INTENT(IN)          :: H
   REAL(kind=Rkind), INTENT(IN)    :: t

   ! Local variables -------------------------------------------------------------
   REAL(kind=Rkind)                   :: Norm, Norm0
   LOGICAL :: debug

   ! --- Begin propagation ---
   debug = propa%debug
   IF (debug) THEN
       WRITE (out_unit, *) 'BEGINNING march_SIL ', t,' Δt = ', propa%delta_t
       FLUSH(out_unit)
   END IF

   ! Initialize psi_dt by copying the starting wavefunction
   psi_dt%CVec(:) = psi%CVec(:)
   
   ! Call the Short Iterative Lanczos (SIL) propagation routine
   CALL march_mSIL(H, psi_dt, psi%Basis, propa%delta_t, propa%max_iter, propa%eps)
   
   ! --- Calculate norms for debug diagnostics ---
   IF (debug) THEN
       ! Calculate the norm of the initial wavefunction
       call Calc_Norm_OF_Psi(psi, Norm0)
       ! Calculate the norm of the propagated wavefunction
       call Calc_Norm_OF_Psi(psi_dt, Norm)
       
       write (out_unit, *) '<psi_dt|psi_dt> = ', Norm, ' |ΔNorm| = ', abs(Norm0 - Norm)
       write (out_unit, *) 'END march_SIL'
       flush(out_unit)
   END IF

END SUBROUTINE


   
END MODULE sub_propa_m
