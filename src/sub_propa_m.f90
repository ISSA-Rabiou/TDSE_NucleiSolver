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
      CHARACTER(len=:), ALLOCATABLE  :: propa_name
      CHARACTER(len=:), ALLOCATABLE  :: propa_name2
      LOGICAL                        :: Beta 
      LOGICAL                        :: P
      LOGICAL                        :: renorm
   END TYPE propa_t

   ! Define public interfaces
   PUBLIC :: march_taylor, march
   PUBLIC :: write_propa, Analyse, creat_file_unit, read_propa
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
      call string_uppercase_TO_lowercase(name)
      select case (name)
      case ('rk1') 
        CALL march_RK1(psi, psi_dt, t, propa, H)
      case ('rk2') 
        CALL march_RK2(psi, psi_dt, t, propa, H)
      case ('rk3') 
        CALL march_RK3(psi, psi_dt, t, propa, H)
      case ('rk4') 
        CALL march_RK4(psi, psi_dt, t, propa, H)
      case ('taylor')
         CALL march_taylor(psi, psi_dt, t, propa, H)
      case ('sil') 
        CALL march_SIL(psi, psi_dt, t, propa,H)
      case ('split')
        CALL march_SPLIT(psi, psi_dt, t, propa, H)
      case ('cn')
        CALL march_CN(psi, psi_dt, t, propa, H)
      case ('cn_gmres')
         write (out_unit, *) 'cn_gmres is not implemented'
         flush(out_unit)
        !CALL march_CN_With_GMRES(psi, psi_dt, t, propa, H)
      case ('cheby')
        CALL march_Chebychev(psi, psi_dt, t, propa, H)
      case default
         write (out_unit, *) 'March name is not in the list'
         flush(out_unit)
      end select      
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

      !allocate variables ----------------------------------------------------
      Ndim = size(psi%Basis%tab_basis)
      ALLOCATE (Pop(psi%Basis%tab_basis(Ndim)%nb))
      ALLOCATE (Qp(psi%Basis%tab_basis(Ndim)%nb))
      ALLOCATE (Qm(Ndim - 1))

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
      
      write (3, *) t, E, Norm, pop
      write (4, *) t, Qm
      write (5, *) t, Qp
      
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
   
      ! Local variables-------------------------------------------------------------------------------
      REAL(KIND=Rkind)                 :: alpha
      REAL(KIND=Rkind)                 :: Norm, Norm0
      INTEGER                          :: kk
      LOGICAL                          :: debug
   
      ! initialize variables --------------------------------------------------
      CALL init_psi(Hpsi, psi%basis, cplx=.true., grid=.FALSE.)
      CALL init_psi(psi0, psi%basis, cplx=.true., grid=.FALSE.)
   
      debug = .true.
      alpha = TEN**10
      psi_dt%CVec = psi%CVec
      psi0%CVec = psi%CVec

      IF (debug) THEN
         WRITE (out_unit, *) 'BEGINNING march_taylor  ', t,'Δt=', propa%delta_t
         FLUSH(out_unit)
      END IF
   
      ! main loop -------------------------------------------------------------
      Do kk = 1, propa%max_iter, 1
         CALL mEyeHPsi(psi0, Hpsi, H)
         psi0%CVec(:) = Hpsi%CVec(:)*(propa%delta_t/kk)
         psi_dt%CVec(:) = psi_dt%CVec(:) + psi0%CVec(:)
         CALL Calc_Norm_OF_Psi(psi0, Norm)
         WRITE (out_unit, *) 'sqrt(<Hpsi|Hpsi>) = ', kk, Norm
    
         ! check conditions ------------------------------------------------------
         IF (Norm >= alpha) THEN
            STOP "wrong choice of delta_t"
         ELSEIF (Norm <= propa%eps) THEN
            WRITE (out_unit, *) 'Taylor condition is fulfilled after', kk, 'iteration'
            EXIT
         END IF
      END DO

      IF (debug) THEN
         CALL Calc_Norm_OF_Psi(psi, Norm0)
         CALL Calc_Norm_OF_Psi(psi_dt, Norm)
         WRITE (out_unit, *) '<psi_dt|psi_dt> = ', Norm, ' |ΔNorm| =', abs(Norm0 - Norm)
         WRITE (out_unit, *) 'END march_taylor'
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
      LOGICAL, PARAMETER :: debug=.TRUE.
      
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
   LOGICAL, PARAMETER :: debug=.TRUE.
   
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
   LOGICAL, PARAMETER :: debug=.TRUE.

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
   LOGICAL, PARAMETER :: debug=.TRUE.
   TYPE(psi_t) :: K1, K2, K3, K4, psi_tmp
   REAL(KIND=Rkind) :: Norm0, Norm

   CALL init_psi(K1, psi%basis, cplx=.TRUE., grid=.FALSE.)
   CALL init_psi(K2, psi%basis, cplx=.TRUE., grid=.FALSE.)
   CALL init_psi(K3, psi%basis, cplx=.true., grid=.FALSE.)
   CALL init_psi(K4, psi%basis, cplx=.true., grid=.FALSE.)
   CALL init_psi(psi_tmp, psi%basis, cplx=.true., grid=.FALSE.)

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
   ! ------------------------------------------------------------
   ! Propagates the wavefunction using a symmetric split-operator
   ! Scheme: psi(t+dt) ≈ e^{-i V dt/2} e^{-i T dt} e^{-i V dt/2} psi(t)
   !
   ! Steps implemented:
   !   1) Apply half-potential + transform -> momentum  (via ApplyExpVdtAndTransformToMomentum)
   !   2) Apply kinetic in p and transform -> position   (via ApplyExpTdtAndTransformToPosition)
   !   3) Apply final half-potential in position (in-place, no transform)
   !
   ! Notes:
   !   - psi%CVec and psi_dt%CVec are 1-D vectors that represent reshaped matrices
   !     (nq, nsurf) or (np, nsurf) depending on routine. We use pointer reshapes
   !     where needed in called routines; here we apply the final V-phase in x
   !     directly using the appropriate buffers.
   ! ------------------------------------------------------------
   IMPLICIT NONE
 
   TYPE(psi_t), INTENT(INOUT)       :: psi_dt
   TYPE(psi_t), INTENT(IN)          :: psi
   TYPE(propa_t), INTENT(IN)        :: propa
   TYPE(Op_t), INTENT(IN)           :: H
   REAL(KIND=Rkind), INTENT(IN)     :: t
 
   ! Local variables
   REAL(KIND=Rkind), ALLOCATABLE :: Q(:,:), WQ(:), P(:,:), WP(:),Q0d(:)
   TYPE(psi_t)                   :: psi_temp1, psi_temp2,psi0
   REAL(KIND=Rkind)              :: Norm0, Norm
   LOGICAL, PARAMETER            :: debug = .TRUE.
   LOGICAL, PARAMETER            :: transform = .FALSE.   ! utiliser .FALSE. => no transform in ApplyExpVdtAndTransformToMomentum
 
 
   IF (debug) THEN
     WRITE(out_unit, *) 'BEGINNING march_SPLIT', ' t=', t, ' Δt=', propa%delta_t
     WRITE(out_unit, *) ''
   END IF
 
   ! Initialize output wavefunction and temps (allocations internes à init_psi)
   CALL init_psi(psi0, psi%basis, cplx=.TRUE., grid=.TRUE.)
   CALL init_psi(psi_temp1, psi%basis, cplx=.TRUE., grid=.TRUE.)
   CALL init_psi(psi_temp2, psi%basis, cplx=.TRUE., grid=.TRUE.)

   IF(psi%Grid) THEN
      psi0%CVec = psi%CVec
   ELSE
      CALL BasisTOGrid_nD_cplx(psi0%CVec, psi%CVec, psi%Basis) 
   END IF     
 
   ! Build NUFFT grids (position Q and momentum P with their weights WQ, WP)
   CALL Build_NUFT_Grids(psi%basis, Q, WQ, P, WP, Q0d)
   ! ------------------------------
   ! Step 1: First half-potential and transform to momentum
   !   - input:  psi (position)
   !   - output: psi_temp (momentum)
   ! ------------------------------
   CALL ApplyExpVdtAndTransformToMomentum(psi0%CVec, psi_temp1%CVec,H%Scal_pot,propa%delta_t,Q, WQ, P, WP)
   ! ------------------------------
   ! Step 2: Full kinetic step in momentum then transform back to position
   !   - input:  psi_temp (momentum)
   !   - output: psi_temp2 (position)
   ! ------------------------------
   CALL ApplyExpTdtAndTransformToPosition(psi_temp1%CVec, psi_temp2%CVec, Q, WQ, P, WP,propa%delta_t)
   ! ------------------------------
   ! Step 3: Second half-potential in position (no transform)
   !   - apply exp(-i V dt/2) to psi_temp2 (position) and store in psi_dt (position)
   !
   ! Note: the routine ApplyExpVdtAndTransformToMomentum supports an optional
   !       'transform' flag; passing transform = .FALSE. means "apply V-phase in x
   !       only, do not transform to momentum". We therefore use:
   !          input  = psi_temp2%CVec (position)
   !          output = psi_dt%CVec   (position after V/2)
   ! ------------------------------
   CALL ApplyExpVdtAndTransformToMomentum(psi_temp2%CVec, psi_temp1%CVec,H%Scal_pot,propa%delta_t,Q, WQ, P, WP,transform)

   IF(psi_dt%Grid) THEN
      psi_dt%CVec = psi_temp1%CVec
   ELSE
      CALL GridTOBasis_nD_cplx(psi_dt%CVec, psi_temp1%CVec, psi%Basis) 
   END IF     
   ! ------------------------------
   ! Diagnostics
   ! ------------------------------
   IF (debug) THEN
     CALL Calc_Norm_OF_Psi(psi, Norm0)
     CALL Calc_Norm_OF_Psi(psi_dt, Norm)
     WRITE(out_unit, *) '<psi_dt|psi_dt> = ', Norm, ' |ΔNorm| = ', abs(Norm0 - Norm)
     WRITE(out_unit, *) 'END march_SPLIT'
     WRITE(out_unit, *) ''
     FLUSH(out_unit)
   END IF
 
   ! Clean up: free temporary psi structures and NUFFT grids
   CALL dealloc_psi(psi_temp1)
   CALL dealloc_psi(psi_temp2)
   CALL dealloc_psi(psi0)
   DEALLOCATE(Q, WQ, P, WP)
 
 END SUBROUTINE march_SPLIT


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
   LOGICAL, PARAMETER                 :: debug = .FALSE.
   INTEGER, ALLOCATABLE               :: Tab_Iq(:, :) 
   TYPE(Op_t)                         :: Op
   ! --- Begin propagation ---------------------------------------------
   IF (debug) THEN
      WRITE(out_unit,*) 'BEGINNING march_CN ', t,'Δt=', propa%delta_t
      FLUSH(out_unit)
   END IF

   !initialize  Op
   !CALL set_op(Op, psi%basis,Tab_Iq,alloc=.TRUE.)
   CALL Make_Mat_H(psi%basis,Op)

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
!>   - Debugging messages are printed if the local variable 'debug' is set to .TRUE.
!----------------------------------------------------------------------------------

!============================================================
! SUBROUTINE: solve_tdse_cn_gmres
! Purpose: propagate TDSE one time step using CN+GMRES
!============================================================
SUBROUTINE solve_tdse_cn_gmres(psi_n, H_op, psi_np1, dt, max_gmres, gmres_tol)
   IMPLICIT NONE
   TYPE(psi_t),  INTENT(IN)  :: psi_n
   TYPE(op_t),   INTENT(IN)  :: H_op
   TYPE(psi_t),  INTENT(INOUT) :: psi_np1
   REAL(KIND=Rkind),INTENT(IN)  :: dt
   INTEGER,      INTENT(IN)  :: max_gmres
   REAL(KIND=Rkind),INTENT(IN)  :: gmres_tol
   LOGICAL, PARAMETER :: debug=.TRUE.
 
   TYPE(psi_t) :: b
   COMPLEX(KIND=Rkind) :: cn_fac
 
   ! Compute CN factor: i*dt/2
   cn_fac = CMPLX(0.0_Rkind, 0.5_Rkind*dt, kind=Rkind)
 
   ! Initialize RHS b = (I - cn_fac*H)*psi_n
   CALL init_psi(b, psi_n%Basis, cplx=.TRUE., grid=.FALSE.)
   CALL matvec_cn_op(H_op, psi_n, b, -cn_fac)  ! b = psi_n - cn_fac*H*psi_n
 
   ! Initialize guess psi_np1 = 0
   CALL init_psi(psi_np1, psi_n%Basis, cplx=.TRUE., grid=.FALSE.)
   psi_np1%CVec = CZERO
 
   ! Solve (I + cn_fac*H) psi_np1 = b using GMRES
   CALL gmres_solver(H_op, b, psi_np1, max_gmres, gmres_tol, cn_fac)
 
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
   LOGICAL, PARAMETER :: debug=.TRUE.
 
   IF(debug) THEN
      WRITE(out_unit,*) '>>> BEGIN march_CN_With_GMRES at time=', t
      IF(.NOT. ASSOCIATED(psi%Basis)) STOP 'psi%Basis not associated'
      IF(.NOT. ASSOCIATED(psi_dt%Basis)) STOP 'psi_dt%Basis not associated'
   END IF
 
   CALL solve_tdse_cn_gmres(psi, H, psi_dt, propa%delta_t, propa%max_iter, propa%eps)
 
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
   LOGICAL, PARAMETER :: debug=.FALSE.
   
   IF (debug) THEN
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

   IF (debug) THEN
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

   IF (debug) THEN
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
 !   - Builds a Krylov subspace using the Arnoldi process
 !   - Uses Modified Gram–Schmidt to orthogonalize basis vectors
 !   - Applies Givens rotations to maintain a small least-squares problem
 !   - Iteratively improves solution until residual < tol or max_iter reached
 !
 ! ARGUMENTS:
 !   H_op      (IN)    : Hamiltonian operator
 !   psi_b     (IN)    : Right-hand side vector b
 !   psi_x    (INOUT)  : On input initial guess (often zero), on output solution
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
SUBROUTINE gmres_solver(H_op, psi_b, psi_x, max_iter, tol, cn_fac)
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
   LOGICAL, PARAMETER :: debug=.TRUE.
 
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
 
      ! Apply Givens rotations to Hessenberg
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
!   - Initializes guess ψ_{n+1}=0
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
!    J_k is the Bessel function of the first kind and `arg` is scaled
!    time step `delta_t` and spectral radius `delta_E`. The routine
!    stops when the magnitude of the Bessel function drops below a
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
  
   IMPLICIT NONE
   
   ! INPUTS
   REAL(KIND=Rkind), INTENT(IN) :: delta_E
   REAL(KIND=Rkind), INTENT(IN) :: delta_t
   REAL(KIND=Rkind), INTENT(IN) :: limit
   
   ! OUTPUTS
   INTEGER, INTENT(OUT) :: num_coeffs
   COMPLEX(KIND=Rkind), ALLOCATABLE, INTENT(OUT) :: a(:)
   
   ! Local variables-------------------------------
   INTEGER :: k
   REAL(KIND=Rkind) :: arg_bessel
   LOGICAL, PARAMETER :: debug = .FALSE.
   
   ! Calculate the argument for the Bessel function J_k
   arg_bessel = delta_E * delta_t
   
   !------------------------------------------
   ! Calculate coefficients and determine num_coeffs
   !------------------------------------------
   ! First pass to determine the number of coefficients
   ! This is the most proper method to avoid multiple allocations/deallocations
   num_coeffs = 1
   IF (ABS(bessel_jn(num_coeffs, arg_bessel) * TWO) > limit) THEN
      DO WHILE (ABS(bessel_jn(num_coeffs, arg_bessel) * TWO) > limit)
         num_coeffs = num_coeffs + 1
         IF (debug) WRITE(out_unit,*) 'abs(bessel_jn(', num_coeffs, ', ', arg_bessel, ') * TWO) = ', ABS(bessel_jn(num_coeffs, arg_bessel) * TWO)
      END DO
   END IF
   
   IF (debug) WRITE(out_unit,*) 'Chebychev limit reached for coeffs num = ', num_coeffs
   
   ! Allocate array with the correct size
   ALLOCATE(a(num_coeffs))
   
   !------------------------------------------
   ! Calculate final coefficients
   !------------------------------------------
   ! Case k = 0
   a(1) = bessel_jn(0, arg_bessel)
   
   ! Loop for coefficients k > 0
   DO k = 1, num_coeffs-1
     ! Calculer le terme (-i)^k et le coefficient J_k
     SELECT CASE (MOD(k, 4))
       CASE (0)
         a(k+1) = CMPLX(TWO, ZERO, KIND=Rkind) * bessel_jn(k, arg_bessel)
       CASE (1)
         a(k+1) = CMPLX(ZERO, -TWO, KIND=Rkind) * bessel_jn(k, arg_bessel)
       CASE (2)
         a(k+1) = CMPLX(-TWO, ZERO, KIND=Rkind) * bessel_jn(k, arg_bessel)
       CASE (3)
         a(k+1) = CMPLX(ZERO, TWO, KIND=Rkind) * bessel_jn(k, arg_bessel)
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
   LOGICAL, PARAMETER :: debug = .FALSE.
   
   ! Outputs
   TYPE(psi_t), INTENT(INOUT) :: Hprime_psi
   
   ! Local variables
   TYPE(psi_t)                 :: Hpsi
 
   ! Initialize temporary vector Hpsi
   CALL init_psi(Hpsi, psi%Basis, cplx=.TRUE., grid=.FALSE.)
   Hpsi%CVec(:) = CZERO
   Hprime_psi%CVec(:) = CZERO
   
   IF (debug) THEN
      WRITE(out_unit,*) 'BEGINNING Apply_Hprime_To_Psi alpha = ', alpha, ' ΔE = ', delta_E
      FLUSH(out_unit)
   END IF
   ! Étape 1 : Calculer H|psi> en utilisant la routine existante
   ! La routine calc_Oppsi doit être accessible
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
   LOGICAL, PARAMETER                 :: debug = .TRUE.

   ! --- Begin propagation ---------------------------------------------
   IF (debug) THEN
      WRITE(out_unit,*) 'BEGINNING march_Chebychev ', t, ' Δt = ', propa%delta_t
      FLUSH(out_unit)
   END IF

   ! 1. Extract propagation parameters and calculate exponential term
   alpha = HALF*(propa%E_min+propa%E_max)
   delta_E = HALF*(propa%E_max-propa%E_min)
   exp_alpha_term = EXP(-EYE*alpha * propa%delta_t)

   ! 2. Calculate Bessel coefficients a_k
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

   ! 5. Calculate phi_1 = H' * phi_0
   CALL Apply_Hprime_To_Psi(H, phi_0, phi_1, alpha, delta_E)
   psi_dt%CVec(:) = psi_dt%CVec(:) + cheby_coeffs(2) * phi_1%CVec(:) ! Terme a_1 * phi_1

   ! Initialiser les vecteurs pour la boucle de récurrence
   phi_k_minus_1%CVec = phi_0%CVec
   phi_k%CVec         = phi_1%CVec

   ! 6. Loop for k=2 to N, calculating phi_k and accumulating the solution
   IF (num_coeffs > 2) THEN
      DO k = 2, num_coeffs-1
         CALL Cheby_Step(H, phi_k, phi_k_minus_1, phi_k_plus_1, alpha, delta_E)
         ! phi_k_plus_1 now contains T_{k}
         psi_dt%CVec(:) = psi_dt%CVec(:) + cheby_coeffs(k+1) * phi_k_plus_1%CVec(:)

         ! rotate
         phi_k_minus_1%CVec = phi_k%CVec
         phi_k%CVec = phi_k_plus_1%CVec
      END DO
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
   LOGICAL, PARAMETER :: debug = .true.
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

      ! Output progress
      IF (MOD(iter, iprint) == 0) THEN
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
   CHARACTER(LEN=40)            :: propa_name, propa_name2
   INTEGER                      :: max_iter
   LOGICAL                      :: Beta, P, renorm

   NAMELIST /prop/ t0, tf, delta_t, max_iter, eps, E_min, E_max, propa_name, propa_name2, &
                   Beta, P, renorm

   ! Default values
   t0          = ZERO
   tf          = TEN
   delta_t     = ONETENTH**3
   eps         = ONETENTH**10
   E_min       = ZERO
   E_max       = TEN
   max_iter    = 500
   propa_name  = 'non_hagedorn'
   propa_name2 = 'rk4'
   Beta        = .TRUE.
   P           = .TRUE.
   renorm      = .TRUE.

   ! Read the namelist
   READ (*, NML=prop)

   ! Assign to structure
   propa%t0          = t0
   propa%tf          = tf
   propa%delta_t     = delta_t
   propa%eps         = eps
   propa%E_min       = E_min
   propa%E_max       = E_max
   propa%max_iter    = max_iter
   propa%propa_name  = propa_name
   propa%propa_name2 = propa_name2
   propa%Beta        = Beta
   propa%P           = P 
   propa%renorm      = renorm

   ! Convert method names to lowercase
   CALL string_uppercase_TO_lowercase(propa%propa_name)
   CALL string_uppercase_TO_lowercase(propa%propa_name2)

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
   WRITE (out_unit, '(A15, " : ", A)')      'propa_name',  TRIM(propa%propa_name)
   WRITE (out_unit, '(A15, " : ", A)')      'propa_name2', TRIM(propa%propa_name2)
   WRITE (out_unit, '(A15, " : ", L1)')     'Beta',       propa%Beta
   WRITE (out_unit, '(A15, " : ", L1)')     'P',          propa%P
   WRITE (out_unit, '(A15, " : ", L1)')     'renorm',     propa%renorm
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

         Print*,"psi  is on Grid"
         CALL init_psi(psi_b, psi%Basis, cplx=.TRUE., grid=.false.)
         call GridTOBasis_nD_cplx(psi_b%CVec, psi%CVec, psi%Basis)
         CALL init_psi(Hpsi, psi%Basis, cplx=.TRUE., grid=.false.)
         call calc_OpPsi(H, psi_b, Hpsi)
         !stop 'verifie le calcul d E'
         write(out_unit,*) 'rappele toi , sur la calcule de l energie'
         E = real(dot_product(Hpsi%CVec, psi_b%CVec), kind=Rkind)

      else
         Print*,"psi is on basis"
         CALL init_psi(Hpsi, psi%Basis, cplx=.TRUE., grid=.false.)
         call calc_OpPsi(H, psi, Hpsi)
         !stop 'verifie le calcul d E'
         write(out_unit,*) 'rappele toi , sur la calcule de l energie'
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
      debug = .FALSE.  ! Set to .TRUE. to print E
   
      ! ===== Expectation value calculation =====
      IF (psi%Grid) THEN
         ! Convert wavefunction from grid to basis representation
         CALL init_psi(psi_b, psi%Basis, cplx=.TRUE., grid=.FALSE.)
         CALL GridTOBasis_nD_cplx(psi_b%CVec, psi%CVec, psi%Basis)
   
         ! Apply operator H to wavefunction in basis
         CALL init_psi(Hpsi, psi%Basis, cplx=.TRUE., grid=.FALSE.)
         CALL calc_OpPsi(H, psi_b, Hpsi)
   
         ! Compute <psi_b | H | psi_b>
         E = REAL(DOT_PRODUCT(Hpsi%CVec, psi_b%CVec), KIND=Rkind)
      ELSE
         ! Directly apply H to wavefunction
         CALL init_psi(Hpsi, psi%Basis, cplx=.TRUE., grid=.FALSE.)
         CALL calc_OpPsi(H, psi, Hpsi)
   
         ! Compute <psi | H | psi>
         E = REAL(DOT_PRODUCT(Hpsi%CVec, psi%CVec), KIND=Rkind)
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
!> @param[in] nio     Numéro d'unité logique (unit) pour le fichier.
!> @param[in] name    Nom de base à utiliser pour créer le fichier.
!> @param[in] propa   Structure contenant les informations sur la
!>                    propagation (nom, type, pas de temps, etc.).
!>
!> @details
!> Le nom final du fichier est construit comme :
!>   name_propaName_propaName2.txt
!> Exemple :
!>   si name = "Qt", propa%propa_name = "STD", propa%propa_name2 = "nb10"
!>   alors le fichier sera : "Qt_STD_nb10.txt"
!>
!> Le fichier est ouvert avec le buffering désactivé (`BUFFERED='NO'`)
!> pour que chaque écriture soit immédiatement visible sur disque.
!===============================================================
SUBROUTINE creat_file_unit(nio, name, propa)

   !==== Arguments ====
   INTEGER,          INTENT(IN) :: nio       ! Numéro d'unité du fichier
   CHARACTER(*),     INTENT(IN) :: name      ! Nom de base du fichier
   TYPE(propa_t),    INTENT(IN) :: propa     ! Données de la propagation

   !==== Variables locales ====
   CHARACTER(100) :: name_tot   ! Nom complet du fichier
   CHARACTER(20)  :: dt         ! Pas de temps converti en chaîne
   CHARACTER(8)   :: fmt        ! Format pour afficher le pas de temps

   !==== Construction du format pour écrire delta_t ====
   fmt = "(E0.1)"
   WRITE(dt, fmt) propa%delta_t

   !==== Construction du nom complet du fichier ====
   name_tot = TRIM(name) // '_' // TRIM(propa%propa_name) // '_' // &
              TRIM(propa%propa_name2) // '.txt'

   !==== Ouverture du fichier ====
   ! BUFFERED='NO' → désactive le buffering pour écriture immédiate
   OPEN(UNIT=nio, FILE=name_tot, STATUS='UNKNOWN', ACTION='WRITE', &
        FORM='FORMATTED')

END SUBROUTINE creat_file_unit

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

      INTEGER :: i, n
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
   LOGICAL, parameter                 :: debug=.TRUE.

   ! --- Begin propagation ---
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