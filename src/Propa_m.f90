module Propa_m
   USE QDUtil_m
   USE psi_m
   USE Op_m
   USE Ana_psi_m
   Use lanczos_m
   USE Auto_corr_m
   Use Hagedorn_m
   USE sub_propa_m
   USE Basis_m

   implicit none

contains



SUBROUTINE propagation(psif, psi0, propa, H)
   USE psi_m
   USE Basis_m
   IMPLICIT NONE

   !====================== Interface ==========================
   TYPE(psi_t),     INTENT(INOUT) :: psif     ! Final wavefunction after propagation
   TYPE(psi_t),     INTENT(IN)    :: psi0     ! Initial wavefunction
   TYPE(propa_t),   INTENT(INOUT) :: propa    ! Propagation parameters
   TYPE(Op_t),      INTENT(INOUT) :: H        ! Hamiltonian initialized by the caller

   !====================== Local constants ====================
   LOGICAL, PARAMETER :: debug = .true. ,ldebug = .false.      ! Debug flag

   !====================== Local variables ====================
   TYPE(Basis_t)              :: Basis           ! Local basis
   REAL(kind=Rkind)           :: t               ! Current time
   REAL(kind=Rkind)           :: t_deltat        ! Time at next step
   REAL(kind=Rkind)           :: Norm            ! Wavefunction norm
   REAL(kind=Rkind)           :: E, E0           ! Energies
   REAL(kind=Rkind)           :: aut_func_arg(1) ! Autocorrelation phase argument
   COMPLEX(kind=Rkind)        :: aut_func(1)     ! Autocorrelation value

   REAL(kind=Rkind), ALLOCATABLE    :: Qt(:)     ! Position centers
   REAL(kind=Rkind), ALLOCATABLE    :: SQt(:)    ! Widths
   REAL(kind=Rkind), ALLOCATABLE    :: Pt(:)     ! Momentum centers
   REAL(kind=Rkind), ALLOCATABLE    :: Qt_ref(:) ! Reference position centers for validation
   REAL(kind=Rkind), ALLOCATABLE    :: Pt_ref(:) ! Reference momentum centers for validation
   COMPLEX(kind=Rkind), ALLOCATABLE :: At(:)     ! Complex widths
   REAL(kind=Rkind), ALLOCATABLE    :: pop(:)    ! Populations per surface
   INTEGER, ALLOCATABLE             :: Tab_iq(:, :)  ! Basis index mapping
   INTEGER :: u_psi, u_Qt, u_E, u_SQt, u_Norm, u_pop, u_Imp_k, u_alpha
   INTEGER :: u_Norm_13, u_E_13, u_auto_cor, u_validation, u_psi_dt

   REAL(kind=Rkind)           :: Norm_ref, E_ref
   REAL(kind=Rkind)           :: max_norm_drift, max_energy_drift
   REAL(kind=Rkind)           :: max_abs_dQ, max_abs_dP, auto_mod_dev
   REAL(kind=Rkind)           :: step_max_abs_dQ, step_max_abs_dP
   REAL(kind=Rkind), PARAMETER :: validation_tol_norm   =  ONETENTH**10
   REAL(kind=Rkind), PARAMETER :: validation_tol_energy = ONETENTH**10

   !TYPE(REDUCED_DENSIRY_t)    :: Rd               ! Reduced density matrix (unused for now)
    
   INTEGER :: Ndim, i, nt, nsurf              ! Dimensions, loop indices
   TYPE(psi_t) :: psi, psi_dt, psi_t0 , Gspsi     ! Wavefunction containers (current, step, and reference)
   character(len=Name_len) :: name

   !====================== Initialization =====================


   name = propa%propa_name2 
   call string_uppercase_TO_lowercase(name)

   IF (debug) THEN
      WRITE(out_unit, *) 'BEGINNING propagation', propa%t0, propa%tf, propa%delta_t
      WRITE(out_unit, *) '-------------propagation parameters---------------'
      CALL write_propa(propa)
   ELSE
      STOP ' check your data!'
      FLUSH(out_unit)
   END IF

   ! Create output files
   if(ldebug) CALL open_propa_file(u_psi, 'psi')
   CALL open_propa_file(u_Qt, 'Qt')
   CALL open_propa_file(u_E, 'E')
   CALL open_propa_file(u_SQt, 'SQt')
   CALL open_propa_file(u_Norm, 'Norm')
   CALL open_propa_file(u_pop, 'pop')
   CALL open_propa_file(u_Imp_k, 'Imp_k')
   CALL open_propa_file(u_alpha, 'alpha')
   CALL open_propa_file(u_Norm_13, 'Norm_13')
   CALL open_propa_file(u_E_13, 'E_13')
   CALL open_propa_file(u_auto_cor, 'auto_cor')
   if(ldebug) CALL open_propa_file(u_psi_dt, 'psi_dt_on_basis0')

   ! Initialize basis from psi0
   CALL init_Basis1_TO_Basis2(Basis, psi0%Basis)
   CALL construct_primitive_basis(Basis)

   ! Get number of degrees of freedom and electronic surfaces
   Ndim = SIZE(psi0%Basis%tab_basis) - 1
   nsurf = psi0%Basis%tab_basis(Ndim + 1)%nb

   ! Allocate Gaussian parameters and surface populations
   ALLOCATE(Qt(Ndim), SQt(Ndim), Pt(Ndim), Qt_ref(Ndim), Pt_ref(Ndim), At(Ndim))
   ALLOCATE(pop(nsurf))

   ! Number of propagation steps
   nt = INT((propa%tf - propa%t0) / propa%delta_t)

   ! Initialize wavefunction containers
   CALL init_psi(psi, psi0%Basis, cplx=.true., grid=.false.)
   CALL init_psi(psi_dt, psi0%Basis, cplx=.true., grid=.false.)
   CALL init_psi(psi_t0, Basis, cplx=.true., grid=.false.)
   CALL init_psi(Gspsi, Basis, cplx=.true., grid=.false.)

   ! Copy initial wavefunction into psi and psi_t0
   psi%CVec(:) = psi0%CVec(:)
   psi_t0%CVec(:) = psi0%CVec(:)

   !call Calc_reduced_density(Rd, psi%CVec, psi%Basis)
   !call Rdensity_Writing(Rd, psi%Basis, nio=21, t=ZERO)
   !STOP 'cc propa'
   IF (propa%propa_name2 == 'vp') THEN
     CALL Get_Basis_Parameters(psi%Basis, Qt_ref, SQt, At, Pt_ref)
   ELSE
     CALL Calc_Basis_parameters(psi, Qt_ref, SQt, At, Pt_ref)
   END IF
   CALL Calc_Av_E(E_ref, psi, H)
   CALL Calc_Norm_OF_Psi(psi, Norm_ref)
   max_norm_drift   = ZERO
   max_energy_drift = ZERO
   max_abs_dQ       = ZERO
   max_abs_dP       = ZERO
   CALL open_propa_file(u_validation, 'validation')
   WRITE(u_validation,'(a)') '# t  norm  dnorm  energy  dE  max_abs_dQ  max_abs_dP  abs_auto_minus_one'

   !====================== Begin Propagation ===================

   !call Hagedorn_temp(psi, psi0, propa)
   !CALL Calc_Av_E(e0, psi, H)
   !STOP 'cc propa'

   IF (name == 'itp') THEN
      CALL relax_propagation(psi_t0, H, propa)
   ELSE IF (name == 'gssil') THEN
      CALL find_ground_state_SIL(E0, Gspsi, H, psi0, psi0%Basis, propa%max_iter, propa%eps) 
      WRITE(out_unit, *) "Ground state energy found: ", E0
      CALL Calc_Av_E(E0, Gspsi, H)
      WRITE(out_unit,*) ': <Gspsi|H|Gspsi>  =: ', E0
   ELSE

   DO i = 0, nt
      t = i * propa%delta_t
      t_deltat = t + propa%delta_t

      !WRITE(out_unit, *) propa%propa_name2, i, t, t_deltat

      IF (propa%propa_name2 == 'vp') THEN
         CALL Get_Basis_Parameters(psi%Basis, Qt, SQt, At, Pt)
      ELSE
         CALL Calc_Basis_parameters(psi, Qt, SQt, At, Pt)
      END IF

      CALL Calc_Av_E(E, psi, H)
      CALL Calc_Norm_OF_psi(psi, Norm)
      CALL Population(psi, pop)

      CALL Calc_Auto_corr(psi_t0, psi, aut_func(1), aut_func_arg(1), propa, t=t,debug=ldebug)
      CALL Write_Complex(aut_func, nio=u_auto_cor, op='t_z', t=t)
      CALL Write_Complex(At, nio=u_alpha, op='t_z', t=t)
      step_max_abs_dQ = MAXVAL(ABS(Qt - Qt_ref))
      step_max_abs_dP = MAXVAL(ABS(Pt - Pt_ref))
      max_norm_drift   = MAX(max_norm_drift,   ABS(Norm - Norm_ref))
      max_energy_drift = MAX(max_energy_drift, ABS(E - E_ref))
      max_abs_dQ       = MAX(max_abs_dQ, step_max_abs_dQ)
      max_abs_dP       = MAX(max_abs_dP, step_max_abs_dP)
      auto_mod_dev     = ABS(ABS(aut_func(1)) - ONE)
      WRITE(u_validation,'(8(1x,es24.16))') t, Norm, ABS(Norm - Norm_ref), E, ABS(E - E_ref), &
                              step_max_abs_dQ, step_max_abs_dP, auto_mod_dev

      WRITE(u_Qt,*) t, Qt
      WRITE(u_E,*) t, E
      WRITE(u_SQt,*) t, SQt
      WRITE(u_Norm,*) t, Norm
      WRITE(u_pop,*) t, pop
      WRITE(u_Imp_k,*) t, Pt
      

      FLUSH(u_Qt)
      FLUSH(u_E)
      FLUSH(u_SQt)
      FLUSH(u_Norm)
      FLUSH(u_pop)
      FLUSH(u_Imp_k)
      FLUSH(u_alpha)
      FLUSH(u_auto_cor)
      FLUSH(u_validation)

      !call eval_pics(psi, ib=28, t=t_deltat)

      if (mod(i, 60) == 0) then
         if(ldebug) then 
            call write_psi(psi=psi, psi_cplx=.true., print_psi_grid=.false., &
               print_basis=.false., t=t_deltat, int_print=10, real_part=.false.)
            !call eval_pics(psi, ib=28, t=t)
            !call Calc_reduced_density(Rd, psi%CVec, psi%Basis)
            !call Rdensity_Writing(Rd, psi%Basis, nio=21, ib=1, t=t_deltat)
            !call test_propa_Hagedorn(psi, t)
         WRITE(u_psi,*)
         FLUSH(u_psi)
         !write(21,*)
         !flush(21)
         end if
      end if
       !write(filename, '("output_", I0, ".dat")') i
       !open(unit=100+i, file=filename) 
      !call write_psi(psi=psi, psi_cplx=.false., print_psi_grid=.true., &
      !  print_basis=.false., t=t_deltat, int_print=i+100, real_part=.false.)
      !  close(100+i)

      CALL march_Global(psi, psi_dt, t_deltat, propa, H, u_Norm_13, u_E_13)

      IF (propa%propa_name == 'hagedorn') THEN
         !> This part reconstructs the potential after the basis has been updated.
         IF (ALLOCATED(H%Scal_pot)) DEALLOCATE(H%Scal_pot)
         CALL Set_Op(H, psi%Basis, Tab_Iq,nio=in_unit)
      END IF
   END DO

   !====================== Finalization ========================

   IF (propa%propa_name == 'hagedorn') THEN
      psif%CVec(:) = psi%CVec(:)
   ELSE
      psif%CVec(:) = psi_dt%CVec(:)
   END IF

   CALL Calc_Norm_OF_Psi(psif, Norm)

   WRITE(u_validation,'(a)') '# summary'
   WRITE(u_validation,'(a,1x,es24.16)') 'reference_norm', Norm_ref
   WRITE(u_validation,'(a,1x,es24.16)') 'reference_energy', E_ref
   WRITE(u_validation,'(a,1x,es24.16)') 'max_norm_drift', max_norm_drift
   WRITE(u_validation,'(a,1x,es24.16)') 'max_energy_drift', max_energy_drift
   WRITE(u_validation,'(a,1x,es24.16)') 'max_abs_dQ', max_abs_dQ
   WRITE(u_validation,'(a,1x,es24.16)') 'max_abs_dP', max_abs_dP
   WRITE(u_validation,'(a,1x,l1)') 'norm_within_tolerance', max_norm_drift <= validation_tol_norm
   WRITE(u_validation,'(a,1x,l1)') 'energy_within_tolerance', max_energy_drift <= validation_tol_energy
   CLOSE(u_validation)

   IF (debug) THEN
      WRITE(out_unit, *) 'END propagation'
      WRITE(out_unit, *) 'norm, psi_dt', Norm
      WRITE(out_unit, *) 'validation max_norm_drift  ', max_norm_drift
      WRITE(out_unit, *) 'validation max_energy_drift', max_energy_drift
      WRITE(out_unit, *) 'validation max_abs_dQ      ', max_abs_dQ
      WRITE(out_unit, *) 'validation max_abs_dP      ', max_abs_dP
      FLUSH(out_unit)
   END IF

   END IF
END SUBROUTINE propagation


 End module Propa_m
