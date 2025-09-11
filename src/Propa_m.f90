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



SUBROUTINE propagation(psif, psi0, propa)
   USE psi_m
   USE Basis_m
   IMPLICIT NONE

   !====================== Interface ==========================
   TYPE(psi_t),     INTENT(INOUT) :: psif     ! Final wavefunction after propagation
   TYPE(psi_t),     INTENT(IN)    :: psi0     ! Initial wavefunction
   TYPE(propa_t),   INTENT(INOUT) :: propa    ! Propagation parameters

   !====================== Local constants ====================
   LOGICAL, PARAMETER :: debug = .true. ,ldebug = .false.      ! Debug flag

   !====================== Local variables ====================
   TYPE(Basis_t)              :: Basis           ! Local basis
   TYPE(Op_t)                 :: H               ! Hamiltonian operator
   REAL(kind=Rkind)           :: t               ! Current time
   REAL(kind=Rkind)           :: t_deltat        ! Time at next step
   REAL(kind=Rkind)           :: Norm            ! Wavefunction norm
   REAL(kind=Rkind)           :: E, E0           ! Energies
   REAL(kind=Rkind)           :: aut_func_arg(1) ! Autocorrelation phase argument
   COMPLEX(kind=Rkind)        :: aut_func(1)     ! Autocorrelation value

   REAL(kind=Rkind), ALLOCATABLE    :: Qt(:)     ! Position centers
   REAL(kind=Rkind), ALLOCATABLE    :: SQt(:)    ! Widths
   REAL(kind=Rkind), ALLOCATABLE    :: Pt(:)     ! Momentum centers
   COMPLEX(kind=Rkind), ALLOCATABLE :: At(:)     ! Complex widths
   REAL(kind=Rkind), ALLOCATABLE    :: pop(:)    ! Populations per surface
   INTEGER, ALLOCATABLE             :: Tab_iq(:, :)  ! Basis index mapping

   TYPE(REDUCED_DENSIRY_t)    :: Rd               ! Reduced density matrix (unused for now)
    
   INTEGER :: Ndim, i, nt, nf, nsurf              ! Dimensions, loop indices
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
   if(ldebug) CALL creat_file_unit(nio=10, name='psi', propa=propa)
   CALL creat_file_unit(nio=11, name='Qt', propa=propa)
   CALL creat_file_unit(nio=12, name='E', propa=propa)
   CALL creat_file_unit(nio=13, name='SQt', propa=propa)
   CALL creat_file_unit(nio=14, name='Norm', propa=propa)
   CALL creat_file_unit(nio=18, name='pop', propa=propa)
   CALL creat_file_unit(nio=19, name='Imp_k', propa=propa)
   CALL creat_file_unit(nio=20, name='alpha', propa=propa)
   !CALL creat_file_unit(nio=21, name='Rd', propa=propa)
   !CALL creat_file_unit(nio=22, name='maxcoeff', propa=propa)
   !CALL creat_file_unit(nio=23, name='psi_int', propa=propa)
   CALL creat_file_unit(nio=24, name='Norm_13', propa=propa)
   CALL creat_file_unit(nio=25, name='E_13', propa=propa)
   CALL creat_file_unit(nio=26, name='auto_cor', propa=propa)
   if(ldebug) CALL creat_file_unit(nio=27, name='psi_dt_on_basis0', propa=propa)
   !CALL creat_file_unit(nio=28, name='file_norm_pics', propa=propa)

   ! Initialize basis from psi0
   CALL init_Basis1_TO_Basis2(Basis, psi0%Basis)
   CALL construct_primitive_basis(Basis)

   ! Get number of degrees of freedom and electronic surfaces
   Ndim = SIZE(psi0%Basis%tab_basis) - 1
   nsurf = psi0%Basis%tab_basis(Ndim + 1)%nb

   ! Allocate Gaussian parameters and surface populations
   ALLOCATE(Qt(Ndim), SQt(Ndim), Pt(Ndim), At(Ndim))
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

   ! Set up operator
   CALL Calc_tab_Iq0(Tab_Iq, psi0%Basis)
   CALL Set_Op(H, psi0%Basis, Tab_Iq)

   !call Calc_reduced_density(Rd, psi%CVec, psi%Basis)
   !call Rdensity_Writing(Rd, psi%Basis, nio=21, t=ZERO)
   !STOP 'cc propa'

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
      CALL Write_Complex(aut_func, nio=26, op='t_z', t=t)
      CALL Write_Complex(At, nio=20, op='t_z', t=t)

      WRITE(11,*) t, Qt
      WRITE(12,*) t, E
      WRITE(13,*) t, SQt
      WRITE(14,*) t, Norm
      WRITE(18,*) t, pop
      WRITE(19,*) t, Pt
      

      FLUSH(11)
      FLUSH(12)
      FLUSH(13)
      FLUSH(14)
      FLUSH(18)
      FLUSH(19)
      FLUSH(20)
      FLUSH(26)

      !call eval_pics(psi, ib=28, t=t_deltat)

      if (mod(i, 60) == 0) then
         if(ldebug) then 
            call write_psi(psi=psi, psi_cplx=.true., print_psi_grid=.false., &
               print_basis=.false., t=t_deltat, int_print=10, real_part=.false.)
            !call eval_pics(psi, ib=28, t=t)
            !call Calc_reduced_density(Rd, psi%CVec, psi%Basis)
            !call Rdensity_Writing(Rd, psi%Basis, nio=21, ib=1, t=t_deltat)
            !call test_propa_Hagedorn(psi, t)
         write(10,*)
         flush(10)
         !write(21,*)
         !flush(21)
         end if
      end if

      CALL march_Global(psi, psi_dt, t_deltat, propa, H)

      IF (propa%propa_name == 'hagedorn') THEN
         !> This part reconstructs the potential after the basis has been updated.
         DEALLOCATE(H%Scal_pot)
         CALL Set_Op(H, psi%Basis, Tab_Iq)
      END IF
   END DO

   !====================== Finalization ========================

   psif%CVec(:) = psi_dt%CVec(:)
   CALL Calc_Norm_OF_Psi(psif, Norm)

   IF (debug) THEN
      WRITE(out_unit, *) 'END propagation'
      WRITE(out_unit, *) 'norm, psi_dt', Norm
      FLUSH(out_unit)
   END IF

   END IF
END SUBROUTINE propagation


 End module Propa_m
