PROGRAM TD_SCHROD
   USE QDUtil_m
   USE Basis_m
   USE psi_m
   USE op_m
   USE InputStaging_m, ONLY: Stage_Standard_Input
   USE Initial_GWP_Basis_m, ONLY: Apply_Initial_GWP_To_Basis
   USE HamiltonianSpectrum_m, ONLY: H_spectrum
   USE TISE_m, ONLY: Solve_TISE, Relax_TISE
   USE Param_WP0_m
   USE Propa_m
   USE Ana_psi_m
   USE lanczos_m
   USE Sub_Vp_m
   IMPLICIT NONE

   !======================= Main Variables ===========================
   TYPE(Basis_t), TARGET         :: Basis        ! Basis object
   TYPE(Op_t)                    :: H !,Op           ! Hamiltonian operator
   TYPE(psi_t)                   :: psi0, psif   ! Initial and final wavefunctions
   TYPE(propa_t)                :: propa        ! Propagation parameters
   TYPE(GWP_t), ALLOCATABLE     :: tab_GWP(:)   ! Gaussian wavepacket parameters
   INTEGER, ALLOCATABLE         :: Tab_iq(:, :) ! Index mapping for operator setup

   ! Variables for initial WP parameters
   REAL(Kind=Rkind) :: E, Norm
   REAL(Kind=Rkind), ALLOCATABLE :: qmodel_q0(:), qmodel_sigma(:)
   INTEGER :: idim

   CALL Reset_Debug_Controls()
   CALL Stage_Standard_Input()

   !======================= QML and Potential Info ====================
    ! Initialize QML potential model
   call InitQModelBasisParam(nio=in_unit, qmodel_q0=qmodel_q0, qmodel_sigma=qmodel_sigma)
   !CALL sub_Init_Qmodel(ndim, nsurf, pot_name, adiabatic, option)
   !===================================================================
   ! Read basis set parameters from input file and construct basis
   CALL Read_Basis(Basis, nio=in_unit)

   ! Read the initial Gaussian once, then optionally derive the primitive
   ! HO/Hagedorn basis parameters from that same definition.
   CALL Read_tab_GWP(tab_GWP=tab_GWP, nb_GWP=1, nio=in_unit)
   IF (TRIM(tab_GWP(1)%source) == 'qmodel') THEN
      IF (.NOT. ALLOCATED(qmodel_q0) .OR. .NOT. ALLOCATED(qmodel_sigma)) &
         ERROR STOP 'QML Gaussian parameters were not returned by InitQModelBasisParam.'
      IF (SIZE(qmodel_q0) /= tab_GWP(1)%ndim .OR. SIZE(qmodel_sigma) /= tab_GWP(1)%ndim) &
         ERROR STOP 'QML and initial GWP dimensions do not match.'
      DO idim = 1, tab_GWP(1)%ndim
         tab_GWP(1)%tab_GWP1D(idim)%Q0 = qmodel_q0(idim)
         tab_GWP(1)%tab_GWP1D(idim)%sigma = qmodel_sigma(idim)
      END DO
      IF (Basis_Debug_Active() .OR. Debug_Enabled()) &
         WRITE(out_unit,*) 'Initial Gaussian centers and widths derived from InitQModelBasisParam.'
   END IF
   ! Read propagation controls now, before basis/operator initialization, so
   ! propa%debug also enables diagnostics in those initialization routines.
   CALL read_propa(propa)
   IF (propa%debug .OR. Basis_Debug_Active()) THEN
      WRITE(out_unit,*) 'Initial Gaussian source/dimension: ', TRIM(tab_GWP(1)%source), tab_GWP(1)%ndim
      IF (ALLOCATED(qmodel_q0)) WRITE(out_unit,*) 'QML centers: ', qmodel_q0
      IF (ALLOCATED(qmodel_sigma)) WRITE(out_unit,*) 'QML widths:  ', qmodel_sigma
      DO idim = 1, tab_GWP(1)%ndim
         WRITE(out_unit,*) 'Initial GWP dim, Q0, sigma, Beta, Imp_k, gamma = ', idim, &
              tab_GWP(1)%tab_GWP1D(idim)%Q0, tab_GWP(1)%tab_GWP1D(idim)%sigma, &
              tab_GWP(1)%tab_GWP1D(idim)%Beta, tab_GWP(1)%tab_GWP1D(idim)%imp_k, &
              tab_GWP(1)%tab_GWP1D(idim)%gamma
      END DO
   END IF
   CALL Apply_Initial_GWP_To_Basis(Basis, tab_GWP(1), propa%debug)
   IF (propa%debug .OR. Basis_Debug_Active()) THEN
      DO idim = 1, SIZE(Basis%tab_basis) - 1
         WRITE(out_unit,*) 'Basis dim, Q0, scaleQ, Imp_k, Alpha = ', idim, &
              Basis%tab_basis(idim)%Q0, Basis%tab_basis(idim)%scaleQ, &
              Basis%tab_basis(idim)%Imp_k, Basis%tab_basis(idim)%Alpha
      END DO
   END IF
   CALL construct_primitive_basis(Basis)
   !nb = Basis%nb
   !CALL Write_Basis(Basis) 
   !STOP 'cc'

   !WRITE(*,*) "Basis is allocated", Basis_IS_allocated(Basis)

   IF (propa%debug .OR. Basis_Debug_Active()) &
      WRITE(out_unit, *) "------------------- Initialization of  psi0 ---------------------"

   ! Initialize wavefunctions
   CALL init_psi(psi0, Basis, cplx=.true., grid=.false.)
   CALL init_psi(psif, Basis, cplx=.true., grid=.false.)

   !CALL test_basitogridgridtobasis(Basis)

   ! Construct operator H
   CALL Calc_tab_Iq0(Tab_Iq,Basis)
   CALL Set_Op(H, Basis, Tab_Iq,alloc=.FALSE.,nio=in_unit)
   ! Always project the requested Gaussian onto the chosen basis. The basis may
   ! intentionally differ from the packet (especially for Hagedorn propagation).
   CALL psi_init_GWP(psi=psi0, Tab_GWP=tab_GWP)
   !CALL Write_VecMat(psi0%CVec, out_unit, 5, info='psi0%CVec')

   !CALL Save_Potential_ToFile(H%Scal_pot,Basis)
   !CALL H_spectrum(H=H,N_states=Basis%nb,psi=psif)
   !CALL Test_Calc_S()
   !CALL Ecrire_psi(psi0, nio=100, t=ZERO)
   !PRINT*, "psi0", psi0%CVec(:)
   !psi0%CVec = CZERO
   !psi0%CVec(1) = CONE

   !CALL Calc_Basis_parameters(psi0, Qt, SQt, At, Pt)
   !CALL Get_Basis_Parameters(Basis, Qt, SQt, At, Pt)
   !CALL Calc_Avg_A_nD(psi0, At)

   ! Compute initial energy and norm
   IF (propa%debug .OR. Basis_Debug_Active()) THEN
      CALL Calc_Av_E(E, psi0, H)
      CALL Calc_Norm_OF_Psi(psi0, Norm)
      WRITE(out_unit, *) '-------------Energy And Norm of initial WP0-----------------------------'
      WRITE(out_unit, *) ' <psi|H|psi> ', E, ' <psi|psi> ', Norm
      WRITE(out_unit, *) '-------------End Initialization of psi0 -------------------------------'
   END IF

   IF (Basis_UsesSmolyakSparseGrid(Basis)) THEN
      CALL Validate_Smolyak_Implementation(H, psi0)
   END IF

   !CALL Test_projection_hagedorn(Basis, psi0)
   !CALL test_op(Basis, psi0)  ! Diagonalization test
   !CALL Calc_Scalar_Pot(V, Basis)
   !CALL test_openmp_op(Basis)
   !CALL calc_nac(Basis)
   !CALL calc_VV(Basis)

   ! Read propagation parameters
   IF (TRIM(propa%dynamics) == 'tise') THEN
      WRITE(out_unit,*) '===== Stationary Schrodinger equation (TISE) ====='
      IF (propa%tise_relax) CALL Relax_TISE(H, psi0, propa%tise_relax_delta_t, &
                                            propa%tise_relax_max_iter, propa%tise_relax_tol, propa%debug)
      CALL Solve_TISE(H, propa%tise_method, propa%tise_nstates, propa%tise_max_iter, &
                     propa%tise_tol, propa%energy_unit, propa%tise_energy_file, &
                     propa%tise_transition_file, psi0, propa%debug)
      WRITE(out_unit,*) '===== End stationary calculation ====='
      CALL dealloc_psi(psi0)
      CALL dealloc_psi(psif)
      STOP
   ELSEIF(propa%writeV)THEN
      call Write_Potential(Basis)
      Write(out_unit,*) 'Write Potential completed'
   END IF
   IF (propa%action_H) CALL H_action_diagnostic(H, psi0)
   !CALL Print_Grids(Basis)
   !CALL Test_NUFFT_Cycle(psi0%CVec,Basis)

   !call Test_Compare_Overlap_ThreeMethods()
   !CALL Test_Calc_M22(Basis)
   !CALL Verify_M22_2D(Basis)
   !CALL Test_All_Derivatives(psi0)
   !CALL Test_Matrices(psi0,H)
   !CALL Test_Psi_Lambda_Transform(psi0)
   !CALL test_taylor(psi0, psif, propa, H)
   !STOP 'calcul de H|psi> est fait'

   !====================== Time Propagation ============================
   CALL propagation(psif, psi0, propa, H)
   WRITE(out_unit, *)  'Fin de la propagation'

   !CALL Write_psi(psif)

   !====================== Finalization ================================
   WRITE(out_unit, *) 'deallocation'
   CALL dealloc_psi(psi0)
   CALL dealloc_psi(psif)
END PROGRAM TD_SCHROD
