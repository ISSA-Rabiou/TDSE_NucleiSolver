PROGRAM TD_SCHROD
   USE QDUtil_m
   USE Basis_m
   USE psi_m
   USE op_m
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

   CALL Stage_Standard_Input()

   !======================= QML and Potential Info ====================
    ! Initialize QML potential model
   call  InitQModelBasisParam(nio=in_unit)
   !CALL sub_Init_Qmodel(ndim, nsurf, pot_name, adiabatic, option)
   !===================================================================
   ! Read basis set parameters from input file and construct basis
   CALL Read_Basis(Basis, nio=in_unit)
   CALL construct_primitive_basis(Basis)
   !nb = Basis%nb
   !CALL Write_Basis(Basis) 
   !STOP 'cc'

   !WRITE(*,*) "Basis is allocated", Basis_IS_allocated(Basis)

   WRITE(out_unit, *) "------------------- Initialization of  psi0 ---------------------"

   ! Initialize wavefunctions
   CALL init_psi(psi0, Basis, cplx=.true., grid=.false.)
   CALL init_psi(psif, Basis, cplx=.true., grid=.false.)

   ! Read Gaussian wavepacket parameters from file
   CALL Read_tab_GWP(tab_GWP=tab_GWP, nb_GWP=1, nio=in_unit)

   !CALL test_basitogridgridtobasis(Basis)

   ! Construct operator H
   CALL Calc_tab_Iq0(Tab_Iq,Basis)
   CALL Set_Op(H, Basis, Tab_Iq,alloc=.TRUE.,nio=in_unit)
   ! Initialize psi0 using the GWP parameters
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
   CALL Calc_Av_E(E, psi0, H)
   CALL Calc_Norm_OF_Psi(psi0, Norm)

   WRITE(out_unit, *) '-------------Energy And Norm of initial WP0-----------------------------'
   WRITE(out_unit, *) ' <psi|H|psi> ', E, ' <psi|psi> ', Norm
   WRITE(out_unit, *) '-------------End Initialization of psi0 -------------------------------'

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
   CALL read_propa(propa)
   IF (propa%diag_H) THEN
      WRITE(out_unit,*) '===== Hamiltonian diagnostic requested ====='
      CALL H_spectrum(H, MIN(10, Get_TotalBasisCoefficientCount(Basis)))
      WRITE(out_unit,*) '===== End Hamiltonian diagnostic ====='
   END IF
   IF (propa%action_H) CALL H_action_diagnostic(H, psi0)
   !CALL Make_Mat_H(Op)
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
