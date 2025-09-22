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
   !REAL(Kind=Rkind), ALLOCATABLE :: V(:, :, :)  ! Scalar potential (optional)

   ! Variables for initial WP parameters
   REAL(Kind=Rkind) :: E, Norm
   !REAL(Kind=Rkind) :: Pt(2), Qt(2), SQt(2)
   !COMPLEX(Kind=Rkind) :: At(2)
   !REAL(Kind=Rkind) :: t_deltat = 0

   !======================= QML and Potential Info ====================
   INTEGER :: ndim, nsurf, option
   LOGICAL :: adiabatic
   CHARACTER(len=16) :: pot_name

   !======================= Initialization ============================
   ndim = -1
   nsurf = -1
   pot_name = 'read_model'
   adiabatic = .false.
   option = 1

   ! Initialize QML potential model
   CALL sub_Init_Qmodel(ndim, nsurf, pot_name, adiabatic, option)
   WRITE(out_unit, *) 'ndim, nsurf', ndim, nsurf
   WRITE(out_unit, *) 'pot_name'
   

   !===================================================================
   ! Read basis set parameters from input file and construct basis
   CALL Read_Basis(Basis, nio=in_unit)
   CALL construct_primitive_basis(Basis)
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
   CALL Calc_tab_Iq0(Tab_Iq, psi0%Basis)
   !stop 'cc'
   CALL Set_Op(H, psi0%Basis, Tab_Iq,alloc=.TRUE.)

   ! Initialize psi0 using the GWP parameters
   CALL psi_init_GWP(psi=psi0, Tab_GWP=tab_GWP)

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

   !CALL Test_projection_hagedorn(Basis, psi0)
   !CALL test_op(Basis, psi0)  ! Diagonalization test
   !CALL Calc_Scalar_Pot(V, Basis)
   !CALL test_openmp_op(Basis)
   !CALL calc_nac(Basis)
   !CALL calc_VV(Basis)

   ! Read propagation parameters
   CALL read_propa(propa)
   !CALL Make_Mat_H(Basis,Op)

   !CALL Print_Grids(Basis)
   !CALL Test_NUFFT_Cycle(psi0%CVec,Basis)

   !CALL Test_Calc_M22(Basis)
   !CALL Verify_M22_2D(Basis)
   !CALL Test_All_Derivatives(psi0)
   !CALL Test_Matrices(psi0,H)
   !CALL Test_Psi_Lambda_Transform(psi0)
   !CALL test_taylor(psi0, psif, propa, H)
   !STOP 'calcul de H|psi> est fait'

   !====================== Time Propagation ============================
   CALL propagation(psif, psi0, propa)
   PRINT*, 'Fin de la propagation'

   !CALL Write_psi(psif)

   !====================== Finalization ================================
   WRITE(out_unit, *) 'deallocation'
   CALL dealloc_psi(psi0)
   CALL dealloc_psi(psif)

END PROGRAM TD_SCHROD