module Op_m
   !$ USE omp_lib
   USE QDUtil_m
   USE Basis_m, ONLY: Basis_t, Basis_UsesSmolyakSparseGrid, Get_TotalBasisCoefficientCount, &
     BasisTOGrid_nD_cplx, GridTOBasis_nD_cplx, Calc_tab_Iq0, Calc_index, Debug_Enabled, Basis_Debug_Active
   USE NDindex_m, ONLY: Init_tab_ind, increase_NDindex
    USE SmolyakSparseGridMetadata_m, ONLY: ExpandSparseBasisToFullBasis, RestrictFullBasisToSparseBasis, &
      PrintSmolyakSparseGridDescriptor
   USE SmolyakSparseBasisAlgebra_m, ONLY: DotSmolyakBasis
   USE Ana_psi_m, ONLY: Calc_Av_imp_k_nD, Calc_AVQ_SQ_nD

   USE Molec_m
   USE mod_Tnum, ONLY: Tnum, CoordType
   USE mod_paramQ, ONLY: sub_QactTOQit
   USE mod_dnSVM, ONLY: Type_dnMat, Type_dnS, alloc_dnSVM, dealloc_dnSVM, Set_ZERO_TO_dnSVM
   USE mod_dnGG_dng, ONLY: get_dng_dnGG
   USE mod_dnDetGG_dnDetg, ONLY: sub3_dnDetGG
   Use Psi_m
   implicit none
   private
   TYPE :: Tdef_t
     ! Coefficients f2^{ij}(Q)
     real(kind=Rkind), allocatable :: F2(:,:,:) ! (nq, ndim, ndim)
  
     ! Coefficients f1^{i}(Q)
     real(kind=Rkind), allocatable :: F1(:,:) ! (nq, ndim)
  
     ! One scalar potential storage: V or V+Vep depending on the method.
     real(kind=Rkind), allocatable :: V(:,:,:) ! (nq, nsurf, nsurf)
  
     real(kind=Rkind), allocatable :: rho(:) ! (nq)
     ! Direct metric data used by the independent KEO kernels.
     real(kind=Rkind), allocatable :: G(:,:,:)       ! G(iq,i,j)
     real(kind=Rkind), allocatable :: J(:) ! volume Jacobian
     ! References retained for on-the-fly similarity metric evaluation.
     TYPE(Tnum), ALLOCATABLE :: para_Tnum_ref
     TYPE(CoordType), ALLOCATABLE :: mole_ref
     Character(len=Name_len)       :: model_type = '', coord_type = '', pot_type = '', sym_type = ''
     ! KEO construction selected in &coordQMLpot:
     !   developed  : expanded F2/F1/Vep form (PDF method 1)
     !   divergence : direct rho-weighted flux form (PDF method 2c)
     !   similarity : direct similarity-transformed flux form (PDF method 2b)
     Character(len=Name_len)       :: construction_type = 'developed'
     ! -1 keeps active coordinates; non-negative values select a Tnum stage.
     INTEGER(KIND=Ikind)           :: pot_itQtransfo = -1_Ikind
     real(kind=Rkind)              :: pot0

  END TYPE Tdef_t
  TYPE :: Op_t
     TYPE(Basis_t), pointer           :: Basis
  
     ! opérateurs matriciels existants
     real(kind=Rkind), allocatable    :: RMat(:,:)
     complex(kind=Rkind), allocatable :: CMat(:,:)
     logical                          :: UsesSmolyakSparseBasis = .FALSE.
  
     ! potentiel déjà existant
     real(kind=Rkind), allocatable    :: Scal_pot(:,:,:)
  
     ! Tdef matrix elements
     TYPE(Tdef_t)                     :: Tdef
   END TYPE Op_t

   PUBLIC :: Op_t, write_Op, Set_Op, dealloc_Op, calc_OpPsi, Calc_Hpsi, Kpsi_nD
   PUBLIC :: Calc_Scalar_Pot, test_openmp_op,calc_tab_Iq,calc_nac,calc_VV
   PUBLIC :: Build_CN_Matrices, Build_CN_Propagator,Save_Potential_ToFile
   PUBLIC :: InitQModelBasisParam, HPsiGridTdef, Tdef_F1_Grad_PsiGrid
   PUBLIC :: Tdef_F2_CrossDeriv_PsiGrid,TdefOnPsiGrid, Tdef_F2_DiagDerivPsiGrid
   PUBLIC :: Validate_Smolyak_Implementation, pdv2psi_nD
   PUBLIC :: H_action_diagnostic,Write_Potential

CONTAINS

   SUBROUTINE alloc_Op(Op, nb)
      TYPE(Op_t), intent(inout) :: Op
      integer, intent(in)    :: nb

      CALL dealloc_Op(Op)

      IF (nb < 1) STOP 'ERROR in init_Op: nb < 1!'

      allocate (Op%RMat(nb, nb))

   END SUBROUTINE alloc_Op

   SUBROUTINE dealloc_Op(Op)
      TYPE(Op_t), intent(inout) :: Op

      nullify (Op%Basis)

      IF (ALLOCATED(Op%RMat)) THEN
         DEALLOCATE (Op%RMat)
      END IF
      IF (ALLOCATED(Op%CMat)) THEN
         DEALLOCATE (Op%CMat)
      END IF
      IF (ALLOCATED(Op%Scal_pot)) THEN
         DEALLOCATE (Op%Scal_pot)
      END IF
      Op%UsesSmolyakSparseBasis = .FALSE.
      IF(ALLOCATED(Op%Tdef%F1)) THEN
         DEALLOCATE(Op%Tdef%F1)
      END IF
      IF(ALLOCATED(Op%Tdef%F2)) THEN
         DEALLOCATE(Op%Tdef%F2)
      END IF   
      IF(ALLOCATED(Op%Tdef%V)) DEALLOCATE(Op%Tdef%V)
      

       IF(ALLOCATED(Op%Tdef%rho)) THEN
         DEALLOCATE(Op%Tdef%rho)
      END IF
      IF(ALLOCATED(Op%Tdef%G)) DEALLOCATE(Op%Tdef%G)
      IF(ALLOCATED(Op%Tdef%J)) DEALLOCATE(Op%Tdef%J)
      IF(ALLOCATED(Op%Tdef%para_Tnum_ref)) DEALLOCATE(Op%Tdef%para_Tnum_ref)
      IF(ALLOCATED(Op%Tdef%mole_ref)) DEALLOCATE(Op%Tdef%mole_ref)
  
   END SUBROUTINE dealloc_Op

   SUBROUTINE write_Op(Op)
      USE QDUtil_m

      TYPE(Op_t), intent(in) :: Op

      ! integer :: i

      IF (associated(Op%Basis)) THEN
         write (out_unit, *) ' The basis is linked to Op.'
      END IF

      IF (allocated(Op%RMat)) THEN
         write (out_unit, *) 'Writing Op (real):'
         write (out_unit, *)
         CALL Write_VecMat(Op%RMat, out_unit, 5,  info='Op%Rmat')
         write (out_unit, *) 'END Writing Op'
      END IF

   END SUBROUTINE write_Op


!--------------------------------------------------------------------
!> @brief  Construct the Crank–Nicolson matrices M_minus and M_plus.
!>
!> Given the Hamiltonian matrix H, this routine builds the two 
!> matrices required for Crank–Nicolson propagation:
!>      H_minus = I - (i dt / 2) H
!>      H_plus  = I + (i dt / 2) H
!>
!> @param[in]  H       Hamiltonian matrix (nb x nb)
!> @param[in]  dt      Time step Δt
!> @param[out] H_minus Matrix (nb x nb): I - (i dt / 2) H
!> @param[out] H_plus  Matrix (nb x nb): I + (i dt / 2) H
!--------------------------------------------------------------------
SUBROUTINE Build_CN_Matrices(H, dt, H_minus, H_plus)
   IMPLICIT NONE
 
   !------------ Arguments --------------------------------------------
   COMPLEX(KIND=Rkind), INTENT(IN)  :: H(:,:)
   REAL(KIND=Rkind),    INTENT(IN)  :: dt
   COMPLEX(KIND=Rkind), ALLOCATABLE, INTENT(OUT) :: H_minus(:,:), H_plus(:,:)
 
   !------------ Local variables ----------------------------------------

   INTEGER                       :: nb, i
   COMPLEX(KIND=Rkind)           :: coef
 
   !------------ Get matrix size ----------------------------------------
   nb = SIZE(H, 1)
 
   !------------ Allocate outputs ---------------------------------------
   ALLOCATE(H_minus(nb, nb), H_plus(nb, nb))
 
   !------------ Compute prefactor coef = i dt / 2 ----------------------
   coef = EYE * dt / TWO
 
   !------------ Build matrices -----------------------------------------
   H_minus(:,:) = -coef * H(:,:)
   H_plus(:,:)  =  coef * H(:,:)
 
   ! Add identity contribution on the diagonal
   DO i = 1, nb
      H_minus(i,i) = H_minus(i,i) + CONE
      H_plus(i,i)  = H_plus(i,i)  + CONE
   END DO
 
   !------------ Debug info ---------------------------------------------
   IF (Debug_Enabled()) THEN
      WRITE(out_unit,*) '>>> Build_CN_Matrices: matrices constructed'
      WRITE(out_unit,*) '    nb   =', nb
      WRITE(out_unit,*) '    dt   =', dt
      FLUSH(out_unit)
   END IF
 
 END SUBROUTINE Build_CN_Matrices


 !---------------------------------------------------------------------
 !> @brief  Construct the Crank–Nicolson propagator U.
 !>
 !> This routine builds the Crank–Nicolson (Cayley) propagator
 !>      U = (I + i dt/2 H)^(-1) (I - i dt/2 H)
 !>
 !> It calls Build_CN_Matrices to generate M_minus and M_plus,
 !> then inverts M_plus with inv_Of_Mat_TO.
 !>
 !> @param[in]  H   Hamiltonian matrix (nb x nb)
 !> @param[in]  dt  Time step Δt
 !> @param[out] U   Propagator matrix (nb x nb)
 !---------------------------------------------------------------------
 SUBROUTINE Build_CN_Propagator(H, dt, U)
    IMPLICIT NONE
 
    !-------- Arguments ------------------------------------------------
    COMPLEX(KIND=Rkind), INTENT(IN)  :: H(:,:)
    REAL(KIND=Rkind),    INTENT(IN)  :: dt
    COMPLEX(KIND=Rkind), ALLOCATABLE, INTENT(OUT) :: U(:,:)
 
    !-------- Local variables -------------------------------------------

    INTEGER            :: nb
    COMPLEX(KIND=Rkind), ALLOCATABLE :: M_minus(:,:), M_plus(:,:)
 
    !-------- Get matrix size -------------------------------------------
    nb = SIZE(H, 1)
 
    !-------- Build M_minus and M_plus ----------------------------------
    CALL Build_CN_Matrices(H, dt, M_minus, M_plus)
    IF (Debug_Enabled()) CALL Write_VecMat(H, out_unit, 5, info='H')
 
    !-------- Invert the left-hand-side matrix I + i*dt*H/2 ------------
    M_plus = inv_OF_Mat(M_plus,inv_type=0,epsi=1.d-20)
 
    !-------- Compute propagator ----------------------------------------
    ALLOCATE(U(nb, nb))
    U(:,:) = MATMUL(M_plus, M_minus)
 
    !-------- Debug info ------------------------------------------------
    IF (Debug_Enabled()) THEN
       WRITE(out_unit,*) '>>> Build_Propagator: propagator constructed'
       WRITE(out_unit,*) '    nb   =', nb
       WRITE(out_unit,*) '    dt   =', dt
       FLUSH(out_unit)
    END IF


    !Print the propagator
    IF (Debug_Enabled()) CALL Write_VecMat(U, out_unit, 5, info='U')
    !Cleanup
    IF (ALLOCATED(M_minus)) DEALLOCATE(M_minus)
    IF (ALLOCATED(M_plus))  DEALLOCATE(M_plus)
 
 END SUBROUTINE Build_CN_Propagator

!===============================================================================
!> SUBROUTINE: write_3D_array_column_major
!>
!> Purpose:
!>   Writes a 3D real array V(:,:,:) in column-major order.
!>   The output is linearized such that the first index varies fastest,
!>   then the second, then the third (Fortran-style storage order).
!>
!> Arguments:
!>   - V         : (IN)  3D real array
!>   - out_unit  : (IN)  Fortran output unit (screen or file)
!>
!> Notes:
!>   - Dimensions are automatically determined from V.
!>   - Each line contains one element with its linearized index and (i,j,k) indices.
!===============================================================================
 SUBROUTINE write_3D_array(V, out_unit)
   IMPLICIT NONE
   !----------------- Arguments -----------------
   REAL(KIND=Rkind), INTENT(IN) :: V(:,:,:)
   INTEGER, INTENT(IN) :: out_unit
 
   !----------------- Local variables -----------------
   INTEGER :: i, j, k, idx
   INTEGER :: dimx, dimy, dimz
 
   !----------------- Determine array size -----------------
   dimx = SIZE(V,1)
   dimy = SIZE(V,2)
   dimz = SIZE(V,3)
 
   WRITE(out_unit,*) '--- Writing 3D array in column-major order ---'
   WRITE(out_unit,*) 'Dimensions: (', dimx, ',', dimy, ',', dimz, ')'
 
   idx = 0
   DO k = 1, dimz
      DO j = 1, dimy
         DO i = 1, dimx
            idx = idx + 1
            WRITE(out_unit,'(A,I6,A,I4,A,I4,A,I4,A,F14.6)') 'Index ', idx, ': V(', i, ',', j, ',', k, ') = ', V(i,j,k)
         END DO
      END DO
   END DO
 
   !----------------- Flush output -----------------
   FLUSH(out_unit)
 END SUBROUTINE write_3D_array
 
!----------------------------------------------------------------------
!> @brief  Initialize or update the operator Op based on the given Basis.
!>
!> This routine sets up the operator `Op` using the information contained
!> in the basis `Basis` and the index table `Tab_Iq`. The behavior can be 
!> controlled with an optional logical argument `alloc`:
!>   - If `alloc=.TRUE.`, the operator structure is fully allocated 
!>     (memory allocation + initialization).
!>   - If `alloc=.FALSE.` (default), only the scalar potential `Scal_pot`
!>     is computed using the provided basis.
!>
!> A local debug flag is available to print important internal information
!> such as dimension sizes, number of surfaces, and allocation choices.
!>
!> @param[inout] Op      Operator object to initialize
!> @param[in]    Basis   Basis object defining the system
!> @param[in]    Tab_Iq  Index table for quadrature points
!> @param[in]    alloc   [optional] logical flag controlling allocation
!----------------------------------------------------------------------
SUBROUTINE Set_Op(Op, Basis, Tab_Iq, alloc, nio)
   IMPLICIT NONE
 
   !--------- Arguments --------------------------------------------------
   TYPE(Op_t),     INTENT(INOUT)         :: Op
   TYPE(Basis_t),  INTENT(IN),   TARGET  :: Basis
   INTEGER,        INTENT(IN)            :: Tab_Iq(:, :)
   LOGICAL,        INTENT(IN), OPTIONAL  :: alloc
   INTEGER,        INTENT(INOUT)             :: nio
 
   !--------- Local variables ------------------------------------------
   INTEGER :: ndim, nsurf, nq,nb,err_io,rewind_stat
   LOGICAL :: do_alloc

   CHARACTER (len=Name_len)  :: coord_type, vcoord_type, construction_type
   CHARACTER (len=Name_len)  :: sym_type, model_type
   INTEGER                   :: effective_basis_size
   INTEGER(KIND=Ikind)       :: pot_itQtransfo
   REAL(kind=Rkind)          :: pot0

   ! vcoord_type is accepted as an alias for coord_type so old and new
   ! input decks can coexist. construction_type selects the KEO form.
   NAMELIST /coordQMLpot/ coord_type,vcoord_type,construction_type,sym_type,model_type,pot0,pot_itQtransfo
 
   !======== Check basis ===============================================
   IF (.NOT. Basis_IS_allocated(Basis)) THEN
      STOP 'ERROR in Set_Op: the Basis is not initialized'
   END IF
 
   !======== Extract dimensions ========================================
   ndim  = SIZE(Basis%tab_basis) - 1
   nq    = Basis%nq
   nsurf = Basis%tab_basis(ndim+1)%nb
   nb = Basis%nb*nsurf
   effective_basis_size = Get_TotalBasisCoefficientCount(Basis)
   Op%Basis => Basis

   coord_type        = 'CC'
   vcoord_type      = ''
    construction_type = 'developed'
   sym_type          = 'SYM'
   model_type        = 'QML'
   pot0              = ZERO
      pot_itQtransfo    = -1_Ikind
 
   !======== Decide allocation mode ====================================
   IF (PRESENT(alloc)) THEN
      do_alloc = alloc
   ELSE
      do_alloc = .FALSE.
   END IF
   ! stdin is a pipe in normal runs and cannot be REWIND'ed.  Read the
   ! coordinate namelist only on the first call and reuse the stored choice
   ! when Hagedorn rebuilds the operator on a new basis.
   IF (LEN_TRIM(Op%Tdef%coord_type) == 0) THEN
      REWIND(nio,IOSTAT=rewind_stat)
      IF (rewind_stat /= 0) THEN
         ERROR STOP 'Set_Op: input unit must be seekable on first curvilinear setup'
      END IF
      READ(nio,NML=coordQMLpot,IOSTAT=err_io)
      IF (err_io /= 0) THEN
         WRITE(*,*) 'Error reading coordQMLpot namelist, IOSTAT =', err_io
         ERROR STOP 'Set_Op: unable to read coordQMLpot'
      END IF
      IF (LEN_TRIM(vcoord_type) > 0) coord_type = vcoord_type
      CALL string_uppercase_TO_lowercase(construction_type)
      Op%Tdef%coord_type = coord_type
      Op%Tdef%construction_type = construction_type
      Op%Tdef%model_type = model_type
      Op%Tdef%sym_type = sym_type
      Op%Tdef%pot0 = pot0
      Op%Tdef%pot_itQtransfo = pot_itQtransfo
   ELSE
      coord_type = Op%Tdef%coord_type
      construction_type = Op%Tdef%construction_type
      model_type = Op%Tdef%model_type
      sym_type = Op%Tdef%sym_type
      pot0 = Op%Tdef%pot0
      pot_itQtransfo = Op%Tdef%pot_itQtransfo
   END IF

    WRITE(out_unit,*)''
    WRITE(out_unit,*) '============== coordQMLpot ========================='
    WRITE(out_unit,*) 'Op%Tdef%coord_type:=',Op%Tdef%coord_type
    WRITE(out_unit,*) 'Op%Tdef%construction_type:=',Op%Tdef%construction_type
    WRITE(out_unit,*) 'Op%Tdef%model_type:=',Op%Tdef%model_type
    WRITE(out_unit,*) 'Op%Tdef%pot_itQtransfo:=',Op%Tdef%pot_itQtransfo
    WRITE(out_unit,*) 'Op%Tdef%pot0:=',Op%Tdef%pot0
    WRITE(out_unit,*)''
   !======== Perform allocation or scalar potential calculation ========
   SELECT CASE (Op%Tdef%coord_type)
   CASE('CC') ! CC=Cartesian Coor
   IF (do_alloc) THEN
      ! Full allocation of the operator structure
      IF(.NOT. ALLOCATED(Op%CMat)) ALLOCATE(Op%CMat(effective_basis_size,effective_basis_size))
      IF(.NOT. ALLOCATED(Op%Scal_pot)) ALLOCATE(Op%Scal_pot(nq, nsurf, nsurf))

     ! Initialize the operator matrix 
      Op%CMat(:,:) = CZERO
      Op%Scal_pot(:,:,:) = ZERO

     CALL Calc_Scalar_Pot(Op%Scal_pot, Basis)


     IF (Debug_Enabled()) CALL Write_VecMat(Op%CMat, out_unit, 5, info='Op%CMat')
     IF (Debug_Enabled()) CALL write_3D_array(Op%Scal_pot, out_unit)
 
      IF (Debug_Enabled()) THEN
         WRITE(out_unit,*) '>>> Set_Op: Operator fully allocated'
         WRITE(out_unit,*) '    ndim   =', ndim
         WRITE(out_unit,*) '    nq     =', nq
         WRITE(out_unit,*) '    nb     =', Basis%nb
         WRITE(out_unit,*) '    nsurf  =', nsurf
      END IF
   ELSE
      ! Only compute scalar potential
      CALL Calc_Scalar_Pot(Op%Scal_pot, Basis)
 
      IF (Debug_Enabled()) THEN
         WRITE(out_unit,*) '>>> Set_Op: Only scalar potential computed'
         WRITE(out_unit,*) '    ndim   =', ndim
         WRITE(out_unit,*) '    nq     =', nq
         WRITE(out_unit,*) '    nb     =', Basis%nb
         WRITE(out_unit,*) '    nsurf  =', nsurf
      END IF
   END IF
   CASE('CURV')
   CALL Set_Op_Tdef(Op, Tab_Iq, alloc)
   CASE DEFAULT
    STOP 'Unknown coordinate type'
   END SELECT
 
 END SUBROUTINE Set_Op



 SUBROUTINE Set_Op_Tdef(Op, Tab_Iq, alloc)
  USE mod_dnSVM
  USE mod_Constant
  USE TnumTana_system_m
  USE mod_Coord_KEO                
  USE mod_PrimOp


   IMPLICIT NONE
 
   !--------- Arguments --------------------------------------------------
   TYPE(Op_t),     INTENT(INOUT)         :: Op
   INTEGER,        INTENT(IN)            :: Tab_Iq(:, :)
   LOGICAL,        INTENT(IN), OPTIONAL  :: alloc
 
   !--------- Local variables ------------------------------------------
   INTEGER :: ndim, nsurf, nq,nb
   INTEGER(KIND=8) :: dense_bytes
   LOGICAL :: do_alloc

     ! Paramètres
   TYPE(constant)  :: const_phys ! Constantes physiques
   TYPE(CoordType) :: mole      ! Description de la molécule
   TYPE(Tnum)      :: para_Tnum ! Paramètres pour le KEO
   TYPE(PrimOp_t)  :: PrimOp

 
   !======== Check basis ===============================================
   IF (.NOT. Basis_IS_allocated(Op%Basis)) THEN
      STOP 'ERROR in Set_Op: the Basis is not initialized'
   END IF
 
   !======== Extract dimensions ========================================
   ndim  = SIZE(Op%Basis%tab_basis) - 1
   nq    = Op%Basis%nq
   nsurf = Op%Basis%tab_basis(ndim+1)%nb
   nb = Op%Basis%nb*nsurf
  ! Tnum coordinate setup is immutable during a propagation.  Reuse the
  ! finalized descriptors when Hagedorn rebuilds only the primitive basis.
  ! This avoids rereading the coordinate namelists and rebuilding all Tnum
  ! transformation metadata at every time step.
  IF (ALLOCATED(Op%Tdef%para_Tnum_ref) .AND. ALLOCATED(Op%Tdef%mole_ref)) THEN
     para_Tnum = Op%Tdef%para_Tnum_ref
     mole      = Op%Tdef%mole_ref
  ELSE
     CALL Read_CoordType(mole, para_Tnum, const_phys)
     CALL read_RefGeom(mole, para_Tnum)
     CALL Finalize_TnumTana_Coord_PrimOp(para_Tnum, mole, PrimOp)
     ALLOCATE(Op%Tdef%para_Tnum_ref, Op%Tdef%mole_ref)
     Op%Tdef%para_Tnum_ref = para_Tnum
     Op%Tdef%mole_ref = mole
  END IF

   !======== Decide allocation mode ====================================
   IF (PRESENT(alloc)) THEN
      do_alloc = alloc
   ELSE
      do_alloc = .FALSE.
   END IF
 
   !======== Perform allocation or scalar potential calculation ========
   IF (do_alloc) THEN
      dense_bytes = 16_8 * INT(nb, KIND=8) * INT(nb, KIND=8)
      IF (dense_bytes > 2_8 * 1024_8**3) THEN
         WRITE(out_unit,*) 'Dense Hamiltonian request refused: basis size = ', nb, &
                           ', matrix storage = ', dense_bytes / 1024_8**3, ' GiB.'
         ERROR STOP 'Select matrix-free propagation or reduce the basis dimension.'
      END IF
      ! Full allocation of the operator structure
     IF(.NOT. ALLOCATED(Op%CMat)) ALLOCATE(Op%CMat(nb,nb))
     IF(.NOT. ALLOCATED(Op%Scal_pot)) ALLOCATE(Op%Scal_pot(nq, nsurf, nsurf))

     ! Initialize the operator matrix 
      Op%CMat(:,:) = CZERO
      Op%Scal_pot(:,:,:) = ZERO
     CALL Init_Op_Tdef_F2_F1_Vep_FromTnum(Op, para_Tnum, mole)

     IF (Debug_Enabled()) CALL Write_VecMat(Op%CMat, out_unit, 5, info='Op%CMat')
     IF (Debug_Enabled()) THEN
        IF (ALLOCATED(Op%Scal_pot)) CALL write_3D_array(Op%Scal_pot, out_unit)
     END IF
 
      IF (Debug_Enabled()) THEN
         WRITE(out_unit,*) '>>> Set_Op: Operator fully allocated'
         WRITE(out_unit,*) '    ndim   =', ndim
         WRITE(out_unit,*) '    nq     =', nq
         WRITE(out_unit,*) '    nb     =', nb
         WRITE(out_unit,*) '    nsurf  =', nsurf
      END IF
   ELSE
      ! Only compute scalar potential
      CALL Init_Op_Tdef_F2_F1_Vep_FromTnum(Op, para_Tnum, mole)
 
      IF (Debug_Enabled()) THEN
         WRITE(out_unit,*) '>>> Set_Op: Only scalar potential computed'
         WRITE(out_unit,*) '    ndim   =', ndim
         WRITE(out_unit,*) '    nq     =', nq
         WRITE(out_unit,*) '    nb     =', nb
         WRITE(out_unit,*) '    nsurf  =', nsurf
      END IF
   END IF
  ! Clean Up
  CALL dealloc_CoordType(mole)

 END SUBROUTINE Set_Op_Tdef

 
 
 SUBROUTINE Kpsi_nD(KPsi_g, Psi_g, Basis)
 
   USE Basis_m
   USE  QDUtil_m
   USE Molec_m
 
      TYPE(Basis_t), INTENT(IN), TARGET               :: Basis
      COMPLEX(KIND=Rkind), INTENT(IN), TARGET         :: psi_g(:)
      COMPLEX(KIND=Rkind), INTENT(INOUT), TARGET      :: Kpsi_g(:)
      COMPLEX(KIND=Rkind), POINTER                    :: psi_ggb(:, :, :)
      COMPLEX(KIND=Rkind), POINTER                    :: d2gg(:, :)
      COMPLEX(KIND=Rkind), POINTER                    :: Kpsi_ggb(:, :, :)
      REAL(KIND=Rkind), ALLOCATABLE                   :: GGdef(:,:)
      LOGICAL                                         :: debug
      INTEGER                                         :: iq, i1, i3, inb, ndim
      INTEGER                                         :: Iq1, Iq2, Iq3
 
      debug = Basis%debug .OR. Debug_Enabled()
      IF (debug) THEN
         !write(out_unit,*) 'BEGINNING Kpsi'
         flush (out_unit)
      END IF
      ndim = size(Basis%tab_basis)-1
      ALLOCATE(GGdef(ndim , ndim))
 
      CALL get_Qmodel_GGdef(GGdef)

      Kpsi_g(:) = CZERO
 
      DO inb = 1, ndim 

         IF (inb == 1) THEN

         Iq1 = 1
         Iq2 = Basis%tab_basis(1)%nq
         Iq3 = Product(Basis%tab_basis(2:ndim)%nq)*Basis%tab_basis(ndim + 1)%nb
         ELSE IF (inb == ndim) THEN

         Iq1  = Product(Basis%tab_basis(1:ndim - 1)%nq)
         Iq2  = Basis%tab_basis(ndim)%nq
         Iq3  = Basis%tab_basis(ndim + 1)%nb

      ELSE
         Iq1 = Product(Basis%tab_basis(1:inb - 1)%nq)
         Iq2 = Basis%tab_basis(inb)%nq
         Iq3 = Product(Basis%tab_basis(inb + 1:ndim)%nq)*Basis%tab_basis(ndim + 1)%nb
      END IF
 
      Kpsi_ggb(1:Iq1, 1:Iq2, 1:Iq3) => Kpsi_g
      psi_ggb(1:Iq1, 1:Iq2, 1:Iq3) => psi_g
      d2gg(1:Iq2, 1:Iq2) => Basis%tab_basis(inb)%d2gg(:,:,1,1)

         DO i3 = 1, ubound(psi_ggb, dim=3)
 
            DO i1 = 1, ubound(psi_ggb, dim=1)

               KPsi_ggb(i1, :, i3) = KPsi_ggb(i1, :, i3) - HALF*GGdef(inb, inb)*matmul(d2gg, psi_ggb(i1, :, i3))
 
            END DO
 
         END DO
 
      END DO
 
      IF (debug) THEN
         write(out_unit,*) 'END KPsi_nD'
         flush (out_unit)
      END IF
 
   END SUBROUTINE Kpsi_nD

   SUBROUTINE Calc_Scalar_Pot(V, Basis, mole, pot_itQtransfo)
      use QDUtil_m
      use Basis_m
      use omp_lib
      use mod_Tnum, only: CoordType
      use mod_paramQ, only: sub_QactTOQit
      implicit none
   
      ! ===== Arguments =====
      type(Basis_t), intent(in), target                :: Basis
      real(kind=Rkind), allocatable, intent(inout)     :: V(:, :, :)
      type(CoordType), intent(in), optional            :: mole
      integer(kind=Ikind), intent(in), optional        :: pot_itQtransfo
   
      ! ===== Local variables =====
      integer(kind=Ikind)                              :: ndim, iq, nq, nsurf
      real(kind=Rkind), allocatable                    :: Q(:,:)
      real(kind=Rkind), allocatable                    :: Qpot(:,:), Qit(:)
      integer(kind=Ikind)                              :: it_pot
   
      ! ===== Initialization =====
      ndim  = size(Basis%tab_basis, KIND=Ikind) - 1_Ikind ! Number of nuclear degrees of freedom
      nq    = Basis%nq                         ! Total number of grid points
      nsurf = Basis%tab_basis(ndim + 1)%nb     ! Number of electronic surfaces
      it_pot = -1_Ikind
      if (present(pot_itQtransfo)) it_pot = pot_itQtransfo
      IF (it_pot < -1_Ikind) ERROR STOP 'Calc_Scalar_Pot: pot_itQtransfo must be -1 or a valid Tnum stage.'
   
      IF (.NOT. allocated(V)) allocate(V(nq, nsurf, nsurf))
      CALL Calc_Q_grid(Q, Basis)

      ! Match ElVibRot's pot_itQtransfo semantics: send each active grid point
      ! through Tnum's own coordinate chain before evaluating the QML PES.
      IF (it_pot >= 0_Ikind) THEN
         IF (.NOT. PRESENT(mole)) ERROR STOP 'Calc_Scalar_Pot: pot_itQtransfo requires Tnum coordinates.'
         IF (it_pot > mole%nb_Qtransfo) ERROR STOP 'Calc_Scalar_Pot: pot_itQtransfo exceeds nb_Qtransfo.'
         IF (SIZE(Q, DIM=2, KIND=Ikind) /= mole%nb_act) &
            ERROR STOP 'Calc_Scalar_Pot: active grid dimension differs from Tnum nb_act.'
         DO iq = 1, nq
            CALL sub_QactTOQit(Q(iq,:), Qit, it_pot, mole)
            IF (iq == 1_Ikind) THEN
               ALLOCATE(Qpot(nq,SIZE(Qit,KIND=Ikind)))
               WRITE(out_unit,*) 'QML potential uses Tnum transformation stage ',it_pot
               WRITE(out_unit,*) 'First grid point dimensions: active=',SIZE(Q,DIM=2,KIND=Ikind), &
                                 ' PES=',SIZE(Qit,KIND=Ikind)
               WRITE(out_unit,*) 'First grid point Q(active) -> Q(PES): ',Q(iq,:),Qit(:)
            ELSE IF (SIZE(Qit,KIND=Ikind) /= SIZE(Qpot,DIM=2,KIND=Ikind)) THEN
               ERROR STOP 'Calc_Scalar_Pot: transformed-coordinate dimension changed on the grid.'
            END IF
            Qpot(iq,:) = Qit(:)
            DEALLOCATE(Qit)
         END DO
      END IF
   
      ! ===== Parallel loop avec variables privées =====
      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP          SHARED(V, Q, Qpot, it_pot, nq, nsurf, ndim) &
      !$OMP          PRIVATE(iq)
      
      !$OMP DO SCHEDULE(DYNAMIC)
      Do iq = 1, nq
         IF (it_pot >= 0) THEN
            call calc_potential_point(V, Qpot, iq, nsurf, SIZE(Qpot,DIM=2,KIND=Ikind))
         ELSE
            call calc_potential_point(V, Q, iq, nsurf, ndim)
         END IF
      end do
      !$OMP END DO
   
      !$OMP END PARALLEL

      !call Save_Potential_ToFile(V,Q)
   
      ! ===== Deallocation =====
      deallocate(Q)
      IF (ALLOCATED(Qpot)) DEALLOCATE(Qpot)
   
   END SUBROUTINE Calc_Scalar_Pot
   
   ! ===== Subroutine auxiliaire pour éviter les tableaux temporaires =====
   SUBROUTINE calc_potential_point(V, Q, iq, nsurf, ndim)
      use QDUtil_m
      implicit none
      
      ! Arguments
      real(kind=Rkind), intent(inout)                  :: V(:, :, :)
      real(kind=Rkind), intent(in)                     :: Q(:, :)
      integer(kind=Ikind), intent(in)                  :: iq, nsurf, ndim
      
      ! Variables locales contiguës
      real(kind=Rkind)                                 :: V_local(nsurf, nsurf)
      real(kind=Rkind)                                 :: Q_local(ndim)
      integer(kind=Ikind)                              :: i, j
      
      ! Copie explicite pour garantir la contiguïté
      do i = 1_Ikind, ndim
         Q_local(i) = Q(iq, i)
      end do
      
      ! Initialisation pour respecter -finit-real=nan
      V_local = ZERO
      
      ! Appel avec tableaux continus
      call sub_Qmodel_V(V_local, Q_local)
      
      ! Diagnostic ciblé : QML ne doit jamais produire NaN/Inf.
      IF (ANY(V_local /= V_local)) THEN
         WRITE(*,*) 'QML non-finite potential at iq=',iq,' Q=',Q_local
         ERROR STOP 'non-finite QML potential'
      END IF

      ! Copie du résultat avec boucles explicites
      do j = 1_Ikind, nsurf
         do i = 1_Ikind, nsurf
            V(iq, i, j) = V_local(i, j)
         end do
      end do
      
   END SUBROUTINE calc_potential_point


   SUBROUTINE Write_Potential(Basis)

   USE QDUtil_m

   IMPLICIT NONE

   ! Argument
   TYPE(Basis_t), INTENT(IN) :: Basis

   ! Tableaux
   REAL(KIND=Rkind), ALLOCATABLE :: V(:, :)
   REAL(KIND=Rkind), ALLOCATABLE :: Q(:, :)
   REAL(KIND=Rkind), ALLOCATABLE :: Qpoint(:)

   ! Indices et dimensions
   INTEGER :: iq
   INTEGER :: i, j
   INTEGER :: nq
   INTEGER :: ndim
   INTEGER :: nsurf

   ! Fichier
   INTEGER :: unit_file
   INTEGER :: ios
   CHARACTER(LEN=200) :: iomsg

   ! ===== Initialisation =====

   ndim  = SIZE(Basis%tab_basis) - 1
   nq    = Basis%nq
   nsurf = Basis%tab_basis(ndim + 1)%nb

   ALLOCATE(V(nsurf, nsurf))
   ALLOCATE(Qpoint(ndim))

   ! Calcul de la grille Q
   CALL Calc_Q_grid(Q, Basis)

   ! Ouverture du fichier avec un numéro automatique
   OPEN(NEWUNIT=unit_file,FILE='potential.dat',STATUS='REPLACE',ACTION='WRITE',FORM='FORMATTED',IOSTAT=ios,IOMSG=iomsg)

   IF (ios /= 0) THEN
      WRITE(*,*) 'Erreur lors de l''ouverture du fichier : ', TRIM(iomsg)

      DEALLOCATE(V)
      DEALLOCATE(Qpoint)

      IF (ALLOCATED(Q)) DEALLOCATE(Q)

      RETURN
   END IF

   ! ===== Calcul et écriture =====

   DO iq = 1, nq

      ! Copie dans un tableau contigu
      Qpoint = Q(iq, :)

      ! Initialisation du potentiel
      V = ZERO

      ! Calcul du potentiel au point Qpoint
      CALL sub_Qmodel_V(V, Qpoint)

      ! Écriture :
      ! Q(1), Q(2), ..., Q(ndim),
      ! V(1,1), V(1,2), ..., V(1,nsurf),
      ! V(2,1), ..., V(nsurf,nsurf)
      WRITE(unit_file, *) Qpoint, &
                          ((V(i, j), j = 1, nsurf), i = 1, nsurf)

   END DO

   ! Fermeture du fichier
   CLOSE(unit_file)

   ! Libération de la mémoire
   DEALLOCATE(V)
   DEALLOCATE(Qpoint)

   IF (ALLOCATED(Q)) DEALLOCATE(Q)

END SUBROUTINE Write_Potential

SUBROUTINE Save_Potential_ToFile(V,Basis)
  USE QDUtil_m
  IMPLICIT NONE

  ! Arguments
  real(kind=Rkind), INTENT(IN)                          :: V(:,:,:)
  TYPE(Basis_t), INTENT(IN), TARGET                     :: Basis

  ! Local variables
  INTEGER :: ndim, nsurf, nq, Ib, Iq,I,II
  real(kind=Rkind), ALLOCATABLE :: Q(:,:)

  CALL Calc_Q_grid(Q, Basis)

  !-----------------------------------
  ! Dimensions
  !-----------------------------------
  
  ndim  = SIZE(Basis%tab_basis)-1
  nq    = Basis%nq
  nsurf = Basis%tab_basis(ndim+1)%nb



  !-----------------------------------
  ! Save to text file
  !-----------------------------------
  OPEN(unit=10, file='potential.txt')

  DO Iq = 1, nq
   ! Comparaison pour ligne vide
   IF(ndim>1)THEN
     IF (Iq > 1) THEN
         IF (Q(Iq-1,ndim) /= Q(Iq,ndim)) WRITE(10,*)
     END IF
   END IF
   WRITE(10,*) (Q(Iq,I), I=1,ndim), V(Iq,1,1)
  END DO
  CLOSE(10)

END SUBROUTINE Save_Potential_ToFile
  !
  ! OBJECTIVE: Build, diagonalize H matrix, display the first N energies, 
  !            and write the first N eigenstates to files.
  ! Energies are assumed to be in Atomic Units (Hartree).
  !


 SUBROUTINE Calc_Hpsi(psi_g, HPsi_g, Basis,V)
    USE Basis_m
    USE Psi_m
    USE Molec_m
    complex(kind=Rkind), intent(in), target            :: psi_g(:)
    complex(kind=Rkind), intent(inout)                 :: HPsi_g(:)
    type(Basis_t), intent(in), target                  :: Basis
    real(kind=Rkind) , intent(in)                      :: V(:, :, :)

    complex(kind=Rkind), allocatable, target           :: VPsi_g(:), KPsi_g(:)
    complex(kind=Rkind), pointer                       :: VPsi_gb(:, :), Psi_gb(:, :)
    Integer                                            :: i, j, ndim,nsurf,nq


    IF (.not. allocated(Basis%tab_basis)) THEN
       STOP 'ERROR in Set_Op: the Basis%tab_bais is not initialized'
    END IF

    ndim = SIZE(Basis%tab_basis) - 1
    nq =Basis%nq
    nsurf=Basis%tab_basis(ndim+1)%nb

    ! action potential V|Psi_g> 
    allocate(VPsi_g(nq*nsurf))
    allocate(KPsi_g(nq*nsurf))

    VPsi_g(:) = CZERO
    KPsi_g(:) = CZERO
    HPsi_g(:) = CZERO

   VPsi_gb(1:nq,1:nsurf) => VPsi_g
   Psi_gb(1:nq, 1:nsurf) => Psi_g

   DO i = 1, nsurf

      DO j = 1, nsurf

         VPsi_gb(:, i) = VPsi_gb(:, i) + V(:, i, j)*Psi_gb(:, j)

      END DO

   END DO

    ! action potential K|Psi_g>

     call Kpsi_nD(KPsi_g, Psi_g, Basis)
     HPsi_g(:) = VPsi_g(:) + KPsi_g(:)

    DEALLOCATE (VPsi_g)
    DEALLOCATE (KPsi_g)

 END SUBROUTINE 



  SUBROUTINE HPsiGridTdef(psi_g, HPsi_g, H)
    USE Basis_m
    USE Psi_m
    USE Molec_m
    COMPLEX(KIND=Rkind), INTENT(IN), TARGET            :: psi_g(:)
    COMPLEX(KIND=Rkind), INTENT(INOUT)                 :: HPsi_g(:)
    TYPE(Op_t), INTENT(IN), TARGET                     :: H

    COMPLEX(KIND=Rkind), ALLOCATABLE, TARGET           :: VeffPsi_g(:), TdefPsi_g(:)
    COMPLEX(KIND=Rkind), POINTER                       :: VeffPsi_gb(:, :), Psi_gb(:, :)
    INTEGER                                            :: i, j, ndim,nsurf,nq


    IF (.NOT. ALLOCATED(H%Basis%tab_basis)) THEN
       STOP 'ERROR in Set_Op: the Basis%tab_bais is not initialized'
    END IF

    ndim = SIZE(H%Basis%tab_basis) - 1
    nq =H%Basis%nq
    nsurf=H%Basis%tab_basis(ndim+1)%nb

    ! action potential V|Psi_g> 
    ALLOCATE(VeffPsi_g(nq*nsurf))
    ALLOCATE(TdefPsi_g(nq*nsurf))

    VeffPsi_g(:) = CZERO
    TdefPsi_g(:) = CZERO
    HPsi_g(:) = CZERO

   VeffPsi_gb(1:nq,1:nsurf) => VeffPsi_g
   Psi_gb(1:nq, 1:nsurf) => Psi_g

   DO i = 1, nsurf

      DO j = 1, nsurf

         IF (TRIM(H%Tdef%construction_type) == 'similarity') THEN
            VeffPsi_gb(:, i) = VeffPsi_gb(:, i) + H%Tdef%V(:, i, j)*Psi_gb(:, j)
         ELSE
            VeffPsi_gb(:, i) = VeffPsi_gb(:, i) + H%Tdef%V(:, i, j)*Psi_gb(:, j)
         END IF

      END DO

   END DO

    ! action potential K|Psi_g>
     CALL TdefOnPsiGrid(TdefPsi_g, Psi_g, H)
     HPsi_g(:) = TdefPsi_g(:) +VeffPsi_g(:)

    !Clean Up
    DEALLOCATE (VeffPsi_g)
    DEALLOCATE (TdefPsi_g)

 END SUBROUTINE 




 SUBROUTINE calc_OpPsi_temp(Op, psi, Oppsi)

   USE psi_m
   TYPE(Op_t), intent(in)           :: Op
   TYPE(psi_t), intent(in)          :: psi
   TYPE(psi_t), intent(inout)       :: Oppsi
   TYPE(psi_t)                      :: psi_g, Oppsi_g

   IF (allocated(Op%RMat)) THEN

      Oppsi%CVec = matmul(Op%RMat, psi%CVec)

   ELSE

      IF (psi%Grid) THEN

         call init_psi(Oppsi_g, Psi%Basis, .true., .true.)
         call Calc_Hpsi(psi%CVec, Oppsi_g%CVec, psi%Basis,Op%Scal_pot)
         call GridTOBasis_nD_cplx(Oppsi%CVec, Oppsi_g%CVec, psi%Basis)
         call dealloc_psi(Oppsi_g)

      ELSE

         call init_psi(psi_g, Psi%Basis, .true., .true.)
         call init_psi(Oppsi_g, psi%Basis, .true., .true.)
         call BasisTOGrid_nD_cplx(psi_g%CVec, psi%CVec, psi%Basis)
         call Calc_Hpsi(psi_g%CVec, Oppsi_g%CVec, psi%Basis,Op%Scal_pot)
         call GridTOBasis_nD_cplx(Oppsi%CVec, Oppsi_g%CVec, psi%Basis)
         call dealloc_psi(psi_g)
         call dealloc_psi(Oppsi_g)

      END IF

   END IF

END SUBROUTINE 


SUBROUTINE calc_OpPsi(Op, psi, Oppsi)
   !---------------------------------------------------------------
   !> Applies the operator Op to the wavefunction psi to compute
   !> Oppsi = Op * psi. This works for both grid and basis representations.
   !---------------------------------------------------------------
 
   USE psi_m
 
   TYPE(Op_t),   INTENT(IN)     :: Op       ! Operator to apply
   TYPE(psi_t),  INTENT(IN)     :: psi      ! Input wavefunction
   TYPE(psi_t),  INTENT(INOUT)  :: Oppsi    ! Resulting wavefunction after application
   TYPE(psi_t)                  :: psi_g    ! Temporary wavefunction on grid
   TYPE(psi_t)                  :: Oppsi_g  ! Temporary result on grid
 
   IF (ALLOCATED(Op%RMat)) THEN
      ! If operator is represented as a matrix, apply it directly
      Oppsi%CVec = MATMUL(Op%RMat, psi%CVec)
 
   ELSE
      IF (psi%Grid) THEN
         ! If psi is on the grid: apply operator directly in grid representation
         CALL init_psi(Oppsi_g, psi%Basis, cplx = .true., grid = .true.)
         IF(Op%Tdef%coord_type == 'CC') THEN
            CALL Calc_Hpsi(psi%CVec, Oppsi_g%CVec, psi%Basis, Op%Scal_pot)
         ELSE
            CALL HPsiGridTdef(psi%CVec, Oppsi_g%CVec, Op)
         END IF

         CALL GridTOBasis_nD_cplx(Oppsi%CVec, Oppsi_g%CVec, psi%Basis)
         CALL dealloc_psi(Oppsi_g)
 
      ELSE
         ! If psi is in the basis: convert to grid, apply operator, and convert back
         CALL init_psi(psi_g,   psi%Basis, cplx = .true., grid = .true.)
         CALL init_psi(Oppsi_g, psi%Basis, cplx = .true., grid = .true.)
 
         CALL BasisTOGrid_nD_cplx(psi_g%CVec,   psi%CVec,   psi%Basis)
         
         IF(Op%Tdef%coord_type == 'CC') THEN
            CALL Calc_Hpsi(psi_g%CVec, Oppsi_g%CVec, psi%Basis, Op%Scal_pot)
         ELSE
            CALL HPsiGridTdef(psi_g%CVec, Oppsi_g%CVec, Op)
         END IF

         CALL GridTOBasis_nD_cplx(Oppsi%CVec, Oppsi_g%CVec, psi%Basis)
 
         CALL dealloc_psi(psi_g)
         CALL dealloc_psi(Oppsi_g)
      END IF
   END IF
 
 END SUBROUTINE calc_OpPsi


 SUBROUTINE Validate_Smolyak_Implementation(H, Psi)

   TYPE(Op_t),  INTENT(IN) :: H
   TYPE(Psi_t), INTENT(IN) :: Psi
   TYPE(Psi_t) :: PsiGrid, HPsiSparse
   COMPLEX(KIND=Rkind), ALLOCATABLE :: FullPsi(:), FullPsiBack(:), SparseCycle(:)
   COMPLEX(KIND=Rkind), ALLOCATABLE :: GridPsi(:), GridHpsi(:), FullHpsi(:), SparseHpsiRef(:)
   REAL(KIND=Rkind) :: NormBasis, NormGrid, EnergySparse, EnergyProjected
   REAL(KIND=Rkind), ALLOCATABLE :: AvgQ(:), WidthQ(:), AvgP(:)
   INTEGER, ALLOCATABLE :: Tab_Iq(:,:)
   INTEGER :: ndim, nsurf, DenseBasisSize, GridSize
   REAL(KIND=Rkind) :: err_expand, err_cycle, err_hpsi
   IF (.NOT. Basis_UsesSmolyakSparseGrid(Psi%Basis)) RETURN
   IF (.NOT. ALLOCATED(Psi%CVec)) RETURN

   ndim = SIZE(Psi%Basis%tab_basis) - 1
   nsurf = Psi%Basis%tab_basis(ndim+1)%nb
   DenseBasisSize = Psi%Basis%nb * nsurf
   GridSize       = Psi%Basis%nq * nsurf

   ALLOCATE(FullPsi(DenseBasisSize), FullPsiBack(SIZE(Psi%CVec)))
   ALLOCATE(GridPsi(GridSize), GridHpsi(GridSize), FullHpsi(DenseBasisSize), SparseHpsiRef(SIZE(Psi%CVec)))
   ALLOCATE(SparseCycle(SIZE(Psi%CVec)))
   ALLOCATE(AvgQ(ndim), WidthQ(ndim), AvgP(ndim))
   CALL ExpandSparseBasisToFullBasis(FullPsi, Psi%CVec, Psi%Basis%SmolyakSparseGrid)
   CALL RestrictFullBasisToSparseBasis(FullPsiBack, FullPsi, Psi%Basis%SmolyakSparseGrid)
   err_expand = MAXVAL(ABS(FullPsiBack - Psi%CVec))
   CALL BasisTOGrid_nD_cplx(GridPsi, Psi%CVec, Psi%Basis)
   CALL GridTOBasis_nD_cplx(SparseCycle, GridPsi, Psi%Basis)
   err_cycle = MAXVAL(ABS(SparseCycle - Psi%CVec))
   CALL init_psi(PsiGrid, Psi%Basis, cplx=.TRUE., grid=.TRUE.)
   PsiGrid%CVec(:) = GridPsi(:)
   CALL Calc_Norm_OF_Psi(Psi, NormBasis)
   CALL Calc_Norm_OF_Psi(PsiGrid, NormGrid)
   CALL init_psi(HPsiSparse, Psi%Basis, cplx=.TRUE., grid=.FALSE.)
   CALL calc_OpPsi(H, Psi, HPsiSparse)

   IF (H%Tdef%coord_type == 'CC') THEN
      CALL Calc_Hpsi(GridPsi, GridHpsi, Psi%Basis, H%Scal_pot)
   ELSE
      CALL HPsiGridTdef(GridPsi, GridHpsi, H)
   END IF

   CALL GridTOBasis_nD_cplx(FullHpsi, GridHpsi, Psi%Basis)
   CALL RestrictFullBasisToSparseBasis(SparseHpsiRef, FullHpsi, Psi%Basis%SmolyakSparseGrid)

   err_hpsi = MAXVAL(ABS(HPsiSparse%CVec - SparseHpsiRef))

    IF (Basis_UsesSmolyakSparseGrid(Psi%Basis)) THEN
       EnergySparse = REAL(DotSmolyakBasis(Psi%CVec,HPsiSparse%CVec,Psi%Basis%SmolyakSparseGrid),KIND=Rkind)
    ELSE
       EnergySparse = REAL(DOT_PRODUCT(Psi%CVec,HPsiSparse%CVec),KIND=Rkind)
    END IF

   EnergyProjected = REAL(DOT_PRODUCT(FullPsi, FullHpsi), KIND=Rkind)

   CALL Calc_tab_Iq0(Tab_Iq, Psi%Basis)
   CALL Calc_AVQ_SQ_nD(Psi, AvgQ, WidthQ, Tab_Iq)
   CALL Calc_Av_imp_k_nD(Psi, AvgP)

   WRITE(out_unit,*) ''

   WRITE(out_unit,*) '================ Validation Smolyak ================'
   CALL PrintSmolyakSparseGridDescriptor(Psi%Basis%SmolyakSparseGrid, out_unit)
   WRITE(out_unit,*) 'Dense basis size            = ', DenseBasisSize
   WRITE(out_unit,*) 'Sparse basis size           = ', SIZE(Psi%CVec)
   WRITE(out_unit,*) 'Dense grid size             = ', GridSize
   WRITE(out_unit,*) '||Restrict(Expand(psi))-psi||inf = ', err_expand
   WRITE(out_unit,*) '||G->B(B->G(psi))-psi||inf       = ', err_cycle
   WRITE(out_unit,*) '||H_sparse psi - P_s H_dense psi||inf = ', err_hpsi
   WRITE(out_unit,*) '<psi|psi> basis             = ', NormBasis
   WRITE(out_unit,*) '<psi|psi> grid              = ', NormGrid
   WRITE(out_unit,*) '<psi|H|psi> sparse          = ', EnergySparse
   WRITE(out_unit,*) '<psi|H|psi> projected dense = ', EnergyProjected
   WRITE(out_unit,*) '<Q>                         = ', AvgQ
   WRITE(out_unit,*) '<P>                         = ', AvgP
   WRITE(out_unit,*) '===================================================='
   WRITE(out_unit,*) ''

   CALL dealloc_psi(PsiGrid)
   CALL dealloc_psi(HPsiSparse)

   DEALLOCATE(FullPsi, FullPsiBack, SparseCycle, GridPsi, GridHpsi, FullHpsi, SparseHpsiRef)
   DEALLOCATE(AvgQ, WidthQ, AvgP, Tab_Iq)

END SUBROUTINE Validate_Smolyak_Implementation

SUBROUTINE Calc_tab_Iq(Tab_Iq,Basis)
   USE  QDUtil_m
   USE Basis_m
   TYPE(Basis_t), intent(in), target             :: Basis
   integer, allocatable ,intent(inout)           :: Tab_iq(:, :)
   integer, allocatable                          :: Tab_iq0(:)
   integer                                       :: ndim, iq,nq
   logical                                       :: Endloop

   ndim = size(Basis%tab_basis) - 1
   nq =Basis%nq

    allocate (Tab_iq(ndim,nq),Tab_iq0(ndim))
   Call Init_tab_ind(Tab_iq0, Basis%NDindexq)
   Iq = 0
   DO
      Iq = Iq + 1
      CALL increase_NDindex(Tab_Iq0, Basis%NDindexq, Endloop)
      IF (Endloop) exit
      Tab_iq(:,Iq) = Tab_Iq0
       !print*,iq,Tab_Iq(:,Iq)
   END DO
   deallocate(Tab_Iq0)
END SUBROUTINE 


SUBROUTINE Calc_Scalar_Pot_openmp(V, Basis,Tab_Iq)
   !$ USE omp_lib
   USE  QDUtil_m
   USE Basis_m
   TYPE(Basis_t), intent(in), target               :: Basis
   real(kind=Rkind), allocatable ,intent(inout)    :: V(:, :, :)
   integer ,intent(in)                             :: Tab_iq(:, :)

   integer                                          :: Ib, ndim, iq,nq,nsurf,maxth
   real(kind=Rkind)               , allocatable     :: Q(:)
   logical                                          :: Endloop

   ndim = size(Basis%tab_basis) - 1
   nq   =Basis%nq
   nsurf=Basis%tab_basis(ndim+1)%nb

    allocate (Q(ndim))
    allocate(V(nq, nsurf, nsurf))
    V(:,:,:) = ZERO
    Q(:) = ZERO

    maxth = 1
    !$ maxth  = omp_get_max_threads()
    write(*,*) 'nbr de coeurs :',maxth
   !$OMP   PARALLEL DEFAULT(NONE) &
   !$OMP   SHARED(Basis,Tab_Iq,maxth,V,nq,ndim,Ib) &
   !$OMP   PRIVATE(Iq,Q) &
   !$OMP   NUM_THREADS(maxth)

   !$OMP BARRIER
   !$OMP   DO SCHEDULE(STATIC)

   DO  Iq = 1,nq

      DO Ib = 1, ndim
          Q(Ib) = Basis%tab_basis(Ib)%x(Tab_Iq(Ib,Iq))
      ENDDO

      !CALL sub_pot(V(Iq, :, :), Q(:), 0)
      CALL sub_Qmodel_V(V(iq, :, :), Q(:))

   END DO
   !$OMP END DO
   !$OMP BARRIER
   !$OMP END PARALLEL

   deallocate(Q)
END SUBROUTINE 



SUBROUTINE test_openmp_op(Basis)
   USE  QDUtil_m
   USE Basis_m
   TYPE(Basis_t), intent(in), target             :: Basis
   real(kind=Rkind), allocatable                 :: V(:, :, :)
   integer, allocatable                          :: Tab_Iq(:, :)
   integer                                        :: nq

  nq = Basis%nq
  print*,"nq=",nq
   call Calc_tab_Iq(Tab_Iq,Basis)
   call  Calc_Scalar_Pot(V, Basis)
   
   !call cpu_time(time=t1)
    !call Calc_Scalar_Pot_openmp(V, Basis,Tab_Iq)
   ! call cpu_time(time=t2)
    !tpsopenmp = t2-t1
    !write(*,*)"temps CPU  openmp :",tpsopenmp
    !deallocate(V)
    !call cpu_time(time=t1)
    !call Calc_Scalar_Pot(V, Basis)
    !call cpu_time(time=t2)
    !tps = t2-t1
    !write(*,*)"temps CPU  non parallele :",tps

END SUBROUTINE 

SUBROUTINE calc_nac(Basis)
   USE  QDUtil_m
   USE Basis_m
   TYPE(Basis_t), intent(in), target                :: Basis
   real(kind=Rkind), allocatable                    :: NAC(:, :, :)

   integer, allocatable                             :: Tab_iq(:)
   integer                                          :: inb, ndim, iq,nq,nsurf,i,nq1,nq2
   real(Kind=Rkind)               , allocatable     :: Q(:),G(:,:,:),V(:,:)
   logical                                          :: Endloop

   ndim = size(Basis%tab_basis) - 1
   nq =Basis%nq
   nsurf=Basis%tab_basis(ndim+1)%nb
   nq1= Basis%tab_basis(1)%nq
   nq2= Basis%tab_basis(2)%nq

   open (unit = 300, file = "N121-ad=t.txt")
   open (unit = 301, file = "N122-ad=t.txt")

   open (unit = 302, file = "N211-ad=t.txt")
   open (unit = 303, file = "N212-ad=t.txt")
   open (unit = 500, file = "x.txt")
   open (unit = 501, file = "y.txt")

    allocate (Q(ndim),Tab_iq(ndim))
    allocate(V(nsurf, nsurf), NAC(nsurf, nsurf,ndim), G(nsurf, nsurf,ndim))

   Call Init_tab_ind(Tab_iq, Basis%NDindexq)
   Iq = 0
   DO
      Iq = Iq + 1
      i =Tab_iq(2)
      CALL increase_NDindex(Tab_iq, Basis%NDindexq, Endloop)
      IF (Endloop) exit
      do inb = 1, ndim
          Q(inb) = Basis%tab_basis(inb)%x(Tab_iq(inb))
      end do
      CALL sub_Qmodel_VG_NAC(V, G, NAC, Q)
      if ( (Tab_iq(2)/= i) )then 
          write(300,*)
          write(301,*)
          write(302,*)
          write(302,*)
      end if  

      write(300,*) Q(:),NAC(1, 2, 1)
      write(301,*) Q(:),NAC(1, 2, 2)

      write(302,*) Q(:),NAC(2, 1, 1)
      write(303,*) Q(:),NAC(2, 1, 2)

   END DO
   do i = 1, nq1
      write(500,*) Basis%tab_basis(1)%x(i)
   end do
   do i = 1, nq2
      write(501,*) Basis%tab_basis(2)%x(i)
   end do

   deallocate(Tab_iq,Q, NAC,G,V)

END SUBROUTINE 




SUBROUTINE calc_VV(Basis)
   USE  QDUtil_m
   USE Basis_m
   TYPE(Basis_t), intent(in), target                :: Basis

   integer, allocatable                             :: Tab_iq(:)
   integer                                          :: inb, ndim, iq,nq,nsurf,i,nq1,nq2
   real(Kind=Rkind)               , allocatable     :: Q(:),V(:,:)
   logical                                          :: Endloop

   ndim = size(Basis%tab_basis) - 1
   nq =Basis%nq
   nsurf=Basis%tab_basis(ndim+1)%nb
   nq1= Basis%tab_basis(1)%nq
   nq2= Basis%tab_basis(2)%nq

   !open (unit = 400, file = "v11-ad=t.txt")
   !open (unit = 401, file = "v22-ad=t.txt")

   !open (unit = 402, file = "v21-ad=t.txt")
   !open (unit = 403, file = "v12-ad=t.txt")
   !open (unit = 404, file = "x.txt")
   !open (unit = 405, file = "y.txt")

    allocate (Q(ndim),Tab_iq(ndim))
    allocate(V(nsurf, nsurf))

   Call Init_tab_ind(Tab_iq, Basis%NDindexq)
   Iq = 0
   DO
      Iq = Iq + 1
      i =Tab_iq(2)
      CALL increase_NDindex(Tab_iq, Basis%NDindexq, Endloop)
      IF (Endloop) exit
      do inb = 1, ndim
          Q(inb) = Basis%tab_basis(inb)%x(Tab_iq(inb))
      end do
      CALL sub_Qmodel_V(V, Q)
      !if ( (Tab_iq(2)/= i) )then 
       !   write(400,*)
        !  write(401,*)
          !write(402,*)
         ! write(403,*)
      !end if  

      !write(400,*) Q(:), V(1, 1)
      !write(401,*) Q(:), V(2, 2)

      !write(402,*) Q(:), V(2, 1)
      !write(403,*) Q(:),V(1, 2)

   END DO
   !do i = 1, nq1
    !  write(404,*) Basis%tab_basis(1)%x(i)
   !end do
   !do i = 1, nq2
    !  write(405,*) Basis%tab_basis(2)%x(i)
   !end do

   deallocate(Tab_iq,Q,V)

END SUBROUTINE 

SUBROUTINE Add_QCoeff_AllSurfaces(AccPsiGrid, CoeffQ, VecPsiGrid, nq, nsurf)
  IMPLICIT NONE
  COMPLEX(KIND=Rkind), INTENT(INOUT), TARGET :: AccPsiGrid(:)
  REAL(KIND=Rkind),    INTENT(IN)            :: CoeffQ(:)
  COMPLEX(KIND=Rkind), INTENT(IN),    TARGET :: VecPsiGrid(:)
  INTEGER,             INTENT(IN)            :: nq, nsurf


  COMPLEX(KIND=Rkind), POINTER :: AccPsiGrid2D(:,:), VecPsiGrid2D(:,:)
  INTEGER :: isurf

  AccPsiGrid2D(1:nq,1:nsurf) => AccPsiGrid
  VecPsiGrid2D(1:nq,1:nsurf) => VecPsiGrid

  DO isurf = 1, nsurf
     AccPsiGrid2D(:,isurf) = AccPsiGrid2D(:,isurf) + CoeffQ(:) * VecPsiGrid2D(:,isurf)
  END DO
END SUBROUTINE Add_QCoeff_AllSurfaces

!============================================================
! Tdef_F2_CrossDeriv_PsiGrid: Apply cross derivative term f2^{i,j}(Q) * (d²psi/dQ_iQ_j)
!
! Computes: f2^{i,j}(Q) * (d²psi/dQ_iQ_j)
!============================================================
SUBROUTINE Tdef_F2_CrossDeriv_PsiGrid(TdefPsiGrid, PsiGrid, H, icoord1, icoord2)

  IMPLICIT NONE
  COMPLEX(KIND=Rkind), INTENT(INOUT), TARGET :: TdefPsiGrid(:)
  COMPLEX(KIND=Rkind), INTENT(IN),    TARGET :: PsiGrid(:)
  TYPE(Op_t),          INTENT(IN)            :: H
  INTEGER,             INTENT(IN)            :: icoord1, icoord2
  COMPLEX(KIND=Rkind), POINTER :: PsiGrid3D(:,:,:)
  COMPLEX(KIND=Rkind), POINTER :: GradPsiGrid3D(:,:,:)
  COMPLEX(KIND=Rkind), POINTER :: GradGradPsiGrid3D(:,:,:)
  COMPLEX(KIND=Rkind), POINTER :: d1gg(:,:)
  COMPLEX(KIND=Rkind), ALLOCATABLE, TARGET :: GradPsiGrid(:)
  COMPLEX(KIND=Rkind), ALLOCATABLE, TARGET :: GradGradPsiGrid(:)

  INTEGER :: iq1, iq2, iq3
  INTEGER :: nq, nq_tot, ndim, nsurf

  IF (Debug_Enabled()) THEN
     WRITE(out_unit,*) 'BEGINNING Tdef_F2_CrossDeriv_PsiGrid'
     FLUSH(out_unit)
  END IF

  ndim  = SIZE(H%Basis%tab_basis) - 1
  nq    = H%Basis%nq
  nsurf = H%Basis%tab_basis(ndim+1)%nb
  nq_tot = nq * nsurf
  ALLOCATE(GradPsiGrid(nq_tot), GradGradPsiGrid(nq_tot))
  GradPsiGrid     = CZERO
  GradGradPsiGrid = CZERO
  ! ---- First derivative with respect to icoord1 ----

  iq1 = 1
  IF (icoord1 > 1) iq1 = PRODUCT(H%Basis%tab_basis(1:icoord1-1)%nq)
  iq2 = H%Basis%tab_basis(icoord1)%nq
  iq3 = nsurf

  IF (icoord1 < ndim) iq3 = iq3 * PRODUCT(H%Basis%tab_basis(icoord1+1:ndim)%nq)

  PsiGrid3D(1:iq1,1:iq2,1:iq3)     => PsiGrid
  GradPsiGrid3D(1:iq1,1:iq2,1:iq3) => GradPsiGrid

  d1gg(1:iq2,1:iq2)                => H%Basis%tab_basis(icoord1)%d1gg(:,:,1)

  DO iq3 = 1, UBOUND(PsiGrid3D,dim=3)
     DO iq1 = 1, UBOUND(PsiGrid3D,dim=1)
        GradPsiGrid3D(iq1,:,iq3) = MATMUL(d1gg, PsiGrid3D(iq1,:,iq3))
     END DO
  END DO

  ! ---- Second derivative with respect to icoord2 ----

  iq1 = 1
  IF (icoord2 > 1) iq1 = PRODUCT(H%Basis%tab_basis(1:icoord2-1)%nq)
  iq2 = H%Basis%tab_basis(icoord2)%nq
  iq3 = nsurf
  IF (icoord2 < ndim) iq3 = iq3 * PRODUCT(H%Basis%tab_basis(icoord2+1:ndim)%nq)

  GradPsiGrid3D(1:iq1,1:iq2,1:iq3)     => GradPsiGrid
  GradGradPsiGrid3D(1:iq1,1:iq2,1:iq3) => GradGradPsiGrid
  d1gg(1:iq2,1:iq2)                    => H%Basis%tab_basis(icoord2)%d1gg(:,:,1)

  DO iq3 = 1, UBOUND(GradPsiGrid3D,dim=3)

     DO iq1 = 1, UBOUND(GradPsiGrid3D,dim=1)

        GradGradPsiGrid3D(iq1,:,iq3) = MATMUL(d1gg, GradPsiGrid3D(iq1,:,iq3))

     END DO

  END DO
  ! ---- Multiply by f2^{i,j}(Q) on EACH electronic surface ----

  CALL Add_QCoeff_AllSurfaces(TdefPsiGrid, H%Tdef%F2(:,icoord1,icoord2), GradGradPsiGrid, nq, nsurf)

  DEALLOCATE(GradPsiGrid, GradGradPsiGrid)

  IF (Debug_Enabled()) THEN
     WRITE(out_unit,*) 'END Tdef_F2_CrossDeriv_PsiGrid'
     FLUSH(out_unit)
  END IF

END SUBROUTINE Tdef_F2_CrossDeriv_PsiGrid

!============================================================
! Tdef_F2_DiagDeriv_PsiGrid: Apply diagonal second derivative d²/dQ_i²
!
! Computes: f2^{i,i}(Q) * (d²psi/dQ_i²)
!============================================================

SUBROUTINE Tdef_F2_DiagDerivPsiGrid(TdefPsiGrid, PsiGrid, H, icoord)

  IMPLICIT NONE
  COMPLEX(KIND=Rkind), INTENT(INOUT), TARGET :: TdefPsiGrid(:)
  COMPLEX(KIND=Rkind), INTENT(IN),    TARGET :: PsiGrid(:)
  TYPE(Op_t),          INTENT(IN)            :: H
  INTEGER,             INTENT(IN)            :: icoord
  COMPLEX(KIND=Rkind), POINTER :: PsiGrid3D(:,:,:)
  COMPLEX(KIND=Rkind), POINTER :: GradGradPsiGrid3D(:,:,:)
  COMPLEX(KIND=Rkind), POINTER :: d2gg(:,:)
  COMPLEX(KIND=Rkind), ALLOCATABLE, TARGET :: GradGradPsiGrid(:)
  INTEGER :: iq1, iq2, iq3
  INTEGER :: nq, nq_tot, ndim, nsurf

  IF (Debug_Enabled()) THEN
     WRITE(out_unit,*) 'BEGINNING Tdef_F2_DiagDerivPsiGrid'
     FLUSH(out_unit)
  END IF

  ndim  = SIZE(H%Basis%tab_basis) - 1
  nq    = H%Basis%nq
  nsurf = H%Basis%tab_basis(ndim+1)%nb
  nq_tot = nq * nsurf
  ALLOCATE(GradGradPsiGrid(nq_tot))

  GradGradPsiGrid = CZERO
  iq1 = 1

  IF (icoord > 1) iq1 = PRODUCT(H%Basis%tab_basis(1:icoord-1)%nq)

  iq2 = H%Basis%tab_basis(icoord)%nq

  iq3 = nsurf

  IF (icoord < ndim) iq3 = iq3 * PRODUCT(H%Basis%tab_basis(icoord+1:ndim)%nq)

  GradGradPsiGrid3D(1:iq1,1:iq2,1:iq3) => GradGradPsiGrid
  PsiGrid3D(1:iq1,1:iq2,1:iq3)         => PsiGrid
  d2gg(1:iq2,1:iq2)                    => H%Basis%tab_basis(icoord)%d2gg(:,:,1,1)

  DO iq3 = 1, UBOUND(PsiGrid3D,dim=3)

     DO iq1 = 1, UBOUND(PsiGrid3D,dim=1)

        GradGradPsiGrid3D(iq1,:,iq3) = MATMUL(d2gg, PsiGrid3D(iq1,:,iq3))

     END DO

  END DO

  CALL Add_QCoeff_AllSurfaces(TdefPsiGrid, H%Tdef%F2(:,icoord,icoord), GradGradPsiGrid, nq, nsurf)

  DEALLOCATE(GradGradPsiGrid)

  IF (Debug_Enabled()) THEN

     WRITE(out_unit,*) 'END Tdef_F2_DiagDerivPsiGrid'

     FLUSH(out_unit)

  END IF

END SUBROUTINE Tdef_F2_DiagDerivPsiGrid



!============================================================
! Tdef_F1_Grad_PsiGrid: Apply first derivative term d/dQ_i
!
! Computes: f1^{i}(Q) * (dpsi/dQ_i)
!============================================================

SUBROUTINE Tdef_F1_Grad_PsiGrid(TdefPsiGrid, PsiGrid, H, icoord)

  IMPLICIT NONE
  COMPLEX(KIND=Rkind), INTENT(INOUT), TARGET :: TdefPsiGrid(:)
  COMPLEX(KIND=Rkind), INTENT(IN),    TARGET :: PsiGrid(:)
  TYPE(Op_t),          INTENT(IN)            :: H
  INTEGER,             INTENT(IN)            :: icoord
  COMPLEX(KIND=Rkind), POINTER :: PsiGrid3D(:,:,:)
  COMPLEX(KIND=Rkind), POINTER :: GradPsiGrid3D(:,:,:)
  COMPLEX(KIND=Rkind), POINTER :: d1gg(:,:)
  COMPLEX(KIND=Rkind), ALLOCATABLE, TARGET :: GradPsiGrid(:)
  INTEGER :: iq1, iq2, iq3
  INTEGER :: nq, nq_tot, ndim, nsurf

  IF (Debug_Enabled()) THEN
     WRITE(out_unit,*) 'BEGINNING Tdef_F1_Grad_PsiGrid'
     FLUSH(out_unit)
  END IF

  ndim  = SIZE(H%Basis%tab_basis) - 1
  nq    = H%Basis%nq
  nsurf = H%Basis%tab_basis(ndim+1)%nb
  nq_tot = nq * nsurf

  ALLOCATE(GradPsiGrid(nq_tot))

  GradPsiGrid = CZERO

  iq1 = 1

  IF (icoord > 1) iq1 = PRODUCT(H%Basis%tab_basis(1:icoord-1)%nq)

  iq2 = H%Basis%tab_basis(icoord)%nq

  iq3 = nsurf

  IF (icoord < ndim) iq3 = iq3 * PRODUCT(H%Basis%tab_basis(icoord+1:ndim)%nq)

  GradPsiGrid3D(1:iq1,1:iq2,1:iq3) => GradPsiGrid

  PsiGrid3D(1:iq1,1:iq2,1:iq3)     => PsiGrid

  d1gg(1:iq2,1:iq2)                => H%Basis%tab_basis(icoord)%d1gg(:,:,1)

  DO iq3 = 1, UBOUND(PsiGrid3D,dim=3)

     DO iq1 = 1, UBOUND(PsiGrid3D,dim=1)

        GradPsiGrid3D(iq1,:,iq3) = MATMUL(d1gg, PsiGrid3D(iq1,:,iq3))

     END DO

  END DO

  CALL Add_QCoeff_AllSurfaces(TdefPsiGrid, H%Tdef%F1(:,icoord), GradPsiGrid, nq, nsurf)

  DEALLOCATE(GradPsiGrid)

  IF (Debug_Enabled()) THEN

     WRITE(out_unit,*) 'END Tdef_F1_Grad_PsiGrid'

     FLUSH(out_unit)

  END IF

END SUBROUTINE Tdef_F1_Grad_PsiGrid


!============================================================
! TdefOnPsiGrid: Main routine for TPsi calculation on grid
!============================================================

SUBROUTINE TdefOnPsiGrid(TdefPsiGrid, PsiGrid, H)

  USE QDUtil_m

  IMPLICIT NONE
  COMPLEX(KIND=Rkind), INTENT(IN),    TARGET :: PsiGrid(:)
  COMPLEX(KIND=Rkind), INTENT(INOUT), TARGET :: TdefPsiGrid(:)
  TYPE(Op_t),          INTENT(IN)            :: H
  INTEGER :: icoord, icoord1, icoord2
  INTEGER :: ndim

  IF (Debug_Enabled()) THEN

     WRITE(out_unit,*) '========================================'

     WRITE(out_unit,*) 'BEGINNING TdefOnPsiGrid'

     WRITE(out_unit,*) 'Computing Tdef Psi on grid'

     WRITE(out_unit,*) '========================================'

     FLUSH(out_unit)

  END IF

  TdefPsiGrid = CZERO

  ndim = SIZE(H%Basis%tab_basis) - 1

  ! Exactly the three formulations documented in the supplied PDF.
  ! f2f1 is the existing expanded F2/F1/Vep implementation.
  SELECT CASE (TRIM(H%Tdef%construction_type))
  CASE ('developed')
     CALL TdefOnPsiGrid_F2F1(TdefPsiGrid, PsiGrid, H, ndim)
  CASE ('divergence')
     CALL TdefOnPsiGrid_DirectFlux(TdefPsiGrid, PsiGrid, H, ndim, 2)
  CASE ('similarity')
     CALL TdefOnPsiGrid_DirectFlux(TdefPsiGrid, PsiGrid, H, ndim, 3)
  CASE DEFAULT
     WRITE(out_unit,*) 'Unknown curvilinear KEO construction: ', TRIM(H%Tdef%construction_type)
     STOP 'TdefOnPsiGrid: invalid construction_type'
  END SELECT

  IF (Debug_Enabled()) THEN

     WRITE(out_unit,*) '======================================='

     WRITE(out_unit,*) 'END TdefOnPsiGrid'

     WRITE(out_unit,*) '======================================='

     FLUSH(out_unit)

  END IF

END SUBROUTINE TdefOnPsiGrid

! Legacy coefficient backend. It is intentionally preserved as a
! regression/reference route and remains based on Tnum-generated F1/F2.
SUBROUTINE TdefOnPsiGrid_F2F1(TdefPsiGrid, PsiGrid, H, ndim)
  COMPLEX(KIND=Rkind), INTENT(IN), TARGET :: PsiGrid(:)
  COMPLEX(KIND=Rkind), INTENT(INOUT), TARGET :: TdefPsiGrid(:)
  TYPE(Op_t), INTENT(IN) :: H
  INTEGER, INTENT(IN) :: ndim
  INTEGER :: icoord, icoord1, icoord2
  DO icoord1=1,ndim
     CALL Tdef_F2_DiagDerivPsiGrid(TdefPsiGrid,PsiGrid,H,icoord1)
     DO icoord2=icoord1+1,ndim
        CALL Tdef_F2_CrossDeriv_PsiGrid(TdefPsiGrid,PsiGrid,H,icoord1,icoord2)
     END DO
  END DO
  DO icoord=1,ndim
     CALL Tdef_F1_Grad_PsiGrid(TdefPsiGrid,PsiGrid,H,icoord)
  END DO
END SUBROUTINE TdefOnPsiGrid_F2F1

! Direct flux implementation. Form=2 is (2c), Form=3 is (2b).
! The derivative operators are applied to complete flux vectors, so mixed
! terms are never counted twice and no hand-expanded F1/F2 is needed.
SUBROUTINE TdefOnPsiGrid_DirectFlux(Tout,Psi,H,ndim,Form)
  COMPLEX(KIND=Rkind), INTENT(IN), TARGET :: Psi(:)
  COMPLEX(KIND=Rkind), INTENT(INOUT), TARGET :: Tout(:)
  TYPE(Op_t), INTENT(IN) :: H
  INTEGER, INTENT(IN) :: ndim,Form
  COMPLEX(KIND=Rkind), ALLOCATABLE :: Work(:,:),Flux(:),State(:),Deriv(:)
  REAL(KIND=Rkind), ALLOCATABLE :: Glocal(:,:)
  INTEGER, ALLOCATABLE :: Tab_iq(:,:)
  REAL(KIND=Rkind), ALLOCATABLE :: Qact(:)
  INTEGER :: nq,nsurf,nall,iq,is,i,j,idx,inb
  LOGICAL :: Endloop
  REAL(KIND=Rkind) :: rho,jac,a,scale,jac_local
  nq=SIZE(H%Tdef%rho); nsurf=H%Basis%tab_basis(ndim+1)%nb
  nall=nq*nsurf
  ALLOCATE(Work(nall,ndim),Flux(nall),State(nall),Deriv(nall))
  ALLOCATE(Tab_iq(ndim,nq), Qact(ndim), Glocal(ndim,ndim))
  CALL Init_tab_ind(Tab_iq(:,1), H%Basis%NDindexq)
  Endloop=.FALSE.
  DO iq=1,nq
     IF (iq > 1) CALL increase_NDindex(Tab_iq(:,iq-1), H%Basis%NDindexq, Endloop)
     IF (Endloop) ERROR STOP 'TdefOnPsiGrid_DirectFlux: inconsistent grid index bounds'
     IF (iq > 1) Tab_iq(:,iq)=Tab_iq(:,iq-1)
  END DO
  ! Retain only the compact integer index table, not a real Q(nq,ndim) array.
  DO iq=1,nq
     DO inb=1,ndim
        Qact(inb)=H%Basis%tab_basis(inb)%x(Tab_iq(inb,iq))
     END DO
     IF (Form==3) THEN
        rho=MAX(H%Tdef%rho(iq),ONETENTH**12)
        CALL SimilarityMetricAtPoint(Glocal,jac_local,H,Qact)
        jac=MAX(jac_local,ONETENTH**12)
        a=SQRT(rho/jac)
        DO is=1,nsurf
           idx=iq+(is-1)*nq; State(idx)=a*Psi(idx)
        END DO
     END IF
  END DO
  IF (Form/=3) THEN
     State=Psi
  END IF
  DO j=1,ndim
     CALL KEO_Grid_FirstDerivative(Work(:,j),State,H,j)
  END DO
  Tout=CZERO
  DO i=1,ndim
     Flux=CZERO
     DO iq=1,nq
        DO inb=1,ndim
           Qact(inb)=H%Basis%tab_basis(inb)%x(Tab_iq(inb,iq))
        END DO
        IF (Form==3) THEN
           CALL SimilarityMetricAtPoint(Glocal,jac_local,H,Qact)
           scale=jac_local
        ELSE
           scale=H%Tdef%rho(iq)
        END IF
        DO is=1,nsurf
           idx=iq+(is-1)*nq
           DO j=1,ndim
              IF (Form==3) THEN
                 Flux(idx)=Flux(idx)+scale*Glocal(i,j)*Work(idx,j)
              ELSE
                 Flux(idx)=Flux(idx)+scale*H%Tdef%G(iq,i,j)*Work(idx,j)
              END IF
           END DO
        END DO
     END DO
     CALL KEO_Grid_FirstDerivative(Deriv,Flux,H,i)
     DO iq=1,nq
        DO inb=1,ndim
           Qact(inb)=H%Basis%tab_basis(inb)%x(Tab_iq(inb,iq))
        END DO
        IF (Form==3) THEN
           rho=MAX(H%Tdef%rho(iq),ONETENTH**12)
           CALL SimilarityMetricAtPoint(Glocal,jac_local,H,Qact)
           jac=MAX(jac_local,ONETENTH**12)
           scale=SQRT(jac/rho)/jac

        ELSE
           scale=ONE/MAX(H%Tdef%rho(iq),ONETENTH**12)
        END IF
        DO is=1,nsurf
           idx=iq+(is-1)*nq
           Tout(idx)=Tout(idx)-HALF*scale*Deriv(idx)
        END DO
     END DO
  END DO
  DEALLOCATE(Work,Flux,State,Deriv,Tab_iq,Qact,Glocal)
END SUBROUTINE TdefOnPsiGrid_DirectFlux

SUBROUTINE SimilarityMetricAtPoint(Glocal,Jlocal,H,Qact)
  REAL(KIND=Rkind), INTENT(OUT) :: Glocal(:,:)
  REAL(KIND=Rkind), INTENT(OUT) :: Jlocal
  TYPE(Op_t), INTENT(IN) :: H
  REAL(KIND=Rkind), INTENT(IN) :: Qact(:)
  TYPE(Type_dnMat) :: dnGG
  TYPE(Type_dnS) :: dnJac
  INTEGER :: ndim, kdet
  ndim=SIZE(Qact)
  CALL alloc_dnSVM(dnGG,H%Tdef%mole_ref%ndimG,H%Tdef%mole_ref%ndimG,H%Tdef%mole_ref%nb_act,0)
  CALL alloc_dnSVM(dnJac,H%Tdef%mole_ref%nb_act,0)
  CALL get_dng_dnGG(Qact,H%Tdef%para_Tnum_ref,H%Tdef%mole_ref,dnGG=dnGG,nderiv=0)
  IF (H%Tdef%mole_ref%ndimG>ndim) THEN
     DO kdet=ndim+1,H%Tdef%mole_ref%ndimG
        dnGG%d0(kdet,:)=ZERO; dnGG%d0(:,kdet)=ZERO; dnGG%d0(kdet,kdet)=ONE
     END DO
  END IF
  Glocal=dnGG%d0(1:ndim,1:ndim)
  CALL Set_ZERO_TO_dnSVM(dnJac)
  CALL sub3_dnDetGG(dnJac,dnGG,0,H%Tdef%mole_ref%masses,H%Tdef%mole_ref%Mtot_inv,H%Tdef%mole_ref%ncart)
  ! sub3_dnDetGG returns the inverse metric volume factor used by Tnum.
  Jlocal=dnJac%d0
  CALL dealloc_dnSVM(dnGG); CALL dealloc_dnSVM(dnJac)
END SUBROUTINE SimilarityMetricAtPoint


! Direct derivative primitive used by the flux kernels.
SUBROUTINE KEO_Grid_FirstDerivative(DPsi, PsiGrid, H, icoord)
  COMPLEX(KIND=Rkind), INTENT(IN), TARGET :: PsiGrid(:)
  COMPLEX(KIND=Rkind), INTENT(OUT), TARGET :: DPsi(:)
  TYPE(Op_t), INTENT(IN) :: H
  INTEGER, INTENT(IN) :: icoord
  COMPLEX(KIND=Rkind), POINTER :: Pin(:,:,:), Pout(:,:,:), D1(:,:)
  INTEGER :: nq, nsurf, n1, n2, n3, k1, k2, k3
  nq=H%Basis%nq
  nsurf=SIZE(PsiGrid)/nq
  n1=1
  IF (icoord>1) THEN
     DO k3=1,icoord-1
        n1=n1*SIZE(H%Basis%tab_basis(k3)%d1gg,1)
     END DO
  END IF
  n2=SIZE(H%Basis%tab_basis(icoord)%d1gg,1)
  ! The electronic surface index is folded into the last dimension,
  n3=(nq/(n1*n2))*nsurf
  Pin(1:n1,1:n2,1:n3)=>PsiGrid; Pout(1:n1,1:n2,1:n3)=>DPsi
  D1(1:n2,1:n2)=>H%Basis%tab_basis(icoord)%d1gg(:,:,1)
  DO k3=1,n3; DO k1=1,n1; Pout(k1,:,k3)=MATMUL(D1,Pin(k1,:,k3)); END DO; END DO
END SUBROUTINE KEO_Grid_FirstDerivative

SUBROUTINE KEO_Grid_SecondDerivative(D2Psi, PsiGrid, H, icoord1, icoord2)
  COMPLEX(KIND=Rkind), INTENT(IN), TARGET :: PsiGrid(:)
  COMPLEX(KIND=Rkind), INTENT(OUT), TARGET :: D2Psi(:)
  TYPE(Op_t), INTENT(IN) :: H
  INTEGER, INTENT(IN) :: icoord1, icoord2
  COMPLEX(KIND=Rkind), ALLOCATABLE, TARGET :: Tmp(:)
  COMPLEX(KIND=Rkind), POINTER :: Pin(:,:,:), Pout(:,:,:), Ptmp(:,:,:), D1(:,:), D2(:,:)
  INTEGER :: nq, nsurf, n1, n2, n3, k1, k2, k3
  nq=H%Basis%nq
  nsurf=SIZE(PsiGrid)/nq
  IF (icoord1==icoord2) THEN
     n1=1
     IF (icoord1>1) THEN
        DO k3=1,icoord1-1
           n1=n1*SIZE(H%Basis%tab_basis(k3)%d1gg,1)
        END DO
     END IF
     n2=SIZE(H%Basis%tab_basis(icoord1)%d1gg,1)
     n3=(nq/(n1*n2))*nsurf
     Pin(1:n1,1:n2,1:n3)=>PsiGrid; Pout(1:n1,1:n2,1:n3)=>D2Psi
     D2(1:n2,1:n2)=>H%Basis%tab_basis(icoord1)%d2gg(:,:,1,1)
     DO k3=1,n3; DO k1=1,n1; Pout(k1,:,k3)=MATMUL(D2,Pin(k1,:,k3)); END DO; END DO
  ELSE
     ALLOCATE(Tmp(SIZE(PsiGrid)))
     CALL KEO_Grid_FirstDerivative(Tmp,PsiGrid,H,icoord1)
     CALL KEO_Grid_FirstDerivative(D2Psi,Tmp,H,icoord2)
     DEALLOCATE(Tmp)
  END IF
END SUBROUTINE KEO_Grid_SecondDerivative

!==============================================================================
! Initialize Hamiltonian operator with kinetic terms (F1, F2) and effective potential
! Computes these quantities on a grid for later use in dynamics calculations
!
! Input:
!   Op        - Hamiltonian operator structure (inout)
!   para_Tnum - KEO parameters (in)
!   mole      - Molecular geometry (in)
!
! Output stored in Op%Tdef:
!   F2   - ndim x ndim kinetic tensor (nq, ndim, ndim)
!   F1   - ndim kinetic vector (nq, ndim)
!   rho  - Metric determinant (nq)
!   V    - Scalar potential, optionally including Vep (nq, nsurf, nsurf)
!==============================================================================
SUBROUTINE Init_Op_Tdef_F2_F1_Vep_FromTnum(H, para_Tnum, mole)

   USE mod_Coord_KEO
   USE mod_PrimOp
   USE mod_dnSVM
   USE mod_dnRho
   USE mod_Constant
   USE TnumTana_system_m
   
   IMPLICIT NONE

   ! Arguments
   TYPE(Op_t),      INTENT(inout) :: H
   TYPE(Tnum),      INTENT(in)    :: para_Tnum
   TYPE(CoordType), INTENT(in)    :: mole

   ! Local variables
   INTEGER :: iq, ndim, nq, nsurf,isurf,kdet,inb
   INTEGER, ALLOCATABLE :: Tab_iq(:)
   LOGICAL :: Endloop
   REAL(kind=Rkind), allocatable :: Qact(:)     ! Current active coordinates

   ! Kinetic energy terms
   REAL(kind=Rkind), allocatable :: Tdef2(:,:), Tdef1(:)
   ! Coriolis and rotational terms (computed but not stored)
   REAL(kind=Rkind), allocatable :: Tcor2(:,:), Tcor1(:), Trot(:,:)
   REAL(kind=Rkind)              :: Vep, rho    ! Extrapotential and metric det
   TYPE(Type_dnMat)               :: dnGG
   TYPE(Type_dnS)                 :: dnJac, dnrho

   !--------------------------------------------------------------------
   ! 1. Get dimensions from basis set
   !    Last basis element = electronic surfaces
   !--------------------------------------------------------------------
   ndim = size(H%Basis%tab_basis) - 1
   nsurf = H%Basis%tab_basis(ndim+1)%nb
   nq = H%Basis%nq

   ! Safety checks
   IF (ndim <= 0) THEN
      WRITE(*,*) 'Error: ndim <= 0 in Init_Op_Tdef_F2_F1_Vep_FromTnum'
      STOP
   END IF

   !--------------------------------------------------------------------
   ! 2. Allocate local workspace
   !--------------------------------------------------------------------
   ALLOCATE(Qact(ndim))
   ALLOCATE(Tab_iq(ndim))
   CALL Init_tab_ind(Tab_iq, H%Basis%NDindexq)
   
   ! Arrays for KEO calculation
   ALLOCATE(Tdef2(ndim,ndim), Tdef1(ndim))
   ALLOCATE(Tcor2(ndim,3), Tcor1(3), Trot(3,3))  ! Not used
   CALL alloc_dnSVM(dnGG,mole%ndimG,mole%ndimG,mole%nb_act,0)
   CALL alloc_dnSVM(dnJac,mole%nb_act,0)
   CALL alloc_dnSVM(dnrho,mole%nb_act,0)

   !--------------------------------------------------------------------
   ! 3. Allocate only the arrays required by the selected formulation.
   ! developed : F1, F2, V
   ! divergence: G, rho, V
   ! similarity: G, rho, J, V
   !--------------------------------------------------------------------
   IF (ALLOCATED(H%Tdef%F1)) DEALLOCATE(H%Tdef%F1)
   IF (ALLOCATED(H%Tdef%F2)) DEALLOCATE(H%Tdef%F2)
   IF (ALLOCATED(H%Tdef%G))  DEALLOCATE(H%Tdef%G)
   IF (ALLOCATED(H%Tdef%rho)) DEALLOCATE(H%Tdef%rho)
   IF (ALLOCATED(H%Tdef%J)) DEALLOCATE(H%Tdef%J)
   IF (ALLOCATED(H%Tdef%V)) DEALLOCATE(H%Tdef%V)

   SELECT CASE (TRIM(H%Tdef%construction_type))
   CASE ('developed')
      ALLOCATE(H%Tdef%F1(nq,ndim), H%Tdef%F2(nq,ndim,ndim))
   CASE ('divergence')
      ALLOCATE(H%Tdef%G(nq,ndim,ndim), H%Tdef%rho(nq))
   CASE ('similarity')
      ! The similarity kernel evaluates G and J on demand.
      ALLOCATE(H%Tdef%rho(nq))
   CASE DEFAULT
      WRITE(out_unit,*) 'Unknown curvilinear KEO construction: ', TRIM(H%Tdef%construction_type)
      STOP 'Init_Op_Tdef: invalid construction_type'
   END SELECT
   ALLOCATE(H%Tdef%V(nq,nsurf,nsurf))
   IF (ALLOCATED(H%Tdef%F2)) H%Tdef%F2=ZERO
   IF (ALLOCATED(H%Tdef%F1)) H%Tdef%F1=ZERO
   IF (ALLOCATED(H%Tdef%rho)) H%Tdef%rho=ZERO
   IF (ALLOCATED(H%Tdef%G)) H%Tdef%G=ZERO
   IF (ALLOCATED(H%Tdef%J)) H%Tdef%J=ZERO
   H%Tdef%V=ZERO
   !--------------------------------------------------------------------
   ! 4. Compute scalar potential on grid
   !--------------------------------------------------------------------
   CALL Calc_Scalar_Pot(H%Scal_pot, H%Basis, mole, H%Tdef%pot_itQtransfo)

   !--------------------------------------------------------------------
   ! 5. Loop over grid points
   !--------------------------------------------------------------------
   DO iq = 1, nq
      ! Recover the current grid point without materialising Q(nq,ndim).
      ! increase_NDindex follows the same Fortran ordering as Basis%x.
      CALL increase_NDindex(Tab_iq, H%Basis%NDindexq, Endloop)
      IF (Endloop) ERROR STOP 'Init_Op_Tdef: inconsistent grid index bounds'
      DO inb = 1, ndim
         Qact(inb) = H%Basis%tab_basis(inb)%x(Tab_iq(inb))
      END DO

      ! Compute only the data required by the selected formulation.
      ! The similarity branch must not evaluate the legacy F1/F2/Vep path:
      ! it obtains G(Q) and J(Q) in the direct-flux kernel at action time.
      SELECT CASE (TRIM(H%Tdef%construction_type))
      CASE ('developed','divergence')
         CALL calc3_f2_f1Q_num( Qact, Tdef2, Tdef1, Vep, rho, &
                                Tcor2, Tcor1, Trot, para_Tnum, mole )
      CASE ('similarity')
         CALL get_dng_dnGG(Qact,para_Tnum,mole,dnGG=dnGG,nderiv=0)
         CALL Set_ZERO_TO_dnSVM(dnJac)
         IF (para_Tnum%nrho == 0 .AND. .NOT. para_Tnum%Gcte) THEN
            CALL sub3_dndetGG(dnJac,dnGG,0,mole%masses,mole%Mtot_inv,mole%ncart)
         END IF
         CALL sub3_dnrho(dnrho,dnJac,Qact,mole,0,para_Tnum%num_x, &
                         para_Tnum%stepT,para_Tnum%nrho)
         rho = dnrho%d0
      END SELECT

      ! Store metric data only for the divergence backend. The similarity
      ! backend deliberately keeps G(Q) and J(Q) transient and recomputes
      ! them in SimilarityMetricAtPoint during each H-application.
      IF (TRIM(H%Tdef%construction_type) == 'divergence') THEN
         ! Tnum stores the full metric; retain only active coordinates.
         CALL get_dng_dnGG(Qact,para_Tnum,mole,dnGG=dnGG,nderiv=0)
         IF (mole%ndimG > ndim) THEN
            DO kdet=ndim+1,mole%ndimG
               dnGG%d0(kdet,:) = ZERO; dnGG%d0(:,kdet) = ZERO
               dnGG%d0(kdet,kdet) = ONE
            END DO
         END IF
         H%Tdef%G(iq,:,:) = dnGG%d0(1:ndim,1:ndim)
      END IF
      !stop 'premier iteration'
      ! Store results
      IF (ALLOCATED(H%Tdef%F2)) H%Tdef%F2(iq,:,:) = Tdef2(:,:)
      IF (ALLOCATED(H%Tdef%F1)) H%Tdef%F1(iq,:)   = Tdef1(:)
      IF (ALLOCATED(H%Tdef%rho)) H%Tdef%rho(iq)   = rho
      
      IF (Debug_Enabled()) THEN
         IF (ALLOCATED(H%Tdef%F1)) THEN
            WRITE(*,*) ' '
            CALL Write_VecMat(H%Tdef%F1(iq,:), out_unit, ndim, info='F1')
            WRITE(*,*) ' '
         END IF
      END IF

      
      ! Effective potential = scalar potential + extrapotential
      ! Note: Vep assumed same for all electronic surfaces
      H%Tdef%V(iq, :, :) = H%Scal_pot(iq,:,:)   ! physical potential
      IF (TRIM(H%Tdef%construction_type) /= 'similarity') THEN
         DO isurf = 1, nsurf
            H%Tdef%V(iq,isurf,isurf) = H%Tdef%V(iq,isurf,isurf) + Vep
         END DO
      END IF

   END DO
     
      IF (Debug_Enabled()) THEN
         IF (ALLOCATED(H%Tdef%rho)) THEN
            WRITE(*,*) ' '
            CALL Write_VecMat(H%Tdef%rho, out_unit, ndim, info='Rho')
            WRITE(*,*) ' '
         END IF
      END IF

   !--------------------------------------------------------------------
   ! 6. Clean up local arrays
   !--------------------------------------------------------------------
   CALL dealloc_dnSVM(dnGG)
   CALL dealloc_dnSVM(dnJac)
   CALL dealloc_dnSVM(dnrho)
   DEALLOCATE(Tab_iq, Qact, Tdef2, Tdef1, Tcor2, Tcor1, Trot,H%Scal_pot)

END SUBROUTINE Init_Op_Tdef_F2_F1_Vep_FromTnum

 SUBROUTINE InitQModelBasisParam(nio, qmodel_q0, qmodel_sigma, debug_output)
      USE  QDUtil_m
      USE mod_Constant
      USE TnumTana_system_m
      USE mod_Coord_KEO
      USE mod_PrimOp
      integer, intent(in)                             ::nio
      REAL(KIND=Rkind), ALLOCATABLE, OPTIONAL, INTENT(OUT) :: qmodel_q0(:), qmodel_sigma(:)
      REAL(KIND=Rkind), ALLOCATABLE                   :: GGdef(:,:)
      LOGICAL                                         :: debug
      LOGICAL, OPTIONAL, INTENT(IN)                 :: debug_output
     Integer                                          :: ndim,nsurf,option,icoord,isurf,jcoord
     Logical                                          :: adiabatic
     Character (len=16)                               :: pot_name
     Real(kind=Rkind), allocatable                    :: V(:,:)
     Real(kind=Rkind), allocatable                    :: Hessian(:,:,:,:),Q0(:),Q0_active(:)
     Real(kind=Rkind), allocatable                    :: G(:,:,:)
     Real(kind=Rkind), allocatable                    :: ScaleQ(:),M(:),sigma(:),a(:)
     REAL(KIND=Rkind), ALLOCATABLE                    :: Gkin(:,:),MapJac(:,:),MapHess(:,:,:),Hess_active(:,:)
     TYPE(constant)                                   :: const_phys
     TYPE(CoordType)                                  :: mole
     TYPE(Tnum)                                       :: para_Tnum
     TYPE(PrimOp_t)                                   :: PrimOp
     TYPE(Type_dnMat)                                 :: dnGG


    ! Paramètres de la base
      INTEGER :: nb_basis
      INTEGER :: nq_basis
      ! Fichier
      INTEGER :: unit_file
      INTEGER :: ios
      CHARACTER(LEN=200) :: iomsg
      INTEGER(KIND=Ikind) :: pot_itQtransfo
      CHARACTER(LEN=Name_len) :: coord_type,vcoord_type,construction_type,sym_type,model_type
      REAL(KIND=Rkind) :: pot0
      NAMELIST /coordQMLpot/ coord_type,vcoord_type,construction_type,sym_type,model_type,pot0,pot_itQtransfo

      
      debug = Debug_Enabled() .OR. Basis_Debug_Active()
      IF (PRESENT(debug_output)) debug = debug .OR. debug_output
      IF (debug) THEN
         !write(out_unit,*) 'BEGINNING GetBasisParameter'
         flush (out_unit)
      END IF

      ! Valeurs initiales demandées
      nb_basis = 2
      nq_basis = 2
      coord_type = 'CC'
      vcoord_type = ''
      construction_type = 'developed'
      sym_type = 'SYM'
      model_type = ''
      pot0 = ZERO
      pot_itQtransfo = -1_Ikind

      REWIND(nio)
      READ(nio,NML=coordQMLpot,IOSTAT=ios)
      IF (ios /= 0) pot_itQtransfo = -1_Ikind
      REWIND(nio)

      ndim = -1
      nsurf = -1
      pot_name = 'read_model'
      adiabatic = .false.
      option = 1
      

      Call sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)

      IF (debug) WRITE(out_unit, *) 'ndim, nsurf', ndim, nsurf
      IF (debug) WRITE(out_unit, *) 'pot_name'
      
      Allocate (GGdef(ndim , ndim))
      Allocate(Hessian(nsurf,nsurf,ndim,ndim))
      Allocate(G(nsurf,nsurf,ndim))
      Allocate(Q0(ndim))
      Allocate(M(ndim))
      Allocate(ScaleQ(ndim),sigma(ndim),a(ndim))
      Allocate(V(nsurf,nsurf))

      Call get_Qmodel_Q0(Q0,option)
      CALL get_Qmodel_GGdef(GGdef)
      IF (debug) CALL Write_VecMat(GGdef, out_unit, ndim, info='<GGdef>')
      IF (debug) write(out_unit,*)''
      IF (debug) CALL Write_VecMat(Q0, out_unit, ndim, info='<Q0>')

      Call sub_Qmodel_VGH(V,G,Hessian,Q0)
      isurf=1
      IF (debug) CALL Write_VecMat(G(isurf,isurf,:), out_unit, ndim, info='<grad V(Q0), surface 1>')
      IF (debug) CALL Write_VecMat(Hessian(isurf,isurf,:,:), out_unit, ndim, info='<Hessian V(Q0), surface 1>')
      IF (MAXVAL(ABS(G(isurf,isurf,:))) > 1.E-8_Rkind) THEN
         WRITE(out_unit,*) 'WARNING InitQModelBasisParam: Q0 has a non-negligible gradient; ', &
                           'the diagonal harmonic parameters are local curvatures, not a minimum fit.'
      END IF

      ALLOCATE(Q0_active(ndim),MapJac(ndim,ndim),MapHess(ndim,ndim,ndim),Hess_active(ndim,ndim))
      ALLOCATE(Gkin(ndim,ndim))
      IF (TRIM(coord_type) == 'CC') THEN
         ! The direct QML/CC path is already expressed in the model coordinates.
         ! Use its own kinetic metric instead of unconditionally requiring a
         ! Tnum molecule/geometry namelist, which is not part of reduced QML decks.
         IF (pot_itQtransfo >= 0_Ikind) &
            ERROR STOP 'InitQModelBasisParam: pot_itQtransfo requires coord_type=CURV.'
         Q0_active = Q0
         Hess_active = Hessian(isurf,isurf,:,:)
         Gkin = GGdef
         IF (debug) THEN
            WRITE(out_unit,*) 'QML GWP parameters and kinetic metric used directly in CC coordinates.'
            CALL Write_VecMat(Gkin,out_unit,ndim,info='<QML CC kinetic metric>')
         END IF
      ELSE
         REWIND(nio)
         CALL Read_CoordType(mole,para_Tnum,const_phys)
         CALL read_RefGeom(mole,para_Tnum)
         CALL Finalize_TnumTana_Coord_PrimOp(para_Tnum,mole,PrimOp)
         IF (mole%nb_act /= ndim .OR. mole%ndimG < ndim) THEN
            WRITE(out_unit,*) 'ERROR InitQModelBasisParam: QML and Tnum active dimensions differ ', &
                              ndim,mole%nb_act,mole%ndimG
            STOP 'InitQModelBasisParam: incompatible QML/Tnum coordinates'
         END IF
         IF (pot_itQtransfo < -1_Ikind .OR. pot_itQtransfo > mole%nb_Qtransfo) THEN
            WRITE(out_unit,*) 'ERROR InitQModelBasisParam: invalid pot_itQtransfo ',pot_itQtransfo, &
                              ' valid range is -1..',mole%nb_Qtransfo
            STOP 'InitQModelBasisParam: invalid QML coordinate stage'
         END IF
         IF (pot_itQtransfo >= 0_Ikind) THEN
            CALL Invert_Tnum_Stage_Map(Q0,pot_itQtransfo,mole,Q0_active,MapJac,MapHess)
            Hess_active = MATMUL(TRANSPOSE(MapJac), &
                                 MATMUL(Hessian(isurf,isurf,:,:),MapJac))
            DO icoord = 1, ndim
               DO jcoord = 1, ndim
                  Hess_active(icoord,jcoord) = Hess_active(icoord,jcoord) + &
                       DOT_PRODUCT(G(isurf,isurf,:),MapHess(:,icoord,jcoord))
               END DO
            END DO
            WRITE(out_unit,*) 'QML GWP parameters mapped from Tnum stage ',pot_itQtransfo,' to active coordinates.'
         ELSE
            Q0_active = Q0
            Hess_active = Hessian(isurf,isurf,:,:)
         END IF
         Q0 = Q0_active
         CALL alloc_dnSVM(dnGG,mole%ndimG,mole%ndimG,mole%nb_act,0)
         CALL get_dng_dnGG(Q0,para_Tnum,mole,dnGG=dnGG,nderiv=0)
         Gkin(:,:) = dnGG%d0(1:ndim,1:ndim)
         CALL Write_VecMat(Gkin,out_unit,ndim,info='<Tnum Gkin(Q0)>')
         CALL dealloc_dnSVM(dnGG)
         CALL dealloc_CoordType(mole)
      END IF

      DO icoord = 1, ndim
       IF (ABS(Gkin(icoord,icoord)) <= TINY(ONE)) THEN
         WRITE(out_unit,*) 'ERROR InitQModelBasisParam: zero Tnum kinetic metric at coordinate ',icoord
         STOP 'InitQModelBasisParam: zero Tnum kinetic metric'
       END IF
       M(icoord)      = ONE/Gkin(icoord,icoord)
       IF (M(icoord) <= ZERO) THEN
         WRITE(out_unit,*) 'ERROR InitQModelBasisParam: non-positive effective mass at coordinate ',icoord
         STOP 'InitQModelBasisParam: non-positive effective mass'
       END IF
       IF (Hess_active(icoord,icoord)/M(icoord) <= ZERO) THEN
         WRITE(out_unit,*) 'ERROR InitQModelBasisParam: non-positive diagonal curvature at coordinate ',icoord
         STOP 'InitQModelBasisParam: non-positive Hessian curvature'
       END IF
       ! Keep a coordinate-aligned product basis: use diagonal elements of the
       ! selected kinetic metric and local Hessian. Off-diagonal couplings are
       ! retained by the Hamiltonian but not by this product-basis estimate.
       ScaleQ(icoord) = SQRT(M(icoord)*SQRT(Hess_active(icoord,icoord)/M(icoord)))
       sigma(icoord) =SQRT(TWO)/ ScaleQ(icoord)
       a(icoord) = ScaleQ(icoord)*ScaleQ(icoord)
      END DO

      IF (PRESENT(qmodel_q0)) THEN
         ALLOCATE(qmodel_q0(ndim))
         qmodel_q0 = Q0
      END IF
      IF (PRESENT(qmodel_sigma)) THEN
         ALLOCATE(qmodel_sigma(ndim))
         qmodel_sigma = sigma
      END IF

      write(out_unit,*)''
      CALL Write_VecMat(M, out_unit, ndim, info='<M>')
      write(out_unit,*)''
      CALL Write_VecMat(ScaleQ, out_unit, ndim, info='<ScaleQ>')
      write(out_unit,*)''
      CALL Write_VecMat(sigma, out_unit, ndim, info=' <sigma>')
      write(out_unit,*)''
      CALL Write_VecMat(a, out_unit, ndim, info='<Re(alpha)>')

   ! ============================================================
   ! Ouverture du fichier de sortie
   ! ============================================================

   OPEN(NEWUNIT=unit_file,FILE='qmodel_basis.txt',STATUS='REPLACE',ACTION='WRITE',FORM='FORMATTED',IOSTAT=ios,IOMSG=iomsg)

   IF (ios /= 0) THEN
      WRITE(out_unit, *) &
         'Erreur lors de l''ouverture de qmodel_basis.txt : ', &
         TRIM(iomsg)
      STOP
   END IF

   ! ============================================================
   ! Écriture dans qmodel_basis.txt
   ! ============================================================
   DO icoord = 1, ndim
      ! Ligne de définition de la gaussienne initiale
      WRITE(unit_file, *) &
         '&defWP0 sigma=', sigma(icoord), &
         ' Beta=0.0 Qeq=', Q0(icoord), &
         ' imp_k=0.0 gamma=0.0 /'
      ! Ligne de définition de la base
      WRITE(unit_file, *) &
         "&basis_nd name='HO' nb=", nb_basis, &
         " nq=", nq_basis, &
         " Q0=", Q0(icoord), &
         " SCALEQ=", ScaleQ(icoord), &
         " Imp_k=0.0 alpha=(", a(icoord), ",0.0) /"
      ! Ligne vide entre deux coordonnées
      WRITE(unit_file, *)
   END DO
   ! Fermeture du fichier
   CLOSE(unit_file)

   DEALLOCATE(Gkin)
      DEALLOCATE(GGdef,Hessian,G,Q0,Q0_active,M,ScaleQ,sigma,a,V,MapJac,MapHess,Hess_active)
   IF (debug) THEN
      ! write(out_unit,*) 'END InitQModelBasisParam'
      flush (out_unit)
   END IF
END SUBROUTINE   InitQModelBasisParam

SUBROUTINE Map_Active_To_Tnum_Stage(Qactive,it_stage,mole,Qstage)
   USE QDUtil_m, ONLY: Rkind
   USE mod_Coord_KEO, ONLY: CoordType
   IMPLICIT NONE
   REAL(KIND=Rkind), INTENT(IN) :: Qactive(:)
   INTEGER(KIND=Ikind), INTENT(IN) :: it_stage
   TYPE(CoordType), INTENT(IN) :: mole
   REAL(KIND=Rkind), ALLOCATABLE, INTENT(OUT) :: Qstage(:)

   CALL sub_QactTOQit(Qactive,Qstage,it_stage,mole,.FALSE.)
   IF (SIZE(Qstage) /= SIZE(Qactive)) THEN
      WRITE(out_unit,*) 'QML/Tnum stage dimensions are not square/invertible: ',SIZE(Qstage),SIZE(Qactive)
      ERROR STOP 'InitQModelBasisParam: cannot map Gaussian through a dimension-changing Tnum stage.'
   END IF
END SUBROUTINE Map_Active_To_Tnum_Stage

SUBROUTINE Invert_Tnum_Stage_Map(Qstage,it_stage,mole,Qactive,MapJac,MapHess)
   USE QDUtil_m, ONLY: Rkind,ONE,TWO,ONETENTH
   USE mod_Coord_KEO, ONLY: CoordType
   IMPLICIT NONE
   REAL(KIND=Rkind), INTENT(IN) :: Qstage(:)
   INTEGER(KIND=Ikind), INTENT(IN) :: it_stage
   TYPE(CoordType), INTENT(IN) :: mole
   REAL(KIND=Rkind), INTENT(OUT) :: Qactive(:),MapJac(:,:),MapHess(:,:,:)
   REAL(KIND=Rkind), ALLOCATABLE :: Qmapped(:),Qtrial(:),Qmapped_trial(:),residual(:),delta(:),Jac(:,:)
   REAL(KIND=Rkind) :: resid_norm,trial_norm,alpha,tolerance
   INTEGER :: ndim,iteration,line_search
   LOGICAL :: converged,accepted,solve_ok

   ndim = SIZE(Qstage)
   IF (SIZE(Qactive) /= ndim .OR. SIZE(MapJac,1) /= ndim .OR. SIZE(MapJac,2) /= ndim .OR. &
       SIZE(MapHess,1) /= ndim .OR. SIZE(MapHess,2) /= ndim .OR. SIZE(MapHess,3) /= ndim) &
      ERROR STOP 'Invert_Tnum_Stage_Map: inconsistent array dimensions.'

   ALLOCATE(Jac(ndim,ndim))
   ! Tnum provides active -> requested-stage mapping. Invert it locally with
   ! a damped Newton solve, valid for any smooth square transformation.
   Qactive = Qstage
   tolerance = ONETENTH**10*MAX(ONE,MAXVAL(ABS(Qstage)))
   converged = .FALSE.
   DO iteration = 1,60
      CALL Map_Active_To_Tnum_Stage(Qactive,it_stage,mole,Qmapped)
      residual = Qstage-Qmapped
      resid_norm = MAXVAL(ABS(residual))
      IF (resid_norm <= tolerance) THEN
         converged = .TRUE.
         EXIT
      END IF
      CALL Calc_Tnum_Stage_Jacobian(Qactive,it_stage,mole,Jac)
      CALL SolveDenseSystem_GWP(Jac,residual,delta,solve_ok)
      IF (.NOT. solve_ok) THEN
         WRITE(out_unit,*) 'Singular Tnum Jacobian while mapping QML minimum; stage=',it_stage
         ERROR STOP 'InitQModelBasisParam: selected Tnum stage is not locally invertible.'
      END IF

      alpha = ONE
      accepted = .FALSE.
      DO line_search = 1,24
         Qtrial = Qactive+alpha*delta
         CALL Map_Active_To_Tnum_Stage(Qtrial,it_stage,mole,Qmapped_trial)
         trial_norm = MAXVAL(ABS(Qstage-Qmapped_trial))
         IF (trial_norm < resid_norm) THEN
            Qactive = Qtrial
            accepted = .TRUE.
            EXIT
         END IF
         alpha = alpha/TWO
      END DO
      IF (.NOT. accepted) THEN
         WRITE(out_unit,*) 'Tnum inverse did not decrease the coordinate residual; stage=',it_stage, &
                           ' residual=',resid_norm
         ERROR STOP 'InitQModelBasisParam: failed to invert the selected Tnum coordinate stage.'
      END IF
   END DO
   IF (.NOT. converged) THEN
      CALL Map_Active_To_Tnum_Stage(Qactive,it_stage,mole,Qmapped)
      resid_norm = MAXVAL(ABS(Qstage-Qmapped))
      IF (resid_norm > tolerance) THEN
         WRITE(out_unit,*) 'Tnum inverse did not converge; stage=',it_stage,' residual=',resid_norm
         ERROR STOP 'InitQModelBasisParam: failed to converge while inverting Tnum coordinates.'
      END IF
   END IF

   CALL Calc_Tnum_Stage_Jacobian(Qactive,it_stage,mole,MapJac)
   CALL Calc_Tnum_Stage_Hessian(Qactive,it_stage,mole,MapHess)
   IF (Debug_Enabled() .OR. Basis_Debug_Active()) THEN
      CALL Map_Active_To_Tnum_Stage(Qactive,it_stage,mole,Qmapped)
      WRITE(out_unit,*) 'Tnum inverse coordinate residual: ',MAXVAL(ABS(Qstage-Qmapped))
   END IF
END SUBROUTINE Invert_Tnum_Stage_Map

SUBROUTINE Calc_Tnum_Stage_Jacobian(Qactive,it_stage,mole,Jac)
   USE QDUtil_m, ONLY: Rkind,ONE,TWO,ONETENTH
   USE mod_Coord_KEO, ONLY: CoordType
   IMPLICIT NONE
   REAL(KIND=Rkind), INTENT(IN) :: Qactive(:)
   INTEGER(KIND=Ikind), INTENT(IN) :: it_stage
   TYPE(CoordType), INTENT(IN) :: mole
   REAL(KIND=Rkind), INTENT(OUT) :: Jac(:,:)
   REAL(KIND=Rkind), ALLOCATABLE :: Qplus(:),Qminus(:),Fplus(:),Fminus(:)
   REAL(KIND=Rkind) :: h
   INTEGER :: i,ndim

   ndim = SIZE(Qactive)
   IF (SIZE(Jac,1) /= ndim .OR. SIZE(Jac,2) /= ndim) &
      ERROR STOP 'Calc_Tnum_Stage_Jacobian: inconsistent array dimensions.'
   DO i = 1,ndim
      h = ONETENTH**5*MAX(ONE,ABS(Qactive(i)))
      Qplus = Qactive; Qminus = Qactive
      Qplus(i) = Qplus(i)+h
      Qminus(i) = Qminus(i)-h
      CALL Map_Active_To_Tnum_Stage(Qplus,it_stage,mole,Fplus)
      CALL Map_Active_To_Tnum_Stage(Qminus,it_stage,mole,Fminus)
      Jac(:,i) = (Fplus-Fminus)/(TWO*h)
   END DO
END SUBROUTINE Calc_Tnum_Stage_Jacobian

SUBROUTINE Calc_Tnum_Stage_Hessian(Qactive,it_stage,mole,MapHess)
   USE QDUtil_m, ONLY: Rkind,ONE,TWO,FOUR,ONETENTH
   USE mod_Coord_KEO, ONLY: CoordType
   IMPLICIT NONE
   REAL(KIND=Rkind), INTENT(IN) :: Qactive(:)
   INTEGER(KIND=Ikind), INTENT(IN) :: it_stage
   TYPE(CoordType), INTENT(IN) :: mole
   REAL(KIND=Rkind), INTENT(OUT) :: MapHess(:,:,:)
   REAL(KIND=Rkind), ALLOCATABLE :: Qp(:),Qm(:),Qpp(:),Qpm(:),Qmp(:),Qmm(:)
   REAL(KIND=Rkind), ALLOCATABLE :: F0(:),Fp(:),Fm(:),Fpp(:),Fpm(:),Fmp(:),Fmm(:)
   REAL(KIND=Rkind), ALLOCATABLE :: h(:)
   INTEGER :: i,j,ndim

   ndim = SIZE(Qactive)
   IF (SIZE(MapHess,1) /= ndim .OR. SIZE(MapHess,2) /= ndim .OR. SIZE(MapHess,3) /= ndim) &
      ERROR STOP 'Calc_Tnum_Stage_Hessian: inconsistent array dimensions.'
   ALLOCATE(h(ndim))
   DO i=1,ndim
      h(i) = ONETENTH**3*MAX(ONE,ABS(Qactive(i)))
   END DO
   CALL Map_Active_To_Tnum_Stage(Qactive,it_stage,mole,F0)
   MapHess = ZERO
   DO i=1,ndim
      Qp=Qactive; Qm=Qactive
      Qp(i)=Qp(i)+h(i); Qm(i)=Qm(i)-h(i)
      CALL Map_Active_To_Tnum_Stage(Qp,it_stage,mole,Fp)
      CALL Map_Active_To_Tnum_Stage(Qm,it_stage,mole,Fm)
      MapHess(:,i,i)=(Fp-TWO*F0+Fm)/(h(i)*h(i))
      DO j=i+1,ndim
         Qpp=Qactive; Qpm=Qactive; Qmp=Qactive; Qmm=Qactive
         Qpp(i)=Qpp(i)+h(i); Qpp(j)=Qpp(j)+h(j)
         Qpm(i)=Qpm(i)+h(i); Qpm(j)=Qpm(j)-h(j)
         Qmp(i)=Qmp(i)-h(i); Qmp(j)=Qmp(j)+h(j)
         Qmm(i)=Qmm(i)-h(i); Qmm(j)=Qmm(j)-h(j)
         CALL Map_Active_To_Tnum_Stage(Qpp,it_stage,mole,Fpp)
         CALL Map_Active_To_Tnum_Stage(Qpm,it_stage,mole,Fpm)
         CALL Map_Active_To_Tnum_Stage(Qmp,it_stage,mole,Fmp)
         CALL Map_Active_To_Tnum_Stage(Qmm,it_stage,mole,Fmm)
         MapHess(:,i,j)=(Fpp-Fpm-Fmp+Fmm)/(FOUR*h(i)*h(j))
         MapHess(:,j,i)=MapHess(:,i,j)
      END DO
   END DO
END SUBROUTINE Calc_Tnum_Stage_Hessian

SUBROUTINE SolveDenseSystem_GWP(Matrix,RHS,Solution,success)
   USE QDUtil_m, ONLY: Rkind,ONE
   IMPLICIT NONE
   REAL(KIND=Rkind), INTENT(IN) :: Matrix(:,:),RHS(:)
   REAL(KIND=Rkind), ALLOCATABLE, INTENT(OUT) :: Solution(:)
   LOGICAL, INTENT(OUT) :: success
   REAL(KIND=Rkind), ALLOCATABLE :: A(:,:),b(:)
   REAL(KIND=Rkind) :: factor,temp,scale
   INTEGER :: n,i,j,k,pivot

   n=SIZE(RHS)
   success=.FALSE.
   IF (SIZE(Matrix,1) /= n .OR. SIZE(Matrix,2) /= n) RETURN
   ALLOCATE(A(n,n),b(n),Solution(n))
   A=Matrix; b=RHS
   scale=MAX(ONE,MAXVAL(ABS(A)))
   DO k=1,n
      pivot=k-1+MAXLOC(ABS(A(k:n,k)),DIM=1)
      IF (ABS(A(pivot,k)) <= 100*EPSILON(ONE)*scale) RETURN
      IF (pivot /= k) THEN
         DO j=k,n
            temp=A(k,j); A(k,j)=A(pivot,j); A(pivot,j)=temp
         END DO
         temp=b(k); b(k)=b(pivot); b(pivot)=temp
      END IF
      DO i=k+1,n
         factor=A(i,k)/A(k,k)
         A(i,k:n)=A(i,k:n)-factor*A(k,k:n)
         b(i)=b(i)-factor*b(k)
      END DO
   END DO
   DO i=n,1,-1
      IF (i < n) b(i)=b(i)-DOT_PRODUCT(A(i,i+1:n),Solution(i+1:n))
      Solution(i)=b(i)/A(i,i)
   END DO
   success=.TRUE.
END SUBROUTINE SolveDenseSystem_GWP


SUBROUTINE InitQModelBasisParam_temp(nio)

   USE QDUtil_m

   IMPLICIT NONE

   ! Argument
   INTEGER, INTENT(IN) :: nio

   ! Tableaux
   REAL(KIND=Rkind), ALLOCATABLE :: GGdef(:, :)
   REAL(KIND=Rkind), ALLOCATABLE :: V(:, :)
   REAL(KIND=Rkind), ALLOCATABLE :: Hessian(:, :, :, :)
   REAL(KIND=Rkind), ALLOCATABLE :: Q0(:)
   REAL(KIND=Rkind), ALLOCATABLE :: G(:, :, :)
   REAL(KIND=Rkind), ALLOCATABLE :: ScaleQ(:)
   REAL(KIND=Rkind), ALLOCATABLE :: M(:)
   REAL(KIND=Rkind), ALLOCATABLE :: sigma(:)
   REAL(KIND=Rkind), ALLOCATABLE :: a(:)

   ! Dimensions et indices
   INTEGER :: ndim
   INTEGER :: nsurf
   INTEGER :: option
   INTEGER :: icoord
   INTEGER :: isurf

   ! Paramètres de la base
   INTEGER :: nb_basis
   INTEGER :: nq_basis

   ! Fichier
   INTEGER :: unit_file
   INTEGER :: ios
   CHARACTER(LEN=200) :: iomsg

   ! Modèle
   LOGICAL :: adiabatic

   CHARACTER(LEN=16) :: pot_name

   ! ============================================================
   ! Initialisation des paramètres du modèle
   ! ============================================================

   ! Valeurs initiales demandées
   nb_basis = 2
   nq_basis = 2

   ndim      = -1
   nsurf     = -1
   pot_name  = 'read_model'
   adiabatic = .FALSE.
   option    = 1

   IF (Debug_Enabled()) THEN
      FLUSH(out_unit)
   END IF

   ! Initialisation du modèle Q
   CALL sub_Init_Qmodel(ndim, nsurf, pot_name, adiabatic, option)

   WRITE(out_unit, *) 'ndim, nsurf = ', ndim, nsurf
   WRITE(out_unit, *) 'pot_name = ', TRIM(pot_name)

   ! ============================================================
   ! Allocation
   ! ============================================================

   ALLOCATE(GGdef(ndim, ndim))
   ALLOCATE(Hessian(nsurf, nsurf, ndim, ndim))
   ALLOCATE(G(nsurf, nsurf, ndim))
   ALLOCATE(Q0(ndim))
   ALLOCATE(M(ndim))
   ALLOCATE(ScaleQ(ndim))
   ALLOCATE(sigma(ndim))
   ALLOCATE(a(ndim))
   ALLOCATE(V(nsurf, nsurf))

   ! ============================================================
   ! Calcul des paramètres
   ! ============================================================

   CALL get_Qmodel_Q0(Q0, option)
   CALL get_Qmodel_GGdef(GGdef)
   CALL Write_VecMat(GGdef, out_unit, ndim, info='<GGdef>')

   WRITE(out_unit, *)
   CALL Write_VecMat(Q0, out_unit, ndim, info='<Q0>')
   ! Calcul du potentiel, du gradient et de la Hessienne
   CALL sub_Qmodel_VGH(V, G, Hessian, Q0)

   isurf = 1

   DO icoord = 1, ndim
      M(icoord) = ONE / GGdef(icoord, icoord)
      ScaleQ(icoord) = SQRT( &
           M(icoord) * SQRT( &
           Hessian(isurf, isurf, icoord, icoord) / M(icoord)))

      sigma(icoord) = SQRT(TWO) / ScaleQ(icoord)

      a(icoord) = ScaleQ(icoord) * ScaleQ(icoord)

   END DO

   ! ============================================================
   ! Affichage des paramètres dans out_unit
   ! ============================================================

   WRITE(out_unit, *)
   CALL Write_VecMat(M, out_unit, ndim, info='<M>')
   WRITE(out_unit, *)
   CALL Write_VecMat(ScaleQ, out_unit, ndim, info='<ScaleQ>')
   WRITE(out_unit, *)
   CALL Write_VecMat(sigma, out_unit, ndim, info='<sigma>')
   WRITE(out_unit, *)
   CALL Write_VecMat(a, out_unit, ndim, info='<Re(alpha)>')

   ! ============================================================
   ! Ouverture du fichier de sortie
   ! ============================================================

   OPEN(NEWUNIT=unit_file,FILE='qmodel_basis.txt',STATUS='REPLACE',ACTION='WRITE',FORM='FORMATTED',IOSTAT=ios,IOMSG=iomsg)

   IF (ios /= 0) THEN

      WRITE(out_unit, *) &
         'Erreur lors de l''ouverture de qmodel_basis.txt : ', &
         TRIM(iomsg)

      DEALLOCATE(GGdef)
      DEALLOCATE(Hessian)
      DEALLOCATE(Q0)
      DEALLOCATE(G)
      DEALLOCATE(M)
      DEALLOCATE(ScaleQ)
      DEALLOCATE(sigma)
      DEALLOCATE(a)
      DEALLOCATE(V)

      RETURN

   END IF

   ! ============================================================
   ! Écriture dans qmodel_basis.txt
   ! ============================================================

   DO icoord = 1, ndim

      ! Ligne de définition de la gaussienne initiale
      WRITE(unit_file, *) &
         '&defWP0 sigma=', sigma(icoord), &
         ' Beta=0.0 Qeq=', Q0(icoord), &
         ' imp_k=0.0 gamma=0.0 /'

      ! Ligne de définition de la base
      WRITE(unit_file, *) &
         "&basis_nd name='HO' nb=", nb_basis, &
         " nq=", nq_basis, &
         " Q0=", Q0(icoord), &
         " SCALEQ=", ScaleQ(icoord), &
         " Imp_k=0.0 alpha=(", a(icoord), ",0.0) /"

      ! Ligne vide entre deux coordonnées
      WRITE(unit_file, *)

   END DO

   ! Fermeture du fichier
   CLOSE(unit_file)

   WRITE(out_unit, *)
   WRITE(out_unit, *) &
      'Le fichier qmodel_basis.txt a ete cree avec succes.'

   ! ============================================================
   ! Libération de la mémoire
   ! ============================================================

   DEALLOCATE(GGdef)
   DEALLOCATE(Hessian)
   DEALLOCATE(Q0)
   DEALLOCATE(G)
   DEALLOCATE(M)
   DEALLOCATE(ScaleQ)
   DEALLOCATE(sigma)
   DEALLOCATE(a)
   DEALLOCATE(V)

   IF (Debug_Enabled()) THEN
      FLUSH(out_unit)
   END IF

END SUBROUTINE InitQModelBasisParam_temp

      SUBROUTINE H_action_diagnostic(H, psi)

      TYPE(Op_t), INTENT(IN) :: H
      TYPE(psi_t), INTENT(IN) :: psi
      TYPE(psi_t) :: Hpsi
      REAL(KIND=Rkind) :: norm_psi, norm_hpsi, energy
      INTEGER :: action_unit

      CALL init_psi(Hpsi, psi%Basis, cplx=.TRUE., grid=.FALSE.)
      CALL calc_OpPsi(H, psi, Hpsi)
      CALL Calc_Norm_OF_Psi(psi, norm_psi)
      CALL Calc_Norm_OF_Psi(Hpsi, norm_hpsi)
      energy = REAL(DOT_PRODUCT(CONJG(psi%CVec), Hpsi%CVec), KIND=Rkind)

      WRITE(out_unit,*) '===== H action diagnostic ====='
      WRITE(out_unit,*) 'norm_psi       = ', norm_psi
      WRITE(out_unit,*) 'norm_Hpsi      = ', norm_hpsi
      WRITE(out_unit,*) 'energy_psi      = ', energy
      WRITE(out_unit,*) '===== End H action diagnostic ====='
      OPEN(NEWUNIT=action_unit, FILE='H_action.dat', STATUS='REPLACE', ACTION='WRITE')
      WRITE(action_unit,'(A)') '# norm_psi norm_Hpsi energy_psi'
      WRITE(action_unit,'(3(1X,ES24.16))') norm_psi, norm_hpsi, energy
      CLOSE(action_unit)

      CALL dealloc_psi(Hpsi)

   END SUBROUTINE H_action_diagnostic 

   ! Tnum reads several namelists by rewinding the input unit.  A shell

   ! here-document may present stdin as a non-seekable pipe, so stage the

   ! complete stream in a regular file once at program start.

 SUBROUTINE pdv2psi_nD(d2psi_g, psi_g, Basis,ib)

    USE  QDUtil_m
    TYPE(Basis_t), intent(in), target               :: Basis
    complex(kind=Rkind), intent(in), target         :: psi_g(:)
    complex(kind=Rkind), intent(inout), target      :: d2psi_g(:)
    integer ,intent(in)                             :: ib
    complex(kind=Rkind), pointer                    :: psi_ggb(:, :, :)
    complex(kind=Rkind), pointer                    :: d2gg(:, :)
    complex(kind=Rkind), pointer                    :: d2psi_ggb(:, :, :)
    logical                                         :: debug
    integer                                         :: i1, i3, Ndim
    integer, allocatable                            :: Iq1(:), Iq2(:), Iq3(:), Ib1(:), Ib2(:), Ib3(:)

    debug = Basis%debug .OR. Debug_Enabled()
    IF (debug) THEN
   
       flush (out_unit)

    END IF

    Ndim = size(Basis%tab_basis)

    call Calc_index(Ib1, Ib2, Ib3, Iq1, Iq2, Iq3, Basis)
    d2psi_g(:) = CZERO

       d2psi_ggb(1:Iq1(ib), 1:Iq2(ib), 1:Iq3(ib)) => d2psi_g
       psi_ggb(1:Iq1(ib), 1:Iq2(ib), 1:Iq3(ib)) => psi_g
       d2gg(1:Iq2(ib), 1:Iq2(ib)) => Basis%tab_basis(ib)%d2gg

       DO i3 = 1, ubound(psi_ggb, dim=3)

          DO i1 = 1, ubound(psi_ggb, dim=1)

             d2psi_ggb(i1, :, i3) = d2psi_ggb(i1, :, i3) + matmul(d2gg, psi_ggb(i1, :, i3))

          END DO

       END DO

    Deallocate (Ib1, Ib2, Ib3, Iq1, Iq2, Iq3)

    IF (debug) THEN

       flush (out_unit)

    END IF

 END SUBROUTINE pdv2psi_nD

end module Op_m
