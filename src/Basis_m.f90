MODULE Basis_m
   USE QDUtil_m
   USE NDindex_m
   USE polyortho_m
   USE SmolyakSparseGridMetadata_m
   USE NUFFT_m, ONLY: NUFFT, INV_NUFFT

   IMPLICIT NONE

   PRIVATE
   LOGICAL, SAVE :: Basis_Debug_Enabled = .FALSE.
   LOGICAL, SAVE :: Propa_Debug_Enabled = .FALSE.
   PUBLIC :: Debug_Enabled, Basis_Debug_Active, Set_Propa_Debug, Reset_Debug_Controls
   PUBLIC :: Basis_t
   PUBLIC :: Read_Basis, Write_Basis
   PUBLIC :: Basis_IS_allocated, Basis_IS_allocatedtot
   PUBLIC :: Deallocate_Basis
   
   PUBLIC :: Calc_dngg_grid
   PUBLIC :: Calc_tranpose_d0gb
   PUBLIC :: test_basitogridgridtobasis
   
   PUBLIC :: GridTOBasis_nD_cplx, BasisTOGrid_nD_cplx
   
   PUBLIC :: Calc_Q_grid, Calc_index
   PUBLIC :: Rdensity_alloc, Rdensity_Writing, Rdensity_Dealloc
   
   PUBLIC :: Scale_Basis, init_Basis1_TO_Basis2
   PUBLIC :: construct_primitive_basis, construct_primitive_basis_temp
   PUBLIC :: Construct_Basis_el,Construct_Basis_Fourier,Construct_Basis_Sin
   PUBLIC :: Construct_Basis_Ho,Construct_Basis_Hag
   PUBLIC :: CheckOrtho_Basis
   
   
   PUBLIC :: Get_Basis_Parameters, Change_Basis_Parameters
   PUBLIC :: Basis_UsesSmolyakSparseGrid, Get_TotalBasisCoefficientCount
   
   PUBLIC :: FinalizeSparseGridRepresentation
   PUBLIC :: REDUCED_DENSIRY_t
   PUBLIC :: calc_tab_Iq0 ,Complete_construct_Basis
   PUBLIC :: Print_Grids, Build_NUFT_Grids, Test_NUFFT_Cycle
   

!----------------------------------------------------------------------
! TYPE: Basis_t
!
! DESCRIPTION:
!   Data structure defining a numerical basis set for quantum dynamics
!   calculations. It contains the number of basis functions, scaling
!   parameters, coordinate grids, weights, and precomputed operator
!   matrices (position, momentum, derivatives, etc.).
!
! FIELDS:
!   nb_basis   : Total number of basis functions.
!   nb         : Number of basis functions along one coordinate.
!   nq         : Number of grid points.
!   A, B       : Scaling or transformation parameters (user-defined meaning).
!   imp_k      : momentum (p₀ or k₀ type quantity).
!   scaleQ     : Scaling factor for coordinate Q.
!   Q0         : Coordinate origin shift.
!   Alpha      : Complex parameter (often Gaussian width or similar).
!   Basis_name : Name of the basis set (dynamically allocated string).
!   x, w       : Grid points and quadrature weights.
!   d0*, d1*, d2* : Precomputed derivative/operator matrices in various forms.
!   dag*, dq0*, dp0* : Precomputed derivative/shift operators.
!   S          : Overlap matrix of the basis functions.
!   NDindexq   : Index mapping for grid representation.
!   NDindexb   : Index mapping for basis representation.
!   tab_basis  : Table of sub-basis structures.
!----------------------------------------------------------------------

   TYPE :: Basis_t
   INTEGER                             :: nb_basis = ZERO
   LOGICAL                             :: debug = .FALSE.
   LOGICAL                             :: use_initial_gwp = .FALSE.
   INTEGER                             :: nb, nq
   REAL(KIND=Rkind)                    :: A, B
   REAL(KIND=Rkind)                    :: imp_k    ! Momentum parameter
   REAL(KIND=Rkind)                    :: scaleQ, Q0
   COMPLEX(KIND=Rkind)                 :: Alpha = CONE
   CHARACTER(LEN=:), ALLOCATABLE       :: Basis_name
   INTEGER                             :: representation_layout = BASIS_LAYOUT_DIRECT_PRODUCT
   ! Smolyak levels: LB controls basis blocks, LG controls grid blocks.
   INTEGER                             :: LB = 0
   INTEGER                             :: LG = 0
   INTEGER                             :: smolyak_rule_shift_B = 0
   INTEGER                             :: smolyak_rule_slope_B = 0
   INTEGER                             :: smolyak_rule_shift_G = 0
   INTEGER                             :: smolyak_rule_slope_G = 0
   TYPE(SmolyakSparseGridDescriptor_t) :: SmolyakSparseGrid
   !<--- Grid points and weights
   REAL(KIND=Rkind),    ALLOCATABLE    :: x(:)      ! Grid points
   REAL(KIND=Rkind),    ALLOCATABLE    :: w(:)      ! Quadrature weights
   !<--- Derivative matrices
   COMPLEX(KIND=Rkind), ALLOCATABLE    :: d0gb(:, :)
   COMPLEX(KIND=Rkind), ALLOCATABLE    :: d0bgw(:, :)
   COMPLEX(KIND=Rkind), ALLOCATABLE    :: d1gb(:, :, :)
   COMPLEX(KIND=Rkind), ALLOCATABLE    :: d1bb(:, :, :)
   COMPLEX(KIND=Rkind), ALLOCATABLE    :: d1gg(:, :, :)
   COMPLEX(KIND=Rkind), ALLOCATABLE    :: d2gb(:, :, :, :)
   COMPLEX(KIND=Rkind), ALLOCATABLE    :: d2gg(:, :, :, :)
   !<--- Overlap matrix
   COMPLEX(KIND=Rkind), ALLOCATABLE    :: S(:, :)   ! Overlap matrix
   !<--- Index mapping
   TYPE(NDindex_t)                     :: NDindexq
   TYPE(NDindex_t)                     :: NDindexb
   !<--- Sub-basis
   TYPE(Basis_t), ALLOCATABLE          :: tab_basis(:)     
END TYPE Basis_t

   TYPE REDUCED_DENSIRY_t
      real(kind=Rkind),    allocatable    :: Norm(:)
      real(kind=Rkind),    allocatable    :: prob(:)
      TYPE(REDUCED_DENSIRY_t),allocatable :: tab_prob(:)
   END TYPE REDUCED_DENSIRY_t

CONTAINS

   LOGICAL FUNCTION Debug_Enabled()
      Debug_Enabled = Propa_Debug_Enabled
   END FUNCTION Debug_Enabled

   LOGICAL FUNCTION Basis_Debug_Active()
      Basis_Debug_Active = Basis_Debug_Enabled
   END FUNCTION Basis_Debug_Active

   SUBROUTINE Set_Propa_Debug(enabled)
      LOGICAL, INTENT(IN) :: enabled
      Propa_Debug_Enabled = enabled
   END SUBROUTINE Set_Propa_Debug

   SUBROUTINE Reset_Debug_Controls()
      Basis_Debug_Enabled = .FALSE.
      Propa_Debug_Enabled = .FALSE.
   END SUBROUTINE Reset_Debug_Controls

   RECURSIVE FUNCTION Basis_IS_allocated(Basis) RESULT(alloc)

      TYPE(Basis_t), intent(in)  :: Basis
      logical                      :: alloc
      integer                      :: i

      IF (Basis%Basis_name == 'el' .AND. Basis%nb > 0) Then
         alloc = .TRUE.
         RETURN
      END IF

      alloc = allocated(Basis%tab_basis)
      IF (allocated(Basis%tab_basis)) THEN
         Do i = 1, size(Basis%tab_basis)
            alloc = alloc .and. Basis_IS_allocated(Basis%tab_basis(i))
         END DO
      ELSE
         alloc = allocated(Basis%x)
         alloc = alloc .AND. allocated(Basis%w)
         alloc = alloc .AND. allocated(Basis%d0gb)
         alloc = alloc .AND. allocated(Basis%d1gb)
         alloc = alloc .AND. allocated(Basis%d2gb)
         !alloc = alloc .AND. allocated(Basis%d0bgw)
      END IF
   END FUNCTION Basis_IS_allocated

   RECURSIVE FUNCTION Basis_IS_allocatedtot(Basis) RESULT(alloc)

      TYPE(Basis_t), intent(in)    :: Basis
      logical                      :: alloc
      integer                      :: i

      alloc = allocated(Basis%tab_basis)
      IF (allocated(Basis%tab_basis)) THEN
         Do i = 1, size(Basis%tab_basis)
            alloc = alloc .and. Basis_IS_allocated(Basis%tab_basis(i))
         END DO
      ELSE
         alloc = allocated(Basis%x)
         alloc = alloc .AND. allocated(Basis%w)
         alloc = alloc .AND. allocated(Basis%d0gb)
         alloc = alloc .AND. allocated(Basis%d1gb)
         alloc = alloc .AND. allocated(Basis%d2gb)
         alloc = alloc .AND. allocated(Basis%d1gg)
         alloc = alloc .AND. allocated(Basis%d2gg)

      END IF

   END FUNCTION Basis_IS_allocatedtot

   RECURSIVE SUBROUTINE Write_Basis(Basis)
      USE QDUtil_m

      TYPE(Basis_t), intent(in)  :: Basis

      !<--- Local variables
      INTEGER                         :: i
      LOGICAL :: debug

      IF(Basis%Basis_name == 'el')THEN
         WRITE(out_unit,*) 'Write_Basis: Basis_name=el'
         WRITE(out_unit,*) 'We cannot write the el basis'
       RETURN
       FLUSH(out_unit)
      END IF

      debug = Basis%debug .OR. Propa_Debug_Enabled

      IF (debug) THEN
         WRITE(out_unit,*) '---------------------------------------------------------------------'
         WRITE(out_unit,*) 'Write_Basis'
         !WRITE(out_unit,*) 'Write_Basis'
         WRITE(out_unit,*) "Basis_name=", Basis%Basis_name
         WRITE(out_unit,*) "n_basis=", Basis%nb_basis
         WRITE(out_unit,*) 'nb,nq',Basis%nb,Basis%nq
         WRITE(out_unit,*) "Q0,scaleQ ", Basis%Q0,Basis%scaleQ
         WRITE(out_unit,*)  'A,B',Basis%A,Basis%B
         WRITE(out_unit,*) 'nb_basis',Basis%nb_basis
         WRITE(out_unit,*) 'imp_k',Basis%imp_k
         WRITE(out_unit,*) 'Alpha',Basis%Alpha
         FLUSH(out_unit)
      END IF 

      IF (.NOT. ALLOCATED(Basis%x)) THEN
         WRITE(out_unit,*)' Basis table x is not allocated.'
      ELSE
         CALL  Write_VecMat(Basis%x, out_unit, 5,  info='Grid x')
      END IF
      WRITE (out_unit, *)
      IF (.NOT. ALLOCATED(Basis%W)) THEN
           WRITE(out_unit,*)' Basis table w is not allocated.'
      ELSE
         CALL  Write_VecMat(Basis%w, out_unit, 5,  info='Weights')
      END IF
      WRITE (out_unit, *)
       IF (.NOT.ALLOCATED(Basis%d0gb)) THEN
        WRITE(out_unit,*)' Basis table d0gb is not allocated.'
      ELSE
         CALL   Write_VecMat(Basis%d0gb, out_unit, 5,  info='d0gb')
       END IF
      WRITE (out_unit, *)
        IF (.NOT.ALLOCATED(Basis%d0bgw)) THEN
           WRITE(out_unit,*)' Basis table d0bgw is not allocated.'
         ELSE
      CALL   Write_VecMat(Basis%d0bgw,out_unit,5, info='d0gbw')
        END IF
      WRITE (out_unit, *)
         IF (.NOT.ALLOCATED(Basis%d1gb)) THEN
           WRITE(out_unit,*)' Basis table d1gb is not allocated.'
         ELSE
      CALL   Write_VecMat(Basis%d1gb(:,:,1),out_unit,5, info='d1gb')
         END IF
      WRITE (out_unit, *)
         IF (.NOT.ALLOCATED(Basis%d1gg)) THEN
           WRITE(out_unit,*)' Basis table d1gb is not allocated.'
         ELSE
      CALL   Write_VecMat(Basis%d1gg(:,:,1),out_unit,5, info='d1gg')
         END IF
      WRITE (out_unit, *)
         IF (.NOT.ALLOCATED(Basis%d2gb)) THEN
           WRITE(out_unit,*)' Basis table d1gb is not allocated.'
         ELSE
      CALL   Write_VecMat(Basis%d2gb(:,:,1,1),out_unit,5, info='d2gb')
         END IF

      WRITE (out_unit, *)
         IF (.NOT.ALLOCATED(Basis%d2gg)) THEN
           WRITE(out_unit,*)' Basis table d2gg is not allocated.'
         ELSE
      CALL   Write_VecMat(Basis%d2gg(:,:,1,1),out_unit,5, info='d2gg')
         END IF
          
      IF(debug) THEN
      WRITE (out_unit, *)
         IF (.NOT.ALLOCATED(Basis%S)) THEN
            WRITE(out_unit,*)' Basis table S is not allocated.'
         ELSE
            CALL  Write_VecMat(Basis%S,out_unit,5, info='S')
         END IF
      END IF

      IF (ALLOCATED(Basis%tab_basis)) THEN
         DO i = 1, size(Basis%tab_basis)
            CALL Write_Basis(Basis%tab_basis(i))
         END DO
      END IF

      IF(debug) THEN
         WRITE (out_unit, *) '--------------------------------------------------------------------------'
         WRITE (out_unit, *) 'END Write_Basis'
         FLUSH(out_unit)
      END IF

   END SUBROUTINE Write_Basis


   !> Recursive subroutine to initialize Basis2 from Basis1
   !! It supports both composite (multi-dimensional) and 1D primitive basis structures.
   RECURSIVE SUBROUTINE init_Basis1_TO_Basis2(Basis2, Basis1)
   USE QDUtil_m
   IMPLICIT NONE
     
     TYPE(Basis_t), intent(in)           :: Basis1   ! Input basis structure
     TYPE(Basis_t), intent(inout)        :: Basis2   ! Output (initialized) basis structure
     integer                             :: ib        ! Loop index
     
     IF (allocated(Basis1%tab_basis)) THEN
        ! Composite basis: recursively initialize each sub-basis
        call Deallocate_Basis(Basis2)
        Basis2%Basis_name  = Basis1%Basis_name
        Basis2%nb_basis    = Basis1%nb_basis
        Basis2%debug       = Basis1%debug
        Basis2%use_initial_gwp = Basis1%use_initial_gwp
        Basis2%representation_layout = Basis1%representation_layout
        Basis2%LB = Basis1%LB
        Basis2%LG = Basis1%LG
         Basis2%smolyak_rule_shift_B = Basis1%smolyak_rule_shift_B
         Basis2%smolyak_rule_slope_B = Basis1%smolyak_rule_slope_B
         Basis2%smolyak_rule_shift_G = Basis1%smolyak_rule_shift_G
         Basis2%smolyak_rule_slope_G = Basis1%smolyak_rule_slope_G
         allocate(Basis2%tab_basis(Basis2%nb_basis))
     
        DO ib = 1, Basis1%nb_basis
           call init_Basis1_TO_Basis2(Basis2%tab_basis(ib), Basis1%tab_basis(ib))
        END DO
     
        ! Recalculate total number of basis functions and quadrature points
        Basis2%nb = 1
        Basis2%nq = 1
        DO ib = 1, Basis1%nb_basis
           IF (Basis2%tab_basis(ib)%Basis_name == 'el') CYCLE
           Basis2%nb = Basis2%nb * Basis2%tab_basis(ib)%nb
           Basis2%nq = Basis2%nq * Basis2%tab_basis(ib)%nq
        END DO
     
     ELSE
        ! Elementary basis: directly copy parameters
        Basis2%Basis_name  = Basis1%Basis_name
        Basis2%nb_basis    = Basis1%nb_basis
        Basis2%debug       = Basis1%debug
        Basis2%use_initial_gwp = Basis1%use_initial_gwp
        Basis2%nb          = Basis1%nb
        Basis2%nq          = Basis1%nq
        Basis2%Q0          = Basis1%Q0
        Basis2%SCALEQ      = Basis1%SCALEQ
        Basis2%Imp_k       = Basis1%Imp_k
        Basis2%Alpha       = Basis1%Alpha
        Basis2%A           = Basis1%A
        Basis2%B           = Basis1%B
        Basis2%S           = Basis1%S
        Basis2%representation_layout = Basis1%representation_layout
        Basis2%LB = Basis1%LB
        Basis2%LG = Basis1%LG
         Basis2%smolyak_rule_shift_B = Basis1%smolyak_rule_shift_B
         Basis2%smolyak_rule_slope_B = Basis1%smolyak_rule_slope_B
         Basis2%smolyak_rule_shift_G = Basis1%smolyak_rule_shift_G
         Basis2%smolyak_rule_slope_G = Basis1%smolyak_rule_slope_G
      END IF
     
     END SUBROUTINE init_Basis1_TO_Basis2

   !> Recursive subroutine to safely deallocate all allocated components
   !! of a Basis_t object, including nested basis structures.
   RECURSIVE SUBROUTINE Deallocate_Basis(Basis)
   USE QDUtil_m
   IMPLICIT NONE
   
   TYPE(Basis_t), intent(inout) :: Basis

   !<--- Local variables
   LOGICAL :: debug
   INTEGER                      :: i

   debug = Basis%debug .OR. Propa_Debug_Enabled

   IF(debug) THEN
      WRITE(out_unit,*) 'BEGINNING Deallocate_Basis'
      FLUSH(out_unit)
   END IF
   
   CALL DestroySmolyakSparseGridDescriptor(Basis%SmolyakSparseGrid)
   Basis%representation_layout = BASIS_LAYOUT_DIRECT_PRODUCT
   Basis%LB = 0
   Basis%LG = 0
   Basis%smolyak_rule_shift_B = 0
   Basis%smolyak_rule_slope_B = 0
   Basis%smolyak_rule_shift_G = 0
   Basis%smolyak_rule_slope_G = 0

   IF(Basis%Basis_name == 'el') THEN
      IF(debug) WRITE(out_unit,*) 'Basis_name=el'
      IF(debug) WRITE(out_unit,*) 'Elementary basis: nothing to deallocate.'
      RETURN
   END IF
   
   ! Deallocate composite basis recursively before releasing the container
   IF (ALLOCATED(Basis%tab_basis)) THEN
      DO i = 1, SIZE(Basis%tab_basis)
         CALL Deallocate_Basis(Basis%tab_basis(i))
      END DO
      DEALLOCATE(Basis%tab_basis)
   END IF

   ! Deallocate ND indexing arrays if present
   IF (ALLOCATED(Basis%NDindexq%Tab0)) DEALLOCATE(Basis%NDindexq%Tab0)
   IF (ALLOCATED(Basis%NDindexb%Tab0)) DEALLOCATE(Basis%NDindexb%Tab0)
   
   ! Deallocate 1D basis arrays if present
   IF (ALLOCATED(Basis%x)) THEN
      DEALLOCATE(Basis%x)
      IF(debug) WRITE(out_unit,*) 'Basis table x is now deallocated.'
   END IF
   
   IF (ALLOCATED(Basis%d0gb)) THEN
      DEALLOCATE(Basis%d0gb)
      IF(debug) WRITE(out_unit,*) 'Basis table d0gb is now deallocated.'
   END IF
   
   IF (ALLOCATED(Basis%d0bgw)) THEN
      DEALLOCATE(Basis%d0bgw)
      IF(debug) WRITE(out_unit,*) 'Basis table d0bgw is now deallocated.'
   END IF
      
   IF (ALLOCATED(Basis%d1gb)) THEN
      DEALLOCATE(Basis%d1gb)
      IF(debug) WRITE(out_unit,*) 'Basis table d1gb is now deallocated.'
   END IF
   
   IF (ALLOCATED(Basis%d1gg)) THEN
      DEALLOCATE(Basis%d1gg)
      IF(debug) WRITE(out_unit,*) 'Basis table d1gg is now deallocated.'
   END IF
   
   IF (ALLOCATED(Basis%d2gb)) THEN
      DEALLOCATE(Basis%d2gb)
      IF(debug) WRITE(out_unit,*) 'Basis table d2gb is now deallocated.'
   END IF
   
   IF (ALLOCATED(Basis%d2gg)) THEN
      DEALLOCATE(Basis%d2gg)
      IF(debug) WRITE(out_unit,*) 'Basis table d2gg is now deallocated.'
   END IF
   

   IF(debug) THEN
      WRITE(out_unit,*) 'END Deallocate_Basis'
      FLUSH(out_unit)
   END IF
   
   END SUBROUTINE Deallocate_Basis

!=======================================================================
!  RECURSIVE SUBROUTINE: Read_Basis
!
!  PURPOSE:
!    Reads a basis set definition from an input unit (NIO) using the 
!    namelist "basis_nD". Supports hierarchical basis definitions 
!    through recursion when `nb_basis > 1`.
!
!  INPUT:
!    nio   : Logical unit number from which the basis namelist is read.
!
!  INPUT/OUTPUT:
!    Basis : Basis_t structure to be filled with parameters and 
!            optional sub-bases.
!
!  NAMELIST FIELDS:
!    name       : Name of the basis (string).
!    nb_basis   : Number of basis components (for composite/multidim bases).
!    nb         : Number of basis functions along one coordinate.
!    nq         : Number of grid points.
!    A, B       : Basis scaling or transformation parameters.
!    scaleQ     : Scaling factor for the coordinate Q.
!    Q0         : Coordinate origin shift.
!    Imp_k      : momentum (e.g., k₀ in Gaussian packets).
!    Alpha      : Complex basis parameter (width, phase, etc.).
!
!  NOTES:
!    - If `nb_basis > 1`, the routine allocates `tab_basis` and calls 
!      itself recursively for each sub-basis.
!    - If `nb_basis <= 1`, it initializes the overlap matrix S to identity.
!=======================================================================
   RECURSIVE SUBROUTINE Read_Basis(Basis, nio)
   USE QDUtil_m
   LOGICAL                              :: debug
   TYPE(Basis_t), INTENT(INOUT)        :: Basis
   INTEGER, INTENT(IN)                 :: nio

   ! Local variables
   INTEGER                              :: err_io, nb, nq, i, nb_basis, ib
   CHARACTER(LEN=Name_len)              :: name
   INTEGER                              :: LB, LG, A_b_smol, B_b_smol, A_g_smol, B_g_smol
   LOGICAL                              :: use_initial_gwp
   LOGICAL                              :: q0_is_set, impk_is_set, alpha_is_set
   CHARACTER(LEN=Name_len)              :: normalized_name
   REAL(KIND=Rkind)                     :: A, B, scaleQ, Q0, Imp_k
   COMPLEX(KIND=Rkind)                  :: Alpha

   ! Namelist declaration
   NAMELIST /basis_nD/ name, nb_basis, nb, nq, A, B, scaleQ, Q0, Imp_k, Alpha, LB, LG, A_b_smol, B_b_smol, A_g_smol, B_g_smol, debug, use_initial_gwp

   ! ===== Default values =====
   nb_basis = 0
   nb       = 0
   nq       = 0
   A        = ZERO
   B        = ZERO
   Q0       = HUGE(ONE)
   scaleQ   = ONE
   name     = '0'
   Imp_k    = HUGE(ONE)
   Alpha    = CMPLX(-HUGE(ONE), -HUGE(ONE), KIND=Rkind)
   LB       = 0
   LG       = 0
   A_b_smol = 0
   B_b_smol = 0
   A_g_smol = 0
   B_g_smol = 0
   debug = .FALSE.
   use_initial_gwp = .FALSE.

   ! ===== Read namelist =====
   READ (nio, NML=basis_nD, IOSTAT=err_io)
   Basis_Debug_Enabled = Basis_Debug_Enabled .OR. debug
   q0_is_set = Q0 /= HUGE(ONE)
   impk_is_set = Imp_k /= HUGE(ONE)
   alpha_is_set = REAL(Alpha, KIND=Rkind) /= -HUGE(ONE)
   IF (.NOT. q0_is_set) Q0 = ZERO
   IF (.NOT. impk_is_set) Imp_k = ZERO
   IF (.NOT. alpha_is_set) Alpha = CONE
   normalized_name = TRIM(ADJUSTL(name))
   CALL string_uppercase_TO_lowercase(normalized_name)
   IF (TRIM(normalized_name) == 'ho' .OR. TRIM(normalized_name) == 'hag') THEN
      ! When any essential Gaussian geometry field is omitted, use the unique
      ! initial-GWP definition instead of silently constructing a mismatched basis.
      use_initial_gwp = use_initial_gwp .OR. (.NOT. q0_is_set .OR. .NOT. impk_is_set .OR. .NOT. alpha_is_set)
   END IF

   ! ===== Debug / formatted output =====
   IF (debug) THEN
      WRITE(out_unit,'(A20,": ",A)')        'Basis name', TRIM(name)
      WRITE(out_unit,'(A20,": ",I5)')       'nb_basis',   nb_basis
      WRITE(out_unit,'(A20,": ",I5)')       'nb',         nb
      WRITE(out_unit,'(A20,": ",I5)')       'nq',         nq
      WRITE(out_unit,'(A20,": ",ES30.15)')  'A',          A
      WRITE(out_unit,'(A20,": ",ES30.15)')  'B',          B
      WRITE(out_unit,'(A20,": ",ES30.15)')  'scaleQ',     scaleQ
      WRITE(out_unit,'(A20,": ",ES30.15)')  'Q0',         Q0
      WRITE(out_unit,'(A20,": ",ES30.15)')  'Imp_k',      Imp_k
      WRITE(out_unit,'(A20,": ",ES30.15)')  'Alpha (real)', REAL(Alpha)
      WRITE(out_unit,'(A20,": ",ES30.15)')  'Alpha (imag)', AIMAG(Alpha)
      WRITE(out_unit,'(A20,": ",I5)')       'LB',         LB
      WRITE(out_unit,'(A20,": ",I5)')       'LG',         LG
      WRITE(out_unit,'(A20,": ",I5)')       'A_b_smol',   A_b_smol
      WRITE(out_unit,'(A20,": ",I5)')       'B_b_smol',   B_b_smol
      WRITE(out_unit,'(A20,": ",I5)')       'A_g_smol',   A_g_smol
      WRITE(out_unit,'(A20,": ",I5)')       'B_g_smol',   B_g_smol
      WRITE(out_unit,*)  ! empty line for spacing
   END IF

   ! ===== Check I/O errors =====
   IF (err_io < 0) THEN
      WRITE (out_unit, *) ' ERROR in Read_Basis: end of file or record.'
      STOP ' ERROR in Read_Basis: namelist reading failed.'
   ELSE IF (err_io > 0) THEN
      WRITE (out_unit, *) ' ERROR in Read_Basis: namelist arguments may be wrong.'
      STOP ' ERROR in Read_Basis: invalid namelist arguments.'
   END IF

   ! ===== Composite basis case =====
   normalized_name = TRIM(ADJUSTL(name))
   CALL string_uppercase_TO_lowercase(normalized_name)

   IF (nb_basis > 1) THEN
      Basis%Basis_name = TRIM(normalized_name)
      Basis%nb_basis   = nb_basis
      Basis%LB = LB
      Basis%LG = LG
      Basis%debug = debug
      Basis%use_initial_gwp = use_initial_gwp
      Basis%A = A
      Basis%B = B
      ! Backward-compatible input: when LG is omitted, use LB for the grid.
      IF (Basis%LG <= 0) Basis%LG = Basis%LB
      SELECT CASE (TRIM(normalized_name))
      CASE ('smolyak', 'smoly')
         Basis%representation_layout = BASIS_LAYOUT_SMOLYAK_SPARSE
      CASE DEFAULT
         Basis%representation_layout = BASIS_LAYOUT_DIRECT_PRODUCT
         Basis%Basis_name = 'dp'
      END SELECT

      ALLOCATE (Basis%tab_basis(nb_basis))
      DO i = 1, nb_basis
         CALL Read_Basis(Basis%tab_basis(i), nio)
         IF (debug) CALL Set_Basis_Debug(Basis%tab_basis(i), .TRUE.)

         ! The Smolyak controls belong to the composite definition, but the
         ! cap rule is evaluated independently in every nuclear dimension.
         ! Therefore propagate composite values to a primitive only when that
         ! primitive did not provide an explicit value of its own.
         IF (TRIM(normalized_name) == 'smolyak' .AND. &
             Basis%tab_basis(i)%Basis_name /= 'el') THEN
            IF (Basis%tab_basis(i)%smolyak_rule_shift_B <= 0) &
                 Basis%tab_basis(i)%smolyak_rule_shift_B = A_b_smol
            IF (Basis%tab_basis(i)%smolyak_rule_slope_B <= 0) &
                 Basis%tab_basis(i)%smolyak_rule_slope_B = B_b_smol
            IF (Basis%tab_basis(i)%smolyak_rule_shift_G <= 0) &
                 Basis%tab_basis(i)%smolyak_rule_shift_G = A_g_smol
            IF (Basis%tab_basis(i)%smolyak_rule_slope_G <= 0) &
                 Basis%tab_basis(i)%smolyak_rule_slope_G = B_g_smol
         END IF
      END DO
      DO i = 1, nb_basis
         IF (Basis%tab_basis(i)%debug) Basis%debug = .TRUE.
      END DO

      ! Keep the composite values available as well, notably for copies of
      ! the basis made by the propagation and validation code.
      Basis%smolyak_rule_shift_B = A_b_smol
      Basis%smolyak_rule_slope_B = B_b_smol
      Basis%smolyak_rule_shift_G = A_g_smol
      Basis%smolyak_rule_slope_G = B_g_smol

      ! Compute total nuclear nb and nq excluding electronic basis blocks.
      Basis%nb = 1
      Basis%nq = 1
      DO i = 1, nb_basis
         IF (Basis%tab_basis(i)%Basis_name == 'el') CYCLE
         Basis%nb = Basis%nb * Basis%tab_basis(i)%nb
         Basis%nq = Basis%nq * Basis%tab_basis(i)%nq
      END DO

   ! ===== Single basis case =====
   ELSE
      Basis%nb_basis   = nb_basis
      Basis%debug      = debug
      Basis%use_initial_gwp = use_initial_gwp
      Basis%nb         = nb
      Basis%nq         = nq
      Basis%Q0         = Q0
      Basis%scaleQ     = scaleQ
      Basis%A          = A
      Basis%B          = B
      Basis%Imp_k      = Imp_k
      Basis%Alpha      = Alpha
      Basis%Basis_name = TRIM(ADJUSTL(name))
      Basis%representation_layout = BASIS_LAYOUT_DIRECT_PRODUCT
      Basis%smolyak_rule_shift_B = A_b_smol
      Basis%smolyak_rule_slope_B = B_b_smol
      Basis%smolyak_rule_shift_G = A_g_smol
      Basis%smolyak_rule_slope_G = B_g_smol

      ! Initialize overlap matrix S as identity
      ALLOCATE (Basis%S(nb, nb))
      Basis%S(:, :) = ZERO
      DO ib = 1, Basis%nb
         Basis%S(ib, ib) = CONE
      END DO

      CALL string_uppercase_TO_lowercase(Basis%Basis_name)
   END IF

END SUBROUTINE Read_Basis

RECURSIVE SUBROUTINE Set_Basis_Debug(Basis, enabled)
   TYPE(Basis_t), INTENT(INOUT) :: Basis
   LOGICAL, INTENT(IN) :: enabled
   INTEGER :: ib
   Basis%debug = enabled
   IF (ALLOCATED(Basis%tab_basis)) THEN
      DO ib = 1, SIZE(Basis%tab_basis)
         CALL Set_Basis_Debug(Basis%tab_basis(ib), enabled)
      END DO
   END IF
END SUBROUTINE Set_Basis_Debug

   
   !> Calculates the index table Tab_iq for multi-dimensional quadrature indexing
   !! Tab_iq(:, iq) contains the index vector for the iq-th point in the total basis grid.
   SUBROUTINE Calc_tab_Iq0(Tab_Iq, Basis)
      USE QDUtil_m
      ! Purpose: Calculate the table of quantum indices (Tab_Iq) for a given basis
      
      ! Input/Output variables:
      TYPE(Basis_t), INTENT(IN), TARGET :: Basis       ! Basis set information
      INTEGER, allocatable, INTENT(INOUT) :: Tab_iq(:,:)  ! Output table of quantum indices
      
      ! Local variables:
      INTEGER, ALLOCATABLE :: Tab_iq0(:)    ! Temporary array for current quantum indices
      INTEGER :: ndim                      ! Number of dimensions
      INTEGER :: iq                        ! Current quantum state index
      INTEGER :: nq                        ! Total number of quantum states
      LOGICAL :: Endloop                   ! Flag to indicate end of loop
      
      ! Get dimensions from the basis set:
      ndim = size(Basis%tab_basis) - 1     ! Number of dimensions (subtracting 1 for 0-based or other reason)
      nq = Basis%nq                        ! Total number of quantum states
      
      ! Allocate memory for the output table and temporary array:
      ALLOCATE(Tab_iq(ndim, nq), Tab_iq0(ndim))
      
      ! Initialize the temporary index array:
      CALL Init_tab_ind(Tab_iq0, Basis%NDindexq)
      
      ! Loop through all quantum states:
      iq = 0
      DO
          iq = iq + 1  ! Increment quantum state counter
          
          ! Get the next set of quantum numbers:
          CALL increase_NDindex(Tab_Iq0, Basis%NDindexq, Endloop)
          
          ! Exit if we've processed all quantum states:
          IF (Endloop) EXIT
          
          ! Store the current quantum numbers in the output table:
          Tab_iq(:, iq) = Tab_Iq0
      END DO
      
      ! Clean up temporary array:
      DEALLOCATE(Tab_Iq0)
  END SUBROUTINE Calc_tab_Iq0

   RECURSIVE SUBROUTINE construct_primitive_basis_temp(Basis)
   USE QDUtil_m
   IMPLICIT NONE

   !==================== Parameters =====================
   LOGICAL :: debug

   !==================== Interface ======================
   TYPE(Basis_t), INTENT(INOUT) :: Basis   ! Basis structure to construct

   !==================== Local variables ================
   INTEGER :: ndim     ! Number of sub-bases
   INTEGER :: ib       ! Index for loop over sub-bases

   !==================== Logic ==========================
   debug = Basis%debug .OR. Propa_Debug_Enabled
   ndim = Basis%nb_basis

   ! If the basis is composite (i.e., tab_basis is allocated),
   ! recursively construct the primitive basis of each subcomponent
   IF (ALLOCATED(Basis%tab_basis)) THEN
      DO ib = 1, ndim
         CALL construct_primitive_basis_temp(Basis%tab_basis(ib))
      END DO

   ELSE
      ! Otherwise, construct the basis based on its name/type
      SELECT CASE (Basis%Basis_name)

      CASE ('el')
         ! Electronic basis: no spatial degree of freedom
         WRITE(6, *) 'Electronic basis. Electronic state number:', Basis%nb
         Basis%nq = 0

      CASE ('boxab')
         ! Particle-in-a-box basis between [A, B]
         CALL Construct_Basis_Sin(Basis)
         Basis%Q0     = Basis%A
         Basis%scaleQ = pi / (Basis%B - Basis%A)
         CALL Complete_construct_Basis(Basis)

      CASE ('fourier')
         ! Plane wave / Fourier basis (unbounded domain)
         Basis%Q0     = ZERO
         Basis%scaleQ = ONE
         CALL Construct_Basis_Fourier(Basis)
         CALL Complete_construct_Basis(Basis)

      CASE ('hag', 'ho')
         ! Harmonic or Hagedorn basis (semiclassical Gaussian functions)
         !CALL Construct_Basis_Ho(Basis)  ! Uncomment if needed
         CALL Construct_Basis_Hag(Basis)
         CALL Complete_construct_Basis(Basis)

      CASE DEFAULT
         ! Unknown or unsupported basis type
         STOP 'ERROR  Nothing to construct'

      END SELECT
   END IF

   ! Optional debug output
   !WRITE(out_unit,*) 'End construct primitive Basis'

END SUBROUTINE construct_primitive_basis_temp

   RECURSIVE SUBROUTINE construct_primitive_basis(Basis)
      USE QDUtil_m
      TYPE(Basis_t), INTENT(INOUT)           :: Basis

      !<==================== Local variables ====================>
      LOGICAL                                :: debug
      INTEGER, ALLOCATABLE                   :: NDend_q(:)
      INTEGER, ALLOCATABLE                   :: NDend_b(:)
      INTEGER                                :: ndim, ib

       debug = Basis%debug .OR. Propa_Debug_Enabled
       ndim = Basis%nb_basis - 1
       ALLOCATE (NDend_q(ndim))
       ALLOCATE (NDend_b(ndim)) 

       IF (debug) THEN
          write (out_unit, *) 'construct_primitive_basis:-------------------------------------------------'
          write (out_unit, *) 'Basis%Basis_name= ', Basis%Basis_name
       END IF

       DO ib = 1, ndim
          NDend_q(ib) = Basis%tab_basis(ib)%nq
          NDend_b(ib) = Basis%tab_basis(ib)%nb
       END DO 

       call Init_NDindex(Basis%NDindexq, NDend_q, ndim)
       call Init_NDindex(Basis%NDindexb, NDend_b, ndim)

       call construct_primitive_basis_temp(Basis)
       call FinalizeSparseGridRepresentation(Basis)

       IF (debug) THEN
          write (out_unit, *) 'End construct_primitive_basis'
       END IF

        
   END SUBROUTINE 


   LOGICAL FUNCTION Basis_UsesSmolyakSparseGrid(Basis)

      TYPE(Basis_t), INTENT(IN) :: Basis

      Basis_UsesSmolyakSparseGrid = Basis%representation_layout == BASIS_LAYOUT_SMOLYAK_SPARSE .AND. Basis%SmolyakSparseGrid%is_initialized

   END FUNCTION Basis_UsesSmolyakSparseGrid



   INTEGER FUNCTION Get_TotalBasisCoefficientCount(Basis) RESULT(ncoeff)

      TYPE(Basis_t), INTENT(IN) :: Basis

      INTEGER :: ndim, nsurf

      IF (Basis_UsesSmolyakSparseGrid(Basis)) THEN

         ncoeff = GetEffectiveTotalBasisSize(Basis%SmolyakSparseGrid)

      ELSE

         IF (ALLOCATED(Basis%tab_basis)) THEN

            ndim = SIZE(Basis%tab_basis) - 1

            IF (ndim >= 1) THEN

               nsurf = Basis%tab_basis(ndim+1)%nb

            ELSE

               nsurf = 1

           END IF

         ELSE

            nsurf = 1

         END IF

         ncoeff = Basis%nb * nsurf

      END IF

   END FUNCTION Get_TotalBasisCoefficientCount



   SUBROUTINE FinalizeSparseGridRepresentation(Basis)

      TYPE(Basis_t), INTENT(INOUT) :: Basis

      INTEGER :: ndim, nsurf, ib

      INTEGER, ALLOCATABLE :: rule_shift(:), rule_slope(:), rule_shift_g(:), rule_slope_g(:), nuclear_nb(:), nuclear_nq(:), ndend_q(:), ndend_b(:)

      TYPE(NDindex_t) :: NuclearGridIndex, NuclearBasisIndex



      IF (.NOT. ALLOCATED(Basis%tab_basis)) RETURN

      IF (Basis%representation_layout /= BASIS_LAYOUT_SMOLYAK_SPARSE) RETURN



      ndim = SIZE(Basis%tab_basis) - 1

      IF (ndim < 1) RETURN

      nsurf = Basis%tab_basis(ndim+1)%nb


      ALLOCATE(rule_shift(ndim), rule_slope(ndim), rule_shift_g(ndim), rule_slope_g(ndim), &
               nuclear_nb(ndim), nuclear_nq(ndim), ndend_q(ndim), ndend_b(ndim))

      DO ib = 1, ndim

         rule_shift(ib) = Basis%tab_basis(ib)%smolyak_rule_shift_B
         rule_slope(ib) = Basis%tab_basis(ib)%smolyak_rule_slope_B

        IF (rule_shift(ib) <= 0) rule_shift(ib) = 1

         IF (rule_slope(ib) <= 0) rule_slope(ib) = 1
         rule_shift_g(ib) = Basis%tab_basis(ib)%smolyak_rule_shift_G
         rule_slope_g(ib) = Basis%tab_basis(ib)%smolyak_rule_slope_G
         IF (rule_shift_g(ib) <= 0) rule_shift_g(ib) = 1
         IF (rule_slope_g(ib) <= 0) rule_slope_g(ib) = 1

         nuclear_nb(ib) = Basis%tab_basis(ib)%nb

        nuclear_nq(ib) = Basis%tab_basis(ib)%nq

         ndend_b(ib)    = Basis%tab_basis(ib)%nb

         ndend_q(ib)    = Basis%tab_basis(ib)%nq

      END DO



     CALL Init_NDindex(NuclearGridIndex, ndend_q, ndim)

     CALL Init_NDindex(NuclearBasisIndex, ndend_b, ndim)

      CALL InitializeSmolyakSparseGridDescriptor(Basis%SmolyakSparseGrid, ndim, nsurf, Basis%LB, Basis%LG,&
      & rule_shift, rule_slope, rule_shift_g, rule_slope_g, NuclearBasisIndex, NuclearGridIndex, nuclear_nb, nuclear_nq, Basis%nb, Basis%nq)



     IF (ALLOCATED(NuclearGridIndex%Tab0)) DEALLOCATE(NuclearGridIndex%Tab0)

      IF (ALLOCATED(NuclearBasisIndex%Tab0)) DEALLOCATE(NuclearBasisIndex%Tab0)

     DEALLOCATE(rule_shift, rule_slope, nuclear_nb, nuclear_nq, ndend_q, ndend_b)

   END SUBROUTINE FinalizeSparseGridRepresentation



    SUBROUTINE Construct_Basis_Sin(Basis) ! sin : boxAB with A=0 and B=pi
      USE QDUtil_m

      TYPE(Basis_t), INTENT(INOUT)        :: Basis

      !<==================== Local variables ====================>
      LOGICAL :: debug
      REAL(KIND=Rkind)                    :: dx
      INTEGER                             :: ib, iq, nb, nq


      debug = Basis%debug .OR. Propa_Debug_Enabled


      IF (debug) THEN
         WRITE (out_unit, *) 'Construct_Basis_Sin:-------------------------------------------------'
         WRITE (out_unit, *) 'Basis%nb= ', Basis%nb
         WRITE (out_unit, *) 'Basis%nq= ', Basis%nq
      END IF

      IF (ALLOCATED(Basis%x)) DEALLOCATE (Basis%x)
      IF (ALLOCATED(Basis%w)) DEALLOCATE (Basis%w)
      IF (ALLOCATED(Basis%d0gb)) DEALLOCATE (Basis%d0gb)
      IF (ALLOCATED(Basis%d1gb)) DEALLOCATE (Basis%d1gb)
      IF (ALLOCATED(Basis%d2gb)) DEALLOCATE (Basis%d2gb)

      nb = Basis%nb
      nq = Basis%nq
      dx = pi/nq

      ! grid and weight
      Basis%x = [(dx*(iq - HALF), iq=1, nq)]
      Basis%w = [(dx, iq=1, nq)]

      ALLOCATE (Basis%d0gb(nq, nb))
      ALLOCATE (Basis%d1gb(nq, nb, 1))
      ALLOCATE (Basis%d2gb(nq, nb, 1, 1))


      DO ib = 1, nb
         DO iq = 1, nq
            CALL d0d1d2box(Basis%x(iq),Basis%d0gb(iq, ib),Basis%d1gb(iq, ib, 1)&
                           &,Basis%d2gb(iq, ib, 1, 1),ib)
         END DO
      END DO

      IF (nb == nq) THEN
         Basis%d0gb(:, nb) = Basis%d0gb(:, nb)/sqrt(TWO)
         Basis%d1gb(:, nb, :) = Basis%d1gb(:, nb, :)/sqrt(TWO)
         Basis%d2gb(:, nb, :, :) = Basis%d2gb(:, nb, :, :)/sqrt(TWO)
      END IF


      IF (debug) THEN
         WRITE (out_unit, *) 'End Construct_Basis_Sin'
      END IF

   END SUBROUTINE Construct_Basis_Sin




   !> Construct the Fourier basis in the interval [-π, π]
   !! This routine allocates the quadrature grid and computes the Fourier
   !! basis functions and their derivatives.
   SUBROUTINE Construct_Basis_Fourier(Basis)
      USE QDUtil_m
      IMPLICIT NONE
      !Input
      TYPE(Basis_t), INTENT(INOUT) :: Basis
      !Local variables
      REAL(KIND=Rkind)             :: dx
      INTEGER                      :: ib, iq, nb, nq
      LOGICAL :: debug


      debug = Basis%debug .OR. Propa_Debug_Enabled


      IF (debug) THEN
         WRITE (out_unit, *) 'Construct_Basis_Fourier:-------------------------------------------------'
         WRITE (out_unit, *) 'Basis%nb= ', Basis%nb
         WRITE (out_unit, *) 'Basis%nq= ', Basis%nq
      END IF   
   
      !---------------------------------------------------------------------------
      ! Deallocate existing arrays if already allocated
      IF (ALLOCATED(Basis%x))       DEALLOCATE(Basis%x)
      IF (ALLOCATED(Basis%w))       DEALLOCATE(Basis%w)
      IF (ALLOCATED(Basis%d0gb))    DEALLOCATE(Basis%d0gb)
      IF (ALLOCATED(Basis%d1gb))    DEALLOCATE(Basis%d1gb)
      IF (ALLOCATED(Basis%d2gb))    DEALLOCATE(Basis%d2gb)
      !---------------------------------------------------------------------------
   
      nb = Basis%nb
      nq = Basis%nq
   
      ! Step size for the uniform grid over [-π, π]
      dx = TWO * pi / nq
   
      !---------------------------------------------------------------------------
      ! Allocate and initialize grid points and weights
      ALLOCATE(Basis%x(nq))
      ALLOCATE(Basis%w(nq))
   
      Basis%x = [(iq * dx - dx / 2 - pi, iq = 1, nq)]
      Basis%w = [(dx, iq = 1, nq)]
      !---------------------------------------------------------------------------
   
      !---------------------------------------------------------------------------
      ! Allocate derivative matrices
      ALLOCATE(Basis%d0gb(nq, nb))
      ALLOCATE(Basis%d1gb(nq, nb, 1))
      ALLOCATE(Basis%d2gb(nq, nb, 1, 1))
      !---------------------------------------------------------------------------
   
      !---------------------------------------------------------------------------
      ! Compute the Fourier basis functions and their derivatives
      DO ib = 1, nb
         DO iq = 1, nq
            CALL d0d1d2fourier(Basis%x(iq), &
                               Basis%d0gb(iq, ib), &
                               Basis%d1gb(iq, ib, 1), &
                               Basis%d2gb(iq, ib, 1, 1), ib)
         END DO
      END DO
      !---------------------------------------------------------------------------
   
      !---------------------------------------------------------------------------
      ! Normalize the last basis function if conditions are met
      IF (nb == nq .AND. MOD(nb, 2) == 0) THEN
         Basis%d0gb(:, nb)        = Basis%d0gb(:, nb) / SQRT(TWO)
         Basis%d1gb(:, nb, :)     = Basis%d1gb(:, nb, :) / SQRT(TWO)
         Basis%d2gb(:, nb, :, :)  = Basis%d2gb(:, nb, :, :) / SQRT(TWO)
      END IF
      !---------------------------------------------------------------------------
   
   END SUBROUTINE Construct_Basis_Fourier





    !> @brief Constructs the 'el' (electronic) type basis
   SUBROUTINE Construct_Basis_el(Basis)
      USE QDUtil_m
      IMPLICIT NONE
    !Input
      TYPE(Basis_t), INTENT(inout) :: Basis
    !Local variables
      LOGICAL :: debug

      debug = Basis%debug .OR. Propa_Debug_Enabled
    
      IF (debug) THEN
         WRITE (out_unit, *) 'Construct_Basis_el:-------------------------------------------------'
         WRITE (out_unit, *) 'Basis%nb= ', Basis%nb
         WRITE (out_unit, *) 'Basis%nq= ', 0
      END IF
    
      ! Set number of grid points to zero for 'el' basis
      Basis%nq = 0
    
      RETURN

      IF (debug) THEN
         WRITE (out_unit, *) 'End Construct_Basis_el'
      END IF
      
    END SUBROUTINE Construct_Basis_el

   !> @brief Constructs the 'ho' (harmonic oscillator) type basis
   !! This subroutine initializes the basis functions and corresponding
   !! grid values for the harmonic oscillator using Hermite polynomials.
   !! It allocates and computes the grid points, weights, and the 0th, 1st, and 2nd
   !! derivatives of the Hermite polynomials at each grid point.
   SUBROUTINE Construct_Basis_Ho(Basis)
      USE QDUtil_m
      IMPLICIT NONE
    !Input
      TYPE(Basis_t), INTENT(inout) :: Basis
    !Local variables
      LOGICAL                      :: debug
      INTEGER                      :: iq, ib
      INTEGER                      :: nb, nq
    
      debug = Basis%debug .OR. Propa_Debug_Enabled
      IF (debug) THEN
         WRITE (out_unit, *) 'Construct_Basis_Ho:-------------------------------------------------'
         WRITE (out_unit, *) 'Basis%nb= ', Basis%nb
         WRITE (out_unit, *) 'Basis%nq= ', Basis%nq
      END IF
    
      ! Retrieve the number of basis functions and grid points
      nb = Basis%nb
      nq = Basis%nq
    
      ! Deallocate previously allocated arrays if necessary
      IF (ALLOCATED(Basis%x))       DEALLOCATE(Basis%x)
      IF (ALLOCATED(Basis%w))       DEALLOCATE(Basis%w)
      IF (ALLOCATED(Basis%d0gb))    DEALLOCATE(Basis%d0gb)
      IF (ALLOCATED(Basis%d1gb))    DEALLOCATE(Basis%d1gb)
      IF (ALLOCATED(Basis%d2gb))    DEALLOCATE(Basis%d2gb)
    
      ! Allocate grid points and weights
      ALLOCATE(Basis%x(nq))
      ALLOCATE(Basis%w(nq))
      CALL hercom(nq, Basis%x(:), Basis%w(:))  ! Hermite quadrature points and weights
    
      ! Allocate arrays for derivatives
      ALLOCATE(Basis%d0gb(nq, nb))
      ALLOCATE(Basis%d1gb(nq, nb, 1))
      ALLOCATE(Basis%d2gb(nq, nb, 1, 1))
    
      ! Compute Hermite polynomials and their derivatives at each grid point
      DO iq = 1, nq
        DO ib = 1, nb
          CALL d0d1d2poly_Hermite_exp(Basis%x(iq), ib - 1,                        &
               Basis%d0gb(iq, ib), Basis%d1gb(iq, ib, 1), Basis%d2gb(iq, ib, 1, 1), .TRUE.)
        END DO
      END DO

      IF (debug) THEN
         WRITE (out_unit, *) 'End Construct_Basis_Ho'
      END IF
     FLUSH(out_unit)
    END SUBROUTINE Construct_Basis_Ho


    SUBROUTINE Construct_Basis_Hag(Basis)
      USE QDUtil_m
      IMPLICIT NONE
      TYPE(Basis_t), INTENT(INOUT) :: Basis
   
      LOGICAL                      :: debug
      INTEGER                      :: iq, ib
      INTEGER                      :: nb, nq
      REAL(KIND=Rkind)             :: q0, sq, p0, a
      COMPLEX(KIND=Rkind)          :: alpha
      COMPLEX(KIND=Rkind)          :: bnm2, bnm1, bn, bnp1, bnp2
   
      q0    = Basis%Q0
      sq    = Basis%scaleQ
      debug = Basis%debug .OR. Propa_Debug_Enabled
      alpha = Basis%alpha
      p0    = Basis%imp_k
      nb    = Basis%nb
      nq    = Basis%nq
      a     = REAL(alpha, KIND=Rkind)  ! Re(alpha)
   
      ! Désallocation éventuelle
      IF (ALLOCATED(Basis%x))    DEALLOCATE(Basis%x)
      IF (ALLOCATED(Basis%w))    DEALLOCATE(Basis%w)
      IF (ALLOCATED(Basis%d0gb)) DEALLOCATE(Basis%d0gb)
      IF (ALLOCATED(Basis%d1gb)) DEALLOCATE(Basis%d1gb)
      IF (ALLOCATED(Basis%d2gb)) DEALLOCATE(Basis%d2gb)
      IF (ALLOCATED(Basis%d1bb)) DEALLOCATE(Basis%d1bb)
   
      ! Grille et poids quadrature
      ALLOCATE (Basis%x(nq))
      ALLOCATE (Basis%w(nq))
      CALL hercom(nq, Basis%x(:), Basis%w(:))
   
      ! Allocation des bases
      ALLOCATE (Basis%d0gb(nq, nb))
      ALLOCATE (Basis%d1bb(nb, nb, 4))
      ALLOCATE (Basis%d1gb(nq, nb, 5))
      ALLOCATE (Basis%d2gb(nq, nb, 1, 1))
   
      ! Construction de la base Hagedorn
      DO iq = 1, nq
         DO ib = 1, nb
            ! ib-1 = n
            CALL d0d1d2poly_Hag_exp(Basis%x(iq), ib-1, &
                 Basis%d0gb(iq, ib), Basis%d1gb(iq, ib, 1), Basis%d2gb(iq, ib, 1, 1), &
                 .TRUE., alpha, p0)
         END DO
      END DO
   
      ! === Construction des dérivées par rapport aux paramètres ===
      DO iq = 1, nq
         DO ib = 1, nb
            ! n = ib-1
            bn   = Basis%d0gb(iq, ib)
            bnm1 = CZERO
            bnm2 = CZERO
            bnp1 = CZERO
            bnp2 = CZERO
   
            IF (ib > 1)  bnm1 = Basis%d0gb(iq, ib-1)
            IF (ib > 2)  bnm2 = Basis%d0gb(iq, ib-2)
            IF (ib < nb) bnp1 = Basis%d0gb(iq, ib+1)
            IF (ib < nb-1) bnp2 = Basis%d0gb(iq, ib+2)
   
            ! ∂_{q0} b_n
            Basis%d1gb(iq,ib,2) = -CONJG(alpha)/SQRT(a) * SQRT(REAL(ib-1,Rkind)/TWO) * bnm1 &
                                  + alpha/SQRT(a) * SQRT(REAL(ib,Rkind)/TWO) * bnp1 &
                                  - EYE * p0 * bn
   
            ! ∂_{p} b_n
            Basis%d1gb(iq,ib,3) = EYE/SQRT(a) * ( SQRT(REAL(ib,Rkind)/TWO) * bnp1 &
                                  + SQRT(REAL(ib-1,Rkind)/TWO) * bnm1 )
   
            ! ∂_{a} b_n
            Basis%d1gb(iq,ib,4) = (ONE/(FOUR*a)) * ( SQRT(REAL((ib-1)*(ib-2),Rkind)) * bnm2 &
                                  - SQRT(REAL(ib*(ib+1),Rkind)) * bnp2 )
   
            ! ∂_{b} b_n
            Basis%d1gb(iq,ib,5) = -EYE/(FOUR*a) * ( SQRT(REAL(ib*(ib+1),Rkind)) * bnp2 &
                                  + REAL(2*ib-1,Rkind) * bn &
                                  + SQRT(REAL((ib-1)*(ib-2),Rkind)) * bnm2 )
         END DO
      END DO
   
      IF (debug) THEN
         WRITE(*,*) 'End Construct_Basis_Hag'
      END IF
   
   END SUBROUTINE Construct_Basis_Hag

   !> Check orthonormality of the basis set
   !! This subroutine checks if the basis functions are orthonormal (⟨b_i | b_j⟩ ≈ δ_ij)
   !! It also checks the overlaps with derivatives if requested (up to second order).

  SUBROUTINE CheckOrtho_Basis(Basis, nderiv)
   USE QDUtil_m
   IMPLICIT NONE
   !Input   
   TYPE(Basis_t), INTENT(IN)          :: Basis
   INTEGER, INTENT(IN)                :: nderiv
   !Local variables
   INTEGER                            :: ib, jb
   COMPLEX(KIND=Rkind), ALLOCATABLE      :: S(:,:)
   REAL(KIND=Rkind)                   :: Sii, Sij
   LOGICAL :: debug

   !-------------------------------------------------------------------------
   ! Check for incompatible basis type
   debug = Basis%debug .OR. Propa_Debug_Enabled
   IF (debug) THEN
      WRITE(out_unit, *) 'CheckOrtho_Basis:-------------------------------------------------'
   END IF
   IF (Basis%Basis_name == 'el') THEN
      WRITE(out_unit, *) 'This routine is not applicable to Basis ''el'''
      RETURN
   END IF
   !-------------------------------------------------------------------------

   ! Proceed only if basis data is allocated
   IF (Basis_IS_allocated(Basis)) THEN
      IF (debug) THEN
      WRITE(out_unit,*) 
      WRITE(out_unit,*) '--------------------------------------------------'
      WRITE(out_unit,*) ' Ortho Basis: ', TRIM(Basis%Basis_name)
      WRITE(out_unit,*) '--------------------------------------------------'
      END IF

      ! Compute overlap matrix: S = <d0b|d0b>
      S = matmul(conjg(Basis%d0bgw), Basis%d0gb)

      ! Check deviation from orthonormality
      Sii = ZERO
      DO ib = 1, Basis%nb
         IF (ABS(S(ib, ib) - ONE) > Sii) Sii = ABS(S(ib, ib) - ONE)
         S(ib, ib) = ZERO
      END DO
      Sij = MAXVAL(ABS(S))

      ! Formatted and aligned output
      WRITE(out_unit,'(A, 2ES12.5)') 'Basis: '//TRIM(Basis%Basis_name)//' | max|Sii-1|, max|Sij| = ', Sii, Sij

      !-------------------------------------------------------------------------
      ! Check overlap with first derivatives if required
      IF (debug) THEN
         IF (nderiv > 0) THEN
            WRITE(out_unit,*) 
            WRITE(out_unit,*) 'First derivative overlaps <d0b|d1b>:'
            S = MATMUL(CONJG(Basis%d0bgw), Basis%d1gb(:,:,1))
            CALL Write_VecMat(S, out_unit, 5, info='<d0b|d1b>')
         END IF
      END IF

      !-------------------------------------------------------------------------
      ! Check overlap with second derivatives if required
      IF (nderiv > 1) THEN
         IF (debug) THEN
            WRITE(out_unit,*) 
            WRITE(out_unit,*) 'Second derivative overlaps <d0b|d2b>:'
            S = MATMUL(CONJG(Basis%d0bgw), Basis%d2gb(:, :, 1, 1))
            CALL Write_VecMat(S, out_unit, 5, info='<d0b|d2b>')
         END IF
      END IF

   ELSE
      WRITE(out_unit,*) 'WARNING in CheckOrtho_Basis: The basis is not allocated. Basis=', TRIM(Basis%Basis_name)
   END IF

   IF (debug) THEN
      WRITE(out_unit,*) 
      WRITE(out_unit,*) '--------------------------------------------------'
      WRITE(out_unit,*) 'Ortho Basis: ', TRIM(Basis%Basis_name)
      WRITE(out_unit,*) '--------------------------------------------------'
   END IF

END SUBROUTINE CheckOrtho_Basis

!----------------------------------------------------------------------
! SUBROUTINE Scale_Basis(Basis, q0, sq)
!  --> Perform the change of variable: x = sqrt(a)*(q - q0)
!      with a = Re(alpha), sq = sqrt(a) provided directly as input.
!
! Inputs:
!   Basis : type(Basis_t) containing grid points x, weights w,
!            basis functions d0gb, first derivatives d1gb, second derivatives d2gb
!   q0    : translation point q0
!   sq    : sqrt(a), scaling factor given as input (do not recompute a)
!
! The subroutine applies:
!   - Translation and rescaling of points x -> q
!   - Adjusts integration weights w
!   - Applies normalization factors using sq directly:
!         a^{1/4}  = sqrt(sq)
!         a^{3/4}  = sq * sqrt(sq)
!         a^{5/4}  = sq^2 * sqrt(sq)
!----------------------------------------------------------------------

  SUBROUTINE Scale_Basis(Basis, q0, sq)
      USE QDUtil_m
      IMPLICIT NONE
      !<--- Input
      TYPE(Basis_t), INTENT(INOUT)  :: Basis
      REAL(KIND=Rkind), INTENT(IN)  :: q0, sq
   
      !<--- Local variables
      LOGICAL                      :: debug
      INTEGER                      :: iq, ib
      INTEGER                      :: nb, nq
      REAL(KIND=Rkind)             :: a   
      COMPLEX(KIND=Rkind)          :: alpha
      REAL(KIND=Rkind) :: sq_half, sq_3half, sq_5half

      debug = Basis%debug .OR. Propa_Debug_Enabled
      IF (debug) THEN
         WRITE(out_unit,*) 
         WRITE(out_unit,*) 'Begin Scale_Basis:-------------------------------------------------'
         WRITE(out_unit,*) 'Basis name: ', TRIM(Basis%Basis_name)
         WRITE(out_unit,*) 'q0: ', q0
         WRITE(out_unit,*) 'sq: ', sq
      END IF
   
   
      IF (Basis%nq == 0) RETURN
      IF (TO_lowercase(Basis%Basis_name) /= 'ho' ) RETURN
   
      IF (abs(sq) > ONETENTH**6 .AND. Basis_IS_allocated(Basis)) THEN
         sq_half  = sqrt(sq)
         sq_3half = sq * sq_half
         sq_5half = sq * sq * sq_half
   
         ! Transform grid points from x to q
         Basis%x(:) = q0 + Basis%x(:)/sq
         ! Adjust integration weights
         Basis%w(:) = Basis%w(:)/sq
   
         ! Apply normalization factor a^{1/4} = sqrt(sq)
         Basis%d0gb(:, :) = Basis%d0gb(:, :) * sq_half
   
         ! Scale derivatives accordingly
         Basis%d1gb(:, :, :) = Basis%d1gb(:, :, :) * sq_3half
         Basis%d2gb(:, :, :, :) = Basis%d2gb(:, :, :, :) * sq_5half
      ELSE
         WRITE(out_unit, *) ' ERROR in Scale_Basis'
         WRITE(out_unit, *) ' sq is too small or the basis is not allocated.'
         WRITE(out_unit,*) 'Basis name: ', TRIM(Basis%Basis_name)
         WRITE(out_unit,*) 'q0: ', q0
         WRITE(out_unit,*) 'sq: ', sq
         FLUSH(out_unit)
         STOP 'ERROR in Scale_Basis'
      END IF

      IF (debug) THEN
         WRITE(out_unit,*) 
         WRITE(out_unit,*) 'End Scale_Basis:-------------------------------------------------'
         FLUSH(out_unit)
      END IF
   END SUBROUTINE Scale_Basis

   SUBROUTINE Complete_construct_Basis(Basis)
      USE QDUtil_m
      IMPLICIT NONE
      !<--- Input
      TYPE(Basis_t), INTENT(INOUT) :: Basis
      !<--- Local variables
      REAL(KIND=Rkind) :: Q0, SQ0, a
      COMPLEX(KIND=Rkind) :: alpha
      LOGICAL :: debug
  
      !--------------------------------------------------------------
      ! Retrieve alpha and compute scaling factor SQ0 = sqrt(Re(alpha))
      !--------------------------------------------------------------
      debug = Basis%debug .OR. Propa_Debug_Enabled
      alpha = Basis%alpha
      a     = REAL(alpha, KIND=Rkind)
      SQ0   = SQRT(a)
  
      !--------------------------------------------------------------
      ! Retrieve the translation center Q0
      !--------------------------------------------------------------
      Q0 = Basis%Q0
  
      !--------------------------------------------------------------
      ! Apply scaling transformation: x = SQ0 * (q - Q0)
      !--------------------------------------------------------------
      call Scale_Basis(Basis, Q0, SQ0)
      !--------------------------------------------------------------
      ! Compute the transpose of d0gb (for efficiency in projections)
      !--------------------------------------------------------------
      call Calc_tranpose_d0gb(Basis)
      !--------------------------------------------------------------
      call Calc_dngg_grid(Basis, debug)
      !--------------------------------------------------------------
      ! Check orthonormality of the constructed basis up to 2nd derivative
      !--------------------------------------------------------------
      call CheckOrtho_Basis(Basis, nderiv=2)

      IF (debug) THEN
         WRITE(out_unit,*) 
         WRITE(out_unit,*) 'End Complete_construct_Basis:-------------------------------------------------'
         FLUSH(out_unit)
      END IF
  
  END SUBROUTINE 
 

    !> Compute the transpose of d0gb and apply the quadrature weights
    !! This subroutine computes the weighted transpose of the d0gb 
    !! matrix and stores it in d0bgw.
    SUBROUTINE Calc_tranpose_d0gb(Basis)
       USE QDUtil_m
       IMPLICIT NONE
       !<--- Input
       TYPE(Basis_t), INTENT(INOUT)      :: Basis
       !<--- Local variables
       INTEGER                           :: ib, nb
       LOGICAL :: debug


       debug = Basis%debug .OR. Propa_Debug_Enabled


       IF (debug) THEN
          WRITE(out_unit,*) 'BEGINNING Calc_tranpose_d0gb'
          WRITE(out_unit,*) 'Basis name: ', TRIM(Basis%Basis_name)
          FLUSH(out_unit)
       END IF
    
       !-------------------------------------------------------------------------
       ! Get the number of basis functions
       nb = Basis%nb
    
       !-------------------------------------------------------------------------
       ! If the output array already exists, deallocate it
       IF (ALLOCATED(Basis%d0bgw)) DEALLOCATE(Basis%d0bgw)
    
       !-------------------------------------------------------------------------
       ! Special case: skip computation if Basis is of type 'el'
       IF (Basis%Basis_name == 'el') RETURN
    
       !-------------------------------------------------------------------------
       ! Compute the transpose and apply weights
       Basis%d0bgw = TRANSPOSE(Basis%d0gb)
    
       DO ib = 1, nb
          Basis%d0bgw(ib, :) = Basis%d0bgw(ib, :) * Basis%w(:)
       END DO

       IF (debug) THEN
          WRITE(out_unit,*) 'End Calc_tranpose_d0gb'
          FLUSH(out_unit)
       END IF
    
    END SUBROUTINE Calc_tranpose_d0gb

    !> Computes several grid-based derivative matrices in the gg (grid-grid) representation
   !! using the precomputed gb (grid-basis) data and the weighted conjugate transpose (d0bgw).
   
    SUBROUTINE Calc_dngg_grid(Basis, analytic_d1bb)
      USE QDUtil_m
      IMPLICIT NONE
      TYPE(Basis_t), INTENT(INOUT) :: Basis
      LOGICAL, OPTIONAL, INTENT(IN) :: analytic_d1bb
      INTEGER :: nb, nq, ib, m, n
      LOGICAL :: use_analytic
      LOGICAL :: debug
      REAL(kind=Rkind) :: a, b, p, q0
      CHARACTER(len=Name_len) :: name
    
      debug = Basis%debug .OR. Propa_Debug_Enabled

      IF (debug) THEN
         WRITE(out_unit,*) 'BEGINNING Calc_dngg_grid'
         WRITE(out_unit,*) 'Basis name: ', TRIM(Basis%Basis_name)
         WRITE(out_unit,*) 'nb: ', Basis%nb
         WRITE(out_unit,*) 'nq: ', Basis%nq
         FLUSH(out_unit)
      END IF
    
   
       name = Basis%Basis_name 
      call string_uppercase_TO_lowercase(name)
      !-------------------------------------------------------------------------
      ! Determine calculation method for d1bb
      IF (PRESENT(analytic_d1bb)) THEN
         use_analytic = analytic_d1bb
         IF (debug) THEN
            WRITE(out_unit,*) 'd1bb calculation method: ', &
                 MERGE('ANALYTIC', 'NUMERIC ', use_analytic)
         END IF
      ELSE
         use_analytic = .FALSE. ! Default: numerical calculation
         IF (debug) THEN
            WRITE(out_unit,*) 'd1bb calculation method: NUMERIC (default)'
         END IF
      END IF
    
      !-------------------------------------------------------------------------
      ! Retrieve dimensions
      nb = Basis%nb
      nq = Basis%nq
    
      !-------------------------------------------------------------------------
      ! Deallocate previous arrays if allocated
      IF (ALLOCATED(Basis%d1gg)) DEALLOCATE(Basis%d1gg)
      IF (ALLOCATED(Basis%d2gg)) DEALLOCATE(Basis%d2gg)
      IF (ALLOCATED(Basis%d1bb)) DEALLOCATE(Basis%d1bb)
    
      !-------------------------------------------------------------------------
      ! Allocate arrays
      ALLOCATE(Basis%d1gg(nq, nq, 1))
      ALLOCATE(Basis%d2gg(nq, nq, 1, 1))
    
      !-------------------------------------------------------------------------
      ! Skip computation for 'el' basis type
      IF (name == 'el') RETURN
    
      !-------------------------------------------------------------------------
      ! Grid-grid matrices
      Basis%d1gg(:,:,1)    = MATMUL(Basis%d1gb(:,:,1), CONJG(Basis%d0bgw))
      Basis%d2gg(:,:,1,1)  = MATMUL(Basis%d2gb(:,:,1,1), CONJG(Basis%d0bgw))
    
      IF(name == 'ho') THEN
         !-------------------------------------------------------------------------
         ! Basis-basis matrices (weighted overlaps)
         ! d1bb(:,:,k) = < b_m | ∂_k b_n >
         ALLOCATE(Basis%d1bb(nb, nb, 5))
         
         IF (use_analytic) THEN
            ! Analytical calculation of d1bb  
            Basis%d1bb(:,:,:) = CZERO
            
            ! Retrieve basis parameters
            a = REAL(Basis%alpha,kind=Rkind)
            b = REAL(AIMAG(Basis%alpha), kind=Rkind)
            p = REAL(Basis%imp_k, kind=Rkind)
            q0 = REAL(Basis%Q0, kind=Rkind)
            
            DO n = 1, nb
               DO m = 1, nb
                  ! d1bb(m,n,1) = < b_m | ∂_x b_n > (not provided in analytical form)
                  ! Will be calculated numerically if needed
                  
                  ! d1bb(m,n,2) = < b_m | ∂_{q0} b_n >
                  IF (m == n) THEN
                     Basis%d1bb(m,n,2) = -EYE * p
                  ELSE IF (m == n+1) THEN
                     Basis%d1bb(m,n,2) = (a + EYE*b) / SQRT(TWO * a * n)
                  ELSE IF (m == n-1) THEN
                     Basis%d1bb(m,n,2) = -(a - EYE*b) / SQRT(TWO * a * (n-1))
                  END IF
                  
                  ! d1bb(m,n,3) = < b_m | ∂_p b_n >
                  IF (m == n-1) THEN
                     Basis%d1bb(m,n,3) = (EYE / SQRT(a)) * SQRT((n-ONE)/TWO)
                  ELSE IF (m == n+1) THEN
                     Basis%d1bb(m,n,3) = (EYE / SQRT(a)) * SQRT(n/TWO)
                  END IF
                  
                  ! d1bb(m,n,4) = < b_m | ∂_a b_n >
                  IF (m == n-2) THEN
                     Basis%d1bb(m,n,4) = SQRT((n-ONE)*(n-TWO)) / (FOUR * a)
                  ELSE IF (m == n+2) THEN
                     Basis%d1bb(m,n,4) = -SQRT(n*(n+ONE)) / (FOUR * a)
                  END IF
                  
                  ! d1bb(m,n,5) = < b_m | ∂_b b_n >
                  IF (m == n-2) THEN
                     Basis%d1bb(m,n,5) = -EYE * SQRT((n-ONE)*(n-TWO)) / (FOUR * a)
                  ELSE IF (m == n) THEN
                     Basis%d1bb(m,n,5) = -EYE *(TWO*n - ONE) / (FOUR * a)
                  ELSE IF (m == n+2) THEN
                     Basis%d1bb(m,n,5) = -EYE * SQRT(n*(n+ONE)) / (FOUR * a)
                  END IF
               END DO
            END DO
            
         ELSE
            ! Numerical calculation (original method)
            Basis%d1bb(:,:,1) = MATMUL(Basis%d0bgw, Basis%d1gb(:,:,1)) ! d/dx
            Basis%d1bb(:,:,2) = MATMUL(Basis%d0bgw, Basis%d1gb(:,:,2)) ! ∂/∂q0
            Basis%d1bb(:,:,3) = MATMUL(Basis%d0bgw, Basis%d1gb(:,:,3)) ! ∂/∂p
            Basis%d1bb(:,:,4) = MATMUL(Basis%d0bgw, Basis%d1gb(:,:,4)) ! ∂/∂a
            Basis%d1bb(:,:,5) = MATMUL(Basis%d0bgw, Basis%d1gb(:,:,5)) ! ∂/∂b
         END IF
      END IF
    
      !-------------------------------------------------------------------------
      ! Debug output
      IF (debug) THEN
         CALL Write_VecMat(Basis%d1gg(:,:,1), out_unit, 5, info='d1gg')
         WRITE(out_unit,*)
         CALL Write_VecMat(Basis%d2gg(:,:,1,1), out_unit, 5, info='d2gg')
         IF(Basis%Basis_name == 'ho') THEN
            DO ib = 2, 5
               WRITE(out_unit,*) 'd1bb(:,:,', ib, ') - ', &
                    MERGE('ANALYTIC', 'NUMERIC ', use_analytic)
               CALL Write_VecMat(Basis%d1bb(:,:,ib), out_unit, 5, info='d1bb')
            END DO
         END IF
         WRITE(out_unit,*) 'END Calc_dngg_grid'
         FLUSH(out_unit)
      END IF
    
    END SUBROUTINE Calc_dngg_grid

   SUBROUTINE test_basitogridgridtobasis(Basis)
      USE QDUtil_m
      !<--- Input
      TYPE(Basis_t), intent(in)       :: Basis
      !<--- Local variables
      LOGICAL :: debug
      COMPLEX(KIND=Rkind), ALLOCATABLE   :: G1(:), B1(:)
      COMPLEX(KIND=Rkind), ALLOCATABLE   :: G2(:), B2(:)
      COMPLEX(KIND=Rkind), ALLOCATABLE   :: B(:)
      REAL(KIND=Rkind), ALLOCATABLE       :: diff_g(:), diff_b(:)
      INTEGER                             :: iq, ndim
      ndim = size(Basis%tab_basis)

      debug = Basis%debug .OR. Propa_Debug_Enabled

      IF (debug) THEN
         WRITE(out_unit, *) 'BEGINNING test_basitogridgridtobasis'
         FLUSH(out_unit)
      END IF

      ALLOCATE (B(Get_TotalBasisCoefficientCount(Basis)))
      ALLOCATE (G1(Basis%nq*Basis%tab_basis(ndim)%nb))
      ALLOCATE (diff_g(Basis%nq*Basis%tab_basis(ndim)%nb))
      ALLOCATE (diff_b(Get_TotalBasisCoefficientCount(Basis)))

      ALLOCATE (G2(Basis%nq*Basis%tab_basis(ndim)%nb))
      ALLOCATE (B1(Get_TotalBasisCoefficientCount(Basis)))
      ALLOCATE (B2(Get_TotalBasisCoefficientCount(Basis)))

      B(:) = CZERO
      B1(:) = CONE
      G1(:) = CONE
      G2(:) = CZERO
      !print *, 'G in', G1
      CALL GridTOBasis_nD_cplx(B, G1, Basis)
      CALL BasisTOGrid_nD_cplx(G2, B, Basis)
      !print *, 'G out', G2

      WRITE(out_unit, *) '---------------------------------------------------------------------'

      diff_g(:) = ABS(G1(:) - G2(:))

      WRITE(out_unit, *) '---------------------------------------------------------------------'
      WRITE(out_unit, *) 'maxval(diff_g(:))=', maxval(diff_g(:))
      WRITE(out_unit, *) 'MINVAL(diff_g(:))=', MINVAL(diff_g(:))
      WRITE(out_unit, *) '---------------------------------------------------------------------'
      G2(:) = CZERO
      B2(:) = CZERO
      B1(:) = CONE
      !print *, 'B in', B1
      CALL BasisTOGrid_nD_cplx(G2, B1, Basis)
      CALL GridTOBasis_nD_cplx(B2, G2, Basis)
      !print *, 'B out', B2
      diff_b(:) = ABS(B2(:) - B1(:))
      WRITE(out_unit, *) '---------------------------------------------------------------------'
      WRITE(out_unit, *) 'maxval(diff_b(:))=', maxval(diff_b(:))
      WRITE(out_unit, *) 'MINVAL(diff_b(:))=', MINVAL(diff_b(:))
      WRITE(out_unit, *) '---------------------------------------------------------------------'

      IF (debug) THEN
         WRITE(out_unit, *) 'END test_basitogridgridtobasis'
         FLUSH(out_unit)
      END IF

   END SUBROUTINE test_basitogridgridtobasis

   !> Transform a function from grid representation (GG) to basis representation (BB)
   !! in one dimension, using complex arithmetic.
   !! The transformation is: BB = <basis|grid> * GG
   SUBROUTINE GridTOBasis_1D_cplx(BB, GG, Basis)
      USE QDUtil_m
      IMPLICIT NONE
      ! Input/Output arguments
      TYPE(Basis_t), intent(in), target        :: Basis
      COMPLEX(kind=Rkind), intent(inout)       :: BB(:, :, :)  ! Output: Basis representation
      COMPLEX(kind=Rkind), intent(in)          :: GG(:, :, :)  ! Input: Grid representation
      ! Optional debug output
      LOGICAL :: debug
      INTEGER                                  :: i1, i3
   
      !-------------------------------------------------------------------------
      ! Optional debug flushing
      debug = Basis%debug .OR. Propa_Debug_Enabled
      IF (debug) THEN
         WRITE(out_unit, *) 'BEGINNING GridTOBasis_1D_cplx'
         WRITE(out_unit, *) 'Size of GG(:, :, :) :: ', SIZE(GG, dim=1), SIZE(GG, dim=2), SIZE(GG, dim=3)
         WRITE(out_unit, *) 'Size of BB(:, :, :) :: ', SIZE(BB, dim=1), SIZE(BB, dim=2), SIZE(BB, dim=3)
         FLUSH(out_unit)
      END IF
   
      !-------------------------------------------------------------------------
      ! Initialize the BB array to zero
      BB = CZERO
   
      !-------------------------------------------------------------------------
      ! Loop over all dimensions (assumed to be 1D data with extra indices)
      DO i3 = 1, UBOUND(GG, dim=3)
         DO i1 = 1, UBOUND(GG, dim=1)
            ! Matrix multiplication: BB(i1,:,i3) = <basis|grid> * GG(i1,:,i3)
            BB(i1, :, i3) = MATMUL(CONJG(Basis%d0bgw), GG(i1, :, i3))
         END DO
      END DO
   
      !-------------------------------------------------------------------------
      ! Optional debug flushing
      IF (debug) THEN
         WRITE(out_unit, *) 'END GridTOBasis_1D_cplx'
         FLUSH(out_unit)
      END IF
   
   END SUBROUTINE GridTOBasis_1D_cplx

   !> Transform a function from basis representation (BB) to grid representation (GB)
   !! in one dimension, using complex arithmetic.
   !! The transformation is: GB = <x|basis> * BB
   SUBROUTINE BasisTOGrid_1D_cplx(GB, BB, Basis)
      USE QDUtil_m
      IMPLICIT NONE
      ! Input/Output arguments
      TYPE(Basis_t), intent(in), target         :: Basis
      COMPLEX(kind=Rkind), intent(inout)        :: GB(:, :, :)  ! Output: Grid representation
      COMPLEX(kind=Rkind), intent(in)           :: BB(:, :, :)  ! Input: Basis representation
      ! Optional debug output
      LOGICAL :: debug
      INTEGER                                   :: i1, i3
   
      !-------------------------------------------------------------------------
      ! Optional debug flushing before processing
      debug = Basis%debug .OR. Propa_Debug_Enabled
      IF (debug) THEN
         WRITE(out_unit, *) 'BEGINNING BasisTOGrid_1D_cplx'
         WRITE(out_unit, *) 'Size of GB(:, :, :) :: ', SIZE(GB, dim=1), SIZE(GB, dim=2), SIZE(GB, dim=3)
         WRITE(out_unit, *) 'Size of BB(:, :, :) :: ', SIZE(BB, dim=1), SIZE(BB, dim=2), SIZE(BB, dim=3)
         FLUSH(out_unit)
      END IF
   
      !-------------------------------------------------------------------------
      ! Initialize the GB array to zero
      GB = CZERO
   
      !-------------------------------------------------------------------------
      ! Loop over all dimensions (i1 and i3 are loop indices for external dimensions)
      DO i3 = 1, UBOUND(BB, dim=3)
         DO i1 = 1, UBOUND(BB, dim=1)
            ! Matrix multiplication:
            ! GB(i1,:,i3) = <x|basis> * BB(i1,:,i3)
            GB(i1, :, i3) = MATMUL(Basis%d0gb(:, :), BB(i1, :, i3))
         END DO
      END DO
   
      !-------------------------------------------------------------------------
      ! Optional debug flushing after processing
      IF (debug) THEN
         WRITE(out_unit, *) 'END BasisTOGrid_1D_cplx'
         FLUSH(out_unit)
      END IF
   
   END SUBROUTINE BasisTOGrid_1D_cplx



   !> Calculates indexing parameters (Ib1, Ib2, Ib3, Iq1, Iq2, Iq3) 
   !! used for reshaping and accessing multidimensional basis/grid data.
   !! This subroutine is typically used in multidimensional simulations (e.g., 3D quantum grids).
   SUBROUTINE Calc_index(Ib1, Ib2, Ib3, Iq1, Iq2, Iq3, Basis)
      ! Input/Output arguments
      IMPLICIT NONE
      ! Input arguments
      TYPE(Basis_t),INTENT(IN), TARGET              :: Basis
      ! Output arguments
      INTEGER, INTENT(INOUT), ALLOCATABLE, OPTIONAL  :: Iq1(:), Iq2(:), Iq3(:)
      INTEGER, INTENT(INOUT), ALLOCATABLE, OPTIONAL  :: Ib1(:), Ib2(:), Ib3(:)
      INTEGER                                        :: Ndim      ! Number of spatial dimensions
      INTEGER                                        :: inb       ! Loop index
      LOGICAL :: debug


      debug = Basis%debug .OR. Propa_Debug_Enabled


      IF(debug) THEN
         WRITE(out_unit, *) 'BEGINNING Calc_index'
         WRITE(out_unit, *) 'Size of Basis%tab_basis :: ', SIZE(Basis%tab_basis)-1
         FLUSH(out_unit)
      END IF
      ! Determine the number of spatial dimensions (tab_basis has Ndim+1 entries)
      Ndim = SIZE(Basis%tab_basis) - 1

      ! Allocate index arrays if they are present in the arguments
      IF (PRESENT(Ib3)) ALLOCATE(Ib3(Ndim))
      IF (PRESENT(Ib2)) ALLOCATE(Ib2(Ndim))
      IF (PRESENT(Ib1)) ALLOCATE(Ib1(Ndim))
      IF (PRESENT(Iq3)) ALLOCATE(Iq3(Ndim))
      IF (PRESENT(Iq2)) ALLOCATE(Iq2(Ndim))
      IF (PRESENT(Iq1)) ALLOCATE(Iq1(Ndim))

   ! Loop over each dimension to calculate indexing terms
   DO inb = 1, Ndim
      IF (inb == 1) THEN
         ! First dimension
         IF (PRESENT(Iq1)) Iq1(1) = 1
         IF (PRESENT(Ib1)) Ib1(1) = 1
         IF (PRESENT(Iq2)) Iq2(1) = Basis%tab_basis(1)%nq
         IF (PRESENT(Ib2)) Ib2(1) = Basis%tab_basis(1)%nb
         IF (PRESENT(Iq3)) Iq3(1) = PRODUCT(Basis%tab_basis(2:Ndim)%nq) * Basis%tab_basis(Ndim + 1)%nb
         IF (PRESENT(Ib3)) Ib3(1) = PRODUCT(Basis%tab_basis(2:Ndim + 1)%nb)

      ELSE IF (inb == Ndim) THEN
         ! Last dimension
         IF (PRESENT(Iq1)) Iq1(inb) = PRODUCT(Basis%tab_basis(1:Ndim - 1)%nq)
         IF (PRESENT(Ib1)) Ib1(inb) = PRODUCT(Basis%tab_basis(1:Ndim - 1)%nb)
         IF (PRESENT(Iq2)) Iq2(inb) = Basis%tab_basis(Ndim)%nq
         IF (PRESENT(Ib2)) Ib2(inb) = Basis%tab_basis(Ndim)%nb
         IF (PRESENT(Iq3)) Iq3(inb) = Basis%tab_basis(Ndim + 1)%nb
         IF (PRESENT(Ib3)) Ib3(inb) = Basis%tab_basis(Ndim + 1)%nb

      ELSE
         ! Intermediate dimensions
         IF (PRESENT(Iq1)) Iq1(inb) = PRODUCT(Basis%tab_basis(1:inb - 1)%nq)
         IF (PRESENT(Ib1)) Ib1(inb) = PRODUCT(Basis%tab_basis(1:inb - 1)%nb)
         IF (PRESENT(Iq2)) Iq2(inb) = Basis%tab_basis(inb)%nq
         IF (PRESENT(Ib2)) Ib2(inb) = Basis%tab_basis(inb)%nb
         IF (PRESENT(Iq3)) Iq3(inb) = PRODUCT(Basis%tab_basis(inb + 1:Ndim)%nq) * Basis%tab_basis(Ndim + 1)%nb
         IF (PRESENT(Ib3)) Ib3(inb) = PRODUCT(Basis%tab_basis(inb + 1:Ndim + 1)%nb)
      END IF
   END DO

   IF(debug) THEN
      WRITE(out_unit, *) 'END Calc_index'
      FLUSH(out_unit)
   END IF

END SUBROUTINE Calc_index



SUBROUTINE BasisTOGrid_nD_cplx_dense(G, B, Basis)
   IMPLICIT NONE

   ! Input/Output arguments
   TYPE(Basis_t), INTENT(IN)                  :: Basis          ! Basis definition and metadata
   COMPLEX(KIND=Rkind), INTENT(IN), TARGET    :: B(:)           ! Input in Basis representation
   COMPLEX(KIND=Rkind), INTENT(INOUT), TARGET :: G(:)           ! Output in Grid representation

   ! Optional debug output
   LOGICAL :: debug
   ! Pointers for reshaped views of arrays
   COMPLEX(KIND=Rkind), POINTER :: BBG(:, :, :), BBB(:, :, :), GGB(:, :, :)
   COMPLEX(KIND=Rkind), POINTER :: GBB(:, :, :)
   ! Indexing variables (flattening multidimensional data)
   INTEGER :: Ib3, Iq1, Iq2, Ib1, Ib2, Iq3
   ! Temporary working arrays
   COMPLEX(KIND=Rkind), ALLOCATABLE, TARGET :: GBB1(:), GGB2(:)
   ! Integer variables
   INTEGER :: inb, ndim

 
    debug = Basis%debug .OR. Propa_Debug_Enabled


    IF(debug) THEN
       WRITE(out_unit, *) 'BEGINNING BasisTOGrid_nD_cplx'
       WRITE(out_unit, *) 'intent(in) :: B(:)', SIZE(B)
       WRITE(out_unit, *) 'intent(inout) :: G(:)', SIZE(G)
       CALL Write_VecMat(B, out_unit, 5, info='B')
       FLUSH(out_unit)
    END IF

   ! Check that basis is allocated
   IF (.NOT. Basis_IS_Allocated(Basis)) THEN
       WRITE(out_unit, *) ' ERROR: the basis is not allocated.'
       STOP 'ERROR BasisTOGrid_Basis: the basis is not allocated.'
   END IF

   ! Check that the sizes of B and G are compatible with the basis structure
   IF (SIZE(B) /= Basis%nb * Basis%tab_basis(SIZE(Basis%tab_basis))%nb) THEN
       WRITE(out_unit, *) ' ERROR: the size of B does not match expected nb.'
       STOP 'ERROR in BasisTOGrid_Basis: wrong B size.'
   END IF
   IF (SIZE(G) /= Basis%nq * Basis%tab_basis(SIZE(Basis%tab_basis))%nb) THEN
       WRITE(out_unit, *) ' ERROR: the size of G does not match expected nq.'
       STOP 'ERROR in BasisTOGrid_Basis: wrong G size.'
   END IF

   ! Determine number of dimensions (excluding electronic dimension)
   ndim = SIZE(Basis%tab_basis) - 1

   ! --- Special Case: 1D Transform ---
   IF (ndim == 1) THEN
       ! Set index factors for single-dimension case
       Iq1 = 1
       Iq2 = Basis%tab_basis(1)%nq
       Iq3 = Basis%tab_basis(ndim + 1)%nb

       Ib1 = 1
       Ib2 = Basis%tab_basis(1)%nb
       Ib3 = PRODUCT(Basis%tab_basis(2:ndim + 1)%nb)

       ! Reshape input B and output G into 3D slices
       BBB(1:Ib1, 1:Ib2, 1:Ib3) => B
       GGB(1:Iq1, 1:Iq2, 1:Iq3) => G

       ! Initialize output
       G = CZERO

       ! Perform transformation for 1D case
       CALL BasisTOGrid_1D_cplx(GGB, BBB, Basis%tab_basis(1))

   ELSE
   ! --- General Case: Multi-Dimensional Transform ---
       ! Initialize first transformation (dimension 1)
       Iq1 = 1
       Iq2 = Basis%tab_basis(1)%nq
       Iq3 = Basis%tab_basis(ndim + 1)%nb

       Ib1 = 1
       Ib2 = Basis%tab_basis(1)%nb
       Ib3 = PRODUCT(Basis%tab_basis(2:ndim + 1)%nb)

       ALLOCATE(GBB1(Iq1 * Iq2 * Ib3))
       BBB(1:Iq1, 1:Ib2, 1:Ib3) => B
       GBB(1:Iq1, 1:Iq2, 1:Ib3) => GBB1
       GBB1 = CZERO

       ! First 1D transformation
       CALL BasisTOGrid_1D_cplx(GBB, BBB, Basis%tab_basis(1))

       ! Loop over intermediate dimensions (2 to ndim-1)
       DO inb = 2, ndim - 1
           Iq1 = PRODUCT(Basis%tab_basis(1:inb - 1)%nq)
           Ib1 = PRODUCT(Basis%tab_basis(1:inb - 1)%nb)
           Iq2 = Basis%tab_basis(inb)%nq
           Ib2 = Basis%tab_basis(inb)%nb
           Ib3 = PRODUCT(Basis%tab_basis(inb + 1:ndim + 1)%nb)
           Iq3 = PRODUCT(Basis%tab_basis(inb + 1:ndim)%nq) * Basis%tab_basis(ndim + 1)%nb

           ! Allocate temporary output
           ALLOCATE(GGB2(Iq1 * Iq2 * Ib3))

           ! Set reshaped views for input/output
           BBB(1:Iq1, 1:Ib2, 1:Ib3) => GBB1
           GBB(1:Iq1, 1:Iq2, 1:Ib3) => GGB2

           ! Perform 1D transformation for dimension `inb`
           CALL BasisTOGrid_1D_cplx(GBB, BBB, Basis%tab_basis(inb))

           ! Move result to input buffer for next iteration
           GBB1 = GGB2
           DEALLOCATE(GGB2)
       END DO

       ! Final transformation for last dimension
       Iq1 = PRODUCT(Basis%tab_basis(1:ndim - 1)%nq)
       Ib1 = PRODUCT(Basis%tab_basis(1:ndim - 1)%nb)
       Ib2 = Basis%tab_basis(ndim)%nb
       Iq2 = Basis%tab_basis(ndim)%nq
       Ib3 = Basis%tab_basis(ndim + 1)%nb
       Iq3 = Basis%tab_basis(ndim + 1)%nb

       BBB(1:Iq1, 1:Ib2, 1:Ib3) => GBB1
       GBB(1:Iq1, 1:Iq2, 1:Ib3) => G

       CALL BasisTOGrid_1D_cplx(GBB, BBB, Basis%tab_basis(ndim))
   END IF

   ! Debugging: output result
   IF (debug) THEN
       WRITE(out_unit, *) 'intent(INOUT) :: G(:)'
       CALL Write_VecMat(G, out_unit, 5, info='G')
       WRITE(out_unit, *) 'END BasisTOGrid_nD_cplx'
       FLUSH(out_unit)
   END IF

END SUBROUTINE

SUBROUTINE GridTOBasis_nD_cplx_dense(B, G, Basis)
   USE QDUtil_m
   !Input/Output arguments
   TYPE(Basis_t), intent(in), target                :: Basis
   COMPLEX(KIND=Rkind), intent(in), target          :: G(:)         ! Input: Grid representation
   COMPLEX(KIND=Rkind), intent(inout), target       :: B(:)         ! Output: Basis representation
   !Pointers for reshaped views of arrays
   COMPLEX(KIND=Rkind), POINTER                     :: BBB(:,:,:), GGB(:,:,:) 
   COMPLEX(KIND=Rkind), ALLOCATABLE, TARGET         :: BGG1(:), BGG2(:)
   COMPLEX(KIND=Rkind), POINTER                     :: GGG(:,:,:) 
   INTEGER                                          :: inb, Ndim
   INTEGER                                          :: Ib1, Ib2, Iq3, Iq1, Iq2, Ib3
   LOGICAL :: debug

   ! Optional debug output
   debug = Basis%debug .OR. Propa_Debug_Enabled
   IF (debug) THEN
      write (out_unit, *) 'BEGINNING GridTOBasis_nD_cplx'
      WRITE(out_unit, *) 'intent(in) :: B(:)', SIZE(B)
      WRITE(out_unit, *) 'intent(in) :: G(:)', SIZE(G)
      CALL Write_VecMat(G, out_unit, 5, info='G')
      flush (out_unit)
   END IF

   ! Check that basis is properly allocated
   IF (.NOT. Basis_IS_Allocated(Basis)) THEN
      write (out_unit, *) ' ERROR in BasisTOGrid_nD_cplx'
      STOP "ERROR BasisTOGrid_Basis: the basis is not Allocated."
   END IF

   ! Check input and output array sizes
   IF (size(B) /= Basis%nb * Basis%tab_basis(size(Basis%tab_basis))%nb) THEN
      write (out_unit, *) ' ERROR: Wrong size for B'
      STOP 'ERROR in GridTOBasis_Basis: wrong B size.'
   END IF

   IF (size(G) /= Basis%nq * Basis%tab_basis(size(Basis%tab_basis))%nb) THEN
      write (out_unit, *) ' ERROR: Wrong size for G'
      STOP 'ERROR in GridTOBasis_Basis: wrong G size'
   END IF

   ndim = size(Basis%tab_basis) - 1   ! Number of spatial dimensions

   ! Special case: 1D
   IF (ndim == 1) THEN
      Iq1 = 1
      Iq2 = Basis%tab_basis(1)%nq
      Iq3 = Basis%tab_basis(ndim + 1)%nb

      Ib1 = 1
      Ib2 = Basis%tab_basis(1)%nb
      Ib3 = Product(Basis%tab_basis(2:ndim+1)%nb)

      B = CZERO
      GGB(1:Iq1, 1:Iq2, 1:Iq3) => G
      BBB(1:Ib1, 1:Ib2, 1:Ib3) => B
      Call GridTOBasis_1D_cplx(BBB, GGB, Basis%tab_basis(1))

   ELSE
      ! Initialization for first dimension
      Iq1 = 1
      Iq2 = Basis%tab_basis(1)%nq
      Iq3 = Product(Basis%tab_basis(2:ndim)%nq) * Basis%tab_basis(ndim + 1)%nb

      Ib1 = 1
      Ib2 = Basis%tab_basis(1)%nb
      Ib3 = Product(Basis%tab_basis(2:ndim+1)%nb)

      allocate (BGG1(Ib1*Ib2*Iq3))
      BGG1 = CZERO
      GGG(1:Ib1, 1:Iq2, 1:Iq3) => G
      GGB(1:Ib1, 1:Ib2, 1:Iq3) => BGG1

      Call GridTOBasis_1D_cplx(GGB, GGG, Basis%tab_basis(1))

      ! Iterative transform for inner dimensions
      DO inb = 2, ndim - 1
         Iq1 = Product(Basis%tab_basis(1:inb - 1)%nq)
         Ib1 = Product(Basis%tab_basis(1:inb - 1)%nb)
         Iq2 = Basis%tab_basis(inb)%nq
         Ib2 = Basis%tab_basis(inb)%nb
         Ib3 = Product(Basis%tab_basis(inb + 1:ndim + 1)%nb)
         Iq3 = Product(Basis%tab_basis(inb + 1:ndim)%nq) * Basis%tab_basis(ndim + 1)%nb 

         allocate (BGG2(Ib1*Ib2*Iq3))
         BGG2 = CZERO
         GGG(1:Ib1, 1:Iq2, 1:Iq3) => BGG1
         GGB(1:Ib1, 1:Ib2, 1:Iq3) => BGG2

         Call GridTOBasis_1D_cplx(GGB, GGG, Basis%tab_basis(inb))
         BGG1 = BGG2
         deallocate (BGG2)
      END DO

      ! Final step for last dimension
      B(:) = CZERO

      Iq1 = Product(Basis%tab_basis(1:ndim - 1)%nq)
      Ib1 = Product(Basis%tab_basis(1:ndim - 1)%nb)
      Ib2 = Basis%tab_basis(ndim)%nb
      Iq2 = Basis%tab_basis(ndim)%nq
      Ib3 = Basis%tab_basis(ndim + 1)%nb
      Iq3 = Basis%tab_basis(ndim + 1)%nb

      GGG(1:Ib1, 1:Iq2, 1:Iq3) => BGG1
      GGB(1:Ib1, 1:Ib2, 1:Iq3) => B
      Call GridTOBasis_1D_cplx(GGB, GGG, Basis%tab_basis(ndim))
   END IF

   IF (debug) THEN
      CALL Write_VecMat(B, out_unit, 5, info='B') 
      WRITE (out_unit, *) 'END GridTOBasis_nD_cplx' 
      flush (out_unit)
   END IF
END SUBROUTINE

 SUBROUTINE BasisTOGrid_nD_cplx_smolyak(G, B, Basis)
   USE QDUtil_m, ONLY: Rkind, CZERO
   TYPE(Basis_t), INTENT(IN) :: Basis
   COMPLEX(KIND=Rkind), INTENT(IN) :: B(:)
   COMPLEX(KIND=Rkind), INTENT(INOUT) :: G(:)
   COMPLEX(KIND=Rkind), ALLOCATABLE :: FullB(:)
   INTEGER :: nfull
   nfull=Basis%nb*Basis%tab_basis(SIZE(Basis%tab_basis))%nb
   IF (SIZE(B) /= GetSmolyakPackedBasisSize(Basis%SmolyakSparseGrid)) &
      STOP 'BasisTOGrid_nD_cplx_smolyak: wrong sparse size'
   IF (SIZE(G) /= Basis%nq*Basis%tab_basis(SIZE(Basis%tab_basis))%nb) &
      STOP 'BasisTOGrid_nD_cplx_smolyak: wrong grid size'
   ALLOCATE(FullB(nfull))
   CALL ExpandSparseBasisToFullBasis(FullB,B,Basis%SmolyakSparseGrid)
   CALL BasisTOGrid_nD_cplx_dense(G,FullB,Basis)
   DEALLOCATE(FullB)
 END SUBROUTINE BasisTOGrid_nD_cplx_smolyak

 SUBROUTINE GridTOBasis_nD_cplx_smolyak(B, G, Basis)
   USE QDUtil_m, ONLY: Rkind
   TYPE(Basis_t), INTENT(IN) :: Basis
   COMPLEX(KIND=Rkind), INTENT(IN) :: G(:)
   COMPLEX(KIND=Rkind), INTENT(INOUT) :: B(:)
   COMPLEX(KIND=Rkind), ALLOCATABLE :: FullB(:)
   INTEGER :: nfull
   nfull=Basis%nb*Basis%tab_basis(SIZE(Basis%tab_basis))%nb
   IF (SIZE(B) /= GetSmolyakPackedBasisSize(Basis%SmolyakSparseGrid)) &
      STOP 'GridTOBasis_nD_cplx_smolyak: wrong sparse size'
   IF (SIZE(G) /= Basis%nq*Basis%tab_basis(SIZE(Basis%tab_basis))%nb) &
      STOP 'GridTOBasis_nD_cplx_smolyak: wrong grid size'
   ALLOCATE(FullB(nfull))
   CALL GridTOBasis_nD_cplx_dense(FullB,G,Basis)
   CALL RestrictFullBasisToSparseBasis(B,FullB,Basis%SmolyakSparseGrid)
   DEALLOCATE(FullB)
 END SUBROUTINE GridTOBasis_nD_cplx_smolyak

 SUBROUTINE decode_block_index(linear, caps, idx)
   INTEGER, INTENT(IN) :: linear, caps(:)
   INTEGER, INTENT(OUT) :: idx(:)
   INTEGER :: r,j
   r=linear-1
   DO j=1,SIZE(caps); idx(j)=MOD(r,caps(j))+1; r=r/caps(j); END DO
 END SUBROUTINE decode_block_index

 SUBROUTINE BasisTOGrid_nD_cplx(G, B, Basis)

   COMPLEX(KIND=Rkind), INTENT(INOUT) :: G(:)

   COMPLEX(KIND=Rkind), INTENT(IN)    :: B(:)

   TYPE(Basis_t),       INTENT(IN), TARGET :: Basis

   COMPLEX(KIND=Rkind), ALLOCATABLE :: FullBasisVector(:)

   INTEGER :: DenseBasisSize



   IF (Basis_UsesSmolyakSparseGrid(Basis) .AND. SIZE(B) == GetSmolyakPackedBasisSize(Basis%SmolyakSparseGrid)) THEN

      DenseBasisSize = Basis%nb * Basis%tab_basis(SIZE(Basis%tab_basis))%nb
      ALLOCATE(FullBasisVector(DenseBasisSize))
      CALL ExpandSparseBasisToFullBasis(FullBasisVector, B, Basis%SmolyakSparseGrid)
      CALL BasisTOGrid_nD_cplx_dense(G, FullBasisVector, Basis)
      DEALLOCATE(FullBasisVector)

   ELSE

      CALL BasisTOGrid_nD_cplx_dense(G, B, Basis)

   END IF

END SUBROUTINE BasisTOGrid_nD_cplx


SUBROUTINE GridTOBasis_nD_cplx(B, G, Basis)

   COMPLEX(KIND=Rkind), INTENT(INOUT) :: B(:)

   COMPLEX(KIND=Rkind), INTENT(IN)    :: G(:)

   TYPE(Basis_t),       INTENT(IN), TARGET :: Basis

   COMPLEX(KIND=Rkind), ALLOCATABLE :: FullBasisVector(:)

   INTEGER :: DenseBasisSize


   IF (Basis_UsesSmolyakSparseGrid(Basis) .AND. SIZE(B) == GetSmolyakPackedBasisSize(Basis%SmolyakSparseGrid)) THEN

      DenseBasisSize = Basis%nb * Basis%tab_basis(SIZE(Basis%tab_basis))%nb
      ALLOCATE(FullBasisVector(DenseBasisSize))
      CALL GridTOBasis_nD_cplx_dense(FullBasisVector, G, Basis)
      CALL RestrictFullBasisToSparseBasis(B, FullBasisVector, Basis%SmolyakSparseGrid)
      DEALLOCATE(FullBasisVector)

   ELSE

      CALL GridTOBasis_nD_cplx_dense(B, G, Basis)

   END IF

END SUBROUTINE GridTOBasis_nD_cplx


SUBROUTINE Calc_Q_grid(Q, Basis)
   IMPLICIT NONE

   ! Input
   TYPE(Basis_t), INTENT(IN)                 :: Basis

   ! Local variables
   INTEGER, ALLOCATABLE                      :: Tab_iq(:)
   INTEGER                                   :: inb, ndim, iq
   REAL(KIND=Rkind), INTENT(INOUT), ALLOCATABLE :: Q(:, :)
   LOGICAL                                   :: Endloop

   ! Determine the number of dimensions excluding the electronic basis
   ndim = SIZE(Basis%tab_basis) - 1

   ! Allocate output array Q and index array Tab_iq.
   ! Q is reused by the direct curvilinear kernels.
   IF (ALLOCATED(Q)) DEALLOCATE(Q)
   ALLOCATE(Q(Basis%nq, ndim), Tab_iq(ndim))

   ! Initialize the index array for multidimensional looping
   CALL Init_tab_ind(Tab_iq, Basis%NDindexq)

   iq = 0
   DO
       iq = iq + 1

       ! Update the multi-dimensional index
       CALL increase_NDindex(Tab_iq, Basis%NDindexq, Endloop)
       IF (Endloop) EXIT

       ! Fill the Q array with the appropriate basis grid values
       DO inb = 1, ndim
           Q(iq, inb) = Basis%tab_basis(inb)%x(Tab_iq(inb))
       END DO

       ! Optional debug print
       ! PRINT*, iq, Q(iq, :)

   END DO

   ! Deallocate local array to avoid memory leak
   DEALLOCATE(Tab_iq)


END SUBROUTINE Calc_Q_grid

 SUBROUTINE Build_NUFT_Grids(Basis, Q, WQ, P, WP, Q0d)
    IMPLICIT NONE
  
    !----------------------------------------------------------------
    ! Build_NUFT_Grids :
    !   - Fill Q(npoints_total, ndim) with position nodes (possibly non-uniform)
    !   - Fill WQ(npoints_total) with product quadrature weights
    !   - Build regular, centered P(npoints_total, ndim) grid suitable for NUFFT:
    !       P component per-dimension: p_k = delta_p_j * (k_index - center_j)
    !   - WP(:) holds the (uniform) product weight in p-space = prod_j delta_p_j
    !   - Q0d(j) returns the center offset Q0 for each physical dimension (if present)
    !
    ! Assumptions:
    !   - Basis%tab_basis(inb)%x(:)     : node coordinates (1..N_j)
    !   - Basis%tab_basis(inb)%w(:)     : quadrature weights (1..N_j)
    !   - Basis%tab_basis(inb)%nq       : number of points along dimension inb (integer)
    !   - Basis%tab_basis(inb)%Q0 (opt) : center offset of the 1D grid (real), optional
    !   - Basis%nq                       : total number of grid points (product of N_j)
    !----------------------------------------------------------------
  
    ! Arguments
    TYPE(Basis_t), INTENT(IN)                       :: Basis
    REAL(KIND=Rkind), INTENT(INOUT), ALLOCATABLE    :: Q(:, :)
    REAL(KIND=Rkind), INTENT(INOUT), ALLOCATABLE    :: WQ(:)
    REAL(KIND=Rkind), INTENT(INOUT), ALLOCATABLE    :: P(:, :)
    REAL(KIND=Rkind), INTENT(INOUT), ALLOCATABLE    :: WP(:)
    REAL(KIND=Rkind), INTENT(OUT), ALLOCATABLE      :: Q0d(:)   ! center per dim (may be 0 if not provided)
  
    ! Locals
    INTEGER, ALLOCATABLE                            :: Tab_iq(:)
    INTEGER                                         :: inb, ndim, iq, ik, nq_total, Nj
    REAL(KIND=Rkind),ALLOCATABLE                    :: L(:)
    REAL(KIND=Rkind),ALLOCATABLE                    :: x_min(:)
    REAL(KIND=Rkind),ALLOCATABLE                    :: x_max(:)
    REAL(KIND=Rkind),ALLOCATABLE                    :: delta_p(:)
    REAL(KIND=Rkind)                                :: prod_dp, p_step
    LOGICAL                                         :: Endloop
    INTEGER                                         :: center_offset, idx1d
    INTEGER                                         :: expected_total, tmpN
  
  
    ! 1) determine dimensionality and expected sizes
    ndim = SIZE(Basis%tab_basis) - 1    ! excluding electronic index
    IF (ndim <= 0) THEN
       WRITE(*,*) 'Build_NUFT_Grids: ndim <= 0, nothing to do.'
       RETURN
    END IF
  
    nq_total = Basis%nq   ! total number of product points expected
    IF (nq_total <= 0) THEN
       WRITE(*,*) 'Build_NUFT_Grids: Basis%nq must be > 0'
       RETURN
    END IF
  
    ! allocate locals
    ALLOCATE(Tab_iq(ndim))
    ALLOCATE(L(ndim))
    ALLOCATE(x_min(ndim))
    ALLOCATE(x_max(ndim))
    ALLOCATE(delta_p(ndim))
  
    ! allocate / reallocate outputs if necessary
    IF (.NOT. ALLOCATED(Q)) THEN
       ALLOCATE(Q(nq_total, ndim))
    END IF
    IF (.NOT. ALLOCATED(WQ)) THEN
       ALLOCATE(WQ(nq_total))
    END IF
    IF (.NOT. ALLOCATED(P)) THEN
       ALLOCATE(P(nq_total, ndim))
    END IF
    IF (.NOT. ALLOCATED(WP)) THEN
       ALLOCATE(WP(nq_total))
    END IF
    IF (.NOT. ALLOCATED(Q0d)) THEN
       ALLOCATE(Q0d(ndim))
    END IF
  
    ! initialize weights and extrema
    WQ(:) = ONE
    x_min(:) = HUGE(ZERO)
    x_max(:) = -HUGE(ZERO)
    Q(:,:) = CZERO
    P(:,:) = CZERO
  
    ! 2) compute expected_total from per-dimension counts (sanity check)
    expected_total = 1
    DO inb = 1, ndim
       Nj = Basis%tab_basis(inb)%nq
       expected_total = expected_total * Nj
    END DO
  
    IF (expected_total /= nq_total) THEN
       WRITE(*,*) 'Warning: product of per-dim nq (', expected_total, ') /= Basis%nq (', nq_total, ')'
       ! proceed but be cautious
    END IF
  
    ! 3) Fill Q and WQ (tensor-product enumeration)
    CALL Init_tab_ind(Tab_iq, Basis%NDindexq)   ! prepare first multi-index according to your convention
    iq = 0
    Endloop = .FALSE.
    DO
       CALL increase_NDindex(Tab_iq, Basis%NDindexq, Endloop)
       IF (Endloop) EXIT
  
       iq = iq + 1
       IF (iq > nq_total) THEN
          WRITE(*,*) 'Build_NUFT_Grids: filled more points than nq_total; breaking.'
          EXIT
       END IF
  
       WQ(iq) = ONE
       DO inb = 1, ndim
          idx1d = Tab_iq(inb)    ! index along this dimension (1..N_j)
          Q(iq, inb) = Basis%tab_basis(inb)%x(idx1d)
          WQ(iq) = WQ(iq) * Basis%tab_basis(inb)%w(idx1d)
          x_min(inb) = MIN(x_min(inb), Q(iq, inb))
          x_max(inb) = MAX(x_max(inb), Q(iq, inb))
       END DO
  
    END DO
  
    ! ensure we've filled expected number
    IF (iq /= nq_total) THEN
       WRITE(*,*) 'Warning: Build_NUFT_Grids filled ', iq, ' points (expected ', nq_total, ').'
       ! Option: trim arrays or leave as is. We proceed using iq as filled count.
    END IF
  
    ! 4) compute extents L_j and delta_p_j (FFT-consistent choice)
    DO inb = 1, ndim
       L(inb) = x_max(inb) - x_min(inb)
       IF (L(inb) <= ZERO) THEN
          WRITE(*,*) 'Error: zero or negative spatial extent in dimension ', inb
          L(inb) = ONE
       END IF
       ! choose delta_p_j consistent with FFT mapping delta_p = 2*pi / L
       delta_p(inb) = (TWO * pi) / L(inb)
    END DO
  
    ! 5) Build centered regular P-grid per dimension
    !    We need per-dimension number of points N_j to determine center offset
    !    Use Basis%tab_basis(inb)%nq for N_j
    CALL Init_tab_ind(Tab_iq, Basis%NDindexq)   ! reset multi-index
    iq = 0
    Endloop = .FALSE.
    DO
       CALL increase_NDindex(Tab_iq, Basis%NDindexq, Endloop)
       IF (Endloop) EXIT
  
       iq = iq + 1
       IF (iq > nq_total) EXIT
  
       DO inb = 1, ndim
          Nj = Basis%tab_basis(inb)%nq
          center_offset = Nj / 2             ! integer division; for even Nj this is Nj/2
          idx1d = Tab_iq(inb)                 ! gives 1..Nj
          ! Build P component centered about zero:
          P(iq, inb) = delta_p(inb) * REAL( (idx1d - 1 - center_offset), Rkind )
       END DO
  
    END DO
  
    ! 6) compute uniform WP weights = product of per-dim delta_p
    prod_dp = ONE
    DO inb = 1, ndim
       prod_dp = prod_dp * delta_p(inb)
    END DO
    WP(:) = prod_dp
  
    ! 7) return Q0d: center of each Q-dimension (use Basis%tab_basis(inb)%Q0 
    DO inb = 1, ndim
       Q0d(inb) = Basis%tab_basis(inb)%Q0
    END DO
  
    ! 8) cleanup
    DEALLOCATE(Tab_iq)
    DEALLOCATE(L)
    DEALLOCATE(x_min)
    DEALLOCATE(x_max)
    DEALLOCATE(delta_p)
  
    RETURN
  END SUBROUTINE Build_NUFT_Grids


  SUBROUTINE Print_Grids(Basis)
   IMPLICIT NONE
   TYPE(Basis_t), INTENT(IN) :: Basis

   !< Locals ----------------------------------------------------
   REAL(KIND=Rkind), ALLOCATABLE :: Q(:,:)   ! nq_total x ndim
   REAL(KIND=Rkind), ALLOCATABLE :: WQ(:)     ! nq_total
   REAL(KIND=Rkind), ALLOCATABLE :: P(:,:)   ! nq_total x ndim
   REAL(KIND=Rkind), ALLOCATABLE :: WP(:)     ! nq_total (constant)
   REAL(KIND=Rkind), ALLOCATABLE :: Q0d(:)    ! ndim
   INTEGER :: ndim, nq, i, j

   !< Build Grids ----------------------------------------------------
   CALL Build_NUFT_Grids(Basis, Q, WQ, P, WP, Q0d)
 
   ndim = SIZE(Q, 2)
   nq   = SIZE(Q, 1)
 
   WRITE(out_unit,*) "==============================================="
   WRITE(out_unit,*) "            NUFT Grid Diagnostic               "
   WRITE(out_unit,*) "==============================================="
 
   ! print dimension and offsets
   WRITE(out_unit,*) "Number of dimensions: ", ndim
   WRITE(out_unit,*) "Total number of points: ", nq
   WRITE(out_unit,*) "Centers Q0d (per dimension):"
   WRITE(out_unit,*) (Q0d(j), j=1, ndim)
   WRITE(out_unit,*)
 
   ! print Q grid and WQ
   WRITE(out_unit,*) "Configuration-space grid Q and weights WQ:"
   DO i = 1, nq
      WRITE(out_unit,*) i, (Q(i,j), j=1, ndim), " | WQ=", WQ(i)
   END DO
   WRITE(out_unit,*)
 
   ! print P grid and WP
   WRITE(out_unit,*) "Momentum-space grid P and weights WP:"
   DO i = 1, nq
      WRITE(out_unit,*) i, (P(i,j), j=1, ndim), " | WP=", WP(i)
   END DO
   WRITE(out_unit,*) "==============================================="
   WRITE(out_unit,*)
 END SUBROUTINE Print_Grids

!>---------------------------------------------------------------------
!> @brief NUFFT: Non-uniform Fast Fourier Transform (approximation)
!> @details Computes the Fourier transform from a non-uniform position grid Q
!>          to a uniform momentum grid P, using multidimensional dot product for phase.
!---------------------------------------------------------------------

!>---------------------------------------------------------------------
!> @brief INV_NUFFT: Inverse Non-uniform Fast Fourier Transform
!> @details Transforms a wavefunction from uniform momentum grid P to non-uniform
!>          position grid Q. Uses multidimensional dot product.
!---------------------------------------------------------------------



!=======================================================================
! TEST_NUFFT_CYCLE: Test routine for NUFFT and INV_NUFFT
!
! Tests the forward and inverse non-uniform Fourier transforms
! using a provided wavefunction structure.
!=======================================================================
SUBROUTINE Test_NUFFT_Cycle(psi_x,Basis)
   IMPLICIT NONE

   ! Input
   TYPE(Basis_t), INTENT(IN)       :: Basis
   COMPLEX(KIND=Rkind), INTENT(IN) :: psi_x(:)

   ! Local variables
   REAL(KIND=Rkind), ALLOCATABLE :: Q(:,:), P(:,:), WQ(:), WP(:),Q0d(:)
   COMPLEX(KIND=Rkind), ALLOCATABLE :: psi_p(:), psi_x_reconstructed(:)
   REAL(KIND=Rkind), PARAMETER :: tolerance = 1.0E-6_Rkind
   REAL(KIND=Rkind) :: sum_diff_sq
   INTEGER :: nq_tot, ndim, iq
   REAL(KIND=Rkind) :: reconstruction_error
   LOGICAL :: test_passed

   ! 1. Build the non-uniform grids for position (Q) and uniform grids for momentum (P)
   CALL Build_NUFT_Grids(Basis, Q, WQ, P, WP,Q0d)

   nq_tot = SIZE(Q, 1)
   ndim   = SIZE(Q, 2)

   ! 2. Check if the input wavefunction has the correct size
   IF (SIZE(psi_x) /= nq_tot) THEN
       WRITE(*,*) "Error: The size of psi_x does not match the grid size."
       RETURN
   END IF

   ! 3. Perform the forward NUFFT
   ALLOCATE(psi_p(nq_tot))
   CALL NUFFT(Q, WQ, P, psi_x, psi_p)

   ! 4. Perform the inverse NUFFT
   ALLOCATE(psi_x_reconstructed(nq_tot))
   CALL INV_NUFFT(Q, WQ, P, WP, psi_p, psi_x_reconstructed)

   ! 5. Compare the reconstructed wavefunction with the initial one (psi_in%CVec)
   ! Calculate the Root Mean Square (RMS) error
   sum_diff_sq = ZERO
   DO iq = 1, nq_tot
       sum_diff_sq = sum_diff_sq + (ABS(psi_x(iq) - psi_x_reconstructed(iq)))**2
   END DO

   reconstruction_error = SQRT(sum_diff_sq / nq_tot)
   test_passed = reconstruction_error < tolerance

   ! 6. Produce a report
   WRITE(*,*) '===================================================='
   WRITE(*,*) 'NUFFT Cycle Test Report'
   WRITE(*,*) '===================================================='
   WRITE(*,*) 'Total grid points (nq_tot): ', nq_tot
   WRITE(*,*) 'Number of dimensions (ndim): ', ndim
   WRITE(*,*) 'Tolerance: ', tolerance
   WRITE(*,*) 'Reconstruction Error: ', reconstruction_error
   IF (test_passed) THEN
       WRITE(*,*) 'Test Result: PASSED'
       WRITE(*,*) 'The reconstruction is accurate within the tolerance.'
   ELSE
       WRITE(*,*) 'Test Result: FAILED'
       WRITE(*,*) 'The reconstruction error exceeds the tolerance.'
   END IF
   WRITE(*,*) '===================================================='



    ! 7. Write data to files for plotting (assuming ndim=2)
   IF (ndim == 2) THEN
      OPEN(UNIT=10, FILE='Q_and_psi_x.txt', STATUS='REPLACE')
      !WRITE(10, *) "x", "y", "psi_x_sq"
      DO iq = 1, nq_tot
          WRITE(10, *) Q(iq,1), Q(iq,2), ABS(psi_x(iq))**2
      END DO
      CLOSE(10)

      OPEN(UNIT=11, FILE='P_and_psi_p.txt', STATUS='REPLACE')
      !WRITE(11, *) "px", "py", "psi_p_sq"
      DO iq = 1, nq_tot
          WRITE(11, *) P(iq,1), P(iq,2), ABS(psi_p(iq))**2
      END DO
      CLOSE(11)
   END IF

   ! 7. Deallocate arrays
   DEALLOCATE(Q, P, WQ, WP)
   DEALLOCATE(psi_p, psi_x_reconstructed)

END SUBROUTINE Test_NUFFT_Cycle




   
    SUBROUTINE d2psi1D_cplx(d2GGB,GGB,Basis)
    IMPLICIT NONE

    COMPLEX(KIND=Rkind), INTENT(IN)                            :: GGB(:,:,:)
    COMPLEX(KIND=Rkind), INTENT(INOUT)                         :: d2GGB(:,:,:)
    TYPE(Basis_t),INTENT(IN)                                   :: Basis
    
    !locals variables---------------------------------------------------------
    LOGICAL                                                    :: debug
  
    INTEGER                                                    :: i1,i3
    
    
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING d2psi1D_cplx'
      flush (out_unit)
    END IF
    
      
      DO i3 = 1, ubound(GGB, dim=3)
         DO i1 = 1, ubound(GGB, dim=1)
           d2GGB(i1, :, i3) = d2GGB(i1, :, i3)+ matmul( Basis%d2gg(:, :, 1, 1) , GGB(i1,:,i3) )
         END DO
      END DO
    
    
    IF (debug) THEN
       flush (out_unit)
    END IF
    
   
 END SUBROUTINE

 SUBROUTINE Rdensity_alloc(Rdensity,Basis)
     
   USE  QDUtil_m  
   TYPE(REDUCED_DENSIRY_t),intent(inout)         :: Rdensity
   TYPE(Basis_t), intent(in), target             :: Basis
   integer                                       :: ib,ndim
   ndim  = SIZE(Basis%tab_basis) - 1
   If(allocated(Rdensity%tab_prob)) deallocate(Rdensity%tab_prob)
   allocate(Rdensity%tab_prob(ndim))
   IF(ndim==1) THEN
     allocate(Rdensity%prob(Basis%tab_basis(1)%nq))
     Rdensity%prob(:) = ZERO
    ! print*,'Rdensity%prob(:)',Rdensity%prob(:)
   ELSE
      DO ib=1,ndim       
        allocate(Rdensity%tab_prob(ib)%prob(Basis%tab_basis(ib)%nq))
        Rdensity%tab_prob(ib)%prob(:) =ZERO
        !print*,' ib,Rdensity%prob(:)',ib,Rdensity%tab_prob(ib)%prob(:) 
      END DO
   END IF        
 END SUBROUTINE



 SUBROUTINE Rdensity_Writing(Rdensity,Basis,nio,ib,t)
    
  USE  QDUtil_m  
  TYPE(REDUCED_DENSIRY_t),intent(inout)         :: Rdensity
  TYPE(Basis_t), intent(in), target             :: Basis
  integer  ,intent(in)                          :: nio
    integer  ,intent(in),optional               :: ib
  real(kind=Rkind),intent(in),optional          :: t 
  integer                                       :: ndim,nq,iq

  ndim  = SIZE(Basis%tab_basis) - 1
   If(present(ib)) then 
     nq= Basis%tab_basis(ib)%nq 
   Else
    nq= Basis%tab_basis(1)%nq 
   End If  

   DO Iq=1,nq
      If(ndim ==1) then
        If(present(t)) then
           write(nio,*) t, Basis%tab_basis(1)%x(Iq), Rdensity%prob(Iq)
        Else
           write(nio,*)  Basis%tab_basis(1)%x(Iq), Rdensity%prob(Iq)
        End If 
      Else 
        If(present(ib)) then        
          If(present(t)) then
            write(nio,*) t, Basis%tab_basis(ib)%X(Iq), Rdensity%tab_prob(ib)%prob(Iq)
          Else
            write(nio,*)  Basis%tab_basis(ib)%X(Iq), Rdensity%tab_prob(ib)%prob(Iq)
          End If 
          Else
          
           If(present(t)) then
             write(nio,*) t, Basis%tab_basis(1)%X(Iq),Basis%tab_basis(2)%X(Iq), &
             & Rdensity%tab_prob(1)%prob(Iq), Rdensity%tab_prob(2)%prob(Iq)
           Else
              write(nio,*) Basis%tab_basis(1)%X(Iq),Basis%tab_basis(2)%X(Iq), &
              & Rdensity%tab_prob(1)%prob(Iq), Rdensity%tab_prob(2)%prob(Iq)
           End If 
         End If
      End If    
   END DO
       
END SUBROUTINE


 SUBROUTINE Rdensity_Dealloc(Rdensity,Basis)
    
  USE  QDUtil_m  
  TYPE(REDUCED_DENSIRY_t),intent(inout)         :: Rdensity
  TYPE(Basis_t), intent(in), target             :: Basis
  integer                                       :: ib,ndim

  ndim  = SIZE(Basis%tab_basis) - 1

   IF (allocated(Rdensity%prob)) Deallocate(Rdensity%prob)

   IF (allocated(Rdensity%tab_prob)) then

   Do ib = 1,ndim
      IF (allocated(Rdensity%tab_prob(ib)%prob)) Deallocate(Rdensity%tab_prob(ib)%prob)
   End  Do

   End If

  IF (allocated(Rdensity%tab_prob)) Deallocate(Rdensity%tab_prob)
       
END SUBROUTINE



 SUBROUTINE Change_Basis_Parameters(Basis, Qt, SQt, At, Pt)
   IMPLICIT NONE
 
   ! === Input/Output ===
   TYPE(Basis_t), INTENT(INOUT)          :: Basis   ! Basis object to modify
   REAL(kind=Rkind), INTENT(IN)          :: Qt(:)    ! New centers (position)
   REAL(kind=Rkind), INTENT(IN)          :: SQt(:)   ! New scaling factors
   REAL(kind=Rkind), INTENT(IN)          :: Pt(:)    ! New impulses (momentum)
   COMPLEX(kind=Rkind), INTENT(IN)       :: At(:)    ! New complex width parameters
 
   ! === Local variables ===
   INTEGER                               :: ndim     ! Number of degrees of freedom
   INTEGER                               :: ib       ! Basis index
 
   ! === Compute number of spatial dimensions (assume last index is for electronic state)
   ndim = SIZE(Basis%tab_basis) - 1
 
   ! === Update basis parameters for each spatial dimension
   DO ib = 1, ndim
     Basis%tab_basis(ib)%Q0      = Qt(ib)
     Basis%tab_basis(ib)%SCALEQ  = SQt(ib)
     Basis%tab_basis(ib)%Imp_k   = Pt(ib)
     Basis%tab_basis(ib)%alpha   = At(ib)
   END DO
 
 END SUBROUTINE Change_Basis_Parameters

 SUBROUTINE Get_Basis_Parameters(Basis, Qt, SQt, At, Pt)
   IMPLICIT NONE
 
   ! === Input ===
   TYPE(Basis_t), INTENT(IN)           :: Basis      ! basis structure
 
   ! === Output ===
   REAL(kind=Rkind), INTENT(OUT)       :: Qt(:)      ! position centers
   REAL(kind=Rkind), INTENT(OUT)       :: SQt(:)     ! scaling (widths)
   REAL(kind=Rkind), INTENT(OUT)       :: Pt(:)      ! momentum components
   COMPLEX(kind=Rkind), INTENT(OUT)    :: At(:)      ! complex parameters alpha
 
   ! === Local variables ===
   INTEGER                             :: ndim       ! Number of spatial dimensions
   INTEGER                             :: ib         ! Loop index
 
   ! Determine the number of spatial dimensions (excluding electronic state)
   ndim = SIZE(Basis%tab_basis) - 1
 
   ! Extract the parameters from the basis
   DO ib = 1, ndim
     Qt(ib)  = Basis%tab_basis(ib)%Q0
     SQt(ib) = Basis%tab_basis(ib)%SCALEQ
     Pt(ib)  = Basis%tab_basis(ib)%Imp_k
     At(ib)  = Basis%tab_basis(ib)%alpha
   END DO
 
 END SUBROUTINE Get_Basis_Parameters


END MODULE Basis_m
