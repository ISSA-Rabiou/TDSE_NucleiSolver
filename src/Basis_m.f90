MODULE Basis_m
   USE QDUtil_m
   USE NDindex_m
   USE polyortho_m

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: Basis_t
   PUBLIC :: Read_Basis, Write_Basis
   PUBLIC :: Basis_IS_allocated, Basis_IS_allocatedtot
   PUBLIC :: Deallocate_Basis
   
   PUBLIC :: Calc_dngg_grid, Calc_dngg_grid_0
   PUBLIC :: Calc_tranpose_d0gb
   PUBLIC :: test_basitogridgridtobasis
   
   PUBLIC :: GridTOBasis_nD_cplx, BasisTOGrid_nD_cplx
   
   PUBLIC :: Calc_Q_grid, Calc_index
   PUBLIC :: Rdensity_alloc, Rdensity_Writing
   
   PUBLIC :: Scale_Basis, init_Basis1_TO_Basis2
   PUBLIC :: construct_primitive_basis, construct_primitive_basis_temp
   PUBLIC :: Construct_Basis_el,Construct_Basis_Fourier,Construct_Basis_Sin
   PUBLIC :: Construct_Basis_Ho,Construct_Basis_Hag
   PUBLIC :: CheckOrtho_Basis
   
   PUBLIC :: pdv2psi_nD
   
   PUBLIC :: Get_Basis_Parameters, Change_Basis_Parameters
   PUBLIC :: Calc_Basis_dPtQtAt
   
   PUBLIC :: Calc_reduced_Density_surf, Calc_reduced_density
   PUBLIC :: REDUCED_DENSIRY_t
   PUBLIC :: calc_tab_Iq0 ,Complete_construct_Basis
   PUBLIC :: Print_Grids, Build_NUFT_Grids, Test_NUFFT_Cycle
   PUBLIC :: ApplyExpVdtAndTransformToMomentum,ApplyExpTdtAndTransformToPosition
   

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
   INTEGER                             :: nb, nq
   REAL(KIND=Rkind)                    :: A, B
   REAL(KIND=Rkind)                    :: imp_k    ! Momentum parameter
   REAL(KIND=Rkind)                    :: scaleQ, Q0
   COMPLEX(KIND=Rkind)                 :: Alpha = CONE
   CHARACTER(LEN=:), ALLOCATABLE       :: Basis_name
   !<--- Grid points and weights
   REAL(KIND=Rkind),    ALLOCATABLE    :: x(:)      ! Grid points
   REAL(KIND=Rkind),    ALLOCATABLE    :: w(:)      ! Quadrature weights
   !<--- Derivative matrices
   COMPLEX(KIND=Rkind), ALLOCATABLE    :: d0gb(:, :)
   COMPLEX(KIND=Rkind), ALLOCATABLE    :: d0bgw(:, :)
   COMPLEX(KIND=Rkind), ALLOCATABLE    :: d1gb(:, :, :)
   COMPLEX(KIND=Rkind), ALLOCATABLE    :: d1gg(:, :, :)
   COMPLEX(KIND=Rkind), ALLOCATABLE    :: d2gb(:, :, :, :)
   COMPLEX(KIND=Rkind), ALLOCATABLE    :: d2gg(:, :, :, :)
   COMPLEX(KIND=Rkind), ALLOCATABLE    :: dagb(:, :)
   COMPLEX(KIND=Rkind), ALLOCATABLE    :: dagg(:, :)
   COMPLEX(KIND=Rkind), ALLOCATABLE    :: dq0gb(:, :)
   COMPLEX(KIND=Rkind), ALLOCATABLE    :: dq0gg(:, :)
   COMPLEX(KIND=Rkind), ALLOCATABLE    :: dp0gb(:, :)
   COMPLEX(KIND=Rkind), ALLOCATABLE    :: dp0gg(:, :)
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
      integer                          :: i

      !write(out_unit,*) '---------------------------------------------------------------------'
      !write(out_unit,*) 'Write_Basis'
      write(out_unit,*) "Basis_name=", Basis%Basis_name
      !write(out_unit,*) "n_basis=", Basis%nb_basis
      !write(out_unit,*) 'nb,nq',Basis%nb,Basis%nq
      !write(out_unit,*) "Q0,scaleQ ", Basis%Q0,Basis%scaleQ
      !write(out_unit,*)  'A,B',Basis%A,Basis%B
      !write(out_unit,*) 'nb_basis',Basis%nb_basis

      IF (.NOT. allocated(Basis%x)) THEN
         write(out_unit,*)' Basis table x is not allocated.'
      ELSE
         call  Write_VecMat(Basis%x, out_unit, 5,  info='x')
      END IF
      write (out_unit, *)
      IF (.NOT. allocated(Basis%W)) THEN
           write(out_unit,*)' Basis table w is not allocated.'
      ELSE
         call  Write_VecMat(Basis%w, out_unit, 5,  info='w')
      END IF
      write (out_unit, *)
       IF (.NOT.allocated(Basis%d0gb)) THEN
        write(out_unit,*)' Basis table d0gb is not allocated.'
      ELSE
      call   Write_VecMat(Basis%d0gb, out_unit, 5,  info='d0gb')
       END IF
          

       write (out_unit, *)
       IF (.NOT.allocated(Basis%dp0gb)) THEN
        write(out_unit,*)' Basis table dp0gb is not allocated.'
      ELSE
      call   Write_VecMat(Basis%dp0gb, out_unit, 5,  info='dp0gb')
       END IF
      
       
       write (out_unit, *)
       IF (.NOT.allocated(Basis%dq0gb)) THEN
        write(out_unit,*)' Basis table dq0gb is not allocated.'
      ELSE
      call   Write_VecMat(Basis%dq0gb, out_unit, 5,  info='dq0gb')
       END IF
       
       write (out_unit, *)
       IF (.NOT.allocated(Basis%dagb)) THEN
        write(out_unit,*)' Basis table dagb is not allocated.'
      ELSE
      call   Write_VecMat(Basis%dagb, out_unit, 5,  info='dagb')
       END IF

       write (out_unit, *)
       IF (.NOT.allocated(Basis%dagg)) THEN
        write(out_unit,*)' Basis table dagg is not allocated.'
      ELSE
      call   Write_VecMat(Basis%dagg, out_unit, 5,  info='dagg')
       END IF

       write (out_unit, *)
       IF (.NOT.allocated(Basis%dp0gg)) THEN
        write(out_unit,*)' Basis table dp0gg is not allocated.'
      ELSE
      call   Write_VecMat(Basis%dp0gg, out_unit, 5,  info='dp0gg')
       END IF
      
       write (out_unit, *)
       IF (.NOT.allocated(Basis%dq0gg)) THEN
        write(out_unit,*)' Basis table dq0gg is not allocated.'
      ELSE
      call   Write_VecMat(Basis%dq0gg, out_unit, 5,  info='dq0gg')
       END IF

      write(out_unit,*)
        IF (.NOT.allocated(Basis%d0bgw)) THEN
           write(out_unit,*)' Basis table d0bgw is not allocated.'
         ELSE
      call   Write_VecMat(Basis%d0bgw,out_unit,5, info='d0gbw')
        END IF

         write(out_unit,*)
         IF (.NOT.allocated(Basis%d1gb)) THEN
           write(out_unit,*)' Basis table d1gb is not allocated.'
         ELSE
      call   Write_VecMat(Basis%d1gb(:,:,1),out_unit,5, info='d1gb')
         END IF
         write(out_unit,*)
         IF (.NOT.allocated(Basis%d1gg)) THEN
           write(out_unit,*)' Basis table d1gb is not allocated.'
         ELSE
      call   Write_VecMat(Basis%d1gg(:,:,1),out_unit,5, info='d1gg')
         END IF

         write(out_unit,*)
         IF (.NOT.allocated(Basis%d2gb)) THEN
           write(out_unit,*)' Basis table d1gb is not allocated.'
         ELSE
      call   Write_VecMat(Basis%d2gb(:,:,1,1),out_unit,5, info='d2gb')
         END IF

         write(out_unit,*)
         IF (.NOT.allocated(Basis%d2gg)) THEN
           write(out_unit,*)' Basis table d2gg is not allocated.'
         ELSE
      call   Write_VecMat(Basis%d2gg(:,:,1,1),out_unit,5, info='d2gg')
         END IF
          

         write(out_unit,*)
         IF (.NOT.allocated(Basis%S)) THEN
            write(out_unit,*)' Basis table S is not allocated.'
         ELSE
      !call  Write_VecMat(Basis%S,out_unit,5, info='S')
         END IF  
         IF (allocated(Basis%tab_basis)) THEN
            DO i = 1, size(Basis%tab_basis)
              CALL Write_Basis(Basis%tab_basis(i))
            END DO
         END IF
      write (out_unit, *) '--------------------------------------------------------------------------'

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
        Basis2%nb          = Basis1%nb
        Basis2%nq          = Basis1%nq
        Basis2%Q0          = Basis1%Q0
        Basis2%SCALEQ      = Basis1%SCALEQ
        Basis2%Imp_k       = Basis1%Imp_k
        Basis2%Alpha       = Basis1%Alpha
        Basis2%A           = Basis1%A
        Basis2%B           = Basis1%B
        Basis2%S           = Basis1%S
     END IF
     
     END SUBROUTINE init_Basis1_TO_Basis2

   !> Recursive subroutine to safely deallocate all allocated components
   !! of a Basis_t object, including nested basis structures.
   RECURSIVE SUBROUTINE Deallocate_Basis(Basis)
   USE QDUtil_m
   IMPLICIT NONE
   
   TYPE(Basis_t), intent(inout) :: Basis
   integer                      :: i
   
   ! Deallocate ND indexing arrays if present
   IF (allocated(Basis%tab_basis)) THEN
      IF (allocated(Basis%NDindexq%Tab0)) deallocate(Basis%NDindexq%Tab0)
      IF (allocated(Basis%NDindexb%Tab0)) deallocate(Basis%NDindexb%Tab0)
      deallocate(Basis%tab_basis)
   END IF
   
   ! Deallocate 1D basis arrays if present
   IF (allocated(Basis%x)) THEN
      deallocate(Basis%x)
      write(out_unit,*) 'Basis table x is now deallocated.'
   END IF
   
   IF (allocated(Basis%d0gb)) THEN
      deallocate(Basis%d0gb)
      write(out_unit,*) 'Basis table d0gb is now deallocated.'
   END IF
   
   IF (allocated(Basis%d0bgw)) THEN
      deallocate(Basis%d0bgw)
      write(out_unit,*) 'Basis table d0bgw is now deallocated.'
   END IF
   
   IF (allocated(Basis%dagb)) THEN
      deallocate(Basis%dagb)
      write(out_unit,*) 'Basis table dagb is now deallocated.'
   END IF
   
   IF (allocated(Basis%dagg)) THEN
      deallocate(Basis%dagg)
      write(out_unit,*) 'Basis table dagg is now deallocated.'
   END IF
   
   IF (allocated(Basis%dq0gb)) THEN
      deallocate(Basis%dq0gb)
      write(out_unit,*) 'Basis table dq0gb is now deallocated.'
   END IF
   
   IF (allocated(Basis%dq0gg)) THEN
      deallocate(Basis%dq0gg)
      write(out_unit,*) 'Basis table dq0gg is now deallocated.'
   END IF
   
   IF (allocated(Basis%dp0gg)) THEN
      deallocate(Basis%dp0gg)
      write(out_unit,*) 'Basis table dp0gg is now deallocated.'
   END IF
   
   IF (allocated(Basis%d1gb)) THEN
      deallocate(Basis%d1gb)
      write(out_unit,*) 'Basis table d1gb is now deallocated.'
   END IF
   
   IF (allocated(Basis%d1gg)) THEN
      deallocate(Basis%d1gg)
      write(out_unit,*) 'Basis table d1gg is now deallocated.'
   END IF
   
   IF (allocated(Basis%d2gb)) THEN
      deallocate(Basis%d2gb)
      write(out_unit,*) 'Basis table d2gb is now deallocated.'
   END IF
   
   IF (allocated(Basis%d2gg)) THEN
      deallocate(Basis%d2gg)
      write(out_unit,*) 'Basis table d2gg is now deallocated.'
   END IF
   
   ! Recursively deallocate sub-basis if composite
   IF (allocated(Basis%tab_basis)) THEN
      DO i = 1, size(Basis%tab_basis)
         CALL Deallocate_Basis(Basis%tab_basis(i))
      END DO
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
   LOGICAL, PARAMETER                  :: debug = .TRUE.
   TYPE(Basis_t), INTENT(INOUT)        :: Basis
   INTEGER, INTENT(IN)                 :: nio

   ! Local variables
   INTEGER                              :: err_io, nb, nq, i, nb_basis, ib
   CHARACTER(LEN=Name_len)              :: name
   REAL(KIND=Rkind)                     :: A, B, scaleQ, Q0, Imp_k
   COMPLEX(KIND=Rkind)                  :: Alpha

   ! Namelist declaration
   NAMELIST /basis_nD/ name, nb_basis, nb, nq, A, B, scaleQ, Q0, Imp_k, Alpha

   ! ===== Default values =====
   nb_basis = 0
   nb       = 0
   nq       = 0
   A        = ZERO
   B        = ZERO
   Q0       = ZERO
   scaleQ   = ONE
   name     = '0'
   Imp_k    = ZERO
   Alpha    = CONE

   ! ===== Read namelist =====
   READ (nio, NML=basis_nD, IOSTAT=err_io)

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
   IF (nb_basis > 1) THEN
      Basis%Basis_name = 'Dp'
      Basis%nb_basis   = nb_basis
      CALL string_uppercase_TO_lowercase(Basis%Basis_name)

      ALLOCATE (Basis%tab_basis(nb_basis))
      DO i = 1, nb_basis
         CALL Read_Basis(Basis%tab_basis(i), nio)
      END DO

      ! Compute total nb and nq excluding "el" bases
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
      Basis%nb         = nb
      Basis%nq         = nq
      Basis%Q0         = Q0
      Basis%scaleQ     = scaleQ
      Basis%A          = A
      Basis%B          = B
      Basis%Imp_k      = Imp_k
      Basis%Alpha      = Alpha
      Basis%Basis_name = TRIM(ADJUSTL(name))

      ! Initialize overlap matrix S as identity
      ALLOCATE (Basis%S(nb, nb))
      Basis%S(:, :) = ZERO
      DO ib = 1, Basis%nb
         Basis%S(ib, ib) = CONE
      END DO

      CALL string_uppercase_TO_lowercase(Basis%Basis_name)
   END IF

END SUBROUTINE Read_Basis

   
   !> Calculates the index table Tab_iq for multi-dimensional quadrature indexing
   !! Tab_iq(:, iq) contains the index vector for the iq-th point in the total basis grid.
   SUBROUTINE Calc_tab_Iq0(Tab_Iq, Basis)
      USE QDUtil_m
      ! Purpose: Calculate the table of quantum indices (Tab_Iq) for a given basis
      
      ! Input/Output variables:
      TYPE(Basis_t), intent(in), target :: Basis       ! Basis set information
      integer, allocatable, intent(inout) :: Tab_iq(:,:)  ! Output table of quantum indices
      
      ! Local variables:
      integer, allocatable :: Tab_iq0(:)    ! Temporary array for current quantum indices
      integer :: ndim                      ! Number of dimensions
      integer :: iq                        ! Current quantum state index
      integer :: nq                        ! Total number of quantum states
      logical :: Endloop                   ! Flag to indicate end of loop
      
      ! Get dimensions from the basis set:
      ndim = size(Basis%tab_basis) - 1     ! Number of dimensions (subtracting 1 for 0-based or other reason)
      nq = Basis%nq                        ! Total number of quantum states
      
      ! Allocate memory for the output table and temporary array:
      allocate(Tab_iq(ndim, nq), Tab_iq0(ndim))
      
      ! Initialize the temporary index array:
      Call Init_tab_ind(Tab_iq0, Basis%NDindexq)
      
      ! Loop through all quantum states:
      Iq = 0
      DO
          Iq = Iq + 1  ! Increment quantum state counter
          
          ! Get the next set of quantum numbers:
          CALL increase_NDindex(Tab_Iq0, Basis%NDindexq, Endloop)
          
          ! Exit if we've processed all quantum states:
          IF (Endloop) exit
          
          ! Store the current quantum numbers in the output table:
          Tab_iq(:, Iq) = Tab_Iq0
      END DO
      
      ! Clean up temporary array:
      deallocate(Tab_Iq0)
  END SUBROUTINE Calc_tab_Iq0

   RECURSIVE SUBROUTINE construct_primitive_basis_temp(Basis)
   USE QDUtil_m
   IMPLICIT NONE

   !==================== Parameters =====================
   LOGICAL, PARAMETER :: debug = .true.
   !LOGICAL, PARAMETER :: debug = .false.  ! Uncomment to disable debug mode

   !==================== Interface ======================
   TYPE(Basis_t), INTENT(INOUT) :: Basis   ! Basis structure to construct

   !==================== Local variables ================
   INTEGER :: ndim     ! Number of sub-bases
   INTEGER :: ib       ! Index for loop over sub-bases

   !==================== Logic ==========================
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
      logical, parameter                     :: debug = .true.
      !logical,             parameter        ::debug = .false.
      TYPE(Basis_t), intent(inout)           :: Basis
      integer, allocatable                   :: NDend_q(:)
      integer, allocatable                   :: NDend_b(:)
      integer                                :: ndim, ib

       ndim = Basis%nb_basis - 1
       allocate (NDend_q(ndim))
       allocate (NDend_b(ndim)) 
       DO ib = 1, ndim
          NDend_q(ib) = Basis%tab_basis(ib)%nq
          NDend_b(ib) = Basis%tab_basis(ib)%nb
       END DO 
       call Init_NDindex(Basis%NDindexq, NDend_q, ndim)
       call Init_NDindex(Basis%NDindexb, NDend_b, ndim)

       call construct_primitive_basis_temp(Basis)

        
   END SUBROUTINE 

    SUBROUTINE Construct_Basis_Sin(Basis) ! sin : boxAB with A=0 and B=pi
      USE QDUtil_m

      TYPE(Basis_t), intent(inout)        :: Basis
      real(kind=Rkind)                    :: dx
      integer                             :: ib, iq, nb, nq


      if (allocated(Basis%x)) deallocate (Basis%x)
      if (allocated(Basis%w)) deallocate (Basis%w)
      if (allocated(Basis%d0gb)) deallocate (Basis%d0gb)
      if (allocated(Basis%d1gb)) deallocate (Basis%d1gb)
      if (allocated(Basis%d2gb)) deallocate (Basis%d2gb)
      if (allocated(Basis%dagb)) deallocate (Basis%dagb)
      if (allocated(Basis%dp0gb)) deallocate (Basis%dp0gb)
      if (allocated(Basis%dq0gb)) deallocate (Basis%dq0gb)

      nb = Basis%nb
      nq = Basis%nq
      dx = pi/nq

      ! grid and weight
      Basis%x = [(dx*(iq - HALF), iq=1, nq)]
      Basis%w = [(dx, iq=1, nq)]

      allocate (Basis%d0gb(nq, nb))
      allocate (Basis%d1gb(nq, nb, 1))
      allocate (Basis%d2gb(nq, nb, 1, 1))


      DO ib = 1, nb
         DO iq = 1, nq
            call d0d1d2box(Basis%x(iq),Basis%d0gb(iq, ib),Basis%d1gb(iq, ib, 1)&
                           &,Basis%d2gb(iq, ib, 1, 1),ib)
         END DO
      END DO

      IF (nb == nq) THEN
         Basis%d0gb(:, nb) = Basis%d0gb(:, nb)/sqrt(TWO)
         Basis%d1gb(:, nb, :) = Basis%d1gb(:, nb, :)/sqrt(TWO)
         Basis%d2gb(:, nb, :, :) = Basis%d2gb(:, nb, :, :)/sqrt(TWO)
      END IF

   END SUBROUTINE Construct_Basis_Sin




   !> Construct the Fourier basis in the interval [-π, π]
   !! This routine allocates the quadrature grid and computes the Fourier
   !! basis functions and their derivatives.
   SUBROUTINE Construct_Basis_Fourier(Basis)
      USE QDUtil_m
      IMPLICIT NONE
   
      TYPE(Basis_t), intent(inout) :: Basis
      real(kind=Rkind)             :: dx
      integer                      :: ib, iq, nb, nq
   
      !---------------------------------------------------------------------------
      ! Deallocate existing arrays if already allocated
      IF (allocated(Basis%x))       deallocate(Basis%x)
      IF (allocated(Basis%w))       deallocate(Basis%w)
      IF (allocated(Basis%d0gb))    deallocate(Basis%d0gb)
      IF (allocated(Basis%d1gb))    deallocate(Basis%d1gb)
      IF (allocated(Basis%d2gb))    deallocate(Basis%d2gb)
      IF (allocated(Basis%dagb))    deallocate(Basis%dagb)
      IF (allocated(Basis%dp0gb))   deallocate(Basis%dp0gb)
      IF (allocated(Basis%dq0gb))   deallocate(Basis%dq0gb)
      !---------------------------------------------------------------------------
   
      nb = Basis%nb
      nq = Basis%nq
   
      ! Step size for the uniform grid over [-π, π]
      dx = TWO * pi / nq
   
      !---------------------------------------------------------------------------
      ! Allocate and initialize grid points and weights
      allocate(Basis%x(nq))
      allocate(Basis%w(nq))
   
      Basis%x = [(iq * dx - dx / 2 - pi, iq = 1, nq)]
      Basis%w = [(dx, iq = 1, nq)]
      !---------------------------------------------------------------------------
   
      !---------------------------------------------------------------------------
      ! Allocate derivative matrices
      allocate(Basis%d0gb(nq, nb))
      allocate(Basis%d1gb(nq, nb, 1))
      allocate(Basis%d2gb(nq, nb, 1, 1))
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
    
      TYPE(Basis_t), INTENT(inout) :: Basis
    
      ! Set number of grid points to zero for 'el' basis
      Basis%nq = 0
    
      RETURN
    END SUBROUTINE Construct_Basis_el

   !> @brief Constructs the 'ho' (harmonic oscillator) type basis
   !! This subroutine initializes the basis functions and corresponding
   !! grid values for the harmonic oscillator using Hermite polynomials.
   !! It allocates and computes the grid points, weights, and the 0th, 1st, and 2nd
   !! derivatives of the Hermite polynomials at each grid point.
   SUBROUTINE Construct_Basis_Ho(Basis)
      USE QDUtil_m
      IMPLICIT NONE
    
      TYPE(Basis_t), INTENT(inout) :: Basis
      INTEGER                      :: iq, ib
      INTEGER                      :: nb, nq
    
      ! Retrieve the number of basis functions and grid points
      nb = Basis%nb
      nq = Basis%nq
    
      ! Deallocate previously allocated arrays if necessary
      IF (ALLOCATED(Basis%x))       DEALLOCATE(Basis%x)
      IF (ALLOCATED(Basis%w))       DEALLOCATE(Basis%w)
      IF (ALLOCATED(Basis%d0gb))    DEALLOCATE(Basis%d0gb)
      IF (ALLOCATED(Basis%d1gb))    DEALLOCATE(Basis%d1gb)
      IF (ALLOCATED(Basis%d2gb))    DEALLOCATE(Basis%d2gb)
      IF (ALLOCATED(Basis%dagb))    DEALLOCATE(Basis%dagb)
      IF (ALLOCATED(Basis%dp0gb))   DEALLOCATE(Basis%dp0gb)
      IF (ALLOCATED(Basis%dq0gb))   DEALLOCATE(Basis%dq0gb)
    
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
    
    END SUBROUTINE Construct_Basis_Ho


   SUBROUTINE Construct_Basis_Hag(Basis)
      USE QDUtil_m
      TYPE(Basis_t), intent(inout) :: Basis
  
      integer :: iq, ib
      integer :: nb, nq
      real(kind=Rkind) :: q0, sq, p0, a, x
      complex(kind=Rkind) :: alpha
  
      !----------------------------------------------------------
      ! Retrieve basis parameters from the Basis structure
      !----------------------------------------------------------
      q0    = Basis%Q0
      sq    = Basis%scaleQ
      alpha = Basis%alpha
      p0    = Basis%imp_k
  
      nb = Basis%nb
      nq = Basis%nq
  
      a = REAL(alpha, kind=Rkind)  ! Re(alpha)
  
      !----------------------------------------------------------
      ! Deallocate previous basis arrays if allocated
      !----------------------------------------------------------
      if (allocated(Basis%x)) deallocate (Basis%x)
      if (allocated(Basis%w)) deallocate (Basis%w)
      if (allocated(Basis%d0gb)) deallocate (Basis%d0gb)
      if (allocated(Basis%d1gb)) deallocate (Basis%d1gb)
      if (allocated(Basis%d2gb)) deallocate (Basis%d2gb)
  
      !----------------------------------------------------------
      ! Allocate grids and quadrature weights (Hermite quadrature)
      !----------------------------------------------------------
      allocate (Basis%x(nq))
      allocate (Basis%w(nq))
      call hercom(nq, Basis%x(:), Basis%w(:))
  
      !----------------------------------------------------------
      ! Allocate basis function arrays
      !----------------------------------------------------------
      allocate (Basis%d0gb(nq, nb))
      allocate (Basis%d1gb(nq, nb, 1))
      allocate (Basis%d2gb(nq, nb, 1, 1))
  
      !----------------------------------------------------------
      ! Construct Hagedorn basis functions directly
      !----------------------------------------------------------
      DO iq = 1, nq
          DO ib = 1, nb
              call d0d1d2poly_Hag_exp(Basis%x(iq), ib - 1, &
                   Basis%d0gb(iq, ib), Basis%d1gb(iq, ib, 1), Basis%d2gb(iq, ib, 1, 1), &
                   .true., alpha, p0)
          END DO
      END DO
  
  END SUBROUTINE Construct_Basis_Hag


   !> Check orthonormality of the basis set
   !! This subroutine checks if the basis functions are orthonormal (⟨b_i | b_j⟩ ≈ δ_ij)
   !! It also checks the overlaps with derivatives if requested (up to second order).

  SUBROUTINE CheckOrtho_Basis(Basis, nderiv)
   USE QDUtil_m
   IMPLICIT NONE

   TYPE(Basis_t), intent(in)          :: Basis
   INTEGER, intent(in)                :: nderiv
   INTEGER                            :: ib, jb
   COMPLEX(KIND=Rkind), ALLOCATABLE      :: S(:,:)
   REAL(KIND=Rkind)                   :: Sii, Sij

   !-------------------------------------------------------------------------
   ! Check for incompatible basis type
   IF (Basis%Basis_name == 'el') THEN
      WRITE(out_unit, '(A)') 'This routine is not applicable to Basis ''el'''
      RETURN
   END IF
   !-------------------------------------------------------------------------

   ! Proceed only if basis data is allocated
   IF (Basis_IS_allocated(Basis)) THEN

      WRITE(out_unit,*) 
     ! WRITE(out_unit,'(A)') '--------------------------------------------------'
      !WRITE(out_unit,'(A,A)') ' Ortho Basis: ', TRIM(Basis%Basis_name)
      WRITE(out_unit,'(A)') '--------------------------------------------------'

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
      IF (nderiv > 0) THEN
         WRITE(out_unit,*) 
         !WRITE(out_unit,'(A)') 'First derivative overlaps <d0b|d1b>:'
         !S = matmul(conjg(Basis%d0bgw), Basis%d1gb(:, :, 1))
         ! Optionally, format and write the matrix here
         ! CALL Write_VecMat(S, out_unit, 5, info='<d0b|d1b>')
      END IF

      !-------------------------------------------------------------------------
      ! Check overlap with second derivatives if required
      IF (nderiv > 1) THEN
         WRITE(out_unit,*) 
         !WRITE(out_unit,'(A)') 'Second derivative overlaps <d0b|d2b>:'
         S = MATMUL(CONJG(Basis%d0bgw), Basis%d2gb(:, :, 1, 1))
         ! Optionally, format and write the matrix here
         ! CALL Write_VecMat(S, out_unit, 5, info='<d0b|d2b>')
      END IF

   ELSE
      WRITE(out_unit,'(A,A)') 'WARNING in CheckOrtho_Basis: The basis is not allocated. Basis=', TRIM(Basis%Basis_name)
   END IF

END SUBROUTINE CheckOrtho_Basis



  SUBROUTINE Scale_Basis(Basis, q0, sq)
      USE QDUtil_m
      TYPE(Basis_t), intent(inout)  :: Basis
      real(kind=Rkind), intent(in)  :: q0, sq
   
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
   
      REAL(KIND=Rkind) :: sq_half, sq_3half, sq_5half
   
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
         write (out_unit, *) ' ERROR in Scale_Basis'
         write (out_unit, *) ' sq is too small or the basis is not allocated.'
         STOP 'ERROR in Scale_Basis'
      END IF
   END SUBROUTINE Scale_Basis

   SUBROUTINE Complete_construct_Basis(Basis)
      USE QDUtil_m
      TYPE(Basis_t), intent(inout) :: Basis
      real(kind=Rkind) :: Q0, SQ0, a
      complex(kind=Rkind) :: alpha
  
      !--------------------------------------------------------------
      ! Retrieve alpha and compute scaling factor SQ0 = sqrt(Re(alpha))
      !--------------------------------------------------------------
      alpha = Basis%alpha
      a     = REAL(alpha, kind=Rkind)
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
      call Calc_dngg_grid_0(Basis)
      !--------------------------------------------------------------
      ! Check orthonormality of the constructed basis up to 2nd derivative
      !--------------------------------------------------------------
      call CheckOrtho_Basis(Basis, nderiv=2)
  
  END SUBROUTINE 
 

    !> Compute the transpose of d0gb and apply the quadrature weights
    !! This subroutine computes the weighted transpose of the d0gb 
    !! matrix and stores it in d0bgw.
    SUBROUTINE Calc_tranpose_d0gb(Basis)
       USE QDUtil_m
       IMPLICIT NONE
    
       TYPE(Basis_t), intent(inout)      :: Basis
       INTEGER                           :: ib, nb
    
       !-------------------------------------------------------------------------
       ! Get the number of basis functions
       nb = Basis%nb
    
       !-------------------------------------------------------------------------
       ! If the output array already exists, deallocate it
       IF (allocated(Basis%d0bgw)) DEALLOCATE(Basis%d0bgw)
    
       !-------------------------------------------------------------------------
       ! Special case: skip computation if Basis is of type 'el'
       IF (Basis%Basis_name == 'el') RETURN
    
       !-------------------------------------------------------------------------
       ! Compute the transpose and apply weights
       Basis%d0bgw = TRANSPOSE(Basis%d0gb)
    
       DO ib = 1, nb
          Basis%d0bgw(ib, :) = Basis%d0bgw(ib, :) * Basis%w(:)
       END DO
    
    END SUBROUTINE Calc_tranpose_d0gb


    SUBROUTINE Calc_Basis_dPtQtAt(Basis)

    USE QDUtil_m
    TYPE(Basis_t), intent(inout)    :: Basis
    integer                         :: ib,iq,nb,nq
    !logical, parameter             :: debug = .true.
    logical, parameter              :: debug = .false.
    real(kind=Rkind)                :: a,Q0 


    Q0=Basis%Q0
    a = Basis%SCALEQ
    nb = Basis%nb
    nq = Basis%nq

    IF (debug) THEN
       write (out_unit, *) 'BEGINNING Calc_dngg_grid'
       CALL Write_Basis(Basis)
       flush (out_unit)
    END IF

    if (allocated(Basis%dagb)) deallocate (Basis%dagb)
    if (allocated(Basis%dq0gb)) deallocate (Basis%dq0gb)
    if (allocated(Basis%dp0gb)) deallocate (Basis%dp0gb)
     
    allocate (Basis%dagb(nq, nb))
    allocate (Basis%dq0gb(nq, nb))
    allocate (Basis%dp0gb(nq, nb))

   if (Basis%Basis_name == 'herm' .or. Basis%Basis_name == 'ho') then

    do ib = 1,nb
       do iq = 1,nq
          Basis%dagb(iq, ib)  = (ONE/a)*(HALF*Basis%d0gb(iq, ib)+(Basis%x(iq)-Q0)*Basis%d1gb(iq, ib, 1))
          Basis%dq0gb(iq, ib) = -Basis%d1gb(iq, ib, 1)
          Basis%dp0gb(iq, ib) = EYE*(Basis%x(iq)-Q0)*Basis%d0gb(iq, ib)
       end do
    end do

    Elseif(Basis%Basis_name == 'boxab') then
   do ib = 1,nb
     do iq = 1,nq
       Basis%dagb(iq, ib)  = (ONE/a)*(HALF*Basis%d0gb(iq, ib)+(Basis%x(iq)-Q0)*Basis%d1gb(iq, ib, 1))
       Basis%dq0gb(iq, ib) = -Basis%d1gb(iq, ib, 1)
       Basis%dp0gb(iq, ib) =CZERO
     end do
   end do
     Elseif(Basis%Basis_name == 'fourier') then

     Basis%dagb(:, :)   = CZERO
     Basis%dq0gb(:,:) = CZERO
     Basis%dp0gb(:,:) = CZERO

     End If

 END SUBROUTINE

    !> Computes several grid-based derivative matrices in the gg (grid-grid) representation
   !! using the precomputed gb (grid-basis) data and the weighted conjugate transpose (d0bgw).
   SUBROUTINE Calc_dngg_grid(Basis)
      USE QDUtil_m
      IMPLICIT NONE
   
      TYPE(Basis_t), intent(inout) :: Basis
      INTEGER                      :: ib, iq, nb, nq
      LOGICAL, PARAMETER           :: debug = .false.  ! Set to true for debugging output
   
      !-------------------------------------------------------------------------
      ! Retrieve number of basis functions and grid points
      nb = Basis%nb
      nq = Basis%nq
   
      !-------------------------------------------------------------------------
      ! Debug output (if enabled)
      IF (debug) THEN
         WRITE(out_unit, *) 'BEGINNING Calc_dngg_grid'
         CALL Write_Basis(Basis)
         FLUSH(out_unit)
      END IF
   
      !-------------------------------------------------------------------------
      ! Deallocate previous gg matrices if already allocated
      IF (allocated(Basis%d1gg)) DEALLOCATE(Basis%d1gg)
      IF (allocated(Basis%d2gg)) DEALLOCATE(Basis%d2gg)
      IF (allocated(Basis%dagg)) DEALLOCATE(Basis%dagg)
      IF (allocated(Basis%dq0gg)) DEALLOCATE(Basis%dq0gg)
      IF (allocated(Basis%dp0gg)) DEALLOCATE(Basis%dp0gg)
   
      !-------------------------------------------------------------------------
      ! Allocate matrices for gg representations
      ALLOCATE(Basis%d1gg(nq, nq, 1))
      ALLOCATE(Basis%d2gg(nq, nq, 1, 1))
      ALLOCATE(Basis%dagg(nq, nq))
      ALLOCATE(Basis%dq0gg(nq, nq))
      ALLOCATE(Basis%dp0gg(nq, nq))
   
      !-------------------------------------------------------------------------
      ! Optional debug output: show d0bgw matrix
      IF (debug) THEN
         CALL Write_VecMat(Basis%d0bgw, out_unit, 5, info='d0bgw')
         WRITE(out_unit, *)
      END IF
   
      !-------------------------------------------------------------------------
      ! Skip computation for 'el' basis type
      IF (Basis%Basis_name == 'el') RETURN
   
      !-------------------------------------------------------------------------
      ! Compute grid-grid matrices via matrix multiplications
      Basis%d1gg(:, :, 1)     = MATMUL(Basis%d1gb(:, :, 1), CONJG(Basis%d0bgw))
      Basis%d2gg(:, :, 1, 1)  = MATMUL(Basis%d2gb(:, :, 1, 1), CONJG(Basis%d0bgw))
      Basis%dagg              = MATMUL(Basis%dagb, CONJG(Basis%d0bgw))
      Basis%dq0gg             = MATMUL(Basis%dq0gb, CONJG(Basis%d0bgw))
      Basis%dp0gg             = MATMUL(Basis%dp0gb, CONJG(Basis%d0bgw))
   
      !-------------------------------------------------------------------------
      ! Debug output: print calculated matrices
      IF (debug) THEN
         CALL Write_VecMat(Basis%d1gg(:, :, 1), out_unit, 5, info='d1gg')
         WRITE(out_unit, *)
         CALL Write_VecMat(Basis%d2gg(:, :, 1, 1), out_unit, 5, info='d2gg')
         WRITE(out_unit, *)
         CALL Write_VecMat(Basis%dp0gg, out_unit, 5, info='dp0gg')
         WRITE(out_unit, *)
         CALL Write_VecMat(Basis%dq0gg, out_unit, 5, info='dq0gg')
         WRITE(out_unit, *)
         CALL Write_VecMat(Basis%dagg, out_unit, 5, info='dagg')
         WRITE(out_unit, *) 'END Calc_dngg_grid'
         FLUSH(out_unit)
      END IF
   
   END SUBROUTINE Calc_dngg_grid
   
   
   !> Computes only the first and second derivative operators (d1gg and d2gg)
   !! in the grid-grid representation, using the basis-grid derivative data.
   SUBROUTINE Calc_dngg_grid_0(Basis)
      USE QDUtil_m
      IMPLICIT NONE
   
      TYPE(Basis_t), intent(inout) :: Basis
      INTEGER                      :: ib, iq, nb, nq
      LOGICAL, PARAMETER           :: debug = .false.  ! Set to true for debug printing
   
      !-------------------------------------------------------------------------
      ! Retrieve the number of basis functions and grid points
      nb = Basis%nb
      nq = Basis%nq
   
      !-------------------------------------------------------------------------
      ! Optional debug print of basis setup
      IF (debug) THEN
         WRITE(out_unit, *) 'BEGINNING Calc_dngg_grid_0'
         CALL Write_Basis(Basis)
         FLUSH(out_unit)
      END IF
   
      !-------------------------------------------------------------------------
      ! Deallocate previous versions of d1gg and d2gg if already allocated
      IF (ALLOCATED(Basis%d1gg)) DEALLOCATE(Basis%d1gg)
      IF (ALLOCATED(Basis%d2gg)) DEALLOCATE(Basis%d2gg)
   
      !-------------------------------------------------------------------------
      ! Allocate arrays for the derivative operators in the grid representation
      ALLOCATE(Basis%d1gg(nq, nq, 1))
      ALLOCATE(Basis%d2gg(nq, nq, 1, 1))
   
      !-------------------------------------------------------------------------
      ! Optional debug print of basis-weighted transpose
      IF (debug) THEN
         CALL Write_VecMat(Basis%d0bgw, out_unit, 5, info='d0bgw')
         WRITE(out_unit, *)
      END IF
   
      !-------------------------------------------------------------------------
      ! Skip computation if basis is of type 'el'
      IF (Basis%Basis_name == 'el') RETURN
   
      !-------------------------------------------------------------------------
      ! Compute derivative matrices in grid-grid representation
      Basis%d1gg(:, :, 1) = MATMUL(Basis%d1gb(:, :, 1), CONJG(Basis%d0bgw))
      Basis%d2gg(:, :, 1, 1) = MATMUL(Basis%d2gb(:, :, 1, 1), CONJG(Basis%d0bgw))
   
      !-------------------------------------------------------------------------
      ! Debug output of resulting matrices
      IF (debug) THEN
         CALL Write_VecMat(Basis%d1gg(:, :, 1), out_unit, 5, info='d1gg')
         WRITE(out_unit, *)
         CALL Write_VecMat(Basis%d2gg(:, :, 1, 1), out_unit, 5, info='d2gg')
         WRITE(out_unit, *)
         WRITE(out_unit, *) 'END Calc_dngg_grid_0'
         FLUSH(out_unit)
      END IF
   
   END SUBROUTINE Calc_dngg_grid_0

   SUBROUTINE test_basitogridgridtobasis(Basis)
      USE QDUtil_m
      TYPE(Basis_t), intent(in)       :: Basis
      logical, parameter                 :: debug = .true.
      complex(kind=Rkind), allocatable   :: G1(:), B1(:)!,G(:), Hpsi(:)
      complex(kind=Rkind), allocatable   :: G2(:), B2(:)
      complex(kind=Rkind), allocatable   :: B(:)
      real(kind=Rkind), allocatable       :: diff_g(:), diff_b(:)
      !REAL(KIND=Rkind)                   :: Norm0,Norm1,min_diff,max_diff
      integer                             :: iq, ndim
      ndim = size(Basis%tab_basis)

      IF (debug) THEN
         write (out_unit, *) 'BEGINNING Test'
         flush (out_unit)
      END IF

      allocate (B(Basis%nb*Basis%tab_basis(ndim)%nb))
      allocate (G1(Basis%nq*Basis%tab_basis(ndim)%nb))
      allocate (diff_g(Basis%nq*Basis%tab_basis(ndim)%nb))
      allocate (diff_b(Basis%nb*Basis%tab_basis(ndim)%nb))
      allocate (G2(Basis%nq*Basis%tab_basis(ndim)%nb))
      allocate (B1(Basis%nb*Basis%tab_basis(ndim)%nb))
      allocate (B2(Basis%nb*Basis%tab_basis(ndim)%nb))

      B(:) = CZERO
      B1(:) = CONE
      G1(:) = CONE
      G2(:) = CZERO
      !print *, 'G in', G1
      Call GridTOBasis_nD_cplx(B, G1, Basis)
      Call BasisTOGrid_nD_cplx(G2, B, Basis)
      !print *, 'G out', G2

      print *, '---------------------------------------------------------------------'

      diff_g(:) = ABS(G1(:) - G2(:))

      print *, '---------------------------------------------------------------------'
      Write (out_unit, *) 'maxval(diff_g(:))=', maxval(diff_g(:))
      Write (out_unit, *) 'MINVAL(diff_g(:))=', MINVAL(diff_g(:))
      print *, '---------------------------------------------------------------------'
      G2(:) = CZERO
      B2(:) = CZERO
      B1(:) = CONE
      !print *, 'B in', B1
      Call BasisTOGrid_nD_cplx(G2, B1, Basis)
      Call GridTOBasis_nD_cplx(B2, G2, Basis)
      !print *, 'B out', B2
      diff_b(:) = ABS(B2(:) - B1(:))
      print *, '---------------------------------------------------------------------'
      Write (out_unit, *) 'maxval(diff_b(:))=', maxval(diff_b(:))
      Write (out_unit, *) 'MINVAL(diff_b(:))=', MINVAL(diff_b(:))
      print *, '---------------------------------------------------------------------'

      IF (debug) THEN
         write (out_unit, *) 'END Test'
         flush (out_unit)
      END IF

   END SUBROUTINE test_basitogridgridtobasis

   !> Transform a function from grid representation (GG) to basis representation (BB)
   !! in one dimension, using complex arithmetic.
   !! The transformation is: BB = <basis|grid> * GG
   SUBROUTINE GridTOBasis_1D_cplx(BB, GG, Basis)
      USE QDUtil_m
      IMPLICIT NONE
   
      TYPE(Basis_t), intent(in), target        :: Basis
      COMPLEX(kind=Rkind), intent(inout)       :: BB(:, :, :)  ! Output: Basis representation
      COMPLEX(kind=Rkind), intent(in)          :: GG(:, :, :)  ! Input: Grid representation
      LOGICAL, PARAMETER                       :: debug = .true.
      INTEGER                                  :: i1, i3
   
      !-------------------------------------------------------------------------
      ! Optional debug flushing
      IF (debug) THEN
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
         FLUSH(out_unit)
      END IF
   
   END SUBROUTINE GridTOBasis_1D_cplx

   !> Transform a function from basis representation (BB) to grid representation (GB)
   !! in one dimension, using complex arithmetic.
   !! The transformation is: GB = <x|basis> * BB
   SUBROUTINE BasisTOGrid_1D_cplx(GB, BB, Basis)
      USE QDUtil_m
      IMPLICIT NONE
   
      TYPE(Basis_t), intent(in), target         :: Basis
      COMPLEX(kind=Rkind), intent(inout)        :: GB(:, :, :)  ! Output: Grid representation
      COMPLEX(kind=Rkind), intent(in)           :: BB(:, :, :)  ! Input: Basis representation
      LOGICAL, PARAMETER                        :: debug = .true.
      INTEGER                                   :: i1, i3
   
      !-------------------------------------------------------------------------
      ! Optional debug flushing before processing
      IF (debug) THEN
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
            ! d0gb contains basis functions evaluated on the grid: φ_b(x_q)
            GB(i1, :, i3) = MATMUL(Basis%d0gb(:, :), BB(i1, :, i3))
         END DO
      END DO
   
      !-------------------------------------------------------------------------
      ! Optional debug flushing after processing
      IF (debug) THEN
         FLUSH(out_unit)
      END IF
   
   END SUBROUTINE BasisTOGrid_1D_cplx



   !> Calculates indexing parameters (Ib1, Ib2, Ib3, Iq1, Iq2, Iq3) 
   !! used for reshaping and accessing multidimensional basis/grid data.
   !! This subroutine is typically used in multidimensional simulations (e.g., 3D quantum grids).
   SUBROUTINE Calc_index(Ib1, Ib2, Ib3, Iq1, Iq2, Iq3, Basis)
   TYPE(Basis_t), intent(in), target              :: Basis
   INTEGER, intent(inout), allocatable, optional  :: Iq1(:), Iq2(:), Iq3(:)
   INTEGER, intent(inout), allocatable, optional  :: Ib1(:), Ib2(:), Ib3(:)
   INTEGER                                        :: Ndim      ! Number of spatial dimensions
   INTEGER                                        :: inb       ! Loop index

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

END SUBROUTINE Calc_index

SUBROUTINE BasisTOGrid_nD_cplx(G, B, Basis)
   USE QDUtil_m
   USE NDindex_m

   ! Input/Output arguments
   TYPE(Basis_t), INTENT(IN)                  :: Basis          ! Basis definition and metadata
   COMPLEX(KIND=Rkind), INTENT(IN), TARGET    :: B(:)           ! Input in Basis representation
   COMPLEX(KIND=Rkind), INTENT(INOUT), TARGET :: G(:)           ! Output in Grid representation

 
   ! Optional debug output
   LOGICAL, PARAMETER :: debug = .false.
   ! Pointers for reshaped views of arrays
   COMPLEX(KIND=Rkind), POINTER :: BBG(:, :, :), BBB(:, :, :), GGB(:, :, :)
   COMPLEX(KIND=Rkind), POINTER :: GBB(:, :, :)

   ! Indexing variables (flattening multidimensional data)
   INTEGER :: Ib3, Iq1, Iq2, Ib1, Ib2, Iq3

   ! Temporary working arrays
   COMPLEX(KIND=Rkind), ALLOCATABLE, TARGET :: GBB1(:), GGB2(:)

   ! Integer variables
   INTEGER :: inb, ndim

   IF (debug) THEN
       WRITE(out_unit, *) 'BEGINNING BasisTOGrid_nD_cplx'
       WRITE(out_unit, *) 'intent(in) :: B(:)', B
       CALL Write_Basis(Basis)
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
       WRITE(out_unit, *) 'intent(INOUT) :: G(:)', G
       WRITE(out_unit, *) 'END BasisTOGrid_nD_cplx'
       FLUSH(out_unit)
   END IF

END SUBROUTINE

SUBROUTINE GridTOBasis_nD_cplx(B, G, Basis)
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
   LOGICAL, PARAMETER :: debug = .FALSE.

   ! Optional debug output
   IF (debug) THEN
      write (out_unit, *) 'BEGINNING GridTOBasis_nD_cplx'
      write (out_unit, *) 'intent(in) :: G(:)', G
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
      write (out_unit, *) 'END GridTOBasis_nD_cplx'
      flush (out_unit)
   END IF
END SUBROUTINE

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

   ! Allocate output array Q and index array Tab_iq
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
SUBROUTINE NUFFT(Q, WQ, P, psi_x, psi_p)
   IMPLICIT NONE 
   !< Input parameters >
   
   REAL(KIND=Rkind), INTENT(IN)          :: Q(:,:)      ! Non-uniform position grid (nq_tot x ndim)
   REAL(KIND=Rkind), INTENT(IN)          :: WQ(:)       ! Position weights
   REAL(KIND=Rkind), INTENT(IN)          :: P(:,:)      ! Uniform momentum grid (nq_tot x ndim)
   COMPLEX(KIND=Rkind), INTENT(IN)       :: psi_x(:)    ! Wavefunction in position space

   !< Output >
   COMPLEX(KIND=Rkind), INTENT(INOUT)    :: psi_p(:)    ! Wavefunction in momentum space

   !< Local variables >
   INTEGER                                :: nq_tot, ndim
   INTEGER                                :: iq, j
   REAL(KIND=Rkind)                       :: phase
   COMPLEX(KIND=Rkind)                    :: sum
   REAL(KIND=Rkind), ALLOCATABLE          :: p_vector(:), q_vector(:)
   REAL(KIND=Rkind)                       :: normalization_factor

   ! Get sizes from Q
   nq_tot = SIZE(Q,1)
   ndim   = SIZE(Q,2)
   ALLOCATE(p_vector(ndim))
   ALLOCATE(q_vector(ndim))
   p_vector = ZERO
   q_vector = ZERO

   normalization_factor = ONE / SQRT((TWO * PI) ** ndim)

   ! Loop over all points in momentum space
   DO iq = 1, nq_tot
       sum = CZERO
       p_vector = P(iq, :)

       ! Loop over all points in position space
       DO j = 1, nq_tot
           q_vector = Q(j, :)

           ! Compute multidimensional phase: dot product of momentum and position vectors
           phase = -DOT_PRODUCT(p_vector, q_vector)

           ! Accumulate contribution with position weight and complex exponential
           sum = sum + psi_x(j) * WQ(j) * CMPLX(COS(phase), SIN(phase), KIND=Rkind)
       END DO

       ! Normalization factor 1/(2*pi)^(ndim)
       ! This factor might need adjustment based on the quadrature scheme and convention.
       
       psi_p(iq) = sum*normalization_factor
   END DO
  
   !clean up
   DEALLOCATE(p_vector)
   DEALLOCATE(q_vector)

END SUBROUTINE NUFFT

!>---------------------------------------------------------------------
!> @brief INV_NUFFT: Inverse Non-uniform Fast Fourier Transform
!> @details Transforms a wavefunction from uniform momentum grid P to non-uniform
!>          position grid Q. Uses multidimensional dot product.
!---------------------------------------------------------------------
SUBROUTINE INV_NUFFT(Q, WQ, P, WP, psi_p, psi_x)
   IMPLICIT NONE

   !< Input parameters >
   REAL(KIND=Rkind), INTENT(IN)       :: Q(:,:)       ! Non-uniform position grid (nq_tot x ndim)
   REAL(KIND=Rkind), INTENT(IN)       :: WQ(:)        ! Position weights
   REAL(KIND=Rkind), INTENT(IN)       :: P(:,:)       ! Uniform momentum grid (nq_tot x ndim)
   REAL(KIND=Rkind), INTENT(IN)       :: WP(:)        ! Momentum weights
   COMPLEX(KIND=Rkind), INTENT(IN)    :: psi_p(:)     ! Wavefunction in momentum space
   
   !< Output >
   COMPLEX(KIND=Rkind), INTENT(OUT)   :: psi_x(:)     ! Wavefunction in position space

   !< Local variables >
   INTEGER                            :: nq_tot, ndim
   INTEGER                            :: iq, j
   REAL(KIND=Rkind)                   :: phase
   COMPLEX(KIND=Rkind)                :: sum
   REAL(KIND=Rkind), ALLOCATABLE      :: p_vector(:), q_vector(:)
   REAL(KIND=Rkind)                   :: normalization_factor

   ! Get sizes from Q
   nq_tot = SIZE(Q,1)
   ndim   = SIZE(Q,2)
   ALLOCATE(p_vector(ndim))
   ALLOCATE(q_vector(ndim))
   p_vector = ZERO
   q_vector = ZERO
   normalization_factor = ONE / SQRT((TWO * PI) ** ndim)

   ! Loop over all points in position space
   DO iq = 1, nq_tot
       sum = CZERO
       q_vector = Q(iq, :)

       ! Loop over all points in momentum space
       DO j = 1, nq_tot
           p_vector = P(j, :)

           phase = DOT_PRODUCT(p_vector, q_vector)

           sum = sum + psi_p(j) * WP(j) * CMPLX(COS(phase), SIN(phase), KIND=Rkind)
       END DO

       ! Normalization factor 1/(2*pi)^(ndim)
       psi_x(iq) = sum * normalization_factor
   END DO

   !clean up
   DEALLOCATE(p_vector)
   DEALLOCATE(q_vector)
END SUBROUTINE INV_NUFFT



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



SUBROUTINE ApplyExpVdtAndTransformToMomentum(psi_x, psi_p, V, dt, Q, WQ, P, WP,transform)
   ! ------------------------------------------------------------
   ! Apply exp(-i V(q) dt/2) on psi_x, then Fourier transform
   ! each surface to momentum space (psi_p).
   !
   ! Input:
   !   Q(nq,ndim) : position quadrature
   !   WQ(nq)     : quadrature weights in position
   !   P(np,ndim) : momentum quadrature
   !   WP(np)     : quadrature weights in momentum
   !   V(nq,nsurf,nsurf) : potential (diagonal in surface index)
   !   dt         : time step
   !
   ! Input/Output:
   !   psi_x(:)   : wavefunction in position space (vectorized: nq*nsurf)
   !   psi_p(:)   : wavefunction in momentum space (vectorized: np*nsurf)
   ! ------------------------------------------------------------
 
   IMPLICIT NONE
   REAL(KIND=Rkind), INTENT(IN) :: Q(:,:), P(:,:), WQ(:), WP(:)
   REAL(KIND=Rkind), INTENT(IN) :: V(:,:,:)
   REAL(KIND=Rkind), INTENT(IN) :: dt
   COMPLEX(KIND=Rkind), INTENT(IN), TARGET :: psi_x(:)
   COMPLEX(KIND=Rkind), INTENT(INOUT),   TARGET :: psi_p(:)
   LOGICAL, INTENT(IN),OPTIONAL :: transform
 
   ! Work arrays and pointers
   COMPLEX(KIND=Rkind), ALLOCATABLE, TARGET :: tmp(:)
   COMPLEX(KIND=Rkind), POINTER :: Psi_gb(:, :), Psi_pb(:, :), Psi_gbtemp(:, :)
 
   INTEGER :: nq, nsurf, np,ndim
   INTEGER :: iq, isurf,jsurf
   COMPLEX(KIND=Rkind) :: phase
 
   nq    = SIZE(Q,1)
   np    = SIZE(P,1)
   nsurf = SIZE(V,2)
   ndim = SIZE(Q,2)
 
   ! --- Map input vectors to 2D views ---
   Psi_gb(1:nq,1:nsurf) => psi_x
   Psi_pb(1:np,1:nsurf) => psi_p

   ! --- Temporary buffer for modified wavefunction ---
   ALLOCATE(tmp(nq*nsurf))
   tmp = CZERO
   Psi_gbtemp(1:nq,1:nsurf) => tmp

 
   ! --- Apply exp(-i V dt/2) (diagonal only) ---
   DO iq = 1, nq
      DO isurf = 1, nsurf
         phase = CZERO
         DO jsurf = 1, nsurf
           phase =  phase  -EYE*HALF * dt * V(iq,isurf,jsurf)
           Psi_gbtemp(iq,isurf) = Psi_gb(iq,isurf) *EXP(phase)
         END DO
      END DO
   END DO
 
   ! --- Transform each surface separately to momentum space ---
   IF(PRESENT(transform)) THEN
      IF(transform) THEN
         DO isurf = 1, nsurf
            CALL NUFFT(Q, WQ, P, Psi_gbtemp(:,isurf), Psi_pb(:,isurf))
         END DO

      ELSE 
         WRITE(out_unit,*) 'ApplyExpVdtAndTransformToMomentum: transform = .FALSE.'
         psi_p(:) = tmp(:)
      END IF
   ELSE
      ! Default to transform
      DO isurf = 1, nsurf
         CALL NUFFT(Q, WQ, P, Psi_gbtemp(:,isurf), Psi_pb(:,isurf))
      END DO
   END IF
 
   ! --- Cleanup ---
   DEALLOCATE(tmp)
   NULLIFY(Psi_gb, Psi_pb, Psi_gbtemp)
 
 END SUBROUTINE ApplyExpVdtAndTransformToMomentum

 SUBROUTINE ApplyExpTdtAndTransformToPosition(psi_p, psi_x, Q, WQ, P, WP, dt)
   ! ------------------------------------------------------------
   ! Applies exp(-i T dt) in momentum space
   ! then transforms the wavefunction back to position space
   ! using an inverse NUFFT.
   !
   ! Steps:
   !   1. Applies the kinetic energy propagator: psi_p = exp(-i*T*dt) * psi_p
   !   2. Performs an inverse NUFFT: psi_x = F^{-1}(psi_p)
   !
   ! Inputs:
   !   P   (np, ndim) : Momentum grid points
   !   WP  (np)       : Quadrature weights for momentum space
   !   Q   (nq, ndim) : Position grid points
   !   WQ  (nq)       : Quadrature weights for position space
   !   dt             : Time step
   !
   ! Input/Output:
   !   psi_p(np*nsurf) : Wavefunction in momentum space (vector input). Overwritten with intermediate result.
   !   psi_x(nq*nsurf) : Wavefunction in position space (output).
   ! ------------------------------------------------------------
   IMPLICIT NONE
 
   ! Arguments
   REAL(KIND=Rkind),    INTENT(IN)    :: P(:,:), WP(:), Q(:,:), WQ(:), dt
   COMPLEX(KIND=Rkind), INTENT(INOUT),TARGET :: psi_p(:) 
   COMPLEX(KIND=Rkind), INTENT(OUT)  ,TARGET :: psi_x(:) 
 
   ! Local variables
   INTEGER :: np, nq, ndim, nsurf, isurf, ip, i
   REAL(KIND=Rkind) :: Tval
   REAL(KIND=Rkind), ALLOCATABLE :: GGdef(:,:)
   COMPLEX(KIND=Rkind), ALLOCATABLE ,TARGET :: work_array(:) ! Work array for the operation
 
   ! Pointers for easier array manipulation
   COMPLEX(KIND=Rkind), POINTER :: Psi_gb(:,:), Psi_pb(:,:), Psi_pb_work(:,:)
 
   ! Get dimensions from input arrays
   np    = SIZE(P, 1)
   nq    = SIZE(Q, 1)  ! Corrected: Get nq from input argument Q
   ndim  = SIZE(P, 2)
   nsurf = SIZE(psi_p) / np ! Calculate number of surfaces/states
 
   ! Check output array size
   IF (SIZE(psi_x) /= nq * nsurf) THEN
     WRITE(*,*) 'ERROR in ApplyExpTdtAndTransformToPosition: size mismatch for psi_x'
     WRITE(*,*) 'Expected size: ', nq * nsurf, ' Actual size: ', SIZE(psi_x)
     STOP
   END IF
 
   ! Allocate and get the inverse mass metric 
   ALLOCATE(GGdef(ndim, ndim))
   CALL get_Qmodel_GGdef(GGdef) 
 
   ! Map 1D vectors to 2D arrays using pointers for clarity
   Psi_gb(1:nq, 1:nsurf)     => psi_x
   Psi_pb(1:np, 1:nsurf)     => psi_p
 
   ! Allocate a work array and copy input psi_p into it
   ALLOCATE(work_array(SIZE(psi_p)))
   work_array(:) = psi_p(:) ! Copy the input wavefunction
   Psi_pb_work(1:np, 1:nsurf) => work_array
 
   ! Apply exp(-i * T * dt) in momentum space
   ! T = 0.5 * sum_{i} GGdef(i,i) * P_i^2
   DO isurf = 1, nsurf
     DO ip = 1, np
       Tval = ZERO
       ! Calculate kinetic energy T for this momentum point P(ip, :)
       DO i = 1, ndim
         Tval = Tval + HALF * GGdef(i, i) * P(ip, i)**2
       END DO
       ! Apply the propagator: Multiply the wavefunction by exp(-i*T*dt)
       Psi_pb_work(ip, isurf) = Psi_pb_work(ip, isurf) * EXP(-EYE * Tval * dt)
     END DO
   END DO
 
   ! Inverse NUFFT: Transform from momentum (p) to position (x) space
   DO isurf = 1, nsurf
     CALL INV_NUFFT(Q, WQ, P, WP, Psi_pb_work(:, isurf), Psi_gb(:, isurf))
   END DO
 
   ! Clean up
   DEALLOCATE(work_array, GGdef)
   ! Nullify pointers to avoid potential dangling references
   NULLIFY(Psi_gb, Psi_pb, Psi_pb_work)
 
 END SUBROUTINE ApplyExpTdtAndTransformToPosition




   
    SUBROUTINE d2psi1D_cplx(d2GGB,GGB,Basis)
    USE QDUtil_m
    complex(kind=Rkind), intent(in)                            :: GGB(:,:,:)
    complex(kind=Rkind), intent(inout)                         :: d2GGB(:,:,:)
    TYPE(Basis_t),intent(in)                                   :: Basis
    
    !locals variables---------------------------------------------------------
    
    
    logical, parameter                                         :: debug = .true.
  
    integer                                                    :: i1,i3
    
    
    !debuging----------------------------------------------------------------
    
    
    IF (debug) THEN
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

 SUBROUTINE pdv2psi_nD(d2psi_g, psi_g, Basis,ib)

    USE  QDUtil_m
    TYPE(Basis_t), intent(in), target               :: Basis
    complex(kind=Rkind), intent(in), target         :: psi_g(:)
    complex(kind=Rkind), intent(inout), target      :: d2psi_g(:)
    integer ,intent(in)                             :: ib
    complex(kind=Rkind), pointer                    :: psi_ggb(:, :, :)
    complex(kind=Rkind), pointer                    :: d2gg(:, :)
    complex(kind=Rkind), pointer                    :: d2psi_ggb(:, :, :)
    logical, parameter                              :: debug = .true.
    integer                                         :: i1, i3, Ndim
    integer, allocatable                            :: Iq1(:), Iq2(:), Iq3(:), Ib1(:), Ib2(:), Ib3(:)

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

 SUBROUTINE psi_per_surf(Gsurf,G,Basis,ibe)

     IMPLICIT NONE
     TYPE(Basis_t), intent(in)                              :: Basis
     complex(kind=Rkind), intent(in) ,target                :: G(:)
     complex(kind=Rkind), intent(inout) ,allocatable        :: Gsurf(:)
     integer  , intent(in)                                  :: ibe
     
     complex(kind=Rkind), pointer                           :: psi_gb(:, :)
     integer                                                :: nq,nsurf

     IF(Allocated(Gsurf)) deallocate(Gsurf)

     nq = Basis%nq
     nsurf = Basis%tab_basis(size(Basis%tab_basis))%nb
     allocate(Gsurf(nq))
     Gsurf = CZERO

      psi_gb(1:nq, 1:nsurf) => G
    
      Gsurf(1:nq) =  psi_gb(1:nq, ibe) 


  END SUBROUTINE


  SUBROUTINE Calc_reduced_Density_surf(Rdensity,G,Basis)

     USE QDUtil_m

    TYPE(Basis_t), intent(in)                                    :: Basis
    complex(Kind=Rkind),intent(in)                               :: G(:)                                         
   TYPE(REDUCED_DENSIRY_t),intent(inout)                         :: Rdensity           

    !Locals variables ----------------------------------------------------------------------------

    integer, allocatable                                         :: Tab_iq(:)
    integer                                                      :: ndim, Iq,inb,ib
    real(Kind=Rkind)                                             :: W
    logical                                                      :: Endloop
      ndim  = SIZE(Basis%tab_basis) - 1
      allocate (Tab_iq(ndim))
      call Init_tab_ind(Tab_iq, Basis%NDindexq) 
         Iq = 0
         DO
            Iq = Iq + 1

            CALL increase_NDindex(Tab_iq, Basis%NDindexq, Endloop)
            IF (Endloop) exit
             W = ONE
            DO inb = 1, ndim  
               W = W*Basis%tab_basis(inb)%w(Tab_iq(inb))
            END DO
            If(ndim ==1) then
               Rdensity%prob(Iq) = Rdensity%prob(Iq) + conjg(G(Iq))*G(Iq)*W
            Else
               DO ib = 1, ndim
                   Rdensity%tab_prob(ib)%prob(Tab_iq(ib)) = Rdensity%tab_prob(ib)%prob(Tab_iq(ib)) + conjg(G(Iq))*G(Iq)*W
               END DO
            End If
         END DO            
      Deallocate(Tab_iq)
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

 SUBROUTINE Calc_Rdensity_Tot(Rdensity,Basis)
    
  USE  QDUtil_m  
  TYPE(REDUCED_DENSIRY_t),intent(inout)         :: Rdensity
  TYPE(Basis_t), intent(in), target             :: Basis

  integer                                       :: ndim,Ib

  ndim  = SIZE(Basis%tab_basis) - 1

  allocate(Rdensity%Norm(ndim))
      

       Rdensity%Norm(:) = ZERO
  
     DO Ib=1,ndim
       
           Rdensity%Norm(Ib) = sum(Rdensity%tab_prob(Ib)%prob(:))
  
     END DO

     Write(*,*) 'RedTot',Rdensity%Norm(:)


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

 SUBROUTINE Calc_reduced_density(Rdensity,B,Basis)
      USE  QDUtil_m
      TYPE(Basis_t), intent(in), target               :: Basis
      complex(kind=Rkind), intent(in), target         :: B(:)
      TYPE(REDUCED_DENSIRY_t),intent(inout)           :: Rdensity

        TYPE(REDUCED_DENSIRY_t)                       :: Rdensitytemp
        complex(kind=Rkind),allocatable               :: G(:),Gsurf(:)
        integer                                       :: nsurf,nq,ibe,ndim,ib

         nsurf = Basis%tab_basis(size(Basis%tab_basis))%nb
         nq = nsurf*Basis%nq
         ndim = size(Basis%tab_basis)-1

         call Rdensity_Dealloc(Rdensity,Basis)
         call Rdensity_alloc(Rdensity,Basis)          
         allocate(G(nq))
         CALL BasisTOGrid_nD_cplx(G,B,Basis)
        
       DO ibe = 1,nsurf

          CALL psi_per_surf(Gsurf,G,Basis,ibe)
          call  Rdensity_alloc(Rdensitytemp,Basis)
          CALL Calc_reduced_Density_surf(Rdensitytemp,Gsurf,Basis)
           If(ndim==1) then
             Rdensity%prob(:)  = Rdensity%prob(:)+Rdensitytemp%prob(:)
          Else
           Do ib =1, ndim
              Rdensity%tab_prob(ib)%prob(:)  = Rdensity%tab_prob(ib)%prob(:)+Rdensitytemp%tab_prob(ib)%prob(:)
           END DO
          End If
          call   Rdensity_Dealloc(Rdensitytemp,Basis)
          

       END DO   
      deallocate(G)

 END SUBROUTINE


 SUBROUTINE Calc_sum_psi(Int,G1,G2,Basis)
   USE QDUtil_m
   logical, parameter                            :: debug = .false.
   TYPE(Basis_t),intent(in)                      :: Basis
   complex(kind=Rkind), intent(in)   ,target     :: G1(:),G2(:)
   complex(kind=Rkind) ,intent(inout)            :: Int
   !Locals variabls ----------------------------------------------------------
   complex(kind=Rkind), pointer                  :: psi_gb1(:, :),psi_gb2(:, :)
   logical                                       :: Endloop_q
   complex(kind=Rkind), allocatable              :: Intel(:)
   real(kind=Rkind), allocatable                 :: Nel(:)
   real(kind=Rkind)                              :: W
   integer, allocatable                          :: Tab_iq(:)
   integer                                       :: iq, inbe,inb,nq,nsurf,ndim
   IF (debug) THEN
      write (out_unit, *) 'Beging Evaluating integral'
      flush (out_unit)
   END IF
   nq = Basis%nq
   nsurf = Basis%tab_basis(size(Basis%tab_basis))%nb
   ndim = size(Basis%tab_basis) - 1
   allocate (Nel(nsurf))
   allocate (Intel(nsurf)) 
   allocate (Tab_iq(ndim))
   psi_gb1(1:nq,1:nsurf) => G1
   psi_gb2(1:nq,1:nsurf) => G2
 Nel(:) = ZERO
 Intel(:)=ZERO
 Int = ZERO
 DO inbe = 1, Basis%tab_basis(size(Basis%tab_basis))%nb !electronic state
    Intel(inbe) = CZERO
    Call Init_tab_ind(Tab_iq, Basis%NDindexq)
    iq = 0
    DO
       iq = iq + 1
       CALL increase_NDindex(Tab_iq, Basis%NDindexq, Endloop_q)
       IF (Endloop_q) exit
      
        W = ONE
        DO inb = 1, ndim
           W = W*Basis%tab_basis(inb)%w(tab_iq(inb))
        END DO
       Nel(inbe) = Nel(inbe) + REAL(CONJG(psi_gb1(iq, inbe))*psi_gb1(iq, inbe)*W,kind=Rkind)
       Intel(inbe) = Intel(inbe) + REAL(CONJG(psi_gb1(iq, inbe))*psi_gb2(iq, inbe)*W,kind=Rkind)
    END DO
 END DO
 Int = sum(Intel)/(Sum(Nel)**2)
 Deallocate (Tab_iq)
 Deallocate(Intel,Nel)
 IF (debug) THEN
    write (out_unit, *) 'END Evaluating integral'
    flush (out_unit)
 END IF
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
