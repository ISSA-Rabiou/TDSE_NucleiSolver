module psi_m

   USE param_WP0_m
   USE Basis_m
   USE NDindex_m
   implicit none

   TYPE :: psi_t
   
      type(Basis_t), POINTER           :: Basis
      REAL(kind= Rkind), ALLOCATABLE    :: RVec(:)
      COMPLEX(kind= Rkind), ALLOCATABLE :: CVec(:)
      LOGICAL                            :: Grid = .TRUE.

   CONTAINS

      PRIVATE
      PROCEDURE, PASS         :: Copy_psi    ! Copy content from other psi1 to psi2 instance,
      GENERIC, PUBLIC         :: ASSIGNMENT(=) => Copy_psi

   END TYPE psi_t

   public :: psi_t, write_psi, init_psi, dealloc_psi, write_psi_grid
   public :: write_psi_basis, Calc_Norm_OF_PsiBasis, Calc_Norm_OF_PsiGrid, Calc_Norm_OF_Psi
   public:: Projection,Projection1,psi_init_GWP,Ecrire_psi
contains


   SUBROUTINE copy_psi(psi_out, psi_in)

      CLASS(psi_t), intent(in)     :: psi_in
      CLASS(psi_t), intent(inout)  :: psi_out

      IF (allocated(psi_in%RVec)) THEN
         !write(out_unit,*) 'Coping psi_in in psi_out (real):'
         psi_out%RVec(:) = psi_in%RVec(:)
         !CALL init_Basis1_TO_Basis2 (psi_in%Basis,psi_out%Basis)
         !CALL  construct_primitive_basis(psi_out%Basis)
         psi_out%Basis => psi_in%Basis
         !write(out_unit,*) 'END Coping psi_in in psi_out'
      END IF
      IF (allocated(psi_in%CVec)) THEN
         !write(out_unit,*) 'Coping psi_in in psi_out (complex):'
         psi_out%CVec(:) = psi_in%CVec(:)
         !CALL init_Basis1_TO_Basis2 (psi_in%Basis,psi_out%Basis)
         !CALL  construct_primitive_basis(psi_out%Basis)
         psi_out%Basis => psi_in%Basis
         !write(out_unit,*) 'END Coping psi_in in psi_out'
      END IF

   END SUBROUTINE copy_psi

   SUBROUTINE init_psi(psi, Basis, cplx, grid)
      USE Basis_m
   
      !======================================================================
      ! PURPOSE:
      !   Initialize a wavefunction (psi_t) structure, allocating
      !   memory for either grid or basis representation.
      !
      ! ARGUMENTS:
      !   psi   : (INOUT) Wavefunction to initialize
      !   Basis : (IN)   Basis set (must be allocated and valid)
      !   cplx  : (IN)   .TRUE. => complex representation, .FALSE. => real
      !   grid  : (IN)   .TRUE. => allocate in grid representation
      !                     .FALSE. => allocate in basis representation
      !
      ! SAFETY:
      !   - Avoids memory leaks by always deallocating previous data first.
      !   - Uses STAT= in ALLOCATE to catch memory allocation errors.
      !   - Checks validity of the Basis before allocation.
      !======================================================================
   
      TYPE(psi_t),     INTENT(INOUT) :: psi
      TYPE(Basis_t),   INTENT(IN), TARGET :: Basis
      LOGICAL,         INTENT(IN)    :: cplx, grid
   
      INTEGER :: alloc_size   ! size of vector to allocate
      INTEGER :: istat        ! allocation status
      LOGICAL , PARAMETER :: debug = .FALSE.
      INTEGER :: i
   
      !-----------------------------------------------------------
      ! Step 1: Deallocate any existing psi data
      !-----------------------------------------------------------
      CALL dealloc_psi(psi)
   
      !-----------------------------------------------------------
      ! Step 2: Validate basis
      !-----------------------------------------------------------
      IF (.NOT. Basis_IS_allocated(Basis)) THEN
         WRITE(out_unit,*) 'ERROR in init_psi: the Basis is not initialized'
         RETURN
      ELSE
         IF (debug) THEN
            WRITE(out_unit,*) 'Basis is allocated in Init_psi'
         END IF
      END IF
   
      IF (Basis%nb < 1) THEN
         WRITE(out_unit,*) 'ERROR in init_psi: Basis%nb < 1!'
         RETURN
      ELSE
         IF (debug) THEN
            WRITE(out_unit,*) 'Basis%nb is ', Basis%nb
         END IF
      END IF
   
      !-----------------------------------------------------------
      ! Step 3: Link psi to basis and set representation flag
      !-----------------------------------------------------------
      psi%Basis => Basis
      psi%Grid  = grid
   
      !-----------------------------------------------------------
      ! Step 4: Determine allocation size
      !-----------------------------------------------------------
      IF (grid) THEN
         ! GRID representation
         IF (ALLOCATED(Basis%tab_basis)) THEN
            IF (debug) THEN
               WRITE(out_unit,*) 'Basis%tab_basis is allocated in Init_psi and size is ', SIZE(Basis%tab_basis)
            END IF
            alloc_size = Basis%nq * Basis%tab_basis(SIZE(Basis%tab_basis))%nb
         ELSE
            IF (debug) THEN
               WRITE(out_unit,*) 'Basis%tab_basis is not allocated in Init_psi'
            END IF
            alloc_size = Basis%nq
         END IF
      ELSE
         ! BASIS representation
         IF (ALLOCATED(Basis%tab_basis)) THEN
            IF (debug) THEN
               WRITE(out_unit,*) 'Basis%tab_basis is allocated in Init_psi and size is ', SIZE(Basis%tab_basis)
            END IF
            alloc_size = Basis%nb * Basis%tab_basis(SIZE(Basis%tab_basis))%nb
         ELSE
            IF (debug) THEN
               WRITE(out_unit,*) 'Basis%tab_basis is not allocated in Init_psi'
            END IF
            alloc_size = Basis%nb
         END IF
      END IF
   
      !-----------------------------------------------------------
      ! Step 5: Allocate the correct coefficient vector and initialize to zero
      ! then write the initial state if debug is true
      !-----------------------------------------------------------
      istat = 0
      IF (cplx) THEN
         ALLOCATE(psi%CVec(alloc_size), STAT=istat)
         psi%CVec = CZERO
         IF (debug) THEN
            CALL Write_VecMat(psi%CVec, out_unit, 5, info='psi%CVec')
         END IF
      ELSE
         ALLOCATE(psi%RVec(alloc_size), STAT=istat)
         psi%RVec = ZERO
         IF (debug) THEN
            CALL Write_VecMat(psi%CVec, out_unit, 5, info='psi%CVec')
         END IF
      END IF
   
      !-----------------------------------------------------------
      ! Step 6: Handle allocation errors
      !-----------------------------------------------------------
      IF (istat /= 0) THEN
         PRINT *, 'ERROR in init_psi: Allocation failed, STAT=', istat
         CALL dealloc_psi(psi)
         RETURN
      END IF
   
   END SUBROUTINE init_psi


   SUBROUTINE dealloc_psi(psi)
      !======================================================================
      ! PURPOSE:
      !   Deallocate the memory associated with a wavefunction (psi_t type)
      !   and nullify its basis pointer.
      !
      ! ARGUMENTS:
      !   psi : (INOUT) The wavefunction to be deallocated.
      !
      ! NOTES:
      !   - Any allocated real or complex coefficient vectors are released.
      !   - The basis pointer is nullified to avoid dangling references.
      !======================================================================
   
      TYPE(psi_t), INTENT(INOUT) :: psi
   
      !-----------------------------------------------------------
      ! Nullify the basis pointer (safe even if already nullified)
      !-----------------------------------------------------------
      NULLIFY(psi%Basis)
   
      !-----------------------------------------------------------
      ! Deallocate real coefficient vector if allocated
      !-----------------------------------------------------------
      IF (ALLOCATED(psi%RVec)) THEN
         DEALLOCATE(psi%RVec)
      END IF
   
      !-----------------------------------------------------------
      ! Deallocate complex coefficient vector if allocated
      !-----------------------------------------------------------
      IF (ALLOCATED(psi%CVec)) THEN
         DEALLOCATE(psi%CVec)
      END IF
   
   END SUBROUTINE dealloc_psi


   SUBROUTINE write_psi(psi, psi_cplx, print_psi_grid, print_basis, t, int_print, real_part)
      !======================================================================
      ! PURPOSE:
      !   Write the wavefunction psi either in grid representation
      !   or basis representation, depending on user request.
      !
      ! ARGUMENTS:
      !   psi            : quantum state (type psi_t)
      !   psi_cplx       : .true. => print complex values, .false. => real only
      !   print_psi_grid : .true. => print in grid representation
      !   print_basis    : .true. => also print the basis
      !   t              : (optional) time
      !   int_print      : (optional) output unit number
      !   real_part      : .true. => print only the real part
      !
      ! NOTES:
      !   - If the representation of psi does not match the desired
      !     output (grid/basis), a temporary conversion is performed.
      !   - Temporary psi objects are always deallocated before return.
      !======================================================================
   
      TYPE(psi_t),        INTENT(IN)           :: psi
      LOGICAL,            INTENT(IN)           :: print_psi_grid, print_basis, psi_cplx, real_part
      INTEGER, OPTIONAL,  INTENT(IN)           :: int_print
      REAL(KIND=Rkind), OPTIONAL, INTENT(IN)   :: t
      LOGICAL, PARAMETER                      :: debug = .true.
   
      ! Local temporary variables for conversion
      TYPE(psi_t) :: psi_g   ! Temporary wavefunction in grid representation
      TYPE(psi_t) :: psi_b   ! Temporary wavefunction in basis representation
      LOGICAL     :: alloc_g, alloc_b
   
      alloc_g = .FALSE.
      alloc_b = .FALSE.
   
      !===============================
      ! Main printing logic
      !===============================
      IF (print_psi_grid) THEN
         !---------------------------------------------------------
         ! CASE 1: User requested grid representation
         !---------------------------------------------------------
         IF (psi%Grid) THEN
            ! Already in grid representation → print directly
            CALL call_write_psi_grid(psi, psi_cplx, t, int_print, real_part)
         ELSE
            ! Convert basis → grid
            CALL init_psi(psi_g, psi%Basis, cplx=.TRUE., grid=.TRUE.)
            alloc_g = .TRUE.
            CALL BasisTOGrid_nD_cplx(psi_g%CVec, psi%CVec, psi%Basis)
   
            CALL call_write_psi_grid(psi_g, psi_cplx, t, int_print, real_part)
         END IF
   
      ELSE
         !---------------------------------------------------------
         ! CASE 2: User requested basis representation
         !---------------------------------------------------------
         IF (psi%Grid) THEN
            ! Convert grid → basis
            CALL init_psi(psi_b, psi%Basis, cplx=.TRUE., grid=.FALSE.)
            alloc_b = .TRUE.
            CALL GridTOBasis_nD_cplx(psi_b%CVec, psi%CVec, psi%Basis)
   
            CALL call_write_psi_basis(psi_b, psi_cplx, t, int_print, real_part)
         ELSE
            ! Already in basis representation → print directly
            CALL call_write_psi_basis(psi, psi_cplx, t, int_print, real_part)
         END IF
      END IF
   
      !===============================
      ! Optional: print the basis
      !===============================
      IF (print_basis) THEN
         CALL Write_Basis(psi%Basis)
      END IF
   
      !===============================
      ! Debug: flush output
      !===============================
      IF (debug) THEN
         FLUSH(out_unit)
      END IF
   
      !===============================
      ! Memory cleanup
      !===============================
      IF (alloc_g) CALL dealloc_psi(psi_g)
      IF (alloc_b) CALL dealloc_psi(psi_b)
   
   CONTAINS
   
      !-------------------------------------------------------------------
      ! Helper routine to write psi in grid representation
      ! This avoids code repetition when checking PRESENT arguments.
      !-------------------------------------------------------------------
      SUBROUTINE call_write_psi_grid(psi_data, cplx_flag, t_opt, nio_opt, real_flag)
         TYPE(psi_t),        INTENT(IN) :: psi_data
         LOGICAL,            INTENT(IN) :: cplx_flag, real_flag
         REAL(KIND=Rkind), OPTIONAL, INTENT(IN) :: t_opt
         INTEGER, OPTIONAL,  INTENT(IN) :: nio_opt
   
         IF (PRESENT(t_opt) .AND. PRESENT(nio_opt)) THEN
            CALL write_psi_grid(psi_data, print_cplx=cplx_flag, t=t_opt, nio=nio_opt, real_part=real_flag)
         ELSEIF (PRESENT(t_opt)) THEN
            CALL write_psi_grid(psi_data, print_cplx=cplx_flag, t=t_opt, real_part=real_flag)
         ELSEIF (PRESENT(nio_opt)) THEN
            CALL write_psi_grid(psi_data, print_cplx=cplx_flag, nio=nio_opt, real_part=real_flag)
         ELSE
            CALL write_psi_grid(psi_data, print_cplx=cplx_flag, real_part=real_flag)
         END IF
      END SUBROUTINE call_write_psi_grid
   
      !-------------------------------------------------------------------
      ! Helper routine to write psi in basis representation
      !-------------------------------------------------------------------
      SUBROUTINE call_write_psi_basis(psi_data, cplx_flag, t_opt, nio_opt, real_flag)
         TYPE(psi_t),        INTENT(IN) :: psi_data
         LOGICAL,            INTENT(IN) :: cplx_flag, real_flag
         REAL(KIND=Rkind), OPTIONAL, INTENT(IN) :: t_opt
         INTEGER, OPTIONAL,  INTENT(IN) :: nio_opt
   
         IF (PRESENT(t_opt) .AND. PRESENT(nio_opt)) THEN
            CALL write_psi_basis(psi_data, print_cplx=cplx_flag, t=t_opt, nio=nio_opt, real_part=real_flag)
         ELSEIF (PRESENT(t_opt)) THEN
            CALL write_psi_basis(psi_data, print_cplx=cplx_flag, t=t_opt, real_part=real_flag)
         ELSEIF (PRESENT(nio_opt)) THEN
            CALL write_psi_basis(psi_data, print_cplx=cplx_flag, nio=nio_opt, real_part=real_flag)
         ELSE
            CALL write_psi_basis(psi_data, print_cplx=cplx_flag, real_part=real_flag)
         END IF
      END SUBROUTINE call_write_psi_basis
   
   END SUBROUTINE write_psi


   !=====================================================================
   !> @brief
   !>   Hagedorn projection of psi_dt_1 onto psi_dt_2.
   !>   Supports N-dimensional basis sets with sequential 1D transformations.
   !>
   !> @details
   !>   - Uses temporary pointers and allocatables for intermediate transformations.
   !>   - Computes initial and final norms for verification.
   !>   - Debug output is local and only printed if the local flag is .true.
   !=====================================================================

   SUBROUTINE Projection(psi_dt_2, psi_dt_1)
      TYPE(psi_t), intent(in), target                  :: psi_dt_1
      TYPE(psi_t), intent(inout), target               :: psi_dt_2
      complex(kind=Rkind), pointer                     :: BBB1(:, :, :), BBB2(:, :, :)
      complex(kind=Rkind), allocatable, target         :: B1(:), B2(:)
      real(kind=Rkind)                                 :: Norm0
      logical, parameter                               :: debug = .true.
      integer                                          :: inb, i1, i3, Ndim
      Integer, allocatable                             :: Ib1(:), Ib2(:), Ib3(:)
   
      ! Compute index arrays for multidimensional basis representation
      call Calc_index(Ib1=Ib1, Ib2=Ib2, Ib3=Ib3, Basis=psi_dt_1%Basis)
   
      ! Determine the number of dimensions (basis functions)
      Ndim = size(psi_dt_1%Basis%tab_basis) - 1
   
      ! Compute initial norm of psi_dt_1 for verification
      call Calc_Norm_OF_psi(psi_dt_1, Norm0)
      write (out_unit, *) 'Begin Hagedorn projection <psi|psi>=', Norm0
   
      ! Perform the Hagedorn projection
      if (Ndim == 1) then
         ! Link the CVec arrays to pointer views for easier access
         BBB1(1:Ib1(1), 1:Ib2(1), 1:Ib3(1)) => psi_dt_1%CVec
         BBB2(1:Ib1(1), 1:Ib2(1), 1:Ib3(1)) => psi_dt_2%CVec
   
         ! Initialize psi_dt_2%CVec to zero
         psi_dt_2%CVec(:) = CZERO
   
         ! Apply the basis transformation along dimension 1
         DO i3 = 1, ubound(BBB1, dim=3)
            DO i1 = 1, ubound(BBB1, dim=1)
               BBB2(i1, :, i3) = matmul(psi_dt_1%Basis%tab_basis(1)%S, BBB1(i1, :, i3))
            END DO
         END DO
      else
         ! Allocate temporary vector B1 to hold intermediate data
         allocate (B1(Ib1(1)*Ib2(1)*Ib3(1)))
   
         ! Link pointers for manipulation
         BBB1(1:Ib1(1), 1:Ib2(1), 1:Ib3(1)) => psi_dt_1%CVec
         BBB2(1:Ib1(1), 1:Ib2(1), 1:Ib3(1)) => B1
   
         ! First transformation along dimension 1
         DO i3 = 1, ubound(BBB1, dim=3)
            DO i1 = 1, ubound(BBB1, dim=1)
               BBB2(i1, :, i3) = matmul(psi_dt_1%Basis%tab_basis(1)%S, BBB1(i1, :, i3))
            END DO
         END DO
   
         ! Loop through remaining dimensions to apply projection
         DO inb = 2, Ndim
            ! Allocate temporary buffer B2
            allocate (B2(Ib1(inb)*Ib2(inb)*Ib3(inb)))
   
            ! Link BBB1 to B1 (previous result), BBB2 to B2 (target of transformation)
            BBB1(1:Ib1(inb), 1:Ib2(inb), 1:Ib3(inb)) => B1
            BBB2(1:Ib1(inb), 1:Ib2(inb), 1:Ib3(inb)) => B2
   
            ! Apply transformation along current dimension
            DO i3 = 1, ubound(BBB1, dim=3)
               DO i1 = 1, ubound(BBB1, dim=1)
                  BBB2(i1, :, i3) = matmul(psi_dt_1%Basis%tab_basis(inb)%S, BBB1(i1, :, i3))
               END DO
            END DO
   
            ! Prepare for next iteration
            B1 = B2
            B2 = CZERO
            deallocate (B2)
   
            ! If last iteration, copy final result to psi_dt_2%CVec
            if (inb == Ndim) psi_dt_2%CVec = B1
         END DO
      END IF
   
      ! Compute final norm for verification
      call Calc_Norm_OF_psi(psi_dt_2, Norm0)
   
      ! Clean up dynamically allocated memory
      deallocate(Ib1, Ib2, Ib3)
      if (allocated(B1)) deallocate(B1)
      if (allocated(B2)) deallocate(B2)
   
      write (out_unit, *) 'END Hagedorn projection <psi|psi>=', Norm0
   END SUBROUTINE Projection

   SUBROUTINE Projection1(psi_dt_2, psi_dt_1)
   TYPE(psi_t), intent(in), target                  :: psi_dt_1
   TYPE(psi_t), intent(inout), target               :: psi_dt_2
   complex(kind= Rkind), pointer                    :: BBB1(:, :, :), BBB2(:, :, :)
   complex(kind= Rkind), allocatable, target        :: B1(:), B2(:)
   real(kind= Rkind)                                :: Norm0
   logical, parameter                               :: debug = .true.
   integer                                          :: inb, i1, i3, Ndim
   Integer, allocatable                             :: Ib1(:), Ib2(:), Ib3(:)
   call Calc_index(Ib1=Ib1, Ib2=Ib2, Ib3=Ib3, Basis=psi_dt_1%Basis)
   Ndim = size(psi_dt_1%Basis%tab_basis) - 1
   call Calc_Norm_OF_psi(psi_dt_1,Norm0)
   write (out_unit, *) 'Begin Hagedorn projection <psi|psi>=',Norm0
   !write (out_unit, *) 'out',psi_dt_2%CVec
   !call  Write_VecMat(psi_dt_1%Basis%tab_basis(1)%S, out_unit, 5,  info='S')
   If (Ndim == 1) then
      BBB1(1:Ib1(1), 1:Ib2(1), 1:Ib3(1)) => psi_dt_1%CVec
      BBB2(1:Ib1(1), 1:Ib2(1), 1:Ib3(1)) => psi_dt_2%CVec
      psi_dt_2%CVec(:) = CZERO
      DO i3 = 1, ubound(BBB1, dim=3)
      DO i1 = 1, ubound(BBB1, dim=1)
         BBB2(i1, :, i3) = matmul(psi_dt_1%Basis%tab_basis(1)%S,BBB1(i1, :, i3))
      END DO
      END DO
   else
      allocate (B1(Ib1(1)*Ib2(1)*Ib3(1)))
      BBB1(1:Ib1(1), 1:Ib2(1), 1:Ib3(1)) => psi_dt_1%CVec
      BBB2(1:Ib1(1), 1:Ib2(1), 1:Ib3(1)) => B1
      DO i3 = 1, ubound(BBB1, dim=3)
      DO i1 = 1, ubound(BBB1, dim=1)
         BBB2(i1, :, i3) = matmul(psi_dt_1%Basis%tab_basis(1)%S,BBB1(i1, :, i3))
      END DO
      END DO
      Do inb = 2, Ndim
         allocate (B2(Ib1(inb)*Ib2(inb)*Ib3(inb)))
         BBB1(1:Ib1(inb), 1:Ib2(inb), 1:Ib3(inb)) => B1
         BBB2(1:Ib1(inb), 1:Ib2(inb), 1:Ib3(inb)) => B2
         DO i3 = 1, ubound(BBB1, dim=3)
         DO i1 = 1, ubound(BBB1, dim=1)
            BBB2(i1, :, i3) = matmul( psi_dt_1%Basis%tab_basis(inb)%S,BBB1(i1, :, i3))
         END DO
         END DO
         B1 = B2
         B2 = CZERO
         deallocate (B2)
         if (inb == Ndim) psi_dt_2%CVec = B1
      END DO
   END IF
     call Calc_Norm_OF_psi(psi_dt_2,Norm0)
     deallocate(Ib1, Ib2, Ib3)
     if(allocated(B1)) deallocate(B1)
     if(allocated(B2)) deallocate(B2)
   write (out_unit, *) 'END Hagedorn projection <psi|psi>=',Norm0
   !write (out_unit, *) 'out',psi_dt_2%CVec
END SUBROUTINE 

SUBROUTINE Calc_Norm_OF_Psi(psi, Norm)
   !===========================================================================
   ! PURPOSE:
   !   Compute the norm of the wavefunction psi, either on the grid or in basis.
   !   Optionally prints the norm if debug is enabled.
   !
   ! ARGUMENTS:
   !   psi   : (IN)  Wavefunction structure
   !   Norm  : (INOUT) Norm value of psi
   !===========================================================================

   IMPLICIT NONE

   TYPE(Psi_t), INTENT(IN)          :: psi
   REAL(KIND=Rkind), INTENT(INOUT)  :: Norm
   LOGICAL, PARAMETER                :: debug = .FALSE.  ! Set to .TRUE. to print norm

   !------------------ Compute norm based on psi representation ----------------
   IF (psi%Grid) THEN
      CALL Calc_Norm_OF_PsiGrid(psi, Norm)
   ELSE
      CALL Calc_Norm_OF_PsiBasis(psi, Norm)
   END IF

   !------------------ Debug output if enabled ---------------------------------
   IF (debug) THEN
      WRITE(out_unit,*) '<psi|psi> =', Norm
      FLUSH(out_unit)
   END IF

END SUBROUTINE Calc_Norm_OF_Psi

   SUBROUTINE Calc_Norm_OF_PsiGrid(psi_g, Norm)
      !======================================================================
      ! PURPOSE:
      !   Compute the total norm of a wavefunction (psi_g) defined on a grid,
      !   summing contributions from all electronic states.
      !
      ! ARGUMENTS:
      !   psi_g : (IN)     Psi_t structure containing the wavefunction data.
      !   Norm  : (OUT)    Computed total norm.
      !
      ! NOTES:
      !   - Assumes psi_g%CVec is allocated and contains the complex values
      !     of the wavefunction in grid representation.
      !   - Basis weights are used to compute the weighted norm.
      !   - No memory leaks: all allocated arrays are deallocated at the end.
      !======================================================================
   
      USE QDUtil_m
      LOGICAL, PARAMETER                  :: debug = .FALSE.
      TYPE(Psi_t), INTENT(IN)             :: psi_g
      REAL(KIND=Rkind), INTENT(OUT)       :: Norm
   
      ! Local variables
      COMPLEX(KIND=Rkind), ALLOCATABLE    :: psi_gb(:, :)
      LOGICAL                             :: Endloop_q
      REAL(KIND=Rkind), ALLOCATABLE       :: Norme(:)
      REAL(KIND=Rkind)                    :: WnD
      INTEGER, ALLOCATABLE                :: Tab_iq(:)
      INTEGER                             :: iq, inb, inbe
      INTEGER                             :: nq, nbasis, nstates
   
      !-----------------------------------------------------------
      ! Cache frequently used values for efficiency
      !-----------------------------------------------------------
      nq       = psi_g%Basis%nq
      nbasis   = SIZE(psi_g%Basis%tab_basis)
      nstates  = psi_g%Basis%tab_basis(nbasis)%nb
   
      !-----------------------------------------------------------
      ! Debug output (optional)
      !-----------------------------------------------------------
      IF (debug) THEN
         WRITE (out_unit, *) 'Beginning NormGrid computation'
         FLUSH (out_unit)
      END IF
   
      !-----------------------------------------------------------
      ! Allocate working arrays
      !-----------------------------------------------------------
      ALLOCATE(psi_gb(nq, nstates))
      ALLOCATE(Tab_iq(nbasis - 1))
      ALLOCATE(Norme(nstates))
   
      ! Reshape psi_g%CVec into 2D form: (grid points, nb_states)
      psi_gb(:, :) = RESHAPE(psi_g%CVec, SHAPE=[nq, nstates])
   
      ! Initialize total norm accumulator
      Norm = ZERO
   
      !-----------------------------------------------------------
      ! Loop over electronic states
      !-----------------------------------------------------------
      DO inbe = 1, nstates
   
         Norme(inbe) = ZERO
   
         CALL Init_tab_ind(Tab_iq, psi_g%Basis%NDindexq)
   
         iq = 0
         DO
            iq = iq + 1
            CALL increase_NDindex(Tab_iq, psi_g%Basis%NDindexq, Endloop_q)
            IF (Endloop_q) EXIT
   
            ! Compute the product of the weights for this grid point
            WnD = ONE
            DO inb = 1, nbasis - 1
               WnD = WnD * psi_g%Basis%tab_basis(inb)%w(Tab_iq(inb))
            END DO
   
            ! Accumulate weighted squared amplitude for this state
            Norme(inbe) = Norme(inbe) + REAL(CONJG(psi_gb(iq, inbe)) * psi_gb(iq, inbe) * WnD,kind=Rkind)
   
         END DO
   
         ! Add contribution from this state to the total norm
         Norm = Norm + Norme(inbe)
   
      END DO
   
      ! Final square root after summing over all states
      Norm = SQRT(Norm)
   
      !-----------------------------------------------------------
      ! Free allocated memory to avoid leaks
      !-----------------------------------------------------------
      DEALLOCATE(Tab_iq, Norme, psi_gb)
   
      !-----------------------------------------------------------
      ! Debug output (optional)
      !-----------------------------------------------------------
      IF (debug) THEN
         WRITE (out_unit, *) 'End of NormGrid computation'
         FLUSH (out_unit)
      END IF
   
   END SUBROUTINE Calc_Norm_OF_PsiGrid

   SUBROUTINE Calc_Norm_OF_PsiBasis(psi, Norm)
      !======================================================================
      ! PURPOSE:
      !   Compute the norm of a wavefunction (psi) represented in the basis.
      !   Works automatically for real (RVec) or complex (CVec) representations.
      !
      ! ARGUMENTS:
      !   psi  : (IN)  Psi_t structure containing the wavefunction coefficients
      !   Norm : (OUT) Computed norm of the wavefunction
      !
      ! NOTES:
      !   - If both RVec and CVec are allocated, CVec takes priority.
      !   - Uses Euclidean norm (dot product for real/complex vectors)
      !======================================================================
   
      TYPE(psi_t), INTENT(IN)       :: psi
      REAL(KIND=Rkind), INTENT(OUT) :: Norm
   
      ! Check which vector is allocated
      IF (ALLOCATED(psi%CVec)) THEN
         ! Compute norm for complex coefficients
         Norm = SQRT(REAL(DOT_PRODUCT(psi%CVec, psi%CVec), KIND=Rkind))
      ELSEIF (ALLOCATED(psi%RVec)) THEN
         ! Compute norm for real coefficients
         Norm = SQRT(DOT_PRODUCT(psi%RVec, psi%RVec))
      ELSE
         ! No vector allocated: return zero and print warning
         Norm = ZERO
         PRINT *, 'WARNING: Calc_Norm_OF_PsiBasis: no coefficients allocated in psi'
      END IF
   
   END SUBROUTINE Calc_Norm_OF_PsiBasis

   SUBROUTINE write_psi_basis(psi, print_cplx, t, nio, real_part)
      !======================================================================
      ! PURPOSE:
      !   Write the wavefunction psi in the basis representation to a file
      !   or standard output.
      !
      ! ARGUMENTS:
      !   psi        : (IN)  Psi_t structure containing the wavefunction
      !   print_cplx : (IN)  Logical flag to print complex coefficients
      !   t          : (IN, OPTIONAL) Time value to include in output
      !   nio        : (IN, OPTIONAL) Output unit number
      !   real_part  : (IN)  Logical flag to print only the real part
      !======================================================================
   
      TYPE(psi_t), INTENT(IN), TARGET          :: psi
      REAL(KIND=Rkind), INTENT(IN), OPTIONAL   :: t
      INTEGER, INTENT(IN), OPTIONAL            :: nio
      LOGICAL, INTENT(IN)                      :: print_cplx, real_part
   
      COMPLEX(KIND=Rkind), POINTER             :: psibe(:, :)
      INTEGER                                  :: ib, Ndim
      INTEGER                                  :: outunit_local
      LOGICAL                                   :: has_time
   
      !-----------------------------------------------------------
      ! Determine dimensions and set pointer to complex vector
      !-----------------------------------------------------------
      Ndim = SIZE(psi%Basis%tab_basis)
      psibe(1:psi%Basis%nb, 1:psi%Basis%tab_basis(Ndim)%nb) => psi%CVec
   
      !-----------------------------------------------------------
      ! Determine output unit
      !-----------------------------------------------------------
      IF (PRESENT(nio)) THEN
         outunit_local = nio
      ELSE
         outunit_local = out_unit
      END IF
   
      has_time = PRESENT(t)
   
      !-----------------------------------------------------------
      ! Loop over basis functions
      !-----------------------------------------------------------
      DO ib = 1, psi%Basis%nb
         IF (print_cplx) THEN
            IF (real_part) THEN
               IF (has_time) THEN
                  WRITE(outunit_local, *) ib, REAL(psibe(ib, :), KIND=Rkind), &
                                         AIMAG(psibe(ib, :)), t
               ELSE
                  WRITE(outunit_local, *) ib, REAL(psibe(ib, :), KIND=Rkind), &
                                         AIMAG(psibe(ib, :))
               END IF
            ELSE
               IF (has_time) THEN
                  WRITE(outunit_local, *) t, ib, psibe(ib, :)
               ELSE
                  WRITE(outunit_local, *) ib, psibe(ib, :)
               END IF
            END IF
         ELSE
            ! Print modulus squared |psi|^2
            IF (has_time) THEN
               WRITE(outunit_local, *) t, ib, ABS(psibe(ib, :))**2
            ELSE
               WRITE(outunit_local, *) ib, ABS(psibe(ib, :))**2
            END IF
         END IF
      END DO
   
   END SUBROUTINE write_psi_basis

   SUBROUTINE write_psi_grid(psi, print_cplx, t, nio, real_part)
      !======================================================================
      ! PURPOSE:
      !   Write the wavefunction psi in grid representation to a file or 
      !   standard output.
      !
      ! ARGUMENTS:
      !   psi        : (IN)  Psi_t structure containing the wavefunction
      !   print_cplx : (IN)  Logical flag to print complex coefficients
      !   t          : (IN, OPTIONAL) Time value to include in output
      !   nio        : (IN, OPTIONAL) Output unit number
      !   real_part  : (IN)  Logical flag to print only the real part
      !
      ! NOTES:
      !   - Uses calc_Q_grid to obtain grid coordinates
      !   - Handles complex, real-part-only, or modulus squared output
      !   - Single loop over grid points avoids code duplication
      !======================================================================
   
      TYPE(psi_t), INTENT(IN), TARGET          :: psi
      REAL(KIND=Rkind), INTENT(IN), OPTIONAL   :: t
      INTEGER, INTENT(IN), OPTIONAL            :: nio
      LOGICAL, INTENT(IN)                      :: print_cplx, real_part
   
      COMPLEX(KIND=Rkind), POINTER             :: psige(:, :)
      INTEGER                                  :: iq, Ndim
      REAL(KIND=Rkind), ALLOCATABLE            :: Q(:, :)
      INTEGER                                  :: outunit_local
      LOGICAL                                  :: has_time
   
      !-----------------------------------------------------------
      ! Determine dimensions and set pointer to complex vector
      !-----------------------------------------------------------
      Ndim = SIZE(psi%Basis%tab_basis)
      psige(1:psi%Basis%nq, 1:psi%Basis%tab_basis(Ndim)%nb) => psi%CVec
   
      ! Compute grid coordinates
      CALL calc_Q_grid(Q, psi%Basis)
   
      ! Determine output unit and if time is present
      IF (PRESENT(nio)) THEN
         outunit_local = nio
      ELSE
         outunit_local = out_unit
      END IF
      has_time = PRESENT(t)
   
      !-----------------------------------------------------------
      ! Loop over grid points
      !-----------------------------------------------------------
      DO iq = 1, psi%Basis%nq
   
         IF (print_cplx) THEN
            IF (real_part) THEN
               IF (has_time) THEN
                  WRITE(outunit_local, *) Q(iq, :), REAL(psige(iq, :), KIND=Rkind), &
                                         AIMAG(psige(iq, :)), t
               ELSE
                  WRITE(outunit_local, *) Q(iq, :), REAL(psige(iq, :), KIND=Rkind), &
                                         AIMAG(psige(iq, :))
               END IF
            ELSE
               IF (has_time) THEN
                  WRITE(outunit_local, *) t, Q(iq, :), psige(iq, :)
               ELSE
                  WRITE(outunit_local, *) Q(iq, :), psige(iq, :)
               END IF
            END IF
         ELSE
            ! Print modulus squared |psi|^2
            IF (has_time) THEN
               WRITE(outunit_local, *) t, Q(iq, :), ABS(psige(iq, :))**2
            ELSE
               WRITE(outunit_local, *) Q(iq, :), ABS(psige(iq, :))**2
            END IF
         END IF
   
         ! Optional: add a blank line at end of each grid block (for readability)
         IF (MOD(iq, psi%Basis%tab_basis(1)%nq) == 0) THEN
            WRITE(outunit_local, *)
         END IF
   
      END DO
   
      !-----------------------------------------------------------
      ! Free allocated memory for grid coordinates
      !-----------------------------------------------------------
      DEALLOCATE(Q)
   
   END SUBROUTINE write_psi_grid

   SUBROUTINE psi_init_GWP(psi, Tab_GWP)
      !======================================================================
      ! PURPOSE:
      !   Initialize a wavefunction psi as a Gaussian Wave Packet (GWP)
      !   in grid or basis space, based on the value of psi%grid.
      !
      ! ARGUMENTS:
      !   psi      : (INOUT) Psi_t structure to be initialized
      !   Tab_GWP  : (IN)    Array of GWP_t structures containing GWP parameters
      !======================================================================
  
      USE QDUtil_m
      IMPLICIT NONE
  
      TYPE(psi_t), INTENT(INOUT), TARGET   :: psi
      TYPE(GWP_t), INTENT(IN)              :: Tab_GWP(:)
  
      !------------------ Local variables ---------------------------------------
      TYPE(psi_t), TARGET                  :: psi_temp           ! Temporary wavefunction on grid
      REAL(KIND=Rkind), ALLOCATABLE        :: Q(:, :)            ! Grid coordinates
      INTEGER                              :: iq
      COMPLEX(KIND=Rkind), POINTER         :: gb(:, :)           ! Pointer to temporary/final coefficients
      REAL(KIND=Rkind)                     :: NormB ,normG  ,Re_normB           ! Norm in final representation
      INTEGER                              :: ndim, n_surf, nq, i_surf
      LOGICAL , PARAMETER                  :: debug = .FALSE.
  
      !------------------ Initialization and setup ---------------------------------
      ndim   = SIZE(psi%Basis%tab_basis) - 1
      n_surf = psi%Basis%tab_basis(ndim+1)%nb
      nq     = psi%Basis%nq
      i_surf = Tab_GWP(1)%Elecindex
      
      IF (debug) THEN
          WRITE(out_unit, *) 'Initializing GWP with the following parameters:'
          WRITE(out_unit, *) 'ndim = ', ndim
          WRITE(out_unit, *) 'n_surf = ', n_surf
          IF(psi%grid) THEN
              WRITE(out_unit, *) 'nq = ', nq
          ELSE
              WRITE(out_unit, *) 'nb = ', psi%Basis%nb
          END IF
          WRITE(out_unit, *) 'i_surf = ', i_surf
      END IF

      IF(.NOT. ASSOCIATED(psi%Basis)) THEN
          WRITE(out_unit, *) 'Error: psi%Basis is not allocated in psi_init_GWP.'
          STOP 'Error: psi%Basis is not allocated in psi_init_GWP.'
      ELSE
          IF (debug) THEN
              WRITE(out_unit, *) 'psi%Basis is allocated in psi_init_GWP.'
          END IF
      END IF
      
      IF (.NOT. ALLOCATED(psi%CVec)) THEN
          WRITE(out_unit, *) 'Error: psi%CVec is not allocated in psi_init_GWP.'
          STOP 'Error: psi%CVec is not allocated in psi_init_GWP.'
      ELSE
          IF (debug) THEN
              WRITE(out_unit, *) 'psi%CVec is allocated in psi_init_GWP.'
          END IF
      END IF
      
      
      !------------------ Build Gaussian Wave Packet -----------------------------
      IF (psi%grid) THEN
          ! Case 1: Initialize directly on the grid (psi%grid is .true.)
          IF(debug) WRITE(out_unit, *) 'Initializing GWP directly on the grid (psi%grid = .TRUE.)=:',psi%grid
          
          CALL calc_Q_grid(Q, psi%Basis)
          
          gb(1:nq, 1:n_surf) => psi%CVec
          
          DO iq = 1, nq
              gb(iq, i_surf) = calc_GWP(Tab_GWP(1), Q(iq, :))
          END DO
          
          DEALLOCATE(Q)
  
      ELSE
          ! Case 2: Initialize in the basis representation (psi%grid is .false.)
          IF(debug) WRITE(out_unit, *) 'Initializing GWP on grid, then transforming to basis (psi%grid = .FALSE.)=:',psi%grid
          IF(debug) CALL Write_VecMat(psi%CVec, out_unit, 5, info='psi%CVec')
          
          CALL init_psi(psi_temp, psi%Basis, cplx=.true., grid=.true.)
          CALL calc_Q_grid(Q, psi%Basis)
          
          psi_temp%CVec(:) = CZERO
          gb(1:nq, 1:n_surf) => psi_temp%CVec
          
          DO iq = 1, nq
              gb(iq, i_surf) = calc_GWP(Tab_GWP(1), Q(iq, :))
          END DO
  
          DEALLOCATE(Q)
  
          CALL GridTOBasis_nD_cplx(psi%CVec, psi_temp%CVec, psi%Basis)
  
          CALL dealloc_psi(psi_temp)
      END IF
  
      !------------------ Norm Checking and Renormalization -------------------------
      WRITE(out_unit, *) '------------------- Norm Checking after GWP initialization ----------------'
      WRITE(out_unit, *)
      IF (psi%grid) THEN
        CALL Calc_Norm_Of_Psi(psi, NormG)
        psi%CVec(:) = psi%CVec(:)/NormG
        CALL Calc_Norm_Of_Psi(psi,NormG)
      ELSE
        CALL Calc_Norm_Of_Psi(psi, NormB)
        psi%CVec(:) = psi%CVec(:)/NormB
        CALL Calc_Norm_Of_Psi(psi,NormB)
        psi%CVec(:) = psi%CVec(:)/NormB
        CALL Calc_Norm_Of_Psi(psi,Re_normB)
      END IF
      
      ! Display final state
      IF (psi%grid) THEN
      WRITE(out_unit, *) 'Norm of the initial GWP set on grid  normG := ', NormG  
      ELSE
         WRITE(out_unit, *) 'Norm of the initial GWP set on basis  normB := ', NormB  
         WRITE(out_unit, *) 'Norm of the  renormalized initial GWP set on basis  Re_normB := ', Re_normB  
      END IF
  
      WRITE(out_unit, *) '------------------------------ End Norm Checking --------------------------'
      WRITE(out_unit, *)
  
  END SUBROUTINE psi_init_GWP
  



   SUBROUTINE psi_init_GWPtemp(psi, Tab_GWP)
      !======================================================================
      ! PURPOSE:
      !   Initialize a wavefunction psi as a Gaussian Wave Packet (GWP)
      !   in grid space, then transform it to the basis representation.
      !
      ! ARGUMENTS:
      !   psi      : (INOUT) Psi_t structure to be initialized
      !   Tab_GWP  : (IN)    Array of GWP_t structures containing GWP parameters
      !
      ! NOTES:
      !   - Norm checking is performed before and after Grid->Basis transformation
      !   - Memory for temporary psi0 and grid coordinates Q is properly managed
      !======================================================================
   
      USE QDUtil_m
   
      TYPE(psi_t), INTENT(INOUT), TARGET  :: psi
      TYPE(GWP_t)                         :: Tab_GWP(:)
   
      !------------------ Local variables ---------------------------------------
      TYPE(psi_t), TARGET                  :: psi0               ! temporary wavefunction on grid
      REAL(KIND=Rkind), ALLOCATABLE        :: Q(:, :)            ! grid coordinates
      INTEGER                              :: iq, iGWP
      COMPLEX(KIND=Rkind), POINTER         :: gb(:, :)           ! pointer to psi0 coefficients
      REAL(KIND=Rkind)                      :: NormG, NormB       ! norms before/after basis transformation
      INTEGER                               :: ndim, n_surf, nq, i_surf
   
      !------------------ Initialization ----------------------------------------
      CALL init_psi(psi0, psi%Basis, cplx=.true., grid=.true.)   ! allocate psi0 on grid
      CALL calc_Q_grid(Q, psi%Basis)                             ! compute grid coordinates
   
      ! Zero target wavefunction coefficients
      psi%CVec(:) = CZERO
   
      ndim   = SIZE(psi%Basis%tab_basis) - 1
      n_surf = psi%Basis%tab_basis(ndim+1)%nb
      nq     = psi%Basis%nq
      i_surf = Tab_GWP(1)%Elecindex
   
      ! Point gb to psi0 grid coefficients and zero it
      gb(1:nq, 1:n_surf) => psi0%CVec
      gb(:, :) = CZERO
   
      !------------------ Build Gaussian Wave Packet on grid --------------------
      DO iq = 1, nq
         gb(iq, i_surf) = calc_GWP(Tab_GWP(1), Q(iq, :))
         ! Optional debugging:
         ! WRITE(25,*) Q(iq,:), ABS(gb(iq, Tab_GWP(1)%Elecindex))**2
         ! IF(MOD(iq,10)==0) WRITE(25,*)
      END DO
   
      !------------------ Compute norm on grid ---------------------------------
      CALL Calc_Norm_OF_Psi(psi0, NormG)
      WRITE(out_unit, *) '------------------- Beginning Norm Checking ---------------------------'
      WRITE(out_unit, *)
      WRITE(out_unit, *) 'NormG = ', NormG
   
      !------------------ Transform Grid -> Basis --------------------------------
      CALL GridTOBasis_nD_cplx(psi%CVec, psi0%CVec, psi%Basis)
   
      ! Compute norm in basis representation
      CALL Calc_Norm_OF_Psi(psi, NormB)
      PRINT *, 'NormB = ', NormB
      PRINT *, 'NormB (renormalized) = ', NormB   ! optionally renormalize if needed
   
      WRITE(out_unit, *) '------------------------------ End Norm Checking --------------------------'
      WRITE(out_unit, *)
   
      !------------------ Clean up memory ---------------------------------------
      DEALLOCATE(Q)
      CALL dealloc_psi(psi0)
   
   END SUBROUTINE 


   subroutine test(Basis)
      implicit none
      TYPE(Basis_t), target, intent(in)     :: Basis
      TYPE(psi_t)                            :: psiB, psiG
      real(kind= Rkind)                           :: NormG, NormB
      CALL init_psi(psiB, Basis, cplx=.TRUE., grid=.false.)
      CALL init_psi(psiG, Basis, cplx=.TRUE., grid=.true.)

      psiB%CVec(:) = CZERO
      psiB%CVec(1) = CONE
      call Calc_Norm_OF_Psi(PsiB, NormB)
      psiB%CVec(:) = psiB%CVec(:)/NormB
      call Calc_Norm_OF_Psi(psiB, NormB)
      call BasisTOGrid_nD_cplx(PsiG%CVec, psiB%CVec, Basis)
      !call BasisTOGrid_nE_cplx(PsiG%CVec,psiB%CVec,Basis)
      call Calc_Norm_OF_Psi(psiG, NormG)
      print *, 'NormB = ', NormB
      print *, 'NormG = ', NormG
      print *, '--------------------------------------------------------'
      psiG%CVec(:) = CONE
      ! psiG%CVec(1) = CONE
      call Calc_Norm_OF_Psi(PsiG, NormG)
      psiG%CVec(:) = psiG%CVec(:)/NormG
      call GridTOBasis_nD_cplx(PsiB%CVec, psiG%CVec, Basis)
      call Calc_Norm_OF_Psi(psiB, NormB)
      call Calc_Norm_OF_Psi(psiG, NormG)
      print *, 'NormG = ', NormG
      print *, 'NormB = ', NormB
   end subroutine test



   subroutine Ecrire_psi(psi,nio,t)
      implicit none
      TYPE(psi_t)           , intent(in)     :: psi
      TYPE(psi_t)    ,target                        :: psiG
      real(kind= Rkind), allocatable         :: Q(:, :)
      complex(Kind= Rkind), pointer          :: gb(:, :)
      real(kind= Rkind),intent(in)           :: t
      integer,intent(in)                     :: nio
      integer                                :: iq,ib,nq,nsurf,ndim,nio_2,nq1
      CHARACTER(len=100)                     :: t_string
      CHARACTER(len=100)                     :: name_1,name_2

      ndim = size(psi%Basis%tab_basis)-1
      nsurf = psi%Basis%tab_basis(ndim+1)%nb
      nq = psi%Basis%nq
      nio_2 =nio+1
      nq1 = psi%Basis%tab_basis(1)%nq

      ! Convertir delta_t en chaîne de caractères
      WRITE(t_string, '(F12.6)') t 
      name_1 = 'density_1_t=' // TRIM(ADJUSTL(t_string)) // '.txt'
      name_2 = 'density_2_t=' // TRIM(ADJUSTL(t_string)) // '.txt'
      
      open (unit=nio, file=name_1)
      open (unit=nio_2, file=name_2)

      call calc_Q_grid(Q, psi%Basis)
      CALL init_psi(psiG, psi%Basis, cplx=.TRUE., grid=.true.)
      call BasisTOGrid_nD_cplx(psiG%CVec, psi%CVec, psi%Basis)

      gb(1:nq, 1:nsurf) => psiG%CVec
      DO iq = 1,nq
       write(nio,*)  Q(iq,:),    abs( gb(iq,1))**2
       write(nio_2,*)  Q(iq,:),   abs( gb(iq,2))**2
       if (mod(iq, nq1 ) == 0) write(nio,*)
       if (mod(iq, nq1 ) == 0)  write(nio_2,*)
      END DO

      
      deallocate (Q)
      CALL dealloc_psi(psiG)
   end subroutine 





end module psi_m
