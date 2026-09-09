MODULE SmolyakSparseGridMetadata_m
  !! ===========================================================================
  !! SmolyakSparseGridMetadata_m
  !! ---------------------------------------------------------------------------
  !! Purpose
  !!   Provide a compact, explicit and academically traceable data model for the
  !!   Smolyak sparse-grid / sparse-basis representation attached to the user's
  !!   multidimensional Basis_t structure.
  !!
  !! Philosophy
  !!   1. Keep the dense direct-product basis as the reference mathematical
  !!      object already implemented in the current code base.
  !!   2. Build a sparse nuclear basis by selecting the direct-product indices
  !!      that belong to at least one Smolyak tensor block.
  !!   3. Keep the electronic-state expansion explicit and untouched: the sparse
  !!      selection is applied only to the nuclear part, then replicated for all
  !!      electronic surfaces.
  !!   4. Allow the rest of the code to convert sparse <-> dense basis vectors
  !!      while reusing the existing robust grid algebra.
  !!
  !! Conventions
  !!      ell_i = s_i + 1
  !!      q     = sum_i (ell_i - 1)
  !!      qmin  = max(0, K - D)
  !!      c(ell)= (-1)^(K-1-q) * C(D-1, K-1-q)
  !! ===========================================================================
  USE QDUtil_m,       ONLY: Rkind, ZERO, CZERO, out_unit, binomial
  USE NDindex_m,      ONLY: NDindex_t, Init_tab_ind, increase_NDindex
  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: BASIS_LAYOUT_DIRECT_PRODUCT    = 0
  INTEGER, PARAMETER, PUBLIC :: BASIS_LAYOUT_SMOLYAK_SPARSE    = 1


! ==========================================================================
! The smolyak term descriptor type
!
! This type contains the information about a smolyak term.
!
! Parameters:
!   term_id : The term id
!   q : The q norm
!   coefficient : The coefficient
!   ell : The level ell
!
! Returns:
!   None
!
! ==========================================================================
  TYPE, PUBLIC :: SmolyakTermDescriptor_t
     INTEGER :: term_id = 0
     INTEGER :: q         = 0
     INTEGER :: coefficient = 0
     INTEGER, ALLOCATABLE :: ell(:)
  END TYPE SmolyakTermDescriptor_t

! ==========================================================================
! The smolyak sparse grid descriptor type
!
! This type contains the information about a smolyak sparse grid.
!
! Parameters:
!   layout_kind : The layout kind
!   is_initialized : Whether the descriptor is initialized
!   ndim : The number of dimensions
!   nsurf : The number of electronic surfaces
!   order_level : The Smolyak level used to build a term family
!   qmin : The minimum q value
!   nterms : The number of terms
!   rule_shift : The rule shift
!   rule_slope : The rule slope
!   active_nuclear_basis_index : The active nuclear basis index
!   active_nuclear_grid_index : The active nuclear grid index
!   nb : Global tensor-product basis size
!   nb_SG : Active Smolyak sparse basis size
!
! Returns:
!   None
!
! ==========================================================================
  TYPE, PUBLIC :: SmolyakSparseGridDescriptor_t
     INTEGER :: layout_kind = BASIS_LAYOUT_DIRECT_PRODUCT   ! To choose between direct product and sparse layout
     LOGICAL :: is_initialized = .FALSE.
     INTEGER :: ndim = 0                                    ! Number of dimensions
     INTEGER :: nsurf = 0                                   ! Number of electronic surfaces
     INTEGER :: order_LB = 0                                ! Smolyak order for the basis
     INTEGER :: order_LG = 0                                ! Smolyak order for the grid
     INTEGER :: qmin_B = 0                                  ! Minimum q for the basis family
     INTEGER :: qmax_B = -1                                 ! Maximum q for the basis family
     INTEGER :: qmin_G = 0                                  ! Minimum q for the grid family
     INTEGER :: qmax_G = -1                                 ! Maximum q for the grid family
     INTEGER :: nterms_b = 0                                ! Number of basis terms (LB)
     INTEGER :: nterms_g = 0                                ! Number of grid terms (LG)
     INTEGER :: nb = 0                 ! Total size of the full direct-product basis
     INTEGER :: nq  = 0                 ! Total size of the full direct-product grid
     INTEGER :: nb_SG = 0               ! Size of the sparse nuclear basis
     INTEGER :: nq_SG  = 0               ! Size of the sparse nuclear grid
     INTEGER, ALLOCATABLE :: rule_shift_B(:)                 ! Basis-rule shift A_b_smol
     INTEGER, ALLOCATABLE :: rule_slope_B(:)                 ! Basis-rule slope B_b_smol
     INTEGER, ALLOCATABLE :: rule_shift_G(:)                 ! Grid-rule shift A_g_smol
     INTEGER, ALLOCATABLE :: rule_slope_G(:)                 ! Grid-rule slope B_g_smol
     INTEGER, ALLOCATABLE :: basis_cap_by_level(:,:)         ! Basis capacity by level
     INTEGER, ALLOCATABLE :: grid_cap_by_level(:,:)          ! Grid capacity by level
     TYPE(SmolyakTermDescriptor_t), ALLOCATABLE :: term_B(:) ! Basis terms (LB)
     TYPE(SmolyakTermDescriptor_t), ALLOCATABLE :: term_G(:) ! Grid terms (LG)
     LOGICAL, ALLOCATABLE :: active_nuclear_basis_mask(:)    ! Active basis mask
     LOGICAL, ALLOCATABLE :: active_nuclear_grid_mask(:)     ! Active grid mask
     INTEGER, ALLOCATABLE :: active_nuclear_basis_index(:)   ! Active basis indices
     INTEGER, ALLOCATABLE :: active_nuclear_grid_index(:)    ! Active grid indices
     ! Canonical SG4 storage: each term owns an independent packed block.
     INTEGER, ALLOCATABLE :: term_basis_size(:)
     INTEGER, ALLOCATABLE :: term_Grid_size(:)
     INTEGER, ALLOCATABLE :: term_basis_offset(:)
     INTEGER, ALLOCATABLE :: term_Grid_offset(:)
  END TYPE SmolyakSparseGridDescriptor_t

  PUBLIC :: InitializeSmolyakSparseGridDescriptor
  PUBLIC :: DestroySmolyakSparseGridDescriptor
  PUBLIC :: GetEffectiveTotalBasisSize
  PUBLIC :: ExpandSparseBasisToFullBasis
  PUBLIC :: RestrictFullBasisToSparseBasis
  PUBLIC :: PrintSmolyakSparseGridDescriptor
  PUBLIC :: GetSmolyakPackedBasisSize, GetSmolyakPackedGridSize
  PUBLIC :: PackFullBasisToSmolyakTerms, UnpackSmolyakTermsToFullBasis
  PUBLIC :: CombineSmolyakBlocksToFull
  PUBLIC :: DotSmolyakBasis

CONTAINS

! ====================================================================================
! The Destructor for the SmolyakSparseGridDescriptor_t type
! This subroutine deallocates all allocatable components of the descriptor
! Parameters:
!   Descriptor: The descriptor to be destroyed
  SUBROUTINE DestroySmolyakSparseGridDescriptor(Descriptor)
    TYPE(SmolyakSparseGridDescriptor_t), INTENT(INOUT) :: Descriptor
    INTEGER :: it
    Descriptor%layout_kind = BASIS_LAYOUT_DIRECT_PRODUCT
    Descriptor%is_initialized = .FALSE.
    Descriptor%ndim = 0.      ! Number of dimensions
    Descriptor%nsurf = 0      ! Number of electronic surfaces
    Descriptor%order_LB = 0
    Descriptor%order_LG = 0
    Descriptor%qmin_B = 0
    Descriptor%qmax_B = -1
    Descriptor%qmin_G = 0
    Descriptor%qmax_G = -1
    Descriptor%nterms_b = 0
    Descriptor%nterms_g = 0
    Descriptor%nb = 0    ! Total size of the full direct-product basis
    Descriptor%nq  = 0    ! Total size of the full direct-product grid
    Descriptor%nb_SG = 0  ! Size of the sparse nuclear basis
    Descriptor%nq_SG  = 0  ! Size of the sparse nuclear grid
    IF (ALLOCATED(Descriptor%rule_shift_B)) DEALLOCATE(Descriptor%rule_shift_B)
    IF (ALLOCATED(Descriptor%rule_slope_B)) DEALLOCATE(Descriptor%rule_slope_B)
    IF (ALLOCATED(Descriptor%rule_shift_G)) DEALLOCATE(Descriptor%rule_shift_G)
    IF (ALLOCATED(Descriptor%rule_slope_G)) DEALLOCATE(Descriptor%rule_slope_G)
    IF (ALLOCATED(Descriptor%basis_cap_by_level)) DEALLOCATE(Descriptor%basis_cap_by_level)
    IF (ALLOCATED(Descriptor%grid_cap_by_level))  DEALLOCATE(Descriptor%grid_cap_by_level)
    IF (ALLOCATED(Descriptor%term_B)) THEN
       DO it = 1, SIZE(Descriptor%term_B)
          IF (ALLOCATED(Descriptor%term_B(it)%ell)) DEALLOCATE(Descriptor%term_B(it)%ell)
       END DO
       DEALLOCATE(Descriptor%term_B)
    END IF
    IF (ALLOCATED(Descriptor%term_G)) THEN
       DO it = 1, SIZE(Descriptor%term_G)
          IF (ALLOCATED(Descriptor%term_G(it)%ell)) DEALLOCATE(Descriptor%term_G(it)%ell)
       END DO
       DEALLOCATE(Descriptor%term_G)
    END IF
    IF (ALLOCATED(Descriptor%active_nuclear_basis_mask))  DEALLOCATE(Descriptor%active_nuclear_basis_mask)
    IF (ALLOCATED(Descriptor%active_nuclear_grid_mask))   DEALLOCATE(Descriptor%active_nuclear_grid_mask)
    IF (ALLOCATED(Descriptor%active_nuclear_basis_index)) DEALLOCATE(Descriptor%active_nuclear_basis_index)
    IF (ALLOCATED(Descriptor%active_nuclear_grid_index))  DEALLOCATE(Descriptor%active_nuclear_grid_index)
    IF (ALLOCATED(Descriptor%term_basis_size)) DEALLOCATE(Descriptor%term_basis_size)
    IF (ALLOCATED(Descriptor%term_grid_size)) DEALLOCATE(Descriptor%term_grid_size)
    IF (ALLOCATED(Descriptor%term_basis_offset)) DEALLOCATE(Descriptor%term_basis_offset)
    IF (ALLOCATED(Descriptor%term_grid_offset)) DEALLOCATE(Descriptor%term_grid_offset)
  END SUBROUTINE DestroySmolyakSparseGridDescriptor



  !=================================================================
  !The Smolyak term count function
  !
  ! This function calculates the number of terms in a Smolyak sparse grid
  ! for a given number of dimensions and order.
  !
  ! Parameters:
  !   ndim : Number of dimensions
  !   order_level : Smolyak level (L)
  !
  ! Returns:
  !   The number of terms in the Smolyak sparse grid
  !
  INTEGER FUNCTION SmolyakTermCount(ndim, order_level) RESULT(nterms)
    INTEGER, INTENT(IN) :: ndim, order_level

    INTEGER :: q
    nterms = 0
    DO q = MAX(0, order_level - ndim), order_level - 1
       nterms = nterms + NINT(binomial(q + ndim - 1, ndim - 1))
    END DO
  END FUNCTION SmolyakTermCount

  !=================================================================
  !The Smolyak coefficient function
  !
  ! This function calculates the Smolyak coefficient(D_l) for a given
  ! dimension, order, and q value.
  !
  ! Parameters:
  !   ndim : Number of dimensions
  !   order_level : Smolyak level (L)
  !   q : Q value
  !
  ! Returns:
  !   The Smolyak coefficient
  !
  INTEGER FUNCTION SmolyakCoefficient(ndim, order_level, q) RESULT(coeff)
    INTEGER, INTENT(IN) :: ndim, order_level, q
    INTEGER :: p
    p = order_level - 1 - q
    IF (p < 0 .OR. p > ndim - 1) THEN
       coeff = 0
    ELSE
       IF (ndim == 1 .AND. p == 0) THEN
         coeff = 1
      ELSE
         coeff = NINT(binomial(ndim - 1, p))
      END IF
       IF (MOD(p, 2) /= 0) coeff = -coeff
    END IF
  END FUNCTION SmolyakCoefficient

!======================================================================================
! The Q composition enumeration function
!
! This function recursively enumerates all compositions of a given sum q
! into ndim parts.
!
! Parameters:
!   ndim : Number of dimensions
!   q : Sum to compose
!   position : Current position in the composition
!   Work : Working array for the current composition
!   Table : Table to store all compositions
!   offset : Current offset in the table
!
! Returns:
!   None (modifies Table and offset in place)
!
  RECURSIVE SUBROUTINE EnumerateQCompositions(ndim, q, position, Work, Table, offset)
    INTEGER, INTENT(IN)    :: ndim, q, position
    INTEGER, INTENT(INOUT) :: Work(:)
    INTEGER, INTENT(INOUT) :: Table(:,:)
    INTEGER, INTENT(INOUT) :: offset
    INTEGER :: v
    IF (position == ndim) THEN
       Work(position) = q
       offset = offset + 1
       Table(:, offset) = Work(:)
       RETURN
    END IF
    DO v = 0, q
       Work(position) = v
       CALL EnumerateQCompositions(ndim, q-v, position+1, Work, Table, offset)
    END DO
  END SUBROUTINE EnumerateQCompositions

!======================================================================================
! The point belongs to any Smolyak block function
!
! This function checks if a given point belongs to any Smolyak block.
!
! Parameters:
!   index_vector : Index vector of the point
!   cap_table : Table of cap values
!   Descriptor : Smolyak sparse grid descriptor
!
! Returns:
!   .TRUE. if the point belongs to any Smolyak block, .FALSE. otherwise
!
  PURE LOGICAL FUNCTION PointBelongsToAnySmolyakBlock(index_vector, cap_table, Descriptor, use_grid) RESULT(is_active)
    INTEGER, INTENT(IN) :: index_vector(:)
    INTEGER, INTENT(IN) :: cap_table(:,:)
    TYPE(SmolyakSparseGridDescriptor_t), INTENT(IN) :: Descriptor
    LOGICAL, INTENT(IN), OPTIONAL :: use_grid
    INTEGER :: it, idim, nterms_local
    LOGICAL :: inside, grid_mode
    grid_mode = .FALSE.; IF (PRESENT(use_grid)) grid_mode = use_grid
    nterms_local = Descriptor%nterms_b
    is_active = .FALSE.
    DO it = 1, MERGE(Descriptor%nterms_g, Descriptor%nterms_b, grid_mode)
       inside = .TRUE.
       DO idim = 1, SIZE(index_vector)
          IF (grid_mode) THEN
             IF (index_vector(idim) > cap_table(idim, Descriptor%term_G(it)%ell(idim))) inside = .FALSE.
          ELSE
             IF (index_vector(idim) > cap_table(idim, Descriptor%term_B(it)%ell(idim))) inside = .FALSE.
          END IF
          IF (.NOT. inside) EXIT
       END DO
       IF (inside) THEN; is_active = .TRUE.; RETURN; END IF
    END DO
  END FUNCTION PointBelongsToAnySmolyakBlock

!======================================================================================
! The build Smolyak term table function
!
! This function builds the Smolyak term table for a given descriptor.
!
! Parameters:
!   Descriptor : Smolyak sparse grid descriptor
!
! Returns:
!   None (modifies Descriptor in place)
!
  SUBROUTINE BuildSmolyakTermTable(Descriptor)
    TYPE(SmolyakSparseGridDescriptor_t), INTENT(INOUT) :: Descriptor
    CALL BuildSmolyakTermTableForOrder(Descriptor%term_B, Descriptor%nterms_b, Descriptor%ndim, Descriptor%order_LB)
    CALL BuildSmolyakTermTableForOrder(Descriptor%term_G, Descriptor%nterms_g, Descriptor%ndim, Descriptor%order_LG)
    Descriptor%qmin_B = MAX(0, Descriptor%order_LB - Descriptor%ndim)
    Descriptor%qmax_B = Descriptor%order_LB - 1
    Descriptor%qmin_G = MAX(0, Descriptor%order_LG - Descriptor%ndim)
    Descriptor%qmax_G = Descriptor%order_LG - 1
  END SUBROUTINE BuildSmolyakTermTable

  SUBROUTINE BuildSmolyakTermTableForOrder(terms, nterms, ndim, order_level)
    TYPE(SmolyakTermDescriptor_t), ALLOCATABLE, INTENT(OUT) :: terms(:)
    INTEGER, INTENT(OUT) :: nterms
    INTEGER, INTENT(IN) :: ndim, order_level
    INTEGER :: q, nq, offset_local, coeff, it
    INTEGER, ALLOCATABLE :: work(:), tmp(:,:)
    nterms = SmolyakTermCount(ndim, order_level)
    ALLOCATE(terms(nterms), work(ndim))
    it = 0
    DO q = MAX(0, order_level-ndim), order_level-1
       nq = NINT(binomial(q+ndim-1, ndim-1))
       ALLOCATE(tmp(ndim,nq)); work=0; offset_local=0
       CALL EnumerateQCompositions(ndim,q,1,work,tmp,offset_local)
       coeff = SmolyakCoefficient(ndim,order_level,q)
       DO offset_local=1,nq
          it=it+1
          terms(it)%term_id=it; terms(it)%q=q; terms(it)%coefficient=coeff
          ALLOCATE(terms(it)%ell(ndim))
          terms(it)%ell=tmp(:,offset_local)+1
       END DO
       DEALLOCATE(tmp)
    END DO
    DEALLOCATE(work)
  END SUBROUTINE BuildSmolyakTermTableForOrder

!======================================================================================
! The build Smolyak cap tables function
!
! This function builds the Smolyak cap tables for a given descriptor.
!
! Parameters:
!   Descriptor : Smolyak sparse grid descriptor
!   nuclear_nb_1d : Nuclear basis size for each dimension
!   nuclear_nq_1d : Nuclear grid size for each dimension
!
! Returns:
!   None (modifies Descriptor in place)
!
  SUBROUTINE BuildSmolyakCapTables(Descriptor, nuclear_nb_1d, nuclear_nq_1d)
    TYPE(SmolyakSparseGridDescriptor_t), INTENT(INOUT) :: Descriptor
    INTEGER, INTENT(IN) :: nuclear_nb_1d(:), nuclear_nq_1d(:)
    INTEGER :: idim, ilevel, maxlevel
    maxlevel = MAX(Descriptor%order_LB, Descriptor%order_LG)
    ALLOCATE(Descriptor%basis_cap_by_level(Descriptor%ndim, Descriptor%order_LB))
    ALLOCATE(Descriptor%grid_cap_by_level (Descriptor%ndim, Descriptor%order_LG))
    DO idim = 1, Descriptor%ndim
       DO ilevel = 1, Descriptor%order_LB
          Descriptor%basis_cap_by_level(idim, ilevel) = MIN(nuclear_nb_1d(idim), Descriptor%rule_shift_B(idim) + Descriptor%rule_slope_B(idim) * ilevel)
          Descriptor%basis_cap_by_level(idim, ilevel) = MAX(1, Descriptor%basis_cap_by_level(idim, ilevel))
       END DO
       DO ilevel = 1, Descriptor%order_LG
          Descriptor%grid_cap_by_level(idim, ilevel) = MIN(nuclear_nq_1d(idim), Descriptor%rule_shift_G(idim) + Descriptor%rule_slope_G(idim) * ilevel)
          Descriptor%grid_cap_by_level(idim, ilevel) = MAX(1, Descriptor%grid_cap_by_level(idim, ilevel))
       END DO
    END DO
  END SUBROUTINE BuildSmolyakCapTables

  SUBROUTINE BuildSmolyakTermStorage(Descriptor)
    TYPE(SmolyakSparseGridDescriptor_t), INTENT(INOUT) :: Descriptor
    INTEGER :: it, i, nb, nq, boff, goff
    ALLOCATE(Descriptor%term_basis_size(Descriptor%nterms_b))
    ALLOCATE(Descriptor%term_grid_size(Descriptor%nterms_g))
    ALLOCATE(Descriptor%term_basis_offset(Descriptor%nterms_b))
    ALLOCATE(Descriptor%term_grid_offset(Descriptor%nterms_g))
    boff=1; goff=1
    DO it=1,Descriptor%nterms_b
       nb=1
       DO i=1,Descriptor%ndim
          nb=nb*Descriptor%basis_cap_by_level(i,Descriptor%term_B(it)%ell(i))
       END DO
       Descriptor%term_basis_size(it)=nb*MAX(1,Descriptor%nsurf)
       Descriptor%term_basis_offset(it)=boff
       boff=boff+Descriptor%term_basis_size(it)
    END DO
    DO it=1,Descriptor%nterms_g
       nq=1
       DO i=1,Descriptor%ndim
          nq=nq*Descriptor%grid_cap_by_level(i,Descriptor%term_G(it)%ell(i))
       END DO
       Descriptor%term_grid_size(it)=nq*MAX(1,Descriptor%nsurf)
       Descriptor%term_grid_offset(it)=goff
       goff=goff+Descriptor%term_grid_size(it)
    END DO
  END SUBROUTINE BuildSmolyakTermStorage

  SUBROUTINE BuildActiveIndicesDirect(indices, mask, nfull, caps_by_level, terms, nterms, ndim)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: indices(:)
    LOGICAL, ALLOCATABLE, INTENT(OUT) :: mask(:)
    INTEGER, INTENT(IN) :: nfull, caps_by_level(:,:), nterms, ndim
    TYPE(SmolyakTermDescriptor_t), INTENT(IN) :: terms(:)
    INTEGER, ALLOCATABLE :: candidates(:), caps(:), dims(:), idx(:)
    INTEGER :: it, idim, iloc, nloc, ncand, iglob, j, nuniq, key, pos
    ALLOCATE(mask(nfull)); mask=.FALSE.
    ncand=0
    DO it=1,nterms
       nloc=1
       DO idim=1,ndim
          nloc=nloc*caps_by_level(idim,terms(it)%ell(idim))
       END DO
       ncand=ncand+nloc
    END DO
    ALLOCATE(candidates(MAX(1,ncand)), caps(ndim), dims(ndim), idx(ndim))
    dims=caps_by_level(:,SIZE(caps_by_level,2)); ncand=0
    DO it=1,nterms
       DO idim=1,ndim; caps(idim)=caps_by_level(idim,terms(it)%ell(idim)); END DO
       nloc=PRODUCT(caps)
       DO iloc=1,nloc
          CALL DecodeSmolyakLinear(iloc,caps,idx)
          iglob=idx(1)
          DO j=2,ndim; iglob=iglob+(idx(j)-1)*PRODUCT(dims(1:j-1)); END DO
          ncand=ncand+1; candidates(ncand)=iglob
       END DO
    END DO
    ! Insertion sort is deterministic and keeps the Fortran linear ordering.
    DO j=2,ncand
       key=candidates(j); pos=j-1
        DO WHILE (pos>=1)
           IF (candidates(pos)<=key) EXIT
           candidates(pos+1)=candidates(pos)
           pos=pos-1
        END DO
        candidates(pos+1)=key
    END DO
    nuniq=0
    DO j=1,ncand
       IF (j==1) THEN
          nuniq=nuniq+1
          candidates(nuniq)=candidates(j)
          mask(candidates(j))=.TRUE.
       ELSEIF (candidates(j)/=candidates(j-1)) THEN
          nuniq=nuniq+1
          candidates(nuniq)=candidates(j)
          mask(candidates(j))=.TRUE.
       END IF
    END DO
    ALLOCATE(indices(nuniq)); IF (nuniq>0) indices=candidates(1:nuniq)
    DEALLOCATE(candidates,caps,dims,idx)
  END SUBROUTINE BuildActiveIndicesDirect

!======================================================================================
! The build Smolyak active masks function
!
! This function builds the Smolyak active masks for a given descriptor.
!
! Parameters:
!   Descriptor : Smolyak sparse grid descriptor
!   NuclearBasisIndex : Nuclear basis index
!   NuclearGridIndex : Nuclear grid index
!
! Returns:
!   None (modifies Descriptor in place)
!
  SUBROUTINE BuildSmolyakActiveMasks(Descriptor, NuclearBasisIndex, NuclearGridIndex)
    TYPE(SmolyakSparseGridDescriptor_t), INTENT(INOUT) :: Descriptor
    TYPE(NDindex_t), INTENT(IN) :: NuclearBasisIndex, NuclearGridIndex
    CALL BuildActiveIndicesDirect(Descriptor%active_nuclear_basis_index, Descriptor%active_nuclear_basis_mask, &
         Descriptor%nb, Descriptor%basis_cap_by_level, Descriptor%term_B, Descriptor%nterms_b, Descriptor%ndim)
    CALL BuildActiveIndicesDirect(Descriptor%active_nuclear_grid_index, Descriptor%active_nuclear_grid_mask, &
         Descriptor%nq, Descriptor%grid_cap_by_level, Descriptor%term_G, Descriptor%nterms_g, Descriptor%ndim)
    Descriptor%nb_SG = SIZE(Descriptor%active_nuclear_basis_index)
    Descriptor%nq_SG = SIZE(Descriptor%active_nuclear_grid_index)
  END SUBROUTINE BuildSmolyakActiveMasks

!======================================================================================
! The initialize Smolyak sparse grid descriptor function
!
! This function initializes the Smolyak sparse grid descriptor.
!
! Parameters:
!   Descriptor : Smolyak sparse grid descriptor
!   ndim : Number of dimensions
!   nsurf : Number of surfaces
!   order_level : Smolyak level (L)
!   rule_shift : Rule shift for each dimension
!   rule_slope : Rule slope for each dimension
!   NuclearBasisIndex : Nuclear basis index
!   NuclearGridIndex : Nuclear grid index
!   nuclear_nb_1d : Nuclear basis size for each dimension
!   nuclear_nq_1d : Nuclear grid size for each dimension
!   nb : Global tensor-product basis size
!   nq : Global tensor-product grid size
!
! Returns:
!   None (modifies Descriptor in place)
!
  SUBROUTINE InitializeSmolyakSparseGridDescriptor(Descriptor, ndim, nsurf, LB, LG, rule_shift_B, rule_slope_B, rule_shift_G, &
   &rule_slope_G, NuclearBasisIndex, NuclearGridIndex, nuclear_nb_1d, nuclear_nq_1d, nb, nq)
    TYPE(SmolyakSparseGridDescriptor_t), INTENT(INOUT) :: Descriptor
    INTEGER, INTENT(IN) :: ndim, nsurf, LB, LG
    INTEGER, INTENT(IN) :: rule_shift_B(:), rule_slope_B(:), rule_shift_G(:), rule_slope_G(:)
    TYPE(NDindex_t), INTENT(IN) :: NuclearBasisIndex, NuclearGridIndex
    INTEGER, INTENT(IN) :: nuclear_nb_1d(:), nuclear_nq_1d(:)
    INTEGER, INTENT(IN) :: nb, nq

    CALL DestroySmolyakSparseGridDescriptor(Descriptor)
    IF (ndim < 1) RETURN
    IF (LB < 1 .OR. LG < 1) RETURN

    Descriptor%layout_kind = BASIS_LAYOUT_SMOLYAK_SPARSE
    Descriptor%is_initialized = .TRUE.
    Descriptor%ndim = ndim
    Descriptor%nsurf = nsurf
    Descriptor%order_LB = LB
    Descriptor%order_LG = LG
    Descriptor%nb = nb
    Descriptor%nq  = nq

    ALLOCATE(Descriptor%rule_shift_B(ndim), Descriptor%rule_slope_B(ndim), &
             Descriptor%rule_shift_G(ndim), Descriptor%rule_slope_G(ndim))
    Descriptor%rule_shift_B = rule_shift_B
    Descriptor%rule_slope_B = rule_slope_B
    Descriptor%rule_shift_G = rule_shift_G
    Descriptor%rule_slope_G = rule_slope_G

    CALL BuildSmolyakTermTable(Descriptor)
    CALL BuildSmolyakCapTables(Descriptor, nuclear_nb_1d, nuclear_nq_1d)
    CALL BuildSmolyakTermStorage(Descriptor)
    CALL BuildSmolyakActiveMasks(Descriptor, NuclearBasisIndex, NuclearGridIndex)
  END SUBROUTINE InitializeSmolyakSparseGridDescriptor

  !=========================================================================
  ! The get effective total basis size function
  !
  ! This function returns the effective total basis size for a given descriptor.
  !
  ! Parameters:
  !   Descriptor : Smolyak sparse grid descriptor
  !
  ! Returns:
  !   The effective total basis size
  !
  PURE INTEGER FUNCTION GetSmolyakPackedBasisSize(Descriptor) RESULT(ncoeff)
    TYPE(SmolyakSparseGridDescriptor_t), INTENT(IN) :: Descriptor
    IF (Descriptor%layout_kind == BASIS_LAYOUT_SMOLYAK_SPARSE .AND. Descriptor%is_initialized .AND. &
        ALLOCATED(Descriptor%active_nuclear_basis_index)) THEN
       ! Canonical sparse representation: one value per active nuclear index
       ! and per electronic surface.  Term blocks remain metadata only.
       ncoeff=SIZE(Descriptor%active_nuclear_basis_index)*MAX(1,Descriptor%nsurf)
    ELSE
       ncoeff=Descriptor%nb*MAX(1,Descriptor%nsurf)
    END IF
  END FUNCTION GetSmolyakPackedBasisSize

  PURE INTEGER FUNCTION GetSmolyakPackedGridSize(Descriptor) RESULT(ncoeff)
    TYPE(SmolyakSparseGridDescriptor_t), INTENT(IN) :: Descriptor
    IF (Descriptor%layout_kind == BASIS_LAYOUT_SMOLYAK_SPARSE .AND. Descriptor%is_initialized .AND. ALLOCATED(Descriptor%term_grid_size)) THEN
       ncoeff=SUM(Descriptor%term_grid_size)
    ELSE
       ncoeff=Descriptor%nq*MAX(1,Descriptor%nsurf)
    END IF
  END FUNCTION GetSmolyakPackedGridSize

  PURE INTEGER FUNCTION GetEffectiveTotalBasisSize(Descriptor) RESULT(ncoeff)
    TYPE(SmolyakSparseGridDescriptor_t), INTENT(IN) :: Descriptor
    ncoeff=GetSmolyakPackedBasisSize(Descriptor)
  END FUNCTION GetEffectiveTotalBasisSize

  SUBROUTINE DecodeSmolyakLinear(linear, caps, idx)
    INTEGER, INTENT(IN) :: linear, caps(:)
    INTEGER, INTENT(OUT) :: idx(:)
    INTEGER :: r, i
    r=linear-1
    DO i=1,SIZE(caps)
       idx(i)=MOD(r,caps(i))+1
       r=r/caps(i)
    END DO
  END SUBROUTINE DecodeSmolyakLinear

  SUBROUTINE PackFullBasisToSmolyakTerms(packed, full, Descriptor)
    COMPLEX(KIND=Rkind), INTENT(OUT) :: packed(:)
    COMPLEX(KIND=Rkind), INTENT(IN) :: full(:)
    TYPE(SmolyakSparseGridDescriptor_t), INTENT(IN) :: Descriptor
    INTEGER :: is, off_dense, off_sparse, nactive
    nactive=SIZE(Descriptor%active_nuclear_basis_index)
    IF (SIZE(packed) /= nactive*MAX(1,Descriptor%nsurf)) STOP 'PackSmolyak: wrong sparse size'
    IF (SIZE(full) /= Descriptor%nb*MAX(1,Descriptor%nsurf)) STOP 'PackSmolyak: wrong full size'
    DO is=1,MAX(1,Descriptor%nsurf)
       off_dense=(is-1)*Descriptor%nb
       off_sparse=(is-1)*nactive
       packed(off_sparse+1:off_sparse+nactive)= &
          full(off_dense+Descriptor%active_nuclear_basis_index(:))
    END DO
  END SUBROUTINE PackFullBasisToSmolyakTerms

  SUBROUTINE CombineSmolyakBlocksToFull(full, packed, Descriptor)
    COMPLEX(KIND=Rkind), INTENT(OUT) :: full(:)
    COMPLEX(KIND=Rkind), INTENT(IN) :: packed(:)
    TYPE(SmolyakSparseGridDescriptor_t), INTENT(IN) :: Descriptor
    INTEGER :: is, off_dense, off_sparse, nactive
    nactive=SIZE(Descriptor%active_nuclear_basis_index)
    IF (SIZE(packed) /= nactive*MAX(1,Descriptor%nsurf)) STOP 'UnpackSmolyak: wrong sparse size'
    IF (SIZE(full) /= Descriptor%nb*MAX(1,Descriptor%nsurf)) STOP 'UnpackSmolyak: wrong full size'
    full=CZERO
    DO is=1,MAX(1,Descriptor%nsurf)
       off_dense=(is-1)*Descriptor%nb
       off_sparse=(is-1)*nactive
       full(off_dense+Descriptor%active_nuclear_basis_index(:))= &
          packed(off_sparse+1:off_sparse+nactive)
    END DO
  END SUBROUTINE CombineSmolyakBlocksToFull

  SUBROUTINE UnpackSmolyakTermsToFullBasis(full, packed, Descriptor)
    COMPLEX(KIND=Rkind), INTENT(OUT) :: full(:)
    COMPLEX(KIND=Rkind), INTENT(IN)  :: packed(:)
    TYPE(SmolyakSparseGridDescriptor_t), INTENT(IN) :: Descriptor
    CALL CombineSmolyakBlocksToFull(full, packed, Descriptor)
  END SUBROUTINE UnpackSmolyakTermsToFullBasis

  FUNCTION DotSmolyakBasis(a, b, Descriptor) RESULT(value)
    COMPLEX(KIND=Rkind), INTENT(IN) :: a(:), b(:)
    TYPE(SmolyakSparseGridDescriptor_t), INTENT(IN) :: Descriptor
    COMPLEX(KIND=Rkind) :: value
    COMPLEX(KIND=Rkind), ALLOCATABLE :: fa(:), fb(:)
    ALLOCATE(fa(Descriptor%nb*MAX(1,Descriptor%nsurf)))
    ALLOCATE(fb(SIZE(fa)))
    CALL UnpackSmolyakTermsToFullBasis(fa,a,Descriptor)
    CALL UnpackSmolyakTermsToFullBasis(fb,b,Descriptor)
    value=DOT_PRODUCT(fa,fb)
    DEALLOCATE(fa,fb)
  END FUNCTION DotSmolyakBasis

  !=========================================================================
  ! The expand sparse basis to full basis function
  !
  ! This function expands a sparse basis to a full basis.
  !
  ! Parameters:
  !   full_vector : Full basis vector
  !   sparse_vector : Sparse basis vector
  !   Descriptor : Smolyak sparse grid descriptor
  !
  ! Returns:
  !   None (modifies full_vector in place)
  !
  SUBROUTINE ExpandSparseBasisToFullBasis(full_vector, sparse_vector, Descriptor)
    COMPLEX(KIND=Rkind), INTENT(OUT) :: full_vector(:)
    COMPLEX(KIND=Rkind), INTENT(IN)  :: sparse_vector(:)
    TYPE(SmolyakSparseGridDescriptor_t), INTENT(IN) :: Descriptor
    INTEGER :: isurf, offset_dense, offset_sparse
    IF (Descriptor%layout_kind /= BASIS_LAYOUT_SMOLYAK_SPARSE .OR. .NOT. Descriptor%is_initialized) THEN
       full_vector = sparse_vector
    ELSEIF (SIZE(sparse_vector) == GetSmolyakPackedBasisSize(Descriptor)) THEN
       CALL UnpackSmolyakTermsToFullBasis(full_vector, sparse_vector, Descriptor)
    ELSE
       full_vector = CZERO
       DO isurf=1,Descriptor%nsurf
          offset_dense=(isurf-1)*Descriptor%nb
          offset_sparse=(isurf-1)*Descriptor%nb_SG
          full_vector(offset_dense+Descriptor%active_nuclear_basis_index(:)) = &
             sparse_vector(offset_sparse+1:offset_sparse+Descriptor%nb_SG)
       END DO
    END IF
  END SUBROUTINE ExpandSparseBasisToFullBasis

  !=========================================================================
  ! The restrict full basis to sparse basis function
  !
  ! This function restricts a full basis to a sparse basis.
  !
  ! Parameters:
  !   sparse_vector : Sparse basis vector
  !   full_vector : Full basis vector
  !   Descriptor : Smolyak sparse grid descriptor
  !
  ! Returns:
  !   None (modifies sparse_vector in place)
  !
  SUBROUTINE RestrictFullBasisToSparseBasis(sparse_vector, full_vector, Descriptor)
    COMPLEX(KIND=Rkind), INTENT(OUT) :: sparse_vector(:)
    COMPLEX(KIND=Rkind), INTENT(IN)  :: full_vector(:)
    TYPE(SmolyakSparseGridDescriptor_t), INTENT(IN) :: Descriptor
    INTEGER :: isurf, offset_dense, offset_sparse
    IF (Descriptor%layout_kind /= BASIS_LAYOUT_SMOLYAK_SPARSE .OR. .NOT. Descriptor%is_initialized) THEN
       sparse_vector = full_vector
    ELSEIF (SIZE(sparse_vector) == GetSmolyakPackedBasisSize(Descriptor)) THEN
       CALL PackFullBasisToSmolyakTerms(sparse_vector, full_vector, Descriptor)
    ELSE
       DO isurf = 1, Descriptor%nsurf
          offset_dense  = (isurf - 1) * Descriptor%nb
          offset_sparse = (isurf - 1) * Descriptor%nb_SG
          sparse_vector(offset_sparse + 1 : offset_sparse + Descriptor%nb_SG) = &
             full_vector(offset_dense + Descriptor%active_nuclear_basis_index(:))
       END DO
    END IF
  END SUBROUTINE RestrictFullBasisToSparseBasis

  !=========================================================================
  ! The print smolyak sparse grid descriptor function
  !
  ! This function prints the smolyak sparse grid descriptor.
  !
  ! Parameters:
  !   Descriptor : Smolyak sparse grid descriptor
  !   unit : Output unit
  !
  ! Returns:
  !   None
  !
  SUBROUTINE PrintSmolyakSparseGridDescriptor(Descriptor, unit)
    TYPE(SmolyakSparseGridDescriptor_t), INTENT(IN) :: Descriptor
    INTEGER, INTENT(IN), OPTIONAL :: unit
    INTEGER :: u, it
    u = out_unit
    IF (PRESENT(unit)) u = unit
    WRITE(u,*) 'Smolyak sparse-grid descriptor:'
    WRITE(u,*) '  ndim                     = ', Descriptor%ndim
    WRITE(u,*) '  nsurf                    = ', Descriptor%nsurf
    WRITE(u,*) '  order_LB                 = ', Descriptor%order_LB
    WRITE(u,*) '  order_LG                 = ', Descriptor%order_LG
    WRITE(u,*) '  qmin_B/qmax_B            = ', Descriptor%qmin_B, Descriptor%qmax_B
    WRITE(u,*) '  qmin_G/qmax_G            = ', Descriptor%qmin_G, Descriptor%qmax_G
    WRITE(u,*) '  rule_shift_B             = ', Descriptor%rule_shift_B
    WRITE(u,*) '  rule_slope_B             = ', Descriptor%rule_slope_B
    WRITE(u,*) '  rule_shift_G             = ', Descriptor%rule_shift_G
    WRITE(u,*) '  rule_slope_G             = ', Descriptor%rule_slope_G
    WRITE(u,*) '  nterms_B/nterms_G        = ', Descriptor%nterms_b, Descriptor%nterms_g
    WRITE(u,*) '  nb = ', Descriptor%nb
    WRITE(u,*) '  nb_SG        = ', Descriptor%nb_SG
    WRITE(u,*) '  nq = ', Descriptor%nq
    WRITE(u,*) '  nq_SG        = ', Descriptor%nq_SG
    WRITE(u,*) '  --- termes de base B (LB) ---'
    DO it = 1, Descriptor%nterms_b
       WRITE(u,'(A,I4,A,I4,A,I4,A,100(I4))') '  B term ', Descriptor%term_B(it)%term_id, ' coeff=', Descriptor%term_B(it)%coefficient, ' q_B=', Descriptor%term_B(it)%q, ' ell_B=', Descriptor%term_B(it)%ell
    END DO
    WRITE(u,*) '  --- termes de grille G (LG) ---'
    DO it = 1, Descriptor%nterms_g
       WRITE(u,'(A,I4,A,I4,A,I4,A,100(I4))') '  G term ', Descriptor%term_G(it)%term_id, ' coeff=', Descriptor%term_G(it)%coefficient, ' q_G=', Descriptor%term_G(it)%q, ' ell_G=', Descriptor%term_G(it)%ell
    END DO
  END SUBROUTINE PrintSmolyakSparseGridDescriptor

END MODULE SmolyakSparseGridMetadata_m
