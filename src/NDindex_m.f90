MODULE NDindex_m
  USE QDUtil_m, ONLY: out_unit
  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: NDINDEX_TENSOR  = 0
  INTEGER, PARAMETER, PUBLIC :: NDINDEX_SMOLYAK = 1

  TYPE, PUBLIC :: NDindex_t
    INTEGER                       :: Ndim = 0
    INTEGER                       :: mode = NDINDEX_TENSOR
    INTEGER                       :: qmin = 0
    INTEGER                       :: qmax = -1
    INTEGER, ALLOCATABLE          :: Tab0(:)
    INTEGER, ALLOCATABLE          :: NDsize(:)
    INTEGER, ALLOCATABLE          :: NDend(:)
  END TYPE NDindex_t

  PUBLIC :: increase_NDindex, Init_tab_ind, Testindex, Init_NDindex
  PUBLIC :: dealloc_NDindex, Write_NDindex

CONTAINS

  SUBROUTINE Init_NDindex(NDindex, NDend, Ndim, index_mode, smolyak_qmin, smolyak_qmax)
    TYPE(NDindex_t), INTENT(INOUT) :: NDindex
    INTEGER, INTENT(IN)            :: NDend(:)
    INTEGER, INTENT(IN)            :: Ndim
    INTEGER, INTENT(IN), OPTIONAL   :: index_mode, smolyak_qmin, smolyak_qmax
    INTEGER :: i

    CALL dealloc_NDindex(NDindex)
    IF (Ndim < 1 .OR. SIZE(NDend) /= Ndim) STOP 'Init_NDindex: invalid dimension or upper bounds'

    NDindex%Ndim = Ndim
    NDindex%NDend = NDend
    ALLOCATE(NDindex%Tab0(Ndim), NDindex%NDsize(Ndim))
    NDindex%NDsize = NDend + 1
    NDindex%mode = NDINDEX_TENSOR
    IF (PRESENT(index_mode)) NDindex%mode = index_mode

    IF (NDindex%mode == NDINDEX_SMOLYAK) THEN
      NDindex%qmin = 0
      NDindex%qmax = SUM(NDend)
      IF (PRESENT(smolyak_qmin)) NDindex%qmin = MAX(0, smolyak_qmin)
      IF (PRESENT(smolyak_qmax)) NDindex%qmax = MIN(SUM(NDend), smolyak_qmax)
      IF (NDindex%qmax < NDindex%qmin) STOP 'Init_NDindex: invalid Smolyak q interval'
      NDindex%Tab0 = 0
    ELSE
      ! Historical convention: the first index starts at zero and all
      ! following components start at one; the first increment yields 1-based
      ! tensor-product indices in the order used by the project.
      NDindex%Tab0 = 1
      NDindex%Tab0(1) = 0
    END IF
    DO i=1,Ndim
      IF (NDindex%Tab0(i) > NDindex%NDend(i)) STOP 'Init_NDindex: empty index range'
    END DO
  END SUBROUTINE Init_NDindex

  SUBROUTINE Init_tab_ind(Tab_ind, NDindex)
    TYPE(NDindex_t), INTENT(IN)    :: NDindex
    INTEGER, INTENT(INOUT)         :: Tab_ind(:)
    IF (SIZE(Tab_ind) /= NDindex%Ndim) STOP 'Init_tab_ind: incompatible size'
    Tab_ind = NDindex%Tab0
  END SUBROUTINE Init_tab_ind

  SUBROUTINE increase_NDindex(Tab_ind, NDindex, Endloop)
    TYPE(NDindex_t), INTENT(IN) :: NDindex
    INTEGER, INTENT(INOUT)      :: Tab_ind(:)
    LOGICAL, INTENT(INOUT)      :: Endloop
    INTEGER :: i

    IF (NDindex%mode == NDINDEX_SMOLYAK) THEN
      DO
        Tab_ind(1) = Tab_ind(1) + 1
        DO i=1,NDindex%Ndim-1
          IF (Tab_ind(i) <= NDindex%NDend(i)) EXIT
          Tab_ind(i) = 0
          Tab_ind(i+1) = Tab_ind(i+1) + 1
        END DO
        IF (Tab_ind(NDindex%Ndim) > NDindex%NDend(NDindex%Ndim)) THEN
          Endloop = .TRUE.
          RETURN
        END IF
        IF (SUM(Tab_ind) >= NDindex%qmin .AND. SUM(Tab_ind) <= NDindex%qmax) THEN
          Endloop = .FALSE.
          RETURN
        END IF
      END DO
    END IF

    Tab_ind(1) = Tab_ind(1) + 1
    DO i=1,NDindex%Ndim-1
      IF (Tab_ind(i) > NDindex%NDend(i)) THEN
        Tab_ind(i+1) = Tab_ind(i+1) + 1
        Tab_ind(i) = 1
      END IF
    END DO
    Endloop = (Tab_ind(NDindex%Ndim) == NDindex%NDend(NDindex%Ndim)+1)
  END SUBROUTINE increase_NDindex

  SUBROUTINE Write_NDindex(NDindex)
    TYPE(NDindex_t), INTENT(IN) :: NDindex
    WRITE(out_unit,*) 'NDindex%Ndim =', NDindex%Ndim
    WRITE(out_unit,*) 'NDindex%mode =', NDindex%mode
    WRITE(out_unit,*) 'NDindex%qmin/qmax =', NDindex%qmin, NDindex%qmax
    WRITE(out_unit,*) 'NDindex%NDend =', NDindex%NDend
    WRITE(out_unit,*) 'NDindex%Tab0 =', NDindex%Tab0
  END SUBROUTINE Write_NDindex

  SUBROUTINE dealloc_NDindex(NDindex)
    TYPE(NDindex_t), INTENT(INOUT) :: NDindex
    IF (ALLOCATED(NDindex%NDend))  DEALLOCATE(NDindex%NDend)
    IF (ALLOCATED(NDindex%NDsize)) DEALLOCATE(NDindex%NDsize)
    IF (ALLOCATED(NDindex%Tab0))   DEALLOCATE(NDindex%Tab0)
    NDindex%Ndim = 0
    NDindex%mode = NDINDEX_TENSOR
    NDindex%qmin = 0
    NDindex%qmax = -1
  END SUBROUTINE dealloc_NDindex

  SUBROUTINE Testindex(NDindex)
    TYPE(NDindex_t), INTENT(IN) :: NDindex
    INTEGER, ALLOCATABLE :: Tab_ind(:)
    LOGICAL :: Endloop
    ALLOCATE(Tab_ind(NDindex%Ndim))
    CALL Init_tab_ind(Tab_ind, NDindex)
    DO
      CALL increase_NDindex(Tab_ind, NDindex, Endloop)
      IF (Endloop) EXIT
      WRITE(out_unit,*) Tab_ind
    END DO
    DEALLOCATE(Tab_ind)
  END SUBROUTINE Testindex

END MODULE NDindex_m
