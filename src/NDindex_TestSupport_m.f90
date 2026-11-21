MODULE NDindex_TestSupport_m
  USE QDUtil_m, ONLY: out_unit
  USE NDindex_m, ONLY: NDindex_t, Init_tab_ind, increase_NDindex
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Testindex
CONTAINS

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


END MODULE NDindex_TestSupport_m
