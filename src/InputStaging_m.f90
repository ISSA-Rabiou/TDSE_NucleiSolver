MODULE InputStaging_m
   USE QDUtil_m, ONLY: in_unit
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: Stage_Standard_Input
CONTAINS

   SUBROUTINE Stage_Standard_Input()

      INTEGER :: cache_unit, ios

      CHARACTER(LEN=8192) :: line



      OPEN(NEWUNIT=cache_unit, FILE='TD_SCHROD.stdin.cache', STATUS='REPLACE', &

           ACTION='WRITE', IOSTAT=ios)

      IF (ios /= 0) ERROR STOP 'TD_SCHROD: cannot create stdin cache'

      DO

         READ(in_unit,'(A)',IOSTAT=ios) line

         IF (ios < 0) EXIT

         IF (ios > 0) ERROR STOP 'TD_SCHROD: cannot stage stdin'

         WRITE(cache_unit,'(A)') TRIM(line)

      END DO

      CLOSE(cache_unit)

      CLOSE(in_unit)

      OPEN(UNIT=in_unit, FILE='TD_SCHROD.stdin.cache', STATUS='OLD', ACTION='READ', &

           IOSTAT=ios)

      IF (ios /= 0) ERROR STOP 'TD_SCHROD: cannot reopen stdin cache'

   END SUBROUTINE Stage_Standard_Input

END MODULE InputStaging_m
