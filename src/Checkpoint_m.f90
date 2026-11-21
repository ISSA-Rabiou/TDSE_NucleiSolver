MODULE Checkpoint_m
   USE QDUtil_m, ONLY: Rkind, Ikind, ILkind, ZERO, ONE, ONETENTH, in_unit
   USE psi_m, ONLY: psi_t
   USE Basis_m, ONLY: Basis_t
   USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_char, c_int, c_null_char
   IMPLICIT NONE
   PRIVATE

   CHARACTER(LEN=*), PARAMETER :: checkpoint_magic = 'TD_SCHROD_CHECKPOINT_V1'
   INTEGER(KIND=ILkind), PARAMETER :: signature_modulus = INT(HUGE(0_Ikind), KIND=ILkind)

   TYPE, PUBLIC :: checkpoint_config_t
      REAL(KIND=Rkind) :: t0 = ZERO
      REAL(KIND=Rkind) :: delta_t = ZERO
      REAL(KIND=Rkind) :: eps = ZERO
      REAL(KIND=Rkind) :: E_min = ZERO
      REAL(KIND=Rkind) :: E_max = ZERO
      INTEGER(KIND=Ikind) :: max_iter = 0_Ikind
      INTEGER(KIND=Ikind) :: split_order = 2_Ikind
      LOGICAL :: Beta = .FALSE.
      LOGICAL :: P = .FALSE.
      LOGICAL :: renorm = .FALSE.
      CHARACTER(LEN=32) :: propa_name = ''
      CHARACTER(LEN=32) :: integrator = ''
      INTEGER(KIND=ILkind) :: input_signature = 0_ILkind
   END TYPE checkpoint_config_t

   PUBLIC :: Input_Signature
   PUBLIC :: Save_Propagation_Checkpoint, Load_Propagation_Checkpoint
   PUBLIC :: Trim_Propagation_Outputs

   INTERFACE
      FUNCTION c_rename(oldpath, newpath) BIND(C, NAME='rename') RESULT(rc)
         IMPORT :: c_char, c_int
         CHARACTER(KIND=c_char), INTENT(IN) :: oldpath(*)
         CHARACTER(KIND=c_char), INTENT(IN) :: newpath(*)
         INTEGER(c_int) :: rc
      END FUNCTION c_rename
   END INTERFACE

CONTAINS

   FUNCTION Input_Signature() RESULT(signature)
      INTEGER(KIND=ILkind) :: signature
      INTEGER(KIND=ILkind) :: hash
      INTEGER(KIND=Ikind) :: ios, i, nchar
      CHARACTER(LEN=8192) :: line
      CHARACTER(LEN=8192) :: lowered
      INTEGER(KIND=Ikind) :: code

      hash = 146959810_ILkind
      REWIND(in_unit, IOSTAT=ios)
      IF (ios /= 0_Ikind) ERROR STOP 'Input_Signature: the staged input is not rewindable.'

      DO
         READ(in_unit, '(A)', IOSTAT=ios) line
         IF (ios /= 0_Ikind) EXIT
         lowered = line
         DO i = 1_Ikind, LEN_TRIM(lowered)
            code = IACHAR(lowered(i:i))
            IF (code >= IACHAR('A') .AND. code <= IACHAR('Z')) lowered(i:i) = ACHAR(code + 32)
         END DO
         IF (INDEX(ADJUSTL(lowered), '&prop') == 1_Ikind) EXIT

         nchar = LEN_TRIM(line)
         DO i = 1_Ikind, nchar
            hash = MODULO(hash * 257_ILkind + INT(IACHAR(line(i:i)), KIND=ILkind), signature_modulus)
         END DO
         hash = MODULO(hash * 257_ILkind + INT(IACHAR(NEW_LINE('A')), KIND=ILkind), signature_modulus)
      END DO

      REWIND(in_unit, IOSTAT=ios)
      IF (ios /= 0_Ikind) ERROR STOP 'Input_Signature: unable to rewind the staged input.'
      signature = hash
   END FUNCTION Input_Signature


   SUBROUTINE Save_Propagation_Checkpoint(file_name, current_time, step, psi, psi_reference, &
                                          stats, qt_reference, pt_reference, config, ierr, message)
      CHARACTER(*), INTENT(IN) :: file_name
      REAL(KIND=Rkind), INTENT(IN) :: current_time
      INTEGER(KIND=Ikind), INTENT(IN) :: step
      TYPE(psi_t), INTENT(IN) :: psi, psi_reference
      REAL(KIND=Rkind), INTENT(IN) :: stats(6), qt_reference(:), pt_reference(:)
      TYPE(checkpoint_config_t), INTENT(IN) :: config
      INTEGER(KIND=Ikind), INTENT(OUT) :: ierr
      CHARACTER(*), INTENT(OUT) :: message

      CHARACTER(LEN=1024) :: temporary_file
      INTEGER(KIND=Ikind) :: unit, ios
      INTEGER(c_int) :: rename_status

      ierr = 0_Ikind
      message = ''
      IF (.NOT. ASSOCIATED(psi%Basis)) THEN
         ierr = 1_Ikind
         message = 'Cannot checkpoint a wavefunction without an associated basis.'
         RETURN
      END IF
      IF (.NOT. ALLOCATED(psi%CVec) .OR. .NOT. ALLOCATED(psi_reference%CVec)) THEN
         ierr = 1_Ikind
         message = 'Cannot checkpoint unallocated wavefunction coefficients.'
         RETURN
      END IF
      IF (SIZE(psi%CVec, KIND=Ikind) /= SIZE(psi_reference%CVec, KIND=Ikind)) THEN
         ierr = 1_Ikind
         message = 'Current and reference coefficient vectors have different sizes.'
         RETURN
      END IF
      IF (LEN_TRIM(file_name) == 0_Ikind .OR. LEN_TRIM(file_name) > 900_Ikind) THEN
         ierr = 1_Ikind
         message = 'Invalid checkpoint_file path (empty or too long).'
         RETURN
      END IF

      temporary_file = TRIM(file_name) // '.tmp'
      OPEN(NEWUNIT=unit, FILE=TRIM(temporary_file), STATUS='REPLACE', ACTION='WRITE', &
           FORM='FORMATTED', IOSTAT=ios)
      IF (ios /= 0_Ikind) THEN
         ierr = ios
         message = 'Unable to create temporary checkpoint file: ' // TRIM(temporary_file)
         RETURN
      END IF

      WRITE(unit, '(A)') checkpoint_magic
      WRITE(unit, *) config%input_signature
      WRITE(unit, '(A)') TRIM(config%propa_name)
      WRITE(unit, '(A)') TRIM(config%integrator)
      WRITE(unit, *) config%t0, config%delta_t, config%eps, config%E_min, config%E_max
      WRITE(unit, *) config%max_iter, config%split_order, config%Beta, config%P, config%renorm
      WRITE(unit, *) current_time, step, SIZE(psi%CVec, KIND=Ikind), SIZE(psi_reference%CVec, KIND=Ikind)
      WRITE(unit, *) stats
      WRITE(unit, *) SIZE(qt_reference, KIND=Ikind), SIZE(pt_reference, KIND=Ikind)
      WRITE(unit, *) qt_reference
      WRITE(unit, *) pt_reference
      CALL Write_Basis_Tree(unit, psi%Basis)
      CALL Write_Coefficients(unit, psi%CVec)
      CALL Write_Coefficients(unit, psi_reference%CVec)
      FLUSH(unit, IOSTAT=ios)
      IF (ios == 0_Ikind) CLOSE(unit, IOSTAT=ios)
      IF (ios /= 0_Ikind) THEN
         ierr = ios
         message = 'Write/close failure while saving checkpoint: ' // TRIM(temporary_file)
         RETURN
      END IF

      rename_status = c_rename(TRIM(temporary_file) // C_NULL_CHAR, TRIM(file_name) // C_NULL_CHAR)
      IF (rename_status /= 0_c_int) THEN
         ierr = INT(rename_status, KIND=Ikind)
         message = 'Atomic replacement of checkpoint failed: ' // TRIM(file_name)
      END IF
   END SUBROUTINE Save_Propagation_Checkpoint


   SUBROUTINE Load_Propagation_Checkpoint(file_name, psi, psi_reference, config, &
                                          current_time, step, stats, qt_reference, pt_reference, ierr, message)
      CHARACTER(*), INTENT(IN) :: file_name
      TYPE(psi_t), INTENT(INOUT) :: psi, psi_reference
      TYPE(checkpoint_config_t), INTENT(IN) :: config
      REAL(KIND=Rkind), INTENT(OUT) :: current_time
      INTEGER(KIND=Ikind), INTENT(OUT) :: step
      REAL(KIND=Rkind), INTENT(OUT) :: stats(6), qt_reference(:), pt_reference(:)
      INTEGER(KIND=Ikind), INTENT(OUT) :: ierr
      CHARACTER(*), INTENT(OUT) :: message

      CHARACTER(LEN=128) :: magic, saved_method, saved_integrator
      INTEGER(KIND=Ikind) :: unit, ios, ncoef, nref_coeff, nref_q, nref_p, saved_max_iter, saved_split_order
      INTEGER(KIND=Ikind) :: saved_step
      INTEGER(KIND=ILkind) :: saved_signature
      REAL(KIND=Rkind) :: saved_t0, saved_dt, saved_eps, saved_emin, saved_emax, saved_time
      LOGICAL :: saved_beta, saved_p, saved_renorm
      REAL(KIND=Rkind) :: tol

      ierr = 0_Ikind
      message = ''
      current_time = ZERO
      step = 0_Ikind
      stats = ZERO
      IF (LEN_TRIM(file_name) == 0_Ikind) THEN
         ierr = 1_Ikind
         message = 'restart=T but checkpoint_file is empty.'
         RETURN
      END IF
      IF (.NOT. ASSOCIATED(psi%Basis) .OR. .NOT. ALLOCATED(psi%CVec) .OR. &
          .NOT. ALLOCATED(psi_reference%CVec)) THEN
         ierr = 1_Ikind
         message = 'Restart destination wavefunction/basis is not initialized.'
         RETURN
      END IF

      OPEN(NEWUNIT=unit, FILE=TRIM(file_name), STATUS='OLD', ACTION='READ', &
           FORM='FORMATTED', IOSTAT=ios)
      IF (ios /= 0_Ikind) THEN
         ierr = ios
         message = 'restart=T but the checkpoint file cannot be opened: ' // TRIM(file_name)
         RETURN
      END IF

      READ(unit, '(A)', IOSTAT=ios) magic
      IF (ios == 0_Ikind) READ(unit, *, IOSTAT=ios) saved_signature
      IF (ios == 0_Ikind) READ(unit, '(A)', IOSTAT=ios) saved_method
      IF (ios == 0_Ikind) READ(unit, '(A)', IOSTAT=ios) saved_integrator
      IF (ios == 0_Ikind) READ(unit, *, IOSTAT=ios) saved_t0, saved_dt, saved_eps, saved_emin, saved_emax
      IF (ios == 0_Ikind) READ(unit, *, IOSTAT=ios) saved_max_iter, saved_split_order, &
           saved_beta, saved_p, saved_renorm
      IF (ios == 0_Ikind) READ(unit, *, IOSTAT=ios) saved_time, saved_step, ncoef, nref_coeff
      IF (ios == 0_Ikind) READ(unit, *, IOSTAT=ios) stats
      IF (ios == 0_Ikind) READ(unit, *, IOSTAT=ios) nref_q, nref_p
      IF (ios /= 0_Ikind) THEN
         CLOSE(unit)
         ierr = ios
         message = 'Checkpoint is truncated or has an invalid header.'
         RETURN
      END IF

      IF (TRIM(magic) /= checkpoint_magic) THEN
         CLOSE(unit)
         ierr = 1_Ikind
         message = 'Unsupported or corrupt checkpoint version.'
         RETURN
      END IF
      IF (saved_signature /= config%input_signature) THEN
         CLOSE(unit)
         ierr = 1_Ikind
         message = 'Input signature mismatch: keep the same pre-&prop input (potential, coordinates, basis, initial packet).'
         RETURN
      END IF
      IF (TRIM(saved_method) /= TRIM(config%propa_name) .OR. &
          TRIM(saved_integrator) /= TRIM(config%integrator)) THEN
         CLOSE(unit)
         ierr = 1_Ikind
         message = 'Propagation method/integrator differs from the checkpoint.'
         RETURN
      END IF
      tol = ONETENTH**12 * MAX(ONE, ABS(config%t0), ABS(config%delta_t))
      IF (ABS(saved_t0-config%t0) > tol .OR. ABS(saved_dt-config%delta_t) > tol .OR. &
          ABS(saved_eps-config%eps) > MAX(TINY(ONE), ONETENTH**12*ABS(config%eps)) .OR. &
          ABS(saved_emin-config%E_min) > ONETENTH**12*MAX(ONE,ABS(config%E_min)) .OR. &
          ABS(saved_emax-config%E_max) > ONETENTH**12*MAX(ONE,ABS(config%E_max))) THEN
         CLOSE(unit)
         ierr = 1_Ikind
         message = 't0, delta_t or eps differs from the checkpoint.'
         RETURN
      END IF
      IF (saved_max_iter /= config%max_iter .OR. saved_split_order /= config%split_order .OR. &
          (saved_beta .NEQV. config%Beta) .OR. (saved_p .NEQV. config%P) .OR. &
          (saved_renorm .NEQV. config%renorm)) THEN
         CLOSE(unit)
         ierr = 1_Ikind
         message = 'Propagation accuracy/physics controls differ from the checkpoint.'
         RETURN
      END IF
      IF (ncoef /= SIZE(psi%CVec, KIND=Ikind) .OR. nref_coeff /= SIZE(psi_reference%CVec, KIND=Ikind) .OR. &
          nref_q /= SIZE(qt_reference, KIND=Ikind) .OR. nref_p /= SIZE(pt_reference, KIND=Ikind) .OR. &
          SIZE(psi_reference%CVec, KIND=Ikind) /= SIZE(psi%CVec, KIND=Ikind)) THEN
         CLOSE(unit)
         ierr = 1_Ikind
         message = 'Checkpoint coefficient or reference-array dimensions do not match this input.'
         RETURN
      END IF
      IF (nref_q /= nref_p) THEN
         CLOSE(unit)
         ierr = 1_Ikind
         message = 'Checkpoint contains inconsistent basis-reference dimensions.'
         RETURN
      END IF

      READ(unit, *, IOSTAT=ios) qt_reference
      IF (ios == 0_Ikind) READ(unit, *, IOSTAT=ios) pt_reference
      IF (ios /= 0_Ikind) THEN
         CLOSE(unit)
         ierr = ios
         message = 'Checkpoint is truncated while reading reference observables.'
         RETURN
      END IF

      current_time = saved_time
      step = saved_step
      CALL Read_Basis_Tree(unit, psi%Basis, ierr, message)
      IF (ierr == 0_Ikind) CALL Read_Coefficients(unit, psi%CVec, ierr, message)
      IF (ierr == 0_Ikind) CALL Read_Coefficients(unit, psi_reference%CVec, ierr, message)
      CLOSE(unit)
      IF (ierr /= 0_Ikind) RETURN

      psi%Grid = .FALSE.
      psi_reference%Grid = .FALSE.
   END SUBROUTINE Load_Propagation_Checkpoint


   RECURSIVE SUBROUTINE Write_Basis_Tree(unit, basis)
      INTEGER(KIND=Ikind), INTENT(IN) :: unit
      TYPE(Basis_t), INTENT(IN) :: basis
      INTEGER(KIND=Ikind) :: nchild, i
      CHARACTER(LEN=128) :: basis_name

      nchild = 0_Ikind
      IF (ALLOCATED(basis%tab_basis)) nchild = SIZE(basis%tab_basis, KIND=Ikind)
      basis_name = ''
      IF (ALLOCATED(basis%Basis_name)) basis_name = TRIM(basis%Basis_name)
      WRITE(unit, *) basis%nb_basis, basis%nb, basis%nq, basis%representation_layout, &
           basis%LB, basis%LG, nchild
      WRITE(unit, '(A)') TRIM(basis_name)
      WRITE(unit, *) basis%Q0, basis%scaleQ, basis%imp_k, REAL(basis%Alpha, KIND=Rkind), AIMAG(basis%Alpha)
      IF (ALLOCATED(basis%tab_basis)) THEN
         DO i = 1_Ikind, SIZE(basis%tab_basis, KIND=Ikind)
            CALL Write_Basis_Tree(unit, basis%tab_basis(i))
         END DO
      END IF
   END SUBROUTINE Write_Basis_Tree


   RECURSIVE SUBROUTINE Read_Basis_Tree(unit, basis, ierr, message)
      INTEGER(KIND=Ikind), INTENT(IN) :: unit
      TYPE(Basis_t), INTENT(INOUT) :: basis
      INTEGER(KIND=Ikind), INTENT(OUT) :: ierr
      CHARACTER(*), INTENT(OUT) :: message
      INTEGER(KIND=Ikind) :: ios, nb_basis, nb, nq, layout, lb, lg, nchild, expected_children, i
      REAL(KIND=Rkind) :: q0, scaleq, impk, alpha_real, alpha_imag
      CHARACTER(LEN=128) :: basis_name, current_name

      ierr = 0_Ikind
      message = ''
      READ(unit, *, IOSTAT=ios) nb_basis, nb, nq, layout, lb, lg, nchild
      IF (ios == 0_Ikind) READ(unit, '(A)', IOSTAT=ios) basis_name
      IF (ios == 0_Ikind) READ(unit, *, IOSTAT=ios) q0, scaleq, impk, alpha_real, alpha_imag
      IF (ios /= 0_Ikind) THEN
         ierr = ios
         message = 'Truncated checkpoint while reading basis metadata.'
         RETURN
      END IF

      expected_children = 0_Ikind
      IF (ALLOCATED(basis%tab_basis)) expected_children = SIZE(basis%tab_basis, KIND=Ikind)
      current_name = ''
      IF (ALLOCATED(basis%Basis_name)) current_name = TRIM(basis%Basis_name)
      IF (nb_basis /= basis%nb_basis .OR. nb /= basis%nb .OR. nq /= basis%nq .OR. &
          layout /= basis%representation_layout .OR. lb /= basis%LB .OR. lg /= basis%LG .OR. &
          nchild /= expected_children .OR. TRIM(basis_name) /= TRIM(current_name)) THEN
         ierr = 1_Ikind
         message = 'Basis layout/name/size mismatch while loading checkpoint.'
         RETURN
      END IF

      basis%Q0 = q0
      basis%scaleQ = scaleq
      basis%imp_k = impk
      basis%Alpha = CMPLX(alpha_real, alpha_imag, KIND=Rkind)
      IF (ALLOCATED(basis%tab_basis)) THEN
         DO i = 1_Ikind, SIZE(basis%tab_basis, KIND=Ikind)
            CALL Read_Basis_Tree(unit, basis%tab_basis(i), ierr, message)
            IF (ierr /= 0_Ikind) RETURN
         END DO
      END IF
   END SUBROUTINE Read_Basis_Tree


   SUBROUTINE Write_Coefficients(unit, coefficients)
      INTEGER(KIND=Ikind), INTENT(IN) :: unit
      COMPLEX(KIND=Rkind), INTENT(IN) :: coefficients(:)
      INTEGER(KIND=Ikind) :: i
      DO i = 1_Ikind, SIZE(coefficients, KIND=Ikind)
         WRITE(unit, '(2(1X,ES26.17E3))') REAL(coefficients(i), KIND=Rkind), AIMAG(coefficients(i))
      END DO
   END SUBROUTINE Write_Coefficients


   SUBROUTINE Read_Coefficients(unit, coefficients, ierr, message)
      INTEGER(KIND=Ikind), INTENT(IN) :: unit
      COMPLEX(KIND=Rkind), INTENT(OUT) :: coefficients(:)
      INTEGER(KIND=Ikind), INTENT(OUT) :: ierr
      CHARACTER(*), INTENT(OUT) :: message
      INTEGER(KIND=Ikind) :: i, ios
      REAL(KIND=Rkind) :: re, im

      ierr = 0_Ikind
      message = ''
      DO i = 1_Ikind, SIZE(coefficients, KIND=Ikind)
         READ(unit, *, IOSTAT=ios) re, im
         IF (ios /= 0_Ikind) THEN
            ierr = ios
            message = 'Checkpoint is truncated while reading wavefunction coefficients.'
            RETURN
         END IF
         coefficients(i) = CMPLX(re, im, KIND=Rkind)
      END DO
   END SUBROUTINE Read_Coefficients


   SUBROUTINE Trim_Propagation_Outputs(checkpoint_time, ierr, message)
      REAL(KIND=Rkind), INTENT(IN) :: checkpoint_time
      INTEGER(KIND=Ikind), INTENT(OUT) :: ierr
      CHARACTER(*), INTENT(OUT) :: message
      CHARACTER(LEN=32), PARAMETER :: output_names(11_Ikind) = [CHARACTER(LEN=32) :: &
           'Qt', 'E', 'SQt', 'Norm', 'pop', 'Imp_k', 'alpha', 'auto_cor', 'E_13', 'Norm_13', 'validation']
      INTEGER(KIND=Ikind) :: i

      ierr = 0_Ikind
      message = ''
      DO i = 1_Ikind, SIZE(output_names, KIND=Ikind)
         CALL Trim_One_Output(TRIM(output_names(i)) // '.txt', checkpoint_time, ierr, message)
         IF (ierr /= 0_Ikind) RETURN
      END DO
   END SUBROUTINE Trim_Propagation_Outputs


   SUBROUTINE Trim_One_Output(file_name, checkpoint_time, ierr, message)
      CHARACTER(*), INTENT(IN) :: file_name
      REAL(KIND=Rkind), INTENT(IN) :: checkpoint_time
      INTEGER(KIND=Ikind), INTENT(OUT) :: ierr
      CHARACTER(*), INTENT(OUT) :: message
      CHARACTER(LEN=1024) :: temp_file
      CHARACTER(LEN=8192) :: line, adjusted
      INTEGER(KIND=Ikind) :: in_file, out_file, ios, write_status, rename_status
      REAL(KIND=Rkind) :: row_time, tol
      LOGICAL :: exists, in_summary

      ierr = 0_Ikind
      message = ''
      INQUIRE(FILE=TRIM(file_name), EXIST=exists)
      IF (.NOT. exists) RETURN
      temp_file = TRIM(file_name) // '.restart.tmp'
      OPEN(NEWUNIT=in_file, FILE=TRIM(file_name), STATUS='OLD', ACTION='READ', &
           FORM='FORMATTED', IOSTAT=ios)
      IF (ios /= 0_Ikind) THEN
         ierr = ios
         message = 'Unable to read observable file during restart: ' // TRIM(file_name)
         RETURN
      END IF
      OPEN(NEWUNIT=out_file, FILE=TRIM(temp_file), STATUS='REPLACE', ACTION='WRITE', &
           FORM='FORMATTED', IOSTAT=ios)
      IF (ios /= 0_Ikind) THEN
         CLOSE(in_file)
         ierr = ios
         message = 'Unable to create temporary observable file: ' // TRIM(temp_file)
         RETURN
      END IF

      tol = ONETENTH**10 * MAX(ONE, ABS(checkpoint_time))
      in_summary = .FALSE.
      DO
         READ(in_file, '(A)', IOSTAT=ios) line
         IF (ios /= 0_Ikind) EXIT
         adjusted = ADJUSTL(line)
         IF (INDEX(adjusted, '# summary') == 1_Ikind) THEN
            in_summary = .TRUE.
            CYCLE
         END IF
         IF (in_summary) CYCLE
         IF (LEN_TRIM(adjusted) == 0_Ikind) THEN
            WRITE(out_file, '(A)') TRIM(line)
            CYCLE
         END IF
         IF (adjusted(1:1) == '#') THEN
            WRITE(out_file, '(A)') TRIM(line)
            CYCLE
         END IF
         READ(line, *, IOSTAT=ios) row_time
         IF (ios == 0_Ikind) THEN
            IF (row_time <= checkpoint_time + tol) WRITE(out_file, '(A)') TRIM(line)
         ELSE
            WRITE(out_file, '(A)') TRIM(line)
         END IF
      END DO
      CLOSE(in_file)
      CLOSE(out_file, IOSTAT=write_status)
      IF (write_status /= 0_Ikind) THEN
         ierr = write_status
         message = 'Unable to close trimmed observable file: ' // TRIM(temp_file)
         RETURN
      END IF

      rename_status = INT(c_rename(TRIM(temp_file) // C_NULL_CHAR, TRIM(file_name) // C_NULL_CHAR), KIND=Ikind)
      IF (rename_status /= 0_Ikind) THEN
         ierr = rename_status
         message = 'Atomic replacement of trimmed observable file failed: ' // TRIM(file_name)
      END IF
   END SUBROUTINE Trim_One_Output

END MODULE Checkpoint_m
