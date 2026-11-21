MODULE Propa_m
   USE QDUtil_m
   USE psi_m
   USE Op_m
   USE Ana_psi_m
   USE lanczos_m
   USE Auto_corr_m
   USE Hagedorn_m
   USE sub_propa_m
   USE Basis_m
   USE Checkpoint_m
   IMPLICIT NONE

CONTAINS

   SUBROUTINE propagation(psif, psi0, propa, H)
      TYPE(psi_t),   INTENT(INOUT) :: psif
      TYPE(psi_t),   INTENT(IN)    :: psi0
      TYPE(propa_t), INTENT(INOUT) :: propa
      TYPE(Op_t),    INTENT(INOUT) :: H

      TYPE(Basis_t), TARGET :: Basis
      TYPE(psi_t) :: psi, psi_dt, psi_t0, Gspsi
      TYPE(checkpoint_config_t) :: checkpoint_config
      REAL(KIND=Rkind), ALLOCATABLE :: Qt(:), SQt(:), Pt(:), Qt_ref(:), Pt_ref(:), pop(:)
      COMPLEX(KIND=Rkind), ALLOCATABLE :: At(:)
      COMPLEX(KIND=Rkind) :: aut_func(1)
      REAL(KIND=Rkind) :: aut_func_arg(1)
      REAL(KIND=Rkind) :: t, current_time, Norm, E, E0, Norm_ref, E_ref
      REAL(KIND=Rkind) :: max_norm_drift, max_energy_drift, max_abs_dQ, max_abs_dP
      REAL(KIND=Rkind) :: step_max_abs_dQ, step_max_abs_dP, auto_mod_dev
      REAL(KIND=Rkind) :: ratio, time_tol
      REAL(KIND=Rkind) :: stats(6)
      INTEGER, ALLOCATABLE :: Tab_iq(:, :)
      INTEGER :: u_Qt, u_E, u_SQt, u_Norm, u_pop, u_Pt, u_alpha
      INTEGER :: u_Norm_13, u_E_13, u_auto_cor, u_validation
      INTEGER(KIND=Ikind) :: ndim, nsurf, total_steps, current_step, i, ierr, file_size
      CHARACTER(LEN=Name_len) :: name
      CHARACTER(LEN=512) :: checkpoint_message
      LOGICAL :: need_basis, need_energy, need_norm, need_autocorrelation
      LOGICAL :: append_outputs, validation_header_needed, file_exists
      LOGICAL :: time_propagation
      REAL(KIND=Rkind), PARAMETER :: validation_tol_norm = ONETENTH**6
      REAL(KIND=Rkind), PARAMETER :: validation_tol_energy = ONETENTH**5

      u_Qt = 0; u_E = 0; u_SQt = 0; u_Norm = 0; u_pop = 0; u_Pt = 0
      u_alpha = 0; u_Norm_13 = 0; u_E_13 = 0; u_auto_cor = 0; u_validation = 0
      name = propa%propa_name2
      CALL string_uppercase_TO_lowercase(name)
      time_propagation = (name /= 'itp' .AND. name /= 'gssil')
      IF (propa%restart .AND. .NOT. time_propagation) &
         ERROR STOP 'propagation: ITP/GSSIL are standalone calculations and cannot be restarted as a time propagation.'

      IF (propa%debug) THEN
         WRITE(out_unit, *) 'BEGINNING propagation', propa%t0, propa%tf, propa%delta_t
         WRITE(out_unit, *) 'restart/checkpoint: ', propa%restart, propa%checkpoint, &
              TRIM(propa%checkpoint_file), ' every ', propa%checkpoint_every, ' step(s)'
         WRITE(out_unit, *) '------------- propagation parameters ---------------'
         CALL write_propa(propa)
      END IF

      CALL init_Basis1_TO_Basis2(Basis, psi0%Basis)
      CALL construct_primitive_basis(Basis)
      ndim = SIZE(psi0%Basis%tab_basis) - 1
      nsurf = psi0%Basis%tab_basis(ndim+1)%nb

      ALLOCATE(Qt(ndim), SQt(ndim), Pt(ndim), Qt_ref(ndim), Pt_ref(ndim), At(ndim))
      ALLOCATE(pop(nsurf))
      Qt = ZERO; SQt = ZERO; Pt = ZERO; Qt_ref = ZERO; Pt_ref = ZERO; At = CZERO
      pop = ZERO; aut_func = CZERO; aut_func_arg = ZERO

      CALL init_psi(psi, psi0%Basis, cplx=.TRUE., grid=.FALSE.)
      CALL init_psi(psi_dt, psi0%Basis, cplx=.TRUE., grid=.FALSE.)
      CALL init_psi(psi_t0, Basis, cplx=.TRUE., grid=.FALSE.)
      CALL init_psi(Gspsi, Basis, cplx=.TRUE., grid=.FALSE.)
      psi%CVec(:) = psi0%CVec(:)
      psi_dt%CVec(:) = psi0%CVec(:)
      psi_t0%CVec(:) = psi0%CVec(:)

      need_basis = propa%write_Qt .OR. propa%write_SQt .OR. propa%write_Pt .OR. &
                   propa%write_alpha .OR. propa%write_validation
      need_energy = propa%write_E .OR. propa%write_validation
      need_norm = propa%write_Norm .OR. propa%write_validation
      need_autocorrelation = propa%write_auto_cor .OR. propa%write_validation

      ratio = (propa%tf - propa%t0) / propa%delta_t
      IF (ABS(ratio) > REAL(HUGE(total_steps), KIND=Rkind)) ERROR STOP 'propagation: requested number of steps exceeds integer range.'
      total_steps = NINT(ratio, KIND=Ikind)
      time_tol = ONETENTH**9 * MAX(ONE, ABS(ratio))
      IF (ABS(ratio - REAL(total_steps,KIND=Rkind)) > time_tol) THEN
         WRITE(out_unit,*) 'tf-t0 must be an integer multiple of delta_t: ratio=', ratio
         ERROR STOP 'propagation: final time is not aligned with the time step.'
      END IF

      checkpoint_config%t0 = propa%t0
      checkpoint_config%delta_t = propa%delta_t
      checkpoint_config%eps = propa%eps
      checkpoint_config%E_min = propa%E_min
      checkpoint_config%E_max = propa%E_max
      checkpoint_config%max_iter = propa%max_iter
      checkpoint_config%split_order = propa%split_order
      checkpoint_config%Beta = propa%Beta
      checkpoint_config%P = propa%P
      checkpoint_config%renorm = propa%renorm
      checkpoint_config%propa_name = propa%propa_name
      checkpoint_config%integrator = name
      checkpoint_config%input_signature = Input_Signature()

      current_step = 0_Ikind
      current_time = propa%t0
      Norm_ref = ZERO; E_ref = ZERO
      max_norm_drift = ZERO; max_energy_drift = ZERO; max_abs_dQ = ZERO; max_abs_dP = ZERO
      IF (propa%restart) THEN
         CALL Load_Propagation_Checkpoint(TRIM(propa%checkpoint_file), psi, psi_t0, checkpoint_config, &
              current_time, current_step, stats, Qt_ref, Pt_ref, ierr, checkpoint_message)
         IF (ierr /= 0_Ikind) THEN
            WRITE(out_unit,'(a)') 'Restart refused: ' // TRIM(checkpoint_message)
            ERROR STOP 'propagation: checkpoint compatibility/read failure.'
         END IF
         IF (current_step < 0_Ikind .OR. current_step > total_steps) &
            ERROR STOP 'propagation: checkpoint step is outside the requested [t0,tf] interval.'
         IF (ABS(current_time - (propa%t0 + REAL(current_step,KIND=Rkind)*propa%delta_t)) > &
             ONETENTH**9*MAX(ONE,ABS(current_time))) &
            ERROR STOP 'propagation: checkpoint time and saved step are inconsistent.'
         IF (current_time > propa%tf + ONETENTH**9*MAX(ONE,ABS(propa%tf))) &
            ERROR STOP 'propagation: requested tf is earlier than the checkpoint time.'

         Norm_ref = stats(1); E_ref = stats(2)
         max_norm_drift = stats(3); max_energy_drift = stats(4)
         max_abs_dQ = stats(5); max_abs_dP = stats(6)
         IF (propa%propa_name == 'hagedorn') THEN
            CALL construct_primitive_basis(psi%Basis)
            psi_dt%Basis => psi%Basis
            psi_dt%CVec(:) = psi%CVec(:)
            IF (ALLOCATED(H%Scal_pot)) DEALLOCATE(H%Scal_pot)
            CALL Calc_tab_Iq(Tab_iq, psi%Basis)
            CALL Set_Op(H, psi%Basis, Tab_iq, nio=in_unit)
         END IF
         WRITE(out_unit,*) 'Restarting at t=', current_time, ' (step ', current_step, ') of ', total_steps
         CALL Trim_Propagation_Outputs(current_time, ierr, checkpoint_message)
         IF (ierr /= 0_Ikind) THEN
            WRITE(out_unit,'(a)') 'Unable to reconcile observable files: ' // TRIM(checkpoint_message)
            ERROR STOP 'propagation: cannot trim output history for restart.'
         END IF
      ELSE IF (propa%write_validation) THEN
         IF (name == 'vp') THEN
            CALL Get_Basis_Parameters(psi%Basis, Qt_ref, SQt, At, Pt_ref)
         ELSE
            CALL Calc_Basis_parameters(psi, Qt_ref, SQt, At, Pt_ref)
         END IF
         CALL Calc_Av_E(E_ref, psi, H)
         CALL Calc_Norm_OF_Psi(psi, Norm_ref)
      END IF

      append_outputs = propa%restart
      IF (propa%write_Qt) CALL open_propa_file(u_Qt, 'Qt', append=append_outputs)
      IF (propa%write_E) CALL open_propa_file(u_E, 'E', append=append_outputs)
      IF (propa%write_SQt) CALL open_propa_file(u_SQt, 'SQt', append=append_outputs)
      IF (propa%write_Norm) CALL open_propa_file(u_Norm, 'Norm', append=append_outputs)
      IF (propa%write_pop) CALL open_propa_file(u_pop, 'pop', append=append_outputs)
      IF (propa%write_Pt) CALL open_propa_file(u_Pt, 'Imp_k', append=append_outputs)
      IF (propa%write_alpha) CALL open_propa_file(u_alpha, 'alpha', append=append_outputs)
      IF (propa%write_auto_cor) CALL open_propa_file(u_auto_cor, 'auto_cor', append=append_outputs)
      IF (propa%write_E_13) CALL open_propa_file(u_E_13, 'E_13', append=append_outputs)
      IF (propa%write_Norm_13) CALL open_propa_file(u_Norm_13, 'Norm_13', append=append_outputs)
      IF (propa%write_validation) THEN
         validation_header_needed = .NOT. append_outputs
         IF (append_outputs) THEN
            INQUIRE(FILE='validation.txt', EXIST=file_exists, SIZE=file_size)
            validation_header_needed = .NOT. file_exists
            IF (file_exists) validation_header_needed = (file_size == 0_Ikind)
         END IF
         CALL open_propa_file(u_validation, 'validation', append=append_outputs)
         IF (validation_header_needed) &
            WRITE(u_validation,'(a)') '# t norm dnorm energy dE max_abs_dQ max_abs_dP abs_auto_minus_one'
      END IF

      IF (name == 'itp') THEN
         CALL relax_propagation(psi_t0, H, propa)
      ELSE IF (name == 'gssil') THEN
         CALL find_ground_state_SIL(E0, Gspsi, H, psi0, psi0%Basis, propa%max_iter, propa%eps)
         IF (propa%debug) WRITE(out_unit, *) 'Ground state energy found: ', E0
         CALL Calc_Av_E(E0, Gspsi, H)
         IF (propa%debug) WRITE(out_unit,*) '<Gspsi|H|Gspsi> = ', E0
      ELSE
         IF (.NOT. propa%restart) THEN
            CALL Record_Observables(current_time)
            IF (propa%checkpoint) CALL Save_Current_Checkpoint(current_time, current_step)
         END IF

         DO i = current_step + 1_Ikind, total_steps
            t = propa%t0 + REAL(i,KIND=Rkind)*propa%delta_t
            CALL march_Global(psi, psi_dt, t, propa, H, u_Norm_13, u_E_13)
            current_step = i
            current_time = t
            CALL Record_Observables(current_time)
            IF (propa%checkpoint) THEN
               IF (MOD(current_step,propa%checkpoint_every) == 0_Ikind .OR. current_step == total_steps) &
                  CALL Save_Current_Checkpoint(current_time, current_step)
            END IF
         END DO

         psif%CVec(:) = psi%CVec(:)
         psif%Basis => psi%Basis
         Norm = ZERO
         IF (propa%debug .OR. propa%write_Norm .OR. propa%write_validation) CALL Calc_Norm_OF_Psi(psif, Norm)

         IF (propa%write_validation) THEN
            WRITE(u_validation,'(a)') '# summary'
            WRITE(u_validation,'(a,1x,es24.16)') 'reference_norm', Norm_ref
            WRITE(u_validation,'(a,1x,es24.16)') 'reference_energy', E_ref
            WRITE(u_validation,'(a,1x,es24.16)') 'max_norm_drift', max_norm_drift
            WRITE(u_validation,'(a,1x,es24.16)') 'max_energy_drift', max_energy_drift
            WRITE(u_validation,'(a,1x,es24.16)') 'tolerance_norm', validation_tol_norm
            WRITE(u_validation,'(a,1x,es24.16)') 'tolerance_energy_hartree', validation_tol_energy
            WRITE(u_validation,'(a,1x,es24.16)') 'max_abs_dQ', max_abs_dQ
            WRITE(u_validation,'(a,1x,es24.16)') 'max_abs_dP', max_abs_dP
            WRITE(u_validation,'(a,1x,l1)') 'norm_within_tolerance', max_norm_drift <= validation_tol_norm
            WRITE(u_validation,'(a,1x,l1)') 'energy_within_tolerance', max_energy_drift <= validation_tol_energy
         END IF
      END IF

      IF (propa%debug) THEN
         WRITE(out_unit, *) 'END propagation; final norm = ', Norm
         FLUSH(out_unit)
      END IF
      IF (propa%write_validation) CLOSE(u_validation)
      IF (propa%write_Qt) CLOSE(u_Qt)
      IF (propa%write_E) CLOSE(u_E)
      IF (propa%write_SQt) CLOSE(u_SQt)
      IF (propa%write_Norm) CLOSE(u_Norm)
      IF (propa%write_pop) CLOSE(u_pop)
      IF (propa%write_Pt) CLOSE(u_Pt)
      IF (propa%write_alpha) CLOSE(u_alpha)
      IF (propa%write_auto_cor) CLOSE(u_auto_cor)
      IF (propa%write_E_13) CLOSE(u_E_13)
      IF (propa%write_Norm_13) CLOSE(u_Norm_13)

      CALL dealloc_psi(psi)
      CALL dealloc_psi(psi_dt)
      CALL dealloc_psi(psi_t0)
      CALL dealloc_psi(Gspsi)
      IF (ALLOCATED(Tab_iq)) DEALLOCATE(Tab_iq)
      DEALLOCATE(Qt, SQt, Pt, Qt_ref, Pt_ref, At, pop)

   CONTAINS

      SUBROUTINE Record_Observables(time_now)
         REAL(KIND=Rkind), INTENT(IN) :: time_now

         E = ZERO
         Norm = ZERO
         IF (need_basis) THEN
            IF (name == 'vp') THEN
               CALL Get_Basis_Parameters(psi%Basis, Qt, SQt, At, Pt)
            ELSE IF (propa%write_validation) THEN
               CALL Calc_Basis_parameters(psi, Qt, SQt, At, Pt)
            ELSE
               CALL Calc_Basis_parameters(psi, Qt, SQt, At, Pt, &
                    need_Pt=propa%write_Pt, need_At=propa%write_alpha)
            END IF
         END IF
         IF (need_energy) CALL Calc_Av_E(E, psi, H)
         IF (need_norm) CALL Calc_Norm_OF_Psi(psi, Norm)
         IF (propa%write_pop) CALL Population(psi, pop)
         IF (need_autocorrelation) &
            CALL Calc_Auto_corr(psi_t0, psi, aut_func(1), aut_func_arg(1), propa, &
                                t=time_now, debug=propa%debug)

         IF (propa%write_auto_cor) CALL Write_Complex(aut_func, nio=u_auto_cor, op='t_z', t=time_now)
         IF (propa%write_alpha) CALL Write_Complex(At, nio=u_alpha, op='t_z', t=time_now)
         IF (propa%write_Qt) WRITE(u_Qt,*) time_now, Qt
         IF (propa%write_E) WRITE(u_E,*) time_now, E
         IF (propa%write_SQt) WRITE(u_SQt,*) time_now, SQt
         IF (propa%write_Norm) WRITE(u_Norm,*) time_now, Norm
         IF (propa%write_pop) WRITE(u_pop,*) time_now, pop
         IF (propa%write_Pt) WRITE(u_Pt,*) time_now, Pt

         IF (propa%write_validation) THEN
            step_max_abs_dQ = MAXVAL(ABS(Qt - Qt_ref))
            step_max_abs_dP = MAXVAL(ABS(Pt - Pt_ref))
            max_norm_drift = MAX(max_norm_drift, ABS(Norm - Norm_ref))
            max_energy_drift = MAX(max_energy_drift, ABS(E - E_ref))
            max_abs_dQ = MAX(max_abs_dQ, step_max_abs_dQ)
            max_abs_dP = MAX(max_abs_dP, step_max_abs_dP)
            auto_mod_dev = ABS(ABS(aut_func(1)) - ONE)
            WRITE(u_validation,'(8(1x,es24.16))') time_now, Norm, ABS(Norm-Norm_ref), E, &
                 ABS(E-E_ref), step_max_abs_dQ, step_max_abs_dP, auto_mod_dev
         END IF

         IF (propa%write_Qt) FLUSH(u_Qt)
         IF (propa%write_E) FLUSH(u_E)
         IF (propa%write_SQt) FLUSH(u_SQt)
         IF (propa%write_Norm) FLUSH(u_Norm)
         IF (propa%write_pop) FLUSH(u_pop)
         IF (propa%write_Pt) FLUSH(u_Pt)
         IF (propa%write_alpha) FLUSH(u_alpha)
         IF (propa%write_auto_cor) FLUSH(u_auto_cor)
         IF (propa%write_validation) FLUSH(u_validation)
      END SUBROUTINE Record_Observables


      SUBROUTINE Save_Current_Checkpoint(time_now, saved_step)
         REAL(KIND=Rkind), INTENT(IN) :: time_now
         INTEGER(KIND=Ikind), INTENT(IN) :: saved_step
         stats = [Norm_ref, E_ref, max_norm_drift, max_energy_drift, max_abs_dQ, max_abs_dP]
         CALL Save_Propagation_Checkpoint(TRIM(propa%checkpoint_file), time_now, saved_step, &
              psi, psi_t0, stats, Qt_ref, Pt_ref, checkpoint_config, ierr, checkpoint_message)
         IF (ierr /= 0_Ikind) THEN
            WRITE(out_unit,'(a)') 'Checkpoint write failed: ' // TRIM(checkpoint_message)
            ERROR STOP 'propagation: failed to save checkpoint.'
         END IF
      END SUBROUTINE Save_Current_Checkpoint

   END SUBROUTINE propagation

END MODULE Propa_m
