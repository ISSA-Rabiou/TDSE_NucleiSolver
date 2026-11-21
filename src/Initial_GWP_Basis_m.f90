MODULE Initial_GWP_Basis_m
   USE QDUtil_m, ONLY: Rkind, ZERO, TWO, out_unit
   USE Basis_m, ONLY: Basis_t
   USE param_WP0_m, ONLY: GWP_t
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: Apply_Initial_GWP_To_Basis
CONTAINS
   SUBROUTINE Apply_Initial_GWP_To_Basis(Basis, GWP, propagation_debug)
      TYPE(Basis_t), INTENT(INOUT) :: Basis
      TYPE(GWP_t), INTENT(IN) :: GWP
      LOGICAL, INTENT(IN) :: propagation_debug
      INTEGER :: ndim, idim
      REAL(KIND=Rkind) :: alpha_real
      COMPLEX(KIND=Rkind) :: alpha

      IF (.NOT. ALLOCATED(Basis%tab_basis)) RETURN
      ndim = SIZE(Basis%tab_basis) - 1
      IF (GWP%ndim /= ndim) ERROR STOP 'Initial GWP dimension does not match the nuclear basis.'

      DO idim = 1, ndim
         IF (.NOT. Basis%tab_basis(idim)%use_initial_gwp) CYCLE
         IF (Basis%tab_basis(idim)%Basis_name /= 'ho' .AND. &
             Basis%tab_basis(idim)%Basis_name /= 'hag') THEN
            ERROR STOP 'use_initial_gwp is only valid for an HO/Hagedorn primitive basis.'
         END IF
         IF (GWP%tab_GWP1D(idim)%sigma <= ZERO) &
            ERROR STOP 'Initial GWP sigma must be positive.'

         alpha_real = TWO / GWP%tab_GWP1D(idim)%sigma**2
         alpha = CMPLX(alpha_real, TWO*GWP%tab_GWP1D(idim)%Beta, KIND=Rkind)
         IF (REAL(alpha, KIND=Rkind) <= ZERO) &
            ERROR STOP 'Derived Hagedorn Alpha must have a positive real part.'

         Basis%tab_basis(idim)%Q0 = GWP%tab_GWP1D(idim)%Q0
         Basis%tab_basis(idim)%Imp_k = GWP%tab_GWP1D(idim)%imp_k
         Basis%tab_basis(idim)%Alpha = alpha
         Basis%tab_basis(idim)%scaleQ = SQRT(REAL(alpha, KIND=Rkind))

         IF (Basis%tab_basis(idim)%debug .OR. propagation_debug) THEN
            WRITE(out_unit,*) 'Derived basis parameters from initial GWP, dimension ', idim
            WRITE(out_unit,*) 'Q0, scaleQ, Imp_k, Alpha = ', &
                 Basis%tab_basis(idim)%Q0, Basis%tab_basis(idim)%scaleQ, &
                 Basis%tab_basis(idim)%Imp_k, Basis%tab_basis(idim)%Alpha
         END IF
      END DO
   END SUBROUTINE Apply_Initial_GWP_To_Basis
END MODULE Initial_GWP_Basis_m
