MODULE NUFFT_m
   USE QDUtil_m
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: NUFFT, INV_NUFFT
CONTAINS
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

END MODULE NUFFT_m
