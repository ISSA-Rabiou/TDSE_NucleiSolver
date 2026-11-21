MODULE SplitOperator_m
   USE QDUtil_m
   USE NUFFT_m, ONLY: NUFFT, INV_NUFFT
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: ApplyExpVdtAndTransformToMomentum, ApplyExpTdtAndTransformToPosition
CONTAINS



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

END MODULE SplitOperator_m
