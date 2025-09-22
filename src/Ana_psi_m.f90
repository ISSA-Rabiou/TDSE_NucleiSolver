module Ana_psi_m
    USE Basis_m
    USE NDindex_m
    USE psi_m
    implicit none
    private
    public:: Population, Qpop
    public::  Calc_Av_imp_k_nD
    public :: Calc_Avg_A_nD,test_analyse_psi,Calc_AVQ_SQ_nD
 
 contains
 
 
    SUBROUTINE test_analyse_psi(psi)
       USE  QDUtil_m
       USE Basis_m
       type(psi_t)    ,intent(in)                    :: psi
 
       real(kind=Rkind), allocatable                 :: V(:, :, :)
       integer, allocatable                          :: Tab_Iq(:, :)
       real(kind=Rkind) ,allocatable                 :: SQt(:),Qt(:)
       complex(kind=Rkind) ,allocatable              :: At(:)
       integer                                       :: ndim 
      
 
       ndim = size(psi%Basis%tab_basis) - 1
       allocate(Qt(ndim), SQt(ndim),At(ndim))
    
       call Calc_tab_Iq0(Tab_Iq,psi%Basis)
       call Calc_Avg_A_nD(psi, At)
 
         write(out_unit,*) 'At',At
         write(out_unit,*) '=============================================================================='
    
    END SUBROUTINE 
    
    subroutine Population(Psi, Pop)
       implicit none
       type(Psi_t), intent(in), target                 :: Psi
       complex(kind=Rkind), pointer                    :: Psi_bb(:, :)
       real(Kind=Rkind), intent(inout), allocatable    ::Pop(:)
       integer                                         :: inb,nsurf,ndim,nb
       real(Kind=Rkind)                                :: Norm
 
 
       ndim = size(Psi%Basis%tab_basis)-1
       nb = Psi%Basis%nb
       nsurf = Psi%Basis%tab_basis(ndim+1)%nb
       Psi_bb(1:nb, 1:nsurf) => Psi%CVec
       !call Calc_Norm_OF_Psi(Psi, Norm)
 
       do inb = 1, nsurf
        ! Norm = DOT_PRODUCT(Psi_bb(:, inb), Psi_bb(:, inb))
         Pop(inb) = DOT_PRODUCT(Psi_bb(:, inb), Psi_bb(:, inb))!/Norm
       end do
 
    end subroutine 
 
    SUBROUTINE Qpop(Psi, Qp)
       USE Basis_m
       USE QDUtil_m
       type(Psi_t), intent(in), target               :: Psi
       type(Psi_t), target                           :: Psi_g
       complex(kind=Rkind), pointer                     :: psi_gb(:, :)
 
       real(kind=Rkind), intent(inout)                  :: Qp(:)
       real(kind=Rkind)                                 :: Norm(2)
       logical, parameter                            :: debug = .true.
       integer                                       :: iq, inb, ndim
 
       IF (debug) THEN
          !write(out_unit,*) 'BEGINNING Qpop'
          flush (out_unit)
       END IF
       Ndim = size(Psi%Basis%tab_basis)
       call init_psi(psi_g, psi%Basis, cplx=.TRUE., grid=.true.)
       Psi_g%CVec(:) = CZERO
       call BasisTOGrid_nD_cplx(Psi_g%CVec, Psi%CVec, Psi%Basis)
       do inb = 1, Psi%Basis%tab_basis(2)%nb
          psi_gb(1:Psi%Basis%tab_basis(1)%nq, 1:Psi%Basis%tab_basis(2)%nb) => psi_g%CVec
          Qp(inb) = dot_product(psi_gb(:, inb), Psi%Basis%tab_basis(1)%w*Psi%Basis%tab_basis(1)%x*psi_gb(:, inb))
          Norm(inb) = dot_product(psi_gb(:, inb), Psi%Basis%tab_basis(1)%w*psi_gb(:, inb))
          if (Norm(inb) /= ZERO) then
             Qp(inb) = Qp(inb)/Norm(inb)
          end if
       end do
       !Do iq = 1,Psi%Basis%tab_basis(1)%nq
       !   write(666,*)   psi_gb(iq,1) , psi_gb(iq,2)
       !End Do
 
       print *, Qp, Norm
       IF (debug) THEN
          !        write(out_unit,*) 'END Qpop
          flush (out_unit)
       END IF
    END SUBROUTINE Qpop

    !=====================================================================
!> @brief
!>   Computes the averaged complex quantity At for all spatial dimensions
!>   of the wavefunction `psi`.
!>
!> @details
!>   - Computes VQ, VP, VQQ, VQP, and SQt using dedicated nD routines.
!>   - Constructs At(Ib) = CA(Ib) - i*CB(Ib) per dimension.
!>   - Local debug flag prints intermediate results if enabled.
!=====================================================================
   SUBROUTINE Calc_Avg_A_nD(psi, At)
      USE QDUtil_m
      IMPLICIT NONE
  
      !----------------- Arguments -----------------
      TYPE(psi_t),          INTENT(IN)    :: psi     ! Wavefunction with basis info
      COMPLEX(KIND=Rkind),  INTENT(INOUT) :: At(:)   ! Output averaged complex array
  
      !----------------- Locals -----------------
      LOGICAL, PARAMETER      :: debug = .FALSE.
      INTEGER                  :: Ndim, Ib
      REAL(KIND=Rkind), ALLOCATABLE :: VQ(:), VP(:), VQQ(:), VQP(:), SQt(:)
      REAL(KIND=Rkind), ALLOCATABLE :: CA(:), CB(:)
      INTEGER, ALLOCATABLE         :: Tab_Iq(:, :)
  
      !----------------- Debug output -----------------
      IF (debug) THEN
          FLUSH(out_unit)
      END IF
  
      !----------------- Determine number of spatial dimensions -----------------
      Ndim = SIZE(psi%Basis%tab_basis) - 1
  
      !----------------- Initialize quadrature indices -----------------
      CALL Calc_tab_Iq0(Tab_Iq, psi%Basis)
  
      !----------------- Allocate local arrays -----------------
      ALLOCATE(VQ(Ndim), VP(Ndim), VQQ(Ndim), VQP(Ndim), SQt(Ndim))
      ALLOCATE(CA(Ndim), CB(Ndim))
  
      VQ(:) = ZERO
      VP(:) = ZERO
      VQQ(:) = ZERO
      VQP(:) = ZERO
      SQt(:) = ONE
      CA(:) = ZERO
      CB(:) = ZERO
  
      !----------------- Compute necessary quantities -----------------
      CALL Calc_AVQ_SQ_nD_scd(psi, VQ, SQt, VQQ, Tab_Iq)  ! Average Q and its square
      CALL Calc_Av_imp_k_nD(psi, VP)                     ! Average momentum
      CALL Calc_VQP_nD(VQP, psi)                         ! VQP per dimension
  
      !----------------- Compute At array -----------------
      At(:) = CZERO
      DO Ib = 1, Ndim
          CA(Ib) = SQt(Ib)**2
          CB(Ib) = (VQP(Ib) - TWO*VP(Ib)*VQ(Ib)) / (TWO*(VQQ(Ib) - VQ(Ib)**2))
          At(Ib) = CMPLX(CA(Ib), -CB(Ib), KIND=Rkind)
      END DO
  
      !----------------- Debug final output -----------------
      IF (debug) THEN
          WRITE(out_unit,*) 'At = ', At
          FLUSH(out_unit)
      END IF
  
      !----------------- Deallocate local arrays -----------------
      DEALLOCATE(VQ, VP, VQQ, VQP, SQt, CA, CB, Tab_Iq)
  
  END SUBROUTINE Calc_Avg_A_nD

 
 
    SUBROUTINE Calc_Av_imp_k_1D(psi0, K, nio)
       USE QDUtil_m
       type(psi_t), intent(in), target                     :: psi0
       real(kind=Rkind), intent(inout)                     :: K
       integer, intent(in)                                 :: nio
       !locals variables---------------------------------------------------------
       type(psi_t), target                                 :: psi, ikpsi
       type(psi_t), target                                 :: psi_b, ikpsi_b
       logical, parameter                                  :: debug = .FALSE.
       complex(kind=Rkind), pointer                        :: GB(:,:,:)
       complex(kind=Rkind), pointer                        :: d1gg(:,:,:)
       complex(kind=Rkind), pointer                        :: ikpsi0(:,:,:)
       Integer, allocatable                                :: Iq1(:), Iq2(:), Iq3(:)
       Integer, allocatable                                :: Ib1(:), Ib2(:), Ib3(:)
       integer                                             :: nq,inb,Ndim,Inbe,i1,i3
       !debuging----------------------------------------------------------------
       
       IF (debug) THEN
          flush (out_unit)
       END IF
       
       Ndim = size(psi0%Basis%tab_basis)
       call init_psi(psi, psi0%Basis, cplx=.TRUE., grid=.TRUE.)
       call init_psi(ikpsi, psi0%Basis, cplx=.TRUE., grid=.TRUE.)
       call init_psi(psi_b, psi0%Basis, cplx=.TRUE., grid=.FALSE.)
       call init_psi(ikpsi_b, psi0%Basis, cplx=.TRUE., grid=.FALSE.)
       Call Calc_index(Ib1=Ib1, Ib2=Ib2, Ib3=Ib3,Iq1=Iq1, Iq2=Iq2,&
        & Iq3=Iq3, Basis=psi0%Basis)

       psi%CVec = CZERO
       ikpsi%CVec = CZERO
       ikpsi_b%CVec = CZERO
       psi_b%CVec = CZERO

       IF (psi0%Grid) THEN
          psi%CVec(:) = psi0%CVec(:)
       ELSE
          call BasisTOGrid_nD_cplx(psi%CVec, psi0%CVec, psi0%Basis)
       END IF
        
         GB(1:Iq1(nio),1:Iq2(nio),1:Iq3(nio))    => psi%CVec
         ikpsi0(1:Iq1(nio),1:Iq2(nio),1:Iq3(nio))=> ikpsi%CVec
         d1gg(1:Iq2(nio),1:Iq2(nio),1:1)=> psi%Basis%tab_basis(nio)%d1gg(:, :, 1)
         
         DO i3 = 1, UBOUND(GB, dim=3)
           DO i1 = 1, UBOUND(GB, dim=1)
             ikpsi0(i1, :, i3) = ikpsi0(i1, :, i3)- EYE*MATMUL(d1gg(:,:,1),GB(i1,:,i3))
           END DO
         END DO
       call GridTOBasis_nD_cplx(psi_b%CVec, psi%CVec, psi0%Basis)
       call GridTOBasis_nD_cplx(ikpsi_b%CVec, ikpsi%CVec, psi0%Basis)
       K = dot_product(psi_b%CVec, ikpsi_b%CVec)/dot_product(psi_b%CVec, psi_b%CVec)

       IF (debug) THEN
          write (out_unit, *) '<psi/-id_x/psi> =', K
          flush (out_unit)
       END IF
       
      call dealloc_psi(psi)
      call dealloc_psi(ikpsi)
      call dealloc_psi(ikpsi_b)
      call dealloc_psi(psi_b)
      Deallocate(Ib1, Ib2, Ib3,Iq1, Iq2, Iq3)
    END SUBROUTINE

   !=====================================================================
   !> @brief
   !>   Computes the nD average impulse <psi|p_i|psi> along each dimension.
   !>
   !> @details
   !>   - Loops over all spatial dimensions and calls 1D impulse routine.
   !>   - Ensures psi is represented on the grid.
   !>   - Optionally prints debug information.
   !>
   !> @param[in]  psi0  Input wavefunction (Basis or Grid representation)
   !> @param[out] K     Array of average impulse values along each dimension
   !> @param[in]  debug [OPTIONAL] enable verbose prints
   !=====================================================================
   SUBROUTINE Calc_Av_imp_k_nD(psi0, K)
      USE QDUtil_m
      IMPLICIT NONE
  
      !----------------- Arguments -----------------
      TYPE(psi_t),       INTENT(IN)           :: psi0
      REAL(KIND=Rkind),  INTENT(INOUT)        :: K(:)
  
      !----------------- Locals --------------------
      TYPE(psi_t), TARGET                 :: psi
      LOGICAL, PARAMETER                  :: debug = .FALSE.  ! local debug
      INTEGER                             :: Ndim, Inb
  
      !----------------- Debug output ----------------
      IF (debug) THEN
          FLUSH(out_unit)
          !WRITE(out_unit,*) '>>> Calc_Av_imp_k_nD start'
      END IF
  
      !----------------- Setup ---------------------
      Ndim = SIZE(psi0%Basis%tab_basis) - 1
      CALL init_psi(psi, psi0%Basis, cplx=.TRUE., grid=.TRUE.)
  
      ! Zero out wavefunction and output
      psi%CVec(:) = CZERO
      K(:)        = ZERO
  
      ! Ensure psi is on the grid
      IF (psi0%Grid) THEN
          psi%CVec(:) = psi0%CVec(:)
      ELSE
          CALL BasisTOGrid_nD_cplx(psi%CVec, psi0%CVec, psi0%Basis)
      END IF
  
      !----------------- Loop over spatial dimensions -----------------
      DO Inb = 1, Ndim
          CALL Calc_Av_imp_k_1D(psi, K(Inb), Inb)
      END DO
  
      ! Optional debug output
      IF (debug) THEN
         WRITE(out_unit,*) ''
          WRITE(out_unit,*) '<psi|-id_xi/psi> =', K
          !WRITE(out_unit,*) '<<< Calc_Av_imp_k_nD done'
          FLUSH(out_unit)
      END IF
  
      !----------------- Cleanup -------------------
      CALL dealloc_psi(psi)
  
  END SUBROUTINE Calc_Av_imp_k_nD 

 
  SUBROUTINE Calc_VQP_1D(psi0, VQP, Ib)
     USE QDUtil_m
     type(psi_t), intent(in), target              :: psi0
     real(kind=Rkind), intent(inout)              :: VQP
     integer, intent(in)                          :: Ib
     !locals variables----------------------      ---------------
     type(psi_t), target                           :: psi, d1psi
     logical, parameter                            :: debug = .true.
     complex(kind=Rkind), pointer                  :: GB(:,:,:)
     complex(kind=Rkind), pointer                  :: d1gg(:,:,:)
     complex(kind=Rkind), pointer                  :: d1psi0(:,:,:)
     Integer, allocatable                          :: Ib1(:), Ib2(:), Ib3(:)
     Integer, allocatable                          :: Iq1(:), Iq2(:), Iq3(:)
     integer                                       :: nq,inb,Ndim,Inbe,i1,i3,Iq,nsurf
 
     complex(kind=Rkind), allocatable              :: psi_gb(:, :),d1psi_gb(:, :)
     logical                                       :: Endloop_q
     real(kind=Rkind), allocatable                 :: VQPEl(:)
     real(kind=Rkind)                              :: W, X
     real(kind=Rkind), allocatable                 :: N(:)
     integer, allocatable                          :: Tab_iq(:)
     !debuging--------------------------------------------------
     
     IF (debug) THEN
        flush (out_unit)
     END IF
     
     Ndim = size(psi0%Basis%tab_basis)-1
     nsurf = psi0%Basis%tab_basis(Ndim+1)%nb
     call init_psi(psi,   psi0%Basis, cplx=.true., grid=.true.)
     call init_psi(d1psi, psi0%Basis, cplx=.true., grid=.true.)
     Call Calc_index(Ib1=Ib1, Ib2=Ib2, Ib3=Ib3,Iq1=Iq1, Iq2=Iq2,&
      & Iq3=Iq3, Basis=psi0%Basis)
 
     psi%CVec = CZERO
     d1psi%CVec = CZERO
 
     IF (psi0%Grid) then
        psi%CVec(:) = psi0%CVec(:)
     ELSE
        call BasisTOGrid_nD_cplx(psi%CVec, psi0%CVec, psi0%Basis)
     END IF
      
       GB(1:Iq1(Ib),1:Iq2(Ib),1:Iq3(Ib))    => psi%CVec
       d1psi0(1:Iq1(Ib),1:Iq2(Ib),1:Iq3(Ib))=> d1psi%CVec
       d1gg(1:Iq2(Ib),1:Iq2(Ib),1:1) => psi%Basis%tab_basis(Ib)%d1gg(:, :, 1)
       
       DO i3 = 1, ubound(GB, dim=3)
         DO i1 = 1, ubound(GB, dim=1)
           d1psi0(i1, :, i3) = d1psi0(i1, :, i3) + matmul(d1gg(:,:,1),GB(i1,:,i3))
         END DO
      END DO
      
      Allocate (psi_gb(psi%Basis%nq, nsurf))
      Allocate (d1psi_gb(d1psi%Basis%nq, nsurf))
      Allocate (VQPEl(nsurf))
      Allocate (Tab_iq(Ndim),N(nsurf))
 
      psi_gb(:, :) =reshape(psi%CVec, shape=[psi%Basis%nq, nsurf])
      d1psi_gb(:, :) =reshape(d1psi%CVec, shape=[d1psi%Basis%nq, d1psi%Basis%tab_basis(Ndim+1)%nb])
      X = ZERO
      N(:) = ZERO
 
     DO inbe = 1, nsurf !electronic state
      VQPel(inbe) = ZERO
      Call Init_tab_ind(Tab_iq, psi%Basis%NDindexq)
       Iq = 0
       DO
          Iq = Iq + 1
          CALL increase_NDindex(Tab_iq, psi%Basis%NDindexq, Endloop_q)
          IF (Endloop_q) exit
          W = ONE
          DO inb = 1,Ndim
             W = W*psi%Basis%tab_basis(inb)%w(tab_iq(inb))
          END DO
 
          X = psi%Basis%tab_basis(ib)%x(tab_iq(ib))
          N(inbe) = N(inbe) + REAL(conjg(psi_gb(iq, inbe))*psi_gb(iq, inbe)*W,kind=Rkind)
          VQPel(inbe) = VQPel(inbe)-REAL(EYE*conjg(psi_gb(iq, inbe))*W*(psi_gb(iq, inbe)+TWO*X*d1psi_gb(iq, inbe)),kind=Rkind)
 
       END DO
     END DO
      VQP = sum(VQPel)/(Sum(N))
     
 
       call dealloc_psi(psi)
       call dealloc_psi(d1psi)
       deallocate(Ib1, Ib2, Ib3,Iq1, Iq2, Iq3)
       deallocate(N,VQPEl,psi_gb,d1psi_gb,Tab_iq)
 
     IF (debug) THEN
        flush (out_unit)
     END IF
    
  END SUBROUTINE  
 
 
  SUBROUTINE Calc_VQP_nD(VQP,psi)
   TYPE(psi_t), intent(in)                        :: psi
   real(kind=Rkind) , intent(inout)               :: VQP(:)
 
   !Locals variabls ----------------------------------------------------------
 
   integer                                         :: Ib,ndim
 
     ndim  = size(psi%Basis%tab_basis) - 1
 
     Do Ib = 1,ndim
       call Calc_VQP_1D(psi, VQP(Ib), Ib)
     End do
     
     !write (out_unit, *) 'VQP',VQP
 
 END SUBROUTINE
 
 SUBROUTINE Calc_AVQ_SQ_nD(psi_in, AVQ, SQ,Tab_Iq)
    USE QDUtil_m
    TYPE(Psi_t), intent(in)                       :: psi_in
    integer,intent(in)                            :: Tab_Iq(:, :)
    real(kind=Rkind), intent(inout)               :: AVQ(:), SQ(:)
 
    logical, parameter                            :: debug = .false.
    TYPE(Psi_t)                                   :: psi
    complex(kind=Rkind), allocatable              :: psi_gb(:, :)
    real(kind=Rkind), allocatable                 :: AVQel(:, :), SQel(:, :),Q(:)
    real(kind=Rkind)                              :: W
    real(kind=Rkind), allocatable                 :: N(:)
    integer                                       :: iq, inbe, inb,nsurf,nq,ndim,Ib
 
 
    IF (debug) THEN
       write (out_unit, *) 'Beging AVQ'
       flush (out_unit)
    END IF
 
    ndim = size(psi_in%Basis%tab_basis) - 1
    nq =psi_in%Basis%nq
    nsurf=psi_in%Basis%tab_basis(ndim+1)%nb
    allocate (N(nsurf),AVQel(ndim, nsurf),SQel(ndim, nsurf),Q(ndim),psi_gb(nq, nsurf))
    CALL init_psi(psi, psi_in%Basis, cplx=.true., grid=.true.)
 
    IF (psi_in%Grid) then
       psi%CVec(:) = psi_in%CVec(:)
    ELSE
       CALL BasisTOGrid_nD_cplx(psi%CVec, psi_in%CVec, psi_in%Basis)
    END IF
 
    psi_gb(:, :) = reshape(psi%CVec, shape=[nq, nsurf])
    N(:) = ZERO
    AVQel(:, :) = ZERO
    SQel(:, :)  = ZERO
 
    DO inbe = 1, nsurf 
 
       DO Iq =  1,nq
          w = ONE
          Q(:)=ZERO
 
          DO inb = 1, ndim
             w = w*psi%Basis%tab_basis(inb)%w(Tab_Iq(inb,Iq))
             Q(inb)= psi%Basis%tab_basis(inb)%x(Tab_Iq(inb,Iq))
          END DO
          
            N(inbe) = N(inbe) + REAL(conjg(psi_gb(iq, inbe))*psi_gb(iq, inbe)*W,kind=Rkind)
            AVQel(:, inbe) = AVQel(:, inbe) + REAL(conjg(psi_gb(iq, inbe))*(Q(:)*w)*psi_gb(iq, inbe),kind=Rkind)
            SQel(:, inbe) = SQel(:, inbe) + REAL(conjg(psi_gb(iq, inbe))*(Q(:)*Q(:)*w)*psi_gb(iq, inbe),kind=Rkind)
 
       END DO
 
 
    END DO
 
    DO inb = 1, ndim
       AVQ(inb) = sum(AVQel(inb,:))/(Sum(N))
       SQ(inb) = sum(SQel(inb, :))/(Sum(N))
       SQ(inb) = TWO*(SQ(inb) - AVQ(inb)*AVQ(inb))
       SQ(inb) = sqrt(ONE/SQ(inb))
    
    END DO   
 
    Deallocate (Psi_gb,N,AVQel,SQel,Q)
    CALL dealloc_psi(psi)
    IF (debug) THEN
       write (out_unit, *) 'END AVQ'
       flush (out_unit)
    END IF
 END SUBROUTINE
 

 !=====================================================================
!> @brief
!>   Computes average Q (AVQ), standard deviation SQ, and VQQ along all dimensions.
!>
!> @details
!>   - Integrates over grid points with quadrature weights.
!>   - Computes <Q>, <Q^2> and variance for each dimension.
!>   - Normalizes by norm squared of wavefunction.
!>
!> @param[in]  psi_in  Input wavefunction (Basis or Grid)
!> @param[in]  Tab_Iq  ND indices for each grid point (ndim × nq)
!> @param[out] AVQ     Average Q per dimension
!> @param[out] SQ      Standard deviation per dimension
!> @param[out] VQQ     <Q^2> per dimension
!=====================================================================
SUBROUTINE Calc_AVQ_SQ_nD_scd(psi_in, AVQ, SQ, VQQ, Tab_Iq)
   USE QDUtil_m
   IMPLICIT NONE

   !----------------- Arguments -----------------
   TYPE(Psi_t), INTENT(IN)                     :: psi_in
   INTEGER, INTENT(IN)                         :: Tab_Iq(:, :)
   REAL(KIND=Rkind), INTENT(INOUT)            :: AVQ(:), SQ(:), VQQ(:)

   !----------------- Locals --------------------
   LOGICAL, PARAMETER                          :: debug = .FALSE.
   TYPE(Psi_t)                                 :: psi
   COMPLEX(KIND=Rkind), ALLOCATABLE           :: psi_gb(:, :)
   REAL(KIND=Rkind), ALLOCATABLE              :: AVQel(:, :), SQel(:, :), Q(:), N(:)
   REAL(KIND=Rkind)                            :: W
   INTEGER                                     :: iq, inbe, inb, nsurf, nq, ndim

   !----------------- Debug ---------------------
   IF (debug) THEN
       !WRITE(out_unit,*) '>>> Begin Calc_AVQ_SQ_nD_scd'
       FLUSH(out_unit)
   END IF

   !----------------- Setup ---------------------
   ndim  = SIZE(psi_in%Basis%tab_basis) - 1
   nq    = psi_in%Basis%nq
   nsurf = psi_in%Basis%tab_basis(ndim+1)%nb

   ALLOCATE(N(nsurf), AVQel(ndim, nsurf), SQel(ndim, nsurf), Q(ndim))
   ALLOCATE(psi_gb(nq, nsurf))
   CALL init_psi(psi, psi_in%Basis, cplx=.TRUE., grid=.TRUE.)

   ! Ensure psi on the grid
   IF (psi_in%Grid) THEN
       psi%CVec(:) = psi_in%CVec(:)
   ELSE
       CALL BasisTOGrid_nD_cplx(psi%CVec, psi_in%CVec, psi_in%Basis)
   END IF

   psi_gb(:, :) = RESHAPE(psi%CVec, SHAPE=[nq, nsurf])
   N(:)         = ZERO
   AVQel(:, :)  = ZERO
   SQel(:, :)   = ZERO

   !----------------- Loop over electronic states and grid points -----------------
   DO inbe = 1, nsurf
       DO iq = 1, nq
           W = ONE
           Q(:) = ZERO
           ! Compute product of weights and Q values along dimensions
           DO inb = 1, ndim
               W = W * psi%Basis%tab_basis(inb)%w(Tab_Iq(inb, iq))
               Q(inb) = psi%Basis%tab_basis(inb)%x(Tab_Iq(inb, iq))
           END DO

           ! Accumulate norm and moments (convert complex to real)
           N(inbe)       = N(inbe) + REAL(CONJG(psi_gb(iq, inbe)) * psi_gb(iq, inbe) * W, KIND=Rkind)
           AVQel(:, inbe) = AVQel(:, inbe) + REAL(CONJG(psi_gb(iq, inbe)) * (Q(:) * W) * psi_gb(iq, inbe), KIND=Rkind)
           SQel(:, inbe)  = SQel(:, inbe)  + REAL(CONJG(psi_gb(iq, inbe)) * (Q(:)**2 * W) * psi_gb(iq, inbe), KIND=Rkind)
       END DO
   END DO

   !----------------- Compute averages and standard deviations -----------------
   DO inb = 1, ndim
       AVQ(inb) = SUM(AVQel(inb, :)) / (SUM(N))
       VQQ(inb) = SUM(SQel(inb, :)) / (SUM(N))
       SQ(inb)  = SQRT(MAX(VQQ(inb) - AVQ(inb)**2, ZERO))   
       SQ(inb)  = ONE / (SQ(inb) * SQRT(TWO))                    
   END DO

   !----------------- Cleanup -------------------
   DEALLOCATE(psi_gb, N, AVQel, SQel, Q)
   CALL dealloc_psi(psi)

   IF (debug) THEN
      write(out_unit,*) 'AVQ= ', AVQ
      write(out_unit,*) 'SQ= ', SQ
      write(out_unit,*) 'VQQ= ', VQQ
      FLUSH(out_unit)
   END IF

END SUBROUTINE Calc_AVQ_SQ_nD_scd


 
 end module