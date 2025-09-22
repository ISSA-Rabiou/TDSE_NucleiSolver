module Op_m
   !$ USE omp_lib
   USE QDUtil_m
   USE Basis_m, only: Basis_t
   USE Molec_m
   Use Psi_m
   implicit none
   private

   TYPE :: Op_t

      TYPE(Basis_t), pointer           :: Basis
      real(kind=Rkind), allocatable    :: RMat(:, :)
      complex(kind=Rkind), allocatable :: CMat(:, :)
      real(kind=Rkind), allocatable    :: Scal_pot(:, :, :)

   END TYPE Op_t

   PUBLIC :: Op_t, write_Op, Set_Op, dealloc_Op, calc_OpPsi, Calc_Hpsi, Kpsi_nD, Make_Mat_H
   PUBLIC :: test_op,Calc_Scalar_Pot, test_openmp_op,calc_tab_Iq,calc_nac,calc_VV,Popu
   PUBLIC :: Build_CN_Matrices, Build_CN_Propagator

CONTAINS
   SUBROUTINE alloc_Op(Op, nb)
      TYPE(Op_t), intent(inout) :: Op
      integer, intent(in)    :: nb

      CALL dealloc_Op(Op)

      IF (nb < 1) STOP 'ERROR in init_Op: nb < 1!'

      allocate (Op%RMat(nb, nb))

   END SUBROUTINE alloc_Op

   SUBROUTINE dealloc_Op(Op)
      TYPE(Op_t), intent(inout) :: Op

      nullify (Op%Basis)

      IF (ALLOCATED(Op%RMat)) THEN
         DEALLOCATE (Op%RMat)
      END IF
      IF (ALLOCATED(Op%CMat)) THEN
         DEALLOCATE (Op%CMat)
      END IF
      IF (ALLOCATED(Op%Scal_pot)) THEN
         DEALLOCATE (Op%Scal_pot)
      END IF

   END SUBROUTINE dealloc_Op

   SUBROUTINE write_Op(Op)
      USE QDUtil_m

      TYPE(Op_t), intent(in) :: Op

      ! integer :: i

      IF (associated(Op%Basis)) THEN
         write (out_unit, *) ' The basis is linked to Op.'
      END IF

      IF (allocated(Op%RMat)) THEN
         write (out_unit, *) 'Writing Op (real):'
         write (out_unit, *)
         CALL Write_VecMat(Op%RMat, out_unit, 5,  info='Op%Rmat')
         write (out_unit, *) 'END Writing Op'
      END IF

   END SUBROUTINE write_Op


!-----------------------------------------------------------------------
!> @brief  Build the Hamiltonian matrix in the chosen basis.
!>
!> This routine constructs the Hamiltonian matrix `H%CMat` using the basis
!> `Basis`. The Hamiltonian is applied to each basis vector via the operator
!> application routine, and the resulting matrix elements are obtained as
!> scalar products.
!>
!> Steps:
!>   1. Initialize temporary wavefunctions `psi` and `Hpsi`.
!>   2. Allocate the Hamiltonian matrix.
!>   3. Loop over basis vectors (columns):
!>        - Apply the operator H to a basis vector.
!>        - Compute the inner product with all other basis vectors.
!>   4. Store the results in `H%CMat`.
!>
!> A local debug flag is available to print information about the basis,
!> allocation, and progression of the computation.
!>
!> @param[in]    Basis   Basis set information
!> @param[inout] H       Hamiltonian operator (matrix will be filled)
!---------------------------------------------------------------------------------
SUBROUTINE Make_Mat_H(Basis,H)
   IMPLICIT NONE
 
   !--------- Arguments --------------------------------------------------
   TYPE(Basis_t), INTENT(IN)     :: Basis
   TYPE(Op_t),    INTENT(INOUT)  :: H
 
   !--------- Local variables ------------------------------------------
   LOGICAL, PARAMETER            :: debug = .FALSE.  ! Debug flag
   INTEGER, ALLOCATABLE          :: Tab_Iq(:, :)
   TYPE(psi_t)                   :: psi, Hpsi
   INTEGER                       :: nb, nsurf, ib, jb, ndim
 
   !--------- Debug: beginning -----------------------------------------
   IF (debug) THEN
      WRITE(out_unit, *) '>>> BEGINNING Make_Mat_H'
      CALL write_basis(Basis)
      FLUSH(out_unit)
   END IF
 
   !--------- Initialize temporary wavefunctions -------------------------
   CALL init_psi(psi,  Basis, cplx=.TRUE., grid=.FALSE.)
   CALL init_psi(Hpsi, Basis, cplx=.TRUE., grid=.FALSE.)
 
   !--------- Extract dimensions -----------------------------------------
   ndim  = SIZE(Basis%tab_basis) - 1
   nsurf = Basis%tab_basis(ndim+1)%nb
   nb    = (Basis%nb) * nsurf
 
   !--------- Set Op ---------------------------------------------------
   CALL Set_Op(H, Basis, Tab_Iq, alloc=.TRUE.)
 
   !--------- Loop over basis vectors ------------------------------------
   DO jb = 1, nb
      ! Set psi to basis vector |jb>
      psi%CVec(:)  = CZERO
      psi%CVec(jb) = CONE
 
      ! Apply Hamiltonian: Hpsi = H |psi>
      CALL calc_Oppsi(H, psi, Hpsi)
 
      ! Compute matrix elements <ib|H|jb>
      DO ib = 1, nb
         psi%CVec(:)  = CZERO
         psi%CVec(ib) = CONE
         H%CMat(ib, jb) = DOT_PRODUCT(psi%CVec, Hpsi%CVec)
      END DO
   END DO

   !Printing the matrix
   !call Write_VecMat(H%CMat, out_unit, 5, info='H')
 
   !--------- Cleanup --------------------------------------------------
   CALL dealloc_psi(psi)
   CALL dealloc_psi(Hpsi)
 
   !--------- Debug: end -----------------------------------------------
   IF (debug) THEN
      WRITE(out_unit, *) '>>> END Make_Mat_H'
      FLUSH(out_unit)
   END IF
 
 END SUBROUTINE Make_Mat_H



!--------------------------------------------------------------------
!> @brief  Construct the Crank–Nicolson matrices M_minus and M_plus.
!>
!> Given the Hamiltonian matrix H, this routine builds the two 
!> matrices required for Crank–Nicolson propagation:
!>      H_minus = I - (i dt / 2) H
!>      H_plus  = I + (i dt / 2) H
!>
!> @param[in]  H       Hamiltonian matrix (nb x nb)
!> @param[in]  dt      Time step Δt
!> @param[out] H_minus Matrix (nb x nb): I - (i dt / 2) H
!> @param[out] H_plus  Matrix (nb x nb): I + (i dt / 2) H
!--------------------------------------------------------------------
SUBROUTINE Build_CN_Matrices(H, dt, H_minus, H_plus)
   IMPLICIT NONE
 
   !------------ Arguments --------------------------------------------
   COMPLEX(KIND=Rkind), INTENT(IN)  :: H(:,:)
   REAL(KIND=Rkind),    INTENT(IN)  :: dt
   COMPLEX(KIND=Rkind), ALLOCATABLE, INTENT(OUT) :: H_minus(:,:), H_plus(:,:)
 
   !------------ Local variables ----------------------------------------
   LOGICAL, PARAMETER            :: debug = .FALSE.
   INTEGER                       :: nb, i
   COMPLEX(KIND=Rkind)           :: coef
 
   !------------ Get matrix size ----------------------------------------
   nb = SIZE(H, 1)
 
   !------------ Allocate outputs ---------------------------------------
   ALLOCATE(H_minus(nb, nb), H_plus(nb, nb))
 
   !------------ Compute prefactor coef = i dt / 2 ----------------------
   coef = EYE * dt / TWO
 
   !------------ Build matrices -----------------------------------------
   H_minus(:,:) = -coef * H(:,:)
   H_plus(:,:)  =  coef * H(:,:)
 
   ! Add identity contribution on the diagonal
   DO i = 1, nb
      H_minus(i,i) = H_minus(i,i) + CONE
      H_plus(i,i)  = H_plus(i,i)  + CONE
   END DO
 
   !------------ Debug info ---------------------------------------------
   IF (debug) THEN
      WRITE(out_unit,*) '>>> Build_CN_Matrices: matrices constructed'
      WRITE(out_unit,*) '    nb   =', nb
      WRITE(out_unit,*) '    dt   =', dt
      FLUSH(out_unit)
   END IF
 
 END SUBROUTINE Build_CN_Matrices


 !---------------------------------------------------------------------
 !> @brief  Construct the Crank–Nicolson propagator U.
 !>
 !> This routine builds the unitary propagator
 !>      U = (I - i dt/2 H)^(-1) (I + i dt/2 H)
 !>
 !> It calls Build_CN_Matrices to generate M_minus and M_plus,
 !> then inverts M_minus with inv_Of_Mat_TO.
 !>
 !> @param[in]  H   Hamiltonian matrix (nb x nb)
 !> @param[in]  dt  Time step Δt
 !> @param[out] U   Propagator matrix (nb x nb)
 !---------------------------------------------------------------------
 SUBROUTINE Build_CN_Propagator(H, dt, U)
    IMPLICIT NONE
 
    !-------- Arguments ------------------------------------------------
    COMPLEX(KIND=Rkind), INTENT(IN)  :: H(:,:)
    REAL(KIND=Rkind),    INTENT(IN)  :: dt
    COMPLEX(KIND=Rkind), ALLOCATABLE, INTENT(OUT) :: U(:,:)
 
    !-------- Local variables -------------------------------------------
    LOGICAL, PARAMETER :: debug = .FALSE.
    INTEGER            :: nb
    COMPLEX(KIND=Rkind), ALLOCATABLE :: M_minus(:,:), M_plus(:,:)
 
    !-------- Get matrix size -------------------------------------------
    nb = SIZE(H, 1)
 
    !-------- Build M_minus and M_plus ----------------------------------
    CALL Build_CN_Matrices(H, dt, M_minus, M_plus)
    IF (debug) CALL Write_VecMat(H, out_unit, 5, info='H')
 
    !-------- Invert M_minus --------------------------------------------
    M_minus = inv_OF_Mat(M_minus,inv_type=0,epsi=1.d-20)   ! M_minus = (M_minus)^(-1)
 
    !-------- Compute propagator ----------------------------------------
    ALLOCATE(U(nb, nb))
    U(:,:) = MATMUL(M_minus, M_plus)
 
    !-------- Debug info ------------------------------------------------
    IF (debug) THEN
       WRITE(out_unit,*) '>>> Build_Propagator: propagator constructed'
       WRITE(out_unit,*) '    nb   =', nb
       WRITE(out_unit,*) '    dt   =', dt
       FLUSH(out_unit)
    END IF


    !Print the propagator
    IF (debug) CALL Write_VecMat(U, out_unit, 5, info='U')
    !Cleanup
    IF (ALLOCATED(M_minus)) DEALLOCATE(M_minus)
    IF (ALLOCATED(M_plus))  DEALLOCATE(M_plus)
 
 END SUBROUTINE Build_CN_Propagator

!===============================================================================
!> SUBROUTINE: write_3D_array_column_major
!>
!> Purpose:
!>   Writes a 3D real array V(:,:,:) in column-major order.
!>   The output is linearized such that the first index varies fastest,
!>   then the second, then the third (Fortran-style storage order).
!>
!> Arguments:
!>   - V         : (IN)  3D real array
!>   - out_unit  : (IN)  Fortran output unit (screen or file)
!>
!> Notes:
!>   - Dimensions are automatically determined from V.
!>   - Each line contains one element with its linearized index and (i,j,k) indices.
!===============================================================================
 SUBROUTINE write_3D_array(V, out_unit)
   IMPLICIT NONE
   !----------------- Arguments -----------------
   REAL(KIND=Rkind), INTENT(IN) :: V(:,:,:)
   INTEGER, INTENT(IN) :: out_unit
 
   !----------------- Local variables -----------------
   INTEGER :: i, j, k, idx
   INTEGER :: dimx, dimy, dimz
 
   !----------------- Determine array size -----------------
   dimx = SIZE(V,1)
   dimy = SIZE(V,2)
   dimz = SIZE(V,3)
 
   WRITE(out_unit,*) '--- Writing 3D array in column-major order ---'
   WRITE(out_unit,*) 'Dimensions: (', dimx, ',', dimy, ',', dimz, ')'
 
   idx = 0
   DO k = 1, dimz
      DO j = 1, dimy
         DO i = 1, dimx
            idx = idx + 1
            WRITE(out_unit,'(A,I6,A,I4,A,I4,A,I4,A,F14.6)') 'Index ', idx, ': V(', i, ',', j, ',', k, ') = ', V(i,j,k)
         END DO
      END DO
   END DO
 
   !----------------- Flush output -----------------
   FLUSH(out_unit)
 END SUBROUTINE write_3D_array
 
!----------------------------------------------------------------------
!> @brief  Initialize or update the operator Op based on the given Basis.
!>
!> This routine sets up the operator `Op` using the information contained
!> in the basis `Basis` and the index table `Tab_Iq`. The behavior can be 
!> controlled with an optional logical argument `alloc`:
!>   - If `alloc=.TRUE.`, the operator structure is fully allocated 
!>     (memory allocation + initialization).
!>   - If `alloc=.FALSE.` (default), only the scalar potential `Scal_pot`
!>     is computed using the provided basis.
!>
!> A local debug flag is available to print important internal information
!> such as dimension sizes, number of surfaces, and allocation choices.
!>
!> @param[inout] Op      Operator object to initialize
!> @param[in]    Basis   Basis object defining the system
!> @param[in]    Tab_Iq  Index table for quadrature points
!> @param[in]    alloc   [optional] logical flag controlling allocation
!----------------------------------------------------------------------
SUBROUTINE Set_Op(Op, Basis, Tab_Iq, alloc)
   IMPLICIT NONE
 
   !--------- Arguments --------------------------------------------------
   TYPE(Op_t),     INTENT(INOUT)         :: Op
   TYPE(Basis_t),  INTENT(IN),   TARGET  :: Basis
   INTEGER,        INTENT(IN)            :: Tab_Iq(:, :)
   LOGICAL,        INTENT(IN), OPTIONAL  :: alloc
 
   !--------- Local variables ------------------------------------------
   INTEGER :: ndim, nsurf, nq,nb
   LOGICAL :: do_alloc
   LOGICAL , PARAMETER :: debug = .FALSE.
 
   !======== Check basis ===============================================
   IF (.NOT. Basis_IS_allocated(Basis)) THEN
      STOP 'ERROR in Set_Op: the Basis is not initialized'
   END IF
 
   !======== Extract dimensions ========================================
   ndim  = SIZE(Basis%tab_basis) - 1
   nq    = Basis%nq
   nsurf = Basis%tab_basis(ndim+1)%nb
   nb = Basis%nb*nsurf
 
   !======== Decide allocation mode ====================================
   IF (PRESENT(alloc)) THEN
      do_alloc = alloc
   ELSE
      do_alloc = .FALSE.
   END IF
 
   !======== Perform allocation or scalar potential calculation ========
   IF (do_alloc) THEN
      ! Full allocation of the operator structure
     IF(.NOT. ALLOCATED(Op%CMat)) ALLOCATE(Op%CMat(nb,nb))
     IF(.NOT. ALLOCATED(Op%Scal_pot)) ALLOCATE(Op%Scal_pot(nq, nsurf, nsurf))
     Op%Basis => Basis


     ! Initialize the operator matrix 
      Op%CMat(:,:) = CZERO
      Op%Scal_pot(:,:,:) = ZERO

     CALL Calc_Scalar_Pot(Op%Scal_pot, Basis)


     IF(debug) CALL Write_VecMat(Op%CMat, out_unit, 5, info='Op%CMat')
     IF(debug) CALL write_3D_array(Op%Scal_pot, out_unit)
 
      IF (debug) THEN
         WRITE(out_unit,*) '>>> Set_Op: Operator fully allocated'
         WRITE(out_unit,*) '    ndim   =', ndim
         WRITE(out_unit,*) '    nq     =', nq
         WRITE(out_unit,*) '    nb     =', Basis%nb
         WRITE(out_unit,*) '    nsurf  =', nsurf
      END IF
   ELSE
      ! Only compute scalar potential
      CALL Calc_Scalar_Pot(Op%Scal_pot, Basis)
 
      IF (debug) THEN
         WRITE(out_unit,*) '>>> Set_Op: Only scalar potential computed'
         WRITE(out_unit,*) '    ndim   =', ndim
         WRITE(out_unit,*) '    nq     =', nq
         WRITE(out_unit,*) '    nb     =', Basis%nb
         WRITE(out_unit,*) '    nsurf  =', nsurf
      END IF
   END IF
 
 END SUBROUTINE Set_Op
 

   SUBROUTINE Kpsi_nD(KPsi_g, Psi_g, Basis)
      USE Basis_m
      USE  QDUtil_m
      USE Molec_m
      TYPE(Basis_t), INTENT(IN), TARGET               :: Basis
      COMPLEX(KIND=Rkind), INTENT(IN), TARGET         :: psi_g(:)
      COMPLEX(KIND=Rkind), INTENT(INOUT), TARGET      :: Kpsi_g(:)
      COMPLEX(KIND=Rkind), POINTER                    :: psi_ggb(:, :, :)
      COMPLEX(KIND=Rkind), POINTER                    :: d2gg(:, :)
      COMPLEX(KIND=Rkind), POINTER                    :: Kpsi_ggb(:, :, :)
      REAL(KIND=Rkind), ALLOCATABLE                   :: GGdef(:,:)
      LOGICAL, PARAMETER                              :: debug = .TRUE.
      INTEGER                                         :: iq, i1, i3, inb, ndim
      INTEGER                                         :: Iq1, Iq2, Iq3

      IF (debug) THEN
         !write(out_unit,*) 'BEGINNING Kpsi'
         flush (out_unit)
      END IF
      
      ndim = size(Basis%tab_basis)-1
      allocate (GGdef(ndim , ndim))
      CALL get_Qmodel_GGdef(GGdef)

      Kpsi_g(:) = CZERO
      DO inb = 1, ndim 
         
      IF (inb == 1) THEN
      
          Iq1 = 1
          Iq2 = Basis%tab_basis(1)%nq
          Iq3 = Product(Basis%tab_basis(2:ndim)%nq)*Basis%tab_basis(ndim + 1)%nb

      ELSE IF (inb == ndim) THEN

         Iq1  = Product(Basis%tab_basis(1:ndim - 1)%nq)
         Iq2  = Basis%tab_basis(ndim)%nq
         Iq3  = Basis%tab_basis(ndim + 1)%nb

      ELSE
          Iq1 = Product(Basis%tab_basis(1:inb - 1)%nq)
          Iq2 = Basis%tab_basis(inb)%nq
          Iq3 = Product(Basis%tab_basis(inb + 1:ndim)%nq)*Basis%tab_basis(ndim + 1)%nb
      END IF

         Kpsi_ggb(1:Iq1, 1:Iq2, 1:Iq3) => Kpsi_g
         psi_ggb(1:Iq1, 1:Iq2, 1:Iq3) => psi_g
         d2gg(1:Iq2, 1:Iq2) => Basis%tab_basis(inb)%d2gg

         DO i3 = 1, ubound(psi_ggb, dim=3)
            DO i1 = 1, ubound(psi_ggb, dim=1)
               KPsi_ggb(i1, :, i3) = KPsi_ggb(i1, :, i3) - HALF*GGdef(inb, inb)*matmul(d2gg, psi_ggb(i1, :, i3))
            END DO
         END DO

      END DO
      IF (debug) THEN
         !        write(out_unit,*) 'END KPsi_nD'
         flush (out_unit)
      END IF
   END SUBROUTINE Kpsi_nD

   SUBROUTINE Calc_Scalar_Pot(V, Basis)
      use QDUtil_m
      use Basis_m
      use omp_lib
      implicit none
   
      ! ===== Arguments =====
      type(Basis_t), intent(in), target                :: Basis
      real(kind=Rkind), allocatable, intent(inout)     :: V(:, :, :)
   
      ! ===== Local variables =====
      integer                                          :: ndim, iq, nq, nsurf
      real(kind=Rkind), allocatable                    :: Q(:,:)
   
      ! ===== Initialization =====
      ndim  = size(Basis%tab_basis) - 1        ! Number of nuclear degrees of freedom
      nq    = Basis%nq                         ! Total number of grid points
      nsurf = Basis%tab_basis(ndim + 1)%nb     ! Number of electronic surfaces
   
      IF (.NOT. allocated(V)) allocate(V(nq, nsurf, nsurf))
      CALL Calc_Q_grid(Q, Basis)
   
      ! ===== Parallel loop avec variables privées =====
      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP          SHARED(V, Q, nq, nsurf, ndim) &
      !$OMP          PRIVATE(iq)
      
      !$OMP DO SCHEDULE(DYNAMIC)
      Do iq = 1, nq
         ! Appel de la version sans tableaux temporaires
         call calc_potential_point(V, Q, iq, nsurf, ndim)
      end do
      !$OMP END DO
   
      !$OMP END PARALLEL
   
      ! ===== Deallocation =====
      deallocate(Q)
   
   END SUBROUTINE Calc_Scalar_Pot
   
   ! ===== Subroutine auxiliaire pour éviter les tableaux temporaires =====
   SUBROUTINE calc_potential_point(V, Q, iq, nsurf, ndim)
      use QDUtil_m
      implicit none
      
      ! Arguments
      real(kind=Rkind), intent(inout)                  :: V(:, :, :)
      real(kind=Rkind), intent(in)                     :: Q(:, :)
      integer, intent(in)                              :: iq, nsurf, ndim
      
      ! Variables locales contiguës
      real(kind=Rkind)                                 :: V_local(nsurf, nsurf)
      real(kind=Rkind)                                 :: Q_local(ndim)
      integer                                          :: i, j
      
      ! Copie explicite pour garantir la contiguïté
      do i = 1, ndim
         Q_local(i) = Q(iq, i)
      end do
      
      ! Initialisation pour respecter -finit-real=nan
      V_local = ZERO
      
      ! Appel avec tableaux continus
      call sub_Qmodel_V(V_local, Q_local)
      
      ! Copie du résultat avec boucles explicites
      do j = 1, nsurf
         do i = 1, nsurf
            V(iq, i, j) = V_local(i, j)
         end do
      end do
      
   END SUBROUTINE calc_potential_point


SUBROUTINE test_op(Basis,psi0)
   USE  QDUtil_m
   USE Basis_m
   TYPE(Basis_t), intent(in), target             :: Basis
   TYPE(psi_t),intent(in)                        :: psi0
    TYPE(Op_t)                                   :: H
    TYPE(psi_t)                                  :: psi
    real(kind=Rkind), allocatable                :: V(:, :, :)
    complex(kind=Rkind) ,allocatable             :: CEigVal(:),CEigVec(:,:)
    real(kind=Rkind) ,allocatable                :: prob(:),vec(:)
    real(Kind=Rkind)            , allocatable    ::Pop(:)
    integer                                      :: nb,ndim,ib,nsurf,jb

    open(unit=50, file = 'proj.txt')
    open(unit=51, file = 'EingVec.txt')
    open(unit=52, file = 'population.txt')

     call  Make_Mat_H(Basis,H)
     ndim = size(Basis%tab_basis) - 1
     nb =Basis%nb*Basis%tab_basis(ndim+1)%nb
     nsurf =Basis%tab_basis(ndim+1)%nb
     allocate(CEigVal(nb),CEigVec(nb,nb),prob(nb),vec(nb),Pop(nsurf))
     call init_psi(psi, psi0%Basis, cplx=.true., grid=.false.)
     !call   Write_VecMat(H%CMat,out_unit,5, info='<psi|H|psi>',Rformat='e13.4')
     call diagonalization(H%CMat,CEigVal,CEigVec)
     !prob = matmul(conjg(transpose(CEigVec)),psi0%CVec)
     !prob = matmul(conjg(CEigVec),psi0%CVec)
     !print*,'E=', dot_product(psi0%CVec,matmul(H%CMat,psi0%CVec))

     !Do ib = 1,nb
     ! psi%CVec(:)=CEigVec(:,ib)
     ! prob(ib)=ZERO 
     ! prob(ib)= dot_product(conjg(psi%CVec), psi0%CVec)
     ! !Do jb = 1,nb
     !   !prob(ib) = prob(ib) + conjg(psi%CVec(jb))*psi0%CVec(jb)
     !   !print*,jb, psi%CVec(jb), psi0%CVec(jb)
     ! !End do
     ! write(50,*),ib,prob(ib)
     !End Do


     Do ib = 1,nb
      psi%CVec(:)=CEigVec(:,ib)
      prob(ib) = sum(CEigVec(:,ib)*psi0%CVec(:))
      call Popu(psi,Pop)
      write(51,*) ib, CEigVal(ib)%re
      write(50,*) ib, prob(ib)
      write(52,*) ib, Pop(:)
     End Do

     !call Calc_Scalar_Pot(V, Basis)
     

     call dealloc_psi(psi)
END SUBROUTINE 


subroutine Popu(Psi, Pop)
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
   call Calc_Norm_OF_Psi(Psi, Norm)
   do inb = 1, nsurf
      Pop(inb) = dot_product(Psi_bb(:, inb), Psi_bb(:, inb))/Norm
   end do
end subroutine 


 SUBROUTINE Calc_Hpsi(psi_g, HPsi_g, Basis,V)
    USE Basis_m
    USE Psi_m
    USE Molec_m
    complex(kind=Rkind), intent(in), target            :: psi_g(:)
    complex(kind=Rkind), intent(inout)                 :: HPsi_g(:)
    type(Basis_t), intent(in), target                  :: Basis
    real(kind=Rkind) , intent(in)                      :: V(:, :, :)

    complex(kind=Rkind), allocatable, target           :: VPsi_g(:), KPsi_g(:)
    complex(kind=Rkind), pointer                       :: VPsi_gb(:, :), Psi_gb(:, :)
    Integer                                            :: i, j, ndim,nsurf,nq


    IF (.not. allocated(Basis%tab_basis)) THEN
       STOP 'ERROR in Set_Op: the Basis%tab_bais is not initialized'
    END IF

    ndim = SIZE(Basis%tab_basis) - 1
    nq =Basis%nq
    nsurf=Basis%tab_basis(ndim+1)%nb

    ! action potential V|Psi_g> 
    allocate(VPsi_g(nq*nsurf))
    allocate(KPsi_g(nq*nsurf))

    VPsi_g(:) = CZERO
    KPsi_g(:) = CZERO
    HPsi_g(:) = CZERO

   VPsi_gb(1:nq,1:nsurf) => VPsi_g
   Psi_gb(1:nq, 1:nsurf) => Psi_g

   DO i = 1, nsurf

      DO j = 1, nsurf

         VPsi_gb(:, i) = VPsi_gb(:, i) + V(:, i, j)*Psi_gb(:, j)

      END DO

   END DO

    ! action potential K|Psi_g>

     call Kpsi_nD(KPsi_g, Psi_g, Basis)
     HPsi_g(:) = VPsi_g(:) + KPsi_g(:)
     !HPsi_g(:) = VPsi_g(:)

    DEALLOCATE (VPsi_g)
    DEALLOCATE (KPsi_g)

 END SUBROUTINE 



 SUBROUTINE calc_OpPsi_temp(Op, psi, Oppsi)

   USE psi_m
   TYPE(Op_t), intent(in)           :: Op
   TYPE(psi_t), intent(in)          :: psi
   TYPE(psi_t), intent(inout)       :: Oppsi
   TYPE(psi_t)                      :: psi_g, Oppsi_g

   IF (allocated(Op%RMat)) THEN

      Oppsi%CVec = matmul(Op%RMat, psi%CVec)

   ELSE

      IF (psi%Grid) THEN

         call init_psi(Oppsi_g, Psi%Basis, .true., .true.)
         call Calc_Hpsi(psi%CVec, Oppsi_g%CVec, psi%Basis,Op%Scal_pot)
         call GridTOBasis_nD_cplx(Oppsi%CVec, Oppsi_g%CVec, psi%Basis)
         call dealloc_psi(Oppsi_g)

      ELSE

         call init_psi(psi_g, Psi%Basis, .true., .true.)
         call init_psi(Oppsi_g, psi%Basis, .true., .true.)
         call BasisTOGrid_nD_cplx(psi_g%CVec, psi%CVec, psi%Basis)
         call Calc_Hpsi(psi_g%CVec, Oppsi_g%CVec, psi%Basis,Op%Scal_pot)
         call GridTOBasis_nD_cplx(Oppsi%CVec, Oppsi_g%CVec, psi%Basis)
         call dealloc_psi(psi_g)
         call dealloc_psi(Oppsi_g)

      END IF

   END IF

END SUBROUTINE 


SUBROUTINE calc_OpPsi(Op, psi, Oppsi)
   !---------------------------------------------------------------
   !> Applies the operator Op to the wavefunction psi to compute
   !> Oppsi = Op * psi. This works for both grid and basis representations.
   !---------------------------------------------------------------
 
   USE psi_m
 
   TYPE(Op_t),   INTENT(IN)     :: Op       ! Operator to apply
   TYPE(psi_t),  INTENT(IN)     :: psi      ! Input wavefunction
   TYPE(psi_t),  INTENT(INOUT)  :: Oppsi    ! Resulting wavefunction after application
   TYPE(psi_t)                  :: psi_g    ! Temporary wavefunction on grid
   TYPE(psi_t)                  :: Oppsi_g  ! Temporary result on grid
 
   IF (ALLOCATED(Op%RMat)) THEN
      ! If operator is represented as a matrix, apply it directly
      Oppsi%CVec = MATMUL(Op%RMat, psi%CVec)
 
   ELSE
      IF (psi%Grid) THEN
         ! If psi is on the grid: apply operator directly in grid representation
         CALL init_psi(Oppsi_g, psi%Basis, cplx = .true., grid = .true.)
         CALL Calc_Hpsi(psi%CVec, Oppsi_g%CVec, psi%Basis, Op%Scal_pot)
         CALL GridTOBasis_nD_cplx(Oppsi%CVec, Oppsi_g%CVec, psi%Basis)
         CALL dealloc_psi(Oppsi_g)
 
      ELSE
         ! If psi is in the basis: convert to grid, apply operator, and convert back
         CALL init_psi(psi_g,   psi%Basis, cplx = .true., grid = .true.)
         CALL init_psi(Oppsi_g, psi%Basis, cplx = .true., grid = .true.)
 
         CALL BasisTOGrid_nD_cplx(psi_g%CVec,   psi%CVec,   psi%Basis)
         CALL Calc_Hpsi(psi_g%CVec, Oppsi_g%CVec, psi%Basis, Op%Scal_pot)
         CALL GridTOBasis_nD_cplx(Oppsi%CVec, Oppsi_g%CVec, psi%Basis)
 
         CALL dealloc_psi(psi_g)
         CALL dealloc_psi(Oppsi_g)
      END IF
   END IF
 
 END SUBROUTINE calc_OpPsi


SUBROUTINE Calc_tab_Iq(Tab_Iq,Basis)
   USE  QDUtil_m
   USE Basis_m
   TYPE(Basis_t), intent(in), target             :: Basis
   integer, allocatable ,intent(inout)           :: Tab_iq(:, :)
   integer, allocatable                          :: Tab_iq0(:)
   integer                                       :: ndim, iq,nq
   logical                                       :: Endloop

   ndim = size(Basis%tab_basis) - 1
   nq =Basis%nq

    allocate (Tab_iq(ndim,nq),Tab_iq0(ndim))
   Call Init_tab_ind(Tab_iq0, Basis%NDindexq)
   Iq = 0
   DO
      Iq = Iq + 1
      CALL increase_NDindex(Tab_Iq0, Basis%NDindexq, Endloop)
      IF (Endloop) exit
      Tab_iq(:,Iq) = Tab_Iq0
       !print*,iq,Tab_Iq(:,Iq)
   END DO
   deallocate(Tab_Iq0)
END SUBROUTINE 


SUBROUTINE Calc_Scalar_Pot_openmp(V, Basis,Tab_Iq)
   !$ USE omp_lib
   USE  QDUtil_m
   USE Basis_m
   TYPE(Basis_t), intent(in), target               :: Basis
   real(kind=Rkind), allocatable ,intent(inout)    :: V(:, :, :)
   integer ,intent(in)                             :: Tab_iq(:, :)

   integer                                          :: Ib, ndim, iq,nq,nsurf,maxth
   real(kind=Rkind)               , allocatable     :: Q(:)
   logical                                          :: Endloop

   ndim = size(Basis%tab_basis) - 1
   nq   =Basis%nq
   nsurf=Basis%tab_basis(ndim+1)%nb

    allocate (Q(ndim))
    allocate(V(nq, nsurf, nsurf))
    V(:,:,:) = ZERO
    Q(:) = ZERO

    maxth = 1
    !$ maxth  = omp_get_max_threads()
    write(*,*) 'nbr de coeurs :',maxth
   !$OMP   PARALLEL DEFAULT(NONE) &
   !$OMP   SHARED(Basis,Tab_Iq,maxth,V,nq,ndim,Ib) &
   !$OMP   PRIVATE(Iq,Q) &
   !$OMP   NUM_THREADS(maxth)

   !$OMP BARRIER
   !$OMP   DO SCHEDULE(STATIC)

   DO  Iq = 1,nq

      DO Ib = 1, ndim
          Q(Ib) = Basis%tab_basis(Ib)%x(Tab_Iq(Ib,Iq))
      ENDDO

      !CALL sub_pot(V(Iq, :, :), Q(:), 0)
      CALL sub_Qmodel_V(V(iq, :, :), Q(:))

   END DO
   !$OMP END DO
   !$OMP BARRIER
   !$OMP END PARALLEL

   deallocate(Q)
END SUBROUTINE 



SUBROUTINE test_openmp_op(Basis)
   USE  QDUtil_m
   USE Basis_m
   TYPE(Basis_t), intent(in), target             :: Basis
   real(kind=Rkind), allocatable                 :: V(:, :, :)
   integer, allocatable                          :: Tab_Iq(:, :)
   integer                                        :: nq

  nq = Basis%nq
  print*,"nq=",nq
   call Calc_tab_Iq(Tab_Iq,Basis)
   call  Calc_Scalar_Pot(V, Basis)
   
     
   !call cpu_time(time=t1)
    !call Calc_Scalar_Pot_openmp(V, Basis,Tab_Iq)
   ! call cpu_time(time=t2)
    !tpsopenmp = t2-t1
    !write(*,*)"temps CPU  openmp :",tpsopenmp
    !deallocate(V)
    !call cpu_time(time=t1)
    !call Calc_Scalar_Pot(V, Basis)
    !call cpu_time(time=t2)
    !tps = t2-t1
    !write(*,*)"temps CPU  non parallele :",tps

END SUBROUTINE 





SUBROUTINE calc_nac(Basis)
   USE  QDUtil_m
   USE Basis_m
   TYPE(Basis_t), intent(in), target                :: Basis
   real(kind=Rkind), allocatable                    :: NAC(:, :, :)

   integer, allocatable                             :: Tab_iq(:)
   integer                                          :: inb, ndim, iq,nq,nsurf,i,nq1,nq2
   real(Kind=Rkind)               , allocatable     :: Q(:),G(:,:,:),V(:,:)
   logical                                          :: Endloop

   ndim = size(Basis%tab_basis) - 1
   nq =Basis%nq
   nsurf=Basis%tab_basis(ndim+1)%nb
   nq1= Basis%tab_basis(1)%nq
   nq2= Basis%tab_basis(2)%nq

   open (unit = 300, file = "N121-ad=t.txt")
   open (unit = 301, file = "N122-ad=t.txt")

   open (unit = 302, file = "N211-ad=t.txt")
   open (unit = 303, file = "N212-ad=t.txt")
   open (unit = 500, file = "x.txt")
   open (unit = 501, file = "y.txt")

    allocate (Q(ndim),Tab_iq(ndim))
    allocate(V(nsurf, nsurf), NAC(nsurf, nsurf,ndim), G(nsurf, nsurf,ndim))

   Call Init_tab_ind(Tab_iq, Basis%NDindexq)
   Iq = 0
   DO
      Iq = Iq + 1
      i =Tab_iq(2)
      CALL increase_NDindex(Tab_iq, Basis%NDindexq, Endloop)
      IF (Endloop) exit
      do inb = 1, ndim
          Q(inb) = Basis%tab_basis(inb)%x(Tab_iq(inb))
      end do
      CALL sub_Qmodel_VG_NAC(V, G, NAC, Q)
      if ( (Tab_iq(2)/= i) )then 
          write(300,*)
          write(301,*)
          write(302,*)
          write(302,*)
      end if  

      write(300,*) Q(:),NAC(1, 2, 1)
      write(301,*) Q(:),NAC(1, 2, 2)

      write(302,*) Q(:),NAC(2, 1, 1)
      write(303,*) Q(:),NAC(2, 1, 2)

   END DO
   do i = 1, nq1
      write(500,*) Basis%tab_basis(1)%x(i)
   end do
   do i = 1, nq2
      write(501,*) Basis%tab_basis(2)%x(i)
   end do

   deallocate(Tab_iq,Q, NAC,G,V)

END SUBROUTINE 




SUBROUTINE calc_VV(Basis)
   USE  QDUtil_m
   USE Basis_m
   TYPE(Basis_t), intent(in), target                :: Basis

   integer, allocatable                             :: Tab_iq(:)
   integer                                          :: inb, ndim, iq,nq,nsurf,i,nq1,nq2
   real(Kind=Rkind)               , allocatable     :: Q(:),V(:,:)
   logical                                          :: Endloop

   ndim = size(Basis%tab_basis) - 1
   nq =Basis%nq
   nsurf=Basis%tab_basis(ndim+1)%nb
   nq1= Basis%tab_basis(1)%nq
   nq2= Basis%tab_basis(2)%nq

   !open (unit = 400, file = "v11-ad=t.txt")
   !open (unit = 401, file = "v22-ad=t.txt")

   !open (unit = 402, file = "v21-ad=t.txt")
   !open (unit = 403, file = "v12-ad=t.txt")
   !open (unit = 404, file = "x.txt")
   !open (unit = 405, file = "y.txt")

    allocate (Q(ndim),Tab_iq(ndim))
    allocate(V(nsurf, nsurf))

   Call Init_tab_ind(Tab_iq, Basis%NDindexq)
   Iq = 0
   DO
      Iq = Iq + 1
      i =Tab_iq(2)
      CALL increase_NDindex(Tab_iq, Basis%NDindexq, Endloop)
      IF (Endloop) exit
      do inb = 1, ndim
          Q(inb) = Basis%tab_basis(inb)%x(Tab_iq(inb))
      end do
      CALL sub_Qmodel_V(V, Q)
      !if ( (Tab_iq(2)/= i) )then 
       !   write(400,*)
        !  write(401,*)
          !write(402,*)
         ! write(403,*)
      !end if  

      !write(400,*) Q(:), V(1, 1)
      !write(401,*) Q(:), V(2, 2)

      !write(402,*) Q(:), V(2, 1)
      !write(403,*) Q(:),V(1, 2)

   END DO
   !do i = 1, nq1
    !  write(404,*) Basis%tab_basis(1)%x(i)
   !end do
   !do i = 1, nq2
    !  write(405,*) Basis%tab_basis(2)%x(i)
   !end do

   deallocate(Tab_iq,Q,V)

END SUBROUTINE 

 
end module Op_m