MODULE HamiltonianSpectrum_m
   USE QDUtil_m
   USE Basis_m
   USE psi_m
  USE Op_m, ONLY: Op_t
  USE TISE_m, ONLY: Build_Hamiltonian_Matrix
   USE Ana_psi_m, ONLY: Popu
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: H_spectrum, test_op, Build_EigVecs_On_Grid, save_EigVecs
CONTAINS


SUBROUTINE test_op(H,psi0)
   USE  QDUtil_m
   USE Basis_m
   TYPE(Op_t), intent(inout)                  :: H
   TYPE(psi_t),intent(in)                        :: psi0
    TYPE(psi_t)                                  :: psi
    real(kind=Rkind), allocatable                :: V(:, :, :)
    complex(kind=Rkind) ,allocatable             :: CEigVal(:),CEigVec(:,:)
    real(kind=Rkind) ,allocatable                :: prob(:),vec(:)
    real(Kind=Rkind)            , allocatable    ::Pop(:)
    integer                                      :: nb,ndim,ib,nsurf,jb

    open(unit=50, file = 'proj.txt')
    open(unit=51, file = 'EingVec.txt')
    open(unit=52, file = 'population.txt')

     call  Build_Hamiltonian_Matrix(H)
     ndim = size(H%Basis%tab_basis) - 1
     nb =H%Basis%nb*H%Basis%tab_basis(ndim+1)%nb
     nsurf =H%Basis%tab_basis(ndim+1)%nb
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

SUBROUTINE H_spectrum(H, N_states,psi)
  USE QDUtil_m
  USE Basis_m
  IMPLICIT NONE

  TYPE(Op_t), INTENT(INOUT)            :: H
  TYPE(psi_t) , INTENT(INOUT),OPTIONAL :: psi
  INTEGER,       INTENT(IN)            :: N_states

  COMPLEX(kind=Rkind), ALLOCATABLE     :: CEigVal(:), CEigVec(:,:)
  REAL(kind=Rkind), ALLOCATABLE        :: Q(:,:)
  COMPLEX(kind=Rkind), ALLOCATABLE     :: Psi_Grid(:,:)  ! (nq, N_states)
  INTEGER                              :: nb, ndim, nsurf, nq, ngrid
  INTEGER                              :: ist
  INTEGER                              :: irow, icol, unit_hmat, unit_evec

  INTEGER, PARAMETER :: Energie_File = 51
  CHARACTER(len=*), PARAMETER :: Energie_Filename = 'EigVal.txt'
  CHARACTER(len=*), PARAMETER :: Etats_Filename   = 'EigVec.txt'

  ! — Construire H —
  CALL Build_Hamiltonian_Matrix(H)

  ndim  = SIZE(H%Basis%tab_basis) - 1
  nsurf = H%Basis%tab_basis(ndim+1)%nb
  nb    = Get_TotalBasisCoefficientCount(H%Basis)
  nq    = H%Basis%nq
  ngrid = nq * nsurf

  IF (N_states > nb) THEN
    WRITE(out_unit,*) "Erreur : N_states > dimension H"
    STOP
  END IF

  ALLOCATE(CEigVal(nb), CEigVec(nb, nb))
  CALL Calc_Q_grid(Q, H%Basis)

  WRITE(out_unit,*) 'Hamiltonian matrix H:'
  CALL Write_VecMat(H%CMat, out_unit, MIN(nb,10), info='H')

  ! Diagonalisation
  CALL diagonalization(H%CMat, CEigVal, CEigVec,sort=1)
  WRITE(out_unit,*) "Diagonalisation terminée."
  IF(present(psi)) psi%CVec= CEigVec(:,1)

  ! Machine-readable exports used by the H2O comparison campaign.
  OPEN(NEWUNIT=unit_hmat, FILE='H_matrix.txt', STATUS='REPLACE', ACTION='WRITE')
  WRITE(unit_hmat,'(A)') '# i j real(Hij) imag(Hij)'
  DO irow = 1, nb
    DO icol = 1, nb
      WRITE(unit_hmat,*) irow, icol, &
           REAL(H%CMat(irow,icol),KIND=Rkind), AIMAG(H%CMat(irow,icol))
    END DO
  END DO
  CLOSE(unit_hmat)
  OPEN(NEWUNIT=unit_evec, FILE='EigVec_basis.txt', STATUS='REPLACE', ACTION='WRITE')
  WRITE(unit_evec,'(A)') '# i state real(v_i) imag(v_i)'
  DO irow = 1, nb
    DO ist = 1, N_states
      WRITE(unit_evec,*) irow, ist, &
           REAL(CEigVec(irow,ist),KIND=Rkind), AIMAG(CEigVec(irow,ist))
    END DO
  END DO
  CLOSE(unit_evec)

  ! Construire les fonctions d’onde sur la grille
  ALLOCATE(Psi_Grid(ngrid, N_states))
  CALL Build_EigVecs_On_Grid(H%Basis, N_states, CEigVec, Psi_Grid)
  
  ! Sauvegarder Q et la partie réelle des états propres
  CALL save_EigVecs(Q, Psi_Grid, nq, N_states, Etats_Filename)

  ! Sauvegarder les énergies
  OPEN(unit = Energie_File, file = Energie_Filename, status = 'replace')
  DO ist = 1, N_states
    WRITE(Energie_File, *) ist, REAL(CEigVal(ist), kind=Rkind)
  END DO
  CLOSE(Energie_File)

  CALL Write_VecMat(REAL(CEigVal), out_unit, nb, Rformat='F15.10', info='CEigVal')
  WRITE(out_unit,*) 'Hamiltonian eigenvectors (first states):'
  CALL Write_VecMat(CEigVec(:,1:N_states), out_unit, MIN(nb,10), info='CEigVec')

  ! Libération mémoire
  DEALLOCATE(Q, CEigVal, CEigVec, Psi_Grid)
  !CALL dealloc_Op(H)

END SUBROUTINE H_spectrum


subroutine save_EigVecs(Q, EigVec, nq, n_states, filename)
  implicit none

  real(kind=Rkind),    intent(in) :: Q(:,:)             ! (nq, dim)
  complex(kind=Rkind), intent(in) :: EigVec(:,:)        ! (nq, n_states)
  integer,     intent(in) :: nq, n_states
  character(*),intent(in) :: filename

  integer :: i, ist, unit

  open(newunit=unit, file=filename, status='replace')

  do i = 1, nq
     write(unit,*) Q(i,:), ( real(EigVec(i,ist),kind=Rkind), ist=1,n_states )
  end do

  close(unit)

end subroutine 




SUBROUTINE Build_EigVecs_On_Grid(Basis, N_states, CEigVec, Psi_Grid)
  USE Basis_m
  IMPLICIT NONE

  TYPE(Basis_t), INTENT(IN), TARGET :: Basis
  INTEGER, INTENT(IN)              :: N_states
  COMPLEX(kind=Rkind), INTENT(IN)  :: CEigVec(:,:)     ! vecteurs propres en base (nb x nb)
  COMPLEX(kind=Rkind), INTENT(OUT) :: Psi_Grid(:,:)    ! vecteurs propres sur la grille (nq x N_states)

  TYPE(psi_t) :: psi_b, psi_g
  INTEGER     :: k

  DO k = 1, N_states

     CALL init_psi(psi_b, Basis, cplx=.true., grid=.false.)
     CALL init_psi(psi_g, Basis, cplx=.true., grid=.true.)

     ! Charger le vecteur propre n°k dans la base
     psi_b%CVec = CEigVec(:, k)

     ! Passage Base → Grille
     CALL BasisTOGrid_nD_cplx(psi_g%CVec, psi_b%CVec, Basis)

     ! Copier dans la sortie
     Psi_Grid(:, k) = psi_g%CVec(:)

     CALL dealloc_psi(psi_b)
     CALL dealloc_psi(psi_g)

  END DO

END SUBROUTINE

END MODULE HamiltonianSpectrum_m
