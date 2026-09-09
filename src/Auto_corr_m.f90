module Auto_corr_m
   USE QDUtil_m
   USE psi_m
   USE Basis_m
   Use Hagedorn_m
   Use sub_propa_m

   implicit none
   public :: Calc_Auto_corr,test_write_cplx,eval_pics

contains

   !==============================================================
   !> @brief
   !>   Compute the autocorrelation function between the reference
   !>   wavefunction (psi0) and the propagated wavefunction (psi_dt).
   !>
   !> @details
   !>   - If the propagation scheme is 'hagedorn', the routine applies
   !>     the inverse Hagedorn transformation before computing the 
   !>     autocorrelation coefficient.
   !>   - Otherwise, the autocorrelation is computed directly.
   !>   - The routine also evaluates the phase (argument) of the
   !>     autocorrelation coefficient.
   !>   - Optional debug output is provided if debug=.true.
   !>
   !> @param[in]     psi0             Reference wavefunction at t=0
   !> @param[in]     psi_dt           Propagated wavefunction at time t
   !> @param[inout]  corre_coeff      Autocorrelation coefficient
   !> @param[inout]  arg_corre_coeff  Phase of autocorrelation coefficient
   !> @param[in]     propa_name       Propagation scheme name ('hagedorn' or else)
   !> @param[in]     renorm           Logical flag for renormalization
   !> @param[in]     t                (Optional) Time value
   !> @param[in]     it               (Optional) Iteration index
   !> @param[in]     debug            (Optional) If .true., write debug information
   !==============================================================
   SUBROUTINE Calc_Auto_corr(psi0, psi_dt, corre_coeff, arg_corre_coeff,propa, t, it, debug)
   USE QDUtil_m
   USE psi_m
   
   !=====================
   ! Arguments
   !=====================
   TYPE(psi_t),      intent(in)            :: psi0, psi_dt
   complex(Rkind),   intent(inout)         :: corre_coeff
   real(Rkind),      intent(inout)         :: arg_corre_coeff
   TYPE(propa_t),    intent(in)            :: propa
   real(Rkind),      intent(in), optional  :: t
   integer,          intent(in), optional  :: it
   logical,          intent(in), optional  :: debug
   
   !=====================
   ! Local variables
   !=====================
   TYPE(psi_t)       :: psi, psi_t0
   TYPE(Basis_t)     :: Basis0, Basis_dt
   real(Rkind)       :: X, Y
   logical           :: ldebug
   
   !=====================
   ! Initialization
   !=====================
   ldebug = .false.
   IF (PRESENT(debug)) ldebug = debug
   
   !-----------------------------------------------------------
   ! Begin autocorrelation calculation
   !-----------------------------------------------------------
   IF (ldebug) WRITE(out_unit,*) ">>> Entering Calc_Auto_corr"
   
   IF (propa%propa_name == 'hagedorn') THEN
   !--- Initialize Basis0 from psi0 and construct primitive basis
   CALL init_Basis1_TO_Basis2(Basis0, psi0%Basis)
   CALL construct_primitive_basis(Basis0)
   
   !--- Initialize Basis_dt from psi_dt and construct primitive basis
   CALL init_Basis1_TO_Basis2(Basis_dt, psi_dt%Basis)
   CALL construct_primitive_basis(Basis_dt)
   
   !--- Initialize psi and psi_t0 (complex vectors, grid = false)
   CALL init_psi(psi,    Basis0, cplx=.TRUE., grid=.false.)
   CALL init_psi(psi_t0, Basis_dt, cplx=.TRUE., grid=.false.)
   
   !--- Set psi to zero and copy psi_dt to psi_t0
   psi%CVec    = CZERO
   psi_t0%CVec = psi_dt%CVec
   
   !--- Apply inverse Hagedorn transformation
   CALL Hagedorn(psi, psi_t0, propa, use_psi0_params=.TRUE.)
   
   !--- Compute autocorrelation coefficient
   corre_coeff     = Dot_Psi(psi0,psi)
   X               = REAL(corre_coeff, kind=Rkind)
   Y               = AIMAG(corre_coeff)
   arg_corre_coeff = ATAN2(Y, X)
   
   !--- Debug output
   IF (ldebug) THEN
   WRITE(out_unit,*) "Hagedorn autocorrelation:"
   WRITE(out_unit,*) "   corre_coeff     =", corre_coeff
   WRITE(out_unit,*) "   arg_corre_coeff =", arg_corre_coeff
   END IF
   
   !--- Optional wavefunction output
   IF (ldebug) THEN
   IF (PRESENT(t)) THEN
   CALL test_write_cplx(psi, ib=27, t=t)
   WRITE(27,*)
   END IF
   
   IF (PRESENT(it)) CALL test_write(psi, ib=27)    
   IF (PRESENT(it)) WRITE(27,*)
   END IF
   
   !--- Deallocate
   CALL dealloc_psi(psi)
   CALL dealloc_psi(psi_t0)
   
   ELSE
   !-----------------------------------------------------------
   ! Direct autocorrelation calculation (non-Hagedorn case)
   !-----------------------------------------------------------
   corre_coeff     = Dot_Psi(psi0,psi_dt)
   X               = REAL(corre_coeff, kind=Rkind)
   Y               = AIMAG(corre_coeff)
   arg_corre_coeff = ATAN2(Y, X)
   
   !--- Debug output
   IF (ldebug) THEN
   WRITE(out_unit,*) "Direct autocorrelation:"
   WRITE(out_unit,*) "   corre_coeff     =", corre_coeff
   WRITE(out_unit,*) "   arg_corre_coeff =", arg_corre_coeff
   END IF
   
   !--- Optional wavefunction output
   IF (ldebug) THEN
   IF (PRESENT(t)) THEN
   CALL test_write_cplx(psi_dt, ib=27, t=t)
   WRITE(27,*)
   END IF
   
   IF (PRESENT(it)) CALL test_write(psi_dt, ib=27)
   IF (PRESENT(it)) WRITE(27,*)
   END IF
   END IF
   
   IF (ldebug) WRITE(out_unit,*) "<<< Leaving Calc_Auto_corr"
   !-----------------------------------------------------------
   
   END SUBROUTINE Calc_Auto_corr


   SUBROUTINE test_write(psi,ib)

   implicit none

   TYPE(psi_t), intent(in)                    :: psi
   integer,intent(in)                          :: ib

  !do i =1,size(psi%CVec)

  write(ib,*) psi%CVec(:)
 ! end do

   END SUBROUTINE


   SUBROUTINE test_write_cplx(psi,ib,t)
      implicit none
      TYPE(psi_t), intent(in)                     :: psi
      integer,intent(in)                          :: ib
      real(kind=Rkind) ,intent(in)                :: t
     integer                                      :: i

        do i =1,size(psi%CVec)
         write(ib,*) t, i, psi%CVec(i)
        end do

   END SUBROUTINE


   SUBROUTINE eval_pics(psi,ib,t)
      implicit none
      TYPE(psi_t), intent(in)                     :: psi
      integer,intent(in)                          :: ib
      real(kind=Rkind) ,intent(in)                :: t
      real(kind=Rkind)                            :: norm,n0
      integer                                     :: i
       
      norm = 0._Rkind
      call Calc_Norm_OF_Psi(psi,n0)
       do i =2,size(psi%CVec)
         norm = norm + abs(psi%CVec(i))**2
       end do

       if ( t==0 ) then
         write(ib,*) t, abs(n0**2-abs(psi%CVec(1))**2) , ZERO  !FMT= "(F20.10,F20.15,F20.18)"
       else
         write(ib,*) t, abs(n0**2-abs(psi%CVec(1))**2) , norm
       end if

   END SUBROUTINE


   !SUBROUTINE eval_pics_temp(psi,ib,t)
   !   implicit none
   !   TYPE(psi_t), intent(in)                     :: psi
   !   integer,intent(in)                          :: ib
   !   real(kind=Rkind) ,intent(in)                :: t
   !   integer                                     :: nb
   !    
   !     nb = size(psi%CVec)
   !     if ( t==0 ) then
   !      write(ib,*) t, abs(ONE-abs(psi%CVec(1))) , 0.0
   !     else
   !      write(ib,*) t, abs(ONE-abs(psi%CVec(1))**2) ,sqrt( sum(abs(psi%CVec(2:nb))**2))
   !     end if
   !      
   !END SUBROUTINE

end module Auto_corr_m

