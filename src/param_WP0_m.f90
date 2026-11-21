
      MODULE param_WP0_m
      USE QDUtil_m
      USE Basis_m, ONLY: Debug_Enabled, Basis_Debug_Active
      IMPLICIT NONE

        TYPE GWP1D_t

        real (kind=Rkind)                  :: sigma = ONETENTH  ! width of WP0
        real (kind=Rkind)                  :: Beta
        real (kind=Rkind)                  :: Q0    = ZERO      ! position of WP0
        real (kind=Rkind)                  :: imp_k = ZERO      ! impultion for WP0
        real (kind=Rkind)                  :: gamma = ZERO      ! phase for WP0


        END TYPE GWP1D_t

        TYPE GWP_t
          integer                     :: ndim           = 1
          complex (kind=Rkind)        :: Coef           = CZERO
          integer                     ::    Elecindex   = 1
          character(len=16)            :: source = 'input'

          TYPE (GWP1D_t), allocatable :: tab_GWP1D(:)

        END TYPE GWP_t

        CONTAINS

   RECURSIVE SUBROUTINE Read_GWP1D(GWP1D)
    TYPE (GWP1D_t), intent(inout)          :: GWP1D

    !------ initial WP definition -----------------------------
    !     GWP(Q)=exp[-((Q-Qeq)/sigma)2+i*imp_k*(Q-Qeq)+i*gamma]
    real (kind=Rkind)                         :: sigma,Beta,imp_k,Qeq,gamma
    integer                                   :: Rerr

    NAMELIST /defWP0/ sigma,Beta,Qeq,imp_k,gamma


    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Read_GWP1D'

    !logical, parameter :: debug =.TRUE.
    !-----------------------------------------------------------

    sigma    = sqrt(TWO)
    Beta     = ZERO
    Qeq      = ZERO
    imp_k    = ZERO
    gamma    = ZERO

    read(in_unit,defWP0,iostat=Rerr)
    if (Rerr /= 0) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) ' problem while reading the namelist "defWP0"'
      write(out_unit,defWP0)
      STOP 'ERROR in Read_GWP1D: problem while reading the namelist "defWP0"'
    end if
          GWP1D = GWP1D_t(sigma,Beta,Qeq,imp_k,gamma)
       

  END SUBROUTINE Read_GWP1D

         SUBROUTINE Write_GWP1D(GWP1D)
           TYPE (GWP1D_t), intent(in) :: GWP1D
        
           !----- for debuging --------------------------------------------------
           character (len=*), parameter :: name_sub='Write_GWP1D'

           !logical, parameter :: debug =.TRUE.
           !-----------------------------------------------------------
        
           IF (Basis_Debug_Active() .OR. Debug_Enabled()) write(out_unit,*) 'Q0,sigma,imp_k,phase', &
                  GWP1D%Q0,GWP1D%sigma, GWP1D%Beta,GWP1D%imp_k,GWP1D%gamma
        
           write(out_unit,*)                                             &
                'sigma=',GWP1D%sigma, 'Beta=',GWP1D%Beta,'Qeq=',GWP1D%Q0 ,'imp_k =',GWP1D%imp_k,'gamma=',GWP1D%gamma
        
         END SUBROUTINE Write_GWP1D

         FUNCTION calc_GWP1D(GWP1D,Q)
           complex (kind=Rkind)             :: calc_GWP1D
           real (kind=Rkind), intent(in)    :: Q
           TYPE (GWP1D_t),    intent(in) :: GWP1D
        
           real (kind=Rkind) :: ze,zk,DQ
           complex(kind=Rkind) :: At
           !----- for debuging --------------------------------------------------
           character (len=*), parameter :: name_sub='calc_GWP1D'

           !logical, parameter :: debug =.TRUE.
           !-----------------------------------------------------------
           At = complex((TWO/GWP1D%sigma**2),TWO*GWP1D%Beta)
           !print*,'at=',At%re
           DQ = Q-GWP1D%Q0
           ze = (DQ/GWP1D%sigma)**2
           zk = mod(-DQ*DQ*GWP1D%Beta+DQ*GWP1D%imp_k+GWP1D%gamma ,TWO*pi)
        
           calc_GWP1D = exp(-ze + EYE*zk) / sqrt(sqrt(pi/TWO)*GWP1D%sigma)
           !calc_GWP1D = exp(-(At/TWO)*DQ*DQ+ EYE*DQ +GWP1D%gamma)*sqrt(sqrt(At%re/pi))
        
         END FUNCTION calc_GWP1D
        
        SUBROUTINE Read_GWP(GWP,nio)
          TYPE (GWP_t), intent(inout) :: GWP
          integer,intent(in)          :: nio
          !------ initial WP definition -----------------------------
          integer                      :: i,Rerr
          integer                      ::  ndim,Elecindex
          character(len=16)             :: source

          complex (kind=Rkind)            :: Coef

          NAMELIST /defGWP/ ndim,Elecindex,Coef,source

          !----- for debuging --------------------------------------------------
          character (len=*), parameter :: name_sub='Read_GWP'

          !logical, parameter :: debug =.TRUE.
          !-----------------------------------------------------------


         ndim            = 1
         Coef            = CONE
         Elecindex       = 1 
         source          = 'input'

          read(nio,nml=defGWP,iostat=Rerr)
          IF (Rerr /= 0) THEN  
            print*,'Rerr',Rerr
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) ' problem while reading the namelist "defGWP"'
            write(out_unit,defGWP)
            STOP 'ERROR in Read_GWP: problem while reading the namelist "defGWP"'
          END IF
           GWP%ndim          = ndim
           GWP%Coef          = Coef
           GWP%Elecindex     = Elecindex
           GWP%source        = TRIM(ADJUSTL(source))
           CALL string_uppercase_TO_lowercase(GWP%source)

          allocate(GWP%tab_GWP1D(ndim))
          SELECT CASE (TRIM(GWP%source))
          CASE ('input')
            DO i=1,ndim
              CALL Read_GWP1D(GWP%tab_GWP1D(i))
            END DO
          CASE ('qmodel')
            ! Consume the usual per-dimension records; the caller replaces Q0 and sigma
            ! with values returned by InitQModelBasisParam. Empty &defWP0 / records suffice.
            DO i=1,ndim
              CALL Read_GWP1D(GWP%tab_GWP1D(i))
            END DO
          CASE DEFAULT
            ERROR STOP "Read_GWP: source must be 'input' or 'qmodel'."
          END SELECT

        END SUBROUTINE Read_GWP

       FUNCTION calc_GWP(GWP,Q)
         complex (kind=Rkind)          :: calc_GWP
         real (kind=Rkind), intent(in) :: Q(:)
         TYPE (GWP_t),      intent(in) :: GWP

         integer :: i
         !----- for debuging --------------------------------------------------
         character (len=*), parameter :: name_sub='calc_GWP'

         !logical, parameter :: debug =.TRUE.
         !-----------------------------------------------------------

         calc_GWP = CONE
         DO i=1,size(GWP%tab_GWP1D)
           calc_GWP = calc_GWP * calc_GWP1D(GWP%tab_GWP1D(i),Q(i))
         END DO
         calc_GWP = calc_GWP * GWP%Coef

       END FUNCTION calc_GWP


      SUBROUTINE Write_GWP(GWP)
        TYPE (GWP_t), intent(in) :: GWP
      
        integer :: i
        !----- for debuging --------------------------------------------------
        character (len=*), parameter :: name_sub='Write_GWP'

        !logical, parameter :: debug =.TRUE.
        !-----------------------------------------------------------
      
        write(out_unit,*) 'ndim,Elecindex,Coef,source',GWP%ndim , GWP%Elecindex , GWP%Coef, TRIM(GWP%source)
        write(out_unit,'(a)') ' Qeq    Beta     sigma       imp_k       phase'
        DO i=1,size(GWP%tab_GWP1D)
          CALL Write_GWP1D(GWP%tab_GWP1D(i))
        END DO
      
      END SUBROUTINE Write_GWP

         SUBROUTINE Read_tab_GWP(tab_GWP,nb_GWP,nio)
        
           TYPE (GWP_t), allocatable, intent(inout) :: Tab_GWP(:)
           integer,                   intent(in)    :: nb_GWP
           integer,                   intent(in)    :: nio
        
        
        
           !------ initial WP definition -----------------------------
             integer                :: i
        
           !----- for debuging --------------------------------------------------
           character (len=*), parameter :: name_sub='Read_tab_GWP'

           !logical, parameter :: debug =.TRUE.
           !-----------------------------------------------------------
        
           allocate(Tab_GWP(nb_GWP))
           DO i=1,nb_GWP
             CALL Read_GWP(Tab_GWP(i),nio)
           END DO        
           IF (Basis_Debug_Active() .OR. Debug_Enabled()) CALL Write_Tab_GWP(Tab_GWP)
         END SUBROUTINE Read_tab_GWP


         SUBROUTINE Write_Tab_GWP(Tab_GWP)
           TYPE (GWP_t), intent(in) :: Tab_GWP(:)
        
           integer :: i
           !----- for debuging --------------------------------------------------
           character (len=*), parameter :: name_sub='Write_Tab_GWP'

           !logical, parameter :: debug =.TRUE.
           !-----------------------------------------------------------
        
           DO i=1,size(Tab_GWP)
             CALL Write_GWP(Tab_GWP(i))
           END DO
        
         END SUBROUTINE Write_Tab_GWP

         SUBROUTINE calc_Tab_GWP(Tab_GWP,Q,psi)
          complex (kind=Rkind) ,intent(inout)   :: psi(:)
          real (kind=Rkind), intent(in)         :: Q(:)
          TYPE (GWP_t),      intent(in)      :: Tab_GWP(:)
          integer :: iGWP0
          !----- for debuging --------------------------------------------------
          character (len=*), parameter :: name_sub='calc_GWP'

          !logical, parameter :: debug =.TRUE.
          !-----------------------------------------------------------
          psi(:) = CZERO     

          DO iGWP0 = 1,size(Tab_GWP),1
               psi(Tab_GWP(iGWP0)%Elecindex)  = calc_GWP(Tab_GWP(iGWP0),Q)
          END DO
          
        
        END SUBROUTINE calc_Tab_GWP


        


      END MODULE param_WP0_m
