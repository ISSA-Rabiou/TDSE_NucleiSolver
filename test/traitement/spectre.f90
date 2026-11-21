!===========================================================================
!===========================================================================
!This file is part of ElVibRot.
!
! MIT License
!
! Copyright (c) 2024 David Lauvergnat
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!
!===========================================================================
!===========================================================================
PROGRAM spectre
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : in_unitp=>INPUT_UNIT, out_unitp=>OUTPUT_UNIT, Rk=>real64
  implicit none
  
  real (kind=Rk),    parameter :: pi = 3.14159265358979323846264338327950288419716939937511_Rk
  real (kind=Rk),    parameter :: ZERO = 0._Rk
  real (kind=Rk),    parameter :: ONE  = 1._Rk
  real (kind=Rk),    parameter :: TWO  = 2._Rk
  complex (kind=Rk), parameter :: EYE = (ZERO,ONE)

  integer :: npts,nptE
  real (kind=Rk) :: t0,tmax,dt
  real (kind=Rk) :: E,Emin,Emax,dE
  real (kind=Rk) :: conv_E,conv_t
  real (kind=Rk),    allocatable :: t(:)
  complex (kind=Rk), allocatable :: auto0(:),auto(:),auto1(:),auto2(:)
  complex (kind=Rk), allocatable :: Expiwt(:)
  complex (kind=Rk) :: funcE0,funcE,funcE1,funcE2
  character (len=50) :: file_auto

  real (kind=Rk) :: a,b,tdum
  integer :: nio,ioerr,i,n_filter,option

  character (len=50) :: name_dum


  namelist / param / dE,Emin,Emax,conv_E,conv_t,file_auto,n_filter,option


  file_auto = 'file_auto'
  conv_E   = ONE
  conv_t   = ONE
  dE       = -ONE
  Emin     = ZERO
  Emax     = -ONE
  n_filter = 1
  option   = 0
  read(in_unitp,param)
  write(out_unitp,param)

  ! read the time function (autocorrelation...)

  open(newunit=nio,file=file_auto)
  IF (option == 1) read(nio,*,IOSTAT=ioerr) name_dum
  i = 0
  DO
    read(nio,*,IOSTAT=ioerr) name_dum
    IF (ioerr /= 0) EXIT
    i = i + 1
  END DO
  close(nio)
  npts = i
  write(out_unitp,*) 'npts',npts
  IF (npts == 0) THEN
    write(out_unitp,*) ' ERROR with file_auto'
    STOP
  END IF

  allocate(auto0(npts))
  allocate(auto(npts))
  allocate(auto1(npts))
  allocate(auto2(npts))
  auto(:) = ZERO
  allocate(t(npts))

  open(newunit=nio,file=file_auto)
  SELECT CASE(option)
  CASE (1)
    read(nio,*,IOSTAT=ioerr) name_dum
    DO i=1,npts
      read(nio,*) t(i),a,b
      auto0(i) = cmplx(a,b,kind=Rk)
   END DO
  CASE (2)
    DO i=1,npts
      read(nio,*) t(i),auto0(i)
    END DO
  CASE DEFAULT
    DO i=1,npts
      read(nio,*) name_dum,t(i),a,b
      auto0(i) = cmplx(a,b,kind=Rk)
    END DO
  END SELECT
  auto0(1) = auto0(1) / TWO
  close(nio)

  t(:) = t(:) * conv_t
  t0   = t(1)
  tmax = t(npts)
  dt   = t(2)-t(1)

  write(out_unitp,*) 'npts t0,tmax,dt: ',npts,t0,tmax,dt
  ! END read the time function

!     Check the funcErgy grid...
  IF (Emax <0) THEN
    Emax=TWO*Pi/dt
  ELSE
    Emax=Emax/conv_E
  END IF
  IF (Emax > TWO*Pi/dt) THEN
    write(out_unitp,*) 'Emax > 2.d0*Pi/dt',Emax,TWO*Pi/tmax
    STOP
  END IF
  Emax = Emax*conv_E

  IF (dE <= ZERO) dE = TWO*Pi/(tmax*TWO**2)*conv_E

  nptE = int((Emax-Emin)/dE)

  auto1(:) = auto0 * cos(pi/(TWO*tmax)*t)**1
  auto2(:) = auto0 * cos(pi/(TWO*tmax)*t)**2
  auto(:)  = auto0 * exp(-t/tmax)**n_filter

  allocate(Expiwt(npts))
  open(newunit=nio,file='file_spectrum.txt')

  DO i=1,nptE
    E = real(i-1,kind=Rk)*dE
    Expiwt(:) = exp(EYE*E/conv_E*t)

    funcE0 = TWO*dt*sum(Expiwt(:)*auto0(:))
    funcE  = TWO*dt*sum(Expiwt(:)*auto(:))
    funcE1 = TWO*dt*sum(Expiwt(:)*auto1(:))
    funcE2 = TWO*dt*sum(Expiwt(:)*auto2(:))
    write(nio,*) E+Emin,funcE0%re,funcE1%re,funcE2%re,funcE%re
  END DO
  close(nio)  


END PROGRAM spectre
