!==============================================================================
! spectre.f90 - Transformee de Fourier d'une fonction d'autocorrelation
!
! Formats d'entree acceptes (une ligne par point, commentaires # ou ! permis):
!   2 colonnes : t  Re[C(t)]          (partie imaginaire supposee nulle)
!   3 colonnes : t  Re[C(t)]  Im[C(t)]
!   4 colonnes : label  t  Re[C(t)]  Im[C(t)]  (compatibilite historique)
!
! La grille temporelle doit etre croissante et regulierement espacee.
! La convention utilisee est S(E) = 2 Re integral C(t) exp(+i E t) dt,
! adaptee au cas ou seule la partie t >= 0 de l'autocorrelation est fournie.
!
! Configuration namelist (stdin):
! &param file_auto='auto_cor.txt', file_spectrum='file_spectrum.txt',
!        dE=-1, Emin=0, Emax=-1, conv_E=1, conv_t=1, n_filter=1 /
!
! Compilation: ./script.sh
!==============================================================================
program spectre
  use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, real64, iostat_end
  implicit none

  integer, parameter :: rk = real64
  real(rk), parameter :: pi = acos(-1.0_rk), zero = 0.0_rk, one = 1.0_rk
  complex(rk), parameter :: eye = (zero, one)

  integer :: npts, nptE, ios, i, n_filter
  integer :: unit_out, unit_param
  real(rk) :: dE, Emin, Emax, conv_E, conv_t, E, dt, t0, tmax
  real(rk) :: duration, max_step, tolerance
  real(rk), allocatable :: t(:), quad(:)
  complex(rk), allocatable :: ac(:), ac_window(:)
  complex(rk) :: ft
  character(len=256) :: file_auto, file_spectrum, file_parameters
  namelist /param/ dE, Emin, Emax, conv_E, conv_t, file_auto, file_spectrum, file_parameters, n_filter

  ! Valeurs par defaut : un namelist complet n'est donc pas obligatoire.
  file_auto = 'auto_cor.txt'
  file_spectrum = 'file_spectrum.txt'
  file_parameters = 'spectrum_parameters.txt'
  dE = -one; Emin = zero; Emax = -one
  conv_E = one; conv_t = one; n_filter = 1

  ! Le namelist est facultatif : EOF ou ligne vide conserve les valeurs par defaut.
  read(input_unit, nml=param, iostat=ios)
  if (ios > 0) then
    write(output_unit,'(a,i0)') 'Erreur de lecture du namelist param, code ', ios
    error stop 1
  end if

  if (conv_t <= zero .or. conv_E <= zero .or. n_filter < 0) then
    write(output_unit,'(a)') 'Erreur: conv_t et conv_E doivent etre > 0, n_filter >= 0.'
    error stop 1
  end if

  call count_data_lines(trim(file_auto), npts)
  if (npts < 2) then
    write(output_unit,'(a,i0)') 'Erreur: au moins 2 points valides sont requis (trouves: ', npts, ').'
    error stop 1
  end if

  allocate(t(npts), ac(npts), ac_window(npts), quad(npts))
  call read_data(trim(file_auto), t, ac)
  t = t * conv_t

  do i = 2, npts
    if (t(i) <= t(i-1)) then
      write(output_unit,'(a,i0)') 'Erreur: les temps doivent etre strictement croissants, ligne ', i
      error stop 1
    end if
  end do
  dt = (t(npts)-t(1)) / real(npts-1, rk)
  max_step = maxval(abs((t(2:npts)-t(1:npts-1))-dt))
  tolerance = 100.0_rk * epsilon(one) * max(one, abs(dt))
  if (max_step > tolerance) then
    write(output_unit,'(a,es12.4)') 'Erreur: grille temporelle non reguliere, ecart maximal = ', max_step
    error stop 1
  end if

  t0 = t(1); tmax = t(npts); duration = tmax-t0
  if (Emax < zero) Emax = 2.0_rk*pi/dt / conv_E
  if (Emax <= Emin) then
    write(output_unit,'(a)') 'Erreur: Emax doit etre strictement superieur a Emin.'
    error stop 1
  end if
  if (Emax*conv_E > 2.0_rk*pi/dt * (one + 100*epsilon(one))) then
    write(output_unit,'(a)') 'Erreur: Emax depasse la frequence de Nyquist angulaire.'
    error stop 1
  end if
  ! Resolution naturelle : une maille 4 fois plus fine que 2*pi/duree.
  ! Cela evite une grille spectrale demesuree avec une transformee directe.
  if (dE <= zero) dE = 8.0_rk*pi/duration * conv_E
  if (dE <= zero) then
    write(output_unit,'(a)') 'Erreur: dE doit etre positif.'
    error stop 1
  end if
  nptE = floor((Emax-Emin)/dE) + 1

  ! Quadrature trapezoidale : correction des deux extremites, pas seulement t(1).
  quad = dt
  quad(1) = 0.5_rk*dt
  quad(npts) = 0.5_rk*dt
  ac_window = ac
  if (n_filter > 0) ac_window = ac * exp(-((t-t0)/duration)**n_filter)

  open(newunit=unit_out, file=trim(file_spectrum), status='replace', action='write', iostat=ios)
  if (ios /= 0) then
    write(output_unit,'(a)') 'Erreur: impossible de creer '//trim(file_spectrum)
    error stop 1
  end if
  write(unit_out,'(a)') '# E  Re[S0(E)]  Re[S_cos(E)]  Re[S_cos2(E)]  Re[S_window(E)]'
  do i = 0, nptE-1
    E = Emin + real(i,rk)*dE
    ft = sum(quad * ac * exp(eye*(E/conv_E)*t))
    write(unit_out,'(5(es24.16,1x))') E, 2.0_rk*real(ft), &
      2.0_rk*real(sum(quad * (ac*cos(pi*(t-t0)/(2.0_rk*duration))) * exp(eye*(E/conv_E)*t))), &
      2.0_rk*real(sum(quad * (ac*cos(pi*(t-t0)/(2.0_rk*duration))**2) * exp(eye*(E/conv_E)*t))), &
      2.0_rk*real(sum(quad * ac_window * exp(eye*(E/conv_E)*t)))
  end do
  close(unit_out)

  ! Ecrire les valeurs effectivement utilisees, notamment dE automatique.
  open(newunit=unit_param, file=trim(file_parameters), status='replace', &
       action='write', iostat=ios)
  if (ios /= 0) then
    write(output_unit,'(a)') 'Erreur: impossible de creer '//trim(file_parameters)
    error stop 1
  end if
  write(unit_param,'(a)') '# Parametres effectivement utilises par spectre'
  write(unit_param,'(a,a)') 'file_auto       = ', trim(file_auto)
  write(unit_param,'(a,a)') 'file_spectrum   = ', trim(file_spectrum)
  write(unit_param,'(a,a)') 'dE              = ', trim(real_to_text(dE))
  write(unit_param,'(a,a)') 'Emin            = ', trim(real_to_text(Emin))
  write(unit_param,'(a,a)') 'Emax            = ', trim(real_to_text(Emax))
  write(unit_param,'(a,a)') 'dt              = ', trim(real_to_text(dt))
  write(unit_param,'(a,a)') 't0              = ', trim(real_to_text(t0))
  write(unit_param,'(a,a)') 'tmax            = ', trim(real_to_text(tmax))
  write(unit_param,'(a,i0)') 'npts            = ', npts
  write(unit_param,'(a,i0)') 'nptE            = ', nptE
  write(unit_param,'(a,i0)') 'n_filter        = ', n_filter
  close(unit_param)

  write(output_unit,'(a)') 'Traitement termine avec succes.'
  write(output_unit,'(a,i0)') '  Points temporels : ', npts
  write(output_unit,'(a,es12.4)') '  Pas temporel dt  : ', dt
  write(output_unit,'(a,i0)') '  Points spectraux : ', nptE
  write(output_unit,'(a)') '  Sortie           : '//trim(file_spectrum)
  write(output_unit,'(a)') '  Parametres      : '//trim(file_parameters)

contains

  function real_to_text(value) result(text)
    real(rk), intent(in) :: value
    character(len=40) :: text
    write(text,'(es24.16)') value
  end function real_to_text

  subroutine count_data_lines(filename, n)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: n
    integer :: u, stat, nw
    character(len=1024) :: s
    real(rk) :: q1, q2, q3, q4
    open(newunit=u, file=filename, status='old', action='read', iostat=stat)
    if (stat /= 0) then
      write(output_unit,'(a)') 'Erreur: fichier introuvable: '//trim(filename); error stop 1
    end if
    n = 0
    do
      read(u,'(A)',iostat=stat) s
      if (stat /= 0) exit
      nw = count_words(s)
      if (nw == 2 .or. nw == 3 .or. nw == 4) then
        read(s,*,iostat=stat) q1,q2,q3,q4
        if (stat == 0 .or. stat == iostat_end) n = n+1
      end if
    end do
    close(u)
  end subroutine count_data_lines

  subroutine read_data(filename, tt, cc)
    character(len=*), intent(in) :: filename
    real(rk), intent(out) :: tt(:)
    complex(rk), intent(out) :: cc(:)
    integer :: u, stat, k, nw
    character(len=1024) :: s
    real(rk) :: q1, q2, q3, q4
    open(newunit=u, file=filename, status='old', action='read')
    k = 0
    do
      read(u,'(A)',iostat=stat) s
      if (stat /= 0) exit
      nw = count_words(s)
      if (nw < 2 .or. nw > 4) cycle
      q1=zero; q2=zero; q3=zero; q4=zero
      read(s,*,iostat=stat) q1,q2,q3,q4
      if (stat /= 0 .and. stat /= iostat_end) cycle
      k = k+1
      select case(nw)
      case(2); tt(k)=q1; cc(k)=cmplx(q2,zero,rk)
      case(3); tt(k)=q1; cc(k)=cmplx(q2,q3,rk)
      case(4); tt(k)=q2; cc(k)=cmplx(q3,q4,rk)
      end select
    end do
    close(u)
  end subroutine read_data

  integer function count_words(s) result(n)
    character(len=*), intent(in) :: s
    integer :: j, L
    logical :: inside
    character :: c
    n=0; inside=.false.; L=len_trim(s)
    do j=1,L
      c=s(j:j)
      if (c=='#' .or. c=='!') exit
      if (c==' ' .or. c==char(9) .or. c==',') then
        inside=.false.
      else if (.not.inside) then
        n=n+1; inside=.true.
      end if
    end do
  end function count_words
end program spectre
