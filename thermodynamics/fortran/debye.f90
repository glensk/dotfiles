program Debye

  implicit none

  external :: getSplineCoefficients, ispline
  real(8) ispline
  real(8), parameter :: kB=0.086173422, GPaTokBar=10, Ang3ToBohr3=6.7483345, pi=3.1415926535898
  real(8), parameter :: convert=67.48    ! conversion factor, Moruzzi Eq. (4)
  real(8), parameter :: empirical=0.617  ! empirical factor, Moruzzi Eq. (6)
  real(8), parameter :: gammalow=1.      ! low T gamma
  real(8), parameter :: gammahigh=2./3.  ! high T gamma

  integer nV,nx,k1,k2,i,j,k,ind

  real(8) V0,B0,Bprime,mass,Vstart,Vstep,Tstart,Tstep,Tmax,V,T,xStart,xStep,r0,Vfix
  real(8) debyeTlow,debyeThigh,debyeT(2),debyeZero,d,Fqh,x
  real(8), allocatable :: debyeIntegral(:), xx(:), bb(:), cc(:), dd(:)
  character(len=20) :: str,dummy

  call getarg(1,str);  read (str,*) V0
  call getarg(2,str);  read (str,*) B0
  call getarg(3,str);  read (str,*) Bprime
  call getarg(4,str);  read (str,*) mass
  call getarg(5,str);  read (str,*) Vstart
  call getarg(6,str);  read (str,*) Vstep
  call getarg(7,str);  read (str,*) nV
  call getarg(8,str);  read (str,*) Tstart
  call getarg(9,str);  read (str,*) Tstep
  call getarg(10,str); read (str,*) Tmax
  call getarg(11,str); read (str,*) Vfix

  ! load tabulated Debye function
  open(unit=10,file="Debye",status="old")
  read(10,*) dummy,xStart,xStep,nx
  allocate(debyeIntegral(nx),xx(nx),bb(nx),cc(nx),dd(nx))
  do i=1,nx
    xx(i) = xStart + (i-1)*xStep
    read(10,*) debyeIntegral(i)
  enddo
  close(10)
  call getSplineCoefficients(xx,debyeIntegral,bb,cc,dd,nx)

  r0 = (3*V0*Ang3ToBohr3/(4*pi))**(1./3.)
  debyeZero = empirical * convert * sqrt(r0*B0*GPaTokBar/mass)

  write (*,*) "debyeZero: ",debyeZero," K"

  open(unit=12,file="DebyeT_vs_V")
  do i=0,nV
    if (i==0) then
      V = Vfix
    else
      V = Vstart + (i-1)*Vstep
    endif
    debyeTlow  = debyeZero * (V0/V)**(-gammalow  + 0.5*(1+Bprime))
    debyeThigh = debyeZero * (V0/V)**(-gammahigh + 0.5*(1+Bprime))
    debyeT = (/ debyeTlow, debyeThigh /)
    if (i.NE.0) write (12,*) V, debyeTlow, debyeThigh

    write (str,'(F20.8)') V
    str = adjustl(str)
    ind = index(str,"0 ")
    do while (ind > 0)
      str = str(1:ind-1)
      ind = index(str,"0 ")
    enddo
    open(unit=10,file="FqhDebye-lowT_"//trim(str))
    open(unit=11,file="FqhDebye-highT_"//trim(str))
    do j=1,ceiling((Tmax-Tstart)/Tstep+1)
     do k=1,2
      T = Tstart + (j-1)*Tstep
      x = debyeT(k)/T
      if (x>xx(nx)) then
        d = 3/x**3 * debyeIntegral(nx)
      else
        d = 3/x**3 * ispline(x,xx,debyeIntegral,bb,cc,dd,nx)
      endif

      ! Debye free energy, corresponds to Eq. (19) from Moruzzi, but NOTE that
      ! there is an error in Moruzzi, the 3ln term must be positive in the end (i.e., when the brackets are eliminated)
      ! so inside the brackets the 3ln term should be negative to cancel the minus in front of the brackets
      Fqh = (9./8.) * kB*debyeT(k) + T*kB * (3*log(1 - exp(-x)) - d )
      write (9+k,*) T,Fqh
     enddo
    enddo
    close(10)
    close(11)
  enddo
  close(12)
end
