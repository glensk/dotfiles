program densityOfStates_weighted
    use qsort_c_module

    implicit none
    integer :: nPoints, n, i, j, c, nc
    character(len=40) :: specFile, str
    real(8) :: scaleSmearing, scaleLimits, dummy, minValue, maxValue
    real(8) :: pivot, smearing,norm,dx
    real(8), allocatable :: w(:),spectrum(:), x(:), dos(:)
    type(list) l(2)
    logical useWeights
 
    nPoints       = 500
    scaleSmearing = 0.05
    scaleLimits   = 0.15
 
    if ( iargc() == 0) then
        write (*,100) scaleSmearing,scaleLimits,scaleLimits,nPoints
        ! getDOS.sh is the shell wrapper to call dos.x
        100 format (/"usage:"/ &
                   & "   1)  getDOS.sh spectrumFile"/ &
                   & "or 2)  getDOS.sh spectrumFile smearing"/ &
                   & "or 3)  getDOS.sh spectrumfile smearing min max nPoints"// &
                   & "if not given defaults are:"/ &
                   & "  smearing = ",F6.3," x (maxOfSpectrum - minOfSpectrum)"/ &
                   & "  min      = minOfSpectrum - ",F6.3," x (maxOfSpectrum - minOfSpectrum)"/ &
                   & "  max      = maxOfSpectrum + ",F6.3," x (maxOfSpectrum - minOfSpectrum)"/ &
                   & "  nPoints  = ",I6/ &
                   & / &
                   & "Option: ""-w""  for weighted dos with weights taken from ""weights"" file")
       call exit
    end if
 
    call getarg(1,specFile)
    if (specFile=="-w") then
      useWeights=.true.
      nc=iargc()-1
      call getarg(2,specFile)
      c=2
    else
      useWeights=.false.
      nc=iargc()
      c=1
    endif

    open (unit=10,file=specFile,status='old',err=99)
    n = 0
 10 read (10,*,end=20) dummy
    n = n+1
    go to 10
 20 continue
    if (n.eq.0) go to 99
    rewind (10)

    if (n>2000000) then
      write (*,*) ""
      write (*,*) "spectrum file contains more than 2.000.000 entries"
      write (*,*) "this is to large to be handled by quicksort routine"
      write (*,*) "EXITING"
      stop
    endif

    allocate (spectrum(n))
    do i=1,n
      read(10,*) spectrum(i)
    end do
    close(10)

    allocate(w(n))
    if (useWeights) then
      specFile='weights'
      open(unit=11,file=specFile,status='old',err=99)
      do i=1,n
        read(11,*) w(i)
      enddo
      close(11)
    else
      w=1.
    endif

    allocate(l(1)%A(n),l(2)%A(n))
    l(1)%A=spectrum
    l(2)%A=w

    ! QsortC routine produces SegFault for very large lists (somewhere about 2.000.000 to 3.000.000 entries)
    call QsortC(l)
    spectrum=l(1)%A
    w=l(2)%A

    if ( nc < 5 ) then
      minValue = spectrum(1) - scaleLimits * (spectrum(n)-spectrum(1))
      maxValue = spectrum(n) + scaleLimits * (spectrum(n)-spectrum(1))
    else
      call getarg(c+2,str); read (str,*) minValue
      call getarg(c+3,str); read (str,*) maxValue
      call getarg(c+4,str); read (str,*) nPoints
    end if

    if ( nc < 2 ) then
      smearing = spectrum(2)-spectrum(1)
      do i=3,n
        if (spectrum(i)-spectrum(i-1)<smearing) smearing = spectrum(i) - spectrum(i-1)
      end do
      smearing = scaleSmearing * (maxValue-minValue)
    else
      call getarg(c+1,str); read (str,*) smearing
    end if

    allocate (dos(nPoints),x(nPoints))
    norm = 0.
    do j=1,nPoints
      dos(j) = 0.;
      do i=1,n
        x(j) = minValue + j*(maxValue-minValue)/nPoints
        dos(j) = dos(j) + w(i)*exp(-(x(j)-spectrum(i))**2/smearing**2)
      end do 
      if (j>1) then
        dx = x(j)-x(j-1)
        norm = norm + dos(j)*dx
      endif
    end do

    open (unit=11,file="DOS")
    do j=1,nPoints
      write (11,*) x(j),dos(j)/norm
    end do
    close(11)
    write (*,*)
    write (*,*) " min           minOfSpectrum      maxOfSpectrum   max"
    write (*,*) minValue,spectrum(1),"   ",spectrum(n),maxValue
    write (*,*)
    write (*,*) " max-min        smearing          nPoints"
    write (*,*) maxValue-minValue,smearing,nPoints
    write (*,*)
    write (*,*) "DOS written"
    write (*,*)
    write (*,*) "Note: if used with getEIGENVALUES.sh scale by (3-s)*NBANDS/NIONS,"
    write (*,*) "      where s is the number of spins, to get unit (states/eV*atom)"
    stop
 
 99 continue
    write (*,*) "error: cannot open or empty file ",specFile
    stop
end
