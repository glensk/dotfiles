program densityOfStates
    use qsort_c_module

    implicit none
    integer :: nPoints, n, i, j, k, prog
    character(len=40) :: specFile, str
    real :: scaleSmearing, scaleLimits, dummy, minValueX, minValueY, maxValueX, maxValueY
    real :: smearingX,smearingY,norm,dx,dy
    real, allocatable :: spectrum(:,:), x(:), y(:), dos(:,:)
    type(list) l(1)
 
    nPoints       = 100
    scaleSmearing = 0.05
    scaleLimits   = 0.15
 
    if ( iargc() == 0) then
        write (*,100) scaleSmearing,scaleSmearing
        ! getDOS2d.sh is the shell wrapper to call dos2d.x
        100 format (/"usage:"/ &
                   & "   1)  getDOS2d.sh spectrumFile"/ &
                   & "or 2)  getDOS2d.sh spectrumFile smearingX smearingY"/ &
                   & ""/ &
                   & "if not given defaults are:"/ &
                   & "  smearingX = ",F6.3," x (maxOfSpectrumX - minOfSpectrumX)"/ &
                   & "  smearingY = ",F6.3," x (maxOfSpectrumY - minOfSpectrumY)")
       call exit
    end if
 
    call getarg(1,specFile)
    open (unit=10,file=specFile,status='old',err=99)
    n = 0
 10 read (10,*,end=20) dummy
    n = n+1
    go to 10
 20 continue
    if (n.eq.0) go to 99
    rewind (10)

    allocate (spectrum(n,2))
    do i=1,n
      read(10,*) spectrum(i,1), spectrum(i,2)
    end do
    close(10)

    allocate(l(1)%A(n))
    l(1)%A=spectrum(:,1)
    call QsortC(l)
    minValueX = l(1)%A(1) - scaleLimits * (l(1)%A(n)-l(1)%A(1))
    maxValueX = l(1)%A(n) + scaleLimits * (l(1)%A(n)-l(1)%A(1))

    if ( iargc()==1 ) then
      smearingX = l(1)%A(2)-l(1)%A(1)
      do i=3,n
        if (l(1)%A(i)-l(1)%A(i-1)<smearingX) smearingX = l(1)%A(i) - l(1)%A(i-1)
      end do
      smearingX = scaleSmearing * (maxValueX-minValueX)
    else
      call getarg(2,str); read (str,*) smearingX
    end if

    l(1)%A=spectrum(:,2)
    call QsortC(l)
    minValueY = l(1)%A(1) - scaleLimits * (l(1)%A(n)-l(1)%A(1))
    maxValueY = l(1)%A(n) + scaleLimits * (l(1)%A(n)-l(1)%A(1))

    if ( iargc()==1 ) then
      smearingY = l(1)%A(2)-l(1)%A(1)
      do i=3,n
        if (l(1)%A(i)-l(1)%A(i-1)<smearingY) smearingY = l(1)%A(i) - l(1)%A(i-1)
      end do
      smearingY = scaleSmearing * (maxValueY-minValueY)
    else
      call getarg(3,str); read (str,*) smearingY
    end if

    allocate (dos(nPoints,nPoints),x(nPoints),y(nPoints))
    norm = 0.
    write (*,*)
    write (*,*) "Progress:"
    prog = 20
    do j=1,nPoints
     do k=1,nPoints
      dos(j,k) = 0.;
      do i=1,n
        x(j) = minValueX + j*(maxValueX-minValueX)/nPoints
        y(k) = minValueY + k*(maxValueY-minValueY)/nPoints
        dos(j,k) = dos(j,k) + exp(-(x(j)-spectrum(i,1))**2/smearingX**2 &
                 &                -(y(k)-spectrum(i,2))**2/smearingY**2)
      end do 
      if (j>1) then
        dx = x(j)-x(j-1)
        dy = y(k)-y(k-1)
        norm = norm + dos(j,k)*dx*dy
      endif
     end do
     if (1.*j/nPoints*100.GE.prog) then
       write (*,*) prog,"%"
       prog = prog+20
     endif
    end do

    open (unit=11,file="DOS2d")
    do j=1,nPoints
     do k=1,nPoints
      write (11,*) x(j),y(k),dos(j,k)/norm
     end do
    end do
    close(11)
    write (*,*)
    write (*,*) " smearingX      smearingY"
    write (*,*)  smearingX,smearingY
    write (*,*)
    write (*,*) " DOS2d written"
    stop
 
 99 continue
    write (*,*) "error: cannot open or empty file ",specFile
    stop
end
