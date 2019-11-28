
PROGRAM getEinstein

      use qsort_c_module

      implicit none
      INTEGER i,j,k,ii,nat,nmd,nPoints
      real(8) :: pivot, smearing,norm,dx,maxDOS(2),halfDOS1(2),halfDOS2(2)
      real(8) :: scaleSmearing, scaleLimits, dummy, minValue, maxValue
      REAL(8),allocatable :: coords(:,:,:), spectrum(:), x(:), dos(:), eqCoords(:,:), einstein(:,:)
      type(list) l(1)
      character(len=20) :: str,filename,si,sk

      nPoints       = 500
      scaleSmearing = 0.05
      scaleLimits   = 0.15

      call getarg(1,str); read (str,*) nat

      open(unit=10,file="pre_POSITIONs_noJumps",status="old")
      nmd = 0
  10  read (10,*,end=20) dummy
      nmd = nmd+1
      go to 10
  20  continue
      if (nmd.eq.0) go to 99
      rewind(10)
      nmd = nmd/nat

      allocate(coords(nmd,nat,3),spectrum(nmd))
      allocate(l(1)%A(nmd))
      allocate(dos(nPoints),x(nPoints))
      allocate(eqCoords(nat,3),einstein(nat,3))

      DO j=1,nmd
       DO i=1,nat
         read(10,*,end=30) coords(j,i,1),coords(j,i,2),coords(j,i,3)
       ENDDO
      ENDDO

  30  continue
      close(10)

      DO i=1,nat
        write (*,*) "atom", i
        DO k=1,3
          spectrum(:) = coords(:,i,k)

          l(1)%A=spectrum
          ! QsortC routine produces SegFault for very large lists (somewhere about 2.000.000 to 3.000.000 entries)
          call QsortC(l)
          spectrum=l(1)%A

          minValue = spectrum(1) - scaleLimits * (spectrum(nmd)-spectrum(1))
          maxValue = spectrum(nmd) + scaleLimits * (spectrum(nmd)-spectrum(1))

          smearing = spectrum(2)-spectrum(1)
          do ii=3,nmd
            if (spectrum(ii)-spectrum(ii-1)<smearing) smearing = spectrum(ii) - spectrum(ii-1)
          end do
          smearing = scaleSmearing * (maxValue-minValue)

          maxDOS = 0.
          do j=1,nPoints
            dos(j) = 0.;
            do ii=1,nmd
              x(j) = minValue + j*(maxValue-minValue)/nPoints
              dos(j) = dos(j) + exp(-(x(j)-spectrum(ii))**2/smearing**2)
            end do 
            if (dos(j)>maxDOS(2)) then
              maxDOS(1)=x(j)
              maxDOS(2)=dos(j)
            endif
          end do

          halfDOS1 = 0.
          halfDOS2 = 0.
          do j=1,nPoints
            if (halfDOS1(2)==0) then
              if (dos(j)>maxDOS(2)/2) then
                halfDOS1(1)=x(j)
                halfDOS1(2)=dos(j)
              endif
            else
              if (halfDOS2(2)==0) then
                if (dos(j)<maxDOS(2)/2) then
                  halfDOS2(1)=x(j)
                  halfDOS2(2)=dos(j)
                  exit
                endif
              endif
            endif
          enddo
          
          eqCoords(i,k) = maxDOS(1)
          einstein(i,k) = (halfDOS2(1) - halfDOS1(1))/2

!         write (si,'(I3.3)') i
!         write (sk,'(I3.3)') k
!         filename='DOS_'//trim(si)//'_'//trim(sk)
!         open (unit=11,file=filename)
!         do j=1,nPoints
!           write (11,*) x(j),dos(j)
!         end do
!         close(11)

        ENDDO
      ENDDO

      open(unit=11,file="eqCoords")
      open(unit=12,file="einstein")
      DO i=1,nat
        write(11,*) eqCoords(i,:)
        write(12,*) einstein(i,:)
      ENDDO
      close(10)
      close(11)
      
      stop


  99  continue
      write (*,*) "ERROR: cannot open or empty file pre_POSITIONs_noJumps"
      stop
END



