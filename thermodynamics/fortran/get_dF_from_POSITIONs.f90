
PROGRAM getDistanceForcesFromPOSITIONs

      implicit none
      INTEGER i,j,nat,nmd
      REAL(8) dummy
      REAL(8) dist,force,cell(3,3),delta(3),forcevec(3),deltaForce(3)
      REAL(8),allocatable :: coords(:,:), forces(:,:)
      integer,allocatable :: map(:,:)

      ! get number of atoms
      open(unit=10,file='mapping',status='old',err=99)
      nat = 0
  10  read (10,*,end=20) dummy
      nat = nat+1
      go to 10
  20  continue
      if (nat.eq.0) go to 99
      rewind(10)

      allocate(coords(3,nat),forces(3,nat),map(nat,nat))

      call readMap(nat,map,'mapping')
      call readCoords(3,cell,'cell')

      open(unit=10,file="POSITIONs",status="old")
      open(unit=11,file="force_vs_distance",status="unknown")

      DO
       DO i=1,nat
         read(10,*,end=30) coords(1,i),coords(2,i),coords(3,i),forces(1,i),forces(2,i),forces(3,i)
         call toReduced(coords(:,i),cell)
       ENDDO

       DO i=1,nat
        DO j=1,nat
         if (map(i,j).NE.-1) then
           call periodicDiffVec(coords(:,map(i,j)),coords(:,i),delta)
           call toCartesian(delta,cell)
           deltaForce(:) = forces(:,map(i,j)) + forces(:,i)
           forcevec(:) = dot_product(delta(:),deltaForce(:))*delta(:)
           dist = sqrt(sum(delta**2.))
           force = sign(sqrt(sum(forcevec**2.)),dot_product(delta(:),deltaForce(:)))
          if (i==1) then
            write (*,*) i,map(i,j)
            write (*,*) delta
            write (*,*) forcevec
            write (*,*) dist,force
          endif
           write(11,*) dist, force
         else
           exit
         endif
        ENDDO
       ENDDO
      ENDDO

      close(10)
      close(11)

  30  STOP

  99  continue
      write (*,*) "ERROR: cannot open or empty file mapping"
      stop
END


! v1 and v2 in reduced coords expected; delta returned in reduced coords
SUBROUTINE periodicDiffVec(v1,v2,delta)
      implicit none
      real(8) v1(3),v2(3),delta(3)
      integer i

      delta = v1 - v2
      delta = delta - int(delta)
      do i=1,3
        if (delta(i)>0.5) then
          delta(i) = delta(i) - 1
        else
          if (delta(i)<-0.5) delta(i) = delta(i) + 1
        endif
      enddo
END


SUBROUTINE readCoords (NIONS,coords,fileName)
      
      implicit none
      INTEGER NIONS,iostatus
      REAL(8) coords(3,NIONS)
      CHARACTER :: fileName*(*)

      LOGICAL exists
      INTEGER i, TIU0
      CHARACTER*1000 :: line
      REAL(8) check(4)
      REAL(8), PARAMETER :: HUGE=100000._8

      TIU0 = 6 ! STDOUT

      open(unit=11,file=fileName,status='old',err=99)

      check=HUGE
      do i=1,NIONS
	READ (11,'(A)',END=99) line
        READ (line,*,IOSTAT=iostatus) check(:)
	if (iostatus>0.OR.check(3).EQ.HUGE.OR.check(4).NE.HUGE) go to 99
	coords(:,i) = check(1:3)
	check=HUGE
      end do
      READ (11,'(A)',IOSTAT=iostatus) line
      if (iostatus.NE.-1) go to 99
      close (11)

      RETURN

  99  continue
      write (TIU0,*) ""
      write (TIU0,*) "ERROR: cannot open or corrupted "//fileName//" file"
      write (TIU0,*) "(corrupt: too many or few lines or columns or nonnumeric entrees)"
      STOP
  
END


SUBROUTINE readMap(N,map,fileName)

      implicit none 
      INTEGER N, TIU0
      integer map(N,N)
      CHARACTER :: fileName*(*)

      INTEGER i, j, iostatus, pos1, pos2
      CHARACTER*200 :: line

      TIU0 = 6 ! STDOUT

      map=-1
      open(unit=10,file=fileName,err=99)
      do i=1,N
	READ (10,'(A)',END=99) line
        pos1 = 1
        j = 0
        do
          pos2 = index(line(pos1:)," ")
          if (pos2==0) then
            j = j+1
            read (line(pos1:),*,IOSTAT=iostatus) map(i,j)
            exit
          endif
          j = j+1
          read (line(pos1:pos1+pos2-2),*,IOSTAT=iostatus) map(i,j)
          pos1 = pos2+pos1
        enddo
      end do
      READ (10,'(A)',IOSTAT=iostatus) line
      if (iostatus.NE.-1) go to 99
      close(10)

      RETURN

  99  continue
      write (TIU0,*) ""
      write (TIU0,*) "ERROR: cannot open or corrupted "//fileName
      write (TIU0,*) "(corrupt: too many or few lines or columns or nonnumeric entrees)"
      STOP
END

SUBROUTINE toCartesian(coord,cell)

      real(8) coord(3)
      real(8) cell(3,3)

      coord = MATMUL(transpose(cell),coord)
      return
END SUBROUTINE


SUBROUTINE toReduced(coord,cell)

      real(8) coord(3)
      real(8) cell(3,3),inv(3,3)

      call invert (transpose(cell),inv)

      coord = MATMUL(inv,coord)
      return
END SUBROUTINE


SUBROUTINE INVERT(MATR,MATR2)
      IMPLICIT NONE
      real(8) :: MATR(3,3)
      real(8) :: MATR2(3,3)
      real(8) :: V1(3),V2(3),V3(3),b1(3),b2(3),b3(3),A1(3),A2(3),A3(3)

      V1 = MATR(1,:)
      V2 = MATR(2,:)
      V3 = MATR(3,:)

      b1(1) = V2(2)*V3(3)-V3(2)*V2(3)
      b1(2) = V2(3)*V3(1)-V3(3)*V2(1)
      b1(3) = V2(1)*V3(2)-V3(1)*V2(2)

      b2(1) = V3(2)*V1(3)-V1(2)*V3(3)
      b2(2) = V3(3)*V1(1)-V1(3)*V3(1)
      b2(3) = V3(1)*V1(2)-V1(1)*V3(2)

      b3(1) = V1(2)*V2(3)-V2(2)*V1(3)
      b3(2) = V1(3)*V2(1)-V2(3)*V1(1)
      b3(3) = V1(1)*V2(2)-V2(1)*V1(2)
      
      A1 = b1/DOT_PRODUCT(b1,V1)
      A2 = b2/DOT_PRODUCT(b2,V2)
      A3 = b3/DOT_PRODUCT(b3,V3)

      MATR2(:,1) = A1
      MATR2(:,2) = A2
      MATR2(:,3) = A3
END SUBROUTINE INVERT

