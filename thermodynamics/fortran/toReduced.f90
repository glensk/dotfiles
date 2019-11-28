program cartesianToReducedCoordinates

     implicit none
     external :: invert
     real(8) :: cell(3,3), inv(3,3)
     real(8) :: dummy
     real(8), allocatable :: coords(:,:)
     integer :: i,n

     open (unit=10,file="cell",status="old", err=99)
     do i=1,3
        read (10,*) cell(i,:)
     end do
     close(10)

     cell = transpose (cell)
     call invert (cell,inv)

     open (unit=11,file="cartesian_coords",status="old",err=99)
     n=0
 10  read (11,*,end=20) dummy
     n=n+1
     go to 10
 20  continue
     if (n.eq.0) go to 99
     rewind (11)
     allocate(coords(n,3))
     do i=1,n
        read (11,*) coords(i,:)
     end do
     close(11)

     open (unit=12,file="reduced_coords",status="replace")
     do i=1,n
        coords(i,:) = MATMUL (inv,coords(i,:))
        write (12,*) coords(i,:)
     end do
     close (12)

     write (*,*)
     write (*,*) "reduced_coords written"
     stop     

 99  continue
     write (*,*) "error: cannot open or empty cell or cartesian_coords file"
end
