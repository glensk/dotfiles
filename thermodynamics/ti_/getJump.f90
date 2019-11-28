program getJump

     implicit none
     external :: invert
     real(8) :: cell(3,3), inv(3,3)
     real(8) :: dummy,vol,dist,radius,pi,dist2
     real(8), allocatable :: coords(:,:,:),eq(:,:),eq2(:,:)
     integer :: nAll,nAtoms,i,j,dr,dd,k

     pi = 3.14159265358979312

     open (unit=10,file="cell",status="old", err=99)
     do i=1,3
        read (10,*) cell(i,:)
     end do
     close(10)

     cell = transpose (cell)
     call invert (cell,inv)

     open (unit=10,file="atoms_volume_steps",status="old",err=99)
     read (10,*) nAtoms,vol,nAll
     close(10)

     allocate(coords(nAll-1,nAtoms,3),eq(nAtoms,3),eq2(nAtoms+1,3))
     open (unit=11,file="POSITIONs",status="old",err=99)
     do i=1,nAtoms
        read (11,*) eq(i,:)
        eq2(i,:)=eq(i,:)
     end do
     eq2(nAtoms+1,1)=0
     eq2(nAtoms+1,2)=0
     eq2(nAtoms+1,3)=0

     do j=1,nAll-1
       do i=1,nAtoms
         read (11,*) coords(j,i,:)
       end do
     end do
     close(11)

     radius = (3*vol/nAtoms/4/pi)**0.3333333333
     !radius =1.3*radius
     !radius = 2*radius/3

     open (unit=12,file="jumps_dist",status="replace",err=100)
     open (unit=13,file="jumps_radius",status="replace",err=100)
     do j=1,nAll-1
       dr=0; dd=0
       do i=1,nAtoms
         ! look here for atoms moving out of the given radius
         call periodicDistance(coords(j,i,:),eq(i,:),cell,inv,dist)
         if (dist>radius) dr=dr+1
         do k=1,nAtoms+1
           if (k.ne.j) then
             ! look here for atoms being closer to other equilibrium positions than their own
             call periodicDistance(coords(j,i,:),eq2(k,:),cell,inv,dist2)
             if (dist2<dist) then; dd=dd+1; exit; end if
           end if
         end do
       end do
       !write (12,*) dist
       write (12,*) dd
       write (13,*) dr
     end do
     close(12)
     close(13)

     stop     

 99  continue
     write (*,*) "error: cannot open or empty cell or POSITIONs file"
     stop
 100 continue
     write (*,*) "error: cannot open file for writing"
     stop
end


subroutine periodicDistance (c1,c2,cell,inv,dist)
  real(8) :: c1(3), c2(3), cell(3,3), inv(3,3), dist, diff(3)
  integer i

  diff = MATMUL(inv,c1) - MATMUL(inv,c2)
  do i=1,3
    if (diff(i)>0.or.diff(i)==0) diff(i) = diff(i) - floor(diff(i))
    if (diff(i)<0)               diff(i) = diff(i) - ceiling(diff(i))
    if (diff(i)> 0.5) diff(i)=1-diff(i)
    if (diff(i)<-0.5) diff(i)=1+diff(i)
  end do
  diff = MATMUL(cell,diff)
  dist = (diff(1)**2 + diff(2)**2 + diff(3)**2)**0.5
end


