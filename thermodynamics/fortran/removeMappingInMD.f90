PROGRAM removeMappingInMD

   implicit none
   INTEGER nAtoms, force, n, i
   REAL(8), ALLOCATABLE :: prev(:,:), atom(:)
   real(8) :: cell(3,3), inv(3,3)
   character(len=200) :: str,fileIn
   external :: invert
   
   call getarg( 1,fileIn)
   call getarg( 2,str); read (str,*) nAtoms
   call getarg( 3,str); read (str,*) cell(1,1)
   call getarg( 4,str); read (str,*) cell(1,2)
   call getarg( 5,str); read (str,*) cell(1,3)
   call getarg( 6,str); read (str,*) cell(2,1)
   call getarg( 7,str); read (str,*) cell(2,2)
   call getarg( 8,str); read (str,*) cell(2,3)
   call getarg( 9,str); read (str,*) cell(3,1)
   call getarg(10,str); read (str,*) cell(3,2)
   call getarg(11,str); read (str,*) cell(3,3)
   call getarg(12,str); read (str,*) force

   cell = transpose (cell)
   call invert (cell,inv)

   if (force.eq.0) then; n=3; else; n=6; endif
   ALLOCATE (prev(nAtoms,n),atom(n))

   open (8,FILE=fileIn)
   open (9,FILE=trim(adjustl(fileIn))//"_noJumps")
   do i=1,nAtoms
     read (8,*) prev(i,:)
					if (n==3) write (9,'(3F)') prev(i,:)
					if (n==6) write (9,'(6F)') prev(i,:)
   end do

   do while (.true.)
     do i=1,nAtoms
        read (8,*,end=11) atom(:)
        ! transform to reduced coords
        prev(i,1:3) = MATMUL (inv,prev(i,1:3))
        atom(1:3) = MATMUL (inv,atom(1:3))
        ! x direction
        if (atom(1)-prev(i,1)> 1./2) atom(1)=atom(1)-1
        if (atom(1)-prev(i,1)<-1./2) atom(1)=atom(1)+1
        ! y direction
        if (atom(2)-prev(i,2)> 1./2) atom(2)=atom(2)-1
        if (atom(2)-prev(i,2)<-1./2) atom(2)=atom(2)+1
        ! z direction
        if (atom(3)-prev(i,3)> 1./2) atom(3)=atom(3)-1
        if (atom(3)-prev(i,3)<-1./2) atom(3)=atom(3)+1
        ! back to cartesian coords
        prev(i,1:3) = MATMUL (cell,prev(i,1:3))
        atom(1:3) = MATMUL (cell,atom(1:3))
        prev(i,:) = atom(:)
        if (n==3) write (9,'(3F)') atom(:)
        if (n==6) write (9,'(6F)') atom(:)
     end do
   end do

11 continue
   close (9)
   close (8)
   
END PROGRAM

