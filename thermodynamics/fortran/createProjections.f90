PROGRAM createProjections
   CHARACTER*1 dummy
   INTEGER n, nAtoms
   REAL, ALLOCATABLE :: EigVec(:,:), eqStr(:), str(:)
   REAL projection
   
   open (7,FILE="eigVecs.dat")
   n=0
   do while (.true.)
      read (7,*,end=10) dummy
      n=n+1
   end do
10 continue
   close (7)

   ALLOCATE (EigVec(1:n,n),eqStr(1:n),str(1:n))

   call readEigVec (n,"eigVecs.dat",EigVec)
   nAtoms = n/3

   open (8,FILE="POSITIONs")
   do i=1,nAtoms
     read (8,*) eqStr(3*(i-1)+1),eqStr(3*(i-1)+2),eqStr(3*(i-1)+3)
   end do

   open (9,FILE="projections.dat")
   do i=1,nAtoms
      read (8,*) str(3*(i-1)+1),str(3*(i-1)+2),str(3*(i-1)+3)
      do j=1,3
         str(3*(i-1)+j)=str(3*(i-1)+j)-eqStr(3*(i-1)+j)
      end do
   end do

   do i=1,n
      projection=0
      do j=1,n
        projection=projection+str(j)*EigVec(i,j)
      end do
      write (9,"(F)",advance="no") projection
   end do
   write (9,*) ""

11 continue
   close (9)
   close (8)
   
END PROGRAM

SUBROUTINE readEigVec (n,f,EigVec)
   CHARACTER*(*) f
   REAL EigVec(n,n)

   open (7,FILE=f)
   DO i=1,n
      write (6,*) n
      read (7,*) EigVec(i,:)
   END DO
   close(7)
END
