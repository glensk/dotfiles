PROGRAM generateLattice

   implicit none
   INTEGER sc,n,i,j,k
   real(8) :: a, R1(3), R2(3)
   character(len=20) :: str,fileIn
   
   call getarg(1,str); read (str,*) a
   call getarg(2,str); read (str,*) sc

   open (9,FILE="coords",status='old',position='append')

   n=1
   do i=0,sc-1
     do j=0,sc-1
       do k=0,sc-1
         R1 = (/ i*a, j*a, k*a /)
         R2 = (/ (i+0.5)*a, (j+0.5)*a, (k+0.5)*a /)
         write (9,'(I,I,3F15.8)') n,   1, R1(:)
         write (9,'(I,I,3F15.8)') n+1, 1, R2(:)
         n=n+2
       end do
     end do
   end do

   close (9)
   
END PROGRAM

