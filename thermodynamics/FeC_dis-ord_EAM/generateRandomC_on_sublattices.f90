PROGRAM generateRandomC

   implicit none
   INTEGER sc,nC,nLast,i,j,k,n,nsublat,ns
   real(8) :: a, rRepulsive, randx, randy, randz, dist(3), sublat(3,3), seedini
   REAL(8), ALLOCATABLE :: rList(:,:)
   character(len=20) :: str,fileIn
   logical notOk

   integer, allocatable :: seed(:)
   integer size
   
   call getarg(1,str); read (str,*) nC           ! number of C atoms
   call getarg(2,str); read (str,*) rRepulsive   ! radius for avoiding C-C interactions
   call getarg(3,str); read (str,*) a            ! lattice constant
   call getarg(4,str); read (str,*) sc           ! supercell size
   call getarg(5,str); read (str,*) seedini      ! seed init
   call getarg(6,str); read (str,*) nsublat      ! 1: put nC/1 C atoms on first sublattic
                                                 ! 2: put nC/2 C atoms on first and second sublattice
                                                 ! 3: put nC/3 C atoms on first and second and third sublattice
   sublat = 0.
   sublat(1,1) = 0.5
   sublat(2,2) = 0.5
   sublat(3,3) = 0.5

   call random_seed(size=size)
   allocate(seed(size))
   seed(:) = seedini
   call random_seed(put=seed)

!   call random_seed
   allocate(rList(nC,3))

   open (9,FILE="coords",status='old',position='append')
   backspace (9)
   read (9,'(I)') nLast

   do ns=1,nsublat
     do i=1,nC/nsublat
       notOk=.true.
       n=1
       do while (notOk)
         n = n+1
         call random_number(randx)
         call random_number(randy)
         call random_number(randz)
         randx = (nint(sc*randx)+sublat(ns,1))/sc
         randy = (nint(sc*randy)+sublat(ns,2))/sc
         randz = (nint(sc*randz)+sublat(ns,3))/sc
         notOk=.false.
         do j=2,i
           call periodicDiffVec(rList(j,:),(/ randx, randy, randz /), dist)
           if (sqrt(sum((sc*dist)**2))<rRepulsive) notOk=.true. 
         enddo
         if (n>sc**3) then
           write (*,*) "cannot place C atoms apart"
           stop
         endif
       enddo
       rList(i,1) = randx
       rList(i,2) = randy
       rList(i,3) = randz
       write (9,'(I,I,3F15.8)') nLast+(ns-1)*nc/nsublat+i, 2, a*randx*sc, a*randy*sc, a*randz*sc
     end do
   end do

   close (9)
   
END PROGRAM

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

