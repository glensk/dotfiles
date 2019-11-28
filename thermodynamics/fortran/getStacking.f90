program getStacking

     implicit none
     external :: invert
     real(8) :: cell(3,3),inv(3,3),diff(3),primCell(3,3),primInv(3,3),enhancement
     real(8) :: dummy,vol,dist,radius,pi,dist2,delta,red(3),rx,ry,diffA(2),diffB(2),diffC(2),wMin,wMax
     real(8), allocatable :: coords(:,:,:),eq(:,:),redLayer(:,:),zLayer(:),wLayerA(:),wLayerB(:),wLayerC(:),mapping(:,:)
     integer :: nAll,nAtoms,i,j,dr,dd,k,n,Ci
     integer, allocatable :: layers(:), nAtomsPerLayer(:),newOrder(:)
     logical alreadyIn,verbose,Cpresent
     character(len=20) :: str, stacking, phase

     delta = 0.001
     enhancement = 0.01
     wMin = 0.5
     wMax = 1.5

     str=""
     verbose=.FALSE.
     if (iargc()==1) call getarg(1,str) ! must be "v" at input for verbose output
     if (trim(str)=="v") verbose=.TRUE.

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

     allocate(coords(nAll-1,nAtoms,3),eq(nAtoms,3),layers(nAtoms),mapping(nAtoms,3))

     ! assume first structure is the equilibrium one
     open (unit=11,file="POSITIONs",status="old",err=99)
     do i=1,nAtoms
        read (11,*) eq(i,:)
     end do

     ! group atoms into close packed layers according to z coordinate
     layers(1) = 1
     n = 1
     do i=2,nAtoms
       alreadyIn = .FALSE.
       do j=1,i-1
         if (abs(eq(j,3)-eq(i,3))<delta) then
           alreadyIn = .TRUE.
           exit
         endif
       enddo
       if (alreadyIn) then
         layers(i) = layers(j)
       else
         n = n+1
         layers(i) = n
       endif
     enddo

     allocate(redLayer(n,3),nAtomsPerLayer(n),zLayer(n),newOrder(n),wLayerA(n),wLayerB(n),wLayerC(n))
     wLayerA = 1.5 ! weight factor to ensure continuity of A layers
     wLayerB = 1.5 ! weight factor to ensure continuity of B layers
     wLayerC = 1.5 ! weight factor to ensure continuity of C layers

     ! resort layers based on z
     do i=1,nAtoms
       zLayer(layers(i)) = eq(i,3)
     enddo
     do i=1,n
       newOrder(i) = 1
       do j=1,n
         if (zLayer(i)>zLayer(j)+delta) newOrder(i) = newOrder(i)+1
       enddo
     enddo
     do i=1,nAtoms
       layers(i) = newOrder(layers(i))
     enddo

     ! calculate atoms per layer
     nAtomsPerLayer = 0
     do i=1,nAtoms
       nAtomsPerLayer(layers(i)) = nAtomsPerLayer(layers(i)) + 1
       zLayer(layers(i)) = eq(i,3)
     enddo

     ! find primitive cell, assume it is a fraction of the first two lattice vectors
     rx=1.
     ry=1.
     do i=1,nAtoms
       if (layers(i)==1) then
         red = MATMUL(inv,eq(i,:))
         if ((abs(red(2))<delta.OR.abs(red(2)-1)<delta).AND.abs(red(1)>delta).AND.red(1)<rx) then
           rx = red(1)
           primCell(1,:) = rx*cell(1,:)
         endif
         if ((abs(red(1))<delta.OR.abs(red(1)-1)<delta).AND.abs(red(2)>delta).AND.red(2)<ry) then
           ry = red(2)
           primCell(2,:) = ry*cell(2,:)
         endif
       endif
     enddo
     primCell(3,:) = (/ 0. , 0. , 1. /)
     call invert (primCell,primInv)

     ! determine mapping into primitive cell for each atom
     do i=1,nAtoms
       mapping(i,:) = nint(MATMUL(primInv,eq(i,:)))
     enddo

     ! read in all MD coordinates
     do j=1,nAll-1
       do i=1,nAtoms
         read (11,*) coords(j,i,:)
       end do
     end do
     close(11)

     ! here the master loop over all MD snapshots
     open(unit=30,file="stacking",status="replace",err=100)
     if (verbose) then
       do i=1,n
         write (str,'(I)') i
         str = adjustl(str)
         open (unit=12+i,file="layer_order_"//trim(str),status="replace",err=100)
         open (unit=32+i,file="layer_POSITIONs_"//trim(str),status="replace",err=100)
       enddo
     endif
     do j=1,nAll-1
       ! calculate reduced coords w.r.t. primCell and average to a single coord per layer (for better statistics)
       redLayer = 0.
       do i=1,nAtoms
         red = MATMUL(primInv,coords(j,i,:))
         red(:) = red(:) - mapping(i,:)
         redLayer(layers(i),:) = redLayer(layers(i),:) + red(:)/nAtomsPerLayer(layers(i))
       end do

       if (verbose) then
         do i=1,n
           write (12+i,'(I,2F)') j,(redLayer(i,1:2)-redLayer(1,1:2))+i
           do k=1,nAtoms
             if (layers(k)==i) write (32+i,*) coords(j,k,:)
           enddo
         enddo
       endif

       ! now get the stacking sequence, first layer is always A and second B
       ! NOTE: the factor 1./6. is hard coded below because it can be uniquely
       !       used to separate layers (coordinates) in fcc, hcp, dhcp which are:
       !       0   0   .
       !       1/3 2/3 .
       !       2/3 1/3 .
       ! 
       stacking = "AB"
       Cpresent = .FALSE.
       do i=3,n
         diffA = redLayer(i,1:2)-redLayer(1,1:2)
         diffB = redLayer(i,1:2)-redLayer(2,1:2)
         if (Cpresent) diffC = redLayer(i,1:2)-redLayer(Ci,1:2)

         ! check whether layer i is in the range of layer A
         if     ((abs(diffA(1))<wLayerA(i)*1./6..OR.abs(diffA(1)+1)<wLayerA(i)*1./6..OR.abs(diffA(1)-1)<wLayerA(i)*1./6.).AND. &
                 (abs(diffA(2))<wLayerA(i)*1./6..OR.abs(diffA(2)+1)<wLayerA(i)*1./6..OR.abs(diffA(2)-1)<wLayerA(i)*1./6.)) then
           stacking = trim(stacking)//"A"
           if (wLayerA(i)<wMax) wLayerA(i) = wLayerA(i) + enhancement
           if (wLayerB(i)>wMin) wLayerB(i) = wLayerB(i) - enhancement
           if (wLayerC(i)>wMin) wLayerC(i) = wLayerC(i) - enhancement

         ! check whether layer i is in the range of layer B
         elseif ((abs(diffB(1))<wLayerB(i)*1./6..OR.abs(diffB(1)+1)<wLayerB(i)*1./6..OR.abs(diffB(1)-1)<wLayerB(i)*1./6.).AND. &
                 (abs(diffB(2))<wLayerB(i)*1./6..OR.abs(diffB(2)+1)<wLayerB(i)*1./6..OR.abs(diffB(2)-1)<wLayerB(i)*1./6.)) then
           stacking = trim(stacking)//"B"
           if (wLayerA(i)>wMin) wLayerA(i) = wLayerA(i) - enhancement
           if (wLayerB(i)<wMax) wLayerB(i) = wLayerB(i) + enhancement
           if (wLayerC(i)>wMin) wLayerC(i) = wLayerC(i) - enhancement

         ! if we got here, layer i is neither A nor B and if C is not yet present we assume that layer i is the new C layer
         elseif (.NOT.Cpresent) then
           stacking = trim(stacking)//"C"
           if (wLayerA(i)>wMin) wLayerA(i) = wLayerA(i) - enhancement
           if (wLayerB(i)>wMin) wLayerB(i) = wLayerB(i) - enhancement
           if (wLayerC(i)<wMax) wLayerC(i) = wLayerC(i) + enhancement
           Cpresent = .TRUE.
           Ci = i

         ! if C is present then check whether layer i is in the range of layer C
         elseif ((abs(diffC(1))<wLayerC(i)*1./6..OR.abs(diffC(1)+1)<wLayerC(i)*1./6..OR.abs(diffC(1)-1)<wLayerC(i)*1./6.).AND. &
                 (abs(diffC(2))<wLayerC(i)*1./6..OR.abs(diffC(2)+1)<wLayerC(i)*1./6..OR.abs(diffC(2)-1)<wLayerC(i)*1./6.)) then
           stacking = trim(stacking)//"C"
           if (wLayerA(i)>wMin) wLayerA(i) = wLayerA(i) - enhancement
           if (wLayerB(i)>wMin) wLayerB(i) = wLayerB(i) - enhancement
           if (wLayerC(i)<wMax) wLayerC(i) = wLayerC(i) + enhancement

         ! if we got here layer i could not be recognized as A,B,C
         else
           stacking = trim(stacking)//"?"
         endif
       enddo

       ! finally find the corresponding phase
       phase="0  unknown_"  ! the underscores are for easy plotting with xmgrace which would otherwise complain because of two columns with letters
       if (index("ABCABCABCABCABCABCABCABCABCABCABCABC",trim(stacking))==1) phase="1  fcc_____" ! fcc   ABC     ccc     (we know AB is in front)
       if (index("ABABABABABABABABABABABABABABABABABAB",trim(stacking))==1) phase="2  hcp_____" ! hcp   AB      hh    
       if (index("ABACABACABACABACABACABACABACABACABAC",trim(stacking))==1) phase="3  dhcp____" ! dhcp  ABAC    chch    (for dhcp two possibilities with AB in front (see below))
       if (index("ABCBABCBABCBABCBABCBABCBABCBABCBABCB",trim(stacking))==1) phase="3  dhcp____" ! dhcp  ABCB
       if (index("ABCBACABCBACABCBACABCBACABCBACABCBAC",trim(stacking))==1) phase="4  thcp____" ! thcp  ABCBAC  cchcch  Ref: McMahan & Young, Phys. Lett. A 105, 129 (1984).
       if (index("ABACBCABACBCABACBCABACBCABACBCABACBC",trim(stacking))==1) phase="4  thcp____" ! thcp  ABACBC          (for thcp four possibilities with AB in front (see below))
       if (index("ABCACBABCACBABCACBABCACBABCACBABCACB",trim(stacking))==1) phase="4  thcp____" ! thcp  ABCACB      
       if (index("ABACBCABACBCABACBCABACBCABACBCABACBC",trim(stacking))==1) phase="4  thcp____" ! thcp  ABACBC      
       if (index("ABABACABABACABABACABABACABABACABABAC",trim(stacking))==1) phase="5  6H______" ! 6H    ABABAC  chhhch  Ref: Yang et al., J. Alloys Compd. 432, 283 (2007).
       if (index("ABACACABACACABACACABACACABACACABACAC",trim(stacking))==1) phase="5  6H______" ! 6H    ABACAC          (for 6H four possibilities with AB in front (see below))
       if (index("ABCBABABCBABABCBABABCBABABCBABABCBAB",trim(stacking))==1) phase="5  6H______" ! 6H    ABCBAB
       if (index("ABCBCBABCBCBABCBCBABCBCBABCBCBABCBCB",trim(stacking))==1) phase="5  6H______" ! 6H    ABCBCB
       ! one could implement detection of phases with longer stacking sequences (>=8 layers) but for the moment I assume this is as large as we get with our supercells
       write (30,'(I6,A2,A11,A)') j,"  ",phase,stacking

     end do
     close(30)
     do i=1,n
       close(12+i)
       close(32+i)
     enddo

     stop     

 99  continue
     write (*,*) "error: cannot open or empty atoms_volume_steps, cell, or POSITIONs file"
     stop
 100 continue
     write (*,*) "error: cannot open file for writing"
     stop
end

!  derivation of symmetrically equivalent stackings for dhcp and thcp that start with AB
!                         dhcp                                 thcp                                     6H
!  permutations of ABC    permutations    AB_start   unique    permutations        AB_start   unique    permutations        AB_start   unique
!  1. ABC                 1. ABAC ABAC    ABAC       ABAC      1. ABCBAC ABCBAC    ABCBAC     ABCBAC    1. ABABAC ABABAC    ABABAC     ABABAC
!  2. ACB                 2. ACAB ACAB    ABAC                 2. ACBCAB ACBCAB    ABACBC     ABACBC    2. ACACAB ACACAB    ABACAC     ABACAC
!  3. BAC                 3. BABC BABC    ABCB       ABCB      3. BACABC BACABC    ABCBAC               3. BABABC BABABC    ABCBAB     ABCBAB
!  4. BCA                 4. BCBA BCBA    ABCB                 4. BCACBA BCACBA    ABCACB     ABCACB    4. BCBCBA BCBCBA    ABCBCB     ABCBCB
!  5. CAB                 5. CACB CACB    x                    5. CABACB CABACB    ABACBC     ABACBC    5. CACACB CACACB    x
!  6. CBA                 6. CBCA CBCA    x                    6. CBABCA CBABCA    ABCACB               6. CBCBCA CBCBCA    x
!

