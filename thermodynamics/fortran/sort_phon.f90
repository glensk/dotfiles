!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Author: Liang Feng Huang                  !! 
!!   From: Institute of Solid State Physics, !!
!!         Chinese Academy of Sciences.      !!
!!   Homepage: http://www.theory.issp.ac.cn  !!
!!   Address: Shushanhu Road 350, Hefei,     !!
!!            Anhui, China. 230031.          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!What can this code do:                      !! 
!! (1) Sort phonon frequencies by Magnitude;  !!
!! (2) Sort phonon frequencies by Band index; !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Input Parameters :                                      !!
!!  &input                                                 !!
!!   atom_No= xx (an integer),    ! required               !! 
!!   q_No= yy (an integer),       ! required               !!
!!   file_modes='zz' (a file name)! default:"matdyn.modes" !!
!!  /                                                      !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! For more, please refer to                                         !!
!!    (1) the Tutorial (in pdf format);                              !!
!!    (2) the Example (graphene phonon dispersions)                  !!
!!   which are both included in this package.                        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      program Phonon_sort
       implicit none
       integer::atom_No,q_No
       logical::partial_ph
       character(len=32)::file_modes
       character(len=32)::char_tmp,char_tmp1
       character(len=32)::char_transform
       real(4)::max_vec_par
       complex(4)::vec_par
       integer::i,j,k,l,m,max_vec_ind,HSP_No
       real(4)::q_sum,q_proj
       character(len=32)::head_char
       integer::head_No
       real(4)::num_tmp,egvec_x,egvec_y,egvec_z
       real(4),dimension(:,:),allocatable::q_path
       real(4),dimension(:),allocatable::q_path_scal
       real(4),dimension(:,:),allocatable::omega
       complex(4),dimension(:,:,:,:),allocatable::egvec
       complex(4),dimension(:,:),allocatable::vec_tmp
       integer::stat

       namelist /input/ atom_No,q_No,file_modes,partial_ph

       atom_No=0
       q_No=0
       partial_ph=.false.

       vec_par=cmplx(0.0,0.0)
       vec_tmp=cmplx(0.0,0.0)
       write(0,*) "Read input file......"
       read(*,nml=input)
       if(atom_No.lt.1 .AND. q_No.lt.1) then
         print*,"ERROR: atom_No or q_No .lt. 1 !!!"
         stop
       endif

       write(*,*) "Total No. of Phononic Bands"
       write(*,'(i10)') 3*atom_No
       write(*,'(a35)') "***********************************"

       write(char_tmp1,*) 3*atom_No
       char_tmp1=adjustl(char_tmp1)
       char_tmp="(f12.4,"//trim(char_tmp1)//"f12.4)"

!!!!!!!!!!!!!!!!!!!!!!!!!!!  For partial phonons  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       IF(partial_ph==.true.) THEN                    
        do i=1,atom_No
         write(char_tmp1,*) i
         char_tmp1=adjustl(char_tmp1)
         m=100+i
         open(m,file="partial_atom_"//trim(char_tmp1)//".dat", &
                                             form='formatted')
         write(m,*) "# q_path, omega_i(q), |e_s^i|^2 (i=1,2,3), &
                    sum{|e_s^i|^2 (i=1,2,3}"
        enddo
       ENDIF                                           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

       allocate(q_path(1:3,1:q_No))
       allocate(q_path_scal(1:q_No))
       allocate(omega(1:3*atom_No,1:q_No))
       allocate(egvec(1:3,1:atom_No,1:3*atom_No,1:q_No))
       allocate(vec_tmp(1:3,1:atom_No))

       write(0,*) "Open file including vibrational modes......"
       open(1,file=trim(file_modes),status='old',form='formatted', &
            iostat=stat)
       if(stat/=0) then
        print*, "Open modes file ERROR!!!"
        stop
       endif
!       write(0,*) "Creat files including phonon dispersions......"
       open(2,file='phonon_magnitude.dat',form='formatted',iostat=stat)
       if(stat/=0) then
        print*, "Cannot creat file!!!"
        stop
       endif
       write(2,*) "# q_path | omega_i (i=1,2...)"

       open(3,file='phonon_band.dat',form='formatted',iostat=stat)
       if(stat/=0) then
        print*, "Cannot creat file!!!"
        stop
       endif
       write(3,*) "# q_path | omega_i (i=1,2...)"
       write(0,*) "......DONE."

       write(0,*) "Reading data of vibrational modes......"

       DO i=1,q_No
        read(1,*)
        read(1,*)
        read(1,'(a5,3f12.4)') head_char,(q_path(j,i),j=1,3)
!                                 print*,head_char,(q_path(j,i),j=1,3)  !!!TEST
        read(1,*)
        do j=1,3*atom_No
!         read(1,'(a11,i2,a6,f12.6,a11,f12.6)') head_char,head_No, &
!                   head_char,num_tmp,head_char,omega(j,i)
         read(1,'(a19,f12.6,a11,f12.6)') head_char, &
                             num_tmp,head_char,omega(j,i)
!                                  print*,omega(j,i)                 !!! TEST
         do k=1,atom_No
           read(1,'(3(a2,2f11.6))') head_char, &
            egvec(1,k,j,i),head_char,egvec(2,k,j,i),head_char, &
            egvec(3,k,j,i)
         enddo
        enddo
        read(1,*)
       ENDDO     
       write(0,*) "......DONE"
       
      HSP_No=0
      write(0,*) "Sort vibrational frequency by magnitude......"
      q_sum=0.0
 
      q_path_scal(1)=0.0
      DO i=1,q_No
       IF(i>1) THEN
        q_sum=sqrt((q_path(1,i)-q_path(1,i-1))**2.0+(q_path(2,i)-   & 
              q_path(2,i-1))**2.0+(q_path(3,i)-q_path(3,i-1))**2.0) &
              +q_sum
        q_path_scal(i)=q_sum
        if(i<=q_No-1) then
          q_proj=    &
              ((q_path(1,i)-q_path(1,i-1))*(q_path(1,i+1)-q_path(1,i)) &
              +(q_path(2,i)-q_path(2,i-1))*(q_path(2,i+1)-q_path(2,i)) &
              +(q_path(3,i)-q_path(3,i-1))*(q_path(3,i+1)-q_path(3,i)))&
              /sqrt((q_path(1,i)-q_path(1,i-1))**2.0+  &
                    (q_path(2,i)-q_path(2,i-1))**2.0+  &
                    (q_path(3,i)-q_path(3,i-1))**2.0)  &
              /sqrt((q_path(1,i+1)-q_path(1,i))**2.0+  &
                    (q_path(2,i+1)-q_path(2,i))**2.0+  &
                    (q_path(3,i+1)-q_path(3,i))**2.0) 
          if(q_proj>1.001) then
            print*,"ERROR: q_proj > 1.0"
            print*,(q_path(j,i),j=1,3)
          endif
          if((1.0-q_proj)>0.001) then
            write(*,'(a19,i3,a2)') "High symmetry Point",HSP_No+1,":"
            write(*,'(a5,3f8.5)') "q = ",(q_path(j,i),j=1,3)
            write(*,'(a10,f8.5)') "q_path = ",q_sum
            write(*,'(a35)') '***********************************'
            HSP_No=HSP_No+1
          endif
        endif
       ENDIF

        write(2,char_tmp) q_sum,(omega(j,i),j=1,3*atom_No)

      ENDDO
      write(0,*) "......DONE."

      write(*,'(a34,i3//)') "Total No. of High-Symmetry Points:",HSP_No


      write(0,*) "Sort vibrational frequency by band......"

      q_sum=0.0     

      DO i=1,q_No
        IF(i>1) THEN
          q_sum=sqrt((q_path(1,i)-q_path(1,i-1))**2.0+(q_path(2,i)-   &
                q_path(2,i-1))**2.0+(q_path(3,i)-q_path(3,i-1))**2.0) &
                +q_sum
         
         Do j=1,3*atom_No-1
          max_vec_par=0.0
          max_vec_ind=j
          do l=j,3*atom_No
            do k=1,atom_No
              do m=1,3
                vec_par=vec_par+ &
                       egvec(m,k,j,i-1)*conjg(egvec(m,k,l,i))
              enddo
            enddo
            if(abs(vec_par)>max_vec_par) then
              max_vec_par=abs(vec_par)
              max_vec_ind=l
            endif
            vec_par=(0.,0.)
          enddo

            if(max_vec_ind>j) then
              do k=1,atom_No
                do m=1,3
                  vec_tmp(m,k)=egvec(m,k,j,i)
                  egvec(m,k,j,i)=egvec(m,k,max_vec_ind,i)
                  egvec(m,k,max_vec_ind,i)=vec_tmp(m,k)
                enddo
              enddo
              num_tmp=omega(j,i)
              omega(j,i)=omega(max_vec_ind,i)
              omega(max_vec_ind,i)=num_tmp
            endif
         EndDo


        ENDIF
        write(3,char_tmp) q_sum,(omega(j,i),j=1,3*atom_No)

      ENDDO

!!!!!!!!!!!!!!!!!!!! Write data on partial phonons into files !!!!!!!!!!!!!!!!!!!
      IF(partial_ph==.true.) THEN
       DO i=1,3*atom_No
         Do j=1,q_No
           do k=1,atom_No
             m=100+k
             egvec_x=abs(egvec(1,k,i,j)*conjg(egvec(1,k,i,j)))
             egvec_y=abs(egvec(2,k,i,j)*conjg(egvec(2,k,i,j)))
             egvec_z=abs(egvec(3,k,i,j)*conjg(egvec(3,k,i,j)))
             write(m,'(f7.3,5f12.4)') q_path_scal(j), &
                   omega(i,j),egvec_x,egvec_y,egvec_z,    &
                   egvec_x+egvec_y+egvec_z
             if(j==q_No) then
              write(m,'(/)')
              write(m,'(/)')
             endif
           enddo
         Enddo    
       ENDDO
      ENDIF      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(0,*) "......DONE"
       
      end program


!      function char_transform(x) result(new_char)
!        implicit none
!        integer::x
!        integer::z,a,len_int
!        character(len=32)::new_char
!        new_char=''
!        len_int=1+floor(log10(real(x)))
!        do a=1,len_int
!         z=mod(x,10)
!         new_char=char(z+48)//new_char
!         x=x/10
!        enddo
!      end function char_transform


