! find neighbors in the supercell, 
! e.g. 1st, 2nd, 3rd ... nearest neighbors
! Format of the input:
!**********************************
!  e1   e2   e3   !fract. coordin. of the reference atom
! len_scale       !length scale
! x_a  y_a  z_a   !lattice vectors a, b, & c in Cartesian
! x_b  y_b  z_b
! x_c  y_c  z_c
! N_atom          !atomic number
! x_1 y_1 z_1     !fract. coordin. of all the atoms
! x_2 y_2 z_2
! ......
! x_N y_N z_N
!**********************************
!Execution command:
! find_neighbor.x < inputfilename
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      PROGRAM find_neighb
        implicit none
        integer::N_atom
        real,dimension(3)::ref_pos,pos_tmp
        real,dimension(3,3)::latt_const
        real,dimension(:,:),allocatable::pos   !dimension(3,N_atom)
        real,dimension(:),allocatable::dist_nb !dimension (N_atom)
        real::len_scale
        real::dist_tmp,dist_cal       
        integer::i,j,k,l,m
        
        read(*,*) (ref_pos(i),i=1,3)
        read(*,*) len_scale
        do i=1,3
          read(*,*) (latt_const(j,i),j=1,3)
        enddo

!        print*,latt_const

        read(*,*) N_atom

        latt_const=latt_const*len_scale
!        ref_pos=ref_pos*len_scale

        allocate(pos(3,N_atom))
        allocate(dist_nb(N_atom))

        write(*,*) "#Index    Distance    Fract. Coord."
        
        do i=1,N_atom
          read(*,*) (pos(j,i),j=1,3)
          dist_nb(i)=dist_cal(pos(1:3,i),ref_pos,latt_const)
        enddo

        do i=1,N_atom
          do j=i+1, N_atom
            if(dist_nb(i).gt.dist_nb(j)) then

              dist_tmp=dist_nb(i)
              dist_nb(i)=dist_nb(j)
              dist_nb(j)=dist_tmp
 
              pos_tmp=pos(1:3,i)
              pos(1:3,i)=pos(1:3,j)
              pos(1:3,j)=pos_tmp

            endif
          enddo
        enddo

        do i=1,N_atom
          write(*,'(i3,4f12.6)') i, dist_nb(i), (pos(j,i),j=1,3)
        enddo

      END PROGRAM find_neighb



      function dist_cal(pos1,pos2,latt)
        implicit none
        real,dimension(3)::pos1,pos2
        real,dimension(3,3)::latt
        real,dimension(3)::pos1_xyz,pos2_xyz,pos2_tmp
        real::dist_cal
        real::dist_tmp
        integer,dimension(3)::l_ind
        integer::l1,l2,l3
        integer::i,j

        pos1_xyz=0.0
        pos2_xyz=0.0
        dist_cal=0.0
        dist_tmp=0.0

        do i=1,3
          do j=1,3
            pos1_xyz(i)=pos1_xyz(i)+pos1(j)*latt(i,j)
            pos2_xyz(i)=pos2_xyz(i)+pos2(j)*latt(i,j)
          enddo
        enddo

        dist_cal=-1.0

        do l3=-1,1
         do l2=-1,1
          do l1=-1,1
            l_ind(1)=l1
            l_ind(2)=l2
            l_ind(3)=l3

            pos2_tmp=pos2_xyz

            do i=1,3

             do j=1,3
              pos2_tmp(i)=pos2_tmp(i)+l_ind(j)*latt(i,j)
             enddo

             dist_tmp=dist_tmp+(pos1_xyz(i)-pos2_tmp(i))**2.0

            enddo

           if((dist_cal.lt.0.0).or.(dist_cal.gt.dist_tmp)) then
             dist_cal=dist_tmp
           endif

           dist_tmp=0.0

          enddo
         enddo
        enddo

        dist_cal=dist_cal**0.5
      end function

