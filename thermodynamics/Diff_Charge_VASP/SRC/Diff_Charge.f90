!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Author: Liangfeng Huang
!! Address: Max-Planck-Institut fur Eisenforschung
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This code is named Diff_Charge.x, because it is 
!! especially useful for the calculation of the
!! differential charge.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Capability of this code:
!!  Calculating the weighted summation of many sets
!!  charge densities output by VASP (the CHG file).
!!+++++
!!For example:
!! diff_CHG=1.0*CHG_1+(-2.0)*CHG_2+3.0*CHG_3
!!+++++
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Execution command in Linux system:
!!
!!  Diff_Charge.x < input_filename > output_filename
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Format of the input file:
!!++++++++++++++++++++++++++
!! &input
!!  no_files=N,     ! N integer, No. of files
!!  CHG_files(1)= , ! names for the N files
!!  CHG_files(2)= ,
!!  ...
!!  CHG_files(N)= , 
!!  weight(1)= ,    ! weights for the N files
!!  weight(2)= ,    ! (positive/negative)
!!  ...
!!  weight(N)= ,
!!  no_type= ,      ! No. of atom types
!! /
!!++++++++++++++++++++++++++
!! Format of the output charge density:
!!  The same as the input CHG files
!!++++++++++++++++++++++++++
!! The input/output CHG files could be converted to 
!! the *.xsf format using v2xsf, and then could be
!! visualized using XCrysDen.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! There is an example in this code package.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       PROGRAM Differential_Charge
        IMPLICIT NONE
        character(len=32),dimension(6)::CHG_files
        real(8),dimension(6)::weight
        integer::no_files
        integer::no_type
        character(len=8),dimension(:),allocatable::atom_names
        integer,dimension(:),allocatable::atom_no
        integer::tot_at_no
        character(len=50)::char_tmp
        real(8)::alat
        real(8),dimension(3,3)::lat_vect
        character(len=20)::coord_format
        logical::fractional
        real(8),dimension(:,:),allocatable::pos_atom
        real(8),dimension(:),allocatable::charge_sum
        real(8),dimension(:,:),allocatable::charge_each
        integer,dimension(3)::grid
        integer::no_data,no_line,resid_no_col
        integer::i,j,k,l,file_ind

        namelist /input/ no_files,CHG_files,weight,no_type,fractional
 
        read(*,nml=input)

        fractional=.true.

        allocate(atom_names(no_type))
        allocate(atom_no(no_type))
        
        DO i=1,no_files

          file_ind=i*100
          char_tmp=trim(adjustl(CHG_files(i)))
          open(file_ind,file=char_tmp)
          
          read(file_ind,*)

          IF(i==1) THEN
           read(file_ind,*) alat
           read(file_ind,*) (lat_vect(j,1),j=1,3)
           read(file_ind,*) (lat_vect(j,2),j=1,3)
           read(file_ind,*) (lat_vect(j,3),j=1,3)

           lat_vect=lat_vect*alat

           write(char_tmp,*) no_type*5
!           char_tmp="'("//trim(adjustl(char_tmp))//"a"//")'"
!           read(file_ind,char_tmp) (atom_names(j),j=1,no_type)
           read(file_ind,*) (atom_names(j),j=1,no_type)
           
           do j=1,no_type
            atom_names(j)=trim(adjustl(atom_names(j)))
           enddo

           read(file_ind,*) (atom_no(j),j=1,no_type)
 
           tot_at_no=0
           do j=1,no_type
            tot_at_no=tot_at_no+atom_no(j)
           enddo           

           read(file_ind,*) coord_format
           
           allocate(pos_atom(3,tot_at_no))

           do j=1,tot_at_no
            read(file_ind,*) (pos_atom(k,j),k=1,3)
           enddo

!           if(fractional==.true.) then
!            pos_atom(1,:)=pos_atom(1,:)*sqrt(lat_vect(1,1)**2+ &
!                                  lat_vect(2,1)**2+lat_vect(3,1)**2)
!            pos_atom(2,:)=pos_atom(2,:)*sqrt(lat_vect(1,2)**2+ &
!                                  lat_vect(2,2)**2+lat_vect(3,2)**2)
!            pos_atom(3,:)=pos_atom(3,:)*sqrt(lat_vect(1,3)**2+ &
!                                  lat_vect(2,3)**2+lat_vect(3,3)**2)
!           endif

           read(file_ind,*)
           read(file_ind,*) (grid(j),j=1,3)

           no_data=grid(1)*grid(2)*grid(3)
           no_line=no_data/10
           resid_no_col=no_data-no_line*10
       
           allocate(charge_sum(no_data))
           allocate(charge_each(no_data,no_files))

           charge_sum=0.0
           charge_each=0.0

          ELSE 
           do l=1,7
            read(file_ind,*)
           enddo
           do l=1,tot_at_no+2
            read(file_ind,*)
           enddo
 
          ENDIF

          Do j=0,no_line-1
           read(file_ind,*) (charge_each(l,i),l=j*10+1,j*10+10)
          EndDo

          If(resid_no_col.gt.0) then
          read(file_ind,*) (charge_each(l,i),l=j*10+1,j*10+resid_no_col)
          Endif

          close(file_ind)
        ENDDO
        
        do i=1,no_files
          charge_sum=charge_sum+weight(i)*charge_each(:,i)
        enddo

        deallocate(charge_each)

        open(11,file="diff_CHG")

        write(11,*) "Differential Charges"
        write(11,'(1x,f18.13)') alat
        do i=1,3
         write(11,'(3(1x,f12.6))') (lat_vect(j,i)/alat,j=1,3)
        enddo

!        write(char_tmp,*) 5
!        char_tmp="'("//"a"//trim(adjustl(char_tmp))//")'"
        write(11,*) (atom_names(i),i=1,no_type)
        write(11,*) (atom_no(i),i=1,no_type)
        write(11,*) trim(adjustl(coord_format))

        do i=1,tot_at_no
         write(11,'(3f10.6)') (pos_atom(j,i),j=1,3)
        enddo
        write(11,*)
        write(11,'(3i5)') (grid(i),i=1,3)

        Do i=0,no_line-1
         write(11,'(10(1x,f10.3,1x))') (charge_sum(l),l=i*10+1,i*10+10)
        EndDo

        If(resid_no_col.gt.0) then
         write(char_tmp,*) resid_no_col
         char_tmp="("//trim(adjustl(char_tmp))//"(1x,f10.3,1x)"//")"
         write(11,char_tmp) (charge_sum(l),l=i*10+1,i*10+resid_no_col)
        Endif



       END PROGRAM Differential_Charge
