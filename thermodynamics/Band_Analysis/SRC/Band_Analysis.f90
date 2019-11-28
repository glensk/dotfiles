!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!This code is a post-processer for the data !!
!!calculated with VASP code                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  Author: Liangfeng Huang
!!  E-mail: l.huang@mpie.de
!! Address: Computational Materials Design Dept.
!!          Max-Planck-Institut fur Eisenforschung GmbH
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Capability of this code:
!!(1) band dispersions (without symmetry)
!!(2) sorting the electronic bands according to 
!!    their symmetry
!!(3) projection of each electronic state to 
!!    each atom and orbital
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Input files needed for each Job:
!!(1) EIGENVAL
!!(2) EIGENVAL+PROCAR(phase included)
!!(2+3) EIGENVAL+PROCAR(phase included)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Inputfile format:
!!************************
!! &input
!!  band_proj =  !.true. or .false. (default)
!!  band_sym  =  !.true. or .false. (default)
!!  n_atom    =  ! Number of atoms  
!!  n_orbit   =  ! Number orbitals in POSCAR
!!  E_fermi   =  ! Fermi Energy level (default: 0.0 eV)
!!  dE        =  ! Band swapping window (default: 0.3 eV)
!! /
!!************************
!! NOTE: two bands with energy difference larger than dE will not 
!!       swapped. When your band dispersions are not well symmetrized,
!!       you'd better increase the k-point number and/or change the 
!!       value of dE (try some values larger/less than the default 0.3 eV).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Execution command:
!! Band_Analysis.x < inputfile > outputfile 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      PROGRAM Bands_analysis
       
       USE utilities

       IMPLICIT NONE              

       logical::band_proj,band_sym  ! switches for the band projection and sorting
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Default:  band_proj = .false.
!!           band_sym  = .false. if band_proj = .false.
!!                       .true.  if band_proj = .true.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       character(len=80)::char_tmp
       integer::num_tmp
       real(8)::value_tmp
       integer::n_band,n_k,n_atom,n_orbit,n_file
       real(8),dimension(:,:),allocatable::k_point
       real(8),dimension(:),allocatable::k_path       ! one-dimension k path
       real(8),dimension(:,:),allocatable::eigen_k    ! eigenvalues at each k point
       real(8),dimension(:,:,:,:),allocatable::proj_mat  ! projection matrix
       real(8),dimension(:,:,:,:),allocatable::phase_R,phase_I ! projection phase
                                                      ! index: (orbital:atom:band:k)
       real(8)::dE,E_fermi
       real(8)::emin,emax
       integer::i,j,k,l

       namelist /input/ band_proj, band_sym, n_atom, n_orbit, &
                        dE,E_fermi,emin,emax


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      print*, "Setting default values"
       band_proj=.false.
       band_sym=.false.
       n_atom=-1
       n_orbit=-1
       dE=0.3
       E_fermi=0.0
       emin=-10001.0
       emax=-10001.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       read(*,nml=input)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Opening prerequisite files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       open(100,file="output_info")
       open(11,file='EIGENVAL')
       open(22,file='bands.dat')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       if(band_proj.eqv..true.) band_sym=.true.
       if(n_orbit.lt.1) then
         write(*,*) 'ERROR: Orbital Number is missing !!!'
         STOP
       endif
       if(n_atom.lt.0) then
         write(*,*) 'ERROR: Atom Number is missing !!!'
         STOP
       endif
       

       do i=1,5
        read(11,*)
       enddo
       read(11,"(3i5)") num_tmp,n_k,n_band


       allocate(k_point(3,n_k))
       allocate(k_path(n_k))
       allocate(eigen_k(n_band,n_k))

       DA:DO i=1,n_k

        read(11,*)
        read(11,*) (k_point(l,i),l=1,3),value_tmp

        DB:Do j=1,n_band
         read(11,*) num_tmp,eigen_k(j,i)
         eigen_k(j,i)=eigen_k(j,i)-E_fermi
        Enddo DB

       ENDDO DA


       call get_k_path(k_path,k_point,n_k)


       IA:IF(band_sym.eqv..true.) THEN

        allocate(proj_mat(n_orbit,n_atom,n_band,n_k))
        allocate(phase_R(n_orbit,n_atom,n_band,n_k))
        allocate(phase_I(n_orbit,n_atom,n_band,n_k))
                
        call get_proj_mat(proj_mat,phase_R,phase_I,  &
                          n_orbit,n_atom,n_band,n_k)


        call symmetrize_band(eigen_k,proj_mat,phase_R,phase_I,  &
                              n_orbit,n_atom,n_band,n_k,dE)

        deallocate(phase_R)
        deallocate(phase_I)

       ENDIF IA


       write(char_tmp,*) n_band
       char_tmp=adjustl(char_tmp)
       char_tmp="(f10.4,"//trim(char_tmp)//"f10.4"//")"

       DO i=1,n_k
         write(22,char_tmp) k_path(i),eigen_k(:,i)
       ENDDO


       IB:If(band_proj.eqv..true.) then

        n_file=n_atom*n_orbit

        do i=1,n_atom
         do j=1,n_orbit
          num_tmp=100*i+j
          write(char_tmp,*) num_tmp
          char_tmp="proj_"//trim(adjustl(char_tmp))
          open(num_tmp,file=char_tmp)
         enddo
        enddo
         

        Do i=1,n_band

         IF(mod(i,2).eq.1) then

          do j=1,n_k
            do k=1,n_atom
              do l=1,n_orbit
              num_tmp=100*k+l
              write(num_tmp,*) k_path(j),eigen_k(i,j),proj_mat(l,k,i,j)
              enddo
            enddo
          enddo

         ELSEIF(mod(i,2).eq.0) then
          do j=n_k,1,-1
            do k=1,n_atom
              do l=1,n_orbit
               num_tmp=100*k+l
               write(num_tmp,*) k_path(j),eigen_k(i,j),proj_mat(l,k,i,j)
              enddo
            enddo
          enddo           

         ENDIF 

       Enddo

       Endif IB

       deallocate(proj_mat)
       do i=1,n_atom
        do j=1,n_orbit
         num_tmp=100*i+j
         close(num_tmp)
        enddo
       enddo

      END PROGRAM Bands_analysis 



