       MODULE parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Input files needed for each Job:
!!(1) EIGENVAL
!!(2) EIGENVAL+PROCAR
!!(2+3) EIGENVAL+PROCAR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
       real(8),dimension(:,:,:,:),allocatable::proj_mat ! projection matrix. 
                                              ! index: (orbital:atom:band:k)
       real(8)::emin,emax
       integer::i,j,k,l
   
       namelist /input/ band_proj, band_sym, n_atom, n_orbit, &
                        emin, emax

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Setting default values
!       band_proj=.false.
!       band_sym=.false.
!       n_atom=-1
!       n_orbit=-1
!       emin=-10001
!       emax=-10001
!       out_dir=""
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      END MODULE parameters
