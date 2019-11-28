MODULE utilities
CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE get_k_path:
!!   Mapping the k points into a 1-D k path.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE get_k_path(k_path,k_point,n_k)
  IMPLICIT NONE
  integer::n_k
  real(8)::k_tmp,k_scal,dk
  real(8),dimension(n_k)::k_path
  real(8),dimension(3,n_k)::k_point
  real(8),dimension(3,2)::k_HSP ! for finding the high-symmetry points
  integer::i,j

  k_tmp=0.0
  k_scal=0.0
  dk=0.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  DO i=1,n_k

   do j=1,3
     k_tmp=k_tmp+k_point(j,i)**2
   enddo

   if(abs(k_tmp).lt.1.0d-3) write(999,'(3f14.7)') k_point(:,i)
   
   IF(i.eq.1) THEN
     k_path(i)=0.0
   ELSEIF(i.gt.1) THEN

     do j=1,3
       dk=dk+(k_point(j,i)-k_point(j,i-1))**2
     enddo
     dk=sqrt(dk)

     k_path(i)=k_path(i-1)+dk
     dk=0.0

     do j=1,2
      k_HSP(:,j)=k_point(:,i+j-1)-k_point(:,i+j-2)
     enddo
     
     do j=1,3
      k_scal=k_scal+k_HSP(j,1)*k_HSP(j,2)
     enddo

     k_scal=k_scal/sqrt(k_HSP(1,1)**2+k_HSP(2,1)**2+k_HSP(3,1)**2) &
                  /sqrt(k_HSP(1,2)**2+k_HSP(2,2)**2+k_HSP(3,2)**2)

   ENDIF   

   k_tmp=0.00
   k_scal=0.00

  ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

END SUBROUTINE get_k_path

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!SUBROUTINE get_proj_mat:
!!  Reading data from PROCAR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_proj_mat(proj_mat,phase_R,phase_I,n_orbit,n_atom,n_band,n_k)
  IMPLICIT NONE
  integer::n_orbit,n_atom,n_band,n_k
  real(8),dimension(n_orbit,n_atom,n_band,n_k)::proj_mat
  real(8),dimension(n_orbit,n_atom,n_band,n_k)::phase_R,phase_I
  real(8)::val_tmp
  integer::no_tmp
  integer::i,j,k,l
  
  open(33,file="PROCAR")

  read(33,*)
  read(33,*)

  DO i=1,n_k
   read(33,*)
   read(33,*)

   Do j=1,n_band
    read(33,*)
    read(33,*)
    read(33,*)
    read(33,*)

    do k=1,n_atom
      read(33,*) no_tmp,(proj_mat(l,k,j,i),l=1,n_orbit),val_tmp
    enddo

    read(33,*)
    read(33,*)

    do k=1,n_atom
      read(33,*) no_tmp,(phase_R(l,k,j,i),l=1,n_orbit)
      read(33,*) no_tmp,(phase_I(l,k,j,i),l=1,n_orbit)
    enddo

   Enddo

   read(33,*)
  ENDDO


END SUBROUTINE get_proj_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!SUBROUTINE symmetrize_band:
!!  Sorting the band dispersion according to their symmetry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE symmetrize_band(eigen_k,proj_mat,phase_R,phase_I,n_orbit,n_atom,n_band,n_k,dE_0)
  IMPLICIT NONE
  integer::n_orbit,n_atom,n_band,n_k
  real(8),dimension(n_band,n_k)::eigen_k
  real(8),dimension(n_orbit,n_atom,n_band,n_k)::proj_mat
  real(8),dimension(n_orbit,n_atom,n_band,n_k)::phase_R,phase_I
  complex(8),dimension(:,:,:,:),allocatable::proj_mat_full
  complex(8),dimension(:,:),allocatable::proj_tmp_full
  real(8),dimension(:,:),allocatable::proj_tmp
  complex(8)::d_proj_tmp
  real(8)::d_proj
  real(8)::val_tmp
  real(8)::len1,len2
  integer::no_tmp,band_index
  real(8)::delta_E,dE_0   ! Energy criterion
  integer::i,j,k,l,s,m,n

  allocate(proj_mat_full(n_orbit,n_atom,n_band,n_k))  
  allocate(proj_tmp_full(n_orbit,n_atom))
  allocate(proj_tmp(n_orbit,n_atom))


  do i=1,n_k
   do j=1,n_band
    do k=1,n_atom
     do l=1,n_orbit
      proj_mat_full(l,k,j,i)=proj_mat(l,k,j,i)*cmplx(phase_R(l,k,j,i),phase_I(l,k,j,i))
     enddo
    enddo
   enddo
  enddo


  DO i=2,n_k

   Do j=1,n_band

    d_proj=-1.0
    band_index=j

    do k=j,n_band

      d_proj_tmp=(0.0,0.0)

      do l=1,n_atom
       do s=1,n_orbit
        d_proj_tmp=d_proj_tmp+proj_mat_full(s,l,j,i-1)*conjg(proj_mat_full(s,l,k,i))
       enddo
      enddo

      delta_E=abs(eigen_k(j,i-1)-eigen_k(k,i))

      if(abs(d_proj_tmp)>d_proj.and.delta_E<dE_0) then
       d_proj=abs(d_proj_tmp)
       band_index=k
      endif

    enddo

    if(band_index.ne.j) then

     do m=i,n_k
      val_tmp=eigen_k(band_index,m)
      eigen_k(band_index,m)=eigen_k(j,m)
      eigen_k(j,m)=val_tmp

      proj_tmp=proj_mat(:,:,band_index,m)
      proj_mat(:,:,band_index,m)=proj_mat(:,:,j,m)
      proj_mat(:,:,j,m)=proj_tmp

      proj_tmp_full=proj_mat_full(:,:,band_index,m)
      proj_mat_full(:,:,band_index,m)=proj_mat_full(:,:,j,m)
      proj_mat_full(:,:,j,m)=proj_tmp_full
     enddo
    endif

   Enddo

  ENDDO 


  DO i=n_k-1,1,-1

   Do j=1,n_band

    d_proj=-1.0
    band_index=j

    do k=j,n_band

      d_proj_tmp=(0.0,0.0)

      do l=1,n_atom
       do s=1,n_orbit
        d_proj_tmp=d_proj_tmp+proj_mat_full(s,l,j,i+1)*conjg(proj_mat_full(s,l,k,i))
       enddo
      enddo

      delta_E=abs(eigen_k(j,i+1)-eigen_k(k,i))

      if(abs(d_proj_tmp)>d_proj.and.delta_E<dE_0) then
       d_proj=abs(d_proj_tmp)
       band_index=k
      endif

    enddo

    if(band_index.ne.j) then

      val_tmp=eigen_k(band_index,i)
      eigen_k(band_index,i)=eigen_k(j,i)
      eigen_k(j,i)=val_tmp

      proj_tmp=proj_mat(:,:,band_index,i)
      proj_mat(:,:,band_index,i)=proj_mat(:,:,j,i)
      proj_mat(:,:,j,i)=proj_tmp

      proj_tmp_full=proj_mat_full(:,:,band_index,i)
      proj_mat_full(:,:,band_index,i)=proj_mat_full(:,:,j,i)
      proj_mat_full(:,:,j,i)=proj_tmp_full
    endif

   Enddo

  ENDDO 


  deallocate(proj_mat_full)
  deallocate(proj_tmp_full)
  deallocate(proj_tmp)

END SUBROUTINE symmetrize_band
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




END MODULE utilities



