! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing

! Made F conformant by Walt Brainerd

module qsort_c_module

implicit none
public :: QsortC
private :: Partition
type list
  real, allocatable :: A(:)
end type

contains

recursive subroutine QsortC(l)
  type(list),intent(in out), dimension(:) :: l
  type(list),allocatable,dimension(:) :: l1,l2
  integer :: iq,i,n,m

  n=size(l)
  m=size(l(1)%A)

  if(m > 1) then
     call Partition(l, iq)
     allocate(l1(n),l2(n))
     do i=1,n
       allocate(l1(i)%A(iq-1),l2(i)%A(m-iq+1))
       l1(i)%A(:)=l(i)%A(:iq-1)
       l2(i)%A(:)=l(i)%A(iq:)
     enddo
     call QsortC(l1)
     call QsortC(l2)
     do i=1,n
       l(i)%A(:iq-1)=l1(i)%A(:)
       l(i)%A(iq:)=l2(i)%A(:)
     enddo
   ! call QsortC(l(:iq-1))
   ! call QsortC(l(iq:))
  endif
end subroutine QsortC

subroutine Partition(l, marker)
  type(list), intent(in out), dimension(:) :: l
  integer, intent(out) :: marker
  integer :: i, j, k
  real :: temp
  real :: x      ! pivot point
  x = l(1)%A(1)
  i= 0
  j= size(l(1)%A) + 1

  do
     j = j-1
     do
        if (l(1)%A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (l(1)%A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        do k=1,size(l)
          temp = l(k)%A(i)
          l(k)%A(i) = l(k)%A(j)
          l(k)%A(j) = temp
        enddo
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition

end module qsort_c_module

