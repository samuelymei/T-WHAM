subroutine sort( n, array ) 
  use precision_m
  implicit none
  integer(kind=4), intent(in) :: n
  real(kind=fp_kind), intent(in out) :: array(n)
  integer(kind=4) :: swapped
  real(kind=fp_kind) :: temp
  integer(kind=4) :: i, j, k
  swapped = 1
  do j = n-1, 1, -1
    swapped = 0
    do i = 1, j
      if(array(i) > array(i+1))then
        temp = array(i)
        array(i) = array(i+1)
        array(i+1) = temp
        swapped = 1
      end if
    end do
    if(swapped == 0) exit
  end do
end subroutine sort

function rmsd(n,a1,a2)
  use precision_m
  implicit none
  integer(kind=4), intent(in) :: n
  real(kind=fp_kind), intent(in) :: a1(n), a2(n)
  real(kind=fp_kind) :: rmsd
  integer(kind=4) :: i
  if(n==0)then
    rmsd = 0 
    return
  end if
  rmsd = 0.d0
  do i = 1, n
    rmsd = rmsd + (a1(i)-a2(i))**2
  end do
  rmsd = sqrt(rmsd/n)
end function rmsd

