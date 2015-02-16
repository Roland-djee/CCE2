module system_basis
  use type
  implicit none

contains

subroutine create_basis
  implicit none
  ! local variables
  integer :: i,j

  do i=1,2
     do j=1,basis(i)%spin_mt
        basis(i)%vector(j) = - basis(i)%spin_mag + dble(j - 1)
     end do
  end do
  
end subroutine create_basis

end module system_basis
