program CCE2
  use type
  use read
  use system_basis
  implicit none
  integer :: i,j

  call read_basis

  print*,basis(1)%spin_type
  print*,basis(2)%spin_type  
  print*,basis(1)%spin_mag
  print*,basis(2)%spin_mag
  print*,basis(1)%spin_mt
  print*,basis(2)%spin_mt

  call create_basis

  do i=1,2
     do j=1,basis(i)%spin_mt
        print*,basis(i)%vector(j)
     end do
  end do

end program CCE2
