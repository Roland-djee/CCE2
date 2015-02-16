!------------------------------------------------------------------------------
! CCE2 main code
!------------------------------------------------------------------------------
!
! MAIN: CCE2
!
!> @author
!> Dr. Roland Guichard University College London
!
! DESCRIPTION: 
!> This part of the CCE2 code is the main. 
!
! REVISION HISTORY:
! 16-02-2015 - Initial Version
! TODO_16_02_2015 - Complete main - TODO_main
! TODO_16_02_2015 - Deallocate basis%vector - TODO_main
!------------------------------------------------------------------------------

program CCE2
  use type
  use read
  use system_basis
  use build_hamiltonian
  implicit none
  integer :: i,j

!> Read the system basis from input
  call read_basis

  print*,basis(1)%spin_type
  print*,basis(2)%spin_type  
  print*,basis(1)%spin_mag
  print*,basis(2)%spin_mag
  print*,basis(1)%spin_mt
  print*,basis(2)%spin_mt

!> Create the basis
  call create_basis

  do i=1,2
     do j=1,basis(i)%spin_mt
        print*,basis(i)%vector(j)
     end do
  end do

!> Build the diagonal free-Hamiltonian
  call build_diag

end program CCE2
