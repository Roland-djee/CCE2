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
  use write
  use math
  use system_basis
  use build_hamiltonian
  use spin_generator
  implicit none
  integer :: i,j
  double precision :: A1(3,3),B1(2,2),K1(6,6)

!> Read the system basis from input
  call read_basis

  print*,B0%x
  print*,B0%y
  print*,B0%z
  print*,B0%ampli
  print*,basis(1)%spin_type
  print*,basis(1)%spin_sp
  print*,basis(2)%spin_type  
  print*,basis(2)%spin_sp

!> Create the basis
  call create_basis(basis(1)%spin_type,basis(1)%spin_sp)
 
  print*,'vector'
  do i=1,2
     do j=1,basis(i)%spin_mt
        print*,basis(i)%vector(j)
     end do
  end do
  !stop

!> Build the diagonal free-Hamiltonian
  call build_diag(basis(1)%spin_type,basis(1)%spin_sp)

  print*,'H0'
  do i=1,20
     print*,H0_diag(i)
  end do

  deallocate(H0_diag)

  call build_hf

  call spin_matrices(0.5d0)
  
  CALL PRINT_MATRIX( 'Sz', 4, 4, Sz, 4 )
  CALL PRINT_MATRIX( 'Sp', 4, 4, Sp, 4 )
  CALL PRINT_MATRIX( 'Sm', 4, 4, Sm, 4 )

  A1=reshape((/ 1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 6.d0, 0.d0, 8.d0, 0.d0 /), shape(A1))
  CALL PRINT_MATRIX( 'A1', 3, 3, A1, 3 )
  B1=reshape((/ 5.d0, 6.d0, 7.d0, 8.d0 /), shape(B1))
  CALL PRINT_MATRIX( 'B1', 2, 2, B1, 2 )
  K1=0.d0

  call kronecker(A1,3,B1,2,K1)
  
  CALL PRINT_MATRIX( 'K1', 6, 6, K1, 6 )



  stop


!> Create the basis
  call create_basis(basis(2)%spin_type,basis(2)%spin_sp)

  print*,'vector'
 
     do j=1,basis(1)%spin_mt
        print*,basis(1)%vector(j)
     end do
 
  
!> Build the diagonal free-Hamiltonian
  call build_diag(basis(2)%spin_type,basis(2)%spin_sp)

  print*,'H0'
  do i=1,2
     print*,H0_diag(i)
  end do

end program CCE2
