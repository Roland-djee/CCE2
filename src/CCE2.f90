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
  double precision :: Sz(19,10),Sp(10,10),Sm(10,10)

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
  call create_basis

  print*,basis(1)%spin_mag
  print*,basis(1)%spin_mt
  print*,basis(2)%spin_mag  
  print*,basis(2)%spin_mt

!> Build the diagonal matrix elements of the Hamiltonian
  call build_diag

  CALL PRINT_MATRIX( 'H0_diag', tot_basis_mt, 1, H0_diag, tot_basis_mt )

!> Build the off-diagonal matrix elements of the Hamiltonian
  
  call build_interaction_matrices

  XJ1 = 1.315044791543825d-02 * 1.d6
  XJ2 = 1.690283813464443d-02 * 1.d6
  C12 = -3.863519669068800d-04* 1.d6  

  allocate (H_tot(16,16))
  H_tot = 0.d0
  
  !H_tot = C12 * H_int good !!!

  !H_tot = C12 * H_int + (XJ1 - XJ2) * H_shf + H_cs
  !H_tot = C12 * H_int
  !H_tot = XJ1 * H_shf1 + XJ2 * H_shf2 + H_cs
  !H_tot = H_shf1 + H_shf2
  H_tot = (XJ1 + XJ2) * H_shf


  CALL PRINT_MATRIX_BLOCK( 'H_tot 1x1', 1, 1, 8, 8, H_tot, 16 )
  CALL PRINT_MATRIX_BLOCK( 'H_tot 1x2', 1, 9, 8, 16, H_tot, 16 )
  CALL PRINT_MATRIX_BLOCK( 'H_tot 2x1', 9, 1, 16, 8, H_tot, 16 )
  CALL PRINT_MATRIX_BLOCK( 'H_tot 2x2', 9, 9, 16, 16, H_tot, 16 )

end program CCE2
