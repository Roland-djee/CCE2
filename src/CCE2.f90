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

  !print*,B0%x
  !print*,B0%y
  !print*,B0%z
  !print*,B0%ampli
  !print*,basis(1)%spin_type
  !print*,basis(1)%spin_sp
  !print*,basis(2)%spin_type  
  !print*,basis(2)%spin_sp

!> Create the basis vectors for basis 1
!> Warning: from here if basis(1)%spin_type = 'CS'
!> basis(2)%spin_mt and basis(2)%spin_mag only are set to an electron
  call create_basis(basis(1)%spin_type,basis(1)%spin_sp)
 
  !CALL PRINT_MATRIX( 'Basis vector 1', 10, 1, basis(1)%vector, 10 )

!> Build the diagonal free-Hamiltonian for basis 1
  call build_diag(basis(1)%spin_type,basis(1)%spin_sp)

  !CALL PRINT_MATRIX( 'H0 diag.', 20, 1, H0_diag, 20 )

!> Build the hyperfine Hamiltonian matrix for basis 1
  call build_hf(basis(1)%spin_mag,basis(1)%spin_mt,basis(1)%spin_sp,&
                basis(2)%spin_mag,basis(2)%spin_mt)
  
  !CALL PRINT_MATRIX_BLOCK( 'Hf 1x1', 1, 1, 10, 10, H_hf, 20 )
  !CALL PRINT_MATRIX_BLOCK( 'Hf 1x2', 1, 11, 10, 20, H_hf, 20 )
  !CALL PRINT_MATRIX_BLOCK( 'Hf 2x1', 11, 1, 20, 10, H_hf, 20 )
  !CALL PRINT_MATRIX_BLOCK( 'Hf 2x2', 11, 11, 20, 20, H_hf, 20 )

!> Build the full Hamiltonian for the basis 1
  allocate (HCS(size(H_hf,1),size(H_hf,1)))
  HCS = H_hf
  forall (i=1:size(HCS,1)) HCS(i,i) = HCS(i,i) + H0_diag(i)

  deallocate(H0_diag,H_hf,HCS)

!> Create the basis vectors for basis 2
!> Warning: Although the input is relative to basis 2
!> the vector is overwritten in basis(1)
  call create_basis(basis(2)%spin_type,basis(2)%spin_sp)

  CALL PRINT_MATRIX( 'Basis vector 2', 2, 1, basis(2)%vector, 2 )
  
!> Build the diagonal free-Hamiltonian for basis 2
  call build_diag(basis(2)%spin_type,basis(2)%spin_sp)

  CALL PRINT_MATRIX( 'H0 diag.', 2, 1, H0_diag, 2 )

!> Build the interaction matrix for basis 2
  call build_int(basis(2)%spin_mag,basis(2)%spin_mt,basis(2)%spin_sp,&
                basis(2)%spin_mag,basis(2)%spin_mt)

  CALL PRINT_MATRIX( 'H_int', 4, 4, H_int, 4 )

  deallocate(H0_diag,H_int)

end program CCE2
