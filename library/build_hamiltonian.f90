!------------------------------------------------------------------------------
! CCE2 code module
!------------------------------------------------------------------------------
!
! MODULE: build_hamiltonian
!
!> @author
!> Dr. Roland Guichard University College London
!
! DESCRIPTION: 
!> This module contains all subroutines necessary to build the Hamiltonian 
!> matrix
!
! REVISION HISTORY:
! 16-02-2015 - Initial Version
! TODO_16_02_2015 - Finish build_diag - TODO_build_hamiltonian
!------------------------------------------------------------------------------

module build_hamiltonian
  use type
  use constant
  implicit none

contains

  !---------------------------------------------------------------------------  
  !> @author 
  !> Dr. Roland Guichard University College London
  !
  ! DESCRIPTION: 
  !> Computes the diagonal part of the free Hamiltonian
  !> @brief
  !> Computes the Zeeman terms of the Hamiltonian.
  !
  ! REVISION HISTORY:
  ! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
  !
  !> @param[in]  basis   
  !> @param[out] --      
  !> @return     H0_diag
  !---------------------------------------------------------------------------  

  subroutine build_diag
    implicit none
    integer :: i,j,k
    double precision :: factor,B0

    factor = 2.d0 * pi * B0
    
    k = 0
    do i=1,basis(1)%spin_mt
       do j=1,basis(2)%spin_mt
          k = k + 1
          H0_diag(k) = 1.d0
       end do
    end do

  end subroutine build_diag

end module build_hamiltonian
