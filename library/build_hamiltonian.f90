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
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------

module build_hamiltonian
  implicit none

contains

  !---------------------------------------------------------------------------  
  !> @author 
  !> Dr. Roland Guichard University College London
  !
  ! DESCRIPTION: 
  !> Computes the diagonal part of the free Hamiltonian
  !> @brief
  !> Flow method (rate of change of position) used by integrator.
  !> Compute \f$ \frac{d\lambda}{dt} , \frac{d\phi}{dt},  \frac{dz}{dt} \f$
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

  

  end subroutine build_diag

end module build_hamiltonian
