!------------------------------------------------------------------------------
! CCE2 module
!------------------------------------------------------------------------------
!
! MODULE: system_basis
!
!> @author
!> Dr. Roland Guichard University College London
!
! DESCRIPTION: 
!> Defines the vector basis for the spin system considered.
!
! REVISION HISTORY:
! 16-02-2015 - Initial Version
! TODO_16_02_2015 - Complete main - TODO_main
!------------------------------------------------------------------------------

module system_basis
  use type
  implicit none

contains

  !---------------------------------------------------------------------------  
  !> @author 
  !> Dr. Roland Guichard University College London
  !
  ! DESCRIPTION: 
  !> Creates the vector basis.
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
