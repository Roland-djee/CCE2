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
  use constant
  implicit none

contains

  !---------------------------------------------------------------------------  
  !> @author 
  !> Dr. Roland Guichard University College London
  !
  ! DESCRIPTION: 
  !> Sets all the parameters for the construction of the basis vectors.
  !> @brief
  !>
  !
  ! REVISION HISTORY:
  ! TODO_19_02_2015 - Finish cases for 'e' and 'n' - TODO_create_basis
  !
  !> @param[in]  basis   
  !> @param[out] --      
  !> @return     --
  !---------------------------------------------------------------------------

  subroutine create_basis(type,specie)
    implicit none
    ! local variables
    character (len=20), intent(in) :: type,specie

    select case (type)
    case ("CS")     
       !> Central spin system {donor + e-}
       !> Sets spin magnitudes
       select case (specie)
       case ("Bi")
          basis(1)%spin_mag = I_209Bi
          basis(1)%spin_mt  = int(2.d0 * basis(1)%spin_mag + 1.d0)
          basis(2)%spin_mt  = 2
          basis(2)%spin_mag = S_e
       case ("P")
          basis(1)%spin_mag = I_31P
          basis(1)%spin_mt  = int(2.d0 * basis(1)%spin_mag + 1.d0)
          basis(2)%spin_mt  = 2
          basis(2)%spin_mag = S_e
       case ("Si")
          basis(1)%spin_mag = I_29Si
          basis(1)%spin_mt  = int(2.d0 * basis(1)%spin_mag + 1.d0)
          basis(2)%spin_mt  = 2
          basis(2)%spin_mag = S_e
       case ("e")
          write(*,*)'To be finished...'
       case ("n")
          write(*,*)'To be finished...'
       end select       
       call create_vectors2
    case ("Bath")
       !> Bath spin system
       !> Sets spin magnitudes
       select case (specie)
       case ("Bi")
          basis(1)%spin_mag = I_209Bi
          basis(1)%spin_mt  = int(2.d0 * basis(1)%spin_mag + 1.d0)
       case ("P")
          basis(1)%spin_mag = I_31P
          basis(1)%spin_mt  = int(2.d0 * basis(1)%spin_mag + 1.d0)
       case ("Si")
          basis(1)%spin_mag = I_29Si
          basis(1)%spin_mt  = int(2.d0 * basis(1)%spin_mag + 1.d0)
       case ("e")
          write(*,*)'To be finished...'
       case ("n")
          write(*,*)'To be finished...'
       end select       
       basis(1)%spin_mt = int(2.d0 * basis(2)%spin_mag + 1.d0)
       basis(2)%spin_mt = 1
       call create_vectors1
    end select
      
  end subroutine create_basis

  subroutine create_vectors2
    implicit none
    ! Local variables
    integer :: i,j

    allocate (basis(1)%vector(basis(1)%spin_mt))
    allocate (basis(2)%vector(basis(2)%spin_mt))

    do i=1,2
       do j=1,basis(i)%spin_mt
          basis(i)%vector(j) = - basis(i)%spin_mag + dble(j - 1)
       end do
    end do

  end subroutine create_vectors2

  subroutine create_vectors1
    implicit none
    ! Local variables
    integer :: j

    allocate (basis(1)%vector(basis(1)%spin_mt))
    allocate (basis(2)%vector(1))
 
    do j=1,basis(1)%spin_mt
       basis(1)%vector(j) = - basis(1)%spin_mag + dble(j - 1)
    end do

    basis(2)%vector(1) = 0.d0

  end subroutine create_vectors1
  
end module system_basis
