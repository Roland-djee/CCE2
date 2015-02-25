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
  !> Sets all the parameters for the construction of the basis.
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

  subroutine create_basis
    implicit none
    ! local variables
    integer :: i
    character (len=20) :: type,specie

    do i=1,2
       type   = basis(i)%spin_type
       specie = basis(i)%spin_sp
       select case (type)
       case ("CS")     
          !> Central spin system {donor + e-}
          !> Sets spin magnitudes and total m.
          select case (specie)
          case ("Bi")
             basis(i)%spin_mag = I_209Bi
             basis(i)%spin_mt  = int(2.d0 * basis(i)%spin_mag + 1.d0)
          case ("P")
             basis(i)%spin_mag = I_31P
             basis(i)%spin_mt  = int(2.d0 * basis(i)%spin_mag + 1.d0)
          case ("Si")
             basis(i)%spin_mag = I_29Si
             basis(i)%spin_mt  = int(2.d0 * basis(i)%spin_mag + 1.d0)
          case ("e")
             write(*,*)'To be finished...'
          case ("n")
             write(*,*)'To be finished...'
          end select
          !> Central spin total spin including the free electron
          electron%spin_mag = S_e
          electron%spin_mt  = int(2.d0 * electron%spin_mag + 1.d0)
          cs_mt = basis(i)%spin_mt * electron%spin_mt
       case ("Bath")
          !> Bath spin system
          !> Sets spin magnitudes and total m.
          select case (specie)
          case ("Bi")
             basis(i)%spin_mag = I_209Bi
             basis(i)%spin_mt  = int(2.d0 * basis(i)%spin_mag + 1.d0)
          case ("P")
             basis(i)%spin_mag = I_31P
             basis(i)%spin_mt  = int(2.d0 * basis(i)%spin_mag + 1.d0)
          case ("Si")
             basis(i)%spin_mag = I_29Si
             basis(i)%spin_mt  = int(2.d0 * basis(i)%spin_mag + 1.d0)
          case ("e")
             write(*,*)'To be finished...'
          case ("n")
             write(*,*)'To be finished...'
          end select
          !> Bath pair total spin
          bath_mt = basis(i)%spin_mt**2
       end select
    end do
    
    !> Toaol basis dimension
    tot_basis_mt = cs_mt * bath_mt
      
  end subroutine create_basis
 
end module system_basis
