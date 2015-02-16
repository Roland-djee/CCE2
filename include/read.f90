!------------------------------------------------------------------------------
! CCE2 module
!------------------------------------------------------------------------------
!
! MODULE: read
!
!> @author
!> Dr. Roland Guichard University College London
!
! DESCRIPTION: 
!> Reads all input files. 
!
! REVISION HISTORY:
! 16-02-2015 - Initial Version
! TODO_16_02_2015 - Finish module - TODO_read
!------------------------------------------------------------------------------

module read
  use type
  implicit none
  character (len=*), parameter :: fmt_str = "(t50, a)"
  character (len=*), parameter :: fmt_dbl = "(t50, es)"

contains

  !---------------------------------------------------------------------------  
  !> @author 
  !> Dr. Roland Guichard University College London
  !
  ! DESCRIPTION: 
  !> Reads input file basis.inp
  !> @brief
  !> Sets all variables for the vector basis of the spin system.
  !
  ! REVISION HISTORY:
  ! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
  !
  !> @param[in]  basis.inp   
  !> @param[out] --      
  !> @return     basis
  !---------------------------------------------------------------------------

  subroutine read_basis
    implicit none
    integer :: basis_inp = 10
    
    open(unit = basis_inp, file = '../input/basis.inp')

    read (basis_inp, fmt_str) basis(1)%spin_type
    if (basis(1)%spin_type == 'CS') then
    else if (trim(basis(1)%spin_type) == 'Bath') then
    else
       write(*,*)'Input Basis',trim(basis(1)%spin_type),' not recognized...'
       stop
    end if
    read (basis_inp, fmt_dbl) basis(1)%spin_mag
    read (basis_inp, fmt_str) basis(2)%spin_type
    read (basis_inp, fmt_dbl) basis(2)%spin_mag
    
    basis(1)%spin_mt = int(2.d0 * basis(1)%spin_mag + 1.d0)
    basis(2)%spin_mt = int(2.d0 * basis(2)%spin_mag + 1.d0)
    
    allocate (basis(1)%vector(basis(1)%spin_mt))
    allocate (basis(2)%vector(basis(2)%spin_mt))
    allocate (H0_diag(basis(1)%spin_mt * basis(2)%spin_mt))
    
  end subroutine read_basis

end module read
