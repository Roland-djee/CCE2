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
  use constant
  implicit none
  character (len=*), parameter :: fmt_str  = "(t50, a)"
  character (len=*), parameter :: fmt_str2 = "(t50, 2a)"
  character (len=*), parameter :: fmt_int  = "(t50, i)"
  character (len=*), parameter :: fmt_int3 = "(t50, 3i)"
  character (len=*), parameter :: fmt_dbl  = "(t50, es)"

contains

  !---------------------------------------------------------------------------  
  !> @author 
  !> Dr. Roland Guichard University College London
  !
  ! DESCRIPTION: 
  !> Reads input file basis.inp
  !> @brief
  !> Sets all variables for the basis of the spin system.
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

    read (basis_inp, fmt_int3) B0%x,B0%y,B0%z
    read (basis_inp, fmt_dbl) B0%ampli
    read (basis_inp, fmt_str) basis(1)%spin_type
    read (basis_inp, fmt_str) basis(1)%spin_sp
    if (basis(1)%spin_type .ne. 'CS' .and. basis(1)%spin_type .ne. 'Bath') then
       write(*,*)'Input Basis',trim(basis(1)%spin_type),' not recognized...'
       stop
    end if
    read (basis_inp, fmt_str) basis(2)%spin_type
    read (basis_inp, fmt_str) basis(2)%spin_sp
    
  end subroutine read_basis

end module read
