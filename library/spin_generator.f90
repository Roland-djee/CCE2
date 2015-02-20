!------------------------------------------------------------------------------
! CCE2 module
!------------------------------------------------------------------------------
!
! MODULE: spin_generator
!
!> @author
!> Dr. Roland Guichard University College London
!
! DESCRIPTION: 
!> Generates the spin operator matrices for any spin.
!
! REVISION HISTORY:
! 19-02-2015 - Initial Version
! TODO_19_02_2015 - Complete  - TODO_spin_generator
!------------------------------------------------------------------------------

module spin_generator
  use type
  implicit none

contains

  !---------------------------------------------------------------------------  
  !> @author 
  !> Dr. Roland Guichard University College London
  !
  ! DESCRIPTION: 
  !> Generates the spin operator matrices
  !> @brief
  !>
  !
  ! REVISION HISTORY:
  ! TODO_19_02_2015 - Finish spin matrices - TODO_spin_matrices
  !
  !> @param[in]  ms   
  !> @param[out] --      
  !> @return     Sz,S+,S-
  !---------------------------------------------------------------------------

  subroutine spin_matrices(ms)
    implicit none
    integer :: mt,i,j
    double precision, intent(in)  :: ms
    double precision :: m,mp,sqrt_fact

    mt = int(2.d0 * abs(ms) + 1.d0)
    allocate(Sz(mt,mt),Sp(mt,mt),Sm(mt,mt))
    
    Sz = 0.d0
    Sp = 0.d0
    Sm = 0.d0

    do i=1,mt
       m = - abs(ms) + dble(i - 1)
       do j=1,mt
          mp = - abs(ms) + dble(j - 1)
          sqrt_fact = sqrt(abs(ms)*(abs(ms)+1.d0) - mp*m)
          Sp(i,j) = delta1(i,j) * sqrt_fact
          Sm(i,j) = delta2(i,j) * sqrt_fact
          if (j==i) Sz(i,j) = m
       end do
    end do

  end subroutine spin_matrices

  function delta1(i,j)
    implicit none
    integer, intent(in) :: i,j
    double precision :: delta1

    if (j == i + 1) then
       delta1 = 1.d0
    else
       delta1 = 0.d0
    end if

  end function delta1

  function delta2(i,j)
    implicit none
    integer, intent(in) :: i,j
    double precision :: delta2

    if (i == j + 1) then
       delta2 = 1.d0
    else
       delta2 = 0.d0
    end if

  end function delta2

end module spin_generator
