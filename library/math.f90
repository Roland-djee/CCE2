!------------------------------------------------------------------------------
! CCE2 module
!------------------------------------------------------------------------------
!
! MODULE: math
!
!> @author
!> Dr. Roland Guichard University College London
!
! DESCRIPTION: 
!> Computes some mathematical operations.
!
! REVISION HISTORY:
! 19-02-2015 - Initial Version
! TODO_19_02_2015 - Complete  - TODO_math
!------------------------------------------------------------------------------

module math
  implicit none

contains

  !---------------------------------------------------------------------------  
  !> @author 
  !> Dr. Roland Guichard University College London
  !
  ! DESCRIPTION: 
  !> Computes the Kronecker product between two matrices.
  !> @brief
  !> Computes the Kronecker product between two square matrices A and B.
  !
  ! REVISION HISTORY:
  ! TODO_19_02_2015 - Complete - TODO_kronecker
  !
  !> @param[in]  A,B
  !> @param[out] --      
  !> @return     K
  !---------------------------------------------------------------------------

  subroutine kronecker(A,nA,B,nB,K)
    implicit none
    integer :: i,j
    integer, intent(in) :: nA,nB
    double precision, intent(in) :: A(nA,nA),B(nB,nB)
    double precision, intent(out):: K(nA*nB,nA*nB) 

    K = 0.d0
    
    do i=1,nA
       do j=1,nA
          K(1+(i-1)*nB:i*nb,1+(j-1)*nB:j*nb) = A(i,j) * B    
       end do
    end do

  end subroutine kronecker

end module math


