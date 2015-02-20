!------------------------------------------------------------------------------
! CCE2 module
!------------------------------------------------------------------------------
!
! MODULE: write
!
!> @author
!> Dr. Roland Guichard University College London
!
! DESCRIPTION: 
!> Write outputs as screen prints or files. 
!
! REVISION HISTORY:
! 19-02-2015 - Initial Version
! TODO_19_02_2015 - Complete  - TODO_write
!------------------------------------------------------------------------------

module write
  implicit none

contains

!---------------------------------------------------------------------------  
  !> @author 
  !> Dr. Roland Guichard University College London
  !
  ! DESCRIPTION: 
  !> Explicitly prints out matrices.
  !> @brief
  !> Format print out for matrices
  !
  ! REVISION HISTORY:
  ! TODO_19_02_2015 - Complete - TODO_PRINT_MATRIX
  !
  !> @param[in]  DESC, M, N, A, LDA   
  !> @param[out] --      
  !> @return     --
  !---------------------------------------------------------------------------

!  =============================================================================
!
!     Auxiliary routine: printing a matrix.
!
      SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA )
      CHARACTER*(*)    DESC
      INTEGER          M, N, LDA
      DOUBLE PRECISION A( LDA, * )
!
      INTEGER          I, J
!
      WRITE(*,*)
      WRITE(*,*) DESC
      DO I = 1, M
         WRITE(*,9998) ( A( I, J ), J = 1, N )
      END DO
!
! 9998 FORMAT( 11(:,1X,F6.2) )
 9998 FORMAT( 11(:,1X,ES13.5E3) )
      RETURN
      END SUBROUTINE PRINT_MATRIX
!
!  =============================================================================
!
!     Auxiliary routine: printing a complex matrix.
!
      SUBROUTINE PRINT_MATRIX_CMPLX( DESC, M, N, A, LDA )
      CHARACTER*(*)    DESC
      INTEGER          M, N, LDA
      DOUBLE COMPLEX   A( LDA, * )
!
      INTEGER          I, J
!
      WRITE(*,*)
      WRITE(*,*) DESC
      DO I = 1, M
         WRITE(*,9998) ( REAL(A( I, J )), AIMAG(A( I, J )), J = 1, N)
      END DO
      !DO I = 1, M
      !   WRITE(*,9998) ( REAL(A( I, J )), J = 1, N)
      !END DO
!
! 9998 FORMAT( 11(:,1X,F6.2) )
 9998 FORMAT( 11(:,1X,ES13.5E3) )
      RETURN
      END SUBROUTINE PRINT_MATRIX_CMPLX

end module write
