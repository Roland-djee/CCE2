!------------------------------------------------------------------------------
! CCE2 module
!------------------------------------------------------------------------------
!
! MODULE: type
!
!> @author
!> Dr. Roland Guichard University College London
!
! DESCRIPTION: 
!> Defines all variable, array etc. types needed for CCE2
!
! REVISION HISTORY:
! 16-02-2015 - Initial Version
! TODO_16_02_2015 - Deallocate basis%vector - TODO_type
! TODO_16_02_2015 - Deallocate H0_diag - TODO_type
! TODO_19_02_2015 - Deallocate H_hf,Sz,Sp,Sm - TODO_type
!------------------------------------------------------------------------------

module type
  implicit none

  integer :: basis_nb

  double precision, allocatable :: H0_diag(:)
  double precision, allocatable :: H_hf(:,:),H_int(:,:)
  double precision, allocatable :: HCS(:,:)

  type vector
     sequence
     integer :: x,y,z
     double precision :: ampli
  end type vector

  type (vector) :: B0

  type basis_def
     sequence
     character (len=20) :: spin_type, spin_sp
     integer :: spin_mt
     double precision   :: spin_mag
     double precision, allocatable :: vector(:)
  end type basis_def
  type (basis_def) :: basis(2)

end module type
