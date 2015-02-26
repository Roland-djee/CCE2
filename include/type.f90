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

  integer :: cs_mt,bath_mt,tot_basis_mt
  double precision :: XJ1,XJ2,C12

  double precision, allocatable :: H0_diag(:)
  double precision, allocatable :: H_cs(:,:)
  double precision, allocatable :: H_int(:,:)
  double precision, allocatable :: H_shf(:,:)
  double precision, allocatable :: H_tot(:,:)

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
  end type basis_def
  type (basis_def) :: basis(2),electron

end module type
