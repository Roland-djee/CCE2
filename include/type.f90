module type
  implicit none

  integer :: basis_nb

  type basis_def
     sequence
     character (len=20) :: spin_type
     integer :: spin_mt
     double precision   :: spin_mag
     double precision, allocatable :: vector(:)
  end type basis_def
  type (basis_def) :: basis(2)

end module type
