module read
  use type
  implicit none
  character (len=*), parameter :: fmt_str = "(t50, a)"
  character (len=*), parameter :: fmt_dbl = "(t50, es)"

contains

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

end subroutine read_basis

end module read
