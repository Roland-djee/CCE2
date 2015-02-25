!------------------------------------------------------------------------------
! CCE2 module
!------------------------------------------------------------------------------
!
! MODULE: constant
!
!> @author
!> Dr. Roland Guichard University College London
!
! DESCRIPTION: 
!> Defines all physical constants needed for CCE2
!
! REVISION HISTORY:
! 16-02-2015 - Initial Version
!------------------------------------------------------------------------------

module constant
  implicit none

  !> Pi
  double precision, parameter :: pi = 3.141592653589793d0

  !> Vacuum permeability [H/m]
  double precision, parameter :: mu_0 = 4.d-7*pi

  !> Reduced Planck constant [Js]
  !double precision, parameter :: hbar = 1.05457172647D-34
  !> For lengths in Angstroms 
  double precision, parameter :: hbar = 1.05457172647d-4

  !> Free electron gyromagnetic ratio [rad/s/T]
  double precision, parameter :: gamma_e = 1.760859708d11
  !> Free electron spin
  double precision, parameter :: S_e = 0.5d0

  !> Moments for 29Si impurities

  !> 29Si nuclear gyromagnetic ratio [rad/s/T]
  double precision, parameter :: gamma_n_29Si = - 53.1903d6
  !> 29Si nuclear spin
  double precision, parameter :: I_29Si = 0.5d0

  !> Moments for donors
  
  !> Bismuth
  !> 209Bi nuclear gyromagnetic ratio [rad/s/T]
  double precision, parameter :: gamma_n_209Bi = - 43.775000006d6
  !> 209Bi hyperfine constant Hz
  double precision, parameter :: A_209Bi = 1.47539815345d9
  !> 209Bi nuclear spin
  double precision, parameter :: I_209Bi = 4.5d0
  !> 209Bi donor ionization energy [eV]
  double precision, parameter :: E_Bi = 0.069d0

  !> Phosphorus
  !> 31P nuclear gyromagnetic ratio [rad/s/T]
  double precision, parameter :: gamma_n_31P = 108.25928d6
  !> 31P nuclear spin
  double precision, parameter :: I_31P = 0.5d0

  !> Charge density
  double precision, parameter :: eta = 186.d0
  !> lengths for the hyperfine KL wavefunction [angstrom]
  double precision, parameter :: a   = 25.09d0
  double precision, parameter :: b   = 14.43d0
  !> Ionization energy [eV]
  double precision, parameter :: E_p = 0.044d0
  !> n parameter entering the KLW
  double precision, parameter :: n   = dsqrt(0.029d0/E_p)
  !> lattice parameter [angstrom]
  double precision, parameter :: a0  = 5.43d0
  !> k0 parameter entering the KLW
  double precision, parameter :: k0  = 2.d0*pi*0.85d0/a0
  !> Prefactors for the dipolar/hyperfine couplings
  double precision, parameter :: pref = mu_0*gamma_n_29Si**2*hbar/(4.d0*pi)
  double precision, parameter :: p = 4.d0*gamma_e*gamma_n_29Si*hbar*Eta*&
       &mu_0/9.D0
  
  !> Scaling factor sites coordinates => angstroms
  double precision, parameter :: scaling = a0/4.d0

end module constant
