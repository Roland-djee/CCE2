!------------------------------------------------------------------------------
! CCE2 code module
!------------------------------------------------------------------------------
!
! MODULE: build_hamiltonian
!
!> @author
!> Dr. Roland Guichard University College London
!
! DESCRIPTION: 
!> This module contains all subroutines necessary to build the Hamiltonian 
!> matrix
!
! REVISION HISTORY:
! 16-02-2015 - Initial Version
! TODO_16_02_2015 - Finish build_diag - TODO_build_hamiltonian
!------------------------------------------------------------------------------

module build_hamiltonian
  use type
  use constant
  implicit none

contains

  !---------------------------------------------------------------------------  
  !> @author 
  !> Dr. Roland Guichard University College London
  !
  ! DESCRIPTION: 
  !> Computes the diagonal part of the free Hamiltonian
  !> @brief
  !> Computes the Zeeman terms of the Hamiltonian gamma.Sz.B0
  !
  ! REVISION HISTORY:
  ! TODO_17_02_2015 - Finish select cases - TODO_build_diag
  !
  !> @param[in]  basis,B0
  !> @param[out] --      
  !> @return     H0_diag
  !---------------------------------------------------------------------------  

  subroutine build_diag(type,specie)
    implicit none
    ! Local variables
    character (len=20), intent(in) :: type,specie
    integer :: i,j,k
    double precision :: Zeeman1,Zeeman2

    select case (type)
    case ("CS")
       select case (specie)
       case ("Bi")
          Zeeman1 = gamma_n_209Bi
       case ("P")
          Zeeman1 = gamma_n_31P 
       case ("Si")
          Zeeman1 = gamma_n_29Si
       end select
       Zeeman2 = gamma_e
    case ("Bath")
       select case (specie)
       case ("Bi")
          Zeeman1 = gamma_n_209Bi
       case ("P")
          Zeeman1 = gamma_n_31P 
       case ("Si")
          Zeeman1 = gamma_n_29Si
       end select
       Zeeman2 = 0.d0
    end select

    allocate (H0_diag(basis(1)%spin_mt * basis(2)%spin_mt))

    k = 0
    do i=1,basis(1)%spin_mt
       do j=1,basis(2)%spin_mt
          k = k + 1
          H0_diag(k) = Zeeman1 * basis(1)%vector(i) + &
                       Zeeman2 * basis(2)%vector(j)
       end do
    end do

    deallocate (basis(1)%vector,basis(2)%vector)

    H0_diag = H0_diag * B0%ampli

  end subroutine build_diag

  !---------------------------------------------------------------------------  
  !> @author 
  !> Dr. Roland Guichard University College London
  !
  ! DESCRIPTION: 
  !> Computes the hyperfine Hamiltonian matrix
  !> @brief
  !> Computes the 1/2[I+S- + I_S+] + IzSz matrix elements of
  !> the hyperfine matrix.
  !
  ! REVISION HISTORY:
  ! TODO_19_02_2015 - Finish matrix elements - TODO_build_hf
  !
  !> @param[in]  basis
  !> @param[out] --      
  !> @return     H_hf
  !--------------------------------------------------------------------------- 

  subroutine build_hf
    implicit none

    
    
        
  end subroutine build_hf

end module build_hamiltonian
