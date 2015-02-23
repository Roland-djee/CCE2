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
  use spin_generator
  use write
  use math
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
          Zeeman2 = gamma_n_209Bi
       case ("P")
          Zeeman2 = gamma_n_31P 
       case ("Si")
          Zeeman2 = gamma_n_29Si
       end select
       Zeeman1 = 0.d0
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
  !> Computes the A * (1/2[I+S- + I-S+] + IzSz) matrix elements of
  !> the hyperfine matrix.
  !
  ! REVISION HISTORY:
  ! TODO_19_02_2015 - Finish set A values - TODO_build_hf
  !
  !> @param[in]  basis
  !> @param[out] --      
  !> @return     H_hf
  !--------------------------------------------------------------------------- 

  subroutine build_hf(bmag1,mt1,specie,bmag2,mt2)
    implicit none
    integer, intent(in) :: mt1,mt2
    double precision, intent(in) :: bmag1,bmag2
    character (len=20), intent(in) :: specie
    double precision :: A
    double precision :: Sz1(mt1,mt1),Sp1(mt1,mt1),Sm1(mt1,mt1)
    double precision :: Sz2(mt2,mt2),Sp2(mt2,mt2),Sm2(mt2,mt2)
    double precision :: K1(mt1*mt2,mt1*mt2)
    double precision :: K2(mt1*mt2,mt1*mt2)
    double precision :: K3(mt1*mt2,mt1*mt2)

    allocate(H_hf(mt1*mt2,mt1*mt2))

    !> Generate spin operator matrices
    call spin_matrices(bmag1,mt1,Sz1,Sp1,Sm1)
    call spin_matrices(bmag2,mt2,Sz2,Sp2,Sm2)

    !> Compute the Kronecker products
    call kronecker(Sp1,mt1,Sm2,mt2,K1) 
    call kronecker(Sm1,mt1,Sp2,mt2,K2)
    call kronecker(Sz1,mt1,Sz2,mt2,K3)

    select case (specie)
    case ("Bi")
       A = A_209Bi
    case ("P")
       write(*,*)'Set A value for P...'
       !A = A_31P
       A = 1.d0
       stop
    case ("Si")
       write(*,*)'Set A value for Si...'
       !A = A_29Si
       A = 1.d0
       stop
    end select

    H_hf = 0.5d0 * (K1 + K2) + K3
    H_hf = 2.d0 * pi * A * H_hf
         
  end subroutine build_hf

  !---------------------------------------------------------------------------  
  !> @author 
  !> Dr. Roland Guichard University College London
  !
  ! DESCRIPTION: 
  !> Computes the interaction matrix between a pair of bath spins.
  !> @brief
  !> Computes the -1/2[I+S- + I-S+] + IzSz matrix elements of
  !> the dipolar coupling interaction matrix.
  !
  ! REVISION HISTORY:
  ! TODO_20_02_2015 - Finish matrix elements - TODO_build_hf
  !
  !> @param[in]  basis
  !> @param[out] --      
  !> @return     H_int
  !--------------------------------------------------------------------------- 

  subroutine build_int(bmag1,mt1,specie,bmag2,mt2)
    implicit none
    integer, intent(in) :: mt1,mt2
    double precision, intent(in) :: bmag1,bmag2
    character (len=20), intent(in) :: specie
    double precision :: A
    double precision :: Sz1(mt1,mt1),Sp1(mt1,mt1),Sm1(mt1,mt1)
    double precision :: Sz2(mt2,mt2),Sp2(mt2,mt2),Sm2(mt2,mt2)
    double precision :: K1(mt1*mt2,mt1*mt2)
    double precision :: K2(mt1*mt2,mt1*mt2)
    double precision :: K3(mt1*mt2,mt1*mt2)

    allocate(H_int(mt1*mt2,mt1*mt2))

    !> Generate spin operator matrices
    call spin_matrices(bmag1,mt1,Sz1,Sp1,Sm1)
    call spin_matrices(bmag2,mt2,Sz2,Sp2,Sm2)

    !> Compute the Kronecker products
    call kronecker(Sp1,mt1,Sm2,mt2,K1) 
    call kronecker(Sm1,mt1,Sp2,mt2,K2)
    call kronecker(Sz1,mt1,Sz2,mt2,K3)

    select case (specie)
    case ("Bi")
       A = A_209Bi
    case ("P")
       write(*,*)'Set A value for P...'
       !A = A_31P
       A = 1.d0
       stop
    case ("Si")
       write(*,*)'Set A value for Si...'
       !A = A_29Si
       A = 1.d0
       !stop
    end select

    H_int = - 0.5d0 * (K1 + K2) + K3
    H_int = A * H_int  
    
  end subroutine build_int

end module build_hamiltonian
