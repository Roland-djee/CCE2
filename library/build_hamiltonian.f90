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
  !> Computes the off-diagonal part of the Hamiltonian
  !> @brief
  !> Computes successively the central spin (cs), impurity pair interaction 
  !> (int) and superhyperfine (shf) Hamiltonian matrices
  !
  ! REVISION HISTORY:
  ! TODO_24_02_2015 - Finish select cases - TODO_build_interaction_matrices
  !
  !> @param[in]  --
  !> @param[out] --      
  !> @return     H_cs,H_int,H_shf
  !---------------------------------------------------------------------------  

  subroutine build_interaction_matrices
    implicit none

    integer :: mt1 = basis(1)%spin_mt
    integer :: mt2 = basis(2)%spin_mt

    double precision :: Sz1(mt1,mt1),Sp1(mt1,mt1),Sm1(mt1,mt1)
    double precision :: Sze(mt1,mt1),Spe(mt1,mt1),Sme(mt1,mt1)
    double precision :: Sz2(mt2,mt2),Sp2(mt2,mt2),Sm2(mt2,mt2)
    double precision :: K1(mt1*mt2,mt1*mt2)
    double precision :: K2(mt1*mt2,mt1*mt2)
    double precision :: K3(mt1*mt2,mt1*mt2)

    !> Generate spin operator matrices
    call spin_matrices(bmag1,mt1,Sz1,Sp1,Sm1)
    !CALL PRINT_MATRIX( 'Sp1', 2, 2, Sp1, 2 )
    !CALL PRINT_MATRIX( 'Sm1', 2, 2, Sm1, 2 )
    !CALL PRINT_MATRIX( 'Sz1', 2, 2, Sz1, 2 )  

    call spin_matrices(bmag2,mt2,Sz2,Sp2,Sm2)
    !CALL PRINT_MATRIX( 'Sp2', 2, 2, Sp2, 2 )
    !CALL PRINT_MATRIX( 'Sm2', 2, 2, Sm2, 2 )
    !CALL PRINT_MATRIX( 'Sz2', 2, 2, Sz2, 2 )  
    

  end subroutine build_interaction_matrices

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
  !> @param[in]  type,specie
  !> @param[out] --      
  !> @return     H0_diag
  !---------------------------------------------------------------------------  

  subroutine build_diag
    implicit none
    ! Local variables
    character (len=20) :: type,specie
    integer :: i,j,k,l,m
    double precision :: vector1,vector2,vector3,vector4
    double precision :: Zeeman1,Zeeman2

    do i=1,2
       type   = basis(i)%spin_type
       specie = basis(i)%spin_sp
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
       case ("Bath")
          select case (specie)
          case ("Bi")
             Zeeman2 = gamma_n_209Bi
          case ("P")
             Zeeman2 = gamma_n_31P 
          case ("Si")
             Zeeman2 = gamma_n_29Si
          end select
       end select
    end do

    allocate (H0_diag(tot_basis_mt))

    m = 0
    !> Central spin loop
    do i=1,basis(1)%spin_mt
       vector1 = - basis(1)%spin_mag + dble(i - 1)
       !> Free electron loop
       do j=1,electron%spin_mt
          vector2 = - electron%spin_mag + dble(j - 1)
          !> First impurity loop
          do k=1,basis(2)%spin_mt
             vector3 = - basis(2)%spin_mag + dble(k - 1)
             !> Second impurity loop
             do l=1,basis(2)%spin_mt
                vector4 = - basis(2)%spin_mag + dble(l - 1)
                m = m + 1
                H0_diag(m) = Zeeman1 * vector1 + &
                             gamma_e * vector2 + &
                             Zeeman2 * vector3 + &
                             Zeeman2 * vector4               
             end do
          end do
       end do
    end do

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
  !> Computes the interaction matrix between members of a bath spin pair.
  !> @brief
  !> Computes the -1/2(1/2[I+S- + I-S+]) + IzSz matrix elements of
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
    !CALL PRINT_MATRIX( 'Sp1', 2, 2, Sp1, 2 )
    !CALL PRINT_MATRIX( 'Sm1', 2, 2, Sm1, 2 )
    !CALL PRINT_MATRIX( 'Sz1', 2, 2, Sz1, 2 )  

    call spin_matrices(bmag2,mt2,Sz2,Sp2,Sm2)
    !CALL PRINT_MATRIX( 'Sp2', 2, 2, Sp2, 2 )
    !CALL PRINT_MATRIX( 'Sm2', 2, 2, Sm2, 2 )
    !CALL PRINT_MATRIX( 'Sz2', 2, 2, Sz2, 2 )  

    !> Compute the Kronecker products
    call kronecker(Sp1,mt1,Sm2,mt2,K1) 
    call kronecker(Sm1,mt1,Sp2,mt2,K2)
    call kronecker(Sz1,mt1,Sz2,mt2,K3)

    !CALL PRINT_MATRIX( 'K1', 4, 4, K1, 4 )
    !CALL PRINT_MATRIX( 'K2', 4, 4, K2, 4 )
    !CALL PRINT_MATRIX( 'K3', 4, 4, K3, 4 )  

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

    H_int = - 0.25d0 * (K1 + K2) + K3
    H_int = A * H_int  
    
  end subroutine build_int

  subroutine build_int_12(bmag1,mt1,bmag2,mt2,bmag3,mt3)
    implicit none
    integer, intent(in) :: mt1,mt2,mt3
    double precision, intent(in) :: bmag1,bmag2,bmag3
    
    double precision :: Sz1(mt1,mt1),Sp1(mt1,mt1),Sm1(mt1,mt1)
    double precision :: Sz2(mt2,mt2),Sp2(mt2,mt2),Sm2(mt2,mt2)
    double precision :: Sz3(mt3,mt3),Sp3(mt3,mt3),Sm3(mt3,mt3)

    double precision :: Sp12(mt1*mt2,mt1*mt2)
    double precision :: Sm12(mt1*mt2,mt1*mt2)
    double precision :: Sz12(mt1*mt2,mt1*mt2)

    double precision :: K1(mt1*mt2*mt3,mt1*mt2*mt3)
    double precision :: K2(mt1*mt2*mt3,mt1*mt2*mt3)
    double precision :: K3(mt1*mt2*mt3,mt1*mt2*mt3)

    allocate(H_int_12(mt1*mt2*mt3,mt1*mt2*mt3))
   
    !> Generate spin operator matrices
    call spin_matrices(bmag1,mt1,Sz1,Sp1,Sm1)
    CALL PRINT_MATRIX( 'Sp1', 2, 2, Sp1, 2 )
    CALL PRINT_MATRIX( 'Sm1', 2, 2, Sm1, 2 )
    CALL PRINT_MATRIX( 'Sz1', 2, 2, Sz1, 2 )  

    call spin_matrices(bmag2,mt2,Sz2,Sp2,Sm2)
    CALL PRINT_MATRIX( 'Sp2', 2, 2, Sp2, 2 )
    CALL PRINT_MATRIX( 'Sm2', 2, 2, Sm2, 2 )
    CALL PRINT_MATRIX( 'Sz2', 2, 2, Sz2, 2 )  

    call spin_matrices(bmag3,mt3,Sz3,Sp3,Sm3)
    CALL PRINT_MATRIX( 'Sp3', 2, 2, Sp3, 2 )
    CALL PRINT_MATRIX( 'Sm3', 2, 2, Sm3, 2 )
    CALL PRINT_MATRIX( 'Sz3', 2, 2, Sz3, 2 )  

    !> Compute the Kronecker products for spins 1 and 2
    !> S+12 = S+1 x S+2 etc.
    call kronecker(Sp1,mt1,Sp2,mt2,Sp12) 
    call kronecker(Sm1,mt1,Sm2,mt2,Sm12)
    call kronecker(Sz1,mt1,Sz2,mt2,Sz12)

    CALL PRINT_MATRIX( 'Sp12', 4, 4, Sp12, 4 )
    CALL PRINT_MATRIX( 'Sm12', 4, 4, Sm12, 4 )
    CALL PRINT_MATRIX( 'Sz12', 4, 4, Sz12, 4 )

    !CALL PRINT_MATRIX_BLOCK( 'Sp12 1x1', 1, 1, 10, 10, Sp12, 20 )
    !CALL PRINT_MATRIX_BLOCK( 'Sp12 1x2', 1, 11, 10, 20, Sp12, 20 )
    !CALL PRINT_MATRIX_BLOCK( 'Sp12 2x1', 11, 1, 20, 10, Sp12, 20 )
    !CALL PRINT_MATRIX_BLOCK( 'Sp12 2x2', 11, 11, 20, 20, Sp12, 20 )
 
    !> Compute the Kronecker products for spins 12 and 3
    call kronecker(Sp12,mt1*mt2,Sm3,mt3,K1) 
    call kronecker(Sm12,mt1*mt2,Sp3,mt3,K2) 
    call kronecker(Sz12,mt1*mt2,Sz3,mt3,K3)

    CALL PRINT_MATRIX( 'K1', 8, 8, K1, 8 )
    CALL PRINT_MATRIX( 'K2', 8, 8, K2, 8 )
    CALL PRINT_MATRIX( 'K3', 8, 8, K3, 8 )
    !CALL PRINT_MATRIX_BLOCK( 'Sp12 1x1', 1, 1, 10, 10, Sp12, 20 )
    !CALL PRINT_MATRIX_BLOCK( 'Sp12 1x2', 1, 11, 10, 20, Sp12, 20 )
    !CALL PRINT_MATRIX_BLOCK( 'Sp12 2x1', 11, 1, 20, 10, Sp12, 20 )
    !CALL PRINT_MATRIX_BLOCK( 'Sp12 2x2', 11, 11, 20, 20, Sp12, 20 )

    H_int_12 = 0.5d0 * (K1 + K2) + K3

  end subroutine build_int_12

end module build_hamiltonian
