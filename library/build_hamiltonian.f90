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
!> matrix.
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

    integer :: mt1,mt2,mte
    double precision :: bmag1,bmag2,bmage
    character (len=20) :: specie

    double precision, allocatable :: Sz1(:,:),Sp1(:,:),Sm1(:,:)
    double precision, allocatable :: Sze(:,:),Spe(:,:),Sme(:,:)
    double precision, allocatable :: Sz2(:,:),Sp2(:,:),Sm2(:,:)
    double precision, allocatable :: K1(:,:)
    double precision, allocatable :: K2(:,:)
    double precision, allocatable :: K3(:,:)

    bmag1 = basis(1)%spin_mag
    mt1   = basis(1)%spin_mt   
    bmage = electron%spin_mag
    mte   = electron%spin_mt
    bmag2 = basis(2)%spin_mag
    mt2   = basis(2)%spin_mt

    allocate (Sz1(mt1,mt1),Sp1(mt1,mt1),Sm1(mt1,mt1))
    allocate (Sze(mte,mte),Spe(mte,mte),Sme(mte,mte))
    allocate (Sz2(mt2,mt2),Sp2(mt2,mt2),Sm2(mt2,mt2))

    !> Generate spin operator matrices for the central spin system

    call spin_matrices(bmag1,mt1,Sz1,Sp1,Sm1)
    CALL PRINT_MATRIX( 'Sp1', 2, 2, Sp1, 2 )
    CALL PRINT_MATRIX( 'Sm1', 2, 2, Sm1, 2 )
    CALL PRINT_MATRIX( 'Sz1', 2, 2, Sz1, 2 )  

    call spin_matrices(bmage,mte,Sze,Spe,Sme)
    CALL PRINT_MATRIX( 'Spe', 2, 2, Spe, 2 )
    CALL PRINT_MATRIX( 'Sme', 2, 2, Sme, 2 )
    CALL PRINT_MATRIX( 'Sze', 2, 2, Sze, 2 )  

    !> Generate the hyperfine matrix for the central spin system
    specie = basis(1)%spin_sp
    call build_cs(bmag1,mt1,specie,Sz1,Sp1,Sm1,bmage,mte,Sze,Spe,Sme,mt2)

    CALL PRINT_MATRIX_BLOCK( 'H_cs 1x1', 1, 1, 8, 8, H_cs, 16 )
    CALL PRINT_MATRIX_BLOCK( 'H_cs 1x2', 1, 9, 8, 16, H_cs, 16 )
    CALL PRINT_MATRIX_BLOCK( 'H_cs 2x1', 9, 1, 16, 8, H_cs, 16 )
    CALL PRINT_MATRIX_BLOCK( 'H_cs 2x2', 9, 9, 16, 16, H_cs, 16 )

    !> Generate spin operator matrices for the impurity
 
    call spin_matrices(bmag2,mt2,Sz2,Sp2,Sm2)

    !> Generate the interaction matrix for a bath pair
    specie = basis(2)%spin_sp
    call build_int(bmag2,mt2,Sz2,Sp2,Sm2,mt1,mte)

    CALL PRINT_MATRIX_BLOCK( 'H_int 1x1', 1, 1, 8, 8, H_int, 16 )
    CALL PRINT_MATRIX_BLOCK( 'H_int 1x2', 1, 9, 8, 16, H_int, 16 )
    CALL PRINT_MATRIX_BLOCK( 'H_int 2x1', 9, 1, 16, 8, H_int, 16 )
    CALL PRINT_MATRIX_BLOCK( 'H_int 2x2', 9, 9, 16, 16, H_int, 16 )

    call build_shf(mt1,bmage,mte,Sze,Spe,Sme,&
                       bmag2,mt2,Sz2,Sp2,Sm2)

    CALL PRINT_MATRIX_BLOCK( 'H_shf 1x1', 1, 1, 8, 8, H_shf, 16 )
    CALL PRINT_MATRIX_BLOCK( 'H_shf 1x2', 1, 9, 8, 16, H_shf, 16 )
    CALL PRINT_MATRIX_BLOCK( 'H_shf 2x1', 9, 1, 16, 8, H_shf, 16 )
    CALL PRINT_MATRIX_BLOCK( 'H_shf 2x2', 9, 9, 16, 16, H_shf, 16 )
    
    deallocate (Sz1,Sp1,Sm1)
    deallocate (Sz2,Sp2,Sm2)
    deallocate (Sze,Spe,Sme)

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
  !> Computes the hyperfine Hamiltonian matrix of the central spin system.
  !> @brief
  !> Computes the A * (1/2[I+S- + I-S+] + IzSz) matrix elements of
  !> the hyperfine matrix for the central spin system.
  !
  ! REVISION HISTORY:
  ! TODO_19_02_2015 - Finish set A values - TODO_build_hf
  !
  !> @param[in]  bmag1,mt1,specie,Sz1,Sp1,Sm1,bmag2,mt2,Sz2,Sp2,Sm2
  !> @param[out] --      
  !> @return     H_cs
  !--------------------------------------------------------------------------- 

  subroutine build_cs(bmag1,mt1,specie,Sz1,Sp1,Sm1,bmage,mte,Sze,Spe,Sme,mt2)
    implicit none
    integer, intent(in) :: mt1,mt2,mte
    double precision,  intent(in) :: bmag1,bmage
    character(len=20), intent(in) :: specie
    double precision,  intent(in) :: Sz1(mt1,mt1),Sp1(mt1,mt1),Sm1(mt1,mt1)
    double precision,  intent(in) :: Sze(mte,mte),Spe(mte,mte),Sme(mte,mte)

    integer :: i
    double precision :: A

    double precision :: H0(mt1*mte,mt1*mte)
    double precision :: K1(mt1*mte,mt1*mte)
    double precision :: K2(mt1*mte,mt1*mte)
    double precision :: K3(mt1*mte,mt1*mte)

    double precision :: Id(mt2*mt2,mt2*mt2)

    !> Compute the Kronecker products
    call kronecker(Sp1,mt1,Sme,mte,K1) 
    call kronecker(Sm1,mt1,Spe,mte,K2)
    call kronecker(Sz1,mt1,Sze,mte,K3)

    select case (specie)
    case ("Bi")
       A = A_209Bi
    case ("P")
       write(*,*)'Set A value for P...'
       !A = A_31P
       A = 1.d0
       !A = A_209Bi
    case ("Si")
       write(*,*)'Set A value for Si...'
       !A = A_29Si
       A = 1.d0
    end select

    H0 = 0.5d0 * (K1 + K2) + K3
    H0 = 2.d0 * pi * A * H0

    !> Raise H_cs to the full dimensions of the system {cs + pair}

    allocate(H_cs(mt1*mte*mt2*mt2,mt1*mte*mt2*mt2))

    !> Identity matrix
    Id = 0.d0
    forall (i=1:mt2*mt2) Id(i,i) = 1.d0
    
    call kronecker(H0,mt1*mte,Id,mt2*mt2,H_cs) 
         
  end subroutine build_cs

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
  !> @param[in]  bmag1,mt1,specie,Sz1,Sp1,Sm1
  !> @param[out] --      
  !> @return     H_int
  !--------------------------------------------------------------------------- 

  subroutine build_int(bmag2,mt2,Sz2,Sp2,Sm2,mt1,mte)
    implicit none
    integer, intent(in) :: mt1,mt2,mte
    double precision, intent(in) :: bmag2
    double precision, intent(in) :: Sz2(mt2,mt2),Sp2(mt2,mt2),Sm2(mt2,mt2)

    integer :: i
    double precision :: H0(mt2*mt2,mt2*mt2)
    double precision :: K1(mt2*mt2,mt2*mt2)
    double precision :: K2(mt2*mt2,mt2*mt2)
    double precision :: K3(mt2*mt2,mt2*mt2)

    double precision :: Id(mt1*mte,mt1*mte)

    !> Compute the Kronecker products
    call kronecker(Sp2,mt2,Sm2,mt2,K1) 
    call kronecker(Sm2,mt2,Sp2,mt2,K2)
    call kronecker(Sz2,mt2,Sz2,mt2,K3)

    !CALL PRINT_MATRIX( 'K1', 4, 4, K1, 4 )
    !CALL PRINT_MATRIX( 'K2', 4, 4, K2, 4 )
    !CALL PRINT_MATRIX( 'K3', 4, 4, K3, 4 )  

    H0 = - 0.25d0 * (K1 + K2) + K3
 
    !> Raise H_int to the full dimensions of the system {cs + pair}

    allocate(H_int(mt1*mte*mt2*mt2,mt1*mte*mt2*mt2))

    !> Identity matrix
    Id = 0.d0
    forall (i=1:mt1*mte) Id(i,i) = 1.d0
    
    call kronecker(Id,mt1*mte,H0,mt2*mt2,H_int)
    
  end subroutine build_int

  !---------------------------------------------------------------------------  
  !> @author 
  !> Dr. Roland Guichard University College London
  !
  ! DESCRIPTION: 
  !> Computes the interaction matrix between the central electron spin and 
  !> one bath spin impurity.
  !> @brief
  !> Computes the 1/2[S_cs+I- + S_cs-I+] + S_cszIz matrix elements of
  !> the superhyperfine coupling interaction matrix.
  !
  ! REVISION HISTORY:
  ! TODO_20_02_2015 - Finish matrix elements - TODO_build_hf
  !
  !> @param[in]  bmag1,mt1,Sz1,Sp1,Sm1,bmage,mte,Sze,Spe,Sme
  !> @param[in]  bmag2,mt2,Sz2,Sp2,Sm2
  !> @param[out] --      
  !> @return     H_shf
  !---------------------------------------------------------------------------

  subroutine build_shf(mt1,bmage,mte,Sze,Spe,Sme,&
                           bmag2,mt2,Sz2,Sp2,Sm2)
    implicit none
    integer, intent(in) :: mt1,mt2,mte
    double precision, intent(in) :: bmag2,bmage
       
    double precision, intent(in) :: Sz2(mt2,mt2),Sp2(mt2,mt2),Sm2(mt2,mt2)
    double precision, intent(in) :: Sze(mte,mte),Spe(mte,mte),Sme(mte,mte)

    integer :: i

    double precision :: H0(mt2*mte,mt2*mte)
    double precision :: K1(mt2*mte,mt2*mte)
    double precision :: K2(mt2*mte,mt2*mte)
    double precision :: K3(mt2*mte,mt2*mte)

    double precision :: Id(mt2,mt2)

    !> Compute the Kronecker products between the central spin and
    !> one impurity
    call kronecker(Spe,mte,Sm2,mt2,K1) 
    call kronecker(Sme,mte,Sp2,mt2,K2) 
    call kronecker(Sze,mte,Sz2,mt2,K3)

    CALL PRINT_MATRIX( 'K1', 4, 4, K1, 4 )
    CALL PRINT_MATRIX( 'K2', 4, 4, K2, 4 )
    CALL PRINT_MATRIX( 'K3', 4, 4, K3, 4 )
    !CALL PRINT_MATRIX_BLOCK( 'Sp12 1x1', 1, 1, 10, 10, Sp12, 20 )
    !CALL PRINT_MATRIX_BLOCK( 'Sp12 1x2', 1, 11, 10, 20, Sp12, 20 )
    !CALL PRINT_MATRIX_BLOCK( 'Sp12 2x1', 11, 1, 20, 10, Sp12, 20 )
    !CALL PRINT_MATRIX_BLOCK( 'Sp12 2x2', 11, 11, 20, 20, Sp12, 20 )

    H0 = 0.5d0 * (K1 + K2) + K3

    !> Raise H_shf to the full dimensions of the system {cs + pair}

    allocate(H_shf(mt1*mte*mt2*mt2,mt1*mte*mt2*mt2))

    !> Identity matrix
    Id = 0.d0
    forall (i=1:mt2) Id(i,i) = 1.d0
    
    call kronecker(H0,mt1*mt2*mte,Id,mt2,H_shf)

  end subroutine build_shf

end module build_hamiltonian
