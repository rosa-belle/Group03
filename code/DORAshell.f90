!
! 
! D.O.R.A. SHELL-MODEL CODE
!
! To compile the main source: gfortran -c ModuleSource
! To compile the main source: gfortran -O2 -o Executable MainSource
!
! A. Forney, R.-B. Gerst, D. Muir, O. Vasseur
! ==============================================================================


! ==============================================================================
! GLOBAL VARIABLES THAT ARE NEITHER INPUT NOR PARAMETERS
! ==============================================================================
module global
  
  use input_DORA

  implicit none

  ! "Include" parameters (maximum values...)
  include 'param_DORA.f90'
  
  integer, parameter :: NN = A-Z
  
  type sp     ! N.B.: 'structure' and 'record' statements are pre-Fortran 90 statements that are equivalent to 'type' statements, but they are obsolete (not portable). Make sure to always use 'type'.
     integer :: ind   ! Index
     integer :: n     ! Principal quantum number.
     integer :: l     ! Orbital momentum.
     real :: j     ! Total angular momentum.
     real :: mj     ! Total angular momentum projection.
     real :: spe   !Single Particle Energy
  end type sp
  type(sp), dimension(mnsp) :: shell     ! One-dimensional array of sets of quantum numbers (single-particle states).

  integer :: dim     ! Dimension of the basis (it is determined in subroutine 'build_basis').
  integer, dimension (:, :), allocatable :: basis     ! Array in which each row is a Slater determinant (therefore its first dimension is the dimension of the basis).
  integer, dimension (:, :), allocatable :: basis_pair     ! Array in which each row is a Slater determinant for the Pairing Problem (therefore its first dimension is the dimension of the basis).

end module global


! ==============================================================================
! MAIN PROGRAM
! ==============================================================================
program DORA

  implicit none
  
  ! Test for input file
  call read_input

  ! Open the interaction file
  !!! TO DO !!!

  ! Build single-particle states
  call build_sp_states

  ! Build the basis of Slater determinants
  call build_basis

  ! Build Hamiltonian
  call build_hamiltonian

  ! Solve eigenproblem
  !call eigensolver

end program DORA


! ==============================================================================
! READ INPUT FILE
! ==============================================================================
subroutine read_input
  
  use input_DORA
  use global

  implicit none

!  inquire(file='', exist=file_exists)

  ! Test
  if ( mod(NN+int(2*M),2) .eq. 1 ) then
     write(6, *) " ERROR in subroutine 'input_DORA': Inconsistent M input. Exiting."
     stop
  end if
  write(6, *) 'Z= ', Z
  write(6, *) 'N= ', NN
  write(6, *) 'M= ', M

end subroutine read_input



! ==============================================================================
! BUILD SINGLE-PARTICLE STATES AND STORE THEM IN A FILE
! ==============================================================================
subroutine build_sp_states

  use global

  implicit none

  logical :: flag
  integer :: kn, kl, spi, im      ! spi: single-particle state index in the 'shell' array.
  real :: kj, km

  open(unit=4, status='unknown', file='sp_states.dat')     ! File where all the single-particle states are stored.

  ! Write header in 'sp_states.dat'
  write(4, *) 'SINGLE-PARTICLE STATES BUILT BY THE D.O.R.A. SHELL-MODEL CODE'
  write(4, *)
  write(4, *) '          #           n          l           j          m_j            SPE'
  
  ! Loop over principal quantum numbers
  kn = -1
  kl = 0 !For Pairing Problem
  spi = 0
  flag = .false.
 
  do
     if ( flag .eqv. .true. ) then     ! If verified, then the sp array is already full.
        exit
     end if
     kn = kn + 1
     !do kl = 0, kn   !Commented for Pairing Problem
        if ( flag .eqv. .true. ) then     ! If verified, then the sp array is already full.
           exit
        end if
        ! First "m" loop for j=l+1/2
        kj = real(kl) + 0.5
        do im = 0, int(2*kj)     ! Introducing a new summation index because it needs to be integer.
           if ( spi .ge. mnsp ) then     ! If verified, then the sp array is already full.
              flag = .true.
              exit
           end if
           km = -kj + real(im)
           ! Fill in the single-particle state array
           spi = spi + 1
           shell(spi)%n = kn
           shell(spi)%l = kl
           shell(spi)%j = kj
           shell(spi)%mj = km
           shell(spi)%spe = kn
           shell(spi)%ind = spi
           write(4, *) shell(spi)%ind, shell(spi)%n, shell(spi)%l, shell(spi)%j, shell(spi)%mj, shell(spi)%spe
        end do
        ! First "m" loop for j=l-1/2
        kj = real(kl) - 0.5
        do im = 0, int(2*kj)     ! Introducing a new summation index because it needs to be integer. N.B.: Because of the way the number of iterations in a do loop is natively calculated in Fortran, the case where 2*kj = -1 (it only happens for kl=0 and kj=kl-1/2) results in a
           if ( spi .ge. mnsp ) then     ! If verified, then the sp array is already full.
              exit
           end if
           km = -kj + real(im)
           ! Fill in the single-particle state array
           spi = spi + 1
           shell(spi)%n = kn
           shell(spi)%l = kl
           shell(spi)%j = kj
           shell(spi)%mj = km
           write(4, *) shell(spi)%n, shell(spi)%l, shell(spi)%j, shell(spi)%mj
        end do
     !end do
  end do

  write(4, *)
  write(4, *) 'Total number of built sp states:', spi
  write(4, *)
  write(4, *) 'END'
  close(4)
  
end subroutine build_sp_states


! ==============================================================================
! BUILD BASIS OF SLATER DETERMINANTS
! ==============================================================================
subroutine build_basis

  use global
  use input_DORA
  
  implicit none

  integer, dimension (NN) :: slater     ! Slater determinant.
  character (*), parameter :: fmt = '(i0' // repeat(', 1x, i0', NN - 1)//')'
  integer :: k,e, ib, jb,row,col,f,b
  e = 1

  ! Check if 'basis' is allocated:
  if ( allocated(basis) ) then
     write(6, *) " ERROR in subroutine 'build_basis': Array 'basis' is already allocated. Exiting"
     stop
  else
     allocate( basis(dim, NN) )
     allocate( basis_pair(dim, NN) )
  end if
  open(unit=1, status='unknown', file='basis.txt')     ! Basis File
  open(unit=2, status='unknown', file='basis_pair.txt')     ! Pairing Basis File
  
  dim = 0
  call gen(1)
  
  close(1)
  close(2)
  !do f=1,6
  !   write(6, '(1000I14)') (int(basis_pair(b,f)), b=1,4)
  !end do

contains
! ------------------------------------------------------------------------------
  recursive subroutine gen (c)

    implicit none
    integer, intent (in) :: c
    integer :: a, b, b2, f
    real :: mtot
    if (c > NN) then
       !write(*, fmt) slater
       !Test for M
       mtot = 0e0
       do b = 1, NN
          mtot = mtot + shell(slater(b))%mj     ! 'slater(b)' is the single-particle state on which the particle labeled 'b' is sitting.
       end do
       ! Check if 'mtot' is compatible with the allowed value 'M'
       if ( mtot .eq. M ) then
          !dim = dim + 1     ! Uncomment for general case (and comment the same statement in the pairing part).
          ! write(6, fmt) slater
          do b = 1, NN
             basis(b, c) = slater(b)
          end do
          write(1, '(1000I14)') (int(basis(b,c)), b=1,NN)
         !Pairing: Check if n is the same
          if (shell(slater(1))%n .eq. shell(slater(2))%n .and. shell(slater(3))%n .eq. shell(slater(4))%n ) then
             dim = dim + 1
             do b = 1, NN
                basis_pair(b, c) = slater(b)
             end do
             write(2, '(1000I14)') (int(basis_pair(b,c)), b=1,NN)
          end if
       end if
    else
       do a = 1, mnsp
          if ((c .eq. 1) .or. (a .gt. slater(c - 1))) then
             slater(c) = a
             call gen(c + 1)
          end if
       end do
    end if
 !   do f=1,d
 !      write(6, '(1000I14)') (int(basis(b,f)), b=1,NN)
 !   end do
  end subroutine gen

end subroutine build_basis


! ==============================================================================
! BUILD HAMILTONIAN
! ==============================================================================
subroutine build_hamiltonian

  use global
  use input_DORA
  implicit none
  integer, dimension(dim, NN) :: slater_ham     ! Slater determinant: this should be removed once we figure out how not to need a file 'basis_pair.txt'.
  integer, dimension(2) :: phasedMatch     ! Array containing the index of a Slater determinant 
  real, dimension(dim,dim) :: Hmat     ! Hamiltonian matrix.
  real, dimension(:, :), allocatable :: ME     ! Matrix elements of the interaction.
  integer :: d, r, nme, matchedSPState, k
  real :: sum, spse, matchedEnergy, s, t, u, v
  
  open(unit=3, status='unknown', file='basis_pair.txt')     ! Read Basis File
  open(unit=4, status='unknown', file='pairing.int')     ! Read Basis File
  open(unit=7, status='unknown', file='ham.txt')     ! Hmat File
  
  ! Read the SDs in the basis
  do d = 1, dim
     read(3, *) ( slater_ham(d, r), r=1,NN )
  end do
  
  write(6,*) dim

  ! Read the interaction files
  read(4, *) nme
  allocate( ME(nme, 5) )
  do r = 1, nme
     read(4, *) ( ME(r, k), k=1,5 )
  end do
  
  ! Initialize hamiltonian matrix
  do d = 1, dim
     do r = 1, dim
        Hmat(d,r) = 0
     end do
!     write(6,*) ( Hmat(d,r), r=1,dim )
  end do

  ! Build hamiltonian
  !! Loop over all the Slater determinants
  do d = 1, dim
     sum = 0e0
     do k = 1, nme
        phasedMatch = tbev(ME(k,1), ME(k,2), ME(k,3), ME(k,4), d)
        Hmat(d, phasedMatch(1)) = Hmat(d, phasedMatch(1)) + phasedMatch(2)*ME(k,5)
     end do
  end do
  do d = 1, dim
     spse = 0e0
     do r = 1, NN
        matchedSPState = slater_ham(d,r)
        matchedEnergy = shell(matchedSPState)%spe
        spse = spse + matchedEnergy
     end do
     write(6, *) spse
     Hmat(d, d) = Hmat(d, d) + spse
  end do



  ! Print the hamiltonian to a file
  do d=1, dim
     write(7, '(1000F14.7)') (real(Hmat(d,r)), r=1,dim)
  end do

  close(3)
  close(4)
  close(7)

contains
! ------------------------------------------------------------------------------
! Auxiliary function to calculate Two-Body Expectation Values
! ------------------------------------------------------------------------------
  function tbev(s, t, u, v, d)

    use global

    integer :: d
    real :: s, t, u, v
    integer, dimension(2) :: tbev

    do r = 1, NN
       if ( u .ne. slater_ham(d,r)) then
          
       else ( v .ne. slater_ham(d,r) ) then
       end if
    end do

  end function tbev

end subroutine build_hamiltonian
