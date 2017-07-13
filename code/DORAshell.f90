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

  type sp     ! N.B.: 'structure' and 'record' statements are pre-Fortran 90 statements that are equivalent to 'type' statements, but they are obsolete (not portable). Make sure to always use 'type'.
     integer :: n     ! Principal quantum number.
     integer :: l     ! Orbital momentum.
     real :: j     ! Total angular momentum.
     real :: mj     ! Total angular momentum projection.
  end type sp
  type(sp), dimension(mnsp) :: shell     ! One-dimensional array of sets of quantum numbers (single-particle states).
!!$  type(sp), dimension(:), allocatable :: shell     ! One-dimensional array of sets of quantum numbers (single-particle states).

!!$  type determinant
!!$     integer, dimension (A-Z) :: det
!!$  end type determinant
  integer, dimension (A-Z) :: slater
!!$  type(determinant) :: slater
  integer, dimension (:, :), allocatable :: basis

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
  !call build_hamiltonian

  ! Solve eigenproblem
  !call eigensolver

end program DORA


! ==============================================================================
! READ INPUT FILE
! ==============================================================================
subroutine read_input
  
  use input_DORA

  implicit none

!  logical :: file_exists
  integer :: N

!  inquire(file='', exist=file_exists)

  N = A-Z
  ! Test
  if ( mod(N+int(2*M),2) .eq. 1 ) then
     write(6, *) " ERROR in subroutine 'input_DORA': Inconsistent M input. Exiting."
     stop
  end if
  write(6, *) 'Z= ', Z
  write(6, *) 'N= ', N
  write(6, *) 'M= ', M

end subroutine read_input



! ==============================================================================
! BUILD SINGLE-PARTICLE STATES AND STORE THEM IN A FILE
! ==============================================================================
subroutine build_sp_states

  use global

  implicit none

  logical :: flag
  integer :: kn, kl, spi, im     ! spi: single-particle state index in the 'shell' array.
  real :: kj, km

  open(unit=4, status='unknown', file='sp_states.dat')     ! File where all the single-particle states are stored.

  ! Write header in 'sp_states.dat'
  write(4, *) 'SINGLE-PARTICLE STATES BUILT BY THE D.O.R.A. SHELL-MODEL CODE'
  write(4, *)
  write(4, *) '          n          l            j          m_j'
  
  ! Loop over principal quantum numbers
  kn = -1
  spi = 0
  flag = .false.
  do
     if ( flag .eqv. .true. ) then     ! If verified, then the sp array is already full.
        exit
     end if
     kn = kn + 1
     do kl = 0, kn
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
           write(4, *) shell(spi)%n, shell(spi)%l, shell(spi)%j, shell(spi)%mj
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
     end do
  end do

  write(4, *)
  write(4, *) 'Total number of built sp states:', spi
  write(4, *)
  write(4, *) 'END'
  close(4)
  
end subroutine build_sp_states


! ==============================================================================
! BUILD SLATER DETERMINANTS
! ==============================================================================
subroutine build_basis

  use global
  use input_DORA
  
  implicit none

  character (*), parameter :: fmt = '(i0' // repeat(', 1x, i0', A-Z - 1)//')'
  integer :: N, k, dim
  
  N = A-Z     ! Simple variable change.

  ! Check if 'basis' is allocated:
  if ( allocated(basis) ) then
     write(6, *) " ERROR in subroutine 'build_basis': Array 'basis' is already allocated. Exiting"
     stop
  else
     allocate( basis(dim, N) )
  end if

  call gen(1)

! ------------------------------------------------------------------------------
contains

  recursive subroutine gen (c)

    implicit none
    integer, intent (in) :: c
    integer :: a, b
    real :: mtot

    if (c > N) then
       !write(*, fmt) slater
       !Test for M
       mtot = 0e0
       do b = 1, N
          mtot = mtot + shell(slater(b))%mj
       end do
       ! Check if 'mtot' is compatible with the allowed value 'M'
       if ( mtot .eq. M ) then
          write(6, fmt) slater
          do b = 1, N
             basis(c, b) = slater(b)
          end do
       end if
    else
       do a = 1, mnsp
          if ((c .eq. 1) .or. (a .gt. slater(c - 1))) then
             slater(c) = a
             call gen(c + 1)
          end if
       end do
    end if
    
  end subroutine gen

end subroutine build_basis
