!
! 
! D.O.R.A. SHELL-MODEL CODE
!
! 
! ==============================================================================


program DORA

  implicit none

  ! Test for input file
  !  TO DO

  ! "Include" parameters (maximum values...)
  include 'parameters.f90'

  ! Open the interaction file
  ! TO DO

  ! Build single-particle states
  call build_sp_states(mnsp)

  ! Build Hamiltonian
  !call build_hamiltonian

  ! Solve eigenproblem
  !call eigensolver

end program DORA


subroutine build_sp_states(mnsp)

  implicit none

  type sp
     integer :: n     ! Principal quantum number.
     integer :: l     ! Orbital momentum.
     double precision :: j     ! Total angular momentum.
     double precision :: m     ! Total angular momentum projection.
  end type sp

  type(sp) :: shell

  integer :: kn, kl, nsp
  double precision :: kj, km

  open(unit=4, status='old', file='sp_states.dat')

  ! Write header in 'sp_states.dat'.
  ! TO dO
  
  ! Loop over principal quantum numbers
  kn = 0
  do while( nsp .eq. mnsp )
     kn = kn + 1
     do kl = 0,kn
        ! First "m" loop for j=l+1/2
        kj = dble(kl) + O.5
        do km = -kj,kj
           ! Fill the single-particle state array
           shell%n = kn
           shell%l = kl
           shell%j = kj
           shell%m = km
           write(4, *) shell%n, shell%l, shell%j, shell%m
        end do
        ! First "m" loop for j=l+1/2
        kj = dble(kl) - O.5
        do km = -kj,kj
           ! Fill the single-particle state array
           shell%n = kn
           shell%l = kl
           shell%j = kj
           shell%m = km
           write(4, *) shell%n, shell%l, shell%j, shell%m
        end do
     end do
  end do

  close(4)
  
end subroutine build_sp_states

