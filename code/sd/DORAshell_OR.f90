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
  !include 'param_DORA.f90'

  integer :: dim     ! Dimension of the basis (it is determined in subroutine 'build_basis').
  integer, dimension (:, :), allocatable :: basis     ! Array in which each row is a Slater determinant (therefore its first dimension is the dimension of the basis).
  real, dimension(:,:),allocatable :: states

end module global


! ==============================================================================
! MAIN PROGRAM
! ==============================================================================
program DORA

  implicit none
  
  ! Test for input file
  call read_input

  ! Build the basis of Slater determinants
  call build_basis

  ! Build Hamiltonian
  call build_hamiltonian

end program DORA


! ==============================================================================
! READ INPUT FILE
! ==============================================================================
subroutine read_input
  
  use input_DORA
  use global

  implicit none

  ! Test
  if ( mod(NN+int(M),2) .eq. 1 ) then
     write(6, *) " ERROR in subroutine 'input_DORA': Inconsistent M input. Exiting."
     stop
  end if
  write(6, *) 'N= ', NN


end subroutine read_input

! ==============================================================================
! BUILD BASIS OF SLATER DETERMINANTS
! ==============================================================================
subroutine build_basis

  use global
  use input_DORA
  
  implicit none

  integer, dimension (NN) :: slater     ! Slater determinant.
  character (*), parameter :: fmt = '(i0' // repeat(', 1x, i0', NN - 1)//')'
  integer :: k,e, ib, jb,row,col,f,b, nsp
  e = 1
  col=6
  !!!Input States
  open(unit=1, status='unknown',file='sdshell_states.int')
  read(1,*) nsp
  allocate( states(nsp,col))
  do k=1,nsp
     read(1,*)(states(k,ib), ib=1,6)
     !write(6,*) states(k,5)
  end do
  
  ! Check if 'basis' is allocated:
  if ( allocated(basis) ) then
     write(6, *) " ERROR in subroutine 'build_basis': Array 'basis' is already allocated. Exiting"
     stop!
  else
     allocate( basis(dim, NN) )
  end if
  open(unit=1, status='unknown', file='basis.txt')     ! Basis File
  
  dim = 0!
  call gen(1)
  
  close(1)

  do f=1,6
     !write(6, '(1000I14)') (int(basis(b,f)), b=1,4)
  end do

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
          mtot = mtot + states(slater(b),5)     ! 'slater(b)' is the single-particle state on which the particle labeled 'b' is sitting.
          !write(6,*) mtot
       end do
       !write(6,*) mtot
       ! Check if 'mtot' is compatible with the allowed value 'M'
       if ( mtot .eq. M ) then
          dim = dim + 1     ! Uncomment for general case (and comment the same statement in the pairing part).
          do b = 1, NN
             basis(b, c) = slater(b)
          end do
          write(1, '(1000I14)') (basis(b,c), b=1,NN)
         
       end if
    else
       do a = 1, nsp
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
  integer :: d, r, nme, matchedSPState, k,dd
  real :: sum, spse, matchedEnergy, s, t, u, v

  !!For Diagonalization
  integer :: diaN, LDA, LWMAX, INFO, LWORK
  real, dimension(dim) :: W
  real, dimension(3*dim) :: WORK
  diaN = dim
  LDA = dim
  LWORK = 5*diaN

  
  open(unit=3, status='unknown', file='basis.txt')     ! Read Basis File
  open(unit=4, status='unknown', file='sdshell.int')     ! Read Basis File
  open(unit=7, status='unknown', file='ham_O.txt')     ! Hmat File
  
  ! Read the SDs in the basis
  do d = 1, dim
     read(3, *) ( slater_ham(d, r), r=1,NN )
     !WRITE(6, *) ( slater_ham(d, r), r=1,NN )
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
        call tbev(ME(k,1), ME(k,2), ME(k,3), ME(k,4),d, slater_ham, phasedMatch)
        !write(6,* ) d, phasedmatch(1)
        !phasedMatch = tbev(
        Hmat(d, phasedMatch(1)) = Hmat(d, phasedMatch(1)) + phasedMatch(2)*ME(k,5)
     end do
  end do
  do d = 1, dim
     do k = d+1, dim
        Hmat(d,k)=Hmat(k,d)
     end do
  end do
  do d = 1, dim
     spse = 0e0
     do r = 1, NN
        matchedSPState = slater_ham(d,r)
        matchedEnergy = states(r,6)
        spse = spse + matchedEnergy
     end do
     !write(6, *) spse
     Hmat(d, d) = Hmat(d, d) + spse
  end do



  ! Print the hamiltonian to a file
  do d=1, dim
     write(7, '(1000F14.7)') (real(Hmat(d,r)), r=1,dim)
  end do

!!!!!!!!!!!!!!!!!!Diagonalize!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call SSYEV('N','U', diaN, Hmat, LDA, W, WORK, LWORK, INFO)

  write(6,*) INFO
  write(6,*) W
  !!Eigenvalues
  !call print_matrix('Eigenvalues',1,dim,W,1) 

  
  close(3)
  close(4)
  close(7)
end subroutine build_hamiltonian


subroutine print_matrix(DESC,M,dim,Hmat,LDA)
  character*(*) :: DESC
  integer :: M, dim, LDA
  real, dimension(dim,dim) :: Hmat
  integer :: i,j

  write(*,*)
  write(*,*) DESC
  do i=1,M
     write(*,*) (Hmat(i,j),j=1, dim)
  end do

  return

end subroutine print_matrix


! ------------------------------------------------------------------------------
! Auxiliary function to calculate Two-Body Expectation Values
! ------------------------------------------------------------------------------
subroutine tbev(s, t, u, v, dd, slater_ham, phasedmatch)

  use global
  implicit none

  integer :: t1, loc, r, phase, r2,dd
  real :: s, t, u, v
  logical :: inslate1, inslate2, find
  integer, dimension(2) :: phasedmatch
  integer, dimension(dim,NN) :: slater_ham
  integer, dimension(NN) :: test, test2
  inslate1 = .false.
  inslate2 = .false.
  find = .true.

  test2 = slater_ham(dd,:)
  !write(6,*) test2
  phase = 1
  !do r=1, dim
  !write(6,*)(slater_ham(r,r2),r2=1,NN)
  !end do
  do r = 1, NN
     if ( u .eq. test2(r)) then
        inslate1 =  .true.
        test2(r) = s
     end if
  end do
  if (inslate1 .eqv. .false.) then
     phasedmatch(1) = 0
     phasedmatch(2) = 0
  else
     do r = 1, NN
        if ( v .eq. test2(r) ) then
           inslate2 = .true.
           !write(6,*) test2(r)
           test2(r)= t
           !write(6,*) test2(r)
        end if
     end do
     !write(6,*) test2
     if (inslate2 .eqv. .false.) then
        phasedmatch(1) = 0
        phasedmatch(2) = 0
     else
        call quicksort(test2,1,NN,phase)
        do r = 1, NN-1
           if (test2(r) .eq. test2(r+1)) then
              phasedmatch(1)=0
              phasedmatch(2)=0
           end if
        end do
     end if
     !write(6,*) test2
     do r = 1, dim
        find = .true.
        !write(6,*)(slater_ham(r,r2),r2=1,NN)
        do r2=1, NN
           if (test2(r2) .ne. slater_ham(r,r2)) then
              !write(6,*) dim
              find = .false.
           end if
        end do
        if (find .eqv. .true.) then
           phasedmatch(1) = r
           phasedmatch(2) = phase
           !write(6,*) r
        end if
     end do
  end if

  ! phasedmatch(1) = 1
  ! phasedmatch(2) = phase 
  !write(6,*) phase
end subroutine tbev


recursive subroutine quicksort(test2, first, last, phase)
  use global
  implicit none 

  integer, dimension(NN) :: test2
  integer :: first, last, i, j, phase, x, t

  x = test2((first+last)/2)
  i = first
  j = last
  do
     do while (test2(i) < x)
        i=i+1
     end do
     do while (x<test2(j))
        j=j-1
     end do
     if (i>=j) exit
     t = test2(i); test2(i)=test2(j); test2(j)=t
     i=i+1
     j=j-1

     !write(6,*)phase
  end do
  phase = phase
  if (first < i-1) call quicksort(test2, first, i-1, phase)
  if (first < i-1) call quicksort(test2, first, i-1, phase)
  !write(6,*)phase
end subroutine quicksort