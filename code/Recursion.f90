program recursion
  implicit none

  integer, parameter :: N = 3
  integer, parameter :: mnsp = 5
  integer, dimension (N) :: comb
  character (*), parameter :: fmt = '(i0' // repeat(', 1x, i0', N - 1)//')'

  call gen(1)

contains

  recursive subroutine gen (m)

    implicit none
    integer, intent (in) :: m
    integer :: a

    if (m > N) then
       write(*, fmt) comb
    else
       do a = 1, mnsp
          if ((m .eq. 1) .or. (a .gt. comb (m - 1))) then
             comb (m) = a
             call gen(m + 1)
          end if
       end do
    end if
  end subroutine gen
  
end program recursion
