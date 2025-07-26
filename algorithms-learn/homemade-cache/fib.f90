module wrapper
    implicit none
    public

contains
    recursive subroutine fib(n, out)
        integer, intent(in) :: n
        integer, intent(out) :: out
        integer :: tmp1, tmp2
        
        !f2py intent(in) n
        !f2py intent(out) out

        if ( n .le. 1 ) then
            out = n
        else
            call fib(n - 1, tmp1)
            call fib(n - 2, tmp2)
            out = tmp1 + tmp2
        end if
    end subroutine 

end module

! program main
!     use wrapper, only: fib
!     implicit none

!     integer :: n = 30
!     integer :: out

!     call fib(n, out)

!     print *, out
! end program main