PROGRAM problem
USE tensolve
IMPLICIT NONE

!     TENSOLVE finds roots of systems of n nonlinear equations in n unknowns,
!     or minimizers of the sum of squares of m > n nonlinear functions in
!     n unknowns, using tensor methods.

!     This example illustrates the use of TENSOLVE to solve a nonlinear
!     equation problem defined by the subroutine frosen (included below).


! From: HP Authorized Customer <@goodnet.com>
! Subject: non-linear system solver
! Date: Sunday, 1 November 1998 7:26

! I have this system:
! (the '1's and '2's are subscripts)
! w1 + w2               = 1/2            w1 => x(1)
! w1*x1 + w2*x2         = 1/3            w2 => x(2)
! w1*(x1^2) + w2*(x2^2) = 1/4            x1 => x(3)
! w1*(x1^3) + w2*(x2^3) = 1/5.           x2 => x(4)

! Thanx, I appreciate any help.

! jeff


INTEGER            :: m, n, msg, termcd
INTEGER, PARAMETER :: maxm = 100, maxn = 30, maxp = 6
REAL (dp)          :: x0(maxn), xp(maxn), fp(maxm), gp(maxn)

INTERFACE
  SUBROUTINE frosen ( x, f, m, n )
    USE tensolve, ONLY: dp
    IMPLICIT NONE
    REAL (dp), INTENT(IN)    :: x(:)
    REAL (dp), INTENT(OUT)   :: f(:)
    INTEGER, INTENT(IN)      :: m
    INTEGER, INTENT(IN)      :: n
  END SUBROUTINE frosen
END INTERFACE

!     Set dimensions of the problem.

m      = 4
n      = 4

!     Set values for the initial point.

x0(1)  = 0.3_dp
x0(2)  = 0.2_dp
x0(3)  = 0.6_dp
x0(4)  = 0.9_dp
msg    = 0

!     Call TENSOLVE.

CALL tsnesi(maxm, maxn, maxp, x0, m, n, frosen, msg, xp, fp, gp, termcd )

WRITE(*, '(a, i5)')     ' Termination code = ', termcd
WRITE(*, '(a, 4g14.6)') ' Solution: ', xp(1:4)
WRITE(*, '(a, 4g14.6)') ' Function values at solution = ', fp(1:4)
WRITE(*, '(a, 4g14.6)') ' Gradients at solution = ', gp(1:4)

!     end of main program.
STOP
END PROGRAM problem


SUBROUTINE frosen ( x, f, m, n )
USE tensolve, ONLY: dp
IMPLICIT NONE
REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(OUT)  :: f(:)
INTEGER, INTENT(IN)     :: m
INTEGER, INTENT(IN)     :: n

!     frosen defines function values for the function.

f(1) = x(1) + x(2) - 0.5_dp
f(2) = x(1)*x(3) + x(2)*x(m) - 0.333333333333333_dp
f(3) = x(1)*x(3)**2 + x(2)*x(m)**2 - 0.25_dp
f(n) = x(1)*x(3)**3 + x(2)*x(m)**3 - 0.2_dp

!     end of frosen.
RETURN
END SUBROUTINE frosen
