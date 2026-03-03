MODULE common_shor
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

! COMMON/shorab/a,b
REAL (dp), SAVE  :: a(5,10), b(10)

END MODULE common_shor



PROGRAM ShorTest
! Sample file for Shor function
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-06-10  Time: 23:36:50

USE Shor_minimization
IMPLICIT NONE

REAL (dp), DIMENSION(:), allocatable :: x
REAL (dp)  :: options(13), f
INTEGER    :: n
LOGICAL    :: flg, flfc, flgc

! EXTERNAL shorf, shorg, null

INTERFACE
  SUBROUTINE shorf(x, f)
    IMPLICIT NONE
    INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)   :: x(:)
    REAL (dp), INTENT(OUT)  :: f
  END SUBROUTINE shorf

  SUBROUTINE shorg(x, g)
    IMPLICIT NONE
    INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)   :: x(:)
    REAL (dp), INTENT(OUT)  :: g(:)
  END SUBROUTINE shorg
END INTERFACE

! Flag for the use of analytically calculated gradients:
flg=.true.

! Flags for constraints are set to 0:
flfc=.false.
flgc=.false.

! Allocate array X:
n=5
allocate (x(n))

! Initialize the problem constants:
CALL shorinit(x)

! Use the default optional parameters:
CALL soptions(options)

! Call the solver:
CALL solvopt(n, x, f, shorf, flg, shorg, options, flfc, flgc)

! Display the results:
WRITE (*, '(//a)') '     Function Value ===== Evaluations == Iterations'
WRITE (*, '(/ " ", g22.15, "  ", i5, " +", i5, "     ", i5)')  &
          f, NINT(options(10)), NINT(options(11)), NINT(options(9))
WRITE (*, '(//" Optimum Point X:" // ((" ", g22.15)))') x

deallocate (x)
STOP


CONTAINS


!  Initialization routine:

SUBROUTINE shorinit(x)
USE common_shor
IMPLICIT NONE

REAL (dp), INTENT(OUT)  :: x(5)

a(1,1)=0.0_dp
a(1,2)=2.0_dp
a(1,3)=1.0_dp
a(1,4)=1.0_dp
a(1,5)=3.0_dp
a(1,6)=0.0_dp
a(1,7)=1.0_dp
a(1,8)=1.0_dp
a(1,9)=0.0_dp
a(1,10)=1.0_dp
a(2,1)=0.0_dp
a(2,2)=1.0_dp
a(2,3)=2.0_dp
a(2,4)=4.0_dp
a(2,5)=2.0_dp
a(2,6)=2.0_dp
a(2,7)=1.0_dp
a(2,8)=0.0_dp
a(2,9)=0.0_dp
a(2,10)=1.0_dp
a(3,1)=0.0_dp
a(3,2)=1.0_dp
a(3,3)=1.0_dp
a(3,4)=1.0_dp
a(3,5)=1.0_dp
a(3,6)=1.0_dp
a(3,7)=1.0_dp
a(3,8)=1.0_dp
a(3,9)=2.0_dp
a(3,10)=2.0_dp
a(4,1)=0.0_dp
a(4,2)=1.0_dp
a(4,3)=1.0_dp
a(4,4)=2.0_dp
a(4,5)=0.0_dp
a(4,6)=0.0_dp
a(4,7)=1.0_dp
a(4,8)=2.0_dp
a(4,9)=1.0_dp
a(4,10)=0.0_dp
a(5,1)=0.0_dp
a(5,2)=3.0_dp
a(5,3)=2.0_dp
a(5,4)=2.0_dp
a(5,5)=1.0_dp
a(5,6)=1.0_dp
a(5,7)=1.0_dp
a(5,8)=1.0_dp
a(5,9)=0.0_dp
a(5,10)=0.0_dp

b(1)=1.0_dp
b(2)=5.0_dp
b(3)=10.0_dp
b(4)=2.0_dp
b(5)=4.0_dp
b(6)=3.0_dp
b(7)=1.7_dp
b(8)=2.5_dp
b(9)=6.0_dp
b(10)=4.5_dp

x(1)=-1.0_dp
x(2)= 1.0_dp
x(3)=-1.0_dp
x(4)= 1.0_dp
x(5)=-1.0_dp

RETURN
END SUBROUTINE shorinit

END PROGRAM ShorTest



SUBROUTINE shorf(x, f)
! Function SHORF returns the function value at a point
! for Shor's piece-wise quadratic function

USE common_shor
IMPLICIT NONE

REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(OUT)  :: f

REAL (dp)  :: d, s
INTEGER    :: i, j

f=0.0_dp
DO i=1,10
  s=0.0_dp
  DO j=1,5
    s=s + (x(j) - a(j,i))**2
  END DO
  d=b(i)*s
  IF (d > f) f=d
END DO

RETURN
END SUBROUTINE shorf



SUBROUTINE shorg(x, g)
! Function SHORG returns the gradient vector
! of Shor's piece-wise quadratic function at a point

USE common_shor
IMPLICIT NONE

REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(OUT)  :: g(:)

REAL (dp) :: f, d, s
INTEGER   :: i, k

f=0.0_dp
DO i=1,10
  s=SUM( (x(1:5) - a(1:5,i))**2 )
  d=b(i)*s
  IF (d > f) THEN
    f=d
    k=i
  END IF
END DO

g(1:5)=b(k)*0.5_dp*(x(1:5) - a(1:5,k))

RETURN
END SUBROUTINE shorg
