PROGRAM Rosenbrock

! Sample file for Rosenbrock function
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-06-10  Time: 23:36:49

USE Shor_minimization
IMPLICIT NONE

REAL (dp), DIMENSION(:), allocatable  :: x
REAL (dp)  :: options(13), f
INTEGER    :: n
LOGICAL    :: flg, flfc, flgc

INTERFACE
  SUBROUTINE rosenbf(x, f)
    IMPLICIT NONE
    INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)   :: x(:)
    REAL (dp), INTENT(OUT)  :: f
  END SUBROUTINE rosenbf

  SUBROUTINE rosenbg(x, g)
    IMPLICIT NONE
    INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)   :: x(:)
    REAL (dp), INTENT(OUT)  :: g(:)
  END SUBROUTINE rosenbg
END INTERFACE

! Flag for the use of analytically calculated gradients:
flg=.true.

! Flags for constraints are set to 0:
flfc=.false.
flgc=.false.

! Allocate array X:
n=2
allocate (x(n))

! Use the default optional parameters:
CALL soptions(options)

! Starting point:
x(1)=-1.2_dp
x(2)=1.0_dp

! Call the solver:
CALL solvopt(n, x, f, rosenbf, flg, rosenbg, options, flfc, flgc)

! Display the results:
WRITE (*, '(//a)') '     Function Value ===== Evaluations == Iterations'
WRITE (*, '(/ " ", g22.15, "  ", i5, " +", i5, "     ", i5)')  &
          f, NINT(options(10)), NINT(options(11)), NINT(options(9))
WRITE (*, '(//" Optimum Point X:" // ((" ", g22.15)))') x

deallocate (x)
STOP
END PROGRAM Rosenbrock



SUBROUTINE rosenbf(x, f)
! ROSENBF returns the function value at a point
! for Rosenbrock Saddle Function

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(OUT)  :: f

REAL (dp), PARAMETER  :: one = 1.0_dp, hundr = 100.0_dp

f=hundr*(x(2)-x(1)**2)**2 + (one-x(1))**2

RETURN
END SUBROUTINE rosenbf



SUBROUTINE rosenbg(x, g)
! ROSENBG returns the gradient at a point
! for Rosenbrock Saddle Function

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(OUT)  :: g(:)

REAL (dp), PARAMETER  :: fhund = 400.0_dp, thund = 200.0_dp, one = 1.0_dp,  &
                         two = 2.0_dp

g(1)=-fhund*(x(2)-x(1)**2)*x(1) - two*(one-x(1))
g(2)= thund*(x(2)-x(1)**2)

RETURN
END SUBROUTINE rosenbg
