MODULE Shell_Data
IMPLICIT NONE
INTEGER, PARAMETER, PRIVATE  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), SAVE  :: a(10,5) = RESHAPE(  &
    (/ -16.0_dp,  0.0_dp, -3.5_dp,  0.0_dp, 0.0_dp,  2.0_dp, -1.0_dp, -1.0_dp, &
         1.0_dp,  1.0_dp,  2.0_dp, -2.0_dp, 0.0_dp, -2.0_dp, -9.0_dp,  0.0_dp, &
        -1.0_dp, -2.0_dp,  2.0_dp,  1.0_dp, 0.0_dp,  0.0_dp,  2.0_dp,  0.0_dp, &
        -2.0_dp, -4.0_dp, -1.0_dp, -3.0_dp, 3.0_dp,  1.0_dp,  1.0_dp,  0.4_dp, &
         0.0_dp, -4.0_dp,  1.0_dp,  0.0_dp,-1.0_dp, -2.0_dp,  4.0_dp,  1.0_dp, &
         0.0_dp,  2.0_dp,  0.0_dp, -1.0_dp,-2.8_dp,  0.0_dp, -1.0_dp, -1.0_dp, &
         5.0_dp,  1.0_dp /), (/ 10, 5 /) ),  &
                    b(10) = (/ -40.0_dp, -2.0_dp,  -0.25_dp, -4.0_dp,  &
                                -4.0_dp, -1.0_dp, -40.0_dp, -60.0_dp,  &
                                 5.0_dp,  1.0_dp /),  &
                    c(5,5) = RESHAPE(     &
    (/ 30.0_dp, -20.0_dp, -10.0_dp,  32.0_dp, -10.0_dp, -20.0_dp,  39.0_dp, &
       -6.0_dp, -31.0_dp,  32.0_dp, -10.0_dp,  -6.0_dp,  10.0_dp,  -6.0_dp, &
      -10.0_dp,  32.0_dp, -31.0_dp, -6.0_dp,   39.0_dp, -20.0_dp, -10.0_dp, &
       32.0_dp, -10.0_dp, -20.0_dp, 30.0_dp /), (/ 5, 5 /) ),  &
                    d(5) = (/ 4.0_dp, 8.0_dp, 10.0_dp, 6.0_dp, 2.0_dp /),  &
                    e(5) = (/ -15.0_dp, -27.0_dp, -36.0_dp, -18.0_dp, -12.0_dp /)
! A priori given penalty coefficient:
REAL (dp), SAVE  :: penalty = 1000.0_dp
INTEGER, SAVE    :: im

END MODULE Shell_Data



PROGRAM Shell_Dual

! Code converted using TO_F90 by Alan Miller
! Date: 2000-06-10  Time: 23:36:49

! Main program for solving Shell Dual Problem
! by use of the exact penalty function method.
! The penalty coefficient is chosen a priori by the user.
! This is an example on calculating gradients by forward differences

USE Shor_minimization
IMPLICIT NONE

REAL (dp), DIMENSION(:), allocatable :: x
REAL (dp)  :: options(13), f
INTEGER    :: n
LOGICAL    :: flg,flfc,flgc

! EXTERNAL dualshel,dualobjf,dualobjg,null

INTERFACE
  SUBROUTINE dualobjf(x, f)
    IMPLICIT NONE
    INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)   :: x(:)
    REAL (dp), INTENT(OUT)  :: f
  END SUBROUTINE dualobjf

  SUBROUTINE dualobjg(x, g)
    IMPLICIT NONE
    INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)   :: x(:)
    REAL (dp), INTENT(OUT)  :: g(:)
  END SUBROUTINE dualobjg

  SUBROUTINE null(x, g)
    IMPLICIT NONE
    INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)   :: x(:)
    REAL (dp), INTENT(OUT)  :: g(:)
  END SUBROUTINE null
END INTERFACE

! Flag for user-supplied gradients is set to 0:
flg=.false.

! Flags for constraints are set to 0:
flfc=.false.
flgc=.false.

! Allocate array X:
n=15
allocate (x(n))

! Initialize the problem constants:
CALL dualshel(x)

! Use the default optional parameters:
CALL soptions(options)

! Call the solver:
CALL solvopt(n, x, f, dualobjf, flg, null, options, flfc, flgc)

! Display the results:
WRITE (*, '(//a)') '     Function Value ===== Evaluations == Iterations'
WRITE (*, '(/ " ", g22.15, "  ", i5, " +", i5, "     ", i5)')  &
          f, NINT(options(10)), NINT(options(11)), NINT(options(9))
WRITE (*, '(//" Optimum Point X:" // ((" ", g22.15)))') x

deallocate (x)
STOP


CONTAINS


!  Initialization routine:

SUBROUTINE dualshel(x)

REAL (dp), INTENT(OUT)  :: x(:)

! Standard starting point:
x(1:15)=1.d-4
x(12)  =60.0_dp

RETURN
END SUBROUTINE dualshel

END PROGRAM Shell_Dual



SUBROUTINE dualobjf(x,f)
! Penalty function:

USE Shell_Data
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(OUT)  :: f

REAL (dp) :: f0, s
INTEGER   :: i, j

! Objective function:
f0=0.0_dp
DO i=1,5
  DO j=1,5
    f0=f0 + c(i,j)*x(i)*x(j)
  END DO
  f0=f0 + d(i)*x(i)**3*2.d0
END DO
DO i=6,15
  f0=f0 - b(i-5)*x(i)
END DO

! Penalty item:
f=0.d0
im=0
DO i=1,20
  IF(i <= 5) THEN
    s=DOT_PRODUCT( a(1:10,i), x(6:15) )
    DO  j=1,5
      s=s - 2.d0*c(j,i)*x(j)
    END DO
    s=s - e(i) - 3.d0*d(i)*x(i)**2
  ELSE
    s=-x(i-5)
  END IF
  IF (s > f) THEN

! The maximal residual:
    f=s

! The constraint with the maximal residual:
    im=i
  END IF
END DO

! The sum of the above two items:
f=f0 + penalty*f

RETURN
END SUBROUTINE dualobjf



SUBROUTINE null(x, g)
IMPLICIT NONE
INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(OUT)  :: g(:)

g = 0.0_dp

RETURN
END SUBROUTINE null

