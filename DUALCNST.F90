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
INTEGER, SAVE    :: im

END MODULE Shell_Data



PROGRAM Shell_Dual

! Code converted using TO_F90 by Alan Miller
! Date: 2000-06-10  Time: 23:36:49

! Main program for solving Shell Dual Problem
! by use of the exact penalty function method.
! The penalty coefficient is chosen automatically
! by the solver

USE Shor_minimization
IMPLICIT NONE

REAL (dp), DIMENSION(:), allocatable :: x
REAL (dp)  :: options(13), f
INTEGER    :: n
LOGICAL    :: flg,flfc,flgc

! EXTERNAL dualshel,dualobjf,dualobjg,dualcntf,dualcntg

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

  SUBROUTINE dualcntf(x, fc)
    IMPLICIT NONE
    INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)   :: x(:)
    REAL (dp), INTENT(OUT)  :: fc
  END SUBROUTINE dualcntf

  SUBROUTINE dualcntg(x, gc)
    IMPLICIT NONE
    INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)   :: x(:)
    REAL (dp), INTENT(OUT)  :: gc(:)
  END SUBROUTINE dualcntg
END INTERFACE

! All flags are set to 1:
flg=.true.
flfc=.true.
flgc=.true.

! Allocate array X:
n=15
allocate (x(n))

! Initialize the problem constants:
CALL dualshel(x)

! Use the default optional parameters:
CALL soptions(options)

! Call the solver:
CALL solvopt(n, x, f, dualobjf, flg, dualobjg, options, flfc, flgc,  &
             func=dualcntf, gradc=dualcntg)

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
x(12)=60.0_dp

RETURN
END SUBROUTINE dualshel

END PROGRAM Shell_Dual



SUBROUTINE dualobjf(x, f)
! Objective function

USE Shell_Data
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(OUT)  :: f

INTEGER  :: i, j

f=0.0_dp
DO i=1,5
  DO j=1,5
    f=f + c(i,j)*x(i)*x(j)
  END DO
  f=f + d(i)*x(i)**3*2.0_dp
END DO
DO i=6,15
  f=f - b(i-5)*x(i)
END DO

RETURN
END SUBROUTINE dualobjf



SUBROUTINE dualobjg(x, g)
! Gradient of the objective function

USE Shell_Data
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(OUT)  :: g(:)

INTEGER    :: i, j
REAL (dp)  :: s

DO i=1,15
  IF(i <= 5) THEN
    g(i)=6.0_dp*d(i)*x(i)**2
    s=0.0_dp
    DO j=1,5
      s=s + c(i,j)*x(j)
    END DO
    g(i)=g(i) + 2.0_dp*s
  ELSE
    g(i)=-b(i-5)
  END IF
END DO

RETURN
END SUBROUTINE dualobjg



SUBROUTINE dualcntf(x, fc)
! The MAXIMAL RESIDUAL for the set of all constraints:

USE Shell_Data
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(OUT)  :: fc

INTEGER    :: i, j
REAL (dp)  :: s

fc=0.0_dp
im=0
DO i=1,20
  IF(i <= 5) THEN
! The first 5 nonlinear inequalities:
    s=0.0_dp
    DO  j=1,10
      s=s + a(j,i)*x(j+5)
    END DO
    DO  j=1,5
      s=s - 2.0_dp*c(j,i)*x(j)
    END DO
    s=s - e(i) - 3.0_dp*d(i)*x(i)**2
  ELSE

! Nonnegative values:
    s=-x(i-5)
  END IF
  IF (s > fc) THEN
    fc=s
    im=i
  END IF
END DO

RETURN
END SUBROUTINE dualcntf



SUBROUTINE dualcntg(x, gc)
! Gradient of the constraint with the MAXIMAL RESIDUAL:

USE Shell_Data
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(OUT)  :: gc(:)

INTEGER    :: j

gc(1:15)=0.0_dp

IF (im > 0) THEN
  IF(im <= 5) THEN
! The first 5 nonlinear inequalities:
    DO j=1,15
      IF(j <= 5)  THEN
        gc(j)=-2.0_dp*c(j,im)
      ELSE
        gc(j)=a(j-5,im)
      END IF
    END DO
    gc(im)=gc(im) - 6.0_dp*d(im)*x(im)
  ELSE
! Nonnegative values:
    gc(im-5)=-1.0_dp
  END IF
END IF

RETURN
END SUBROUTINE dualcntg
