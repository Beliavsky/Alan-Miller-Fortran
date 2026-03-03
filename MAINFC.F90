MODULE Common_findif
IMPLICIT NONE
INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)

! COMMON /findif/   accrcy
REAL (dp), SAVE  :: accrcy

END MODULE Common_findif



PROGRAM Difference_CNB

!***********************************************************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-14  Time: 16:13:30
 
! FINITE-DIFFERENCING, CUSTOMIZED, NO BOUNDS
!***********************************************************************
! MAIN PROGRAM TO MINIMIZE A FUNCTION (REPRESENTED BY THE ROUTINE SFUN)
! OF N VARIABLES X - CUSTOMIZED VERSION, WITH FINITE DIFFERENCING

USE Truncated_Newton
USE Common_findif, ONLY: accrcy
IMPLICIT NONE

REAL (dp)  :: x(50), f, g(50)
REAL (dp)  :: eta, xtol, stepmx, faccy
INTEGER    :: i, ierror, maxfun, maxit, msglvl, n
! COMMON /findif/   accrcy
! EXTERNAL          sfun

INTERFACE
  SUBROUTINE sfun(n, x, f, g)
    IMPLICIT NONE
    INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(IN)   :: x(:)
    REAL (dp), INTENT(OUT)  :: f, g(:)
  END SUBROUTINE sfun
END INTERFACE

! SET UP FUNCTION AND VARIABLE INFORMATION
! N -  NUMBER OF VARIABLES
! X -  INITIAL ESTIMATE OF THE SOLUTION
! F -  ROUGH ESTIMATE OF FUNCTION VALUE AT SOLUTION

n  = 10
DO  i = 1,n
  x(i) = i / REAL(n+1, KIND=dp)
END DO

! SET UP CUSTOMIZING PARAMETERS
! ETA     - SEVERITY OF THE LINESEARCH
! MAXFUN  - MAXIMUM ALLOWABLE NUMBER OF FUNCTION EVALUATIONS
! XTOL    - DESIRED ACCURACY FOR THE SOLUTION X*
! STEPMX  - MAXIMUM ALLOWABLE STEP IN THE LINESEARCH
! FACCY   - ACCURACY OF COMPUTED FUNCTION VALUE
! ACCRCY  - ACCURACY OF COMPUTED FUNCTION AND GRADIENT VALUES
!           (VALUES ADJUSTED BECAUSE OF FINITE DIFFERENCING OF GRADIENT)
! MSGLVL  - DETERMINES QUANTITY OF PRINTED OUTPUT
!           0 = NONE, 1 = ONE LINE PER MAJOR ITERATION.
! MAXIT   - MAXIMUM NUMBER OF INNER ITERATIONS PER STEP

maxit  = n/2
maxfun = 150*n
eta    = .25D0
stepmx = 1.d1
faccy  = 1.d-15
accrcy = SQRT(faccy)
xtol   = SQRT(accrcy)
msglvl = 1

! MINIMIZE THE FUNCTION

CALL lmqn (ierror, n, x, f, g, sfun,  &
           msglvl, maxit, maxfun, eta, stepmx, accrcy, xtol)

! PRINT THE RESULTS

IF (ierror /= 0) WRITE(*,800) ierror
IF (msglvl >= 1) WRITE(*,810)
IF (msglvl >= 1) WRITE(*,820) (i,x(i),i=1,n)
STOP

800 FORMAT(//' ERROR CODE =', i3,/)
810 FORMAT(t11, 'CURRENT SOLUTION IS '/ t15, 'I', t27, 'X(I)')
820 FORMAT(t11, i5, '  ', g22.15)
END PROGRAM Difference_CNB



SUBROUTINE sfun (n, x, f, g)
USE Common_findif
IMPLICIT NONE

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(OUT)  :: f
REAL (dp), INTENT(OUT)  :: g(:)

REAL (dp)  :: xh, fh, hinv, y(n)
INTEGER    :: i
! COMMON /findif/   h

! ROUTINE TO COMPUTE FUNCTION (F) AND GRADIENT (G) OF THE OBJECTIVE
! FUNCTION AT THE POINT X
! GRADIENT OBTAINED VIA FINITE-DIFFERENCING
! FUNCTION OBTAINED FROM SUBROUTINE GETF

hinv = 1.d0 / accrcy
CALL getf (n, x, f)
y = x(1:n)
DO  i = 1,n
  xh   = x(i)
  y(i) = x(i) + accrcy
  CALL getf (n, y, fh)
  g(i) = (fh - f) * hinv
  y(i) = xh
END DO
RETURN

CONTAINS

SUBROUTINE getf (n, x, f)
INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(OUT)  :: f

REAL (dp) :: t
INTEGER   :: i

! ROUTINE TO EVALUATE FUNCTION (F)

f = 0.d0
DO  i = 1,n
  t = x(i) - i
  f = f + t*t
END DO
RETURN
END SUBROUTINE getf

END SUBROUTINE sfun


