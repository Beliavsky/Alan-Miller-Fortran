PROGRAM Customized_Bounds

!***********************************************************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-14  Time: 16:13:24
 
! CUSTOMIZED, BOUNDS ON VARIABLES
!***********************************************************************
! MAIN PROGRAM TO MINIMIZE A FUNCTION (REPRESENTED BY THE ROUTINE SFUN)
! SUBJECT TO BOUNDS ON THE VARIABLES X - CUSTOMIZED VERSION

USE Truncated_Newton
IMPLICIT NONE

REAL (dp)  :: x(50), f, g(50), low(50), up(50)
REAL (dp)  :: eta, accrcy, xtol, stepmx
INTEGER    :: i, ierror, ipivot(50), maxfun, maxit, msglvl, n
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
! N   - NUMBER OF VARIABLES
! X   - INITIAL ESTIMATE OF THE SOLUTION
! LOW - LOWER BOUNDS
! UP  - UPPER BOUNDS
! F   - ROUGH ESTIMATE OF FUNCTION VALUE AT SOLUTION

n  = 10
DO  i = 1,n
  x(i)   = i / REAL(n+1, KIND=dp)
  low(i) = 0.d0
  up(i)  = 6.d0
END DO
f  = 1.d0

! SET UP CUSTOMIZING PARAMETERS
! ETA    - SEVERITY OF THE LINESEARCH
! MAXFUN - MAXIMUM ALLOWABLE NUMBER OF FUNCTION EVALUATIONS
! XTOL   - DESIRED ACCURACY FOR THE SOLUTION X*
! STEPMX - MAXIMUM ALLOWABLE STEP IN THE LINESEARCH
! ACCRCY - ACCURACY OF COMPUTED FUNCTION VALUES
! MSGLVL - DETERMINES QUANTITY OF PRINTED OUTPUT
!          0 = NONE, 1 = ONE LINE PER MAJOR ITERATION.
! MAXIT  - MAXIMUM NUMBER OF INNER ITERATIONS PER STEP

maxit  = n/2
maxfun = 150*n
eta    = .25D0
stepmx = 1.d1
accrcy = 1.d-15
xtol   = SQRT(accrcy)
msglvl = 1

! MINIMIZE THE FUNCTION

CALL lmqnbc (ierror, n, x, f, g, sfun, low, up, ipivot,  &
             msglvl, maxit, maxfun, eta, stepmx, accrcy, xtol)

! PRINT THE RESULTS

IF (ierror /= 0) WRITE(*,800) ierror
IF (msglvl >= 1) WRITE(*,810)
IF (msglvl >= 1) WRITE(*,820) (i,x(i),i=1,n)
STOP

800 FORMAT(//' ERROR CODE =', i3/)
810 FORMAT(t11, 'CURRENT SOLUTION IS '/ t15, 'I', t27, 'X(I)')
820 FORMAT(t11, i5, '  ', g22.15)
END PROGRAM Customized_Bounds



SUBROUTINE sfun (n, x, f, g)
IMPLICIT NONE
INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(OUT)  :: f
REAL (dp), INTENT(OUT)  :: g(:)

REAL (dp)  :: t
INTEGER    :: i

! ROUTINE TO EVALUATE FUNCTION (F) AND GRADIENT (G) OF THE OBJECTIVE
! FUNCTION AT THE POINT X

f = 0.d0
DO  i = 1,n
  t    = x(i) - i
  f    = f + t*t
  g(i) = 2 * t
END DO
RETURN
END SUBROUTINE sfun
