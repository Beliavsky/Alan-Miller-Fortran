MODULE common_xpk2
IMPLICIT NONE

! COMMON/DATA/ xdat(200),t(200),ndata
REAL, SAVE     :: xdat(200), ydat(200)
INTEGER, SAVE  :: ndata

END MODULE common_xpk2



! Code converted using TO_F90 by Alan Miller
! Date: 2001-07-09  Time: 22:28:59
 
PROGRAM xpk2
!==================================================================
!     Driver program for ellipse fitting problem (Sect. 5.4)
!==================================================================

USE Genetic_Algorithm
IMPLICIT NONE
INTEGER             :: seed, STATUS
INTEGER, PARAMETER  :: n=5
REAL                :: ctrl(12), x(n), f

! EXTERNAL  fit2

INTERFACE
  FUNCTION fit2(n, x) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: n
    REAL, INTENT(IN)     :: x(:)
    REAL                 :: fn_val
  END FUNCTION fit2
END INTERFACE


!     First, initialize the random-number generator

seed=654321
CALL rninit(seed)

!     Initializations

CALL finit()

!     Set control variables (50 individuals for 200 generations
!     under Select-random-delete-worst reproduction plan)

ctrl(1:12) = -1
ctrl(1)=50
ctrl(2)=200
ctrl(9)=0.0
ctrl(10)=3

!     Now call pikaia
CALL pikaia(fit2,n,ctrl,x,f,STATUS)

!     Print the results
WRITE(*,*) ' status: ', STATUS
WRITE(*,*) '      x: ', x
WRITE(*,*) '      f: ', f
WRITE(*,20) ctrl
STOP
20 FORMAT('    ctrl: ', 6F11.6 / t11, 6F11.6)


CONTAINS

!***************************************************************

SUBROUTINE finit()
!===============================================================
!     Reads in synthetic dataset (see Figure 5.9)
!===============================================================

! COMMON/DATA/ xdat(200),ydat(200),ndata
USE common_xpk2

OPEN(1, FILE='syndat2.txt', STATUS='OLD')

READ(1, *) ndata
READ(1, *) xdat(1:ndata)
READ(1, *) ydat(1:ndata)

RETURN
END SUBROUTINE finit

END PROGRAM xpk2


!***********************************************************


FUNCTION fit2(n,x) RESULT(fn_val)
!=============================================================
!     Fitness function for ellipse fitting problem (Sect. 5.4)

!     Ellipse parameters are:
!     x(1)=x_0, x(2)=y_0, x(3)=a, x(4)=b, x(5)=theta_0
!=============================================================

! COMMON/DATA/ xdat(200),ydat(200),ndata
USE common_xpk2

IMPLICIT NONE
INTEGER, INTENT(IN)  :: n
REAL, INTENT(IN)     :: x(:)
REAL                 :: fn_val

INTEGER  :: i
REAL     :: x0, y0, a2, b2, distdat, angdat, rthj, sum
REAL, PARAMETER  :: pi=3.1415926536

!---------- 1. rescale input variables:
!           0 <= x(1...4) <= 2, 0 <= x(5) <= pi
x0    = x(1)*2.
y0    = x(2)*2.
a2    =(x(3)*2.)**2
b2    =(x(4)*2.)**2

!---------- 2. compute merit function
sum=0.
DO  i=1,ndata
!     (a) compute distance d_j from center to data point
  distdat=SQRT((x0-xdat(i))**2 + (y0-ydat(i))**2)
!     (b) compute angle theta_j of segment center---data point
  angdat=fatan((ydat(i)-y0),(xdat(i)-x0)) - x(5)*pi
!     (c) compute radius r(\theta_j) of ellipse at that angle
  rthj=SQRT(a2*b2/(a2*SIN(angdat)**2 + b2*COS(angdat)**2))
!     (d) increment DR merit function
  sum=sum + (distdat-rthj)**2
END DO
!---------- 3. equate fitness to inverse of DR merit function
fn_val=1./sum

RETURN


CONTAINS


FUNCTION fatan(yy,xx) RESULT(fn_val)
!===========================================================
!     Returns arctangent in full circle
!===========================================================

REAL, INTENT(IN)  :: yy
REAL, INTENT(IN)  :: xx
REAL              :: fn_val

REAL :: a1
REAL, PARAMETER  :: pi = 3.1415926536

a1=ATAN(yy/xx)
IF(xx < 0.) a1=a1 + pi
IF(a1 < 0.) a1=2.*pi + a1
fn_val=a1

RETURN
END FUNCTION fatan

END FUNCTION fit2
