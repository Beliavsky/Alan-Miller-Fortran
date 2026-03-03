MODULE common_xpk3
IMPLICIT NONE

! COMMON/DATA/ xdat(200),xdat(200),sigma,ndata
REAL, SAVE     :: xdat(200), ydat(200), sigma
INTEGER, SAVE  :: ndata

END MODULE common_xpk3



! Code converted using TO_F90 by Alan Miller
! Date: 2001-07-09  Time: 22:29:04
 
PROGRAM xpk3
!============================================================
!     Driver program for circle fitting problem using
!     a robust estimator based on the Hough transform
!     (Sect. 5.5)
!============================================================

USE Genetic_Algorithm
IMPLICIT NONE
INTEGER             :: seed, STATUS
INTEGER, PARAMETER  :: n=2
REAL                :: ctrl(12), x(n), f

! EXTERNAL  fit3

INTERFACE
  FUNCTION fit3(n, x) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: n
    REAL, INTENT(IN)     :: x(:)
    REAL                 :: fn_val
  END FUNCTION fit3
END INTERFACE

!     First, initialize the random-number generator

seed = 123456
CALL rninit(seed)

!     Read in synthetic data

CALL finit()

!     Set control variables (use

ctrl(1:12) = -1

!     Now call pikaia
CALL pikaia(fit3,n,ctrl,x,f,STATUS)

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
!     Reads in synthetic dataset
!===============================================================

USE common_xpk3
! COMMON/DATA/ xdat(200),ydat(200),sigma,ndata

OPEN(1, FILE='syndat3.txt', STATUS='OLD')

READ(1, *) ndata
READ(1, *) xdat(1:ndata)
READ(1, *) ydat(1:ndata)

!     Same error estimate in x and y for all data points
sigma=0.05

RETURN
END SUBROUTINE finit

END PROGRAM xpk3



FUNCTION fit3(n,x) RESULT(fn_val)
!=============================================================
!     Fitness function for circle fitting problem using
!     a robust estimator based on the Hough transform
!     (Sect. 5.5)
!=============================================================

USE common_xpk3
! COMMON/DATA/ xdat(200),ydat(200),sigma,ndata

IMPLICIT NONE
INTEGER, INTENT(IN)  :: n
REAL, INTENT(IN)     :: x(:)
REAL                 :: fn_val

INTEGER  :: i
REAL     :: x0, y0, distdat, r1, r2, sum

!---------- 1. rescale input variables:
!           0 <= x(1,2) <= 3
x0    = x(1)*3.
y0    = x(2)*3.
!---------- 2. compute merit function
r1=1.-sigma/2.
r2=1.+sigma/2.
sum=0.
DO  i=1,ndata
!     (a) compute distance d_j from center to data point
  distdat=SQRT((x0-xdat(i))**2 + (y0-ydat(i))**2)
!     (b) increment Hough merit function
  IF(distdat >= r1 .AND. distdat <= r2) sum=sum + 1.
END DO
!---------- 3. equate fitness to Hough merit function
fn_val=sum

RETURN
END FUNCTION fit3


