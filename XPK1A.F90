MODULE common_xpk1a
IMPLICIT NONE

! COMMON/DATA/ f(200),t(200),sigma,ndata
REAL, SAVE     :: f(200), t(200), sigma
INTEGER, SAVE  :: ndata

END MODULE common_xpk1a



PROGRAM xpk1a
 
! Code converted using TO_F90 by Alan Miller
! Date: 2001-07-09  Time: 22:28:29
 
!================================================================
!     Driver program for linear least-squares problem
!     (Sect. 5.2)
!================================================================

USE Genetic_Algorithm
IMPLICIT  NONE
INTEGER             :: seed, STATUS
INTEGER, PARAMETER  :: n=2
REAL                :: ctrl(12), x(n), f

! EXTERNAL  fit1a

INTERFACE
  FUNCTION fit1a(n, x) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: n
    REAL, INTENT(IN)     :: x(:)
    REAL                 :: fn_val
  END FUNCTION fit1a
END INTERFACE

!     First, initialize the random-number generator

seed = 123456
CALL rninit(seed)

!     Initializations

!     open(1, file='fake.i3e', form='unformatted')
CALL finit()

!     Set control variables (evolve 50 individuals over 100
!     generations, use defaults values for other input parameters)

ctrl(1:12) = -1
ctrl(1)=50
ctrl(2)=100

!     Now call pikaia

CALL pikaia(fit1a,n,ctrl,x,f,STATUS)

!     Print the results
WRITE(*,*) ' status: ', STATUS
WRITE(*,*) '      x: ', x
WRITE(*,*) '      f: ', f
WRITE(*,20) ctrl
STOP

20 FORMAT('    ctrl: ', 6F11.6 / t11, 6F11.6)

CONTAINS

!*********************************************************************

SUBROUTINE finit()

!     Reads in synthetic dataset (see Figure 5.4)

! COMMON/DATA/ f(200),t(200),sigma,ndata
USE common_xpk1a

REAL     :: f0(200), vdum(11)

OPEN(1, FILE='syndat1.txt', STATUS='OLD')

READ(1, *) ndata
READ(1, *) vdum(1:11)
READ(1, *) t(1:ndata)
READ(1, *) f0(1:ndata)
READ(1, *) f(1:ndata)
sigma=5.

RETURN
END SUBROUTINE finit

END PROGRAM xpk1a



FUNCTION fit1a(n,x) RESULT(fn_val)
!========================================================
!     Fitness function for linear least squares problem
!     (Sect. 5.2)
!========================================================

USE common_xpk1a
IMPLICIT NONE
INTEGER, INTENT(IN)  :: n
REAL, INTENT(IN)     :: x(:)
REAL                 :: fn_val

! COMMON/DATA/ f(200),t(200),sigma,ndata

INTEGER  :: i
REAL     :: a, b, sum
!---------- 1. rescale input variables:
a=x(1)*10.
b=x(2)*100.
!---------- 2. compute chi**2
sum=0.
DO  i=1,ndata
  sum=sum + ( (a*t(i)+b-f(i))/sigma )**2
END DO
!---------- 3. define fitness
fn_val=1./sum

RETURN
END FUNCTION fit1a


