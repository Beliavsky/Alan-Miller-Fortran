MODULE common_xpk1b
IMPLICIT NONE

! COMMON/DATA/ f(200),t(200),sigma,ndata,m
REAL, SAVE     :: f(200), t(200), sigma
INTEGER, SAVE  :: ndata, m

END MODULE common_xpk1b



PROGRAM xpk1b
 
! Code converted using TO_F90 by Alan Miller
! Date: 2001-07-09  Time: 22:28:53
 
!==================================================================
!     Driver program for non-linear least-squares problem
!     (Sect. 5.3)
!==================================================================

USE Genetic_Algorithm
IMPLICIT  NONE
INTEGER             :: seed, STATUS
INTEGER, PARAMETER  :: n=17
REAL                :: ctrl(12), x(n), f

! EXTERNAL  fit1b

INTERFACE
  FUNCTION fit1b(n, x) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: n
    REAL, INTENT(IN)     :: x(:)
    REAL                 :: fn_val
  END FUNCTION fit1b
END INTERFACE

!     First, initialize the random-number generator

seed=13579
CALL rninit(seed)

!     Initializations

CALL finit()

!     Set control variables (use defaults except for population size)

ctrl(1:12) = -1
ctrl(1)=50

!     Now call pikaia
CALL pikaia(fit1b,n,ctrl,x,f,STATUS)

!     Print the results
WRITE(*,*) ' status: ', STATUS
WRITE(*,*) '      x: ', x
WRITE(*,*) '      f: ', f
WRITE(*,20) ctrl
STOP

20 FORMAT(   '    ctrl: ', 6F11.6/ t11, 6F11.6)


CONTAINS

!***************************************************************

SUBROUTINE finit()
!===============================================================
!     Reads in synthetic dataset (see Figure 5.4)
!===============================================================

USE common_xpk1b
! COMMON/DATA/ f(200),t(200),sigma,ndata,m

REAL     :: f0(200), vdum(11)

OPEN(1, FILE='syndat1.txt', STATUS='OLD')

READ(1, *) ndata
READ(1, *) vdum(1:11)
READ(1, *) t(1:ndata)
READ(1, *) f0(1:ndata)
READ(1, *) f(1:ndata)

!     Use 5 Fourier modes for the fit
m=5

!     same error bar for all point
sigma=5.

RETURN
END SUBROUTINE finit

END PROGRAM xpk1b


 
FUNCTION fit1b(n,x) RESULT(fn_val)
!==========================================================
!     Fitness function for non-linear least squares problem
!     (Sect. 5.3)
!==========================================================

! COMMON/DATA/ f(200),t(200),sigma,ndata,m
USE common_xpk1b
IMPLICIT NONE
INTEGER, INTENT(IN)  :: n
REAL, INTENT(IN)     :: x(:)
REAL                 :: fn_val

INTEGER, PARAMETER :: mmax=10
INTEGER            :: i, j
REAL               :: amp(mmax), per(mmax), a, b, sum, sum2, nyp
REAL, PARAMETER    :: pi=3.1516926536

!---------- 1. rescale input variables:
a=x(1)*10.
b=x(2)*100.
nyp=2.*(t(2)-t(1))
DO  j=1,m
  amp(j)=x(3*j)*100.
  per(j)=x(3*j+1)*(50. - nyp) + nyp
END DO
!---------- 2. compute chi**2
sum=0.
DO  i=1,ndata
  sum2=0.
  DO  j=1,m
    sum2=sum2 + amp(j)*SIN(2.*pi*(t(i)/per(j) + x(3*j+2)))
  END DO
  sum=sum + ( (a*t(i)+b+sum2-f(i))/sigma)**2
END DO
!---------- 3. define fitness
fn_val=1./sum

RETURN
END FUNCTION fit1b
