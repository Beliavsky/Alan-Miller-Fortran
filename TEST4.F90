MODULE common_test4
USE orthpol, ONLY: dp
IMPLICIT NONE

! COMMON/s/c, a, b, e, epsma

REAL (dp), SAVE :: a(81), b(81), e(81), c, epsma

END MODULE common_test4



PROGRAM test4
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-02-04  Time: 20:39:50
! Latest revision - 7 June 2001

! This is a test of the routine  mcdis, which is applied to generate
! the first 80 recurrence coefficients for the orthogonal polynomials
! belonging to the weight function

!       (1-x**2)**(-1/2) + c   on (-1, 1),  c = 1, 10, 100.

! The corresponding inner product is discretized by applying the Gauss-
! Chebyshev quadrature rule to the first term of the weight function, and
! the Gauss-Legendre rule to the second term.  In addition to the
! beta-coefficients (all alpha's are zero), the routine prints the variables
! ncap  and  kount  to confirm convergence after one iteration.

USE orthpol
USE common_test4
IMPLICIT NONE

REAL (dp) :: xp(1), yp(1), endl(1) = 0.0_dp, endr(1) = 0.0_dp,  &
             alpha(80), beta(80), be(80), x(81), w(81), betap(80, 3)
LOGICAL   :: finl = .FALSE., finr = .FALSE.
INTEGER   :: ic, idelta, ie, iem, ierr, iq, irout, k, km1, kount, mc, mp,  &
             n, ncap, ncapm
REAL (dp) :: eps

INTERFACE
  SUBROUTINE qchle(n, x, w, i, ierr)
    USE orthpol
    USE common_test4
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(OUT)  :: x(:)
    REAL (dp), INTENT(OUT)  :: w(:)
    INTEGER, INTENT(IN)     :: i
    INTEGER, INTENT(OUT)    :: ierr
  END SUBROUTINE qchle

  FUNCTION dwf(dx, i) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: dx
    INTEGER, INTENT(IN)    :: i
    REAL (dp)              :: fn_val
  END FUNCTION dwf
END INTERFACE

WRITE(*, 1)
1 FORMAT(/)

! epsma is the machine double precision.

epsma = EPSILON(0.0_dp)
iq = 1
idelta = 2
irout = 1
n = 80
mc = 2
mp = 0
ncapm = 81
eps = 5000.*epsma
iem = 0
c = 0.1_dp
xp = 0.0_dp
yp = 0.0_dp
DO  ic=1, 3
  c = 10.*c
  
! Compute the desired recursion coefficients.  On machines with limited
! exponent range, harmless underflow may occur in the routine  gauss
! used in  qchle  to generate the Gauss-Legendre quadrature rule.
  
  CALL dmcdis(n, ncapm, mc, mp, xp, yp, qchle, eps, iq, idelta, irout, finl,  &
              finr, endl, endr, alpha, beta, ncap, kount, ierr,  &
              ie, be, x, w, dwf)
  WRITE(*, 2) ncap, kount, ierr, ie, c
  2 FORMAT('   ncap = ', i3, ' kount = ', i2, ' ierr = ', i3,  &
           ' ie = ', i3, ' for c = ', f5.1)
  IF(ABS(ie) > iem) iem = ABS(ie)
  IF(ie /= 0 .AND. ABS(ie) <= n) THEN
    CALL dmcdis(ABS(ie)-1, ncapm, mc, mp, xp, yp, qchle, eps, iq, idelta,  &
                irout, finl, finr, endl, endr, alpha, beta, ncap, &
                kount, ierr, ie, be, x, w, dwf)
    WRITE(*, 2) ncap, kount, ierr, ie, c
    WRITE(*, 1)
  END IF
  
! Assemble the results in an array.
  
  DO  k=1, n
    IF(ie == 0 .OR. k-1 < ABS(ie)) THEN
      betap(k, ic) = beta(k)
    ELSE
      betap(k, ic) = 0.
    END IF
  END DO
END DO

! Print the results.

WRITE(*, 3)
3 FORMAT(//'   k    beta(k), c=1      beta(k), c=10    beta(k), c=100'/)
DO  k=1, n
  km1 = k-1
  IF(iem == 0 .OR. km1 < iem) WRITE(*, 4) km1, betap(k,1), betap(k,2), betap(k,3)
  4 FORMAT(' ', i3, 3E18.10)
END DO

STOP
END PROGRAM test4



SUBROUTINE qchle(n, x, w, i, ierr)

USE orthpol
USE common_test4
IMPLICIT NONE

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(OUT)  :: x(:)
REAL (dp), INTENT(OUT)  :: w(:)
INTEGER, INTENT(IN)     :: i
INTEGER, INTENT(OUT)    :: ierr

REAL (dp)  :: fk, fn, pi
INTEGER    :: k

fn = n
pi = 4.*ATAN(1.0_dp)
IF(i == 1) THEN
  DO  k=1, n
    fk = k
    x(k) = COS((2.*fk - 1.0_dp)*pi/(2.0_dp*fn))
    w(k) = pi/fn
  END DO
ELSE
  CALL drecur(n, 1, 0._dp, 0._dp, a, b, ierr)
  CALL dgauss(n, a, b, epsma, x, w, ierr)
  w(1:n) = c*w(1:n)
END IF

RETURN
END SUBROUTINE qchle



FUNCTION dwf(x, i) RESULT(fn_val)

! This is a dummy function. It is never called, since the routine
! qgp  in  mcdis  which requires  wf  is not activated in this test.

USE orthpol, ONLY: dp
IMPLICIT NONE
REAL (dp), INTENT(IN)  :: x
INTEGER, INTENT(IN)    :: i
REAL (dp)              :: fn_val

fn_val = 0.

RETURN
END FUNCTION dwf
