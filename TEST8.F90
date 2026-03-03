MODULE common_test8
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

! COMMON/s/om2
REAL (dp), SAVE  :: om2

END MODULE common_test8



PROGRAM test8

! Code converted using TO_F90 by Alan Miller
! Date: 2000-02-04  Time: 20:39:50

! This test reproduces the results of  test1  in double precision,
! for n=40, using the routine  mccheb  in place of  cheb  and
! Gauss-Chebyshev quadrature to discretize the modified moments.

USE common_test8, ONLY: om2
USE orthpol
IMPLICIT NONE

! EXTERNAL qcheb
INTERFACE
  SUBROUTINE qcheb(n, x, w, i, ierr)
    USE common_test8
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(OUT)  :: x(:)
    REAL (dp), INTENT(OUT)  :: w(:)
    INTEGER, INTENT(IN)     :: i
    INTEGER, INTENT(OUT)    :: ierr
  END SUBROUTINE qcheb

  FUNCTION dwf(dx, i) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: dx
    INTEGER, INTENT(IN)    :: i
    REAL (dp)              :: fn_val
  END FUNCTION dwf
END INTERFACE

REAL (dp)  :: xp(1), yp(1), endl(1), endr(1), xfer(1), wfer(1), a(79), b(79), &
              fnu(80), alpha(40), beta(40), x(500), w(500), s(40), eps, epsma
LOGICAL    :: finl, finr
REAL (dp), PARAMETER  :: oom2(7) = (/ .1_dp, .3_dp, .5_dp, .7_dp, .9_dp,  &
                                      .99_dp, .999_dp /)
INTEGER    :: idelta, ierr, iom, iq, k, km1, kount, mc, mp, n, ncap, ncapm, &
              ndm1

WRITE(*,*)
epsma = EPSILON(0.0_dp)

! epsma is the machine double precision.

iq = 1
idelta = 1
n = 40
ndm1 = 2*n-1
mc = 1
mp = 0
ncapm = 500
eps = 100.*epsma
xp(1) = 0.0_dp
yp(1) = 0.0_dp

! Generate the recurrence coefficients for the (Chebyshev) polynomials
! defining the modified moments.

CALL drecur(ndm1, 3, 0._dp, 0._dp, a, b, ierr)

! Compute the desired recursion coefficients by the discretized
! Chebyshev algorithm.

DO  iom = 1,7
  om2 = oom2(iom)
  CALL dmcheb(n, ncapm, mc, mp, xp, yp, qcheb, eps, iq, idelta, finl, finr,  &
              endl, endr, xfer, wfer, a, b, fnu, alpha, beta, ncap, kount,  &
              ierr, x, w, s, dwf)
  
! On machines with limited single-precision exponent range, the routine
! cheb  may have generated an underflow exception, which however is
! harmless and can be ignored.
  
  WRITE(*,2) ncap, kount, ierr
  2 FORMAT(/'ncap =', i3, '  kount =', i3, '  ierr =', i3/)
  
! Print the results.
  
  WRITE(*,3)
  3 FORMAT(t6, 'k      beta(k)'/)
  DO  k=1,n
    km1 = k-1
    IF(k == 1) THEN
      WRITE(*,4) km1, beta(k), om2
      4 FORMAT(' ', i5, e18.10, '   om2 =', f6.3)
    ELSE
      WRITE(*,5) km1, beta(k)
      5 FORMAT(' ', i5, e18.10)
    END IF
  END DO
  WRITE(*, *)
END DO

STOP
END PROGRAM test8



SUBROUTINE qcheb(n, x, w, i, ierr)

USE common_test8
IMPLICIT NONE

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(OUT)  :: x(:)
REAL (dp), INTENT(OUT)  :: w(:)
INTEGER, INTENT(IN)     :: i
INTEGER, INTENT(OUT)    :: ierr

! COMMON/s/om2

REAL (dp)  :: pi
INTEGER    :: k

pi = 4.*ATAN(1.0_dp)
DO  k=1,n
  x(k) = COS((2.0_dp*k - 1.0_dp)*pi / (2.0_dp*n))
  w(k) = pi/(n*SQRT(1.0_dp - om2*x(k)**2))
END DO
ierr = 0

RETURN
END SUBROUTINE qcheb



FUNCTION dwf(x, i) RESULT(fn_val)

! This is a dummy function. It is never called, since the routine
! qgp  in  mccheb  which requires  wf  is not activated in this test.

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)  :: x
INTEGER, INTENT(IN)    :: i
REAL (dp)              :: fn_val

fn_val = 0.0_dp

RETURN
END FUNCTION dwf
