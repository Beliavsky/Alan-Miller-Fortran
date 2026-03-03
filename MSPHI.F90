MODULE sphi_func
 
! From the book "Computation of Special Functions"
!      by Shanjie Zhang and Jianming Jin
!   Copyright 1996 by John Wiley & Sons, Inc.
! The authors state:
!   "However, we give permission to the reader who purchases this book
!    to incorporate any of these programs into his or her programs
!    provided that the copyright is acknowledged."
 
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
 
CONTAINS
 

SUBROUTINE sphi(n, x, nm, si, di)

!    ========================================================
!    Purpose: Compute modified spherical Bessel functions
!             of the first kind, in(x) and in'(x)
!    Input :  x --- Argument of in(x)
!             n --- Order of in(x) ( n = 0,1,2,... )
!    Output:  SI(n) --- in(x)
!             DI(n) --- in'(x)
!             NM --- Highest order computed
!    Routines called:
!             MSTA1 and MSTA2 for computing the starting
!             point for backward recurrence
!    ========================================================

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x
INTEGER, INTENT(OUT)    :: nm
REAL (dp), INTENT(OUT)  :: si(0:n)
REAL (dp), INTENT(OUT)  :: di(0:n)

REAL (dp)  :: cs, f, f0, f1, si0
INTEGER    :: k, m

nm = n
IF (ABS(x) < 1.0D-100) THEN
  DO  k = 0, n
    si(k) = 0.0_dp
    di(k) = 0.0_dp
  END DO
  si(0) = 1.0_dp
  di(1) = 0.333333333333333_dp
  RETURN
END IF
si(0) = SINH(x) / x
si(1) = -(SINH(x)/x - COSH(x)) / x
si0 = si(0)
IF (n >= 2) THEN
  m = msta1(x, 200)
  IF (m < n) THEN
    nm = m
  ELSE
    m = msta2(x, n, 15)
  END IF
  f0 = 0.0_dp
  f1 = 1.0_dp - 100
  DO  k = m, 0, -1
    f = (2*k+3) * f1 / x + f0
    IF (k <= nm) si(k) = f
    f0 = f1
    f1 = f
  END DO
  cs = si0 / f
  si(0:nm) = cs * si(0:nm)
END IF
di(0) = si(1)
DO  k = 1, nm
  di(k) = si(k-1) - (k+1)/x * si(k)
END DO
RETURN
END SUBROUTINE sphi


FUNCTION msta1(x, mp) RESULT(fn_val)

!    ===================================================
!    Purpose: Determine the starting point for backward
!             recurrence such that the magnitude of
!             Jn(x) at that point is about 10^(-MP)
!    Input :  x     --- Argument of Jn(x)
!             MP    --- Value of magnitude
!    Output:  MSTA1 --- Starting point
!    ===================================================

REAL (dp), INTENT(IN)  :: x
INTEGER, INTENT(IN)    :: mp
INTEGER                :: fn_val

REAL (dp)  :: a0, f, f0, f1
INTEGER    :: it, n0, n1, nn

a0 = ABS(x)
n0 = INT(1.1*a0) + 1
f0 = envj(n0,a0) - mp
n1 = n0 + 5
f1 = envj(n1,a0) - mp
DO  it = 1, 20
  nn = n1 - (n1-n0) / (1.0_dp - f0/f1)
  f = envj(nn,a0) - mp
  IF (ABS(nn-n1) < 1) EXIT
  n0 = n1
  f0 = f1
  n1 = nn
  f1 = f
END DO

fn_val = nn
RETURN
END FUNCTION msta1



FUNCTION msta2(x, n, mp) RESULT(fn_val)

!    ===================================================
!    Purpose: Determine the starting point for backward
!             recurrence such that all Jn(x) has MP
!             significant digits
!    Input :  x  --- Argument of Jn(x)
!             n  --- Order of Jn(x)
!             MP --- Significant digit
!    Output:  MSTA2 --- Starting point
!    ===================================================

REAL (dp), INTENT(IN)  :: x
INTEGER, INTENT(IN)    :: n
INTEGER, INTENT(IN)    :: mp
INTEGER                :: fn_val

REAL (dp)  :: a0, ejn, f, f0, f1, hmp, obj
INTEGER    :: it, n0, n1, nn

a0 = ABS(x)
hmp = 0.5_dp * mp
ejn = envj(n, a0)
IF (ejn <= hmp) THEN
  obj = mp
  n0 = INT(1.1*a0)
ELSE
  obj = hmp + ejn
  n0 = n
END IF
f0 = envj(n0,a0) - obj
n1 = n0 + 5
f1 = envj(n1,a0) - obj
DO  it = 1, 20
  nn = n1 - (n1-n0) / (1.0_dp - f0/f1)
  f = envj(nn, a0) - obj
  IF (ABS(nn-n1) < 1) EXIT
  n0 = n1
  f0 = f1
  n1 = nn
  f1 = f
END DO

fn_val = nn + 10
RETURN
END FUNCTION msta2



FUNCTION envj(n, x) RESULT(fn_val)

INTEGER, INTENT(IN)    :: n
REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val

fn_val = 0.5_dp * LOG10(6.28_dp*n) - n * LOG10(1.36_dp*x/n)
RETURN
END FUNCTION envj

END MODULE sphi_func
 
 
 
PROGRAM msphi
USE sphi_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:47

!    ======================================================
!    Purpose: This program computes the modified spherical
!             Bessel functions of the first kind in(x) and
!             in'(x) using subroutine SPHI
!    Input :  x --- Argument of in(x)
!             n --- Order of in(x) ( 0 ó n ó 250 )
!    Output:  SI(n) --- in(x)
!             DI(n) --- in'(x)
!    Example: x = 10.0
!               n          in(x)               in'(x)
!             --------------------------------------------
!               0     .1101323287D+04     .9911909633D+03
!               1     .9911909633D+03     .9030850948D+03
!               2     .8039659985D+03     .7500011637D+03
!               3     .5892079640D+03     .5682828129D+03
!               4     .3915204237D+03     .3934477522D+03
!               5     .2368395827D+03     .2494166741D+03
!    ======================================================

REAL (dp)  :: si(0:250), di(0:250), x
INTEGER    :: k, n, nm, ns

WRITE (*,*) 'Please enter n and x '
READ (*,*) n, x
WRITE (*,5100) n, x
IF (n <= 10) THEN
  ns = 1
ELSE
  WRITE (*,*) 'Please enter order step Ns '
  READ (*,*) ns
END IF
CALL sphi(n, x, nm, si, di)
WRITE (*,*)
WRITE (*,*) '  n          in(x)               in''(X)'
WRITE (*,*) '--------------------------------------------'
DO  k = 0, nm, ns
  WRITE (*,5000) k, si(k), di(k)
END DO
STOP

5000 FORMAT (' ', i3, 2g20.10)
5100 FORMAT ('   Nmax =', i3, ',     x =', f6.1)
END PROGRAM msphi
