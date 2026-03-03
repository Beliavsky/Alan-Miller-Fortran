MODULE sphj_func
 
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
 

SUBROUTINE sphj(n, x, nm, sj, dj)

!    =======================================================
!    Purpose: Compute spherical Bessel functions jn(x) and
!             their derivatives
!    Input :  x --- Argument of jn(x)
!             n --- Order of jn(x)  ( n = 0,1,תתת )
!    Output:  SJ(n) --- jn(x)
!             DJ(n) --- jn'(x)
!             NM --- Highest order computed
!    Routines called:
!             MSTA1 and MSTA2 for computing the starting
!             point for backward recurrence
!    =======================================================

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x
INTEGER, INTENT(OUT)    :: nm
REAL (dp), INTENT(OUT)  :: sj(0:n)
REAL (dp), INTENT(OUT)  :: dj(0:n)

REAL (dp)  :: cs, f, f0, f1, sa, sb
INTEGER    :: k, m

nm = n
IF (ABS(x) == 1.0D-100) THEN
  DO  k = 0, n
    sj(k) = 0.0_dp
    dj(k) = 0.0_dp
  END DO
  sj(0) = 1.0_dp
  dj(1) = .3333333333333333_dp
  RETURN
END IF
sj(0) = SIN(x) / x
sj(1) = (sj(0)-COS(x)) / x
IF (n >= 2) THEN
  sa = sj(0)
  sb = sj(1)
  m = msta1(x, 200)
  IF (m < n) THEN
    nm = m
  ELSE
    m = msta2(x, n, 15)
  END IF
  f0 = 0.0_dp
  f1 = 1.0_dp - 100
  DO  k = m, 0, -1
    f = (2*k+3) * f1 / x - f0
    IF (k <= nm) sj(k) = f
    f0 = f1
    f1 = f
  END DO
  IF (ABS(sa) > ABS(sb)) cs = sa / f
  IF (ABS(sa) <= ABS(sb)) cs = sb / f0
  sj(0:nm) = cs * sj(0:nm)
END IF
dj(0) = (COS(x)-SIN(x)/x) / x
DO  k = 1, nm
  dj(k) = sj(k-1) - (k+1) * sj(k) / x
END DO
RETURN
END SUBROUTINE sphj


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

END MODULE sphj_func
 
 
 
PROGRAM msphj
USE sphj_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:47

!    =========================================================
!    Purpose: This program computes the spherical Bessel
!             functions jn(x) and jn'(x) using subroutine SPHJ
!    Input :  x --- Argument of jn(x)
!             n --- Order of jn(x)  ( n = 0,1,תתת,ף 250 )
!    Output:  SJ(n) --- jn(x)
!             DJ(n) --- jn'(x)
!    Example:   x =10.0
!               n          jn(x)               jn(x)
!             --------------------------------------------
!               0    -.5440211109D-01    -.7846694180D-01
!               1     .7846694180D-01    -.7009549945D-01
!               2     .7794219363D-01     .5508428371D-01
!               3    -.3949584498D-01     .9374053162D-01
!               4    -.1055892851D+00     .1329879757D-01
!               5    -.5553451162D-01    -.7226857814D-01
!    =========================================================

REAL (dp)  :: sj(0:250), dj(0:250), x
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
CALL sphj(n, x, nm, sj, dj)
WRITE (*,*)
WRITE (*,*) '  n          jn(x)               jn''(X)'
WRITE (*,*) '--------------------------------------------'
DO  k = 0, nm, ns
  WRITE (*,5000) k, sj(k), dj(k)
END DO
STOP

5000 FORMAT (' ', i3, 2g20.10)
5100 FORMAT ('   Nmax =', i3, ',     x=', f5.1)
END PROGRAM msphj
