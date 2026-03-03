MODULE iknb_func
 
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
 

SUBROUTINE iknb(n,x,nm,bi,di,bk,dk)

!    ============================================================
!    Purpose: Compute modified Bessel functions In(x) and Kn(x),
!             and their derivatives
!    Input:   x --- Argument of In(x) and Kn(x) ( 0 ó x ó 700 )
!             n --- Order of In(x) and Kn(x)
!    Output:  BI(n) --- In(x)
!             DI(n) --- In'(x)
!             BK(n) --- Kn(x)
!             DK(n) --- Kn'(x)
!             NM --- Highest order computed
!    Routines called:
!             MSTA1 and MSTA2 for computing the starting point
!             for backward recurrence
!    ===========================================================

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x
INTEGER, INTENT(OUT)    :: nm
REAL (dp), INTENT(OUT)  :: bi(0:n)
REAL (dp), INTENT(OUT)  :: di(0:n)
REAL (dp), INTENT(OUT)  :: bk(0:n)
REAL (dp), INTENT(OUT)  :: dk(0:n)

REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp, el = 0.5772156649015329_dp
REAL (dp)  :: a0, bkl, bs, f, f0, f1, g, g0, g1, r, s0, sk0, vt
INTEGER    :: k, k0, l, m

nm = n
IF (x <= 1.0D-100) THEN
  DO  k = 0, n
    bi(k) = 0.0D0
    di(k) = 0.0D0
    bk(k) = 1.0D+300
    dk(k) = -1.0D+300
  END DO
  bi(0) = 1.0D0
  di(1) = 0.5D0
  RETURN
END IF
IF (n == 0) nm = 1
m = msta1(x, 200)
IF (m < nm) THEN
  nm = m
ELSE
  m = msta2(x, nm, 15)
END IF
bs = 0.0D0
sk0 = 0.0D0
f0 = 0.0D0
f1 = 1.0D-100
DO  k = m, 0, -1
  f = 2*(k+1)/x*f1 + f0
  IF (k <= nm) bi(k) = f
  IF (k /= 0 .AND. k == 2*INT(k/2)) sk0 = sk0 + 4.0D0 * f / k
  bs = bs + 2.0D0 * f
  f0 = f1
  f1 = f
END DO
s0 = EXP(x) / (bs-f)
bi(0:nm) = s0 * bi(0:nm)
IF (x <= 8.0D0) THEN
  bk(0) = -(LOG(0.5D0*x)+el) * bi(0) + s0 * sk0
  bk(1) = (1.0D0/x-bi(1)*bk(0)) / bi(0)
ELSE
  a0 = SQRT(pi/(2.0D0*x)) * EXP(-x)
  k0 = 16
  IF (x >= 25.0) k0 = 10
  IF (x >= 80.0) k0 = 8
  IF (x >= 200.0) k0 = 6
  DO  l = 0, 1
    bkl = 1.0D0
    vt = 4 * l
    r = 1.0D0
    DO  k = 1, k0
      r = 0.125D0 * r * (vt - (2*k-1)**2) / (k*x)
      bkl = bkl + r
    END DO
    bk(l) = a0 * bkl
  END DO
END IF
g0 = bk(0)
g1 = bk(1)
DO  k = 2, nm
  g = 2*(k-1)/x*g1 + g0
  bk(k) = g
  g0 = g1
  g1 = g
END DO
di(0) = bi(1)
dk(0) = -bk(1)
DO  k = 1, nm
  di(k) = bi(k-1) - k / x * bi(k)
  dk(k) = -bk(k-1) - k / x * bk(k)
END DO
RETURN
END SUBROUTINE iknb


FUNCTION msta1(x, mp) RESULT(fn_val)

!       ===================================================
!       Purpose: Determine the starting point for backward
!                recurrence such that the magnitude of
!                Jn(x) at that point is about 10^(-MP)
!       Input :  x     --- Argument of Jn(x)
!                MP    --- Value of magnitude
!       Output:  MSTA1 --- Starting point
!       ===================================================

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

!       ===================================================
!       Purpose: Determine the starting point for backward
!                recurrence such that all Jn(x) has MP
!                significant digits
!       Input :  x  --- Argument of Jn(x)
!                n  --- Order of Jn(x)
!                MP --- Significant digit
!       Output:  MSTA2 --- Starting point
!       ===================================================

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

END MODULE iknb_func
 
 
 
PROGRAM miknb
USE iknb_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:40

!       =============================================================
!       Purpose: This program computes modified Bessel functions
!                In(x) and Kn(x), and their derivatives using
!                subroutine IKNB
!       Input:   x --- Argument of In(x) and Kn(x) ( 0 ó x ó 700 )
!                n --- Order of In(x) and Kn(x)
!                      ( n = 0,1,..., n ó 250 )
!       Output:  BI(n) --- In(x)
!                DI(n) --- In'(x)
!                BK(n) --- Kn(x)
!                DK(n) --- Kn'(x)
!       Example: Nmax = 5,    x = 10.0

!     n      In(x)          In'(x)         Kn(x)         Kn'(x)
!    ---------------------------------------------------------------
!     0   .2815717D+04   .2670988D+04   .1778006D-04  -.1864877D-04
!     1   .2670988D+04   .2548618D+04   .1864877D-04  -.1964494D-04
!     2   .2281519D+04   .2214685D+04   .2150982D-04  -.2295074D-04
!     3   .1758381D+04   .1754005D+04   .2725270D-04  -.2968563D-04
!     4   .1226491D+04   .1267785D+04   .3786144D-04  -.4239728D-04
!     5   .7771883D+03   .8378964D+03   .5754185D-04  -.6663236D-04
!       =============================================================

REAL (dp)  :: bi(0:250), di(0:250), bk(0:250), dk(0:250), x
INTEGER    :: k, n, nm, ns

WRITE (*,*) '  Please enter n, x '
READ (*,*) n, x
WRITE (*,5100) n, x
WRITE (*,*)
IF (n <= 10) THEN
  ns = 1
ELSE
  WRITE (*,*) '  Please enter order step Ns '
  READ (*,*) ns
END IF
CALL iknb(n, x, nm, bi, di, bk, dk)
WRITE (*,*) '  n      In(x)          In''(X)         Kn(x)         Kn''(X) '
WRITE (*,*) ' ---------------------------------------------------------------'
DO  k = 0, nm, ns
  WRITE (*,5000) k, bi(k), di(k), bk(k), dk(k)
END DO
STOP

5000 FORMAT (' ', i3, 4g15.7)
5100 FORMAT (t4, 'Nmax =', i3, ',    x =', f5.1)
END PROGRAM miknb
