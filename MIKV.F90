MODULE ikv_func
 
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


SUBROUTINE ikv(v, x, vm, bi, di, bk, dk)

!    =======================================================
!    Purpose: Compute modified Bessel functions Iv(x) and
!             Kv(x), and their derivatives
!    Input :  x --- Argument ( x ע 0 )
!             v --- Order of Iv(x) and Kv(x)
!                   ( v = n+v0, n = 0,1,2,..., 0 ף v0 < 1 )
!    Output:  BI(n) --- In+v0(x)
!             DI(n) --- In+v0'(x)
!             BK(n) --- Kn+v0(x)
!             DK(n) --- Kn+v0'(x)
!             VM --- Highest order computed
!    Routines called:
!         (1) GAMMA for computing the gamma function
!         (2) MSTA1 and MSTA2 to compute the starting
!             point for backward recurrence
!    =======================================================

REAL (dp), INTENT(IN)   :: v
REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: vm
REAL (dp), INTENT(OUT)  :: bi(0:)
REAL (dp), INTENT(OUT)  :: di(0:)
REAL (dp), INTENT(OUT)  :: bk(0:)
REAL (dp), INTENT(OUT)  :: dk(0:)

REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp
REAL (dp)  :: a1, a2, bi0, bk0, bk1, bk2, ca, cb, cs, ct, f, f1, f2, gan, gap,  &
              piv, r, r1, r2, sum, v0, v0n, v0p, vt, w0, wa, ww, x2
INTEGER    :: k, k0, m, n

x2 = x * x
n = INT(v)
v0 = v - n
IF (n == 0) n = 1
IF (x < 1.0D-100) THEN
  DO  k = 0, n
    bi(k) = 0.0D0
    di(k) = 0.0D0
    bk(k) = -1.0D+300
    dk(k) = 1.0D+300
  END DO
  IF (v == 0.0) THEN
    bi(0) = 1.0D0
    di(1) = 0.5D0
  END IF
  vm = v
  RETURN
END IF
piv = pi * v0
vt = 4.0D0 * v0 * v0
IF (v0 == 0.0D0) THEN
  a1 = 1.0D0
ELSE
  v0p = 1.0D0 + v0
  CALL gamma(v0p, gap)
  a1 = (0.5D0*x) ** v0 / gap
END IF
k0 = 14
IF (x >= 35.0) k0 = 10
IF (x >= 50.0) k0 = 8
IF (x <= 18.0) THEN
  bi0 = 1.0D0
  r = 1.0D0
  DO  k = 1, 30
    r = 0.25D0 * r * x2 / (k*(k+v0))
    bi0 = bi0 + r
    IF (ABS(r/bi0) < 1.0D-15) EXIT
  END DO
  bi0 = bi0 * a1
ELSE
  ca = EXP(x) / SQRT(2.0D0*pi*x)
  sum = 1.0D0
  r = 1.0D0
  DO  k = 1, k0
    r = -0.125D0 * r * (vt - (2*k-1)**2) / (k*x)
    sum = sum + r
  END DO
  bi0 = ca * sum
END IF
m = msta1(x,200)
IF (m < n) THEN
  n = m
ELSE
  m = msta2(x, n, 15)
END IF
f2 = 0.0D0
f1 = 1.0D-100
DO  k = m, 0, -1
  f = 2.0D0 * (v0+k+1) / x * f1 + f2
  IF (k <= n) bi(k) = f
  f2 = f1
  f1 = f
END DO
cs = bi0 / f
DO  k = 0, n
  bi(k) = cs * bi(k)
END DO
di(0) = v0 / x * bi(0) + bi(1)
DO  k = 1, n
  di(k) = -(k+v0) / x * bi(k) + bi(k-1)
END DO
IF (x <= 9.0D0) THEN
  IF (v0 == 0.0D0) THEN
    ct = -LOG(0.5D0*x) - 0.5772156649015329_dp
    cs = 0.0D0
    w0 = 0.0D0
    r = 1.0D0
    DO  k = 1, 50
      w0 = w0 + 1.0D0 / k
      r = 0.25D0 * r / (k*k) * x2
      cs = cs + r * (w0+ct)
      wa = ABS(cs)
      IF (ABS((wa-ww)/wa) < 1.0D-15) EXIT
      ww = wa
    END DO
    bk0 = ct + cs
  ELSE
    v0n = 1.0D0 - v0
    CALL gamma(v0n,gan)
    a2 = 1.0D0 / (gan*(0.5D0*x)**v0)
    a1 = (0.5D0*x) ** v0 / gap
    sum = a2 - a1
    r1 = 1.0D0
    r2 = 1.0D0
    DO  k = 1, 120
      r1 = 0.25D0 * r1 * x2 / (k*(k-v0))
      r2 = 0.25D0 * r2 * x2 / (k*(k+v0))
      sum = sum + a2 * r1 - a1 * r2
      wa = ABS(sum)
      IF (ABS((wa-ww)/wa) < 1.0D-15) EXIT
      ww = wa
    END DO
    bk0 = 0.5D0 * pi * sum / SIN(piv)
  END IF
ELSE
  cb = EXP(-x) * SQRT(0.5D0*pi/x)
  sum = 1.0D0
  r = 1.0D0
  DO  k = 1, k0
    r = 0.125D0 * r * (vt-(2.0*k-1.0)**2.0) / (k*x)
    sum = sum + r
  END DO
  bk0 = cb * sum
END IF
bk1 = (1.0D0/x-bi(1)*bk0) / bi(0)
bk(0) = bk0
bk(1) = bk1
DO  k = 2, n
  bk2 = 2.0D0 * (v0+k-1.0D0) / x * bk1 + bk0
  bk(k) = bk2
  bk0 = bk1
  bk1 = bk2
END DO
dk(0) = v0 / x * bk(0) - bk(1)
DO  k = 1, n
  dk(k) = -(k+v0) / x * bk(k) - bk(k-1)
END DO
vm = n + v0
RETURN
END SUBROUTINE ikv


SUBROUTINE gamma(x, ga)

!    ==================================================
!    Purpose: Compute gamma function ג(x)
!    Input :  x  --- Argument of ג(x)
!                    ( x is not equal to 0,-1,-2,תתת)
!    Output:  GA --- ג(x)
!    ==================================================

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: ga

REAL (dp), PARAMETER  :: g(26) = (/ 1.0D0, 0.5772156649015329D0,  &
     -0.6558780715202538D0, -0.420026350340952D-1, 0.1665386113822915D0,   &
     -0.421977345555443D-1, -.96219715278770D-2, .72189432466630D-2,  &
     -0.11651675918591D-2, -.2152416741149D-3, .1280502823882D-3,  &
     -0.201348547807D-4, -.12504934821D-5, .11330272320D-5, -.2056338417D-6, &
      0.61160950D-8, .50020075D-8, -.11812746D-8, .1043427D-9, .77823D-11,  &
      -.36968D-11, .51D-12, -.206D-13, -.54D-14, .14D-14, .1D-15 /)
REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp
REAL (dp)  :: gr, r, z
INTEGER    :: k, m, m1

IF (x == INT(x)) THEN
  IF (x > 0.0D0) THEN
    ga = 1.0D0
    m1 = x - 1
    DO  k = 2, m1
      ga = ga * k
    END DO
  ELSE
    ga = 1.0D+300
  END IF
ELSE
  IF (ABS(x) > 1.0D0) THEN
    z = ABS(x)
    m = z
    r = 1.0_dp
    DO  k = 1, m
      r = r * (z-k)
    END DO
    z = z - m
  ELSE
    z = x
  END IF
  gr = g(26)
  DO  k = 25, 1, -1
    gr = gr * z + g(k)
  END DO
  ga = 1.0D0 / (gr*z)
  IF (ABS(x) > 1.0D0) THEN
    ga = ga * r
    IF (x < 0.0D0) ga = -pi / (x*ga*SIN(pi*x))
  END IF
END IF
RETURN
END SUBROUTINE gamma


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

END MODULE ikv_func
 
 
 
PROGRAM mikv
USE ikv_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:40

!     ============================================================
!     Purpose: This program computes modified Bessel functions
!              Iv(x) and Kv(x) with an arbitrary order, and
!              their derivatives using subroutine IKV
!     Input :  x --- Argument ( x ע 0 )
!              v --- Order of Iv(x) and Kv(x)
!                    ( v = n+v0, 0 ף n ף 250 , 0 ף v0 < 1 )
!     Output:  BI(n) --- In+v0(x)
!              DI(n) --- In+v0'(x)
!              BK(n) --- Kn+v0(x)
!              DK(n) --- Kn+v0'(x)
!     Example: v = n+v0,   v0 = .25,   x = 10.0

!   n       Iv(x)          Iv'(x)         Kv(x)          Kv'(x)
!  ----------------------------------------------------------------
!   0   .28064359D+04  .26631677D+04  .17833184D-04 -.18709581D-04
!   1   .25930068D+04  .24823101D+04  .19155411D-04 -.20227611D-04
!   2   .21581842D+04  .21074153D+04  .22622037D-04 -.24245369D-04
!   3   .16218239D+04  .16310915D+04  .29335327D-04 -.32156018D-04
!   4   .11039987D+04  .11526244D+04  .41690000D-04 -.47053577D-04
!   5   .68342498D+03  .74520058D+03  .64771827D-04 -.75695209D-04
!     =============================================================

REAL (dp)  :: bi(0:250), di(0:250), bk(0:250), dk(0:250), v, v0, vm, x
INTEGER    :: k, n, nm, ns

WRITE (*,*) '  Please enter v, x '
READ (*,*) v, x
n = v
v0 = v - n
WRITE (*,5100) v0, x
IF (n <= 10) THEN
  ns = 1
ELSE
  WRITE (*,*) '  Please enter order step Ns '
  READ (*,*) ns
END IF
CALL ikv(v, x, vm, bi, di, bk, dk)
nm = vm
WRITE(*,*)
WRITE(*,*)'    n       Iv(x)           Iv''(X)          Kv(x)           Kv''(X)'
WRITE(*,*)'  ---------------------------------------------------------------------'
DO  k = 0, nm, ns
  WRITE (*,5000) k, bi(k), di(k), bk(k), dk(k)
END DO
STOP

5000 FORMAT ('  ', i4, ' ', 4g16.8)
5100 FORMAT (t9, 'v = n+v0', ',   v0 =', f7.5, ',   x =', f5.1)
END PROGRAM mikv
