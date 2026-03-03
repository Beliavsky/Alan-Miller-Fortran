MODULE lamv_func
 
! From the book "Computation of Special Functions"
!      by Shanjie Zhang and Jianming Jin
!   Copyright 1996 by John Wiley & Sons, Inc.
! The authors state:
!   "However, we give permission to the reader who purchases this book
!    to incorporate any of these programs into his or her programs
!    provided that the copyright is acknowledged."

! Latest revision - 5 February 2002
! Corrections by Alan Miller (amiller @ bigpond.net.au)
! (25) in routine GAM0 changed to g(25)
 
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
 
CONTAINS
 

SUBROUTINE lamv(v, x, vm, vl, dl)

!    =========================================================
!    Purpose: Compute lambda function with arbitrary order v,
!             and their derivative
!    Input :  x --- Argument of lambda function
!             v --- Order of lambda function
!    Output:  VL(n) --- Lambda function of order n+v0
!             DL(n) --- Derivative of lambda function
!             VM --- Highest order computed
!    Routines called:
!         (1) MSTA1 and MSTA2 for computing the starting
!             point for backward recurrence
!         (2) GAM0 for computing gamma function (|x| ó 1)
!    =========================================================

REAL (dp), INTENT(IN)      :: v
REAL (dp), INTENT(IN OUT)  :: x
REAL (dp), INTENT(OUT)     :: vm
REAL (dp), INTENT(OUT)     :: vl(0:)
REAL (dp), INTENT(OUT)     :: dl(0:)

REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp, rp2 = 0.63661977236758_dp
REAL (dp)  :: a0, bjv0, bjv1, bk, ck, cs, f, f0, f1, f2, fac, ga,   &
              px, r, r0, rc, rp, rq, qx, sk, uk, v0, vk, vv, x2, xk
INTEGER    :: i, j, k, k0, m, n

x = ABS(x)
x2 = x * x
n = v
v0 = v - n
vm = v
IF (x <= 12.0_dp) THEN
  DO  k = 0, n
    vk = v0 + k
    bk = 1.0_dp
    r = 1.0_dp
    DO  i = 1, 50
      r = -0.25_dp * r * x2 / (i*(i+vk))
      bk = bk + r
      IF (ABS(r) < ABS(bk)*1.0D-15) EXIT
    END DO
    vl(k) = bk
    uk = 1.0_dp
    r = 1.0_dp
    DO  i = 1, 50
      r = -0.25_dp * r * x2 / (i*(i+vk+1.0_dp))
      uk = uk + r
      IF (ABS(r) < ABS(uk)*1.0D-15) EXIT
    END DO
    dl(k) = -0.5_dp * x / (vk+1.0_dp) * uk
  END DO
  RETURN
END IF
k0 = 11
IF (x >= 35.0_dp) k0 = 10
IF (x >= 50.0_dp) k0 = 8
DO  j = 0, 1
  vv = 4.0_dp * (j+v0) * (j+v0)
  px = 1.0_dp
  rp = 1.0_dp
  DO  k = 1, k0
    rp = -0.78125D-2 * rp * (vv-(4*k-3)**2) * (vv-(4*k-1)**2) / (k*(2*k-1)*x2)
    px = px + rp
  END DO
  qx = 1.0_dp
  rq = 1.0_dp
  DO  k = 1, k0
    rq = -0.78125D-2 * rq * (vv-(4*k-1)**2) * (vv-(4*k+1)**2) / (k*(2*k+1)*x2)
    qx = qx + rq
  END DO
  qx = 0.125_dp * (vv-1.0_dp) * qx / x
  xk = x - (0.5_dp*(j+v0) + 0.25_dp) * pi
  a0 = SQRT(rp2/x)
  ck = COS(xk)
  sk = SIN(xk)
  IF (j == 0) bjv0 = a0 * (px*ck - qx*sk)
  IF (j == 1) bjv1 = a0 * (px*ck - qx*sk)
END DO
IF (v0 == 0.0_dp) THEN
  ga = 1.0_dp
ELSE
  CALL gam0(v0, ga)
  ga = v0 * ga
END IF
fac = (2.0_dp/x) ** v0 * ga
vl(0) = bjv0
dl(0) = -bjv1 + v0 / x * bjv0
vl(1) = bjv1
dl(1) = bjv0 - (1.0_dp+v0) / x * bjv1
r0 = 2.0_dp * (1.0_dp+v0) / x
IF (n <= 1) THEN
  vl(0) = fac * vl(0)
  dl(0) = fac * dl(0) - v0 / x * vl(0)
  vl(1) = fac * r0 * vl(1)
  dl(1) = fac * r0 * dl(1) - (1.0_dp+v0) / x * vl(1)
  RETURN
END IF
IF (n >= 2 .AND. n <= INT(0.9*x)) THEN
  f0 = bjv0
  f1 = bjv1
  DO  k = 2, n
    f = 2.0_dp * (k+v0-1.0_dp) / x * f1 - f0
    f0 = f1
    f1 = f
    vl(k) = f
  END DO
ELSE IF (n >= 2) THEN
  m = msta1(x, 200)
  IF (m < n) THEN
    n = m
  ELSE
    m = msta2(x, n, 15)
  END IF
  f2 = 0.0_dp
  f1 = 1.0D-100
  DO  k = m, 0, -1
    f = 2.0_dp * (v0+k+1.0_dp) / x * f1 - f2
    IF (k <= n) vl(k) = f
    f2 = f1
    f1 = f
  END DO
  IF (ABS(bjv0) > ABS(bjv1)) cs = bjv0 / f
  IF (ABS(bjv0) <= ABS(bjv1)) cs = bjv1 / f2
  vl(0:n) = cs * vl(0:n)
END IF
vl(0) = fac * vl(0)
DO  j = 1, n
  rc = fac * r0
  vl(j) = rc * vl(j)
  dl(j-1) = -0.5_dp * x / (j+v0) * vl(j)
  r0 = 2.0_dp * (j+v0+1) / x * r0
END DO
dl(n) = 2.0_dp * (v0+n) * (vl(n-1)-vl(n)) / x
vm = n + v0
RETURN
END SUBROUTINE lamv



SUBROUTINE gam0(x, ga)

!    ================================================
!    Purpose: Compute gamma function â(x)
!    Input :  x  --- Argument of â(x)  ( |x| ó 1 )
!    Output:  GA --- â(x)
!    ================================================

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: ga

REAL (dp), PARAMETER  :: g(25) = (/  &
       1.0_dp, 0.5772156649015329_dp, -0.6558780715202538_dp,  &
      -0.420026350340952D-1, 0.1665386113822915_dp,   &
      -0.421977345555443D-1, -.96219715278770D-2,  &
      .72189432466630D-2, -.11651675918591D-2, -.2152416741149D-3,  &
      .1280502823882D-3, -.201348547807D-4, -.12504934821D-5,  &
      .11330272320D-5, -.2056338417D-6, .61160950D-8,  &
      .50020075D-8, -.11812746D-8, .1043427D-9, .77823D-11,  &
      -.36968D-11, .51D-12, -.206D-13, -.54D-14, .14D-14 /)
REAL (dp)  :: gr
INTEGER    :: k

gr = g(25)
DO  k = 24, 1, -1
  gr = gr * x + g(k)
END DO
ga = 1.0_dp / (gr*x)
RETURN
END SUBROUTINE gam0


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

END MODULE lamv_func
 
 
 
PROGRAM mlamv
USE lamv_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:43

!    =======================================================
!    Purpose: This program computes the lambda functions
!             for an arbitrary order, and their derivative
!             using subroutine LAMV
!    Input :  x --- Argument of lambda function
!             v --- Order of lambda function
!                   ( v = n+v0, 0 ó n ó 250, 0 ó v0 < 1 )
!    Output:  VL(n) --- Lambda function of order n+v0
!             DL(n) --- Derivative of lambda function
!    Example: x = 10.0

!                v         Lambda(x)        Lambda'(x)
!              ------------------------------------------
!               0.25    -.12510515D+00    -.78558916D-01
!               0.50    -.54402111D-01    -.78466942D-01
!               0.75    -.13657787D-01    -.66234027D-01
!               1.00     .86945492D-02    -.50926063D-01
!               1.25     .19639729D-01    -.36186221D-01
!               1.50     .23540083D-01    -.23382658D-01
!               1.75     .23181910D-01    -.12893894D-01
!               2.00     .20370425D-01    -.46703503D-02
!               2.25     .16283799D-01     .15101684D-02
!               2.50     .11691329D-01     .59243767D-02
!    =======================================================

REAL (dp)  :: vl(0:250), dl(0:250), v, v0, vk, vm, x
INTEGER    :: k, nm, ns

WRITE (*,*) '  Please enter v and x '
READ (*,*) v, x
WRITE (*,5100) v, x
IF (v <= 8) THEN
  ns = 1
ELSE
  WRITE (*,*) '  Please enter order step Ns '
  READ (*,*) ns
END IF
WRITE (*,*)
WRITE (*,*) '   v         Lambda(x)        Lambda''(X)'
WRITE (*,*) '-------------------------------------------'
CALL lamv(v, x, vm, vl, dl)
nm = vm
v0 = vm - nm
DO  k = 0, nm, ns
  vk = k + v0
  WRITE (*,5000) vk, vl(k), dl(k)
END DO
STOP

5000 FORMAT (' ', f6.2, 2g18.8)
5100 FORMAT (' v =', f6.2, '    x =', f8.2)
END PROGRAM mlamv
