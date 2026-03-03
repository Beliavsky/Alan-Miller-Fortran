MODULE jyna_func
 
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


SUBROUTINE jyna(n, x, nm, bj, dj, by, dy)

!       ==========================================================
!       Purpose: Compute Bessel functions Jn(x) & Yn(x) and
!                their derivatives
!       Input :  x --- Argument of Jn(x) & Yn(x)  ( x ע 0 )
!                n --- Order of Jn(x) & Yn(x)
!       Output:  BJ(n) --- Jn(x)
!                DJ(n) --- Jn'(x)
!                BY(n) --- Yn(x)
!                DY(n) --- Yn'(x)
!                NM --- Highest order computed
!       Routines called:
!            (1) JY01B to calculate J0(x), J1(x), Y0(x) & Y1(x)
!            (2) MSTA1 and MSTA2 to calculate the starting
!                point for backward recurrence
!       =========================================================

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x
INTEGER, INTENT(OUT)    :: nm
REAL (dp), INTENT(OUT)  :: bj(0:n)
REAL (dp), INTENT(OUT)  :: dj(0:n)
REAL (dp), INTENT(OUT)  :: by(0:n)
REAL (dp), INTENT(OUT)  :: dy(0:n)

REAL (dp)  :: bj0, bj1, bjk, by0, by1, cs, dj0, dj1, dy0, dy1, f, f0, f1, f2
INTEGER    :: k, m

nm = n
IF (x < 1.0D-100) THEN
  DO  k = 0, n
    bj(k) = 0.0_dp
    dj(k) = 0.0_dp
    by(k) = -1.0D+300
    dy(k) = 1.0D+300
  END DO
  bj(0) = 1.0_dp
  dj(1) = 0.5_dp
  RETURN
END IF
CALL jy01b(x, bj0, dj0, bj1, dj1, by0, dy0, by1, dy1)
bj(0) = bj0
bj(1) = bj1
by(0) = by0
by(1) = by1
dj(0) = dj0
dj(1) = dj1
dy(0) = dy0
dy(1) = dy1
IF (n <= 1) RETURN
IF (n < INT(0.9*x)) THEN
  DO  k = 2, n
    bjk = 2.0_dp * (k-1) / x * bj1 - bj0
    bj(k) = bjk
    bj0 = bj1
    bj1 = bjk
  END DO
ELSE
  m = msta1(x, 200)
  IF (m < n) THEN
    nm = m
  ELSE
    m = msta2(x, n, 15)
  END IF
  f2 = 0.0_dp
  f1 = 1.0D-100
  DO  k = m, 0, -1
    f = 2.0_dp * (k+1) / x * f1 - f2
    IF (k <= nm) bj(k) = f
    f2 = f1
    f1 = f
  END DO
  IF (ABS(bj0) > ABS(bj1)) THEN
    cs = bj0 / f
  ELSE
    cs = bj1 / f2
  END IF
  bj(0:nm) = cs * bj(0:nm)
END IF
DO  k = 2, nm
  dj(k) = bj(k-1) - k / x * bj(k)
END DO
f0 = by(0)
f1 = by(1)
DO  k = 2, nm
  f = 2 * (k-1) / x * f1 - f0
  by(k) = f
  f0 = f1
  f1 = f
END DO
DO  k = 2, nm
  dy(k) = by(k-1) - k * by(k) / x
END DO
RETURN
END SUBROUTINE jyna


SUBROUTINE jy01b(x, bj0, dj0, bj1, dj1, by0, dy0, by1, dy1)

!    =======================================================
!    Purpose: Compute Bessel functions J0(x), J1(x), Y0(x),
!             Y1(x), and their derivatives
!    Input :  x   --- Argument of Jn(x) & Yn(x) ( x ע 0 )
!    Output:  BJ0 --- J0(x)
!             DJ0 --- J0'(x)
!             BJ1 --- J1(x)
!             DJ1 --- J1'(x)
!             BY0 --- Y0(x)
!             DY0 --- Y0'(x)
!             BY1 --- Y1(x)
!             DY1 --- Y1'(x)
!    =======================================================

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: bj0
REAL (dp), INTENT(OUT)  :: dj0
REAL (dp), INTENT(OUT)  :: bj1
REAL (dp), INTENT(OUT)  :: dj1
REAL (dp), INTENT(OUT)  :: by0
REAL (dp), INTENT(OUT)  :: dy0
REAL (dp), INTENT(OUT)  :: by1
REAL (dp), INTENT(OUT)  :: dy1

REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp
REAL (dp)  :: a0, p0, p1, q0, q1, t, t2, ta0, ta1

IF (x == 0.0_dp) THEN
  bj0 = 1.0_dp
  bj1 = 0.0_dp
  dj0 = 0.0_dp
  dj1 = 0.5_dp
  by0 = -1.0D+300
  by1 = -1.0D+300
  dy0 = 1.0D+300
  dy1 = 1.0D+300
  RETURN
ELSE IF (x <= 4.0_dp) THEN
  t = x / 4.0_dp
  t2 = t * t
  bj0 = ((((((-.5014415D-3*t2 + .76771853D-2)*t2 - .0709253492_dp)*t2 +   &
        .4443584263_dp)*t2 - 1.7777560599_dp)*t2 + 3.9999973021_dp)*t2 -  &
        3.9999998721_dp)*t2 + 1.0_dp
  bj1 = t * (((((((-.1289769D-3*t2 + .22069155D-2)*t2 - .0236616773_dp)*t2  &
        + .1777582922_dp)*t2 - .8888839649_dp)*t2 + 2.6666660544_dp)*t2 -  &
        3.9999999710_dp)*t2 + 1.9999999998_dp)
  by0 = (((((((-.567433D-4*t2 + .859977D-3)*t2 - .94855882D-2)*t2 +   &
        .0772975809_dp)*t2 - .4261737419_dp)*t2 + 1.4216421221_dp)*t2 -  &
        2.3498519931_dp)*t2 + 1.0766115157)*t2 + .3674669052_dp
  by0 = 2.0_dp / pi * LOG(x/2.0_dp) * bj0 + by0
  by1 = ((((((((.6535773D-3*t2 - .0108175626_dp)*t2 + .107657606_dp)*t2 -  &
        .7268945577_dp)*t2 + 3.1261399273_dp)*t2 - 7.3980241381_dp)*t2 +   &
        6.8529236342_dp)*t2 + .3932562018_dp)*t2 - .6366197726_dp) / x
  by1 = 2.0_dp / pi * LOG(x/2.0_dp) * bj1 + by1
ELSE
  t = 4.0_dp / x
  t2 = t * t
  a0 = SQRT(2.0_dp/(pi*x))
  p0 = ((((-.9285D-5*t2 + .43506D-4)*t2 - .122226D-3)*t2 + .434725D-3)*t2  &
       - .4394275D-2)*t2 + .999999997_dp
  q0 = t * (((((.8099D-5*t2 - .35614D-4)*t2 + .85844D-4)*t2 - .218024D-3)*t2  &
       + .1144106D-2)*t2 - .031249995_dp)
  ta0 = x - .25_dp * pi
  bj0 = a0 * (p0*COS(ta0) - q0*SIN(ta0))
  by0 = a0 * (p0*SIN(ta0) + q0*COS(ta0))
  p1 = ((((.10632D-4*t2 - .50363D-4)*t2 + .145575D-3)*t2 - .559487D-3)*t2  &
       + .7323931D-2)*t2 + 1.000000004_dp
  q1 = t * (((((-.9173D-5*t2 + .40658D-4)*t2 - .99941D-4)*t2 +   &
      .266891D-3)*t2 - .1601836D-2)*t2 + .093749994_dp)
  ta1 = x - .75_dp * pi
  bj1 = a0 * (p1*COS(ta1) - q1*SIN(ta1))
  by1 = a0 * (p1*SIN(ta1) + q1*COS(ta1))
END IF
dj0 = -bj1
dj1 = bj0 - bj1 / x
dy0 = -by1
dy1 = by0 - by1 / x
RETURN
END SUBROUTINE jy01b


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

END MODULE jyna_func
 
 
 
PROGRAM mjyna
USE jyna_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:42

!    ====================================================
!    Purpose: This program computes Bessel functions
!             Jn(x) and Yn(x), and their derivatives
!             using subroutine JYNA
!    Input :  x --- Argument of Jn(x) & Yn(x)  ( x ע 0 )
!             n --- Order of Jn(x) & Yn(x)
!                   ( n = 0,1,2,תתת, n ף 250 )
!    Output:  BJ(n) --- Jn(x)
!             DJ(n) --- Jn'(x)
!             BY(n) --- Yn(x)
!             DY(n) --- Yn'(x)
!    Example:
!             x = 10.0

!             n        Jn(x)           Jn'(x)
!           -------------------------------------
!             0    -.2459358D+00   -.4347275D-01
!            10     .2074861D+00    .8436958D-01
!            20     .1151337D-04    .2011954D-04
!            30     .1551096D-11    .4396479D-11

!             n        Yn(x)           Yn'(x)
!           -------------------------------------
!             0     .5567117D-01   -.2490154D+00
!            10    -.3598142D+00    .1605149D+00
!            20    -.1597484D+04    .2737803D+04
!            30    -.7256142D+10    .2047617D+11
!    ====================================================

REAL (dp)  :: bj(0:250), by(0:250), dj(0:250), dy(0:250), x
INTEGER    :: k, n, nm, ns

WRITE (*,*) '  Please enter n, x '
READ (*,*) n, x
WRITE (*,5000) n, x
IF (n <= 8) THEN
  ns = 1
ELSE
  WRITE (*,*) '  Please enter order step Ns '
  READ (*,*) ns
END IF
WRITE (*,*)
CALL jyna(n, x, nm, bj, dj, by, dy)
WRITE (*,*) '  n        Jn(x)           Jn''(X)'
WRITE (*,*) '--------------------------------------'
DO  k = 0, nm, ns
  WRITE (*,5100) k, bj(k), dj(k)
END DO
WRITE (*,*)
WRITE (*,*) '  n        Yn(x)           Yn''(X)'
WRITE (*,*) '--------------------------------------'
DO  k = 0, nm, ns
  WRITE (*,5100) k, by(k), dy(k)
END DO
STOP

5000 FORMAT ('   Nmax =', i3, ',    x =', f6.1)
5100 FORMAT (' ', i3, ' ', 2g16.7)
END PROGRAM mjyna
