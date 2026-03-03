MODULE ikna_func
 
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
 

SUBROUTINE ikna(n, x, nm, bi, di, bk, dk)

!    ========================================================
!    Purpose: Compute modified Bessel functions In(x) and
!             Kn(x), and their derivatives
!    Input:   x --- Argument of In(x) and Kn(x) ( x ע 0 )
!             n --- Order of In(x) and Kn(x)
!    Output:  BI(n) --- In(x)
!             DI(n) --- In'(x)
!             BK(n) --- Kn(x)
!             DK(n) --- Kn'(x)
!             NM --- Highest order computed
!    Routines called:
!         (1) IK01A for computing I0(x),I1(x),K0(x) & K1(x)
!         (2) MSTA1 and MSTA2 for computing the starting
!             point for backward recurrence
!    ========================================================

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x
INTEGER, INTENT(OUT)    :: nm
REAL (dp), INTENT(OUT)  :: bi(0:n)
REAL (dp), INTENT(OUT)  :: di(0:n)
REAL (dp), INTENT(OUT)  :: bk(0:n)
REAL (dp), INTENT(OUT)  :: dk(0:n)

INTEGER    :: k, m
REAL (dp)  :: bi0, bi1, bk0, bk1, di0, di1, dk0, dk1, f, f0, f1, g, g0, g1,  &
              h, h0, h1, s0

nm = n
IF (x <= 1.0D-100) THEN
  DO  k = 0, n
    bi(k) = 0.0_dp
    di(k) = 0.0_dp
    bk(k) = 1.0D+300
    dk(k) = -1.0D+300
  END DO
  bi(0) = 1.0_dp
  di(1) = 0.5_dp
  RETURN
END IF
CALL ik01a(x, bi0, di0, bi1, di1, bk0, dk0, bk1, dk1)
bi(0) = bi0
bi(1) = bi1
bk(0) = bk0
bk(1) = bk1
di(0) = di0
di(1) = di1
dk(0) = dk0
dk(1) = dk1
IF (n <= 1) RETURN
IF (x > 40.0 .AND. n < 0.25*x) THEN
  h0 = bi0
  h1 = bi1
  DO  k = 2, n
    h = -2*(k-1)/x*h1 + h0
    bi(k) = h
    h0 = h1
    h1 = h
  END DO
ELSE
  m = msta1(x, 200)
  IF (m < n) THEN
    nm = m
  ELSE
    m = msta2(x, n, 15)
  END IF
  f0 = 0.0_dp
  f1 = 1.0D-100
  DO  k = m, 0, -1
    f = 2*(k+1)*f1/x + f0
    IF (k <= nm) bi(k) = f
    f0 = f1
    f1 = f
  END DO
  s0 = bi0 / f
  bi(0:nm) = s0 * bi(0:nm)
END IF
g0 = bk0
g1 = bk1
DO  k = 2, nm
  g = 2*(k-1)/x*g1 + g0
  bk(k) = g
  g0 = g1
  g1 = g
END DO
DO  k = 2, nm
  di(k) = bi(k-1) - k / x * bi(k)
  dk(k) = -bk(k-1) - k / x * bk(k)
END DO
RETURN
END SUBROUTINE ikna


SUBROUTINE ik01a(x, bi0, di0, bi1, di1, bk0, dk0, bk1, dk1)

!    =========================================================
!    Purpose: Compute modified Bessel functions I0(x), I1(1),
!             K0(x) and K1(x), and their derivatives
!    Input :  x   --- Argument ( x ע 0 )
!    Output:  BI0 --- I0(x)
!             DI0 --- I0'(x)
!             BI1 --- I1(x)
!             DI1 --- I1'(x)
!             BK0 --- K0(x)
!             DK0 --- K0'(x)
!             BK1 --- K1(x)
!             DK1 --- K1'(x)
!    =========================================================

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: bi0
REAL (dp), INTENT(OUT)  :: di0
REAL (dp), INTENT(OUT)  :: bi1
REAL (dp), INTENT(OUT)  :: di1
REAL (dp), INTENT(OUT)  :: bk0
REAL (dp), INTENT(OUT)  :: dk0
REAL (dp), INTENT(OUT)  :: bk1
REAL (dp), INTENT(OUT)  :: dk1

REAL (dp), PARAMETER  :: a(12) = (/ 0.125_dp, 7.03125D-2, 7.32421875D-2,  &
      1.1215209960938D-1, 2.2710800170898D-1, 5.7250142097473D-1,  &
      1.7277275025845_dp, 6.0740420012735_dp, 2.4380529699556D01,  &
      1.1001714026925D02, 5.5133589612202D02, 3.0380905109224D03 /)
REAL (dp), PARAMETER  :: b(12) = (/ -0.375_dp, -1.171875D-1, -1.025390625D-1, &
      -1.4419555664063D-1, -2.7757644653320D-1, -6.7659258842468D-1,  &
      -1.9935317337513_dp, -6.8839142681099_dp, -2.7248827311269D01,  &
      -1.2159789187654D02, -6.0384407670507D02, -3.3022722944809D03 /)
REAL (dp), PARAMETER  :: a1(8) = (/ 0.125_dp, 0.2109375_dp, 1.0986328125_dp,  &
      1.1775970458984D01, 2.1461706161499D02, 5.9511522710323D03,  &
      2.3347645606175D05, 1.2312234987631D07 /)
REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp, el = 0.5772156649015329_dp
REAL (dp)  :: ca, cb, ct, r, w0, ww, x2, xr, xr2
INTEGER    :: k, k0

x2 = x * x
IF (x == 0.0_dp) THEN
  bi0 = 1.0_dp
  bi1 = 0.0_dp
  bk0 = 1.0D+300
  bk1 = 1.0D+300
  di0 = 0.0_dp
  di1 = 0.5_dp
  dk0 = -1.0D+300
  dk1 = -1.0D+300
  RETURN
ELSE IF (x <= 18.0) THEN
  bi0 = 1.0_dp
  r = 1.0_dp
  DO  k = 1, 50
    r = 0.25_dp * r * x2 / (k*k)
    bi0 = bi0 + r
    IF (ABS(r/bi0) < 1.0D-15) EXIT
  END DO
  bi1 = 1.0_dp
  r = 1.0_dp
  DO  k = 1, 50
    r = 0.25_dp * r * x2 / (k*(k+1))
    bi1 = bi1 + r
    IF (ABS(r/bi1) < 1.0D-15) EXIT
  END DO
  bi1 = 0.5_dp * x * bi1
ELSE
  k0 = 12
  IF (x >= 35.0) k0 = 9
  IF (x >= 50.0) k0 = 7
  ca = EXP(x) / SQRT(2.0_dp*pi*x)
  bi0 = 1.0_dp
  xr = 1.0_dp / x
  DO  k = 1, k0
    bi0 = bi0 + a(k) * xr ** k
  END DO
  bi0 = ca * bi0
  bi1 = 1.0_dp
  DO  k = 1, k0
    bi1 = bi1 + b(k) * xr ** k
  END DO
  bi1 = ca * bi1
END IF
IF (x <= 9.0_dp) THEN
  ct = -(LOG(x/2.0_dp)+el)
  bk0 = 0.0_dp
  w0 = 0.0_dp
  r = 1.0_dp
  DO  k = 1, 50
    w0 = w0 + 1.0_dp / k
    r = 0.25_dp * r / (k*k) * x2
    bk0 = bk0 + r * (w0+ct)
    IF (ABS((bk0-ww)/bk0) < 1.0D-15) EXIT
    ww = bk0
  END DO
  bk0 = bk0 + ct
ELSE
  cb = 0.5_dp / x
  xr2 = 1.0_dp / x2
  bk0 = 1.0_dp
  DO  k = 1, 8
    bk0 = bk0 + a1(k) * xr2 ** k
  END DO
  bk0 = cb * bk0 / bi0
END IF
bk1 = (1.0_dp/x - bi1*bk0) / bi0
di0 = bi1
di1 = bi0 - bi1 / x
dk0 = -bk1
dk1 = -bk0 - bk1 / x
RETURN
END SUBROUTINE ik01a


FUNCTION msta1(x, mp) RESULT(fn_val)

!       ===================================================
!       Purpose: Determine the starting point for backward
!                recurrence such that the magnitude of
!                Jn(x) at that point is about 10^(-MP)
!       Input :  x     --- Argument of Jn(x)
!                MP    --- Value of magnitude
!       Output:  MSTA1 --- Starting point
!       ===================================================

REAL (dp), INTENT(IN)      :: x
INTEGER, INTENT(IN)        :: mp
INTEGER                    :: fn_val

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

REAL (dp), INTENT(IN)      :: x
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: mp
INTEGER                    :: fn_val

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

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: x
REAL (dp)                  :: fn_val

fn_val = 0.5_dp * LOG10(6.28_dp*n) - n * LOG10(1.36_dp*x/n)
RETURN
END FUNCTION envj

END MODULE ikna_func
 
 
 
PROGRAM mikna
USE ikna_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:40

!     ===============================================================
!     Purpose: This program computes modified Bessel functions In(x)
!              and Kn(x), and their derivatives using subroutine IKNA
!     Input:   x --- Argument of In(x) and Kn(x) ( x ע 0 )
!              n --- Order of In(x) and Kn(x)
!                    ( n = 0,1,תתת, n ף 250 )
!     Output:  BI(n) --- In(x)
!              DI(n) --- In'(x)
!              BK(n) --- Kn(x)
!              DK(n) --- Kn'(x)
!     Example: Nmax = 5,    x = 10.0

!   n      In(x)          In'(x)         Kn(x)         Kn'(x)
!  ---------------------------------------------------------------
!   0   .2815717D+04   .2670988D+04   .1778006D-04  -.1864877D-04
!   1   .2670988D+04   .2548618D+04   .1864877D-04  -.1964494D-04
!   2   .2281519D+04   .2214685D+04   .2150982D-04  -.2295074D-04
!   3   .1758381D+04   .1754005D+04   .2725270D-04  -.2968563D-04
!   4   .1226491D+04   .1267785D+04   .3786144D-04  -.4239728D-04
!   5   .7771883D+03   .8378964D+03   .5754185D-04  -.6663236D-04
!     ===============================================================

REAL (dp)  :: bi(0:250), di(0:250), bk(0:250), dk(0:250), x
INTEGER    :: k, n, nm, ns

WRITE (*,*) '  Please enter n, x '
READ (*,*) n, x
WRITE (*,5100) n, x
WRITE (*,*)
IF (n <= 10) THEN
  ns = 1
ELSE
  WRITE (*,*) '  Please enter order step Ns: '
  READ (*,*) ns
END IF
CALL ikna(n, x, nm, bi, di, bk, dk)
WRITE (*,*) '  n      In(x)          In''(X)         Kn(x)         Kn''(X) '
WRITE (*,*) ' ---------------------------------------------------------------'
DO  k = 0, nm, ns
  WRITE (*,5000) k, bi(k), di(k), bk(k), dk(k)
END DO
STOP

5000 FORMAT (' ', i3, 4g15.7)
5100 FORMAT (t4, 'Nmax =', i3, ',    x =', f5.1)
END PROGRAM mikna
