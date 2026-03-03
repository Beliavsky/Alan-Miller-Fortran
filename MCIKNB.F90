MODULE ciknb_func
 
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


SUBROUTINE ciknb(n, z, nm, cbi, cdi, cbk, cdk)

!       ============================================================
!       Purpose: Compute modified Bessel functions In(z) and Kn(z),
!                and their derivatives for a complex argument
!       Input:   z --- Complex argument
!                n --- Order of In(z) and Kn(z)
!       Output:  CBI(n) --- In(z)
!                CDI(n) --- In'(z)
!                CBK(n) --- Kn(z)
!                CDK(n) --- Kn'(z)
!                NM --- Highest order computed
!       Routones called:
!                MSTA1 and MSTA2 to compute the starting point for
!                backward recurrence
!       ===========================================================

INTEGER, INTENT(IN)        :: n
COMPLEX (dp), INTENT(IN)   :: z
INTEGER, INTENT(OUT)       :: nm
COMPLEX (dp), INTENT(OUT)  :: cbi(0:n)
COMPLEX (dp), INTENT(OUT)  :: cdi(0:n)
COMPLEX (dp), INTENT(OUT)  :: cbk(0:n)
COMPLEX (dp), INTENT(OUT)  :: cdk(0:n)

REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp, el = 0.57721566490153_dp
REAL (dp)     :: a0, fac, vt
INTEGER       :: k, k0, l, m
COMPLEX (dp)  :: ca0, cbkl, cbs, cf, cf0, cf1, cg, cg0, cg1, ci, cr, cs0, csk0, z1

a0 = ABS(z)
nm = n
IF (a0 < 1.0D-100) THEN
  DO  k = 0, n
    cbi(k) = (0.0_dp,0.0_dp)
    cbk(k) = (1.0D+300,0.0_dp)
    cdi(k) = (0.0_dp,0.0_dp)
    cdk(k) = -(1.0D+300,0.0_dp)
  END DO
  cbi(0) = (1.0_dp,0.0_dp)
  cdi(1) = (0.5_dp,0.0_dp)
  RETURN
END IF
z1 = z
ci = (0.0_dp,1.0_dp)
IF (REAL(z) < 0.0) z1 = -z
IF (n == 0) nm = 1
m = msta1(a0, 200)
IF (m < nm) THEN
  nm = m
ELSE
  m = msta2(a0, nm, 15)
END IF
cbs = 0.0_dp
csk0 = 0.0_dp
cf0 = 0.0_dp
cf1 = 1.0D-100
DO  k = m, 0, -1
  cf = 2.0_dp * (k+1.0_dp) * cf1 / z1 + cf0
  IF (k <= nm) cbi(k) = cf
  IF (k /= 0 .AND. k == 2*INT(k/2)) csk0 = csk0 + 4.0_dp * cf / k
  cbs = cbs + 2.0_dp * cf
  cf0 = cf1
  cf1 = cf
END DO
cs0 = EXP(z1) / (cbs-cf)
DO  k = 0, nm
  cbi(k) = cs0 * cbi(k)
END DO
IF (a0 <= 9.0) THEN
  cbk(0) = -(LOG(0.5_dp*z1)+el) * cbi(0) + cs0 * csk0
  cbk(1) = (1.0_dp/z1-cbi(1)*cbk(0)) / cbi(0)
ELSE
  ca0 = SQRT(pi/(2.0_dp*z1)) * EXP(-z1)
  k0 = 16
  IF (a0 >= 25.0) k0 = 10
  IF (a0 >= 80.0) k0 = 8
  IF (a0 >= 200.0) k0 = 6
  DO  l = 0, 1
    cbkl = 1.0_dp
    vt = 4.0_dp * l
    cr = (1.0_dp,0.0_dp)
    DO  k = 1, k0
      cr = 0.125_dp * cr * (vt-(2.0*k-1.0)**2) / (k*z1)
      cbkl = cbkl + cr
    END DO
    cbk(l) = ca0 * cbkl
  END DO
END IF
cg0 = cbk(0)
cg1 = cbk(1)
DO  k = 2, nm
  cg = 2.0_dp * (k-1.0_dp) / z1 * cg1 + cg0
  cbk(k) = cg
  cg0 = cg1
  cg1 = cg
END DO
IF (REAL(z) < 0.0) THEN
  fac = 1.0_dp
  DO  k = 0, nm
    IF (AIMAG(z) < 0.0) THEN
      cbk(k) = fac * cbk(k) + ci * pi * cbi(k)
    ELSE
      cbk(k) = fac * cbk(k) - ci * pi * cbi(k)
    END IF
    cbi(k) = fac * cbi(k)
    fac = -fac
  END DO
END IF
cdi(0) = cbi(1)
cdk(0) = -cbk(1)
DO  k = 1, nm
  cdi(k) = cbi(k-1) - k / z * cbi(k)
  cdk(k) = -cbk(k-1) - k / z * cbk(k)
END DO
RETURN
END SUBROUTINE ciknb



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
  nn = n1 - (n1-n0) / (1.0D0 - f0/f1)
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
hmp = 0.5D0 * mp
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
  nn = n1 - (n1-n0) / (1.0D0 - f0/f1)
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

fn_val = 0.5D0 * LOG10(6.28D0*n) - n * LOG10(1.36D0*x/n)
RETURN
END FUNCTION envj

END MODULE ciknb_func
 
 
 
PROGRAM mciknb
USE ciknb_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:35

!       =============================================================
!       Purpose: This program computes the modified Bessel functions
!                In(z) and Kn(z), and their derivatives for a
!                complex argument using subroutine CIKNB
!       Input:   z --- Complex argument
!                n --- Order of In(z) and Kn(z)
!                      ( n = 0,1,תתת, n ף 250 )
!       Output:  CBI(n) --- In(z)
!                CDI(n) --- In'(z)
!                CBK(n) --- Kn(z)
!                CDK(n) --- Kn'(z)
!       Example: Nmax = 5,   z = 4.0 + i 2.0

!     n     Re[In(z)]      Im[In(z)]      Re[In'(z)]     Im[In'(z)]
!   -----------------------------------------------------------------
!     0  -.19056142D+01  .10403505D+02 -.23059657D+01  .92222463D+01
!     1  -.23059657D+01  .92222463D+01 -.23666457D+01  .83284588D+01
!     2  -.28276772D+01  .62534130D+01 -.24255774D+01  .61553456D+01
!     3  -.25451891D+01  .30884450D+01 -.22270972D+01  .36367893D+01
!     4  -.16265172D+01  .10201656D+01 -.16520416D+01  .16217056D+01
!     5  -.75889410D+00  .15496632D+00 -.94510625D+00  .48575220D+00

!     n     Re[Kn(z)]      Im[Kn(z)]      Re[Kn'(z)]     Im[Kn'(z)]
!   -----------------------------------------------------------------
!     0  -.64221754D-02 -.84393648D-02  .74307276D-02  .89585853D-02
!     1  -.74307276D-02 -.89585853D-02  .88041795D-02  .94880091D-02
!     2  -.11186184D-01 -.10536653D-01  .14012532D-01  .10936010D-01
!     3  -.20594336D-01 -.12913435D-01  .27416815D-01  .12106413D-01
!     4  -.43647447D-01 -.13676173D-01  .60982763D-01  .63953943D-02
!     5  -.10137119D+00  .12264588D-03  .14495731D+00 -.37132068D-01
!       =============================================================

COMPLEX (dp)  :: cbi(0:250), cdi(0:250), cbk(0:250), cdk(0:250)
REAL (dp)     :: x, y
COMPLEX (dp)  :: z
INTEGER       :: k, n, nm, ns

WRITE (*,*) '  Please enter n, x and y (z = x+iy) '
READ (*,*) n, x, y
WRITE (*,5100) n, x, y
z = CMPLX(x, y, KIND=dp)
WRITE (*,*)
IF (n <= 8) THEN
  ns = 1
ELSE
  WRITE (*,*) '  Please enter order step Ns: '
  READ (*,*) ns
END IF
CALL ciknb(n, z, nm, cbi, cdi, cbk, cdk)
WRITE (*,*)'   n      Re[In(z)]       Im[In(z)]      Re[In''(Z)]       IM[IN''(Z)]'
WRITE (*,*)' ---------------------------------------------------------------------'
DO  k = 0, nm, ns
  WRITE (*,5000) k, cbi(k), cdi(k)
END DO
WRITE(*,*)
WRITE(*,*) '   n      Re[Kn(z)]       Im[Kn(z)]      Re[Kn''(Z)]       IM[KN''(Z)]'
WRITE(*,*) ' ---------------------------------------------------------------------'
DO  k = 0, nm, ns
  WRITE (*,5000) k, cbk(k), cdk(k)
END DO
STOP

5000 FORMAT (t3, i3, ' ', 4g16.8)
5100 FORMAT (t4, 'Nmaz =', i3, ',    z =', f6.1, ' + i', f6.1)
END PROGRAM mciknb
