MODULE ch12n_func
 
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
 

SUBROUTINE ch12n(n, z, nm, chf1, chd1, chf2, chd2)

!    ====================================================
!    Purpose: Compute Hankel functions of the first and
!             second kinds and their derivatives for a
!             complex argument
!    Input :  z --- Complex argument
!             n --- Order of Hn(1)(z) and Hn(2)(z)
!    Output:  CHF1(n) --- Hn(1)(z)
!             CHD1(n) --- Hn(1)'(z)
!             CHF2(n) --- Hn(2)(z)
!             CHD2(n) --- Hn(2)'(z)
!             NM --- Highest order computed
!    Routines called:
!          (1) CJYNB for computing Jn(z) and Yn(z)
!          (2) CIKNB for computing In(z) and Kn(z)
!    ====================================================

INTEGER, INTENT(IN)        :: n
COMPLEX (dp), INTENT(IN)   :: z
INTEGER, INTENT(OUT)       :: nm
COMPLEX (dp), INTENT(OUT)  :: chf1(0:n)
COMPLEX (dp), INTENT(OUT)  :: chd1(0:n)
COMPLEX (dp), INTENT(OUT)  :: chf2(0:n)
COMPLEX (dp), INTENT(OUT)  :: chd2(0:n)

COMPLEX (dp)  :: cbj(0:250), cdj(0:250), cby(0:250), cdy(0:250),  &
                 cbi(0:250), cdi(0:250), cbk(0:250), cdk(0:250)
COMPLEX (dp)  :: cfac, cf1, ci, zi
INTEGER       :: k
REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp

ci = (0.0_dp, 1.0_dp)
IF (AIMAG(z) < 0.0_dp) THEN
  CALL cjynb(n, z, nm, cbj, cdj, cby, cdy)
  DO  k = 0, nm
    chf1(k) = cbj(k) + ci * cby(k)
    chd1(k) = cdj(k) + ci * cdy(k)
  END DO
  zi = ci * z
  CALL ciknb(n, zi, nm, cbi, cdi, cbk, cdk)
  cfac = -2.0_dp / (pi*ci)
  DO  k = 0, nm
    chf2(k) = cfac * cbk(k)
    chd2(k) = cfac * ci * cdk(k)
    cfac = cfac * ci
  END DO
ELSE IF (AIMAG(z) > 0.0_dp) THEN
  zi = -ci * z
  CALL ciknb(n, zi, nm, cbi, cdi, cbk, cdk)
  cf1 = -ci
  cfac = 2.0_dp / (pi*ci)
  DO  k = 0, nm
    chf1(k) = cfac * cbk(k)
    chd1(k) = -cfac * ci * cdk(k)
    cfac = cfac * cf1
  END DO
  CALL cjynb(n, z, nm, cbj, cdj, cby, cdy)
  DO  k = 0, nm
    chf2(k) = cbj(k) - ci * cby(k)
    chd2(k) = cdj(k) - ci * cdy(k)
  END DO
ELSE
  CALL cjynb(n, z, nm, cbj, cdj, cby, cdy)
  DO  k = 0, nm
    chf1(k) = cbj(k) + ci * cby(k)
    chd1(k) = cdj(k) + ci * cdy(k)
    chf2(k) = cbj(k) - ci * cby(k)
    chd2(k) = cdj(k) - ci * cdy(k)
  END DO
END IF
RETURN
END SUBROUTINE ch12n


SUBROUTINE cjynb(n, z, nm, cbj, cdj, cby, cdy)

!    =======================================================
!    Purpose: Compute Bessel functions Jn(z), Yn(z) and
!             their derivatives for a complex argument
!    Input :  z --- Complex argument of Jn(z) and Yn(z)
!             n --- Order of Jn(z) and Yn(z)
!    Output:  CBJ(n) --- Jn(z)
!             CDJ(n) --- Jn'(z)
!             CBY(n) --- Yn(z)
!             CDY(n) --- Yn'(z)
!             NM --- Highest order computed
!    Routines called:
!             MSTA1 and MSTA2 to calculate the starting
!             point for backward recurrence
!    =======================================================

INTEGER, INTENT(IN)        :: n
COMPLEX (dp), INTENT(IN)   :: z
INTEGER, INTENT(OUT)       :: nm
COMPLEX (dp), INTENT(OUT)  :: cbj(0:n)
COMPLEX (dp), INTENT(OUT)  :: cdj(0:n)
COMPLEX (dp), INTENT(OUT)  :: cby(0:n)
COMPLEX (dp), INTENT(OUT)  :: cdy(0:n)

REAL (dp)     :: a0, y0
COMPLEX (dp)  :: cbj0, cbj1, cbjk, cbs, cby0, cby1, ce, cf, cf1, cf2,  &
                 cp0, cp1, cq0, cq1, cs0, csu, csv, ct1, ct2, cu, cyy
INTEGER       :: k, m
REAL (dp), PARAMETER  :: el = 0.5772156649015329_dp, pi = 3.141592653589793_dp
REAL (dp), PARAMETER  :: a(4) = (/ -.7031250000000000D-01,   &
     0.1121520996093750_dp, -0.5725014209747314_dp, 0.6074042001273483D+01 /)
REAL (dp), PARAMETER  :: b(4) = (/ 0.7324218750000000D-01,   &
    -0.2271080017089844_dp, 0.1727727502584457D+01, -.2438052969955606D+02 /)
REAL (dp), PARAMETER  :: a1(4) = (/ 0.1171875000000000_dp,  &
    -0.1441955566406250_dp, 0.6765925884246826_dp, -0.6883914268109947D+01 /)
REAL (dp), PARAMETER  :: b1(4) = (/ -0.1025390625000000_dp,  &
     0.2775764465332031_dp, -0.1993531733751297D+01, 0.2724882731126854D+02 /)
REAL (dp), PARAMETER  :: r2p = .63661977236758_dp

y0 = ABS(AIMAG(z))
a0 = ABS(z)
nm = n
IF (a0 < 1.0D-100) THEN
  DO  k = 0, n
    cbj(k) = (0.0_dp,0.0_dp)
    cdj(k) = (0.0_dp,0.0_dp)
    cby(k) = -(1.0D+300,0.0_dp)
    cdy(k) = (1.0D+300,0.0_dp)
  END DO
  cbj(0) = (1.0_dp,0.0_dp)
  cdj(1) = (0.5_dp,0.0_dp)
  RETURN
END IF
IF (a0 <= 300._dp .OR. n > 80) THEN
  IF (n == 0) nm = 1
  m = msta1(a0,200)
  IF (m < nm) THEN
    nm = m
  ELSE
    m = msta2(a0, nm, 15)
  END IF
  cbs = (0.0_dp,0.0_dp)
  csu = (0.0_dp,0.0_dp)
  csv = (0.0_dp,0.0_dp)
  cf2 = (0.0_dp,0.0_dp)
  cf1 = (1.0D-100,0.0_dp)
  DO  k = m, 0, -1
    cf = 2.0_dp * (k+1.0_dp) / z * cf1 - cf2
    IF (k <= nm) cbj(k) = cf
    IF (k == 2*INT(k/2) .AND. k /= 0) THEN
      IF (y0 <= 1.0_dp) THEN
        cbs = cbs + 2.0_dp * cf
      ELSE
        cbs = cbs + (-1) ** (k/2) * 2.0_dp * cf
      END IF
      csu = csu + (-1) ** (k/2) * cf / k
    ELSE IF (k > 1) THEN
      csv = csv + (-1) ** (k/2) * k / (k*k-1.0_dp) * cf
    END IF
    cf2 = cf1
    cf1 = cf
  END DO
  IF (y0 <= 1.0_dp) THEN
    cs0 = cbs + cf
  ELSE
    cs0 = (cbs+cf) / COS(z)
  END IF
  DO  k = 0, nm
    cbj(k) = cbj(k) / cs0
  END DO
  ce = LOG(z/2.0_dp) + el
  cby(0) = r2p * (ce*cbj(0) - 4.0_dp*csu/cs0)
  cby(1) = r2p * (-cbj(0)/z + (ce-1.0_dp)*cbj(1) - 4.0_dp*csv/cs0)
ELSE
  ct1 = z - 0.25_dp * pi
  cp0 = (1.0_dp,0.0_dp)
  DO  k = 1, 4
    cp0 = cp0 + a(k) * z ** (-2*k)
  END DO
  cq0 = -0.125_dp / z
  DO  k = 1, 4
    cq0 = cq0 + b(k) * z ** (-2*k-1)
  END DO
  cu = SQRT(r2p/z)
  cbj0 = cu * (cp0*COS(ct1) - cq0*SIN(ct1))
  cby0 = cu * (cp0*SIN(ct1) + cq0*COS(ct1))
  cbj(0) = cbj0
  cby(0) = cby0
  ct2 = z - 0.75_dp * pi
  cp1 = (1.0_dp,0.0_dp)
  DO  k = 1, 4
    cp1 = cp1 + a1(k) * z ** (-2*k)
  END DO
  cq1 = 0.375_dp / z
  DO  k = 1, 4
    cq1 = cq1 + b1(k) * z ** (-2*k-1)
  END DO
  cbj1 = cu * (cp1*COS(ct2) - cq1*SIN(ct2))
  cby1 = cu * (cp1*SIN(ct2) + cq1*COS(ct2))
  cbj(1) = cbj1
  cby(1) = cby1
  DO  k = 2, nm
    cbjk = 2.0_dp * (k-1) / z * cbj1 - cbj0
    cbj(k) = cbjk
    cbj0 = cbj1
    cbj1 = cbjk
  END DO
END IF
cdj(0) = -cbj(1)
DO  k = 1, nm
  cdj(k) = cbj(k-1) - k / z * cbj(k)
END DO
IF (ABS(cbj(0)) > 1.0_dp) THEN
  cby(1) = (cbj(1)*cby(0) - 2.0_dp/(pi*z)) / cbj(0)
END IF
DO  k = 2, nm
  IF (ABS(cbj(k-1)) >= ABS(cbj(k-2))) THEN
    cyy = (cbj(k)*cby(k-1) - 2.0_dp/(pi*z)) / cbj(k-1)
  ELSE
    cyy = (cbj(k)*cby(k-2) - 4.0_dp*(k-1.0_dp)/(pi*z*z)) / cbj(k-2)
  END IF
  cby(k) = cyy
END DO
cdy(0) = -cby(1)
DO  k = 1, nm
  cdy(k) = cby(k-1) - k / z * cby(k)
END DO
RETURN
END SUBROUTINE cjynb



SUBROUTINE ciknb(n, z, nm, cbi, cdi, cbk, cdk)

!    ============================================================
!    Purpose: Compute modified Bessel functions In(z) and Kn(z),
!             and their derivatives for a complex argument
!    Input:   z --- Complex argument
!             n --- Order of In(z) and Kn(z)
!    Output:  CBI(n) --- In(z)
!             CDI(n) --- In'(z)
!             CBK(n) --- Kn(z)
!             CDK(n) --- Kn'(z)
!             NM --- Highest order computed
!    Routones called:
!             MSTA1 and MSTA2 to compute the starting point for
!             backward recurrence
!    ===========================================================

INTEGER, INTENT(IN)        :: n
COMPLEX (dp), INTENT(IN)   :: z
INTEGER, INTENT(OUT)       :: nm
COMPLEX (dp), INTENT(OUT)  :: cbi(0:n)
COMPLEX (dp), INTENT(OUT)  :: cdi(0:n)
COMPLEX (dp), INTENT(OUT)  :: cbk(0:n)
COMPLEX (dp), INTENT(OUT)  :: cdk(0:n)

REAL (dp)     :: a0, fac, vt
COMPLEX (dp)  :: ca0, cbkl, cbs, cf, cf0, cf1, cg, cg0, cg1, ci,  &
                 cr, cs0, csk0, z1
INTEGER       :: k, k0, l, m
REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp, el = 0.57721566490153_dp

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
m = msta1(a0,200)
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
  cbk(0) = -(LOG(0.5_dp*z1) + el) * cbi(0) + cs0 * csk0
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

!    ===================================================
!    Purpose: Determine the starting point for backward
!             recurrence such that the magnitude of
!             Jn(x) at that point is about 10^(-MP)
!    Input :  x     --- Argument of Jn(x)
!             MP    --- Value of magnitude
!    Output:  MSTA1 --- Starting point
!    ===================================================

REAL (dp), INTENT(IN)   :: x
INTEGER, INTENT(IN)     :: mp
INTEGER                 :: fn_val

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

REAL (dp), INTENT(IN)   :: x
INTEGER, INTENT(IN)     :: n
INTEGER, INTENT(IN)     :: mp
INTEGER                 :: fn_val

REAL (dp)  :: a0, ejn, f, f0, f1, hmp, obj
INTEGER    :: it, n0, n1, nn

a0 = ABS(x)
hmp = 0.5_dp * mp
ejn = envj(n,a0)
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
  nn = n1 - (n1-n0) / (1.0_dp-f0/f1)
  f = envj(nn,a0) - obj
  IF (ABS(nn-n1) < 1) EXIT
  n0 = n1
  f0 = f1
  n1 = nn
  f1 = f
END DO

fn_val = nn + 10
RETURN
END FUNCTION msta2



FUNCTION envj(n,x) RESULT(fn_val)

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x
REAL (dp)               :: fn_val

fn_val = 0.5_dp * LOG10(6.28_dp*n) - n * LOG10(1.36_dp*x/n)
RETURN
END FUNCTION envj
 
END MODULE ch12n_func
 
 
 
PROGRAM mch12n
USE ch12n_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:34

!    =====================================================
!    Purpose: This program computes Hankel functions of
!             the first and second kinds and their derivatives
!             for a complex argument using subroutine CH12N
!    Input :  z --- Complex argument
!             n --- Order of Hn(1)(z) and Hn(2)(z)
!                   ( n = 0,1,תתת, n ף 250 )
!    Output:  CHF1(n) --- Hn(1)(z)
!             CHD1(n) --- Hn(1)'(z)
!             CHF2(n) --- Hn(2)(z)
!             CHD2(n) --- Hn(2)'(z)
!    =====================================================

COMPLEX (dp)  :: chf1(0:250), chd1(0:250), chf2(0:250), chd2(0:250), z
REAL (dp)     :: x, y
INTEGER       :: k, n, nm, ns

WRITE (*,*) '  Please enter n, x and y (z=x+iy): '
READ (*,*) n, x, y
WRITE (*,5100) x, y, n
z = CMPLX(x, y, KIND=dp)
IF (n <= 8) THEN
  ns = 1
ELSE
  WRITE (*,*) '  Please enter order step Ns: '
  READ (*,*) ns
END IF
CALL ch12n(n, z, nm, chf1, chd1, chf2, chd2)
WRITE(*,*)
WRITE(*,*) '   n     Re[Hn(1)(z)]     Im[Hn(1)(z)]',  &
           '      Re[Hn(1)''(Z)]     IM[HN(1)''(Z)]'
WRITE(*,*) ' ----------------------------------------------------------------------------'
DO  k = 0, nm, ns
  WRITE (*,5000) k, chf1(k), chd1(k)
END DO
WRITE(*,*)
WRITE(*,*) '   n     Re[Hn(2)(z)]     Im[Hn(2)(z)]',  &
           '      Re[Hn(2)''(Z)]     IM[HN(2)''(Z)]'
WRITE(*,*) ' ----------------------------------------------------------------------------'
DO  k = 0, nm, ns
  WRITE (*,5000) k, chf2(k), chd2(k)
END DO
STOP

5000 FORMAT (' ', i4, 4g18.10)
5100 FORMAT ('   z =', f8.3, ' + i ', f8.3, ' ,      Nmax =', i4)
END PROGRAM mch12n
