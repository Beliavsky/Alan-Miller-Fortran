MODULE rswfo_func
 
! From the book "Computation of Special Functions"
!      by Shanjie Zhang and Jianming Jin
!   Copyright 1996 by John Wiley & Sons, Inc.
! The authors state:
!   "However, we give permission to the reader who purchases this book
!    to incorporate any of these programs into his or her programs
!    provided that the copyright is acknowledged."

! Latest revision - 30 December 2001
! Corrections by Alan Miller (amiller @ bigpond.net.au)
! Variable sw was used without a value assigned to it.
! There was a call to a routine with the same argument variable repeated.
! This is only permitted when the called routine does not alter either
! argument - i.e. when they are both INTENT(IN) arguments.  The duplicated
! argument was nm2 and the second one has been replaced with nm3.

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
 
CONTAINS
 

SUBROUTINE rswfo(m, n, c, x, cv, kf, r1f, r1d, r2f, r2d)

!    ==========================================================
!    Purpose: Compute oblate radial functions of the first
!             and second kinds, and their derivatives
!    Input :  m  --- Mode parameter,  m = 0,1,2,...
!             n  --- Mode parameter,  n = m,m+1,m+2,...
!             c  --- Spheroidal parameter
!             x  --- Argument (x ע 0)
!             cv --- Characteristic value
!             KF --- Function code
!                    KF=1 for the first kind
!                    KF=2 for the second kind
!                    KF=3 for both the first and second kinds
!    Output:  R1F --- Radial function of the first kind
!             R1D --- Derivative of the radial function of
!                     the first kind
!             R2F --- Radial function of the second kind
!             R2D --- Derivative of the radial function of
!                     the second kind
!    Routines called:
!         (1) SDMN for computing expansion coefficients dk
!         (2) RMN1 for computing prolate or oblate radial
!             function of the first kind
!         (3) RMN2L for computing prolate or oblate radial
!             function of the second kind for a large argument
!         (4) RMN2SO for computing oblate radial functions of
!             the second kind for a small argument
!    ==========================================================

INTEGER, INTENT(IN OUT)    :: m
INTEGER, INTENT(IN OUT)    :: n
REAL (dp), INTENT(IN OUT)  :: c
REAL (dp), INTENT(IN)      :: x
REAL (dp), INTENT(IN OUT)  :: cv
INTEGER, INTENT(IN)        :: kf
REAL (dp), INTENT(IN OUT)  :: r1f
REAL (dp), INTENT(IN OUT)  :: r1d
REAL (dp), INTENT(IN OUT)  :: r2f
REAL (dp), INTENT(IN OUT)  :: r2d

REAL (dp)  :: df(200)
INTEGER    :: id, kd

kd = -1
CALL sdmn(m, n, c, cv, kd, df)
IF (kf /= 2) THEN
  CALL rmn1(m, n, c, x, df, kd, r1f, r1d)
END IF
IF (kf > 1) THEN
  id = 10
  IF (x > 1.0D-8) THEN
    CALL rmn2l(m, n, c, x, df, kd, r2f, r2d, id)
  END IF
  IF (id > -1) THEN
    CALL rmn2so(m, n, c, x, cv, df, kd, r2f, r2d)
  END IF
END IF
RETURN
END SUBROUTINE rswfo



SUBROUTINE sdmn(m, n, c, cv, kd, df)

!    =====================================================
!    Purpose: Compute the expansion coefficients of the
!             prolate and oblate spheroidal functions, dk
!    Input :  m  --- Mode parameter
!             n  --- Mode parameter
!             c  --- Spheroidal parameter
!             cv --- Characteristic value
!             KD --- Function code
!                    KD=1 for prolate; KD=-1 for oblate
!    Output:  DF(k) --- Expansion coefficients dk;
!                       DF(1), DF(2), ... correspond to
!                       d0, d2, ... for even n-m and d1,
!                       d3, ... for odd n-m
!    =====================================================

INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: c
REAL (dp), INTENT(IN OUT)  :: cv
INTEGER, INTENT(IN)        :: kd
REAL (dp), INTENT(OUT)     :: df(200)

REAL (dp)  :: a(200), d(200), g(200)
REAL (dp)  :: cs, d2k, dk0, dk1, dk2, f, f0, f1, f2, fl, fs, r1, r3, r4,  &
              s0, su1, su2, sw
INTEGER    :: i, ip, j, k, k1, kb, nm

nm = 25 + INT(0.5*(n-m)+c)
IF (c < 1.0D-10) THEN
  DO  i = 1, nm
    df(i) = 0_dp
  END DO
  df((n-m)/2+1) = 1.0_dp
  RETURN
END IF
cs = c * c * kd
ip = 1
IF (n-m == 2*INT((n-m)/2)) ip = 0
DO  i = 1, nm + 2
  IF (ip == 0) k = 2 * (i-1)
  IF (ip == 1) k = 2 * i - 1
  dk0 = m + k
  dk1 = m + k + 1
  dk2 = 2 * (m+k)
  d2k = 2 * m + k
  a(i) = (d2k+2.0) * (d2k+1.0) / ((dk2+3.0)*(dk2+5.0)) * cs
  d(i) = dk0 * dk1 + (2.0*dk0*dk1 - 2.0*m*m - 1.0) /   &
         ((dk2 - 1.0)*(dk2 + 3.0)) * cs
  g(i) = k * (k-1.0) / ((dk2-3.0)*(dk2-1.0)) * cs
END DO
fs = 1.0_dp
f1 = 0.0_dp
f0 = 1.0D-100
kb = 0
df(nm+1) = 0.0_dp
DO  k = nm, 1, -1
  f = -((d(k+1)-cv)*f0 + a(k+1)*f1) / g(k+1)
  IF (ABS(f) > ABS(df(k+1))) THEN
    df(k) = f
    fl = df(k+1)
    f1 = f0
    f0 = f
    IF (ABS(f) > 1.0D+100) THEN
      DO  k1 = k, nm
        df(k1) = df(k1) * 1.0D-100
      END DO
      f1 = f1 * 1.0D-100
      f0 = f0 * 1.0D-100
    END IF
  ELSE
    kb = k
    fl = df(k+1)
    f1 = 1.0D-100
    f2 = -(d(1)-cv) / a(1) * f1
    df(1) = f1
    IF (kb == 1) THEN
      fs = f2
    ELSE IF (kb == 2) THEN
      df(2) = f2
      fs = -((d(2)-cv)*f2+g(2)*f1) / a(2)
    ELSE
      df(2) = f2
      DO  j = 3, kb + 1
        f = -((d(j-1)-cv)*f2+g(j-1)*f1) / a(j-1)
        IF (j <= kb) df(j) = f
        IF (ABS(f) > 1.0D+100) THEN
          DO  k1 = 1, j
            df(k1) = df(k1) * 1.0D-100
          END DO
          f = f * 1.0D-100
          f2 = f2 * 1.0D-100
        END IF
        f1 = f2
        f2 = f
      END DO
      fs = f
    END IF
    EXIT
  END IF
END DO

su1 = 0.0_dp
r1 = 1.0_dp
DO  j = m + ip + 1, 2 * (m+ip)
  r1 = r1 * j
END DO
su1 = df(1) * r1
DO  k = 2, kb
  r1 = -r1 * (k+m+ip-1.5_dp) / (k-1.0_dp)
  su1 = su1 + r1 * df(k)
END DO
su2 = 0.0_dp
sw = su2
DO  k = kb + 1, nm
  IF (k /= 1) r1 = -r1 * (k+m+ip-1.5_dp) / (k-1)
  su2 = su2 + r1 * df(k)
  IF (ABS(sw-su2) < ABS(su2)*1.0D-14) EXIT
  sw = su2
END DO

r3 = 1.0_dp
DO  j = 1, (m+n+ip) / 2
  r3 = r3 * (j + 0.5_dp*(n+m+ip))
END DO
r4 = 1.0_dp
DO  j = 1, (n-m-ip) / 2
  r4 = -4.0_dp * r4 * j
END DO
s0 = r3 / (fl*(su1/fs) + su2) / r4
DO  k = 1, kb
  df(k) = fl / fs * s0 * df(k)
END DO
DO  k = kb + 1, nm
  df(k) = s0 * df(k)
END DO
RETURN
END SUBROUTINE sdmn



SUBROUTINE rmn1(m, n, c, x, df, kd, r1f, r1d)

!    =======================================================
!    Purpose: Compute prolate and oblate spheroidal radial
!             functions of the first kind for given m, n,
!             c and x
!    Routines called:
!         (1) SCKB for computing expansion coefficients c2k
!         (2) SPHJ for computing the spherical Bessel
!             functions of the first kind
!    =======================================================

INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: c
REAL (dp), INTENT(IN)      :: x
REAL (dp), INTENT(IN)      :: df(200)
INTEGER, INTENT(IN)        :: kd
REAL (dp), INTENT(OUT)     :: r1f
REAL (dp), INTENT(OUT)     :: r1d

REAL (dp)  :: ck(200), sj(0:251), dj(0:251)
REAL (dp)  :: a0, b0, cx, eps, r, r0, r1, r2, r3, reg,  &
              sa0, suc, sud, sum, sw, sw1
INTEGER    :: ip, j, k, l, lg, nm, nm1, nm2, np

eps = 1.0D-14
ip = 1
nm1 = INT((n-m)/2)
IF (n-m == 2*nm1) ip = 0
nm = 25 + nm1 + INT(c)
reg = 1.0_dp
IF (m+nm > 80) reg = 1.0D-200
r0 = reg
DO  j = 1, 2 * m + ip
  r0 = r0 * j
END DO
r = r0
suc = r * df(1)
DO  k = 2, nm
  r = r * (m+k-1.0) * (m+k+ip-1.5_dp) / (k-1.0_dp) / (k+ip-1.5_dp)
  suc = suc + r * df(k)
  IF (k > nm1 .AND. ABS(suc-sw) < ABS(suc)*eps) EXIT
  sw = suc
END DO

IF (x == 0.0) THEN
  CALL sckb(m, n, c, df, ck)
  sum = 0.0_dp
  DO  j = 1, nm
    sum = sum + ck(j)
    IF (ABS(sum-sw1) < ABS(sum)*eps) EXIT
    sw1 = sum
  END DO

  r1 = 1.0_dp
  DO  j = 1, (n+m+ip) / 2
    r1 = r1 * (j+0.5_dp*(n+m+ip))
  END DO
  r2 = 1.0_dp
  DO  j = 1, m
    r2 = 2.0_dp * c * r2 * j
  END DO
  r3 = 1.0_dp
  DO  j = 1, (n-m-ip) / 2
    r3 = r3 * j
  END DO
  sa0 = (2*(m+ip)+1) * r1 / (2.0**n*c**ip*r2*r3)
  IF (ip == 0) THEN
    r1f = sum / (sa0*suc) * df(1) * reg
    r1d = 0.0_dp
  ELSE IF (ip == 1) THEN
    r1f = 0.0_dp
    r1d = sum / (sa0*suc) * df(1) * reg
  END IF
  RETURN
END IF
cx = c * x
nm2 = 2 * nm + m
CALL sphj(nm2, cx, nm2, sj, dj)
a0 = (1.0_dp-kd/(x*x)) ** (0.5_dp*m) / suc
r1f = 0.0_dp
DO  k = 1, nm
  l = 2 * k + m - n - 2 + ip
  IF (l == 4*INT(l/4)) lg = 1
  IF (l /= 4*INT(l/4)) lg = -1
  IF (k == 1) THEN
    r = r0
  ELSE
    r = r * (m+k-1.0) * (m+k+ip-1.5_dp) / (k-1.0_dp) / (k+ip-1.5_dp)
  END IF
  np = m + 2 * k - 2 + ip
  r1f = r1f + lg * r * df(k) * sj(np)
  IF (k > nm1 .AND. ABS(r1f-sw) < ABS(r1f)*eps) EXIT
  sw = r1f
END DO

r1f = r1f * a0
b0 = kd * m / x ** 3 / (1.0 - kd/(x*x)) * r1f
sud = 0.0_dp
DO  k = 1, nm
  l = 2 * k + m - n - 2 + ip
  IF (l == 4*INT(l/4)) lg = 1
  IF (l /= 4*INT(l/4)) lg = -1
  IF (k == 1) THEN
    r = r0
  ELSE
    r = r * (m+k-1.0) * (m+k+ip-1.5_dp) / (k-1.0_dp) / (k+ip-1.5_dp)
  END IF
  np = m + 2 * k - 2 + ip
  sud = sud + lg * r * df(k) * dj(np)
  IF (k > nm1 .AND. ABS(sud-sw) < ABS(sud)*eps) EXIT
  sw = sud
END DO

r1d = b0 + a0 * c * sud
RETURN
END SUBROUTINE rmn1



SUBROUTINE sckb(m, n, c, df, ck)

!       ======================================================
!       Purpose: Compute the expansion coefficients of the
!                prolate and oblate spheroidal functions, c2k
!       Input :  m  --- Mode parameter
!                n  --- Mode parameter
!                c  --- Spheroidal parameter
!                DF(k) --- Expansion coefficients dk
!       Output:  CK(k) --- Expansion coefficients ck;
!                          CK(1), CK(2), ... correspond to
!                          c0, c2, ...
!       ======================================================

INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: c
REAL (dp), INTENT(IN)      :: df(200)
REAL (dp), INTENT(OUT)     :: ck(200)

REAL (dp)  :: d1, d2, d3, fac, r, r1, reg, sum, sw
INTEGER    :: i, i1, i2, ip, k, nm

IF (c <= 1.0D-10) c = 1.0D-10
nm = 25 + INT(0.5*(n-m)+c)
ip = 1
IF (n-m == 2*INT((n-m)/2)) ip = 0
reg = 1.0_dp
IF (m+nm > 80) reg = 1.0D-200
fac = -0.5_dp ** m
DO  k = 0, nm - 1
  fac = -fac
  i1 = 2 * k + ip + 1
  r = reg
  DO  i = i1, i1 + 2 * m - 1
    r = r * i
  END DO
  i2 = k + m + ip
  DO  i = i2, i2 + k - 1
    r = r * (i+0.5_dp)
  END DO
  sum = r * df(k+1)
  DO  i = k + 1, nm
    d1 = 2.0_dp * i + ip
    d2 = 2.0_dp * m + d1
    d3 = i + m + ip - 0.5_dp
    r = r * d2 * (d2-1.0_dp) * i * (d3+k) / (d1*(d1-1.0_dp)*(i-k)*d3 )
    sum = sum + r * df(i+1)
    IF (ABS(sw-sum) < ABS(sum)*1.0D-14) EXIT
    sw = sum
  END DO
  r1 = reg
  DO  i = 2, m + k
    r1 = r1 * i
  END DO
  ck(k+1) = fac * sum / r1
END DO
RETURN
END SUBROUTINE sckb



SUBROUTINE sphj(n, x, nm, sj, dj)

!       =======================================================
!       Purpose: Compute spherical Bessel functions jn(x) and
!                their derivatives
!       Input :  x --- Argument of jn(x)
!                n --- Order of jn(x)  ( n = 0,1,תתת )
!       Output:  SJ(n) --- jn(x)
!                DJ(n) --- jn'(x)
!                NM --- Highest order computed
!       Routines called:
!                MSTA1 and MSTA2 for computing the starting
!                point for backward recurrence
!       =======================================================

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: x
INTEGER, INTENT(OUT)       :: nm
REAL (dp), INTENT(OUT)     :: sj(0:n)
REAL (dp), INTENT(OUT)     :: dj(0:n)

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
  m = msta1(x,200)
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
  DO  k = 0, nm
    sj(k) = cs * sj(k)
  END DO
END IF
dj(0) = (COS(x)-SIN(x)/x) / x
DO  k = 1, nm
  dj(k) = sj(k-1) - (k+1) * sj(k) / x
END DO
RETURN
END SUBROUTINE sphj



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



SUBROUTINE rmn2l(m, n, c, x, df, kd, r2f, r2d, id)

!    ========================================================
!    Purpose: Compute prolate and oblate spheroidal radial
!             functions of the second kind for given m, n,
!             c and a large cx
!    Routine called:
!             SPHY for computing the spherical Bessel
!             functions of the second kind
!    ========================================================

INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: c
REAL (dp), INTENT(IN)      :: x
REAL (dp), INTENT(IN)      :: df(200)
INTEGER, INTENT(IN)        :: kd
REAL (dp), INTENT(OUT)     :: r2f
REAL (dp), INTENT(OUT)     :: r2d
INTEGER, INTENT(OUT)       :: id

REAL (dp)  :: sy(0:251), dy(0:251)
REAL (dp)  :: a0, b0, cx, eps, eps1, eps2, r, r0, reg, suc, sud, sw
INTEGER    :: id1, id2, ip, j, k, l, lg, nm, nm1, nm2, nm3, np

eps = 1.0D-14
ip = 1
nm1 = INT((n-m)/2)
IF (n-m == 2*nm1) ip = 0
nm = 25 + nm1 + INT(c)
reg = 1.0_dp
IF (m+nm > 80) reg = 1.0D-200
nm2 = 2 * nm + m
cx = c * x
CALL sphy(nm2, cx, nm3, sy, dy)
r0 = reg
DO  j = 1, 2 * m + ip
  r0 = r0 * j
END DO
r = r0
suc = r * df(1)
sw = suc
DO  k = 2, nm
  r = r * (m+k-1.0) * (m+k+ip-1.5_dp) / (k-1.0_dp) / (k+ip-1.5_dp)
  suc = suc + r * df(k)
  IF (k > nm1 .AND. ABS(suc-sw) < ABS(suc)*eps) EXIT
  sw = suc
END DO

a0 = (1.0_dp-kd/(x*x)) ** (0.5_dp*m) / suc
r2f = 0.0
DO  k = 1, nm
  l = 2 * k + m - n - 2 + ip
  IF (l == 4*INT(l/4)) lg = 1
  IF (l /= 4*INT(l/4)) lg = -1
  IF (k == 1) THEN
    r = r0
  ELSE
    r = r * (m+k-1.0) * (m+k+ip-1.5_dp) / (k-1.0_dp) / (k+ip-1.5_dp)
  END IF
  np = m + 2 * k - 2 + ip
  r2f = r2f + lg * r * (df(k)*sy(np))
  eps1 = ABS(r2f-sw)
  IF (k > nm1 .AND. eps1 < ABS(r2f)*eps) EXIT
  sw = r2f
END DO

id1 = INT(LOG10(eps1/ABS(r2f)+eps))
r2f = r2f * a0
IF (np >= nm3) THEN
  id = 10
  RETURN
END IF
b0 = kd * m / x ** 3 / (1.0 - kd/(x*x)) * r2f
sud = 0.0_dp
DO  k = 1, nm
  l = 2 * k + m - n - 2 + ip
  IF (l == 4*INT(l/4)) lg = 1
  IF (l /= 4*INT(l/4)) lg = -1
  IF (k == 1) THEN
    r = r0
  ELSE
    r = r * (m+k-1.0) * (m+k+ip-1.5_dp) / (k-1.0_dp) / (k+ip-1.5_dp)
  END IF
  np = m + 2 * k - 2 + ip
  sud = sud + lg * r * (df(k)*dy(np))
  eps2 = ABS(sud-sw)
  IF (k > nm1 .AND. eps2 < ABS(sud)*eps) EXIT
  sw = sud
END DO

r2d = b0 + a0 * c * sud
id2 = INT(LOG10(eps2/ABS(sud)+eps))
id = MAX(id1,id2)
RETURN
END SUBROUTINE rmn2l



SUBROUTINE sphy(n, x, nm, sy, dy)

!    ======================================================
!    Purpose: Compute spherical Bessel functions yn(x) and
!             their derivatives
!    Input :  x --- Argument of yn(x) ( x ע 0 )
!             n --- Order of yn(x) ( n = 0,1,תתת )
!    Output:  SY(n) --- yn(x)
!             DY(n) --- yn'(x)
!             NM --- Highest order computed
!    ======================================================

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: x
INTEGER, INTENT(OUT)       :: nm
REAL (dp), INTENT(OUT)     :: sy(0:n)
REAL (dp), INTENT(OUT)     :: dy(0:n)

REAL (dp)  :: f, f0, f1
INTEGER    :: k

nm = n
IF (x < 1.0D-60) THEN
  DO  k = 0, n
    sy(k) = -1.0D+300
    dy(k) = 1.0D+300
  END DO
  RETURN
END IF
sy(0) = -COS(x) / x
sy(1) = (sy(0)-SIN(x)) / x
f0 = sy(0)
f1 = sy(1)
DO  k = 2, n
  f = (2*k-1) * f1 / x - f0
  sy(k) = f
  IF (ABS(f) >= 1.0D+300) EXIT
  f0 = f1
  f1 = f
END DO

nm = k - 1
dy(0) = (SIN(x)+COS(x)/x) / x
DO  k = 1, nm
  dy(k) = sy(k-1) - (k+1) * sy(k) / x
END DO
RETURN
END SUBROUTINE sphy



SUBROUTINE rmn2so(m, n, c, x, cv, df, kd, r2f, r2d)

!    =============================================================
!    Purpose: Compute oblate radial functions of the second kind
!             with a small argument, Rmn(-ic,ix) & Rmn'(-ic,ix)
!    Routines called:
!         (1) SCKB for computing the expansion coefficients c2k
!         (2) KMN for computing the joining factors
!         (3) QSTAR for computing the factor defined in (15.7.3)
!         (4) CBK for computing the the expansion coefficient
!             defined in (15.7.6)
!         (5) GMN for computing the function defined in (15.7.4)
!         (6) RMN1 for computing the radial function of the first
!             kind
!    =============================================================

INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: c
REAL (dp), INTENT(IN)      :: x
REAL (dp), INTENT(IN OUT)  :: cv
REAL (dp), INTENT(IN OUT)  :: df(200)
INTEGER, INTENT(IN OUT)    :: kd
REAL (dp), INTENT(OUT)     :: r2f
REAL (dp), INTENT(OUT)     :: r2d

REAL (dp)  :: bk(200), ck(200), dn(200)
REAL (dp)  :: ck1, ck2, eps, gd, gf, h0, pi, qs, qt, r1d, r1f, sum, sw
INTEGER    :: ip, j, nm

IF (ABS(df(1)) <= 1.0D-280) THEN
  r2f = 1.0D+300
  r2d = 1.0D+300
  RETURN
END IF
eps = 1.0D-14
pi = 3.141592653589793_dp
nm = 25 + INT((n-m)/2+c)
ip = 1
IF (n-m == 2*INT((n-m)/2)) ip = 0
CALL sckb(m, n, c, df, ck)
CALL kmn(m, n, c, cv, kd, df, dn, ck1, ck2)
CALL qstar(m, n, c, ck, ck1, qs, qt)
CALL cbk(m, n, c, cv, qt, ck, bk)
IF (x == 0.0_dp) THEN
  sum = 0.0_dp
  DO  j = 1, nm
    sum = sum + ck(j)
    IF (ABS(sum-sw) < ABS(sum)*eps) EXIT
    sw = sum
  END DO

  IF (ip == 0) THEN
    r1f = sum / ck1
    r2f = -0.5_dp * pi * qs * r1f
    r2d = qs * r1f + bk(1)
  ELSE IF (ip == 1) THEN
    r1d = sum / ck1
    r2f = bk(1)
    r2d = -0.5_dp * pi * qs * r1d
  END IF
  RETURN
ELSE
  CALL gmn(m, n, c, x, bk, gf, gd)
  CALL rmn1(m, n, c, x, df, kd, r1f, r1d)
  h0 = ATAN(x) - 0.5_dp * pi
  r2f = qs * r1f * h0 + gf
  r2d = qs * (r1d*h0 + r1f/(1.0_dp + x*x)) + gd
END IF
RETURN
END SUBROUTINE rmn2so



SUBROUTINE qstar(m, n, c, ck, ck1, qs, qt)

!    =========================================================
!    Purpose: Compute Q*mn(-ic) for oblate radial functions
!             with a small argument
!    =========================================================

INTEGER, INTENT(IN)     :: m
INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: c
REAL (dp), INTENT(IN)   :: ck(200)
REAL (dp), INTENT(IN)   :: ck1
REAL (dp), INTENT(OUT)  :: qs
REAL (dp), INTENT(OUT)  :: qt

REAL (dp)  :: ap(200), qs0, r, s, sk
INTEGER    :: i, ip, k, l

ip = 1
IF (n-m == 2*INT((n-m)/2)) ip = 0
r = 1.0_dp / ck(1) ** 2
ap(1) = r
DO  i = 1, m
  s = 0.0_dp
  DO  l = 1, i
    sk = 0.0_dp
    DO  k = 0, l
      sk = sk + ck(k+1) * ck(l-k+1)
    END DO
    s = s + sk * ap(i-l+1)
  END DO
  ap(i+1) = -r * s
END DO
qs0 = ap(m+1)
DO  l = 1, m
  r = 1.0_dp
  DO  k = 1, l
    r = r * (2*k+ip) * (2*k-1+ip) / (2.0_dp*k) ** 2
  END DO
  qs0 = qs0 + ap(m-l+1) * r
END DO
qs = (-1) ** ip * ck1 * (ck1*qs0) / c
qt = -2.0_dp / ck1 * qs
RETURN
END SUBROUTINE qstar



SUBROUTINE cbk(m, n, c, cv, qt, ck, bk)

!    =====================================================
!    Purpose: Compute coefficient Bk's for oblate radial
!             functions with a small argument
!    =====================================================

INTEGER, INTENT(IN)     :: m
INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: c
REAL (dp), INTENT(IN)   :: cv
REAL (dp), INTENT(IN)   :: qt
REAL (dp), INTENT(IN)   :: ck(200)
REAL (dp), INTENT(OUT)  :: bk(200)

REAL (dp)  :: u(200), v(200), w(200)
REAL (dp)  :: eps, r1, s1, sw, t
INTEGER    :: i, i1, ip, j, k, n2, nm

eps = 1.0D-14
ip = 1
IF (n-m == 2*INT((n-m)/2)) ip = 0
nm = 25 + INT(0.5*(n-m)+c)
u(1) = 0.0_dp
n2 = nm - 2
DO  j = 2, n2
  u(j) = c * c
END DO
DO  j = 1, n2
  v(j) = (2*j-1-ip) * (2*(j-m)-ip) + m * (m-1) - cv
END DO
DO  j = 1, nm - 1
  w(j) = (2*j-ip) * (2*j+1-ip)
END DO
IF (ip == 0) THEN
  DO  k = 0, n2 - 1
    s1 = 0.0_dp
    i1 = k - m + 1
    DO  i = i1, nm
      IF (i >= 0) THEN
        r1 = 1.0_dp
        DO  j = 1, k
          r1 = r1 * (i+m-j) / j
        END DO
        s1 = s1 + ck(i+1) * (2*i+m) * r1
        IF (ABS(s1-sw) < ABS(s1)*eps) EXIT
        sw = s1
      END IF
    END DO
    bk(k+1) = qt * s1
  END DO
ELSE IF (ip == 1) THEN
  DO  k = 0, n2 - 1
    s1 = 0.0_dp
    i1 = k - m + 1
    DO  i = i1, nm
      IF (i >= 0) THEN
        r1 = 1.0_dp
        DO  j = 1, k
          r1 = r1 * (i+m-j) / j
        END DO
        IF (i > 0) s1 = s1 + ck(i) * (2*i+m-1) * r1
        s1 = s1 - ck(i+1) * (2*i+m) * r1
        IF (ABS(s1-sw) < ABS(s1)*eps) EXIT
        sw = s1
      END IF
    END DO
    bk(k+1) = qt * s1
  END DO
END IF
w(1) = w(1) / v(1)
bk(1) = bk(1) / v(1)
DO  k = 2, n2
  t = v(k) - w(k-1) * u(k)
  w(k) = w(k) / t
  bk(k) = (bk(k)-bk(k-1)*u(k)) / t
END DO
DO  k = n2 - 1, 1, -1
  bk(k) = bk(k) - w(k) * bk(k+1)
END DO
RETURN
END SUBROUTINE cbk



SUBROUTINE gmn(m, n, c, x, bk, gf, gd)

!    ===========================================================
!    Purpose: Compute gmn(-ic,ix) and its derivative for oblate
!             radial functions with a small argument
!    ===========================================================

INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: c
REAL (dp), INTENT(IN)      :: x
REAL (dp), INTENT(IN)      :: bk(200)
REAL (dp), INTENT(OUT)     :: gf
REAL (dp), INTENT(OUT)     :: gd

REAL (dp)  :: eps, gd0, gd1, gf0, gw, xm
INTEGER    :: ip, k, nm

eps = 1.0D-14
ip = 1
IF (n-m == 2*INT((n-m)/2)) ip = 0
nm = 25 + INT(0.5*(n-m)+c)
xm = (1.0_dp + x*x) ** (-0.5_dp*m)
gf0 = 0.0_dp
DO  k = 1, nm
  gf0 = gf0 + bk(k) * x ** (2*k-2)
  IF (ABS((gf0-gw)/gf0) < eps .AND. k >= 10) EXIT
  gw = gf0
END DO

gf = xm * gf0 * x ** (1-ip)
gd1 = -m * x / (1.0_dp + x*x) * gf
gd0 = 0.0_dp
DO  k = 1, nm
  IF (ip == 0) THEN
    gd0 = gd0 + (2*k-1) * bk(k) * x ** (2*k-2)
  ELSE
    gd0 = gd0 + 2 * k * bk(k+1) * x ** (2*k-1)
  END IF
  IF (ABS((gd0-gw)/gd0) < eps .AND. k >= 10) EXIT
  gw = gd0
END DO

gd = gd1 + xm * gd0
RETURN
END SUBROUTINE gmn



SUBROUTINE kmn(m, n, c, cv, kd, df, dn, ck1, ck2)

!    ===================================================
!    Purpose: Compute the expansion coefficients of the
!             prolate and oblate spheroidal functions
!             and joining factors
!    ===================================================

INTEGER, INTENT(IN)     :: m
INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: c
REAL (dp), INTENT(IN)   :: cv
INTEGER, INTENT(IN)     :: kd
REAL (dp), INTENT(IN)   :: df(200)
REAL (dp), INTENT(OUT)  :: dn(200)
REAL (dp), INTENT(OUT)  :: ck1
REAL (dp), INTENT(OUT)  :: ck2

REAL (dp)  :: u(200), v(200), w(200), tp(200), rk(200)
REAL (dp)  :: cs, dnp, g0, gk0, gk1, gk2, gk3, r, r1, r2, r3, r4, r5,  &
              sa0, sb0, su0, sw, t
INTEGER    :: i, ip, j, k, l, nm, nm1, nn

nm = 25 + INT(0.5*(n-m) + c)
nn = nm + m
cs = c * c * kd
ip = 1
IF (n-m == 2*INT((n-m)/2)) ip = 0
DO  i = 1, nn + 3
  IF (ip == 0) k = -2 * (i-1)
  IF (ip == 1) k = -(2*i-3)
  gk0 = 2 * m + k
  gk1 = (m+k) * (m+k+1)
  gk2 = 2 * (m+k) - 1
  gk3 = 2 * (m+k) + 3
  u(i) = gk0 * (gk0-1.0_dp) * cs / (gk2*(gk2+2.0_dp))
  v(i) = gk1 - cv + (2.0_dp*(gk1-m*m)-1.0_dp) * cs / (gk2*gk3)
  w(i) = (k+1) * (k+2) * cs / ((gk2+2.0_dp)*gk3)
END DO
DO  k = 1, m
  t = v(m+1)
  DO  l = 0, m - k - 1
    t = v(m-l) - w(m-l+1) * u(m-l) / t
  END DO
  rk(k) = -u(k) / t
END DO
r = 1.0_dp
DO  k = 1, m
  r = r * rk(k)
  dn(k) = df(1) * r
END DO
tp(nn) = v(nn+1)
DO  k = nn - 1, m + 1, -1
  tp(k) = v(k+1) - w(k+2) * u(k+1) / tp(k+1)
  IF (k > m+1) rk(k) = -u(k) / tp(k)
END DO
IF (m == 0) dnp = df(1)
IF (m /= 0) dnp = dn(m)
dn(m+1) = (-1) ** ip * dnp * cs / ((2*m-1)*(2*m+1 - 4*ip)* tp(m+1))
DO  k = m + 2, nn
  dn(k) = rk(k) * dn(k-1)
END DO
r1 = 1.0_dp
DO  j = 1, (n+m+ip) / 2
  r1 = r1 * (j + 0.5_dp*(n+m+ip))
END DO
nm1 = (n-m) / 2
r = 1.0_dp
DO  j = 1, 2 * m + ip
  r = r * j
END DO
su0 = r * df(1)
DO  k = 2, nm
  r = r * (m+k-1.0) * (m+k+ip-1.5_dp) / (k-1.0_dp) / (k+ip-1.5_dp)
  su0 = su0 + r * df(k)
  IF (k > nm1 .AND. ABS((su0-sw)/su0) < 1.0D-14) EXIT
  sw = su0
END DO

IF (kd /= 1) THEN
  r2 = 1.0_dp
  DO  j = 1, m
    r2 = 2.0_dp * c * r2 * j
  END DO
  r3 = 1.0_dp
  DO  j = 1, (n-m-ip) / 2
    r3 = r3 * j
  END DO
  sa0 = (2*(m+ip)+1) * r1 / (2.0**n*c**ip*r2*r3*df(1))
  ck1 = sa0 * su0
  IF (kd == -1) RETURN
END IF
r4 = 1.0_dp
DO  j = 1, (n-m-ip) / 2
  r4 = 4.0_dp * r4 * j
END DO
r5 = 1.0_dp
DO  j = 1, m
  r5 = r5 * (j+m) / c
END DO
g0 = dn(m)
IF (m == 0) g0 = df(1)
sb0 = (ip+1) * c ** (ip+1) / (2.0*ip*(m-2.0)+1.0) / (2.0*m-1.0)
ck2 = (-1) ** ip * sb0 * r4 * r5 * g0 / r1 * su0
RETURN
END SUBROUTINE kmn



SUBROUTINE segv(m, n, c, kd, cv, eg)

!    =========================================================
!    Purpose: Compute the characteristic values of spheroidal
!             wave functions
!    Input :  m  --- Mode parameter
!             n  --- Mode parameter
!             c  --- Spheroidal parameter
!             KD --- Function code
!                    KD=1 for Prolate; KD=-1 for Oblate
!    Output:  CV --- Characteristic value for given m, n and c
!             EG(L) --- Characteristic value for mode m and n'
!                       ( L = n' - m + 1 )
!    =========================================================

INTEGER, INTENT(IN)     :: m
INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: c
INTEGER, INTENT(IN)     :: kd
REAL (dp), INTENT(OUT)  :: cv
REAL (dp), INTENT(OUT)  :: eg(200)

REAL (dp)  :: b(100), h(100), d(300), e(300), f(300), cv0(100), a(300), g(300)
INTEGER    :: i, icm, j, k, k1, l, nm, nm1
REAL (dp)  :: cs, d2k, dk0, dk1, dk2, s, t, t1, x1, xa, xb

IF (c < 1.0D-10) THEN
  DO  i = 1, n
    eg(i) = (i+m) * (i+m-1)
  END DO
  GO TO 120
END IF
icm = (n-m+2) / 2
nm = 10 + INT(0.5*(n-m) + c)
cs = c * c * kd
DO  l = 0, 1
  DO  i = 1, nm
    IF (l == 0) k = 2 * (i-1)
    IF (l == 1) k = 2 * i - 1
    dk0 = m + k
    dk1 = m + k + 1.0_dp
    dk2 = 2 * (m+k)
    d2k = 2 * m + k
    a(i) = (d2k+2.0) * (d2k+1.0) / ((dk2+3.0)*(dk2+5.0)) * cs
    d(i) = dk0 * dk1 + (2.0*dk0*dk1-2.0*m*m-1.0) / ((dk2-1.0)*(dk2 +3.0)) * cs
    g(i) = k * (k-1.0) / ((dk2-3.0)*(dk2-1.0)) * cs
  END DO
  DO  k = 2, nm
    e(k) = SQRT(a(k-1)*g(k))
    f(k) = e(k) * e(k)
  END DO
  f(1) = 0.0_dp
  e(1) = 0.0_dp
  xa = d(nm) + ABS(e(nm))
  xb = d(nm) - ABS(e(nm))
  nm1 = nm - 1
  DO  i = 1, nm1
    t = ABS(e(i)) + ABS(e(i+1))
    t1 = d(i) + t
    IF (xa < t1) xa = t1
    t1 = d(i) - t
    IF (t1 < xb) xb = t1
  END DO
  DO  i = 1, icm
    b(i) = xa
    h(i) = xb
  END DO
  DO  k = 1, icm
    DO  k1 = k, icm
      IF (b(k1) < b(k)) THEN
        b(k) = b(k1)
        EXIT
      END IF
    END DO
    IF (k /= 1 .AND. h(k) < h(k-1)) h(k) = h(k-1)

    80 x1 = (b(k)+h(k)) / 2.0_dp
    cv0(k) = x1
    IF (ABS((b(k)-h(k))/x1) >= 1.0D-14) THEN
      j = 0
      s = 1.0_dp
      DO  i = 1, nm
        IF (s == 0.0_dp) s = s + 1.0D-30
        t = f(i) / s
        s = d(i) - t - x1
        IF (s < 0.0_dp) j = j + 1
      END DO
      IF (j < k) THEN
        h(k) = x1
      ELSE
        b(k) = x1
        IF (j >= icm) THEN
          b(icm) = x1
        ELSE
          IF (h(j+1) < x1) h(j+1) = x1
          IF (x1 < b(j)) b(j) = x1
        END IF
      END IF
      GO TO 80
    END IF
    cv0(k) = x1
    IF (l == 0) eg(2*k-1) = cv0(k)
    IF (l == 1) eg(2*k) = cv0(k)
  END DO
END DO

120 cv = eg(n-m+1)
RETURN
END SUBROUTINE segv
 
END MODULE rswfo_func
 
 
 
PROGRAM mrswfo
USE rswfo_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:45

!       =============================================================
!       Purpose: This program computes the radial oblate spheriodal
!                functions of the first and second kinds, and their
!                derivatives using subroutine RSWFO
!       Input :  m  --- Mode parameter, m = 0,1,2,...
!                n  --- Mode parameter, n = m,m+1,m+2,...
!                c  --- Spheroidal parameter
!                x  --- Argument (x ע 0)
!                cv --- Characteristic value
!                KF --- Function code
!                       KF=1 for the first kind
!                       KF=2 for the second kind
!                       KF=3 for both the first and second kinds
!       Output:  R1F --- Radial function of the first kind
!                R1D --- Derivative of the radial function of
!                        the first kind
!                R2F --- Radial function of the second kind
!                R2D --- Derivative of the radial function of
!                        the second kind
!       Example:
!                KD=-1, m = 2, n = 3, c = 5.0 and cv = 2.10980581604

!    x  R23(1)(-ic,ix) R23(1)'(-ic,ix) R23(2)(-ic,ix) R23(2)'(-ic,ix)
!  ------------------------------------------------------------------
!   0.0      0.0000000 (-1) 4.9911346  (-1)-4.0071049  (-1) 3.3065682
!   0.5 (-1) 2.0215966 (-1) 1.9559661  (-1)-1.5154835  (-1) 6.4482526
!   1.0 (-1) 1.0526737 (-1)-5.4010624  (-1) 1.3065731  (-1) 2.7958491
!   1.5 (-1)-1.1611246 (-2)-8.9414690  (-2) 3.7405936  (-1)-5.0118500
!   5.0 (-2) 3.6463251 (-2)-8.3283065  (-2) 1.5549615  (-1) 1.7544481
!       ==============================================================

REAL (dp)  :: c, cv, eg(200), r1f, r1d, r2f, r2d, x
INTEGER    :: kf, m, n

WRITE (*,*) 'Please enter KF: '
READ (*,*) kf
WRITE (*,5000) kf
WRITE (*,*) 'Please enter m, n, c and x ( x ע 0 ): '
READ (*,*) m, n, c, x
CALL segv(m, n, c, -1, cv, eg)
WRITE (*,5100) m, n, c, cv, x
WRITE (*,*)
CALL rswfo(m, n, c, x, cv, kf, r1f, r1d, r2f, r2d)
WRITE (*,*) '   x    Rmn(1)(-ic,ix)  Rmn''(1)(-IC,IX)  ',  &
    ' Rmn(2)(-ic,ix)  Rmn''(2)(-IC,IX)'
WRITE (*,*) '  ---------------------------------------------',  &
    '---------------------------'
IF (kf == 1) THEN
  WRITE (*,5200) x, r1f, r1d
ELSE IF (kf == 2) THEN
  WRITE (*,5300) x, r2f, r2d
ELSE IF (kf == 3) THEN
  WRITE (*,5200) x, r1f, r1d, r2f, r2d
END IF
IF (kf == 3) THEN
  WRITE (*,5400) r1f * r2d - r2f * r1d, 1.0_dp / (c*(x*x + 1.0_dp))
  WRITE (*,5500)
END IF
STOP

5000 FORMAT (' KF=',i3)
5100 FORMAT (' m=', i2, ',   n=', i2, ',   c=', f5.1, ',    cv =',  &
             f18.10, ',  x=', f5.2)
5200 FORMAT (' ', f5.1, 4g17.8)
5300 FORMAT (' ', f5.1, t41, 4g17.8)
5400 FORMAT (/' Wronskian check:'/  &
             ' Computed value =', g17.8, '     Exact value =', g17.8)
5500 FORMAT (/' Caution: This check is not accurate if it involves'/  &
             '          the subtraction of two similar numbers')
END PROGRAM mrswfo
