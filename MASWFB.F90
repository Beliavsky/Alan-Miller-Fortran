MODULE aswfb_func
 
! From the book "Computation of Special Functions"
!      by Shanjie Zhang and Jianming Jin
!   Copyright 1996 by John Wiley & Sons, Inc.
! The authors state:
!   "However, we give permission to the reader who purchases this book
!    to incorporate any of these programs into his or her programs
!    provided that the copyright is acknowledged."

! Latest revision - 27 December 2001
! Corrections by Alan Miller (amiller @ bigpond.net.au)
! Variables sw & fl were used without values assigned to them

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
 
CONTAINS


SUBROUTINE aswfb(m, n, c, x, kd, cv, s1f, s1d)

!    ===========================================================
!    Purpose: Compute the prolate and oblate spheroidal angular
!             functions of the first kind and their derivatives
!    Input :  m  --- Mode parameter,  m = 0,1,2,...
!             n  --- Mode parameter,  n = m,m+1,...
!             c  --- Spheroidal parameter
!             x  --- Argument of angular function, |x| < 1.0
!             KD --- Function code
!                    KD=1 for prolate;  KD=-1 for oblate
!             cv --- Characteristic value
!    Output:  S1F --- Angular function of the first kind
!             S1D --- Derivative of the angular function of
!                     the first kind
!    Routines called:
!         (1) SDMN for computing expansion coefficients dk
!         (2) LPMNS for computing associated Legendre function
!             of the first kind Pmn(x)
!    ===========================================================

INTEGER, INTENT(IN)     :: m
INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: c
REAL (dp), INTENT(IN)   :: x
INTEGER, INTENT(IN)     :: kd
REAL (dp), INTENT(IN)   :: cv
REAL (dp), INTENT(OUT)  :: s1f
REAL (dp), INTENT(OUT)  :: s1d

REAL (dp)  :: df(200), pm(0:251), pd(0:251)
REAL (dp)  :: eps, su1, sw
INTEGER    :: ip, k, mk, nm, nm2

eps = 1.0D-14
ip = 1
IF (n-m == 2*INT((n-m)/2)) ip = 0
nm = 25 + INT((n-m)/2+c)
nm2 = 2 * nm + m
CALL sdmn(m, n, c, cv, kd, df)
CALL lpmns(m, nm2, x, pm, pd)
su1 = 0.0_dp
sw = 0.0_dp
DO  k = 1, nm
  mk = m + 2 * (k-1) + ip
  su1 = su1 + df(k) * pm(mk)
  IF (ABS(sw-su1) < ABS(su1)*eps) EXIT
  sw = su1
END DO

s1f = (-1) ** m * su1
su1 = 0.0_dp
DO  k = 1, nm
  mk = m + 2 * (k-1) + ip
  su1 = su1 + df(k) * pd(mk)
  IF (ABS(sw-su1) < ABS(su1)*eps) EXIT
  sw = su1
END DO

s1d = (-1) ** m * su1
RETURN
END SUBROUTINE aswfb



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
!                       _dp, d2, ... for even n-m and d1,
!                       d3, ... for odd n-m
!    =====================================================

INTEGER, INTENT(IN)     :: m
INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: c
REAL (dp), INTENT(IN)   :: cv
INTEGER, INTENT(IN)     :: kd
REAL (dp), INTENT(OUT)  :: df(200)

REAL (dp)  :: a(200), d(200), g(200)
REAL (dp)  :: cs , d2k, dk0, dk1, dk2, f, f0, f1, f2, fl, fs, r1, r3, r4,  &
              s0, su1, su2, sw
INTEGER    :: i, ip, j, k, kb, nm

nm = 25 + INT(0.5*(n-m)+c)
IF (c < 1.0D-10) THEN
  df(1:nm) = 0.0_dp
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
  d(i) = dk0 * dk1 + (2.0*dk0*dk1-2.0*m*m-1.0) / ((dk2-1.0)*(dk2+ 3.0)) * cs
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
      df(k:nm) = df(k:nm) * 1.0D-100
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
      fs = -((d(2)-cv)*f2 + g(2)*f1) / a(2)
    ELSE
      df(2) = f2
      DO  j = 3, kb + 1
        f = -((d(j-1)-cv)*f2 + g(j-1)*f1) / a(j-1)
        IF (j <= kb) df(j) = f
        IF (ABS(f) > 1.0D+100) THEN
          df(1:j) = df(1:j) * 1.0D-100
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
sw = 0.0_dp
DO  k = kb + 1, nm
  IF (k /= 1) r1 = -r1 * (k+m+ip-1.5_dp) / (k-1.0_dp)
  su2 = su2 + r1 * df(k)
  IF (ABS(sw-su2) < ABS(su2)*1.0D-14) EXIT
  sw = su2
END DO

r3 = 1.0_dp
DO  j = 1, (m+n+ip) / 2
  r3 = r3 * (j+0.5_dp*(n+m+ip))
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
REAL (dp)  :: cs, d2k, dk0, dk1, dk2, s, t, t1, x1, xa, xb
INTEGER    :: i, icm, j, k, k1, l, nm, nm1

IF (c < 1.0D-10) THEN
  DO  i = 1, n
    eg(i) = (i+m) * (i+m-1.0_dp)
  END DO
  GO TO 120
END IF
icm = (n-m+2) / 2
nm = 10 + INT(0.5*(n-m)+c)
cs = c * c * kd
DO  l = 0, 1
  DO  i = 1, nm
    IF (l == 0) k = 2 * (i-1)
    IF (l == 1) k = 2 * i - 1
    dk0 = m + k
    dk1 = m + k + 1
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



SUBROUTINE lpmns(m, n, x, pm, pd)

!    ========================================================
!    Purpose: Compute associated Legendre functions Pmn(x)
!             and Pmn'(x) for a given order
!    Input :  x --- Argument of Pmn(x)
!             m --- Order of Pmn(x),  m = 0,1,2,...,n
!             n --- Degree of Pmn(x), n = 0,1,2,...,N
!    Output:  PM(n) --- Pmn(x)
!             PD(n) --- Pmn'(x)
!    ========================================================

INTEGER, INTENT(IN)     :: m
INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: pm(0:n)
REAL (dp), INTENT(OUT)  :: pd(0:n)

REAL (dp)  :: pm0, pm1, pm2, pmk, x0
INTEGER    :: k

DO  k = 0, n
  pm(k) = 0.0_dp
  pd(k) = 0.0_dp
END DO
IF (ABS(x) == 1.0_dp) THEN
  DO  k = 0, n
    IF (m == 0) THEN
      pm(k) = 1.0_dp
      pd(k) = 0.5_dp * k * (k+1.0)
      IF (x < 0.0) THEN
        pm(k) = (-1) ** k * pm(k)
        pd(k) = (-1) ** (k+1) * pd(k)
      END IF
    ELSE IF (m == 1) THEN
      pd(k) = 1.0D+300
    ELSE IF (m == 2) THEN
      pd(k) = -0.25_dp * (k+2.0) * (k+1.0) * k * (k-1.0)
      IF (x < 0.0) pd(k) = (-1) ** (k+1) * pd(k)
    END IF
  END DO
  RETURN
END IF
x0 = ABS(1.0_dp-x*x)
pm0 = 1.0_dp
pmk = pm0
DO  k = 1, m
  pmk = (2.0_dp*k-1.0_dp) * SQRT(x0) * pm0
  pm0 = pmk
END DO
pm1 = (2.0_dp*m+1.0_dp) * x * pm0
pm(m) = pmk
pm(m+1) = pm1
DO  k = m + 2, n
  pm2 = ((2.0_dp*k-1.0_dp)*x*pm1-(k+m-1.0_dp)*pmk) / (k-m)
  pm(k) = pm2
  pmk = pm1
  pm1 = pm2
END DO
pd(0) = ((1.0_dp-m)*pm(1)-x*pm(0)) / (x*x-1.0)
DO  k = 1, n
  pd(k) = (k*x*pm(k)-(k+m)*pm(k-1)) / (x*x-1.0_dp)
END DO
RETURN
END SUBROUTINE lpmns
 
END MODULE aswfb_func
 
 
 
PROGRAM maswfb
USE aswfb_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:33

!    ============================================================
!    Purpose: This program computes the prolate and oblate
!             spheroidal angular functions of the first kind
!             and their derivatives using subroutine ASWFB
!    Input :  m  --- Mode parameter,  m = 0,1,2,...
!             n  --- Mode parameter,  n = m,m+1,...
!             c  --- Spheroidal parameter
!             x  --- Argument of angular function, |x| ó 1.0
!             KD --- Function code
!                    KD=1 for prolate;  KD=-1 for oblate
!             cv --- Characteristic value
!    Output:  S1F --- Angular function of the first kind
!             S1D --- Derivative of the angular function of
!                     the first kind
!    Examples:
!            KD = 1, m = 2, n = 3, c = 3.0 and cv = 14.8277782138
!               x         Smn(c,x)            Smn'(c,x)
!             --------------------------------------------
!              0.2      .28261309D+01       .12418631D+02
!              0.5      .49938554D+01       .92761604D+00
!              0.8      .31693975D+01      -.12646552D+02

!            KD =-1, m = 2, n = 3, c = 3.0 and cv = 8.8093939208
!               x         Smn(-ic,x)         Smn'(-ic,x)
!             --------------------------------------------
!              0.2      .29417848D+01       .14106305D+02
!              0.5      .64138827D+01       .76007194D+01
!              0.8      .60069873D+01      -.14387479D+02
!    ============================================================

REAL (dp)  :: c, cv, eg(200), s1f, s1d, x
INTEGER    :: i, kd, m, n

WRITE (*,*) 'Please enter KD, m, n and c: '
READ (*,*) kd, m, n, c
WRITE (*,5000) kd, m, n, c
CALL segv(m, n, c, kd, cv, eg)
WRITE (*,5100) cv
WRITE (*,*)
IF (kd == 1) THEN
  WRITE (*,*) '    x         Smn(c,x)            Smn''(C,X)'
ELSE IF (kd == -1) THEN
  WRITE (*,*) '    x         Smn(-ic,x)         Smn''(-IC,X)'
END IF
WRITE (*,*) '  --------------------------------------------'
DO  i = 0, 20
  x = -1.0_dp + 0.1_dp * i
  CALL aswfb(m, n, c, x, kd, cv, s1f, s1d)
  WRITE (*,5200) x, s1f, s1d
END DO
STOP

5000 FORMAT (' KD =', i2, ', ', 'm =', i2, ', ', 'n =', i2, ', ', 'c =', f5.1)
5100 FORMAT ('  cv =', f18.10)
5200 FORMAT (' ', f5.1, 2g20.8)
END PROGRAM maswfb
