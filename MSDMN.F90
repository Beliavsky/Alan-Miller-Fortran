MODULE sdmn_func
 
! From the book "Computation of Special Functions"
!      by Shanjie Zhang and Jianming Jin
!   Copyright 1996 by John Wiley & Sons, Inc.
! The authors state:
!   "However, we give permission to the reader who purchases this book
!    to incorporate any of these programs into his or her programs
!    provided that the copyright is acknowledged."

! Latest revision - 14 January 2002
! Corrections by Alan Miller (amiller @ bigpond.net.au)
! Variable SW was undefined in routine SDMN.

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
 
CONTAINS


SUBROUTINE sdmn(m, n, c, cv, kd, df)

!       =====================================================
!       Purpose: Compute the expansion coefficients of the
!                prolate and oblate spheroidal functions, dk
!       Input :  m  --- Mode parameter
!                n  --- Mode parameter
!                c  --- Spheroidal parameter
!                cv --- Characteristic value
!                KD --- Function code
!                       KD=1 for prolate; KD=-1 for oblate
!       Output:  DF(k) --- Expansion coefficients dk;
!                          DF(1), DF(2),... correspond to
!                          d0, d2,... for even n-m and d1,
!                          d3,... for odd n-m
!       =====================================================

INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: c
REAL (dp), INTENT(IN)      :: cv
INTEGER, INTENT(IN)        :: kd
REAL (dp), INTENT(OUT)     :: df(200)

REAL (dp)  :: a(200), d(200), g(200)
REAL (dp)  :: cs, d2k, dk0, dk1, dk2, f, f0, f1, f2, fl, fs,  &
              r1, r3, r4, s0, su1, su2, sw
INTEGER    :: i, ip, j, k, k1, kb, nm

nm = 25 + INT(0.5*(n-m)+c)
fl = 0.0_dp
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
      fs = -((d(2)-cv)*f2 + g(2)*f1) / a(2)
    ELSE
      df(2) = f2
      DO  j = 3, kb + 1
        f = -((d(j-1)-cv)*f2 + g(j-1)*f1) / a(j-1)
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
  r1 = -r1 * (k+m+ip-1.5_dp) / (k-1)
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

!       =========================================================
!       Purpose: Compute the characteristic values of spheroidal
!                wave functions
!       Input :  m  --- Mode parameter
!                n  --- Mode parameter
!                c  --- Spheroidal parameter
!                KD --- Function code
!                       KD=1 for Prolate; KD=-1 for Oblate
!       Output:  CV --- Characteristic value for given m, n and c
!                EG(L) --- Characteristic value for mode m and n'
!                          ( L = n' - m + 1 )
!       =========================================================

INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: c
INTEGER, INTENT(IN)        :: kd
REAL (dp), INTENT(OUT)     :: cv
REAL (dp), INTENT(OUT)     :: eg(200)

REAL (dp)  :: b(100), h(100), d(300), e(300), f(300), cv0(100), a(300), g(300)
REAL (dp)  :: cs, d2k, dk0, dk1, dk2, s, t, t1, x1, xa, xb
INTEGER    :: i, icm, j, k, k1, l, nm, nm1

IF (c < 1.0D-10) THEN
  DO  i = 1, n
    eg(i) = (i+m) * (i+m-1)
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
 
END MODULE sdmn_func
 
 
 
PROGRAM msdmn
USE sdmn_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:46

!       ===========================================================
!       Purpose: This program computes the expansion coefficients
!                of the prolate and oblate spheroidal functions,
!                dk, using subroutine SDMN
!       Input :  m  --- Mode parameter
!                n  --- Mode parameter
!                c  --- Spheroidal parameter
!                cv --- Characteristic value
!                KD --- Function code
!                       KD=1 for prolate; KD=-1 for oblate
!       Output:  DF(k) --- Expansion coefficients dk;
!                          DF(1), DF(2),... correspond to
!                          d0, d2,... for even n-m and d1,
!                          d3,... for odd n-m
!       Example: Compute the first 12 expansion coefficients for
!                KD= 1, m=2, n=2, c=3.0 and cv=7.1511005241; and
!                KD=-1, m=2, n=2, c=3.0 and cv=4.5264604622

!                Coefficients of Prolate and oblate functions

!                  r          dr(c)             dr(-ic)
!                -------------------------------------------
!                  0     .9237882817D+00    .1115434000D+01
!                  2    -.2901607696D-01    .4888489020D-01
!                  4     .8142246173D-03    .1600845667D-02
!                  6    -.1632270292D-04    .3509183384D-04
!                  8     .2376699010D-06    .5416293446D-06
!                 10    -.2601391701D-08    .6176624069D-08
!                 12     .2209142844D-10    .5407431236D-10
!                 14    -.1494812074D-12    .3745889118D-12
!                 16     .8239302207D-15    .2103624480D-14
!                 18    -.3768260778D-17    .9768323113D-17
!                 20     .1452384658D-19    .3812753620D-19
!                 22    -.4780280430D-22    .1268321726D-21
!       ===========================================================

REAL (dp)  :: c, cv, df(200), eg(200)
INTEGER    :: j, k, kd, m, n, nm

WRITE (*,*) 'Please enter KD, m, n and c '
READ (*,*) kd, m, n, c
CALL segv(m, n, c, kd, cv, eg)
WRITE (*,5100) kd, m, n, c, cv
CALL sdmn(m, n, c, cv, kd, df)
WRITE (*,*)
IF (kd == 1) THEN
  WRITE (*,*) 'Coefficients of Prolate function'
  WRITE (*,*)
  WRITE (*,*) '   r             dr(c)'
ELSE
  WRITE (*,*) 'Coefficients of Oblate function'
  WRITE (*,*)
  WRITE (*,*) '   r            dr(-ic)'
END IF
WRITE (*,*) '----------------------------'
nm = 25 + INT(0.5*(n-m)+c)
DO  k = 1, nm
  IF (n-m == 2*INT((n-m)/2)) THEN
    j = 2 * (k-1)
  ELSE
    j = 2 * k - 1
  END IF
  WRITE (*,5000) j, df(k)
END DO
STOP

5000 FORMAT(t3, i3, '    ', g18.10)
5100 FORMAT(' KD=', i3, ',  m=', i3, ',  n=', i3, ',  c=', f5.1, ',  cv =', f18.10)
END PROGRAM msdmn
