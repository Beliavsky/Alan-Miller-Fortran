MODULE scka_func
 
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


SUBROUTINE scka(m, n, c, cv, kd, ck)

!    ======================================================
!    Purpose: Compute the expansion coefficients of the
!             prolate and oblate spheroidal functions, c2k
!    Input :  m  --- Mode parameter
!             n  --- Mode parameter
!             c  --- Spheroidal parameter
!             cv --- Characteristic value
!             KD --- Function code
!                    KD=1 for prolate; KD=-1 for oblate
!    Output:  CK(k) --- Expansion coefficients ck;
!                       CK(1), CK(2),... correspond to
!                       c0, c2,...
!    ======================================================

INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: c
REAL (dp), INTENT(IN)      :: cv
INTEGER, INTENT(IN)        :: kd
REAL (dp), INTENT(OUT)     :: ck(200)

REAL (dp)  :: cs, f, f0, f1, f2, fl, fs, r1, r2, s0, su1, su2
INTEGER    :: ip, j, k, k1, kb, nm

IF (c <= 1.0D-10) c = 1.0D-10
nm = 25 + INT((n-m)/2+c)
cs = c * c * kd
ip = 1
IF (n-m == 2*INT((n-m)/2)) ip = 0
fs = 1.0D0
f1 = 0.0D0
f0 = 1.0D-100
kb = 0
ck(nm+1) = 0.0D0
DO  k = nm, 1, -1
  f = (((2*k+m+ip)*(2*k+m+1+ip)-cv+cs)*f0 - 4*(k+1)*(k+m+1)*f1) / cs
  IF (ABS(f) > ABS(ck(k+1))) THEN
    ck(k) = f
    f1 = f0
    f0 = f
    IF (ABS(f) > 1.0D+100) THEN
      DO  k1 = nm, k, -1
        ck(k1) = ck(k1) * 1.0D-100
      END DO
      f1 = f1 * 1.0D-100
      f0 = f0 * 1.0D-100
    END IF
  ELSE
    kb = k
    fl = ck(k+1)
    f1 = 1.0D0
    f2 = 0.25D0 * ((m+ip)*(m+ip+1.0)-cv+cs) / (m+1.0) * f1
    ck(1) = f1
    IF (kb == 1) THEN
      fs = f2
    ELSE IF (kb == 2) THEN
      ck(2) = f2
      fs = 0.125D0 * (((m+ip+2)*(m+ip+3)-cv+cs)*f2 - cs*f1) / (m+2)
    ELSE
      ck(2) = f2
      DO  j = 3, kb + 1
        f = 0.25D0 * (((2*j+m+ip-4)*(2*j+m+ip-3)-cv+cs)*f2 - cs*f1) / ((j-1)*(j+m-1))
        IF (j <= kb) ck(j) = f
        f1 = f2
        f2 = f
      END DO
      fs = f
    END IF
    EXIT
  END IF
END DO
su1 = SUM( ck(1:kb) )
su2 = SUM( ck(kb+1:nm) )
r1 = 1.0D0
DO  j = 1, (n+m+ip) / 2
  r1 = r1 * (j + 0.5D0*(n+m+ip))
END DO
r2 = 1.0D0
DO  j = 1, (n-m-ip) / 2
  r2 = -r2 * j
END DO
IF (kb == 0) THEN
  s0 = r1 / (2**n*r2*su2)
ELSE
  s0 = r1 / (2**n*r2*(fl/fs*su1 + su2))
END IF
DO  k = 1, kb
  ck(k) = fl / fs * s0 * ck(k)
END DO
DO  k = kb + 1, nm
  ck(k) = s0 * ck(k)
END DO
RETURN
END SUBROUTINE scka


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
    eg(i) = (i+m) * (i+m-1.0D0)
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
    dk1 = m + k + 1
    dk2 = 2 * (m+k)
    d2k = 2 * m + k
    a(i) = (d2k+2.0) * (d2k+1.0) / ((dk2+3.0)*(dk2+5.0)) * cs
    d(i) = dk0 * dk1 + (2.0*dk0*dk1-2.0*m*m-1.0) / ((dk2-1.0)*(dk2 +3.0)) * cs
    g(i) = k * (k-1) / ((dk2-3.0)*(dk2-1.0)) * cs
  END DO
  DO  k = 2, nm
    e(k) = SQRT(a(k-1)*g(k))
    f(k) = e(k) * e(k)
  END DO
  f(1) = 0.0D0
  e(1) = 0.0D0
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

    80 x1 = (b(k)+h(k)) / 2.0D0
    cv0(k) = x1
    IF (ABS((b(k)-h(k))/x1) >= 1.0D-14) THEN
      j = 0
      s = 1.0D0
      DO  i = 1, nm
        IF (s == 0.0D0) s = s + 1.0D-30
        t = f(i) / s
        s = d(i) - t - x1
        IF (s < 0.0D0) j = j + 1
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
 
END MODULE scka_func
 
 
 
PROGRAM mscka
USE scka_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:46

!    ============================================================
!    Purpose: This program computes the expansion coefficients
!             of the prolate and oblate spheroidal functions,
!             c2k, using subroutine SCKA
!    Input :  m  --- Mode parameter
!             n  --- Mode parameter
!             c  --- Spheroidal parameter
!             cv --- Characteristic value
!             KD --- Function code
!                    KD=1 for prolate; KD=-1 for oblate
!    Output:  CK(k) --- Expansion coefficients ck;
!                       CK(1), CK(2),... correspond to
!                       c0, c2,...
!    Example: Compute the first 13 expansion coefficients C2k for
!             KD= 1, m=2, n=3, c=3.0 and cv=14.8277782138; and
!             KD=-1, m=2, n=3, c=3.0 and cv=8.80939392077

!             Coefficients of Prolate and oblate functions

!               k         C2k(c)              C2k(-ic)
!             ---------------------------------------------
!               0     .9173213327D+01     .2489664942D+02
!               1     .4718258929D+01    -.1205287032D+02
!               2     .9841212916D+00     .2410564082D+01
!               3     .1151870224D+00    -.2735821590D+00
!               4     .8733916403D-02     .2026057157D-01
!               5     .4663888254D-03    -.1061946315D-02
!               6     .1853910398D-04     .4158091152D-04
!               7     .5708084895D-06    -.1264400411D-05
!               8     .1402786472D-07     .3074963448D-07
!               9     .2817194508D-09    -.6120579463D-09
!              10     .4712094447D-11     .1015900041D-10
!              11     .6667838485D-13    -.1427953361D-12
!              12     .8087995432D-15     .1721924955D-14
!    ============================================================

REAL (dp)  :: c, ck(200), cv, eg(200)
INTEGER    :: k, kd, m, n, nm

WRITE (*,*) 'Please KD, m, n and c '
READ (*,*) kd, m, n, c
CALL segv(m, n, c, kd, cv, eg)
WRITE (*,5100) kd, m, n, c, cv
CALL scka(m, n, c, cv, kd, ck)
WRITE (*,*)
IF (kd == 1) THEN
  WRITE (*,*) 'Coefficients of Prolate function'
  WRITE (*,*)
  WRITE (*,*) '   k            C2k(c)'
ELSE
  WRITE (*,*) 'Coefficients of Oblate function'
  WRITE (*,*)
  WRITE (*,*) '   k           C2k(-ic)'
END IF
WRITE (*,*) '----------------------------'
nm = 25 + INT((n-m)/2 + c)
DO  k = 1, nm
  WRITE (*,5000) k-1, ck(k)
END DO
STOP

5000 FORMAT ('  ', i3, '    ', g18.10)
5100 FORMAT (' KD=', i3, ',  m=', i3, ',  n=', i3, ',  c=', f5.1,  &
             ',  cv =', f18.10)
END PROGRAM mscka
