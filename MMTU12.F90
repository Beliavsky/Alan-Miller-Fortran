MODULE mtu12_func
 
! From the book "Computation of Special Functions"
!      by Shanjie Zhang and Jianming Jin
!   Copyright 1996 by John Wiley & Sons, Inc.
! The authors state:
!   "However, we give permission to the reader who purchases this book
!    to incorporate any of these programs into his or her programs
!    provided that the copyright is acknowledged."

! Zhang & Jin's errata applied to initialize array FC in routine FCOEF,
! and substantial reductions in routines REFINE & CVA2.
! Latest revision - 16 December 2002
! Corrections by Alan Miller (amiller @ bigpond.net.au)

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
 
CONTAINS
 

SUBROUTINE mtu12(kf, kc, m, q, x, f1r, d1r, f2r, d2r)

!       ==============================================================
!       Purpose: Compute modified Mathieu functions of the first and
!                second kinds, Mcm(1)(2)(x,q) and Msm(1)(2)(x,q),
!                and their derivatives
!       Input:   KF --- Function code
!                       KF=1 for computing Mcm(x,q)
!                       KF=2 for computing Msm(x,q)
!                KC --- Function Code
!                       KC=1 for computing the first kind
!                       KC=2 for computing the second kind
!                            or Msm(2)(x,q) and Msm(2)'(x,q)
!                       KC=3 for computing both the first and second kinds
!                m  --- Order of Mathieu functions
!                q  --- Parameter of Mathieu functions ( q ò 0 )
!                x  --- Argument of Mathieu functions
!       Output:  F1R --- Mcm(1)(x,q) or Msm(1)(x,q)
!                D1R --- Derivative of Mcm(1)(x,q) or Msm(1)(x,q)
!                F2R --- Mcm(2)(x,q) or Msm(2)(x,q)
!                D2R --- Derivative of Mcm(2)(x,q) or Msm(2)(x,q)
!       Routines called:
!            (1) CVA2 for computing the characteristic values
!            (2) FCOEF for computing expansion coefficients
!            (3) JYNB for computing Jn(x), Yn(x) and their derivatives
!       ==============================================================

INTEGER, INTENT(IN)     :: kf
INTEGER, INTENT(IN)     :: kc
INTEGER, INTENT(IN)     :: m
REAL (dp), INTENT(IN)   :: q
REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: f1r
REAL (dp), INTENT(OUT)  :: d1r
REAL (dp), INTENT(OUT)  :: f2r
REAL (dp), INTENT(OUT)  :: d2r

REAL (dp)  :: fg(251), bj1(0:251), dj1(0:251), bj2(0:251), dj2(0:251),  &
              by1(0:251), dy1(0:251), by2(0:251), dy2(0:251)
REAL (dp)  :: a, c1, c2, eps, qm, u1, u2, w1, w2
INTEGER    :: ic, k, kd, km, nm

eps = 1.0D-14
IF (kf == 1 .AND. m == 2*INT(m/2)) kd = 1
IF (kf == 1 .AND. m /= 2*INT(m/2)) kd = 2
IF (kf == 2 .AND. m /= 2*INT(m/2)) kd = 3
IF (kf == 2 .AND. m == 2*INT(m/2)) kd = 4
CALL cva2(kd, m, q, a)
IF (q <= 1.0_dp) THEN
  qm = 7.5 + 56.1 * SQRT(q) - 134.7 * q + 90.7 * SQRT(q) * q
ELSE
  qm = 17.0 + 3.1 * SQRT(q) - .126 * q + .0037 * SQRT(q) * q
END IF
km = INT(qm+0.5*m)
CALL fcoef(kd, m, q, a, fg)
ic = INT(m/2) + 1
IF (kd == 4) ic = m / 2
c1 = EXP(-x)
c2 = EXP(x)
u1 = SQRT(q) * c1
u2 = SQRT(q) * c2
CALL jynb(km, u1, nm, bj1, dj1, by1, dy1)
CALL jynb(km, u2, nm, bj2, dj2, by2, dy2)
IF (kc /= 2) THEN
  f1r = 0.0_dp
  DO  k = 1, km
    IF (kd == 1) THEN
      f1r = f1r + (-1) ** (ic+k) * fg(k) * bj1(k-1) * bj2(k-1)
    ELSE IF (kd == 2 .OR. kd == 3) THEN
      f1r = f1r + (-1) ** (ic+k) * fg(k) * (bj1(k-1)*bj2(k) +  &
            (-1)**kd*bj1(k)*bj2(k-1))
    ELSE
      f1r = f1r + (-1) ** (ic+k) * fg(k) * (bj1(k-1)*bj2(k+1) -  &
            bj1(k+1)*bj2(k-1))
    END IF
    IF (k >= 5 .AND. ABS(f1r-w1) < ABS(f1r)*eps) EXIT
    w1 = f1r
  END DO

  f1r = f1r / fg(1)
  d1r = 0.0_dp
  DO  k = 1, km
    IF (kd == 1) THEN
      d1r = d1r + (-1) ** (ic+k) * fg(k) * (c2*bj1(k-1)*dj2(k-1)-  &
          c1*dj1(k-1)*bj2(k-1))
    ELSE IF (kd == 2 .OR. kd == 3) THEN
      d1r = d1r + (-1) ** (ic+k) * fg(k) * (c2*(bj1(k-1)*dj2(k)+(  &
          -1)**kd*bj1(k)*dj2(k-1))-c1*(dj1(k-1)*bj2(k)+(-1)**kd* dj1(k)*bj2(k-1)))
    ELSE
      d1r = d1r + (-1)**(ic+k)*fg(k)*(c2*(bj1(k-1)*dj2(k+1) - bj1(k+1)*dj2(k-1)) - &
                  c1*(dj1(k-1)*bj2(k+1) - dj1(k+1)*bj2(k-1)))
    END IF
    IF (k >= 5 .AND. ABS(d1r-w2) < ABS(d1r)*eps) EXIT
    w2 = d1r
  END DO

  d1r = d1r * SQRT(q) / fg(1)
  IF (kc == 1) RETURN
END IF
f2r = 0.0_dp
DO  k = 1, km
  IF (kd == 1) THEN
    f2r = f2r + (-1) ** (ic+k) * fg(k) * bj1(k-1) * by2(k-1)
  ELSE IF (kd == 2 .OR. kd == 3) THEN
    f2r = f2r + (-1)**(ic+k)*fg(k)*(bj1(k-1)*by2(k) + (-1)**kd*bj1(k)*by2(k-1))
  ELSE
    f2r = f2r + (-1)**(ic+k)*fg(k)*(bj1(k-1)*by2(k+1) - bj1(k+1)*by2(k-1))
  END IF
  IF (k >= 5 .AND. ABS(f2r-w1) < ABS(f2r)*eps) EXIT
  w1 = f2r
END DO

f2r = f2r / fg(1)
d2r = 0.0_dp
DO  k = 1, km
  IF (kd == 1) THEN
    d2r = d2r + (-1) ** (ic+k) * fg(k) * (c2*bj1(k-1)*dy2(k-1) - c1*  &
        dj1(k-1)*by2(k-1))
  ELSE IF (kd == 2 .OR. kd == 3) THEN
    d2r = d2r + (-1) ** (ic+k) * fg(k) * (c2*(bj1(k-1)*dy2(k) +   &
          (-1)**kd*bj1(k)*dy2(k-1)) - c1*(dj1(k-1)*by2(k) +   &
          (-1)**kd*dj1(k)*by2(k-1)))
  ELSE
    d2r = d2r + (-1) ** (ic+k) * fg(k) * (c2*(bj1(k-1)*dy2(k+1) -  &
          bj1(k+1)*dy2(k-1)) - c1*(dj1(k-1)*by2(k+1) - dj1(k+1)*by2(k-1)))
  END IF
  IF (k >= 5 .AND. ABS(d2r-w2) < ABS(d2r)*eps) EXIT
  w2 = d2r
END DO

d2r = d2r * SQRT(q) / fg(1)
RETURN
END SUBROUTINE mtu12



SUBROUTINE fcoef(kd, m, q, a, fc)

!       =====================================================
!       Purpose: Compute expansion coefficients for Mathieu
!                functions and modified Mathieu functions
!       Input :  m  --- Order of Mathieu functions
!                q  --- Parameter of Mathieu functions
!                KD --- Case code
!                       KD=1 for cem(x,q)  ( m = 0,2,4,...)
!                       KD=2 for cem(x,q)  ( m = 1,3,5,...)
!                       KD=3 for sem(x,q)  ( m = 1,3,5,...)
!                       KD=4 for sem(x,q)  ( m = 2,4,6,...)
!                A  --- Characteristic value of Mathieu
!                       functions for given m and q
!       Output:  FC(k) --- Expansion coefficients of Mathieu
!                       functions ( k= 1,2,...,KM )
!                       FC(1),FC(2),FC(3),... correspond to
!                       A0,A2,A4,... for KD=1 case, A1,A3,
!                       A5,... for KD=2 case, B1,B3,B5,...
!                       for KD=3 case and B2,B4,B6,... for
!                       KD=4 case
!       =====================================================

INTEGER, INTENT(IN)     :: kd
INTEGER, INTENT(IN)     :: m
REAL (dp), INTENT(IN)   :: q
REAL (dp), INTENT(IN)   :: a
REAL (dp), INTENT(OUT)  :: fc(251)

REAL (dp)  :: f, f1, f2, f3, qm, s, s0, sp, ss, u, v
INTEGER    :: i, j, k, kb, km

fc = 0.0_dp
IF (q <= 1.0_dp) THEN
  qm = 7.5 + 56.1 * SQRT(q) - 134.7 * q + 90.7 * SQRT(q) * q
ELSE
  qm = 17.0 + 3.1 * SQRT(q) - .126 * q + .0037 * SQRT(q) * q
END IF
km = INT(qm+0.5*m)
IF (q == 0.0_dp) THEN
  IF (kd == 1) THEN
    fc((m+2)/2) = 1.0_dp
    IF (m == 0) fc(1) = 1.0_dp / SQRT(2.0_dp)
  ELSE IF (kd == 4) THEN
    fc(m/2) = 1.0_dp
  ELSE
    fc((m+1)/2) = 1.0_dp
  END IF
  RETURN
END IF

kb = 0
s = 0.0_dp
f = 1.0D-100
u = 0.0_dp
fc(km) = 0.0_dp
IF (kd == 1) THEN
  DO  k = km, 3, -1
    v = u
    u = f
    f = (a-4.0_dp*k*k) * u / q - v
    IF (ABS(f) < ABS(fc(k+1))) THEN
      kb = k
      fc(1) = 1.0D-100
      sp = 0.0_dp
      f3 = fc(k+1)
      fc(2) = a / q * fc(1)
      fc(3) = (a - 4.0_dp) * fc(2) / q - 2.0_dp * fc(1)
      u = fc(2)
      f1 = fc(3)
      DO  i = 3, kb
        v = u
        u = f1
        f1 = (a - 4*(i-1)**2) * u / q - v
        fc(i+1) = f1
        IF (i == kb) f2 = f1
        IF (i /= kb) sp = sp + f1 * f1
      END DO
      sp = sp + 2.0_dp * fc(1) ** 2 + fc(2) ** 2 + fc(3) ** 2
      ss = s + sp * (f3/f2) ** 2
      s0 = SQRT(1.0_dp/ss)
      DO  j = 1, km
        IF (j <= kb+1) THEN
          fc(j) = s0 * fc(j) * f3 / f2
        ELSE
          fc(j) = s0 * fc(j)
        END IF
      END DO
      GO TO 160
    ELSE
      fc(k) = f
      s = s + f * f
    END IF
  END DO
  fc(2) = q * fc(3) / (a - 4.0_dp - 2.0_dp*q*q/a)
  fc(1) = q / a * fc(2)
  s = s + 2.0_dp * fc(1) ** 2 + fc(2) ** 2
  s0 = SQRT(1.0_dp/s)
  DO  k = 1, km
    fc(k) = s0 * fc(k)
  END DO
ELSE IF (kd == 2 .OR. kd == 3) THEN
  DO  k = km, 3, -1
    v = u
    u = f
    f = (a - (2*k-1)**2) * u / q - v
    IF (ABS(f) >= ABS(fc(k))) THEN
      fc(k-1) = f
      s = s + f * f
    ELSE
      kb = k
      f3 = fc(k)
      GO TO 80
    END IF
  END DO
  fc(1) = q / (a - 1.0_dp - (-1)**kd*q) * fc(2)
  s = s + fc(1) * fc(1)
  s0 = SQRT(1.0_dp/s)
  DO  k = 1, km
    fc(k) = s0 * fc(k)
  END DO
  GO TO 160

  80 fc(1) = 1.0D-100
  fc(2) = (a - 1.0_dp - (-1)**kd*q) / q * fc(1)
  sp = 0.0_dp
  u = fc(1)
  f1 = fc(2)
  DO  i = 2, kb - 1
    v = u
    u = f1
    f1 = (a - (2*i-1)**2) * u / q - v
    IF (i /= kb-1) THEN
      fc(i+1) = f1
      sp = sp + f1 * f1
    ELSE
      f2 = f1
    END IF
  END DO
  sp = sp + fc(1) ** 2 + fc(2) ** 2
  ss = s + sp * (f3/f2) ** 2
  s0 = 1.0_dp / SQRT(ss)
  DO  j = 1, km
    IF (j < kb) fc(j) = s0 * fc(j) * f3 / f2
    IF (j >= kb) fc(j) = s0 * fc(j)
  END DO
ELSE IF (kd == 4) THEN
  DO  k = km, 3, -1
    v = u
    u = f
    f = (a - 4*k*k) * u / q - v
    IF (ABS(f) >= ABS(fc(k))) THEN
      fc(k-1) = f
      s = s + f * f
    ELSE
      kb = k
      f3 = fc(k)
      GO TO 130
    END IF
  END DO
  fc(1) = q / (a-4.0_dp) * fc(2)
  s = s + fc(1) * fc(1)
  s0 = SQRT(1.0_dp/s)
  DO  k = 1, km
    fc(k) = s0 * fc(k)
  END DO
  GO TO 160

  130 fc(1) = 1.0D-100
  fc(2) = (a-4.0_dp) / q * fc(1)
  sp = 0.0_dp
  u = fc(1)
  f1 = fc(2)
  DO  i = 2, kb - 1
    v = u
    u = f1
    f1 = (a - 4*i*i) * u / q - v
    IF (i /= kb-1) THEN
      fc(i+1) = f1
      sp = sp + f1 * f1
    ELSE
      f2 = f1
    END IF
  END DO
  sp = sp + fc(1) ** 2 + fc(2) ** 2
  ss = s + sp * (f3/f2) ** 2
  s0 = 1.0_dp / SQRT(ss)
  DO  j = 1, km
    IF (j < kb) fc(j) = s0 * fc(j) * f3 / f2
    IF (j >= kb) fc(j) = s0 * fc(j)
  END DO
END IF

160 IF (fc(1) < 0.0_dp) THEN
  DO  j = 1, km
    fc(j) = -fc(j)
  END DO
END IF
RETURN
END SUBROUTINE fcoef



SUBROUTINE cva2(kd, m, q, a)

!       ======================================================
!       Purpose: Calculate a specific characteristic value of
!                Mathieu functions
!       Input :  m  --- Order of Mathieu functions
!                q  --- Parameter of Mathieu functions
!                KD --- Case code
!                       KD=1 for cem(x,q)  ( m = 0,2,4,...)
!                       KD=2 for cem(x,q)  ( m = 1,3,5,...)
!                       KD=3 for sem(x,q)  ( m = 1,3,5,...)
!                       KD=4 for sem(x,q)  ( m = 2,4,6,...)
!       Output:  A  --- Characteristic value
!       Routines called:
!             (1) REFINE for finding accurate characteristic
!                 value using an iteration method
!             (2) CV0 for finding initial characteristic
!                 values using polynomial approximation
!             (3) CVQM for computing initial characteristic
!                 values for q ó 3*m
!             (3) CVQL for computing initial characteristic
!                 values for q ò m*m
!       ======================================================

INTEGER, INTENT(IN)     :: kd
INTEGER, INTENT(IN)     :: m
REAL (dp), INTENT(IN)   :: q
REAL (dp), INTENT(OUT)  :: a

REAL (dp)  :: a1, a2, delta, q1, q2, qq
INTEGER    :: i, iflag, ndiv, nn

IF (m <= 12 .OR. q <= 3.0*m.OR.q > m*m) THEN
  CALL cv0(kd, m, q, a)
  iflag = 1
  IF (q /= 0.0_dp) CALL refine(kd, m, q, a)
ELSE
  ndiv = 10
  delta = (m-3.0) * m / ndiv
  IF ((q-3.0*m) <= (m*m-q)) THEN
    10 nn = INT((q-3.0*m)/delta) + 1
    delta = (q-3.0*m) / nn
    q1 = 2.0 * m
    CALL cvqm(m, q1, a1)
    q2 = 3.0 * m
    CALL cvqm(m, q2, a2)
    qq = 3.0 * m
    DO  i = 1, nn
      qq = qq + delta
      a = (a1*q2 - a2*q1 + (a2-a1)*qq) / (q2-q1)
      iflag = 1
      IF (i == nn) iflag = -1
      CALL refine(kd, m, qq, a)
      q1 = q2
      q2 = qq
      a1 = a2
      a2 = a
    END DO
    IF (iflag == -10) THEN
      ndiv = ndiv * 2
      delta = (m-3.0) * m / ndiv
      GO TO 10
    END IF
  ELSE
    30 nn = INT((m*m-q)/delta) + 1
    delta = (m*m-q) / nn
    q1 = m * (m-1.0)
    CALL cvql(kd, m, q1, a1)
    q2 = m * m
    CALL cvql(kd, m, q2, a2)
    qq = m * m
    DO  i = 1, nn
      qq = qq - delta
      a = (a1*q2 - a2*q1 + (a2-a1)*qq) / (q2-q1)
      iflag = 1
      IF (i == nn) iflag = -1
      CALL refine(kd, m, qq, a)
      q1 = q2
      q2 = qq
      a1 = a2
      a2 = a
    END DO
    IF (iflag == -10) THEN
      ndiv = ndiv * 2
      delta = (m-3.0) * m / ndiv
      GO TO 30
    END IF
  END IF
END IF
RETURN
END SUBROUTINE cva2




SUBROUTINE refine(kd, m, q, a)

!       =====================================================
!       Purpose: calculate the accurate characteristic value
!                by the secant method
!       Input :  m --- Order of Mathieu functions
!                q --- Parameter of Mathieu functions
!                A --- Initial characteristic value
!       Output:  A --- Refineed characteristic value
!       Routine called:  CVF for computing the value of F for
!                        characteristic equation
! N.B. Argument IFLAG has been removed.
!       ========================================================

INTEGER, INTENT(IN)        :: kd
INTEGER, INTENT(IN)        :: m
REAL (dp), INTENT(IN)      :: q
REAL (dp), INTENT(IN OUT)  :: a

REAL (dp)  :: eps, f, f0, f1, x, x0, x1
INTEGER    :: it, mj

eps = 1.0D-14
mj = 10 + m
x0 = a
CALL cvf(kd, m, q, x0, mj, f0)
x1 = 1.002 * a
CALL cvf(kd, m, q, x1, mj, f1)
DO  it = 1, 100
  mj = mj + 1
  x = x1 - (x1-x0) / (1.0_dp - f0/f1)
  CALL cvf(kd, m, q, x, mj, f)
  IF (ABS(1.0-x1/x) < eps .OR. f == 0.0) EXIT
  x0 = x1
  f0 = f1
  x1 = x
  f1 = f
END DO

a = x
RETURN
END SUBROUTINE refine



SUBROUTINE cvf(kd, m, q, a, mj, f)

!       ======================================================
!       Purpose: Compute the value of F for characteristic
!                equation of Mathieu functions
!       Input :  m --- Order of Mathieu functions
!                q --- Parameter of Mathieu functions
!                A --- Characteristic value
!       Output:  F --- Value of F for characteristic equation
!       ======================================================

INTEGER, INTENT(IN)        :: kd
INTEGER, INTENT(IN)        :: m
REAL (dp), INTENT(IN)      :: q
REAL (dp), INTENT(IN)      :: a
INTEGER, INTENT(IN)        :: mj
REAL (dp), INTENT(OUT)     :: f

REAL (dp)  :: b, t0, t1, t2
INTEGER    :: ic, j, j0, jf, l, l0

b = a
ic = INT(m/2)
l = 0
l0 = 0
j0 = 2
jf = ic
IF (kd == 1) l0 = 2
IF (kd == 1) j0 = 3
IF (kd == 2 .OR. kd == 3) l = 1
IF (kd == 4) jf = ic - 1
t1 = 0.0_dp
DO  j = mj, ic + 1, -1
  t1 = -q * q / ((2.0_dp*j+l)**2-b+t1)
END DO
IF (m <= 2) THEN
  t2 = 0.0_dp
  IF (kd == 1 .AND. m == 0) t1 = t1 + t1
  IF (kd == 1 .AND. m == 2) t1 = -2.0 * q * q / (4.0-b+t1) - 4.0
  IF (kd == 2 .AND. m == 1) t1 = t1 + q
  IF (kd == 3 .AND. m == 1) t1 = t1 - q
ELSE
  IF (kd == 1) t0 = 4.0_dp - b + 2.0_dp * q * q / b
  IF (kd == 2) t0 = 1.0_dp - b + q
  IF (kd == 3) t0 = 1.0_dp - b - q
  IF (kd == 4) t0 = 4.0_dp - b
  t2 = -q * q / t0
  DO  j = j0, jf
    t2 = -q * q / ((2.0_dp*j-l-l0)**2-b+t2)
  END DO
END IF
f = (2.0_dp*ic+l) ** 2 + t1 + t2 - b
RETURN
END SUBROUTINE cvf



SUBROUTINE cv0(kd, m, q, a0)

!       =====================================================
!       Purpose: Compute the initial characteristic value of
!                Mathieu functions for m ó 12  or q ó 300 or
!                q ò m*m
!       Input :  m  --- Order of Mathieu functions
!                q  --- Parameter of Mathieu functions
!       Output:  A0 --- Characteristic value
!       Routines called:
!             (1) CVQM for computing initial characteristic
!                 value for q ó 3*m
!             (2) CVQL for computing initial characteristic
!                 value for q ò m*m
!       ====================================================

INTEGER, INTENT(IN)     :: kd
INTEGER, INTENT(IN)     :: m
REAL (dp), INTENT(IN)   :: q
REAL (dp), INTENT(OUT)  :: a0

REAL (dp)  :: q2

q2 = q * q
IF (m == 0) THEN
  IF (q <= 1.0) THEN
    a0 = (((.0036392*q2-.0125868)*q2+.0546875)*q2-.5) * q2
  ELSE IF (q <= 10.0) THEN
    a0 = ((3.999267D-3*q-9.638957D-2)*q-.88297) * q + .5542818
  ELSE
    CALL cvql(kd, m, q, a0)
  END IF
ELSE IF (m == 1) THEN
  IF (q <= 1.0 .AND. kd == 2) THEN
    a0 = (((-6.51E-4*q-.015625)*q-.125)*q+1.0) * q + 1.0
  ELSE IF (q <= 1.0 .AND. kd == 3) THEN
    a0 = (((-6.51E-4*q+.015625)*q-.125)*q-1.0) * q + 1.0
  ELSE IF (q <= 10.0 .AND. kd == 2) THEN
    a0 = (((-4.94603D-4*q+1.92917D-2)*q-.3089229)*q+1.33372) * q + .811752
  ELSE IF (q <= 10.0 .AND. kd == 3) THEN
    a0 = ((1.971096D-3*q-5.482465D-2)*q-1.152218) * q + 1.10427
  ELSE
    CALL cvql(kd, m, q, a0)
  END IF
ELSE IF (m == 2) THEN
  IF (q <= 1.0 .AND. kd == 1) THEN
    a0 = (((-.0036391*q2+.0125888)*q2-.0551939)*q2+.416667) * q2 + 4.0
  ELSE IF (q <= 1.0 .AND. kd == 4) THEN
    a0 = (.0003617*q2-.0833333) * q2 + 4.0
  ELSE IF (q <= 15 .AND. kd == 1) THEN
    a0 = (((3.200972D-4*q-8.667445D-3)*q-1.829032D-4)*q+.9919999)  &
        * q + 3.3290504
  ELSE IF (q <= 10.0 .AND. kd == 4) THEN
    a0 = ((2.38446D-3*q-.08725329)*q-4.732542D-3) * q + 4.00909
  ELSE
    CALL cvql(kd, m, q, a0)
  END IF
ELSE IF (m == 3) THEN
  IF (q <= 1.0 .AND. kd == 2) THEN
    a0 = ((6.348E-4*q+.015625)*q+.0625) * q2 + 9.0
  ELSE IF (q <= 1.0 .AND. kd == 3) THEN
    a0 = ((6.348E-4*q-.015625)*q+.0625) * q2 + 9.0
  ELSE IF (q <= 20.0 .AND. kd == 2) THEN
    a0 = (((3.035731D-4*q-1.453021D-2)*q+.19069602)*q-.1039356) *  &
        q + 8.9449274
  ELSE IF (q <= 15.0 .AND. kd == 3) THEN
    a0 = ((9.369364D-5*q-.03569325)*q+.2689874) * q + 8.771735
  ELSE
    CALL cvql(kd, m, q, a0)
  END IF
ELSE IF (m == 4) THEN
  IF (q <= 1.0 .AND. kd == 1) THEN
    a0 = ((-2.1E-6*q2+5.012E-4)*q2+.0333333) * q2 + 16.0
  ELSE IF (q <= 1.0 .AND. kd == 4) THEN
    a0 = ((3.7E-6*q2-3.669E-4)*q2+.0333333) * q2 + 16.0
  ELSE IF (q <= 25.0 .AND. kd == 1) THEN
    a0 = (((1.076676D-4*q-7.9684875D-3)*q+.17344854)*q-.5924058) *  &
        q + 16.620847
  ELSE IF (q <= 20.0 .AND. kd == 4) THEN
    a0 = ((-7.08719D-4*q+3.8216144D-3)*q+.1907493) * q + 15.744
  ELSE
    CALL cvql(kd, m, q, a0)
  END IF
ELSE IF (m == 5) THEN
  IF (q <= 1.0 .AND. kd == 2) THEN
    a0 = ((6.8E-6*q+1.42E-5)*q2+.0208333) * q2 + 25.0
  ELSE IF (q <= 1.0 .AND. kd == 3) THEN
    a0 = ((-6.8E-6*q+1.42E-5)*q2+.0208333) * q2 + 25.0
  ELSE IF (q <= 35.0 .AND. kd == 2) THEN
    a0 = (((2.238231D-5*q-2.983416D-3)*q+.10706975)*q-.600205) * q + 25.93515
  ELSE IF (q <= 25.0 .AND. kd == 3) THEN
    a0 = ((-7.425364D-4*q+2.18225D-2)*q+4.16399D-2) * q + 24.897
  ELSE
    CALL cvql(kd, m, q, a0)
  END IF
ELSE IF (m == 6) THEN
  IF (q <= 1.0) THEN
    a0 = (.4D-6*q2+.0142857) * q2 + 36.0
  ELSE IF (q <= 40.0 .AND. kd == 1) THEN
    a0 = (((-1.66846D-5*q+4.80263D-4)*q+2.53998D-2)*q-.181233) * q + 36.423
  ELSE IF (q <= 35.0 .AND. kd == 4) THEN
    a0 = ((-4.57146D-4*q+2.16609D-2)*q-2.349616D-2) * q + 35.99251
  ELSE
    CALL cvql(kd, m, q, a0)
  END IF
ELSE IF (m == 7) THEN
  IF (q <= 10.0) THEN
    CALL cvqm(m, q, a0)
  ELSE IF (q <= 50.0 .AND. kd == 2) THEN
    a0 = (((-1.411114D-5*q+9.730514D-4)*q-3.097887D-3)*q+  &
        3.533597D-2) * q + 49.0547
  ELSE IF (q <= 40.0 .AND. kd == 3) THEN
    a0 = ((-3.043872D-4*q+2.05511D-2)*q-9.16292D-2) * q + 49.19035
  ELSE
    CALL cvql(kd, m, q, a0)
  END IF
ELSE IF (m >= 8) THEN
  IF (q <= 3.*m) THEN
    CALL cvqm(m, q, a0)
  ELSE IF (q > m*m) THEN
    CALL cvql(kd, m, q, a0)
  ELSE
    IF (m == 8 .AND. kd == 1) THEN
      a0 = (((8.634308D-6*q-2.100289D-3)*q+.169072)*q-4.64336) * q + 109.4211
    ELSE IF (m == 8 .AND. kd == 4) THEN
      a0 = ((-6.7842D-5*q+2.2057D-3)*q+.48296) * q + 56.59
    ELSE IF (m == 9 .AND. kd == 2) THEN
      a0 = (((2.906435D-6*q-1.019893D-3)*q+.1101965)*q-3.821851) *  &
          q + 127.6098
    ELSE IF (m == 9 .AND. kd == 3) THEN
      a0 = ((-9.577289D-5*q+.01043839)*q+.06588934) * q + 78.0198
    ELSE IF (m == 10 .AND. kd == 1) THEN
      a0 = (((5.44927D-7*q-3.926119D-4)*q+.0612099)*q-2.600805) * q + 138.1923
    ELSE IF (m == 10 .AND. kd == 4) THEN
      a0 = ((-7.660143D-5*q+.01132506)*q-.09746023) * q + 99.29494
    ELSE IF (m == 11 .AND. kd == 2) THEN
      a0 = (((-5.67615D-7*q+7.152722D-6)*q+.01920291)*q-1.081583) * q + 140.88
    ELSE IF (m == 11 .AND. kd == 3) THEN
      a0 = ((-6.310551D-5*q+.0119247)*q-.2681195) * q + 123.667
    ELSE IF (m == 12 .AND. kd == 1) THEN
      a0 = (((-2.38351D-7*q-2.90139D-5)*q+.02023088)*q-1.289) * q + 171.2723
    ELSE IF (m == 12 .AND. kd == 4) THEN
      a0 = (((3.08902D-7*q-1.577869D-4)*q+.0247911)*q-1.05454) * q + 161.471
    END IF
  END IF
END IF
RETURN
END SUBROUTINE cv0



SUBROUTINE cvql(kd, m, q, a0)

!       ========================================================
!       Purpose: Compute the characteristic value of Mathieu
!                functions  for q ò 3m
!       Input :  m  --- Order of Mathieu functions
!                q  --- Parameter of Mathieu functions
!       Output:  A0 --- Initial characteristic value
!       ========================================================

INTEGER, INTENT(IN)     :: kd
INTEGER, INTENT(IN)     :: m
REAL (dp), INTENT(IN)   :: q
REAL (dp), INTENT(OUT)  :: a0

REAL (dp)  :: c1, cv1, cv2, d1, d2, d3, d4, p1, p2, w, w2, w3, w4, w6

IF (kd == 1 .OR. kd == 2) w = 2 * m + 1
IF (kd == 3 .OR. kd == 4) w = 2 * m - 1
w2 = w * w
w3 = w * w2
w4 = w2 * w2
w6 = w2 * w4
d1 = 5.0 + 34.0 / w2 + 9.0 / w4
d2 = (33.0 + 410.0/w2 + 405.0/w4) / w
d3 = (63.0 + 1260.0/w2 + 2943.0/w4 + 486.0/w6) / w2
d4 = (527.0 + 15617.0/w2 + 69001.0/w4 + 41607.0/w6) / w3
c1 = 128.0
p2 = q / w4
p1 = SQRT(p2)
cv1 = -2.0 * q + 2.0 * w * SQRT(q) - (w2+1.0) / 8.0
cv2 = (w+3.0/w) + d1 / (32.0*p1) + d2 / (8.0*c1*p2)
cv2 = cv2 + d3 / (64.0*c1*p1*p2) + d4 / (16.0*c1*c1*p2*p2)
a0 = cv1 - cv2 / (c1*p1)
RETURN
END SUBROUTINE cvql



SUBROUTINE cvqm(m, q, a0)

!       =====================================================
!       Purpose: Compute the characteristic value of Mathieu
!                functions for q ó m*m
!       Input :  m  --- Order of Mathieu functions
!                q  --- Parameter of Mathieu functions
!       Output:  A0 --- Initial characteristic value
!       =====================================================

INTEGER, INTENT(IN)     :: m
REAL (dp), INTENT(IN)   :: q
REAL (dp), INTENT(OUT)  :: a0

REAL (dp)  :: hm1, hm3, hm5
INTEGER    :: m2

m2 = m * m
hm1 = .5 * q / (m2-1)
hm3 = .25 * hm1 ** 3 / (m2-4)
hm5 = hm1 * hm3 * q / ((m2-1)*(m2-9))
a0 = m2 + q * (hm1 + (5.0*m2 + 7.0)*hm3 + (9.0*m**4 + 58.0*m2 + 29.0)*hm5)
RETURN
END SUBROUTINE cvqm



SUBROUTINE jynb(n, x, nm, bj, dj, by, dy)

!       =====================================================
!       Purpose: Compute Bessel functions Jn(x), Yn(x) and
!                their derivatives
!       Input :  x --- Argument of Jn(x) and Yn(x) ( x ò 0 )
!                n --- Order of Jn(x) and Yn(x)
!       Output:  BJ(n) --- Jn(x)
!                DJ(n) --- Jn'(x)
!                BY(n) --- Yn(x)
!                DY(n) --- Yn'(x)
!                NM --- Highest order computed
!       Routines called:
!                MSTA1 and MSTA2 to calculate the starting
!                point for backward recurrence
!       =====================================================

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: x
INTEGER, INTENT(OUT)       :: nm
REAL (dp), INTENT(OUT)     :: bj(0:n)
REAL (dp), INTENT(OUT)     :: dj(0:n)
REAL (dp), INTENT(OUT)     :: by(0:n)
REAL (dp), INTENT(OUT)     :: dy(0:n)

REAL (dp), PARAMETER :: a(4) = (/ -.7031250000000000D-01, .1121520996093750D+00,  &
      -.5725014209747314D+00, .6074042001273483D+01 /)
REAL (dp), PARAMETER :: b(4) = (/ .7324218750000000D-01, -.2271080017089844D+00,  &
      .1727727502584457D+01, -.2438052969955606D+02 /)
REAL (dp), PARAMETER :: a1(4) = (/ .1171875000000000D+00, -.1441955566406250D+00,  &
      .6765925884246826D+00, -.6883914268109947D+01 /)
REAL (dp), PARAMETER :: b1(4) = (/ -.1025390625000000D+00, .2775764465332031D+00,  &
      -.1993531733751297D+01, .2724882731126854D+02 /)
REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp, r2p = .63661977236758_dp
REAL (dp)  :: bj0, bj1, bjk, bs, by0, by1, byk, cu, ec, f, f1, f2,   &
              p0, p1, q0, q1, s0, su, sv, t1, t2
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
IF (x <= 300.0 .OR. n > INT(0.9*x)) THEN
  IF (n == 0) nm = 1
  m = msta1(x,200)
  IF (m < nm) THEN
    nm = m
  ELSE
    m = msta2(x,nm,15)
  END IF
  bs = 0.0_dp
  su = 0.0_dp
  sv = 0.0_dp
  f2 = 0.0_dp
  f1 = 1.0D-100
  DO  k = m, 0, -1
    f = 2.0_dp * (k+1.0_dp) / x * f1 - f2
    IF (k <= nm) bj(k) = f
    IF (k == 2*INT(k/2) .AND. k /= 0) THEN
      bs = bs + 2.0_dp * f
      su = su + (-1) ** (k/2) * f / k
    ELSE IF (k > 1) THEN
      sv = sv + (-1) ** (k/2) * k / (k*k-1.0) * f
    END IF
    f2 = f1
    f1 = f
  END DO
  s0 = bs + f
  DO  k = 0, nm
    bj(k) = bj(k) / s0
  END DO
  ec = LOG(x/2.0_dp) + 0.5772156649015329_dp
  by0 = r2p * (ec*bj(0)-4.0_dp*su/s0)
  by(0) = by0
  by1 = r2p * ((ec-1.0_dp)*bj(1)-bj(0)/x-4.0_dp*sv/s0)
  by(1) = by1
ELSE
  t1 = x - 0.25_dp * pi
  p0 = 1.0_dp
  q0 = -0.125_dp / x
  DO  k = 1, 4
    p0 = p0 + a(k) * x ** (-2*k)
    q0 = q0 + b(k) * x ** (-2*k-1)
  END DO
  cu = SQRT(r2p/x)
  bj0 = cu * (p0*COS(t1)-q0*SIN(t1))
  by0 = cu * (p0*SIN(t1)+q0*COS(t1))
  bj(0) = bj0
  by(0) = by0
  t2 = x - 0.75_dp * pi
  p1 = 1.0_dp
  q1 = 0.375_dp / x
  DO  k = 1, 4
    p1 = p1 + a1(k) * x ** (-2*k)
    q1 = q1 + b1(k) * x ** (-2*k-1)
  END DO
  bj1 = cu * (p1*COS(t2)-q1*SIN(t2))
  by1 = cu * (p1*SIN(t2)+q1*COS(t2))
  bj(1) = bj1
  by(1) = by1
  DO  k = 2, nm
    bjk = 2.0_dp * (k-1.0_dp) / x * bj1 - bj0
    bj(k) = bjk
    bj0 = bj1
    bj1 = bjk
  END DO
END IF
dj(0) = -bj(1)
DO  k = 1, nm
  dj(k) = bj(k-1) - k / x * bj(k)
END DO
DO  k = 2, nm
  byk = 2.0_dp * (k-1.0_dp) * by1 / x - by0
  by(k) = byk
  by0 = by1
  by1 = byk
END DO
dy(0) = -by(1)
DO  k = 1, nm
  dy(k) = by(k-1) - k * by(k) / x
END DO
RETURN
END SUBROUTINE jynb



FUNCTION msta1(x, mp) RESULT(fn_val)

!       ===================================================
!       Purpose: Determine the starting point for backward
!                recurrence such that the magnitude of
!                Jn(x) at that point is about 10^(-MP)
!       Input :  x     --- Argument of Jn(x)
!                MP    --- Value of magnitude
!       Output:  MSTA1 --- Starting point
!       ===================================================

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

!       ===================================================
!       Purpose: Determine the starting point for backward
!                recurrence such that all Jn(x) has MP
!                significant digits
!       Input :  x  --- Argument of Jn(x)
!                n  --- Order of Jn(x)
!                MP --- Significant digit
!       Output:  MSTA2 --- Starting point
!       ===================================================

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

END MODULE mtu12_func
 
 
 
PROGRAM mmtu12
USE mtu12_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:44

!       ===============================================================
!       Purpose: This program computes the modified Mathieu functions
!                of the first and second kinds, Mcm(1)(2)(x,q) and
!                Msm(1)(2)(x,q), and their derivatives using
!                subroutine MTU12
!       Input:   KF --- Function code
!                       KF=1 for computing Mcm(x,q)
!                       KF=2 for computing Msm(x,q)
!                KC --- Function Code
!                       KC=1 for computing Mcm(1)(x,q) and Mcm(1)'(x,q)
!                            or Msm(1)(x,q) and Msm(1)'(x,q)
!                       KC=2 for computing Mcm(2)(x,q) and Mcm(2)'(x,q)
!                            or Msm(2)(x,q) and Msm(2)'(x,q)
!                       KC=3 for both modified Mathieu functions of the
!                            first and second kinds, and their
!                            derivatives
!                m  --- Order of Mathieu functions
!                q  --- Parameter of Mathieu functions
!                x  --- Argument of Mathieu functions
!       Output:  F1R --- Mcm(1)(x,q) or Msm(1)(x,q)
!                D1R --- Derivative of Mcm(1)(x,q) or Msm(1)(x,q)
!                F2R --- Mcm(2)(x,q) or Msm(2)(x,q)
!                D2R --- Derivative of Mcm(2)(x,q) or Msm(2)(x,q)
!       ===============================================================

REAL (dp)  :: q, x, f1r, d1r, f2r, d2r
INTEGER    :: kc, kf, m

WRITE (*,*) 'Please enter KF, m, q and x '
READ (*,*) kf, m, q, x
WRITE (*,5000) kf, m, q, x
kc = 3
CALL mtu12(kf, kc, m, q, x, f1r, d1r, f2r, d2r)
WRITE (*,*)
IF (kf == 1) THEN
  WRITE (*,*) '   x      Mcm(1)(x,q)    Mcm(1)''(X,Q)',  &
      '    Mcm(2)(x,q)     Mcm(2)''(X,Q)'
ELSE
  WRITE (*,*) '   x      Msm(1)(x,q)    Msm(1)''(X,Q)',  &
      '    Msm(2)(x,q)     Msm(2)''(X,Q)'
END IF
WRITE (*,*) ' ---------------------------------------------------------------------'
WRITE (*,5100) x, f1r, d1r, f2r, d2r
WRITE (*,*)
WRITE (*,5200) f1r * d2r - f2r * d1r, .63661977236758_dp
WRITE (*,5300)
STOP

5000 FORMAT (' KF =', i2, ',  m =', i3, ',  q =', f5.1, ',  x =', f5.1)
5100 FORMAT (' ', f5.1, 4g16.8)
5200 FORMAT (' WRONSKIAN=', e16.8, '   should equal   2/PI=', e16.8)
5300 FORMAT (/' Caution: This check is not accurate if it involves'/  &
             '          the subtraction of two similar numbers')
END PROGRAM mmtu12
