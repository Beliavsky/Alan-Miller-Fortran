MODULE airya_func
 
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


SUBROUTINE airya(x, ai, bi, ad, bd)

!    ======================================================
!    Purpose: Compute Airy functions and their derivatives
!    Input:   x  --- Argument of Airy function
!    Output:  AI --- Ai(x)
!             BI --- Bi(x)
!             AD --- Ai'(x)
!             BD --- Bi'(x)
!    Routine called:
!             AJYIK for computing Jv(x), Yv(x), Iv(x) and
!             Kv(x) with v=1/3 and 2/3
!    ======================================================

REAL (dp), INTENT(IN)      :: x
REAL (dp), INTENT(OUT)     :: ai
REAL (dp), INTENT(OUT)     :: bi
REAL (dp), INTENT(OUT)     :: ad
REAL (dp), INTENT(OUT)     :: bd

REAL (dp)  :: xa, pir, c1, c2, sr3, z, xq, vj1, vj2, vy1, vy2, vi1, vi2,  &
              vk1, vk2

xa = ABS(x)
pir = 0.318309886183891_dp
c1 = 0.355028053887817_dp
c2 = 0.258819403792807_dp
sr3 = 1.732050807568877_dp
z = xa ** 1.5 / 1.5_dp
xq = SQRT(xa)
CALL ajyik(z, vj1, vj2, vy1, vy2, vi1, vi2, vk1, vk2)
IF (x == 0.0_dp) THEN
  ai = c1
  bi = sr3 * c1
  ad = -c2
  bd = sr3 * c2
ELSE IF (x > 0.0_dp) THEN
  ai = pir * xq / sr3 * vk1
  bi = xq * (pir*vk1 + 2.0_dp/sr3*vi1)
  ad = -xa / sr3 * pir * vk2
  bd = xa * (pir*vk2 + 2.0_dp/sr3*vi2)
ELSE
  ai = 0.5_dp * xq * (vj1 - vy1/sr3)
  bi = -0.5_dp * xq * (vj1/sr3 + vy1)
  ad = 0.5_dp * xa * (vj2 + vy2/sr3)
  bd = 0.5_dp * xa * (vj2/sr3-vy2)
END IF
RETURN
END SUBROUTINE airya


SUBROUTINE ajyik(x, vj1, vj2, vy1, vy2, vi1, vi2, vk1, vk2)

!       =======================================================
!       Purpose: Compute Bessel functions Jv(x) and Yv(x),
!                and modified Bessel functions Iv(x) and
!                Kv(x), and their derivatives with v=1/3,2/3
!       Input :  x --- Argument of Jv(x),Yv(x),Iv(x) and
!                      Kv(x) ( x ò 0 )
!       Output:  VJ1 --- J1/3(x)
!                VJ2 --- J2/3(x)
!                VY1 --- Y1/3(x)
!                VY2 --- Y2/3(x)
!                VI1 --- I1/3(x)
!                VI2 --- I2/3(x)
!                VK1 --- K1/3(x)
!                VK2 --- K2/3(x)
!       =======================================================

REAL (dp), INTENT(IN)      :: x
REAL (dp), INTENT(OUT)     :: vj1
REAL (dp), INTENT(OUT)     :: vj2
REAL (dp), INTENT(OUT)     :: vy1
REAL (dp), INTENT(OUT)     :: vy2
REAL (dp), INTENT(OUT)     :: vi1
REAL (dp), INTENT(OUT)     :: vi2
REAL (dp), INTENT(OUT)     :: vk1
REAL (dp), INTENT(OUT)     :: vk2

INTEGER    :: k, k0, l
REAL (dp)  :: pi, rp2, gp1, gp2, gn1, gn2, vv0, uu0, x2, vl, vjl, r, vv, px, &
              rp, qx, rq, a0, c0, vsl, gn, sum, xk, ck, sk, uj1, uj2, vil,  &
              pv1, pv2, b0

IF (x == 0.0_dp) THEN
  vj1 = 0.0_dp
  vj2 = 0.0_dp
  vy1 = -1.0D+300
  vy2 = 1.0D+300
  vi1 = 0.0_dp
  vi2 = 0.0_dp
  vk1 = -1.0D+300
  vk2 = -1.0D+300
  RETURN
END IF
pi = 3.141592653589793_dp
rp2 = .63661977236758_dp
gp1 = .892979511569249_dp
gp2 = .902745292950934_dp
gn1 = 1.3541179394264_dp
gn2 = 2.678938534707747_dp
vv0 = 0.444444444444444_dp
uu0 = 1.1547005383793_dp
x2 = x * x
k0 = 12
IF (x >= 35.0) k0 = 10
IF (x >= 50.0) k0 = 8
IF (x <= 12.0) THEN
  DO  l = 1, 2
    vl = l / 3.0_dp
    vjl = 1.0_dp
    r = 1.0_dp
    DO  k = 1, 40
      r = -0.25_dp * r * x2 / (k*(k+vl))
      vjl = vjl + r
      IF (ABS(r) < 1.0D-15) EXIT
    END DO

    a0 = (0.5_dp*x) ** vl
    IF (l == 1) vj1 = a0 / gp1 * vjl
    IF (l == 2) vj2 = a0 / gp2 * vjl
  END DO
ELSE
  DO  l = 1, 2
    vv = vv0 * l * l
    px = 1.0_dp
    rp = 1.0_dp
    DO  k = 1, k0
      rp = -0.78125D-2 * rp * (vv - (4*k - 3)**2) * (vv-(4*k - 1)**2) / &
           (k*(2*k - 1)*x2)
      px = px + rp
    END DO
    qx = 1.0_dp
    rq = 1.0_dp
    DO  k = 1, k0
      rq = -0.78125D-2 * rq * (vv - (4*k - 1)**2) * (vv - (4*k + 1)**2) / &
            (k*(2*k + 1)*x2)
      qx = qx + rq
    END DO
    qx = 0.125_dp * (vv-1.0) * qx / x
    xk = x - (0.5_dp*l/3.0_dp + 0.25_dp) * pi
    a0 = SQRT(rp2/x)
    ck = COS(xk)
    sk = SIN(xk)
    IF (l == 1) THEN
      vj1 = a0 * (px*ck - qx*sk)
      vy1 = a0 * (px*sk + qx*ck)
    ELSE IF (l == 2) THEN
      vj2 = a0 * (px*ck-qx*sk)
      vy2 = a0 * (px*sk+qx*ck)
    END IF
  END DO
END IF
IF (x <= 12.0_dp) THEN
  DO  l = 1, 2
    vl = l / 3.0_dp
    vjl = 1.0_dp
    r = 1.0_dp
    DO  k = 1, 40
      r = -0.25_dp * r * x2 / (k*(k-vl))
      vjl = vjl + r
      IF (ABS(r) < 1.0D-15) EXIT
    END DO

    b0 = (2.0_dp/x) ** vl
    IF (l == 1) uj1 = b0 * vjl / gn1
    IF (l == 2) uj2 = b0 * vjl / gn2
  END DO
  pv1 = pi / 3.0_dp
  pv2 = pi / 1.5_dp
  vy1 = uu0 * (vj1*COS(pv1)-uj1)
  vy2 = uu0 * (vj2*COS(pv2)-uj2)
END IF
IF (x <= 18.0) THEN
  DO  l = 1, 2
    vl = l / 3.0_dp
    vil = 1.0_dp
    r = 1.0_dp
    DO  k = 1, 40
      r = 0.25_dp * r * x2 / (k*(k+vl))
      vil = vil + r
      IF (ABS(r) < 1.0D-15) EXIT
    END DO

    a0 = (0.5_dp*x) ** vl
    IF (l == 1) vi1 = a0 / gp1 * vil
    IF (l == 2) vi2 = a0 / gp2 * vil
  END DO
ELSE
  c0 = EXP(x) / SQRT(2.0_dp*pi*x)
  DO  l = 1, 2
    vv = vv0 * l * l
    vsl = 1.0_dp
    r = 1.0_dp
    DO  k = 1, k0
      r = -0.125_dp * r * (vv-(2.0_dp*k-1.0_dp)**2.0) / (k*x)
      vsl = vsl + r
    END DO
    IF (l == 1) vi1 = c0 * vsl
    IF (l == 2) vi2 = c0 * vsl
  END DO
END IF
IF (x <= 9.0_dp) THEN
  DO  l = 1, 2
    vl = l / 3.0_dp
    IF (l == 1) gn = gn1
    IF (l == 2) gn = gn2
    a0 = (2.0_dp/x) ** vl / gn
    sum = 1.0_dp
    r = 1.0_dp
    DO  k = 1, 60
      r = 0.25_dp * r * x2 / (k*(k-vl))
      sum = sum + r
      IF (ABS(r) < 1.0D-15) EXIT
    END DO

    IF (l == 1) vk1 = 0.5_dp * uu0 * pi * (sum*a0 - vi1)
    IF (l == 2) vk2 = 0.5_dp * uu0 * pi * (sum*a0 - vi2)
  END DO
ELSE
  c0 = EXP(-x) * SQRT(0.5_dp*pi/x)
  DO  l = 1, 2
    vv = vv0 * l * l
    sum = 1.0_dp
    r = 1.0_dp
    DO  k = 1, k0
      r = 0.125_dp * r * (vv-(2.0*k-1.0)**2.0) / (k*x)
      sum = sum + r
    END DO
    IF (l == 1) vk1 = c0 * sum
    IF (l == 2) vk2 = c0 * sum
  END DO
END IF
RETURN
END SUBROUTINE ajyik

END MODULE airya_func
 
 
 
PROGRAM mairya
USE airya_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:33

!       ============================================================
!       Purpose: This program computes Airy functions and their
!                derivatives using subroutine AIRYA
!       Input:   x  --- Argument of Airy function
!       Output:  AI --- Ai(x)
!                BI --- Bi(x)
!                AD --- Ai'(x)
!                BD --- Bi'(x)
!       Example:

!   x       Ai(x)          Bi(x)          Ai'(x)         Bi'(x)
!  ----------------------------------------------------------------
!   0   .35502805D+00  .61492663D+00 -.25881940D+00  .44828836D+00
!  10   .11047533D-09  .45564115D+09 -.35206337D-09  .14292361D+10
!  20   .16916729D-26  .21037650D+26 -.75863916D-26  .93818393D+26
!  30   .32082176D-48  .90572885D+47 -.17598766D-47  .49533045D+48

!   x       Ai(-x)         Bi(-x)         Ai'(-x)        Bi'(-x)
!  ----------------------------------------------------------------
!   0       .35502805      .61492663     -.25881940      .44828836
!  10       .04024124     -.31467983      .99626504      .11941411
!  20      -.17640613     -.20013931      .89286286     -.79142903
!  30      -.08796819     -.22444694     1.22862060     -.48369473
!       ============================================================

REAL (dp)  :: x, ai, bi, ad, bd

DO
  WRITE (*,*) 'Please enter x: '
  READ (*,*) x
  CALL airya(x, ai, bi, ad, bd)
  WRITE (*,5200)
  WRITE (*,5300)
  WRITE (*,5000) x, ai, bi, ad, bd
  WRITE (*,*)
  CALL airya(-x, ai, bi, ad, bd)
  WRITE (*,5400)
  WRITE (*,5300)
  WRITE (*,5100) x, ai, bi, ad, bd
END DO
STOP

5000 FORMAT (' ', f5.1, 4g16.8)
5100 FORMAT (' ', f5.1, 4g16.8)
5200 FORMAT (t5, 'x        Ai(x)           Bi(x)           Ai''(X)          BI''(X)')
5300 FORMAT (t3, '---------------------------------------------------------------------')
5400 FORMAT (t5, 'x        Ai(-x)          Bi(-x)          Ai''(-X)         Bi''(-X)')
END PROGRAM mairya
