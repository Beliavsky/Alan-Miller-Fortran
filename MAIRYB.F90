MODULE airyb_func
 
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
 

SUBROUTINE airyb(x, ai, bi, ad, bd)

!    =======================================================
!    Purpose: Compute Airy functions and their derivatives
!    Input:   x  --- Argument of Airy function
!    Output:  AI --- Ai(x)
!             BI --- Bi(x)
!             AD --- Ai'(x)
!             BD --- Bi'(x)
!    =======================================================

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: ai
REAL (dp), INTENT(OUT)  :: bi
REAL (dp), INTENT(OUT)  :: ad
REAL (dp), INTENT(OUT)  :: bd

REAL (dp)  :: ck(41), dk(41)
REAL (dp)  :: eps, pi, c1, c2, sr3, xa, xq, xm, fx, r, gx, df, dg, xe, xr1, &
              xar, xf, rp, sai, sad, sbi, sbd, xp1, xcs, xss, ssa, sda,  &
              xr2, ssb, sdb
INTEGER    :: k, km

eps = 1.0D-15
pi = 3.141592653589793_dp
c1 = 0.355028053887817_dp
c2 = 0.258819403792807_dp
sr3 = 1.732050807568877_dp
xa = ABS(x)
xq = SQRT(xa)
IF (x > 0.0_dp) xm = 5.0
IF (x <= 0.0_dp) xm = 8.0
IF (x == 0.0_dp) THEN
  ai = c1
  bi = sr3 * c1
  ad = -c2
  bd = sr3 * c2
  RETURN
END IF
IF (xa <= xm) THEN
  fx = 1.0_dp
  r = 1.0_dp
  DO  k = 1, 40
    r = r * x / (3.0_dp*k) * x / (3.0_dp*k - 1.0_dp) * x
    fx = fx + r
    IF (ABS(r) < ABS(fx)*eps) EXIT
  END DO

  gx = x
  r = x
  DO  k = 1, 40
    r = r * x / (3.0_dp*k) * x / (3.0_dp*k + 1.0_dp) * x
    gx = gx + r
    IF (ABS(r) < ABS(gx)*eps) EXIT
  END DO

  ai = c1 * fx - c2 * gx
  bi = sr3 * (c1*fx + c2*gx)
  df = 0.5_dp * x * x
  r = df
  DO  k = 1, 40
    r = r * x / (3.0_dp*k) * x / (3.0_dp*k + 2.0_dp) * x
    df = df + r
    IF (ABS(r) < ABS(df)*eps) EXIT
  END DO

  dg = 1.0_dp
  r = 1.0_dp
  DO  k = 1, 40
    r = r * x / (3.0_dp*k) * x / (3.0_dp*k - 2.0_dp) * x
    dg = dg + r
    IF (ABS(r) < ABS(dg)*eps) EXIT
  END DO

  ad = c1 * df - c2 * dg
  bd = sr3 * (c1*df + c2*dg)
ELSE
  xe = xa * xq / 1.5_dp
  xr1 = 1.0_dp / xe
  xar = 1.0_dp / xq
  xf = SQRT(xar)
  rp = 0.5641895835477563_dp
  r = 1.0_dp
  DO  k = 1, 40
    r = r * (6.0_dp*k - 1.0_dp) / 216.0_dp * (6.0_dp*k - 3.0_dp) /  &
        k * ( 6.0_dp*k - 5.0_dp) / (2.0_dp*k - 1.0_dp)
    ck(k) = r
    dk(k) = -(6.0_dp*k + 1.0_dp) / (6.0_dp*k - 1.0_dp) * ck(k)
  END DO
  km = INT(24.5-xa)
  IF (xa < 6.0) km = 14
  IF (xa > 15.0) km = 10
  IF (x > 0.0_dp) THEN
    sai = 1.0_dp
    sad = 1.0_dp
    r = 1.0_dp
    DO  k = 1, km
      r = -r * xr1
      sai = sai + ck(k) * r
      sad = sad + dk(k) * r
    END DO
    sbi = 1.0_dp
    sbd = 1.0_dp
    r = 1.0_dp
    DO  k = 1, km
      r = r * xr1
      sbi = sbi + ck(k) * r
      sbd = sbd + dk(k) * r
    END DO
    xp1 = EXP(-xe)
    ai = 0.5_dp * rp * xf * xp1 * sai
    bi = rp * xf / xp1 * sbi
    ad = -.5_dp * rp / xf * xp1 * sad
    bd = rp / xf / xp1 * sbd
  ELSE
    xcs = COS(xe + pi/4.0_dp)
    xss = SIN(xe + pi/4.0_dp)
    ssa = 1.0_dp
    sda = 1.0_dp
    r = 1.0_dp
    xr2 = 1.0_dp / (xe*xe)
    DO  k = 1, km
      r = -r * xr2
      ssa = ssa + ck(2*k) * r
      sda = sda + dk(2*k) * r
    END DO
    ssb = ck(1) * xr1
    sdb = dk(1) * xr1
    r = xr1
    DO  k = 1, km
      r = -r * xr2
      ssb = ssb + ck(2*k + 1) * r
      sdb = sdb + dk(2*k + 1) * r
    END DO
    ai = rp * xf * (xss*ssa - xcs*ssb)
    bi = rp * xf * (xcs*ssa + xss*ssb)
    ad = -rp / xf * (xcs*sda + xss*sdb)
    bd = rp / xf * (xss*sda - xcs*sdb)
  END IF
END IF
RETURN
END SUBROUTINE airyb
 
END MODULE airyb_func
 
 
 
PROGRAM mairyb
USE airyb_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:33

!       ============================================================
!       Purpose: This program computes Airy functions and their
!                derivatives using subroutine AIRYB
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
  CALL airyb(x,ai,bi,ad,bd)
  WRITE (*,5200)
  WRITE (*,5300)
  WRITE (*,5000) x, ai, bi, ad, bd
  WRITE (*,*)
  CALL airyb(-x,ai,bi,ad,bd)
  WRITE (*,5400)
  WRITE (*,5300)
  WRITE (*,5100) x, ai, bi, ad, bd
END DO
STOP

5000 FORMAT(' ', f5.1, 4g16.8)
5100 FORMAT(' ', f5.1, 4g16.8)
5200 FORMAT(t5, 'x        Ai(x)           Bi(x)           Ai''(X)          BI''(X)')
5300 FORMAT(t3, '---------------------------------------------------------------------')
5400 FORMAT(t5, 'x        Ai(-x)          Bi(-x)          Ai''(-X)         Bi''(-X)')
END PROGRAM mairyb
