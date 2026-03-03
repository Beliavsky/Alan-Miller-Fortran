MODULE airyzo_func
 
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
 

SUBROUTINE airyzo(nt, kf, xa, xb, xc, xd)

!    =========================================================
!    Purpose: Compute the first NT zeros of Airy functions
!             Ai(x) and Ai'(x), a and a', and the associated
!             values of Ai(a') and Ai'(a); and the first NT
!             zeros of Airy functions Bi(x) and Bi'(x), b and
!             b', and the associated values of Bi(b') and Bi'(b)
!    Input :  NT    --- Total number of zeros
!             KF    --- Function code
!                       KF=1 for Ai(x) and Ai'(x)
!                       KF=2 for Bi(x) and Bi'(x)
!    Output:  XA(m) --- a, the m-th zero of Ai(x) or
!                       b, the m-th zero of Bi(x)
!             XB(m) --- a', the m-th zero of Ai'(x) or
!                       b', the m-th zero of Bi'(x)
!             XC(m) --- Ai(a') or Bi(b')
!             XD(m) --- Ai'(a) or Bi'(b)
!                       ( m --- Serial number of zeros )
!    Routine called: AIRYB for computing Airy functions and
!                    their derivatives
!    ========================================================

INTEGER, INTENT(IN)     :: nt
INTEGER, INTENT(IN)     :: kf
REAL (dp), INTENT(OUT)  :: xa(nt)
REAL (dp), INTENT(OUT)  :: xb(nt)
REAL (dp), INTENT(OUT)  :: xc(nt)
REAL (dp), INTENT(OUT)  :: xd(nt)

REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp
REAL (dp)  :: u, u1, rt, rt0, ad, bd, x, ai, bi
INTEGER    :: i

DO  i = 1, nt
  IF (kf == 1) THEN
    u = 3.0 * pi * (4*i - 1) / 8.0_dp
    u1 = 1 / (u*u)
    rt0 = -(u*u) ** (1.0/3.0) * ((((-15.5902*u1 + .929844)*u1 -  &
          .138889)*u1 + .10416667_dp)*u1 + 1.0_dp)
  ELSE IF (kf == 2) THEN
    IF (i == 1) THEN
      rt0 = -1.17371
    ELSE
      u = 3.0 * pi * (4*i - 3) / 8.0
      u1 = 1.0_dp / (u*u)
      rt0 = -(u*u) ** (1.0/3.0) * ((((-15.5902*u1 + .929844)*u1 -  &
            .138889)*u1 + .10416667)*u1 + 1.0)
    END IF
  END IF

  10 x = rt0
  CALL airyb(x,ai,bi,ad,bd)
  IF (kf == 1) rt = rt0 - ai / ad
  IF (kf == 2) rt = rt0 - bi / bd
  IF (ABS((rt-rt0)/rt) > 1.d-9) THEN
    rt0 = rt
    GO TO 10
  ELSE
    xa(i) = rt
    IF (kf == 1) xd(i) = ad
    IF (kf == 2) xd(i) = bd
  END IF
END DO

DO  i = 1, nt
  IF (kf == 1) THEN
    IF (i == 1) THEN
      rt0 = -1.01879
    ELSE
      u = 3.0 * pi * (4*i - 3) / 8.0
      u1 = 1 / (u*u)
      rt0 = -(u*u) ** (1.0/3.0) * ((((15.0168*u1 - .873954)*u1 +  &
            .121528)*u1 - .145833_dp)*u1 + 1.0_dp)
    END IF
  ELSE IF (kf == 2) THEN
    IF (i == 1) THEN
      rt0 = -2.29444
    ELSE
      u = 3.0 * pi * (4*i - 1) / 8.0
      u1 = 1.0 / (u*u)
      rt0 = -(u*u) ** (1.0/3.0) * ((((15.0168*u1 -.873954)*u1 +  &
            .121528)*u1 - .145833)*u1 + 1.0)
    END IF
  END IF

  30 x = rt0
  CALL airyb(x, ai, bi, ad, bd)
  IF (kf == 1) rt = rt0 - ad / (ai*x)
  IF (kf == 2) rt = rt0 - bd / (bi*x)
  IF (ABS((rt-rt0)/rt) > 1.0D-9) THEN
    rt0 = rt
    GO TO 30
  ELSE
    xb(i) = rt
    IF (kf == 1) xc(i) = ai
    IF (kf == 2) xc(i) = bi
  END IF
END DO
RETURN
END SUBROUTINE airyzo


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
    r = r * x / (3*k) * x / (3*k - 1) * x
    fx = fx + r
    IF (ABS(r/fx) < eps) EXIT
  END DO

  gx = x
  r = x
  DO  k = 1, 40
    r = r * x / (3*k) * x / (3*k + 1) * x
    gx = gx + r
    IF (ABS(r/gx) < eps) EXIT
  END DO

  ai = c1 * fx - c2 * gx
  bi = sr3 * (c1*fx + c2*gx)
  df = .5_dp * x * x
  r = df
  DO  k = 1, 40
    r = r * x / (3*k) * x / (3*k + 2) * x
    df = df + r
    IF (ABS(r/df) < eps) EXIT
  END DO

  dg = 1.0_dp
  r = 1.0_dp
  DO  k = 1, 40
    r = r * x / (3*k) * x / (3*k - 2) * x
    dg = dg + r
    IF (ABS(r/dg) < eps) EXIT
  END DO

  ad = c1 * df - c2 * dg
  bd = sr3 * (c1*df + c2*dg)
ELSE
  xe = xa * xq / 1.5_dp
  xr1 = 1.0_dp / xe
  xar = 1.0_dp / xq
  xf = SQRT(xar)
  rp = .5641895835477563_dp
  r = 1.0_dp
  DO  k = 1, 40
    r = r * (6.0_dp*k - 1.0_dp) / 216.0_dp * (6.0_dp*k - 3.0_dp) /  &
        k * (6.0_dp*k - 5.0_dp) / (2.0_dp*k - 1.0_dp)
    ck(k) = r
    dk(k) = -(6.0_dp*k + 1.0_dp) / (6.0_dp*k - 1.0_dp) * ck(k)
  END DO
  km = INT(24.5 - xa)
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
    ai = .5_dp * rp * xf * xp1 * sai
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
 
END MODULE airyzo_func
 
 
 
PROGRAM mairyzo
USE airyzo_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:33

!    =========================================================
!    Purpose: This program computes the first NT zeros of Airy
!             functions Ai(x) and Ai'(x), and the associated
!             values of Ai(a') and Ai'(a), and the first NT
!             zeros of Airy functions Bi(x) and Bi'(x), and
!             the associated values of Bi(b') and Bi'(b) using
!             subroutine AIRYZO
!    Input :  NT    --- Total number of zeros
!             KF    --- Function code
!                       KF=1 for Ai(x) and Ai'(x)
!                       KF=2 for Bi(x) and Bi'(x)
!    Output:  XA(m) --- a, the m-th zero of Ai(x) or
!                       b, the m-th zero of Bi(x)
!             XB(m) --- a', the m-th zero of Ai'(x) or
!                       b', the m-th zero of Bi'(x)
!             XC(m) --- Ai(a') or Bi(b')
!             XD(m) --- Ai'(a) or Bi'(b)
!                       ( m --- Serial number of zeros )
!    Example: NT=5

!    m         a            Ai'(a)         a'          Ai(a')
!   -----------------------------------------------------------
!    1    -2.33810741     .70121082   -1.01879297    .53565666
!    2    -4.08794944    -.80311137   -3.24819758   -.41901548
!    3    -5.52055983     .86520403   -4.82009921    .38040647
!    4    -6.78670809    -.91085074   -6.16330736   -.35790794
!    5    -7.94413359     .94733571   -7.37217726    .34230124

!    m         b            Bi'(b)         b'          Bi(b')
!   -----------------------------------------------------------
!    1    -1.17371322     .60195789   -2.29443968   -.45494438
!    2    -3.27109330    -.76031014   -4.07315509    .39652284
!    3    -4.83073784     .83699101   -5.51239573   -.36796916
!    4    -6.16985213    -.88947990   -6.78129445    .34949912
!    5    -7.37676208     .92998364   -7.94017869   -.33602624
!    ==========================================================

REAL (dp)  :: xa(50), xb(50), xc(50), xd(50)
INTEGER    :: k, kf, nt

WRITE (*,5200)
WRITE (*,5300)
WRITE (*,*) 'Please enter KF,NT: '
READ (*,*) kf, nt
WRITE (*,5100) kf, nt
IF (kf == 1) THEN
  WRITE (*,*) '  m        a             Ai''(A)        A''           Ai(a'')'
ELSE IF (kf == 2) THEN
  WRITE (*,*) '  m        b             Bi''(B)        B''           Bi(b'')'
END IF
WRITE (*,*) '------------------------------------------------------------'
CALL airyzo(nt, kf, xa, xb, xc, xd)
DO  k = 1, nt
  WRITE (*,5000) k, xa(k), xd(k), xb(k), xc(k)
END DO
STOP

5000 FORMAT(' ', i3, ' ', 3F14.8, f13.8)
5100 FORMAT(' KF=', i2, ',     NT=', i3)
5200 FORMAT(t11, 'KF=1 for Ai(x) and Ai''(X); KF=2 FOR BI(X) and Bi''(X)')
5300 FORMAT(t11, 'NT is the number of the zeros')
END PROGRAM mairyzo
