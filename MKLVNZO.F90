MODULE klvnzo_func
 
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
 

SUBROUTINE klvnzo(nt, kd, zo)

!    ====================================================
!    Purpose: Compute the zeros of Kelvin functions
!    Input :  NT  --- Total number of zeros
!             KD  --- Function code
!             KD=1 to 8 for ber x, bei x, ker x, kei x,
!                       ber'x, bei'x, ker'x and kei'x,
!                       respectively.
!    Output:  ZO(M) --- the M-th zero of Kelvin function
!                       for code KD
!    Routine called:
!             KLVNA for computing Kelvin functions and
!             their derivatives
!    ====================================================

INTEGER, INTENT(IN)     :: nt
INTEGER, INTENT(IN)     :: kd
REAL (dp), INTENT(OUT)  :: zo(nt)

REAL (dp)  :: rt0(8) = (/ 2.84891_dp, 5.02622_dp, 1.71854_dp,  &
              3.91467_dp, 6.03871_dp, 3.77268_dp, 2.66584_dp, 4.93181_dp /)
REAL (dp)  :: rt, ber, bei, ger, gei, der, dei, her, hei, ddi, ddr, gdi, gdr
INTEGER    :: m

rt = rt0(kd)
DO  m = 1, nt
  10 CALL klvna(rt, ber, bei, ger, gei, der, dei, her, hei)
  IF (kd == 1) THEN
    rt = rt - ber / der
  ELSE IF (kd == 2) THEN
    rt = rt - bei / dei
  ELSE IF (kd == 3) THEN
    rt = rt - ger / her
  ELSE IF (kd == 4) THEN
    rt = rt - gei / hei
  ELSE IF (kd == 5) THEN
    ddr = -bei - der / rt
    rt = rt - der / ddr
  ELSE IF (kd == 6) THEN
    ddi = ber - dei / rt
    rt = rt - dei / ddi
  ELSE IF (kd == 7) THEN
    gdr = -gei - her / rt
    rt = rt - her / gdr
  ELSE
    gdi = ger - hei / rt
    rt = rt - hei / gdi
  END IF
  IF (ABS(rt-rt0(kd)) > 5.0D-10) THEN
    rt0(kd) = rt
    GO TO 10
  END IF
  zo(m) = rt
  rt = rt + 4.44_dp
END DO
RETURN
END SUBROUTINE klvnzo


SUBROUTINE klvna(x, ber, bei, ger, gei, der, dei, her, hei)

!    ======================================================
!    Purpose: Compute Kelvin functions ber x, bei x, ker x
!             and kei x, and their derivatives  ( x > 0 )
!    Input :  x   --- Argument of Kelvin functions
!    Output:  BER --- ber x
!             BEI --- bei x
!             GER --- ker x
!             GEI --- kei x
!             DER --- ber'x
!             DEI --- bei'x
!             HER --- ker'x
!             HEI --- kei'x
!    ================================================

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: ber
REAL (dp), INTENT(OUT)  :: bei
REAL (dp), INTENT(OUT)  :: ger
REAL (dp), INTENT(OUT)  :: gei
REAL (dp), INTENT(OUT)  :: der
REAL (dp), INTENT(OUT)  :: dei
REAL (dp), INTENT(OUT)  :: her
REAL (dp), INTENT(OUT)  :: hei

REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp, el = .5772156649015329_dp
REAL (dp)  :: cn0, cp0, cs, eps, fac, gs,  &
              pn0, pn1, pp0, pp1, qn0, qn1, qp0, qp1, r, r0, r1, rc, rs, sn0, sp0, ss,  &
              x2, x4, xc1, xc2, xd, xe1, xe2, xt
INTEGER    :: k, km, m

eps = 1.0D-15
IF (x == 0.0D0) THEN
  ber = 1.0D0
  bei = 0.0D0
  ger = 1.0D+300
  gei = -.25D0 * pi
  der = 0.0D0
  dei = 0.0D0
  her = -1.0D+300
  hei = 0.0D0
  RETURN
END IF
x2 = .25D0 * x * x
x4 = x2 * x2
IF (ABS(x) < 10.0D0) THEN
  ber = 1.0D0
  r = 1.0D0
  DO  m = 1, 60
    r = -.25D0 * r / (m*m) / (2.0D0*m-1.0D0) ** 2 * x4
    ber = ber + r
    IF (ABS(r/ber) < eps) EXIT
  END DO
  bei = x2
  r = x2
  DO  m = 1, 60
    r = -.25D0 * r / (m*m) / (2.0D0*m+1.0D0) ** 2 * x4
    bei = bei + r
    IF (ABS(r/bei) < eps) EXIT
  END DO
  ger = -(LOG(x/2.0D0)+el) * ber + .25D0 * pi * bei
  r = 1.0D0
  gs = 0.0D0
  DO  m = 1, 60
    r = -.25D0 * r / (m*m) / (2.0D0*m-1.0D0) ** 2 * x4
    gs = gs + 1.0D0 / (2.0D0*m-1.0D0) + 1.0D0 / (2.0D0*m)
    ger = ger + r * gs
    IF (ABS(r*gs/ger) < eps) EXIT
  END DO
  gei = x2 - (LOG(x/2.0D0)+el) * bei - .25D0 * pi * ber
  r = x2
  gs = 1.0D0
  DO  m = 1, 60
    r = -.25D0 * r / (m*m) / (2.0D0*m+1.0D0) ** 2 * x4
    gs = gs + 1.0D0 / (2.0D0*m) + 1.0D0 / (2.0D0*m+1.0D0)
    gei = gei + r * gs
    IF (ABS(r*gs/gei) < eps) EXIT
  END DO
  der = -.25D0 * x * x2
  r = der
  DO  m = 1, 60
    r = -.25D0 * r / m / (m+1.0D0) / (2.0D0*m+1.0D0) ** 2 * x4
    der = der + r
    IF (ABS(r/der) < eps) EXIT
  END DO
  dei = .5D0 * x
  r = dei
  DO  m = 1, 60
    r = -.25D0 * r / (m*m) / (2.d0*m-1.d0) / (2.d0*m+1.d0) * x4
    dei = dei + r
    IF (ABS(r/dei) < eps) EXIT
  END DO
  r = -.25D0 * x * x2
  gs = 1.5D0
  her = 1.5D0 * r - ber / x - (LOG(x/2.d0)+el) * der + .25 * pi * dei
  DO  m = 1, 60
    r = -.25D0 * r / m / (m+1.0D0) / (2.0D0*m+1.0D0) ** 2 * x4
    gs = gs + 1.0D0 / (2*m+1.0D0) + 1.0D0 / (2*m+2.0D0)
    her = her + r * gs
    IF (ABS(r*gs/her) < eps) EXIT
  END DO
  r = .5D0 * x
  gs = 1.0D0
  hei = .5D0 * x - bei / x - (LOG(x/2.d0)+el) * dei - .25 * pi * der
  DO  m = 1, 60
    r = -.25D0 * r / (m*m) / (2*m-1.0D0) / (2*m+1.0D0) * x4
    gs = gs + 1.0D0 / (2.0D0*m) + 1.0D0 / (2*m+1.0D0)
    hei = hei + r * gs
    IF (ABS(r*gs/hei) < eps) RETURN
  END DO
ELSE
  pp0 = 1.0D0
  pn0 = 1.0D0
  qp0 = 0.0D0
  qn0 = 0.0D0
  r0 = 1.0D0
  km = 18
  IF (ABS(x) >= 40.0) km = 10
  fac = 1.0D0
  DO  k = 1, km
    fac = -fac
    xt = .25D0 * k * pi - INT(.125D0*k) * 2.0D0 * pi
    cs = COS(xt)
    ss = SIN(xt)
    r0 = .125D0 * r0 * (2.0D0*k-1.0D0) ** 2 / k / x
    rc = r0 * cs
    rs = r0 * ss
    pp0 = pp0 + rc
    pn0 = pn0 + fac * rc
    qp0 = qp0 + rs
    qn0 = qn0 + fac * rs
  END DO
  xd = x / SQRT(2.0D0)
  xe1 = EXP(xd)
  xe2 = EXP(-xd)
  xc1 = 1.d0 / SQRT(2.0D0*pi*x)
  xc2 = SQRT(.5D0*pi/x)
  cp0 = COS(xd+.125D0*pi)
  cn0 = COS(xd-.125D0*pi)
  sp0 = SIN(xd+.125D0*pi)
  sn0 = SIN(xd-.125D0*pi)
  ger = xc2 * xe2 * (pn0*cp0 - qn0*sp0)
  gei = xc2 * xe2 * (-pn0*sp0 - qn0*cp0)
  ber = xc1 * xe1 * (pp0*cn0 + qp0*sn0) - gei / pi
  bei = xc1 * xe1 * (pp0*sn0 - qp0*cn0) + ger / pi
  pp1 = 1.0D0
  pn1 = 1.0D0
  qp1 = 0.0D0
  qn1 = 0.0D0
  r1 = 1.0D0
  fac = 1.0D0
  DO  k = 1, km
    fac = -fac
    xt = .25D0 * k * pi - INT(.125D0*k) * 2.0D0 * pi
    cs = COS(xt)
    ss = SIN(xt)
    r1 = .125D0 * r1 * (4.d0-(2.0D0*k-1.0D0)**2) / k / x
    rc = r1 * cs
    rs = r1 * ss
    pp1 = pp1 + fac * rc
    pn1 = pn1 + rc
    qp1 = qp1 + fac * rs
    qn1 = qn1 + rs
  END DO
  her = xc2 * xe2 * (-pn1*cn0 + qn1*sn0)
  hei = xc2 * xe2 * (pn1*sn0 + qn1*cn0)
  der = xc1 * xe1 * (pp1*cp0 + qp1*sp0) - hei / pi
  dei = xc1 * xe1 * (pp1*sp0 - qp1*cp0) + her / pi
END IF
RETURN
END SUBROUTINE klvna
 
END MODULE klvnzo_func
 
 
 
PROGRAM mklvnzo
USE klvnzo_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:42

!    ==============================================================
!    Purpose: This program computes the first NT zeros of Kelvin
!             functions and their derivatives using subroutine
!             KLVNZO
!    Input :  NT --- Total number of zeros
!    Example: NT = 5

!        Zeros of Kelvin functions ber x, bei x, ker x and kei x

!     m       ber x          bei x          ker x          kei x
!    ---------------------------------------------------------------
!     1     2.84891782     5.02622395     1.71854296     3.91466761
!     2     7.23882945     9.45540630     6.12727913     8.34422506
!     3    11.67396355    13.89348785    10.56294271    12.78255715
!     4    16.11356383    18.33398346    15.00268812    17.22314372
!     5    20.55463158    22.77543929    19.44381663    21.66464214

!       Zeros of Kelvin Functions ber'x, bei'x, ker'x and kei'x

!     m       ber'x          bei'x          ker'x          kei'x
!    ---------------------------------------------------------------
!     1     6.03871081     3.77267330     2.66583979     4.93181194
!     2    10.51364251     8.28098785     7.17212212     9.40405458
!     3    14.96844542    12.74214752    11.63218639    13.85826916
!     4    19.41757493    17.19343175    16.08312025    18.30717294
!     5    23.86430432    21.64114394    20.53067845    22.75379258
!    ==============================================================


REAL (dp)  :: r1(50), r2(50), r3(50), r4(50), r5(50), r6(50), r7(50), r8(50)
INTEGER    :: l, nt

WRITE (*,5300)
WRITE (*,*) 'Please enter NT '
READ (*,*) nt
WRITE (*,5100)
WRITE (*,*)
WRITE (*,*) '   m       ber x          bei x          ker x          kei x'
WRITE (*,*) ' ---------------------------------------------------------------'
CALL klvnzo(nt,1,r1)
CALL klvnzo(nt,2,r2)
CALL klvnzo(nt,3,r3)
CALL klvnzo(nt,4,r4)
DO  l = 1, nt
  WRITE (*,5000) l, r1(l), r2(l), r3(l), r4(l)
END DO
CALL klvnzo(nt,5,r5)
CALL klvnzo(nt,6,r6)
CALL klvnzo(nt,7,r7)
CALL klvnzo(nt,8,r8)
WRITE (*,*)
WRITE (*,5200)
WRITE (*,*)
WRITE (*,*) '   m       ber''X          BEI''X          KER''X          kei''X'
WRITE (*,*) ' ---------------------------------------------------------------'
DO  l = 1, nt
  WRITE (*,5000) l, r5(l), r6(l), r7(l), r8(l)
END DO
CLOSE (1)
STOP

5000 FORMAT (' ', i3, ' ', f14.8, ' ', f14.8, ' ', f14.8, ' ', f14.8)
5100 FORMAT (t5, 'Zeros of Kelvin functions ber x, bei x, ker x and kei x')
5200 FORMAT (t5, 'Zeros of Kelvin functions ber''X, BEI''X, ker''X AND KEI''X')
5300 FORMAT (' NT is the number of the zeros')
END PROGRAM mklvnzo
