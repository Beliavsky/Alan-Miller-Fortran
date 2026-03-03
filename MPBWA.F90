MODULE pbwa_func
 
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


SUBROUTINE pbwa(a, x, w1f, w1d, w2f, w2d)

!    ======================================================
!    Purpose: Compute parabolic cylinder functions W(a,ñx)
!             and their derivatives
!    Input  : a --- Parameter  ( 0 ó |a| ó 5 )
!             x --- Argument of W(a,ñx)  ( 0 ó |x| ó 5 )
!    Output : W1F --- W(a,x)
!             W1D --- W'(a,x)
!             W2F --- W(a,-x)
!             W2D --- W'(a,-x)
!    Routine called:
!            CGAMA for computing complex gamma function
!    ======================================================

REAL (dp), INTENT(IN)   :: a
REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: w1f
REAL (dp), INTENT(OUT)  :: w1d
REAL (dp), INTENT(OUT)  :: w2f
REAL (dp), INTENT(OUT)  :: w2d

REAL (dp)  :: d1, d2, dl, eps, f1, f2, g1, g2, h(100), h0, h1, hl, d(100),  &
              p0, r, r1, ugi, ugr, vgi, vgr, x1, x2, y1, y1d, y1f, y2d, y2f
INTEGER    :: k, l1, l2, m

eps = 1.0D-15
p0 = 0.59460355750136D0
IF (a == 0.0D0) THEN
  g1 = 3.625609908222D0
  g2 = 1.225416702465D0
ELSE
  x1 = 0.25D0
  y1 = 0.5D0 * a
  CALL cgama(x1, y1, 1, ugr, ugi)
  g1 = SQRT(ugr*ugr + ugi*ugi)
  x2 = 0.75D0
  CALL cgama(x2, y1, 1, vgr, vgi)
  g2 = SQRT(vgr*vgr + vgi*vgi)
END IF
f1 = SQRT(g1/g2)
f2 = SQRT(2.0D0*g2/g1)
h0 = 1.0D0
h1 = a
h(1) = a
DO  l1 = 4, 200, 2
  m = l1 / 2
  hl = a * h1 - 0.25D0 * (l1-2.0D0) * (l1-3.0D0) * h0
  h(m) = hl
  h0 = h1
  h1 = hl
END DO
y1f = 1.0D0
r = 1.0D0
DO  k = 1, 100
  r = 0.5D0 * r * x * x / (k*(2*k-1))
  r1 = h(k) * r
  y1f = y1f + r1
  IF (ABS(r1/y1f) <= eps .AND. k > 30) EXIT
END DO
y1d = a
r = 1.0D0
DO  k = 1, 100
  r = 0.5D0 * r * x * x / (k*(2.0D0*k+1.0D0))
  r1 = h(k+1) * r
  y1d = y1d + r1
  IF (ABS(r1/y1d) <= eps .AND. k > 30) EXIT
END DO
y1d = x * y1d
d1 = 1.0D0
d2 = a
d(1) = 1.0D0
d(2) = a
DO  l2 = 5, 160, 2
  m = (l2+1) / 2
  dl = a * d2 - 0.25D0 * (l2-2.0D0) * (l2-3.0D0) * d1
  d(m) = dl
  d1 = d2
  d2 = dl
END DO
y2f = 1.0D0
r = 1.0D0
DO  k = 1, 100
  r = 0.5D0 * r * x * x / (k*(2*k+1))
  r1 = d(k+1) * r
  y2f = y2f + r1
  IF (ABS(r1/y2f) <= eps .AND. k > 30) EXIT
END DO
y2f = x * y2f
y2d = 1.0D0
r = 1.0D0
DO  k = 1, 100
  r = 0.5D0 * r * x * x / (k*(2*k-1))
  r1 = d(k+1) * r
  y2d = y2d + r1
  IF (ABS(r1/y2d) <= eps .AND. k > 30) EXIT
END DO
w1f = p0 * (f1*y1f - f2*y2f)
w2f = p0 * (f1*y1f + f2*y2f)
w1d = p0 * (f1*y1d - f2*y2d)
w2d = p0 * (f1*y1d + f2*y2d)
RETURN
END SUBROUTINE pbwa


SUBROUTINE cgama(x, y, kf, gr, gi)

!    =========================================================
!    Purpose: Compute complex gamma function â(z) or Ln[â(z)]
!    Input :  x  --- Real part of z
!             y  --- Imaginary part of z
!             KF --- Function code
!                    KF=0 for Ln[â(z)]
!                    KF=1 for â(z)
!    Output:  GR --- Real part of Ln[â(z)] or â(z)
!             GI --- Imaginary part of Ln[â(z)] or â(z)
!    ========================================================

REAL (dp), INTENT(IN OUT)  :: x
REAL (dp), INTENT(IN OUT)  :: y
INTEGER, INTENT(IN)        :: kf
REAL (dp), INTENT(OUT)     :: gr
REAL (dp), INTENT(OUT)     :: gi

REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp
REAL (dp), PARAMETER  :: a(10) = (/ 8.333333333333333D-02,  &
   -2.777777777777778D-03, 7.936507936507937D-04, -5.952380952380952D-04,  &
    8.417508417508418D-04, -1.917526917526918D-03, 6.410256410256410D-03,  &
   -2.955065359477124D-02, 1.796443723688307D-01, -1.39243221690590D+00 /)
REAL (dp)  :: g0, gi1, gr1, si, sr, t, th, th1, th2, x0, x1, y1, z1, z2
INTEGER    :: j, k, na

IF (y == 0.0D0 .AND. x == INT(x) .AND. x <= 0.0D0) THEN
  gr = 1.0D+300
  gi = 0.0D0
  RETURN
ELSE IF (x < 0.0D0) THEN
  x1 = x
  y1 = y
  x = -x
  y = -y
END IF
x0 = x
IF (x <= 7.0) THEN
  na = 7-x
  x0 = x + na
END IF
z1 = SQRT(x0*x0 + y*y)
th = ATAN(y/x0)
gr = (x0-.5D0) * LOG(z1) - th * y - x0 + 0.5D0 * LOG(2.0D0*pi)
gi = th * (x0-0.5D0) + y * LOG(z1) - y
DO  k = 1, 10
  t = z1 ** (1-2*k)
  gr = gr + a(k) * t * COS((2*k-1)*th)
  gi = gi - a(k) * t * SIN((2*k-1)*th)
END DO
IF (x <= 7.0) THEN
  gr1 = 0.0D0
  gi1 = 0.0D0
  DO  j = 0, na - 1
    gr1 = gr1 + .5D0 * LOG((x+j)**2 + y*y)
    gi1 = gi1 + ATAN(y/(x+j))
  END DO
  gr = gr - gr1
  gi = gi - gi1
END IF
IF (x1 < 0.0D0) THEN
  z1 = SQRT(x*x + y*y)
  th1 = ATAN(y/x)
  sr = -SIN(pi*x) * COSH(pi*y)
  si = -COS(pi*x) * SINH(pi*y)
  z2 = SQRT(sr*sr + si*si)
  th2 = ATAN(si/sr)
  IF (sr < 0.0D0) th2 = pi + th2
  gr = LOG(pi/(z1*z2)) - gr
  gi = -th1 - th2 - gi
  x = x1
  y = y1
END IF
IF (kf == 1) THEN
  g0 = EXP(gr)
  gr = g0 * COS(gi)
  gi = g0 * SIN(gi)
END IF
RETURN
END SUBROUTINE cgama
 
END MODULE pbwa_func
 
 
 
PROGRAM mpbwa
USE pbwa_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:45

!    ===============================================================
!    Purpose: This program computes the parabolic cylinder functions
!             W(a,ñx) and their derivatives using subroutine PBWA
!    Input  : a --- Parameter  ( 0 ó |a| ó 5 )
!             x --- Argument of W(a,ñx)  ( 0 ó |x| ó 5 )
!    Output : W1F --- W(a,x)
!             W1D --- W'(a,x)
!             W2F --- W(a,-x)
!             W2D --- W'(a,-x)
!    Example: x = 5.0
!              a      W(a,x)     W'(a,x)    W(a,-x)   W'(a,-x)
!           ----------------------------------------------------
!             0.5   .1871153    .1915744  -.8556585   4.4682493
!             1.5  -.0215853    .0899870 -8.8586002  -9.3971967
!             0.0   .3009549   -.7148233   .6599634   1.7552224
!            -0.5  -.1934088  -1.3474400   .6448148   -.6781011
!            -1.5  -.5266539    .8219516  -.2822774  -1.4582283
!            -5.0   .0893618  -1.8118641   .5386084    .2698553
!    ===============================================================

REAL (dp)  :: a, w1f, w1d, w2f, w2d, x

WRITE (*,*) 'Please enter a and x '
READ (*,*) a, x
WRITE(*,5000) a, x
WRITE(*,*)
WRITE(*,*)'   a       W(a,x)          W''(A,X)         W(a,-x)         W''(A,-X)'
WRITE(*,*)' ---------------------------------------------------------------------'
CALL pbwa(a, x, w1f, w1d, w2f, w2d)
WRITE (*,5100) a, w1f, w1d, w2f, w2d
STOP

5000 FORMAT (' a=', f5.1, '   x=', f5.1)
5100 FORMAT (' ', f5.1, 4g16.8)
END PROGRAM mpbwa
