MODULE cgama_func
 
! From the book "Computation of Special Functions"
!      by Shanjie Zhang and Jianming Jin
!   Copyright 1996 by John Wiley & Sons, Inc.
! The authors state:
!   "However, we give permission to the reader who purchases this book
!    to incorporate any of these programs into his or her programs
!    provided that the copyright is acknowledged."

! Latest revision - 27 December 2001
! Corrections by Alan Miller (amiller @ bigpond.net.au)
! Variables x1 & y1 were used without values assigned to them

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
 
CONTAINS


SUBROUTINE cgama(x, y, kf, gr, gi)

!    ======================================================
!    Purpose: Compute the gamma function â(z) or ln[â(z)]
!             for a complex argument
!    Input :  x  --- Real part of z
!             y  --- Imaginary part of z
!             KF --- Function code
!                    KF=0 for ln[â(z)]
!                    KF=1 for â(z)
!    Output:  GR --- Real part of ln[â(z)] or â(z)
!             GI --- Imaginary part of ln[â(z)] or â(z)
!    ======================================================

REAL (dp), INTENT(IN OUT)  :: x
REAL (dp), INTENT(IN OUT)  :: y
INTEGER, INTENT(IN)        :: kf
REAL (dp), INTENT(OUT)     :: gr
REAL (dp), INTENT(OUT)     :: gi

REAL (dp), PARAMETER  :: a(10) = (/ 8.333333333333333D-02,  &
   -2.777777777777778D-03,  7.936507936507937D-04, -5.952380952380952D-04,  &
    8.417508417508418D-04, -1.917526917526918D-03,  6.410256410256410D-03,  &
   -2.955065359477124D-02,  1.796443723688307D-01, -1.39243221690590_dp /)
REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp
REAL (dp)  :: g0, gi1, gr1, sr, si, t, th, th1, th2, x0, x1, y1, z1, z2
INTEGER    :: j, k, na

x1 = x
y1 = y
IF (y == 0.0_dp .AND. x == INT(x) .AND. x <= 0.0_dp) THEN
  gr = 1.0D+300
  gi = 0.0_dp
  RETURN
ELSE IF (x < 0.0_dp) THEN
  x1 = x
  y1 = y
  x = -x
  y = -y
END IF
x0 = x
IF (x <= 7.0) THEN
  na = INT(7-x)
  x0 = x + na
END IF
z1 = SQRT(x0*x0 + y*y)
th = ATAN(y/x0)
gr = (x0-.5_dp) * LOG(z1) - th * y - x0 + 0.5_dp * LOG(2.0_dp*pi)
gi = th * (x0 - 0.5_dp) + y * LOG(z1) - y
DO  k = 1, 10
  t = z1 ** (1-2*k)
  gr = gr + a(k) * t * COS((2*k-1)*th)
  gi = gi - a(k) * t * SIN((2*k-1)*th)
END DO
IF (x <= 7.0) THEN
  gr1 = 0.0_dp
  gi1 = 0.0_dp
  DO  j = 0, na - 1
    gr1 = gr1 + .5_dp * LOG((x+j)**2 + y*y)
    gi1 = gi1 + ATAN(y/(x+j))
  END DO
  gr = gr - gr1
  gi = gi - gi1
END IF
IF (x1 < 0.0_dp) THEN
  z1 = SQRT(x*x + y*y)
  th1 = ATAN(y/x)
  sr = -SIN(pi*x) * COSH(pi*y)
  si = -COS(pi*x) * SINH(pi*y)
  z2 = SQRT(sr*sr + si*si)
  th2 = ATAN(si/sr)
  IF (sr < 0.0_dp) th2 = pi + th2
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
 
END MODULE cgama_func
 
 
 
PROGRAM mcgama
USE cgama_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:34

!    =============================================================
!    Purpose: This program computes the gamma function â(z) or
!             ln[â(z)] for a complex argument using subroutine CGAMA
!    Input :  x  --- Real part of z
!             y  --- Imaginary part of z
!             KF --- Function code
!                    KF=0 for ln[â(z)]
!                    KF=1 for â(z)
!    Output:  GR --- Real part of ln[â(z)] or â(z)
!             GI --- Imaginary part of ln[â(z)] or â(z)
!    Examples:

!      x         y           Re[â(z)]           Im[â(z)]
!    --------------------------------------------------------
!     2.50      5.00     .2267360319D-01    -.1172284404D-01
!     5.00     10.00     .1327696517D-01     .3639011746D-02
!     2.50     -5.00     .2267360319D-01     .1172284404D-01
!     5.00    -10.00     .1327696517D-01    -.3639011746D-02

!      x         y          Re[lnâ(z)]         Im[lnâ(z)]
!   ---------------------------------------------------------
!     2.50      5.00    -.3668103262D+01     .5806009801D+01
!     5.00     10.00    -.4285507444D+01     .1911707090D+02
!     2.50     -5.00    -.3668103262D+01    -.5806009801D+01
!     5.00    -10.00    -.4285507444D+01    -.1911707090D+02
!    =============================================================

REAL (dp)  :: x, y, gr, gi
INTEGER    :: kf

DO
  WRITE (*,*) '  Please enter KF, x and y: '
  READ (*,*) kf, x, y
  WRITE (*,*)
  IF (kf == 1) THEN
    WRITE (*,*) '       x         y           Re[â(z)]', '           Im[â(z)]'
  ELSE
    WRITE (*,*) '       x         y          Re[lnâ(z)]', '         Im[lnâ(z)]'
  END IF
  WRITE (*,*) '    ---------------------------------------------------------'
  CALL cgama(x, y, kf, gr, gi)
  WRITE (*,5000) x, y, gr, gi
END DO
STOP

5000 FORMAT (' ', 2F10.2, 2g20.10)
END PROGRAM mcgama
