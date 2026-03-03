MODULE ittjya_func
 
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


SUBROUTINE ittjya(x, ttj, tty)

!       =========================================================
!       Purpose: Integrate [1-J0(t)]/t with respect to t from 0
!                to x, and Y0(t)/t with respect to t from x to ì
!       Input :  x   --- Variable in the limits  ( x ò 0 )
!       Output:  TTJ --- Integration of [1-J0(t)]/t from 0 to x
!                TTY --- Integration of Y0(t)/t from x to ì
!       =========================================================

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: ttj
REAL (dp), INTENT(OUT)  :: tty

REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp, el = .5772156649015329_dp
REAL (dp)  :: a0, b1, bj0, bj1, by0, by1, e0, g0, g1, px, qx,   &
              r, r0, r1, r2, rs, t, vt, xk
INTEGER    :: k, l

IF (x == 0.0D0) THEN
  ttj = 0.0D0
  tty = -1.0D+300
ELSE IF (x <= 20.0D0) THEN
  ttj = 1.0D0
  r = 1.0D0
  DO  k = 2, 100
    r = -.25D0 * r * (k-1.0D0) / (k*k*k) * x * x
    ttj = ttj + r
    IF (ABS(r) < ABS(ttj)*1.0D-12) EXIT
  END DO
  ttj = ttj * .125D0 * x * x
  e0 = .5D0 * (pi*pi/6.0D0-el*el) - (.5D0*LOG(x/2.0D0) + el) * LOG(x/2.0D0)
  b1 = el + LOG(x/2.0D0) - 1.5D0
  rs = 1.0D0
  r = -1.0D0
  DO  k = 2, 100
    r = -.25D0 * r * (k-1.0D0) / (k*k*k) * x * x
    rs = rs + 1.0D0 / k
    r2 = r * (rs + 1.0D0/(2*k) - (el + LOG(x/2.0D0)))
    b1 = b1 + r2
    IF (ABS(r2) < ABS(b1)*1.0D-12) EXIT
  END DO
  tty = 2.0D0 / pi * (e0+.125D0*x*x*b1)
ELSE
  a0 = SQRT(2.0D0/(pi*x))
  DO  l = 0, 1
    vt = 4.0D0 * l * l
    px = 1.0D0
    r = 1.0D0
    DO  k = 1, 14
      r = -.0078125D0*r*(vt - (4*k-3)**2) / (x*k) * (vt - (4*k-1)**2) / ((2*k-1)*x)
      px = px + r
      IF (ABS(r) < ABS(px)*1.0D-12) EXIT
    END DO
    qx = 1.0D0
    r = 1.0D0
    DO  k = 1, 14
      r = -.0078125D0*r*(vt-(4*k-1)**2) / (x*k) * (vt - (4*k+1)**2) / (2*k+1) / x
      qx = qx + r
      IF (ABS(r) < ABS(qx)*1.0D-12) EXIT
    END DO
    qx = .125D0 * (vt-1.0D0) / x * qx
    xk = x - (.25D0+.5D0*l) * pi
    bj1 = a0 * (px*COS(xk) - qx*SIN(xk))
    by1 = a0 * (px*SIN(xk) + qx*COS(xk))
    IF (l == 0) THEN
      bj0 = bj1
      by0 = by1
    END IF
  END DO
  t = 2.0D0 / x
  g0 = 1.0D0
  r0 = 1.0D0
  DO  k = 1, 10
    r0 = -k * k * t * t * r0
    g0 = g0 + r0
  END DO
  g1 = 1.0D0
  r1 = 1.0D0
  DO  k = 1, 10
    r1 = -k * (k+1.0D0) * t * t * r1
    g1 = g1 + r1
  END DO
  ttj = 2.0D0 * g1 * bj0 / (x*x) - g0 * bj1 / x + el + LOG(x/2.0D0)
  tty = 2.0D0 * g1 * by0 / (x*x) - g0 * by1 / x
END IF
RETURN
END SUBROUTINE ittjya
 
END MODULE ittjya_func
 
 
 
PROGRAM mittjya
USE ittjya_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:41

!    ===========================================================
!    Purpose: This program computes the integral of [1-J0(t)]/t
!             with respect to t from 0 to x and Y0(t)/t with
!             respect to t from x to ì using subroutine ITTJYA
!    Input :  x   --- Variable in the limits  ( x ò 0 )
!    Output:  TTJ --- Integration of [1-J0(t)]/t from 0 to x
!             TTY --- Integration of Y0(t)/t from x to ì
!    Example:
!               x       [1-J0(t)]/tdt       Y0(t)/tdt
!            -------------------------------------------
!              5.0     .15403472D+01    -.46322055D-01
!             10.0     .21778664D+01    -.22987934D-01
!             15.0     .25785507D+01     .38573574D-03
!             20.0     .28773106D+01     .85031527D-02
!             25.0     .31082313D+01     .35263393D-02
!    ===========================================================

REAL (dp)  :: x, ttj, tty

WRITE (*,*) 'Please enter x '
READ (*,*) x
WRITE (*,*) '   x      [1-J0(t)]/tdt       Y0(t)/tdt'
WRITE (*,*) '-------------------------------------------'
CALL ittjya(x, ttj, tty)
WRITE (*,5000) x, ttj, tty
STOP

5000 FORMAT (' ', f5.1, 2g18.8)
END PROGRAM mittjya
