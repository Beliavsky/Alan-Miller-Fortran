MODULE cisia_func
 
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
 


SUBROUTINE cisia(x,ci,si)

!       =============================================
!       Purpose: Compute cosine and sine integrals
!                Si(x) and Ci(x)  ( x ò 0 )
!       Input :  x  --- Argument of Ci(x) and Si(x)
!       Output:  CI --- Ci(x)
!                SI --- Si(x)
!       =============================================

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: ci
REAL (dp), INTENT(OUT)  :: si

REAL (dp)  :: bj(101), eps, x2, xa, xa0, xa1, xcs, xf, xg, xg1, xg2, xr, xs, xss
INTEGER    :: k, m
REAL (dp), PARAMETER  :: p2 = 1.570796326794897_dp, el = .5772156649015329_dp

eps = 1.0D-15
x2 = x * x
IF (x == 0.0_dp) THEN
  ci = -1.0D+300
  si = 0.0_dp
ELSE IF (x <= 16.0_dp) THEN
  xr = -.25_dp * x2
  ci = el + LOG(x) + xr
  DO  k = 2, 40
    xr = -.5_dp * xr * (k-1) / (k*k*(2*k-1)) * x2
    ci = ci + xr
    IF (ABS(xr) < ABS(ci)*eps) EXIT
  END DO
  xr = x
  si = x
  DO  k = 1, 40
    xr = -.5_dp * xr * (2*k-1) / k / (4*k*k+4*k+1) * x2
    si = si + xr
    IF (ABS(xr) < ABS(si)*eps) RETURN
  END DO
ELSE IF (x <= 32.0_dp) THEN
  m = INT(47.2 + .82*x)
  xa1 = 0.0_dp
  xa0 = 1.0D-100
  DO  k = m, 1, -1
    xa = 4.0_dp * k * xa0 / x - xa1
    bj(k) = xa
    xa1 = xa0
    xa0 = xa
  END DO
  xs = bj(1)
  DO  k = 3, m, 2
    xs = xs + 2.0_dp * bj(k)
  END DO
  bj(1) = bj(1) / xs
  bj(2:m) = bj(2:m) / xs
  xr = 1.0_dp
  xg1 = bj(1)
  DO  k = 2, m
    xr = .25_dp * xr * (2.0*k-3.0) ** 2 / ((k-1.0)*(2.0*k-1.0)**2) * x
    xg1 = xg1 + bj(k) * xr
  END DO
  xr = 1.0_dp
  xg2 = bj(1)
  DO  k = 2, m
    xr = .25_dp * xr * (2.0*k-5.0) ** 2 / ((k-1.0)*(2.0*k-3.0)**2) * x
    xg2 = xg2 + bj(k) * xr
  END DO
  xcs = COS(x/2.0_dp)
  xss = SIN(x/2.0_dp)
  ci = el + LOG(x) - x * xss * xg1 + 2 * xcs * xg2 - 2 * xcs * xcs
  si = x * xcs * xg1 + 2 * xss * xg2 - SIN(x)
ELSE
  xr = 1.0_dp
  xf = 1.0_dp
  DO  k = 1, 9
    xr = -2.0_dp * xr * k * (2*k-1) / x2
    xf = xf + xr
  END DO
  xr = 1.0_dp / x
  xg = xr
  DO  k = 1, 8
    xr = -2.0_dp * xr * (2*k+1) * k / x2
    xg = xg + xr
  END DO
  ci = xf * SIN(x) / x - xg * COS(x) / x
  si = p2 - xf * COS(x) / x - xg * SIN(x) / x
END IF
RETURN
END SUBROUTINE cisia
 
END MODULE cisia_func
 
 
 
PROGRAM mcisia
USE cisia_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:35

!      =====================================================
!      Purpose: This program computes the cosine and sine
!               integrals using subroutine CISIA
!      Input :  x  --- Argument of Ci(x) and Si(x)
!      Output:  CI --- Ci(x)
!               SI --- Si(x)
!      Example:
!                   x         Ci(x)          Si(x)
!                ------------------------------------
!                  0.0     - ì             .00000000
!                  5.0     -.19002975     1.54993124
!                 10.0     -.04545643     1.65834759
!                 20.0      .04441982     1.54824170
!                 30.0     -.03303242     1.56675654
!                 40.0      .01902001     1.58698512
!      =====================================================

REAL (dp)  :: ci, si, x

DO
  WRITE (*,*) 'Please enter x '
  READ (*,*) x
  WRITE (*,*) '   x         Ci(x)          Si(x)'
  WRITE (*,*) '------------------------------------'
  CALL cisia(x, ci, si)
  IF (x /= 0.0_dp) WRITE (*,5000) x, ci, si
  IF (x == 0.0_dp) WRITE (*,5100)
END DO
STOP

5000 FORMAT (' ', f5.1, 2F15.8)
5100 FORMAT (t4, ' .0     - ì', t28, '.00000000')
END PROGRAM mcisia
