MODULE itsl0_func
 
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
 

SUBROUTINE itsl0(x, tl0)

!    ===========================================================
!    Purpose: Evaluate the integral of modified Struve function
!             L0(t) with respect to t from 0 to x
!    Input :  x   --- Upper limit  ( x ò 0 )
!    Output:  TL0 --- Integration of L0(t) from 0 to x
!    ===========================================================

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: tl0

REAL (dp)  :: a(18), a0, a1, af, r, rd, s, s0, ti
INTEGER    :: k

REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp, el = .57721566490153_dp, &
                         one = 1.0_dp, two = 2.0_dp, half = 0.5_dp

r = one
IF (x <= 20.0) THEN
  s = half
  DO  k = 1, 100
    rd = one
    IF (k == 1) rd = half
    r = r * rd * k / (k+1) * (x/(2*k+1)) ** 2
    s = s + r
    IF (ABS(r/s) < 1.0D-12) EXIT
  END DO
  tl0 = two / pi * x * x * s
ELSE
  s = one
  DO  k = 1, 10
    r = r * k / (k+1) * ((2*k+1)/x) ** 2
    s = s + r
    IF (ABS(r/s) < 1.0D-12) EXIT
  END DO
  s0 = -s / (pi*x*x) + two / pi * (LOG(two*x)+el)
  a0 = one
  a1 = 0.625_dp
  a(1) = a1
  DO  k = 1, 10
    af = ((1.5D0*(k+half)*(k+5.0D0/6.0D0)*a1 - half*(k+half)**2 *  &
         (k-half)*a0)) / (k+1)
    a(k+1) = af
    a0 = a1
    a1 = af
  END DO
  ti = one
  r = one
  DO  k = 1, 11
    r = r / x
    ti = ti + a(k) * r
  END DO
  tl0 = ti / SQRT(2*pi*x) * EXP(x) + s0
END IF
RETURN
END SUBROUTINE itsl0
 
END MODULE itsl0_func
 
 
 
PROGRAM mitsl0
USE itsl0_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:41

!    ========================================================
!    Purpose: This program evaluates the integral of modified
!             Struve function L0(t) with respect to t from 0
!             to x using subroutine ITSL0
!    Input :  x   --- Upper limit  ( x ò 0 )
!    Output:  TL0 --- Integration of L0(t) from 0 to x
!    Example:
!                   x        L0(t)dt
!                -----------------------
!                  0.0    .0000000D+00
!                  5.0    .3003079D+02
!                 10.0    .2990773D+04
!                 15.0    .3526179D+06
!                 20.0    .4475860D+08
!                 30.0    .7955389D+12
!                 40.0    .1508972D+17
!                 50.0    .2962966D+21
!    ========================================================

REAL (dp)  :: tl0, x

WRITE (*,*) 'Please enter x '
READ (*,*) x
WRITE (*,*) '   x        L0(t)dt'
WRITE (*,*) '-----------------------'
CALL itsl0(x, tl0)
WRITE (*,5000) x, tl0
STOP

5000 FORMAT (' ', f5.1, g16.7)
END PROGRAM mitsl0
