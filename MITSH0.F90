MODULE itsh0_func
 
! From the book "Computation of Special Functions"
!      by Shanjie Zhang and Jianming Jin
!   Copyright 1996 by John Wiley & Sons, Inc.
! The authors state:
!   "However, we give permission to the reader who purchases this book
!    to incorporate any of these programs into his or her programs
!    provided that the copyright is acknowledged."

! The results for x >= 30. agree with those in the program only to about
! 6 decimal places.

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
 
CONTAINS


SUBROUTINE itsh0(x, th0)

!       ===================================================
!       Purpose: Evaluate the integral of Struve function
!                H0(t) with respect to t from 0 and x
!       Input :  x   --- Upper limit  ( x ò 0 )
!       Output:  TH0 --- Integration of H0(t) from 0 and x
!       ===================================================

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: th0

REAL (dp)  :: a(25), a0, a1, af, bf, bg, r, rd, s, s0, ty, xp
INTEGER    :: k
REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp, el = .57721566490153_dp

r = 1.0_dp
IF (x <= 30.0_dp) THEN
  s = 0.5_dp
  DO  k = 1, 100
    rd = 1.0_dp
    IF (k == 1) rd = 0.5_dp
    r = -r * rd * k / (k+1) * (x/(2*k+1)) ** 2
    s = s + r
    IF (ABS(r) < ABS(s)*1.0D-12) EXIT
  END DO
  th0 = 2.0_dp / pi * x * x * s
ELSE
  s = 1.0_dp
  DO  k = 1, 12
    r = -r * k / (k+1) * ((2*k+1)/x) ** 2
    s = s + r
    IF (ABS(r) < ABS(s)*1.0D-12) EXIT
  END DO

  s0 = s / (pi*x*x) + 2.0_dp / pi * (LOG(2.0_dp*x) + el)
  a0 = 1.0_dp
  a1 = 0.625_dp
  a(1) = a1
  DO  k = 1, 20
    af = ((1.5_dp*(k + .5_dp)*(k + 5.0_dp/6.0_dp)*a1 - .5_dp*(k + .5_dp)*(k + &
         .5_dp)*(k - .5_dp)*a0)) / (k+1)
    a(k+1) = af
    a0 = a1
    a1 = af
  END DO
  bf = 1.0_dp
  r = 1.0_dp
  DO  k = 1, 10
    r = -r / (x*x)
    bf = bf + a(2*k) * r
  END DO
  bg = a(1) / x
  r = 1.0_dp / x
  DO  k = 1, 10
    r = -r / (x*x)
    bg = bg + a(2*k+1) * r
  END DO
  xp = x + .25_dp * pi
  ty = SQRT(2.0_dp/(pi*x)) * (bg*COS(xp) - bf*SIN(xp))
  th0 = ty + s0
END IF
RETURN
END SUBROUTINE itsh0
 
END MODULE itsh0_func
 
 
 
PROGRAM mitsh0
USE itsh0_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:41

!    ====================================================
!    Purpose: This program evaluates the integral of
!             Struve function H0(t) with respect to t
!             from 0 and x using subroutine ITSH0
!    Input :  x   --- Upper limit  ( x ò 0 )
!    Output:  TH0 --- Integration of H0(t) from 0 and x
!    Example:
!                 x        H0(t)dt
!              ----------------------
!                0.0       .0000000
!                5.0      2.0442437
!               10.0      2.5189577
!               15.0      2.5415824
!               20.0      2.5484517
!               30.0      3.0625848
!               40.0      3.1484123
!               50.0      3.2445168
!    ====================================================

REAL (dp)  :: th0, x

WRITE (*,*) 'Please enter x '
READ (*,*) x
WRITE (*,*) '   x        H0(t)dt'
WRITE (*,*) '----------------------'
CALL itsh0(x, th0)
WRITE (*,5000) x, th0
STOP

5000 FORMAT (' ', f5.1, g16.7)
END PROGRAM mitsh0
