MODULE itth0_func
 
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
 

SUBROUTINE itth0(x, tth)

!    ===========================================================
!    Purpose: Evaluate the integral H0(t)/t with respect to t
!             from x to infinity
!    Input :  x   --- Lower limit  ( x ò 0 )
!    Output:  TTH --- Integration of H0(t)/t from x to infinity
!    ===========================================================

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: tth

REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp, one = 1.0_dp, two = 2.0_dp
REAL (dp)  :: f0, g0, r, s, t, tty, xt
INTEGER    :: k

s = one
r = one
IF (x < 24.5D0) THEN
  DO  k = 1, 60
    r = -r * x * x * (2*k - one) / (2*k+1) ** 3
    s = s + r
    IF (ABS(r) < ABS(s)*1.0D-12) EXIT
  END DO
  tth = pi / two - two / pi * x * s
ELSE
  DO  k = 1, 10
    r = -r * (2*k-1) ** 3 / ((2*k+1)*x*x)
    s = s + r
    IF (ABS(r) < ABS(s)*1.0D-12) EXIT
  END DO
  tth = two / (pi*x) * s
  t = 8.0D0 / x
  xt = x + .25D0 * pi
  f0 = (((((.18118D-2*t - .91909D-2)*t + .017033D0)*t - .9394D-3)*t -   &
       .051445D0)*t - .11D-5) * t + .7978846D0
  g0 = (((((-.23731D-2*t + .59842D-2)*t + .24437D-2)*t - .0233178D0)*t +   &
       .595D-4)*t + .1620695D0) * t
  tty = (f0*SIN(xt) - g0*COS(xt)) / (SQRT(x)*x)
  tth = tth + tty
END IF
RETURN
END SUBROUTINE itth0
 
END MODULE itth0_func
 
 
 
PROGRAM mitth0
USE itth0_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:41

!    ===========================================================
!    Purpose: This program evaluates the integral of H0(t)/t
!             with respect to t from x to infinity using
!             subroutine ITTH0
!    Input :  x   --- Lower limit  ( x ò 0 )
!    Output:  TTH --- Integration of H0(t)/t from x to infinity
!    Example:
!                 x        H0(t)/t dt
!              -----------------------
!                0.0      1.57079633
!                5.0       .07954575
!               10.0       .04047175
!               15.0       .04276558
!               20.0       .04030796
!               30.0       .01815256
!               40.0       .01621331
!               50.0       .01378661
!    =======================================================

REAL (dp)  :: tth, x

WRITE (*,*) 'Please enter x '
READ (*,*) x
WRITE (*,*) '   x       H0(t)/t dt'
WRITE (*,*) '-----------------------'
CALL itth0(x, tth)
WRITE (*,5000) x, tth
STOP

5000 FORMAT (' ', f5.1, ' ', g16.8)
END PROGRAM mitth0
