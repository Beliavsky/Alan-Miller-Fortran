MODULE stvh1_func
 
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


SUBROUTINE stvh1(x, sh1)

!    =============================================
!    Purpose: Compute Struve function H1(x)
!    Input :  x   --- Argument of H1(x) ( x ò 0 )
!    Output:  SH1 --- H1(x)
!    =============================================

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: sh1

REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp
REAL (dp)  :: a0, by1, p1, q1, r, s, t, t2, ta1
INTEGER    :: k, km

r = 1.0_dp
IF (x <= 20.0_dp) THEN
  s = 0.0_dp
  a0 = -2.0_dp / pi
  DO  k = 1, 60
    r = -r * x * x / (4*k*k-1)
    s = s + r
    IF (ABS(r) < ABS(s)*1.0D-12) EXIT
  END DO
  sh1 = a0 * s
ELSE
  s = 1.0_dp
  km = INT(.5*x)
  IF (x > 50._dp) km = 25
  DO  k = 1, km
    r = -r * (4*k*k-1) / (x*x)
    s = s + r
    IF (ABS(r) < ABS(s)*1.0D-12) EXIT
  END DO
  t = 4.0_dp / x
  t2 = t * t
  p1 = ((((.42414D-5*t2 - .20092D-4)*t2 + .580759D-4)*t2 - .223203D-3)*  &
       t2 + .29218256D-2)*t2 + .3989422819_dp
  q1 = t * (((((-.36594D-5*t2 + .1622D-4)*t2 - .398708D-4)*t2 +  &
       .1064741D-3)*t2 - .63904D-3)*t2 + .0374008364_dp)
  ta1 = x - .75_dp * pi
  by1 = 2.0_dp / SQRT(x) * (p1*SIN(ta1) + q1*COS(ta1))
  sh1 = 2.0 / pi * (1.0_dp + s/(x*x)) + by1
END IF
RETURN
END SUBROUTINE stvh1
 
END MODULE stvh1_func
 
 
 
PROGRAM mstvh1
USE stvh1_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:47

!    ===============================================
!    Purpose: This program computes Struve function
!             H1(x) using subroutine STVH1
!    Input :  x   --- Argument of H1(x) ( x ò 0 )
!    Output:  SH1 --- H1(x)
!    Example:
!                x          H1(x)
!             -----------------------
!               0.0       .00000000
!               5.0       .80781195
!              10.0       .89183249
!              15.0       .66048730
!              20.0       .47268818
!              25.0       .53880362
!    ===============================================

REAL (dp)  :: sh1, x

WRITE (*,*) 'Please enter x '
READ (*,*) x
WRITE (*,*) '   x          H1(x)'
WRITE (*,*) '-----------------------'
CALL stvh1(x, sh1)
WRITE (*,5000) x, sh1
STOP

5000 FORMAT (' ', f5.1, g16.8)
END PROGRAM mstvh1
