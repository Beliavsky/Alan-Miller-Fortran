MODULE e1xb_func
 
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


SUBROUTINE e1xb(x, e1)

!       ============================================
!       Purpose: Compute exponential integral E1(x)
!       Input :  x  --- Argument of E1(x)
!       Output:  E1 --- E1(x)  ( x > 0 )
!       ============================================

REAL (dp), INTENT(IN)      :: x
REAL (dp), INTENT(OUT)     :: e1

REAL (dp)  :: ga, r, t, t0
INTEGER    :: k, m

IF (x == 0.0) THEN
  e1 = 1.0D+300
ELSE IF (x <= 1.0) THEN
  e1 = 1.0_dp
  r = 1.0_dp
  DO  k = 1, 25
    r = -r * k * x / (k+1) ** 2
    e1 = e1 + r
    IF (ABS(r) <= ABS(e1)*1.0D-15) EXIT
  END DO
  ga = 0.5772156649015328_dp
  e1 = -ga - LOG(x) + x * e1
ELSE
  m = 20 + INT(80.0/x)
  t0 = 0.0_dp
  DO  k = m, 1, -1
    t0 = k / (1.0_dp + k/(x+t0))
  END DO
  t = 1.0_dp / (x+t0)
  e1 = EXP(-x) * t
END IF
RETURN
END SUBROUTINE e1xb
 
END MODULE e1xb_func
 
 
 
PROGRAM me1xb
USE e1xb_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:38

!       =========================================================
!       Purpose: This program computes the exponential integral
!                E1(x) using subroutine E1XB
!       Input :  x  --- Argument of E1(x)  ( x > 0 )
!       Output:  E1 --- E1(x)
!       Example:
!                  x          E1(x)
!                -------------------------
!                 0.0     .1000000000+301
!                 1.0     .2193839344E+00
!                 2.0     .4890051071E-01
!                 3.0     .1304838109E-01
!                 4.0     .3779352410E-02
!                 5.0     .1148295591E-02
!       =========================================================

REAL (dp)  :: e1, x

WRITE (*,*) 'Please enter x '
READ (*,*) x
WRITE (*,*) '   x          E1(x)'
WRITE (*,*) ' -------------------------'
CALL e1xb(x, e1)
WRITE (*,5000) x, e1
STOP

5000 FORMAT (' ', f5.1, e20.10)
END PROGRAM me1xb
