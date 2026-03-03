MODULE stvl1_func
 
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
 

SUBROUTINE stvl1(x, sl1)

!    ================================================
!    Purpose: Compute modified Struve function L1(x)
!    Input :  x   --- Argument of L1(x) ( x ò 0 )
!    Output:  SL1 --- L1(x)
!    ================================================

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: sl1

REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp
REAL (dp)  :: a1, bi1, r, s
INTEGER    :: k, km

r = 1.0_dp
IF (x <= 20.0_dp) THEN
  s = 0.0_dp
  DO  k = 1, 60
    r = r * x * x / (4*k*k-1)
    s = s + r
    IF (ABS(r/s) < 1.0D-12) EXIT
  END DO
  sl1 = 2.0_dp / pi * s
ELSE
  s = 1.0_dp
  km = .5*x
  IF (x > 50) km = 25
  DO  k = 1, km
    r = r * (2*k+3) * (2*k+1) / (x*x)
    s = s + r
    IF (ABS(r/s) < 1.0D-12) EXIT
  END DO
  sl1 = 2.0_dp / pi * (-1.0_dp + 1.0_dp/(x*x) + 3.0_dp*s/x**4)
  a1 = EXP(x) / SQRT(2.0_dp*pi*x)
  r = 1.0_dp
  bi1 = 1.0_dp
  DO  k = 1, 16
    r = -0.125_dp * r * (4 - (2*k-1)**2) / (k*x)
    bi1 = bi1 + r
    IF (ABS(r/bi1) < 1.0D-12) EXIT
  END DO
  sl1 = sl1 + a1 * bi1
END IF
RETURN
END SUBROUTINE stvl1
 
END MODULE stvl1_func
 
 
 
PROGRAM mstvl1
USE stvl1_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:47

!    =====================================================
!    Purpose: This program computes the modified Struve
!             function L1(x) using subroutine STVL1
!    Input :  x   --- Argument of L1(x) ( x ò 0 )
!    Output:  SL1 --- L1(x)
!    Example:
!                  x        L1(x)
!              -----------------------
!                0.0   .00000000D+00
!                5.0   .23728216D+02
!               10.0   .26703583D+04
!               15.0   .32812429D+06
!               20.0   .42454973D+08
!               30.0   .76853204D+12
!               40.0   .14707396D+17
!               50.0   .29030786D+21
!    =====================================================

REAL (dp)  :: sl1, x

WRITE (*,*) 'Please enter x '
READ (*,*) x
WRITE (*,*) '   x        L1(x)'
WRITE (*,*) '-----------------------'
CALL stvl1(x,sl1)
WRITE (*,5000) x, sl1
STOP

5000 FORMAT (' ', f5.1, g16.8)
END PROGRAM mstvl1
