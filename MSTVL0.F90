MODULE stvl0_func
 
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


SUBROUTINE stvl0(x, sl0)

!    ================================================
!    Purpose: Compute modified Struve function L0(x)
!    Input :  x   --- Argument of L0(x) ( x ò 0 )
!    Output:  SL0 --- L0(x)
!    ================================================

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: sl0

REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp
REAL (dp)  :: a0, a1, bi0, r, s
INTEGER    :: k, km

s = 1.0_dp
r = 1.0_dp
IF (x <= 20.0_dp) THEN
  a0 = 2.0_dp * x / pi
  DO  k = 1, 60
    r = r * (x/(2*k+1)) ** 2
    s = s + r
    IF (ABS(r/s) < 1.0D-12) EXIT
  END DO
  sl0 = a0 * s
ELSE
  km = .5*(x+1.0)
  IF (x >= 50.0) km = 25
  DO  k = 1, km
    r = r * ((2*k-1)/x) ** 2
    s = s + r
    IF (ABS(r/s) < 1.0D-12) EXIT
  END DO
  a1 = EXP(x) / SQRT(2.0_dp*pi*x)
  r = 1.0_dp
  bi0 = 1.0_dp
  DO  k = 1, 16
    r = 0.125_dp * r * (2*k-1) ** 2 / (k*x)
    bi0 = bi0 + r
    IF (ABS(r/bi0) < 1.0D-12) EXIT
  END DO
  bi0 = a1 * bi0
  sl0 = -2.0_dp / (pi*x) * s + bi0
END IF
RETURN
END SUBROUTINE stvl0
 
END MODULE stvl0_func
 
 
 
PROGRAM mstvl0
USE stvl0_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:47

!    =====================================================
!    Purpose: This program computes modified Struve
!             function L0(x) using subroutine STVL0
!    Input :  x   --- Argument of L0(x) ( x ò 0 )
!    Output:  SL0 --- L0(x)
!    Example:
!                x        L0(x)
!            ------------------------
!               0.0   .00000000D+00
!               5.0   .27105917D+02
!              10.0   .28156522D+04
!              15.0   .33964933D+06
!              20.0   .43558283D+08
!              30.0   .78167230D+12
!              40.0   .14894775D+17
!              50.0   .29325538D+21
!    =====================================================

REAL (dp)  :: sl0, x

WRITE (*,*) 'Please enter x '
READ (*,*) x
WRITE (*,*) '   x        L0(x)'
WRITE (*,*) '-----------------------'
CALL stvl0(x, sl0)
WRITE (*,5000) x, sl0
STOP

5000 FORMAT (' ', f5.1, g16.8)
END PROGRAM mstvl0
