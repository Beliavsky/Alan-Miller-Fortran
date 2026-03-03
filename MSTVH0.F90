MODULE stvh0_func
 
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
 

SUBROUTINE stvh0(x, sh0)

!    =============================================
!    Purpose: Compute Struve function H0(x)
!    Input :  x   --- Argument of H0(x) ( x ò 0 )
!    Output:  SH0 --- H0(x)
!    =============================================

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: sh0

REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp
REAL (dp)  :: a0, by0, p0, q0, r, s, t, t2, ta0
INTEGER    :: k, km

s = 1.0D0
r = 1.0D0
IF (x <= 20.0D0) THEN
  a0 = 2.0 * x / pi
  DO  k = 1, 60
    r = -r * x / (2*k+1) * x / (2*k+1)
    s = s + r
    IF (ABS(r) < ABS(s)*1.0D-12) EXIT
  END DO
  sh0 = a0 * s
ELSE
  km = INT(.5*(x+1.0))
  IF (x >= 50.0) km = 25
  DO  k = 1, km
    r = -r * ((2*k-1)/x) ** 2
    s = s + r
    IF (ABS(r) < ABS(s)*1.0D-12) EXIT
  END DO
  t = 4.0D0 / x
  t2 = t * t
  p0 = ((((-.37043D-5*t2 + .173565D-4)*t2 - .487613D-4)*t2 + .17343D-3)*  &
       t2 - .1753062D-2)*t2 + .3989422793_dp
  q0 = t * (((((.32312D-5*t2 - .142078D-4)*t2 + .342468D-4)*t2 -  &
       .869791D-4)*t2 + .4564324D-3)*t2 - .0124669441_dp)
  ta0 = x - .25D0 * pi
  by0 = 2.0D0 / SQRT(x) * (p0*SIN(ta0) + q0*COS(ta0))
  sh0 = 2.0D0 / (pi*x) * s + by0
END IF
RETURN
END SUBROUTINE stvh0
 
END MODULE stvh0_func
 
 
 
PROGRAM mstvh0
USE stvh0_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:47

!    ====================================================
!    Purpose: This program computes Struve function
!             H0(x) using subroutine STVH0
!    Input :  x   --- Argument of H0(x) ( x ò 0 )
!    Output:  SH0 --- H0(x)
!    Example:
!                x          H0(x)
!             ----------------------
!               0.0       .00000000
!               5.0      -.18521682
!              10.0       .11874368
!              15.0       .24772383
!              20.0       .09439370
!              25.0      -.10182519
!    ====================================================

REAL (dp)  :: sh0, x

WRITE (*,*) 'Please enter x '
READ (*,*) x
WRITE (*,*) '   x          H0(x)'
WRITE (*,*) '----------------------'
CALL stvh0(x, sh0)
WRITE (*,5000) x, sh0
STOP

5000 FORMAT (' ', f5.1, g16.8)
END PROGRAM mstvh0
