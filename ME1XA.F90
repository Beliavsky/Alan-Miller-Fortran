MODULE e1xa_func
 
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


SUBROUTINE e1xa(x, e1)

!       ============================================
!       Purpose: Compute exponential integral E1(x)
!       Input :  x  --- Argument of E1(x)
!       Output:  E1 --- E1(x) ( x > 0 )
!       ============================================

REAL (dp), INTENT(IN)      :: x
REAL (dp), INTENT(OUT)     :: e1

REAL (dp)  :: es1, es2

IF (x == 0.0) THEN
  e1 = 1.0D+300
ELSE IF (x <= 1.0) THEN
  e1 = -LOG(x) + ((((1.07857D-3*x - 9.76004D-3)*x + 5.519968D-2)*x -  &
        0.24991055D0)*x + 0.99999193D0) * x - 0.57721566D0
ELSE
  es1 = (((x + 8.5733287401D0)*x + 18.059016973D0)*x + 8.6347608925D0) *  &
        x + 0.2677737343D0
  es2 = (((x + 9.5733223454D0)*x + 25.6329561486D0)*x + 21.0996530827D0)  &
        * x + 3.9584969228D0
  e1 = EXP(-x) / x * es1 / es2
END IF
RETURN
END SUBROUTINE e1xa
 
END MODULE e1xa_func
 
 
 
PROGRAM me1xa
USE e1xa_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:38

!       =========================================================
!       Purpose: This program computes the exponential integral
!                E1(x) using subroutine E1XA
!       Input :  x  --- Argument of E1(x)  ( x > 0 )
!       Output:  E1 --- E1(x)
!       Example:
!                  x        E1(x)
!                ----------------------
!                 0.0     .1000000+301
!                 1.0     .2193839E+00
!                 2.0     .4890051E-01
!                 3.0     .1304838E-01
!                 4.0     .3779352E-02
!                 5.0     .1148296E-02
!       =========================================================

REAL (dp)  :: e1, x

WRITE (*,*) 'Please enter x '
READ (*,*) x
WRITE (*,*) '   x        E1(x)'
WRITE (*,*) ' ----------------------'
CALL e1xa(x, e1)
WRITE (*,5000) x, e1
STOP

5000 FORMAT (' ', f5.1, g17.7)
END PROGRAM me1xa
