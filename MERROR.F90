MODULE error_func
 
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


SUBROUTINE error(x, ERR)

!       =========================================
!       Purpose: Compute error function erf(x)
!       Input:   x   --- Argument of erf(x)
!       Output:  ERR --- erf(x)
!       =========================================

REAL (dp), INTENT(IN)      :: x
REAL (dp), INTENT(OUT)     :: ERR

REAL (dp)  :: c0, eps, er, r, x2
INTEGER    :: k
REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp

eps = 1.0D-15
x2 = x * x
IF (ABS(x) < 3.5_dp) THEN
  er = 1.0_dp
  r = 1.0_dp
  DO  k = 1, 50
    r = r * x2 / (k+0.5_dp)
    er = er + r
    IF (ABS(r) <= ABS(er)*eps) EXIT
  END DO
  c0 = 2.0_dp / SQRT(pi) * x * EXP(-x2)
  ERR = c0 * er
ELSE
  er = 1.0_dp
  r = 1.0_dp
  DO  k = 1, 12
    r = -r * (k-0.5_dp) / x2
    er = er + r
  END DO
  c0 = EXP(-x2) / (ABS(x)*SQRT(pi))
  ERR = 1.0_dp - c0 * er
  IF (x < 0.0) ERR = -ERR
END IF
RETURN
END SUBROUTINE error
 
END MODULE error_func
 
 
 
PROGRAM merror
USE error_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:38

!       ===================================================
!       Purpose: This program computes the error function
!                erf(x) using subroutine ERROR
!       Input:   x   --- Argument of erf(x)
!       Output:  ERR --- erf(x)
!       Example:
!                  x         erf(x)
!                ---------------------
!                 1.0       .84270079
!                 2.0       .99532227
!                 3.0       .99997791
!                 4.0       .99999998
!                 5.0      1.00000000
!       ===================================================

REAL (dp)  :: er, x

DO
  WRITE (*,*) 'Please enter x '
  READ (*,*) x
  WRITE (*,*) '   x         erf(x)'
  WRITE (*,*) '---------------------'
  CALL error(x, er)
  WRITE (*,5000) x, er
END DO
STOP

5000 FORMAT (' ', f5.2, f15.8)
END PROGRAM merror
