MODULE eix_func
 
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


SUBROUTINE eix(x, ei)

!     ============================================
!     Purpose: Compute exponential integral Ei(x)
!     Input :  x  --- Argument of Ei(x)
!     Output:  EI --- Ei(x) ( x > 0 )
!     ============================================

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: ei

REAL (dp)  :: r
INTEGER    :: k

IF (x == 0.0) THEN
  ei = -1.0D+300
ELSE IF (x <= 40.0) THEN
  ei = 1.0_dp
  r = 1.0_dp
  DO  k = 1, 100
    r = r * k * x / (k+1) ** 2
    ei = ei + r
    IF (ABS(r/ei) <= 1.0D-15) EXIT
  END DO
  ei = 0.5772156649015328_dp + LOG(x) + x * ei
ELSE
  ei = 1.0_dp
  r = 1.0_dp
  DO  k = 1, 20
    r = r * k / x
    ei = ei + r
  END DO
  ei = EXP(x) / x * ei
END IF
RETURN
END SUBROUTINE eix
 
END MODULE eix_func
 
 
 
PROGRAM meix
USE eix_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:38

!     =======================================================
!     Purpose: This program computes the exponential integral
!              Ei(x) using subroutine EIX
!     Example:
!                x        Ei(x)
!              -----------------------
!                0    -.10000000+301
!                1     .18951178E+01
!                2     .49542344E+01
!                3     .99338326E+01
!                4     .19630874E+02
!                5     .40185275E+02
!     =======================================================

REAL (dp)  :: ei, x

WRITE(*,*) 'Please enter x '
READ(*,*) x
WRITE(*,*)
WRITE(*,*) '   x         Ei(x)'
WRITE(*,*) '------------------------'
CALL eix(x, ei)
WRITE(*,5000) x, ei
STOP

5000 FORMAT(' ', f5.1, e18.8)
END PROGRAM meix
