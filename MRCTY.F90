MODULE rcty_func
 
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


SUBROUTINE rcty(n, x, nm, ry, dy)

!    ========================================================
!    Purpose: Compute Riccati-Bessel functions of the second
!             kind and their derivatives
!    Input:   x --- Argument of Riccati-Bessel function
!             n --- Order of yn(x)
!    Output:  RY(n) --- xúyn(x)
!             DY(n) --- [xúyn(x)]'
!             NM --- Highest order computed
!    ========================================================

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x
INTEGER, INTENT(OUT)    :: nm
REAL (dp), INTENT(OUT)  :: ry(0:n)
REAL (dp), INTENT(OUT)  :: dy(0:n)

REAL (dp)  :: rf0, rf1, rf2
INTEGER    :: k

nm = n
IF (x < 1.0D-60) THEN
  DO  k = 0, n
    ry(k) = -1.0D+300
    dy(k) = 1.0D+300
  END DO
  ry(0) = -1.0D0
  dy(0) = 0.0D0
  RETURN
END IF
ry(0) = -COS(x)
ry(1) = ry(0) / x - SIN(x)
rf0 = ry(0)
rf1 = ry(1)
DO  k = 2, n
  rf2 = (2*k-1) * rf1 / x - rf0
  IF (ABS(rf2) > 1.0D+300) EXIT
  ry(k) = rf2
  rf0 = rf1
  rf1 = rf2
END DO
nm = k - 1
dy(0) = SIN(x)
DO  k = 1, nm
  dy(k) = -k * ry(k) / x + ry(k-1)
END DO
RETURN
END SUBROUTINE rcty
 
END MODULE rcty_func
 
 
 
PROGRAM mrcty
USE rcty_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:45

!    =======================================================
!    Purpose: This program computes the Riccati-Bessel
!             functions of the second kind and their
!             derivatives using subroutine RCTY
!    Input:   x --- Argument of Riccati-Bessel function
!             n --- Order of yn(x)
!    Output:  RY(n) --- xúyn(x)
!             DY(n) --- [xúyn(x)]'
!    Example: x = 10.0
!               n        xúyn(x)             [xúyn(x)]'
!             --------------------------------------------
!               0     .8390715291D+00    -.5440211109D+00
!               1     .6279282638D+00     .7762787027D+00
!               2    -.6506930499D+00     .7580668738D+00
!               3    -.9532747888D+00    -.3647106133D+00
!               4    -.1659930220D-01    -.9466350679D+00
!               5     .9383354168D+00    -.4857670106D+00
!    =======================================================

REAL (dp)  :: ry(0:250), dy(0:250), x
INTEGER    :: k, n, nm, ns

WRITE (*,*) '  Please enter n and x '
READ (*,*) n, x
WRITE (*,5100) n, x
IF (n <= 10) THEN
  ns = 1
ELSE
  WRITE (*,*) '  Please enter order step Ns '
  READ (*,*) ns
END IF
WRITE (*,*)
CALL rcty(n, x, nm, ry, dy)
WRITE (*,*)
WRITE (*,*) '  n        xúyn(x)             [xúyn(x)]'''
WRITE (*,*) '--------------------------------------------'
DO  k = 0, nm, ns
  WRITE (*,5000) k, ry(k), dy(k)
END DO
STOP

5000 FORMAT (' ', i3, 2g20.10)
5100 FORMAT ('   Nmax =', i3, ',    x =', f6.2)
END PROGRAM mrcty
