MODULE sphy_func
 
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
 


SUBROUTINE sphy(n, x, nm, sy, dy)

!    ======================================================
!    Purpose: Compute spherical Bessel functions yn(x) and
!             their derivatives
!    Input :  x --- Argument of yn(x) ( x ע 0 )
!             n --- Order of yn(x) ( n = 0,1,תתת )
!    Output:  SY(n) --- yn(x)
!             DY(n) --- yn'(x)
!             NM --- Highest order computed
!    ======================================================

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x
INTEGER, INTENT(OUT)    :: nm
REAL (dp), INTENT(OUT)  :: sy(0:n)
REAL (dp), INTENT(OUT)  :: dy(0:n)

REAL (dp)  :: f, f0, f1
INTEGER    :: k

nm = n
IF (x < 1.0D-60) THEN
  DO  k = 0, n
    sy(k) = -1.0D+300
    dy(k) = 1.0D+300
  END DO
  RETURN
END IF
sy(0) = -COS(x) / x
sy(1) = (sy(0)-SIN(x)) / x
f0 = sy(0)
f1 = sy(1)
DO  k = 2, n
  f = (2*k-1)*f1/x - f0
  sy(k) = f
  IF (ABS(f) >= 1.0D+300) EXIT
  f0 = f1
  f1 = f
END DO
nm = k - 1
dy(0) = (SIN(x) + COS(x)/x) / x
DO  k = 1, nm
  dy(k) = sy(k-1) - (k+1)*sy(k)/x
END DO
RETURN
END SUBROUTINE sphy
 
END MODULE sphy_func
 
 
 
PROGRAM msphy
USE sphy_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:47

!    ======================================================
!    Purpose: This program computes the spherical Bessel
!             functions yn(x) and yn'(x) using subroutine
!             SPHY
!    Input :  x --- Argument of yn(x) ( x ע 0 )
!             n --- Order of yn(x) ( n = 0,1,תתת, ף 250 )
!    Output:  SY(n) --- yn(x)
!             DY(n) --- yn'(x)
!    Example:   x = 10.0
!               n          yn(x)               yn'(x)
!             --------------------------------------------
!               0     .8390715291D-01    -.6279282638D-01
!               1     .6279282638D-01     .7134858763D-01
!               2    -.6506930499D-01     .8231361788D-01
!               3    -.9532747888D-01    -.2693831344D-01
!               4    -.1659930220D-02    -.9449751377D-01
!               5     .9383354168D-01    -.5796005523D-01
!    ======================================================

REAL (dp)  :: sy(0:250), dy(0:250), x
INTEGER    :: k, n, nm, ns

WRITE(*,*) 'Please enter n and x '
READ(*,*) n, x
WRITE(*,5100) n, x
IF (n <= 10) THEN
  ns = 1
ELSE
  WRITE(*,*) 'Please enter order step Ns '
  READ(*,*) ns
END IF
CALL sphy(n, x, nm, sy, dy)
WRITE(*,*)
WRITE(*,*) '  n          yn(x)               yn''(X)'
WRITE(*,*) '--------------------------------------------'
DO  k = 0, nm, ns
  WRITE (*,5000) k, sy(k), dy(k)
END DO
STOP

5000 FORMAT(' ', i3, 2g20.10)
5100 FORMAT('   Nmax =', i3, ',     x=', f6.1)
END PROGRAM msphy
