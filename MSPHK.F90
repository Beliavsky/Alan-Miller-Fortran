MODULE sphk_func
 
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


SUBROUTINE sphk(n, x, nm, sk, dk)

!    =====================================================
!    Purpose: Compute modified spherical Bessel functions
!             of the second kind, kn(x) and kn'(x)
!    Input :  x --- Argument of kn(x)  ( x ò 0 )
!             n --- Order of kn(x) ( n = 0,1,2,... )
!    Output:  SK(n) --- kn(x)
!             DK(n) --- kn'(x)
!             NM --- Highest order computed
!    =====================================================

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x
INTEGER, INTENT(OUT)    :: nm
REAL (dp), INTENT(OUT)  :: sk(0:n)
REAL (dp), INTENT(OUT)  :: dk(0:n)

REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp
REAL (dp)  :: f, f0, f1
INTEGER    :: k

nm = n
IF (x < 1.0D-60) THEN
  DO  k = 0, n
    sk(k) = 1.0D+300
    dk(k) = -1.0D+300
  END DO
  RETURN
END IF
sk(0) = 0.5_dp * pi / x * EXP(-x)
sk(1) = sk(0) * (1.0_dp + 1.0_dp/x)
f0 = sk(0)
f1 = sk(1)
DO  k = 2, n
  f = (2*k-1) * f1 / x + f0
  sk(k) = f
  IF (ABS(f) > 1.0D+300) EXIT
  f0 = f1
  f1 = f
END DO
nm = k - 1
dk(0) = -sk(1)
DO  k = 1, nm
  dk(k) = -sk(k-1) - (k+1) / x * sk(k)
END DO
RETURN
END SUBROUTINE sphk
 
END MODULE sphk_func
 
 
 
PROGRAM msphk
USE sphk_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:47

!    ======================================================
!    Purpose: This program computes the modified spherical
!             Bessel functions kn(x) and kn'(x) using
!             subroutine SPHK
!    Input :  x --- Argument of kn(x)  ( x ò 0 )
!             n --- Order of kn(x) ( n ó 250 )
!    Output:  SK(n) --- kn(x)
!             DK(n) --- kn'(x)
!    Example: x= 10.0
!               n          kn(x)               kn'(x)
!             --------------------------------------------
!               0     .7131404291D-05    -.7844544720D-05
!               1     .7844544720D-05    -.8700313235D-05
!               2     .9484767707D-05    -.1068997503D-04
!               3     .1258692857D-04    -.1451953914D-04
!               4     .1829561771D-04    -.2173473743D-04
!               5     .2905298451D-04    -.3572740841D-04
!    ======================================================


REAL (dp)  :: sk(0:250), dk(0:250), x
INTEGER    :: k, n, nm, ns

WRITE (*,*) 'Please enter n and x '
READ (*,*) n, x
WRITE (*,5100) n, x
IF (n <= 10) THEN
  ns = 1
ELSE
  WRITE (*,*) 'Please enter order step Ns '
  READ (*,*) ns
END IF
CALL sphk(n, x, nm, sk, dk)
WRITE (*,*)
WRITE (*,*) '  n          kn(x)               kn''(X)'
WRITE (*,*) '--------------------------------------------'
DO  k = 0, nm, ns
  WRITE (*,5000) k, sk(k), dk(k)
END DO
STOP

5000 FORMAT (' ', i3, 2g20.10)
5100 FORMAT ('   Nmax =', i3, ',     x =', f6.1)
END PROGRAM msphk
