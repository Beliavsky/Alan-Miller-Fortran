MODULE othpl_func
 
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
 

SUBROUTINE othpl(kf, n, x, pl, dpl)

!    ==========================================================
!    Purpose: Compute orthogonal polynomials: Tn(x) or Un(x),
!             or Ln(x) or Hn(x), and their derivatives
!    Input :  KF --- Function code
!                    KF=1 for Chebyshev polynomial Tn(x)
!                    KF=2 for Chebyshev polynomial Un(x)
!                    KF=3 for Laguerre polynomial Ln(x)
!                    KF=4 for Hermite polynomial Hn(x)
!             n ---  Order of orthogonal polynomials
!             x ---  Argument of orthogonal polynomials
!    Output:  PL(n) --- Tn(x) or Un(x) or Ln(x) or Hn(x)
!             DPL(n)--- Tn'(x) or Un'(x) or Ln'(x) or Hn'(x)
!    =========================================================

INTEGER, INTENT(IN)     :: kf
INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: pl(0:n)
REAL (dp), INTENT(OUT)  :: dpl(0:n)

REAL (dp)  :: a, b, c, dy0, dy1, dyn, y0, y1, yn
INTEGER    :: k

a = 2.0_dp
b = 0.0_dp
c = 1.0_dp
y0 = 1.0_dp
y1 = 2.0_dp * x
dy0 = 0.0_dp
dy1 = 2.0_dp
pl(0) = 1.0_dp
pl(1) = 2.0_dp * x
dpl(0) = 0.0_dp
dpl(1) = 2.0_dp
IF (kf == 1) THEN
  y1 = x
  dy1 = 1.0_dp
  pl(1) = x
  dpl(1) = 1.0_dp
ELSE IF (kf == 3) THEN
  y1 = 1.0_dp - x
  dy1 = -1.0_dp
  pl(1) = 1.0_dp - x
  dpl(1) = -1.0_dp
END IF
DO  k = 2, n
  IF (kf == 3) THEN
    a = -1.0_dp / k
    b = 2.0_dp + a
    c = 1.0_dp + a
  ELSE IF (kf == 4) THEN
    c = 2 * (k-1)
  END IF
  yn = (a*x+b) * y1 - c * y0
  dyn = a * y1 + (a*x + b) * dy1 - c * dy0
  pl(k) = yn
  dpl(k) = dyn
  y0 = y1
  y1 = yn
  dy0 = dy1
  dy1 = dyn
END DO
RETURN
END SUBROUTINE othpl
 
END MODULE othpl_func
 
 
 
PROGRAM mothpl
USE othpl_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:44

!    =========================================================
!    Purpose: This program computes orthogonal polynomials:
!             Tn(x) or Un(x) or Ln(x) or Hn(x), and their
!             derivatives using subroutine OTHPL
!    Input :  KF --- Function code
!                    KF=1 for Chebyshev polynomial Tn(x)
!                    KF=2 for Chebyshev polynomial Un(x)
!                    KF=3 for Laguerre polynomial Ln(x)
!                    KF=4 for Hermite polynomial Hn(x)
!             n ---  Order of orthogonal polynomials
!             x ---  Argument
!    Output:  PL(n) --- Tn(x) or Un(x) or Ln(x) or Hn(x)
!             DPL(n)--- Tn'(x) or Un'(x) or Ln'(x) or Hn'(x)
!                       n = 0,1,2,...,N ( N ó 100 )
!    =========================================================

REAL (dp)  :: pl(0:100), dpl(0:100), x
INTEGER    :: k, kf, n

WRITE (*,*) 'KF,N,x = ? '
READ (*,*) kf, n, x
WRITE (*,5000) kf, n, x
WRITE (*,*)
CALL othpl(kf, n, x, pl, dpl)
IF (kf == 1) WRITE (*,*) '  n          Tn(x)            Tn''(X)'
IF (kf == 2) WRITE (*,*) '  n          Un(x)            Un''(X)'
IF (kf == 3) WRITE (*,*) '  n          Ln(x)            Ln''(X)'
IF (kf == 4) WRITE (*,*) '  n          Hn(x)            Hn''(X)'
WRITE (*,*) '-----------------------------------------'
DO  k = 0, n
  WRITE (*,5100) k, pl(k), dpl(k)
END DO
STOP

5000 FORMAT(' KF=', i3, '     Nmax=', i3, '     x=', f6.3)
5100 FORMAT(' ', i3, 2g18.8)
END PROGRAM mothpl
