MODULE lagzo_func
 
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
 

SUBROUTINE lagzo(n, x, w)

!    =========================================================
!    Purpose : Compute the zeros of Laguerre polynomial Ln(x)
!              in the interval [0,ì], and the corresponding
!              weighting coefficients for Gauss-Laguerre
!              integration
!    Input :   n    --- Order of the Laguerre polynomial
!              X(n) --- Zeros of the Laguerre polynomial
!              W(n) --- Corresponding weighting coefficients
!    =========================================================

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: x(n)
REAL (dp), INTENT(OUT)     :: w(n)

REAL (dp)  :: f0, f1, fd, gd, hn, p, pd, pf, q, wp, z, z0
INTEGER    :: i, it, j, k, nr

hn = 1.0D0 / n
DO  nr = 1, n
  IF (nr == 1) z = hn
  IF (nr > 1) z = x(nr-1) + hn * nr ** 1.27
  it = 0

  10 it = it + 1
  z0 = z
  p = 1.0D0
  DO  i = 1, nr - 1
    p = p * (z-x(i))
  END DO
  f0 = 1.0D0
  f1 = 1.0D0 - z
  DO  k = 2, n
    pf = ((2*k-1-z)*f1 - (k-1)*f0) / k
    pd = k / z * (pf-f1)
    f0 = f1
    f1 = pf
  END DO
  fd = pf / p
  q = 0.0D0
  DO  i = 1, nr - 1
    wp = 1.0D0
    DO  j = 1, nr - 1
      IF (j /= i) THEN
        wp = wp * (z-x(j))
      END IF
    END DO
    q = q + wp
  END DO
  gd = (pd-q*fd) / p
  z = z - fd / gd
  IF (it <= 40 .AND. ABS((z-z0)/z) > 1.0D-15) GO TO 10

  x(nr) = z
  w(nr) = 1.0D0 / (z*pd*pd)
END DO

RETURN
END SUBROUTINE lagzo
 
END MODULE lagzo_func
 
 
 
PROGRAM mlagzo
USE lagzo_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:43

!    ===========================================================
!    Purpose : This program computes the zeros of Laguerre
!              polynomial Ln(x) in the interval [0,ì] and the
!              corresponding weighting coefficients for Gauss-
!              Laguerre integration using subroutine LAGZO
!    Input :   n    --- Order of the Laguerre polynomial
!              X(n) --- Zeros of the Laguerre polynomial
!              W(n) --- Corresponding weighting coefficients
!    ===========================================================

REAL (dp)  :: x(100), w(100)
INTEGER    :: j, n

WRITE (*,*) 'Please enter the order of Ln(x), n: '
READ (*,*) n
WRITE (*,5000) n
CALL lagzo(n, x, w)
WRITE (*,*) '  Nodes and weights for Gauss-Laguerre integration'
WRITE (*,*)
WRITE (*,*) '  i             xi                      Wi'
WRITE (*,*) ' -----------------------------------------------------'
DO  j = 1, n
  WRITE (*,5100) j, x(j), w(j)
END DO
STOP

5000 FORMAT (' n =', i3)
5100 FORMAT (' ', i3, '   ', g22.13, '   ', g22.13)
END PROGRAM mlagzo
