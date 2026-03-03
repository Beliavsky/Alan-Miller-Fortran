MODULE herzo_func
 
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
 

SUBROUTINE herzo(n, x, w)

!       ========================================================
!       Purpose : Compute the zeros of Hermite polynomial Ln(x)
!                 in the interval [-ì,ì], and the corresponding
!                 weighting coefficients for Gauss-Hermite
!                 integration
!       Input :   n    --- Order of the Hermite polynomial
!       Output:   X(n) --- Zeros of the Hermite polynomial
!                 W(n) --- Corresponding weighting coefficients
!       ========================================================

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(OUT)     :: x(n)
REAL (dp), INTENT(OUT)     :: w(n)

REAL (dp)  :: f0, f1, fd, gd, hd, hf, hn, p, q, r, r1, r2, wp, z, z0, zl
INTEGER    :: i, it, j, k, nr

hn = 1.0_dp / n
zl = -1.1611_dp + 1.46_dp * n ** 0.5
DO  nr = 1, n / 2
  IF (nr == 1) z = zl
  IF (nr /= 1) z = z - hn * (n/2+1-nr)
  it = 0

  10 it = it + 1
  z0 = z
  f0 = 1.0_dp
  f1 = 2.0_dp * z
  DO  k = 2, n
    hf = 2.0_dp * z * f1 - 2 * (k-1) * f0
    hd = 2.0_dp * k * f1
    f0 = f1
    f1 = hf
  END DO
  p = 1.0_dp
  DO  i = 1, nr - 1
    p = p * (z-x(i))
  END DO
  fd = hf / p
  q = 0.0_dp
  DO  i = 1, nr - 1
    wp = 1.0_dp
    DO  j = 1, nr - 1
      IF (j /= i) THEN
        wp = wp * (z-x(j))
      END IF
    END DO
    q = q + wp
  END DO
  gd = (hd - q*fd) / p
  z = z - fd / gd
  IF (it <= 40 .AND. ABS((z-z0)/z) > 1.0D-15) GO TO 10
  x(nr) = z
  x(n+1-nr) = -z
  r = 1.0_dp
  DO  k = 1, n
    r = 2.0_dp * r * k
  END DO
  w(nr) = 3.544907701811_dp * r / (hd*hd)
  w(n+1-nr) = w(nr)
END DO
IF (n /= 2*INT(n/2)) THEN
  r1 = 1.0_dp
  r2 = 1.0_dp
  DO  j = 1, n
    r1 = 2.0_dp * r1 * j
    IF (j >= (n+1)/2) r2 = r2 * j
  END DO
  w(n/2+1) = 0.88622692545276_dp * r1 / (r2*r2)
  x(n/2+1) = 0.0_dp
END IF
RETURN
END SUBROUTINE herzo
 
END MODULE herzo_func
 
 
 
PROGRAM mherzo
USE herzo_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:39

!       ===========================================================
!       Purpose : This program computes the zeros of Hermite
!                 polynomial Ln(x) in the interval [-ì,ì] and the
!                 corresponding weighting coefficients for Gauss-
!                 Hermite integration using subroutine HERZO
!       Input :   n    --- Order of the Hermite polynomial
!       Output:   X(n) --- Zeros of the Hermite polynomial
!                 W(n) --- Corresponding weighting coefficients
!       ===========================================================

REAL (dp)  :: x(100), w(100)
INTEGER    :: j, n

WRITE (*,*) 'Please enter the order of Hn(x), n '
READ (*,*) n
WRITE (*,5000) n
CALL herzo(n, x, w)
WRITE (*,*) '  Nodes and weights for Gauss-Hermite integration'
WRITE (*,*)
WRITE (*,*) '  i             xi                      Wi'
WRITE (*,*) ' -----------------------------------------------------'
DO  j = 1, n
  WRITE (*,5100) j, x(j), w(j)
END DO
STOP

5000 FORMAT (' n =', i3)
5100 FORMAT (' ', i3, '   ', g22.13, '   ', g22.13)
END PROGRAM mherzo
