MODULE legzo_func
 
! From the book "Computation of Special Functions"
!      by Shanjie Zhang and Jianming Jin
!   Copyright 1996 by John Wiley & Sons, Inc.
! The authors state:
!   "However, we give permission to the reader who purchases this book
!    to incorporate any of these programs into his or her programs
!    provided that the copyright is acknowledged."

! Zhang & Jin's errata applied to the upper limits of 2 DO loops.
! Latest revision - 5 January 2002
! Corrections by Alan Miller (amiller @ bigpond.net.au)
! Elements of array X were used before being defined.

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
 
CONTAINS
 

SUBROUTINE legzo(n, x, w)

!       =========================================================
!       Purpose : Compute the zeros of Legendre polynomial Pn(x)
!                 in the interval [-1,1], and the corresponding
!                 weighting coefficients for Gauss-Legendre
!                 integration
!       Input :   n    --- Order of the Legendre polynomial
!       Output:   X(n) --- Zeros of the Legendre polynomial
!                 W(n) --- Corresponding weighting coefficients
!       =========================================================

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(OUT)     :: x(n)
REAL (dp), INTENT(OUT)     :: w(n)

REAL (dp)  :: f0, f1, fd, gd, p, pd, pf, q, wp, z, z0
INTEGER    :: i, j, k, n0, nr

n0 = (n+1) / 2
x(1:n0) = 0.0_dp
DO  nr = 1, n0
  z = COS(3.1415926_dp*(nr-0.25_dp)/n)

  10 z0 = z
  p = 1.0_dp
  DO  i = 1, nr - 1
    p = p * (z-x(i))
  END DO
  f0 = 1.0_dp
  IF (nr == n0 .AND. n /= 2*INT(n/2)) z = 0.0_dp
  f1 = z
  DO  k = 2, n
    pf = (2.0_dp - 1.0_dp/k) * z * f1 - (1.0_dp - 1.0_dp/k) * f0
    pd = k * (f1 - z*pf) / (1.0_dp - z*z)
    f0 = f1
    f1 = pf
  END DO
  IF (z /= 0.0) THEN
    fd = pf / p
    q = 0.0_dp
    DO  i = 1, nr
      wp = 1.0_dp
      DO  j = 1, nr
        IF (j /= i) wp = wp * (z-x(j))
      END DO
      q = q + wp
    END DO
    gd = (pd - q*fd) / p
    z = z - fd / gd
    IF (ABS(z-z0) > ABS(z)*1.0D-15) GO TO 10
  END IF
  x(nr) = z
  x(n+1-nr) = -z
  w(nr) = 2.0_dp / ((1.0_dp - z*z)*pd*pd)
  w(n+1-nr) = w(nr)
END DO
RETURN
END SUBROUTINE legzo
 
END MODULE legzo_func
 
 
 
PROGRAM mlegzo
USE legzo_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:43

!       ============================================================
!       Purpose : This program computes the zeros of Legendre
!                 polynomial Pn(x) in the interval [-1,1] and the
!                 corresponding weighting coefficients for Gauss-
!                 Legendre integration using subroutine LEGZO
!       Input :   n    --- Order of the Legendre polynomial
!       Output:   X(n) --- Zeros of the Legendre polynomial
!                 W(n) --- Corresponding weighting coefficients
!       ============================================================

REAL (dp)  :: x(120), w(120)
INTEGER    :: i, n

WRITE (*,*) 'Please enter the order of Pn(x), n '
READ (*,*) n
WRITE (*,5000) n
CALL legzo(n, x, w)
WRITE (*,*) '  Nodes and weights for Gauss-Legendre integration'
WRITE (*,*)
WRITE (*,*) '  i              xi                   Wi'
WRITE (*,*) ' ------------------------------------------------'
DO  i = 1, n
  WRITE (*,5100) i, x(i), w(i)
END DO
STOP

5000 FORMAT (' n =', i3)
5100 FORMAT (' ', i3, ' ', f22.13, g22.13)
END PROGRAM mlegzo
