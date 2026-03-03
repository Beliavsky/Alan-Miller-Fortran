MODULE bernob_func
 
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


SUBROUTINE bernob(n, bn)

!    ======================================
!    Purpose: Compute Bernoulli number Bn
!    Input :  n --- Serial number
!    Output:  BN(n) --- Bn
!    ======================================

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(OUT)  :: bn(0:n)

REAL (dp), PARAMETER  :: tpi = 6.283185307179586_dp
REAL (dp)  :: r1, r2, s
INTEGER    :: k, m

bn(0) = 1.0_dp
bn(1) = -0.5_dp
bn(2) = 1.0_dp / 6.0_dp
r1 = (2.0_dp/tpi) ** 2
DO  m = 4, n, 2
  r1 = -r1 * (m-1) * m / (tpi*tpi)
  r2 = 1.0_dp
  DO  k = 2, 10000
    s = (1.0_dp/k) ** m
    r2 = r2 + s
    IF (s < 1.0D-15) EXIT
  END DO
  bn(m) = r1 * r2
END DO
RETURN
END SUBROUTINE bernob
 
END MODULE bernob_func
 
 
 
PROGRAM mbernob
USE bernob_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:34

!    ===========================================================
!    Purpose: This program computes Bernoulli number Bn using
!             subroutine BERNOB
!    Example: Compute Bernouli number Bn for n = 0,1,...,10
!             Computed results:

!                n            Bn
!              --------------------------
!                0     .100000000000D+01
!                1    -.500000000000D+00
!                2     .166666666667D+00
!                4    -.333333333333D-01
!                6     .238095238095D-01
!                8    -.333333333333D-01
!               10     .757575757576D-01
!    ===========================================================

REAL (dp) :: b(0:200)
INTEGER   :: k, n

WRITE (*,*) '  Please enter Nmax: '
READ (*,*) n
CALL bernob(n,b)
WRITE (*,*) '   n            Bn'
WRITE (*,*) ' --------------------------'
WRITE (*,5000) 0, b(0)
WRITE (*,5000) 1, b(1)
DO  k = 2, n, 2
  WRITE (*,5000) k, b(k)
END DO
STOP

5000 FORMAT (t3, i3, g22.12)
END PROGRAM mbernob
