MODULE lpn_func
 
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


SUBROUTINE lpn(n, x, pn, pd)

!    ===============================================
!    Purpose: Compute Legendre polynomials Pn(x)
!             and their derivatives Pn'(x)
!    Input :  x --- Argument of Pn(x)
!             n --- Degree of Pn(x) ( n = 0,1,...)
!    Output:  PN(n) --- Pn(x)
!             PD(n) --- Pn'(x)
!    ===============================================

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: pn(0:n)
REAL (dp), INTENT(OUT)  :: pd(0:n)

REAL (dp)  :: p0, p1, pf
INTEGER    :: k

pn(0) = 1.0_dp
pn(1) = x
pd(0) = 0.0_dp
pd(1) = 1.0_dp
p0 = 1.0_dp
p1 = x
DO  k = 2, n
  pf = (2.0_dp*k-1.0_dp) / k * x * p1 - (k-1.0_dp) / k * p0
  pn(k) = pf
  IF (ABS(x) == 1.0_dp) THEN
    pd(k) = 0.5_dp * x ** (k+1) * k * (k+1)
  ELSE
    pd(k) = k * (p1 - x*pf) / (1.0_dp - x*x)
  END IF
  p0 = p1
  p1 = pf
END DO
RETURN
END SUBROUTINE lpn
 
END MODULE lpn_func
 
 
 
PROGRAM mlpn
USE lpn_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:43

!    ========================================================
!    Purpose: This program computes the Legendre polynomials
!             Pn(x) and their derivatives Pn'(x) using
!             subroutine LPN
!    Input :  x --- Argument of Pn(x)
!             n --- Degree of Pn(x) ( n = 0,1,...)
!    Output:  PN(n) --- Pn(x)
!             PD(n) --- Pn'(x)
!    Example:    x = 0.5
!               n          Pn(x)            Pn'(x)
!             ---------------------------------------
!               0       1.00000000        .00000000
!               1        .50000000       1.00000000
!               2       -.12500000       1.50000000
!               3       -.43750000        .37500000
!               4       -.28906250      -1.56250000
!               5        .08984375      -2.22656250
!    ========================================================

REAL (dp)  :: pn(0:100), pd(0:100), x
INTEGER    :: k, n

WRITE (*,*) '  Please enter Nmax and x '
READ (*,*) n, x
WRITE (*,5100) x
WRITE (*,*)
CALL lpn(n, x, pn, pd)
WRITE (*,*) '  n         Pn(x)           Pn''(X)'
WRITE (*,*) '---------------------------------------'
DO  k = 0, n
  WRITE (*,5000) k, pn(k), pd(k)
END DO
STOP

5000 FORMAT (' ', i3, 2g17.8)
5100 FORMAT ('   x =', f5.1)
END PROGRAM mlpn
