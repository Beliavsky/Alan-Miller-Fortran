MODULE lpni_func
 
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


SUBROUTINE lpni(n, x, pn, pd, pl)

!    =====================================================
!    Purpose: Compute Legendre polynomials Pn(x), Pn'(x)
!             and the integral of Pn(t) from 0 to x
!    Input :  x --- Argument of Pn(x)
!             n --- Degree of Pn(x) ( n = 0,1,... )
!    Output:  PN(n) --- Pn(x)
!             PD(n) --- Pn'(x)
!             PL(n) --- Integral of Pn(t) from 0 to x
!    =====================================================

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: pn(0:n)
REAL (dp), INTENT(OUT)  :: pd(0:n)
REAL (dp), INTENT(OUT)  :: pl(0:n)

REAL (dp)  :: p0, p1, pf, r
INTEGER    :: j, k, n1

pn(0) = 1.0_dp
pn(1) = x
pd(0) = 0.0_dp
pd(1) = 1.0_dp
pl(0) = x
pl(1) = 0.5_dp * x * x
p0 = 1.0_dp
p1 = x
DO  k = 2, n
  pf = (2.0_dp*k-1.0_dp) / k * x * p1 - (k-1.0_dp) / k * p0
  pn(k) = pf
  IF (ABS(x) == 1.0_dp) THEN
    pd(k) = 0.5_dp * x ** (k+1) * k * (k+1)
  ELSE
    pd(k) = k * (p1 - x*pf) / (1.0_dp-x*x)
  END IF
  pl(k) = (x*pn(k) - pn(k-1)) / (k+1)
  p0 = p1
  p1 = pf
  IF (k /= 2*INT(k/2)) THEN
    r = 1.0_dp / (k+1)
    n1 = (k-1) / 2
    DO  j = 1, n1
      r = (0.5_dp/j - 1.0_dp) * r
    END DO
    pl(k) = pl(k) + r
  END IF
END DO
RETURN
END SUBROUTINE lpni
 
END MODULE lpni_func
 
 
 
PROGRAM mlpni
USE lpni_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:43

!    ========================================================
!    Purpose: This program computes the Legendre polynomials
!             Pn(x), Pn'(x) and the integral of Pn(t) from 0
!             to x using subroutine LPNI
!    Input :  x --- Argument of Pn(x)
!             n --- Degree of Pn(x) ( n = 0,1,... )
!    Output:  PN(n) --- Pn(x)
!             PD(n) --- Pn'(x)
!             PL(n) --- Integral of Pn(t) from 0 to x
!    Example: x = 0.50
!             n       Pn(x)         Pn'(x)        Pn(t)dt
!            ---------------------------------------------
!             0    1.00000000     .00000000     .50000000
!             1     .50000000    1.00000000     .12500000
!             2    -.12500000    1.50000000    -.18750000
!             3    -.43750000     .37500000    -.14843750
!             4    -.28906250   -1.56250000     .05859375
!             5     .08984375   -2.22656250     .11816406
!    ========================================================

REAL (dp)  :: pn(0:100), pd(0:100), pl(0:100), x
INTEGER    :: k, n

WRITE (*,*) '  Please enter Nmax and x '
READ (*,*) n, x
WRITE (*,5100) x
WRITE (*,*)
WRITE (*,*) '  n        Pn(x)          Pn''(X)         Pn(t)dt'
WRITE (*,*) ' ---------------------------------------------------'
CALL lpni(n, x, pn, pd, pl)
DO  k = 0, n
  WRITE (*,5000) k, pn(k), pd(k), pl(k)
END DO
STOP

5000 FORMAT (' ', i3, 3g16.8)
5100 FORMAT ('   x =', f5.2)
END PROGRAM mlpni
