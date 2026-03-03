MODULE cjk_func
 
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
 

SUBROUTINE cjk(km, a)

!       ==============================================================
!       Purpose: Compute the expansion coefficients for the asymptotic
!                expansion of Bessel functions with large orders
!       Input :  Km   --- Maximum k
!       Output:  A(L) --- Cj(k) where j and k are related to L
!                         by L=j+1+[k*(k+1)]/2; j,k=0,1,...,Km
!       ==============================================================

INTEGER, INTENT(IN)        :: km
REAL (dp), INTENT(OUT)     :: a(:)

REAL (dp)  :: f, f0, g, g0
INTEGER    :: j, k, l1, l2, l3, l4

a(1) = 1.0_dp
f0 = 1.0_dp
g0 = 1.0_dp
DO  k = 0, km - 1
  l1 = (k+1) * (k+2) / 2 + 1
  l2 = (k+1) * (k+2) / 2 + k + 2
  f = (0.5_dp*k + 0.125_dp/(k+1)) * f0
  g = -(1.5_dp*k + 0.625_dp/(3.0*(k+1))) * g0
  a(l1) = f
  a(l2) = g
  f0 = f
  g0 = g
END DO
DO  k = 1, km - 1
  DO  j = 1, k
    l3 = k * (k+1) / 2 + j + 1
    l4 = (k+1) * (k+2) / 2 + j + 1
    a(l4) = (j + 0.5_dp*k + 0.125_dp/(2.0*j+k+1.0)) * a(l3) - (j + 0.5_dp*k  &
            - 1.0 + 0.625_dp/(2*j+k+1)) * a(l3-1)
  END DO
END DO
RETURN
END SUBROUTINE cjk
 
END MODULE cjk_func
 
 
 
PROGRAM mcjk
USE cjk_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:36

!       ============================================================
!       Purpose: This program computes the expansion coefficients
!                for the asymptotic expansion of Bessel functions
!                with large orders using subroutine CJK
!       Input :  Km   --- Maximum k
!       Output:  A(L) --- Cj(k) where j and k are related to L by
!                         L=j+1+[k*(k+1)]/2; j,k=0,1,2,...,Km
!       ============================================================

REAL (dp)  :: a(231)
INTEGER    :: k, km, lm

WRITE (*,*) 'Please enter Km ( ó 20 ): '
READ (*,*) km
lm = km + 1 + (km*(km+1)) / 2
CALL cjk(km, a)
DO  k = 1, lm
  WRITE (*,5000) k, a(k)
END DO
STOP

5000 FORMAT (' ', i3, g25.14)
END PROGRAM mcjk
