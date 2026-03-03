MODULE clpn_func
 
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


SUBROUTINE clpn(n, x, y, cpn, cpd)

!       ==================================================
!       Purpose: Compute Legendre polynomials Pn(z) and
!                their derivatives Pn'(z) for a complex
!                argument
!       Input :  x --- Real part of z
!                y --- Imaginary part of z
!                n --- Degree of Pn(z), n = 0,1,2,...
!       Output:  CPN(n) --- Pn(z)
!                CPD(n) --- Pn'(z)
!       ==================================================

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: x
REAL (dp), INTENT(IN)      :: y
COMPLEX (dp), INTENT(OUT)  :: cpn(0:n)
COMPLEX (dp), INTENT(OUT)  :: cpd(0:n)

COMPLEX (dp)  :: cp0, cp1, cpf, z
INTEGER       :: k

z = CMPLX(x, y, KIND=dp)
cpn(0) = (1.0_dp,0.0_dp)
cpn(1) = z
cpd(0) = (0.0_dp,0.0_dp)
cpd(1) = (1.0_dp,0.0_dp)
cp0 = (1.0_dp,0.0_dp)
cp1 = z
DO  k = 2, n
  cpf = (2.0_dp*k-1.0_dp) / k * z * cp1 - (k-1.0_dp) / k * cp0
  cpn(k) = cpf
  IF (ABS(x) == 1.0_dp .AND. y == 0.0_dp) THEN
    cpd(k) = 0.5_dp * x ** (k+1) * k * (k+1)
  ELSE
    cpd(k) = k * (cp1 - z*cpf) / (1.0_dp - z*z)
  END IF
  cp0 = cp1
  cp1 = cpf
END DO
RETURN
END SUBROUTINE clpn
 
END MODULE clpn_func
 
 
 
PROGRAM mclpn
USE clpn_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:37

!       ==========================================================
!       Purpose: This program computes the Legendre polynomials
!                Pn(z) and Pn'(z) for a complex argument using
!                subroutine CLPN
!       Input :  x --- Real part of z
!                y --- Imaginary part of z
!                n --- Degree of Pn(z), n = 0,1,...,N
!       Output:  CPN(n) --- Pn(z)
!                CPD(n) --- Pn'(z)
!       Example: z = 3.0 +2.0 i

!       n    Re[Pn(z)]     Im[Pn(z)]     Re[Pn'(z)]   Im[Pn'(z)]
!      -----------------------------------------------------------
!       0   .100000D+01   .000000D+00   .000000D+00   .000000D+00
!       1   .300000D+01   .200000D+01   .100000D+01   .000000D+00
!       2   .700000D+01   .180000D+02   .900000D+01   .600000D+01
!       3  -.270000D+02   .112000D+03   .360000D+02   .900000D+02
!       4  -.539000D+03   .480000D+03  -.180000D+03   .790000D+03
!       5  -.461700D+04   .562000D+03  -.481500D+04   .441000D+04
!       ==========================================================

COMPLEX (dp)  :: cpn(0:100), cpd(0:100)
REAL (dp)     :: x, y
INTEGER       :: k, n

WRITE (*,*) '  Please enter Nmax, x and y (z=x+iy): '
READ (*,*) n, x, y
WRITE (*,5100) x, y
WRITE (*,*)
CALL clpn(n, x, y, cpn, cpd)
WRITE (*,*) '  n    Re[Pn(z)]     Im[Pn(z)]     Re[Pn''(Z)]   Im[Pn''(Z)]'
WRITE (*,*) ' -----------------------------------------------------------'
DO  k = 0, n
  WRITE (*,5000) k, cpn(k), cpd(k)
END DO
STOP

5000 FORMAT (' ', i3, 4g14.6)
5100 FORMAT ('   x =', f5.1, ',  y =', f5.1)
END PROGRAM mclpn
