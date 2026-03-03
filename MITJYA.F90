MODULE itjya_func
 
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


SUBROUTINE itjya(x, tj, ty)

!    ==========================================================
!    Purpose: Integrate Bessel functions J0(t) & Y0(t) with
!             respect to t from 0 to x
!    Input :  x  --- Upper limit of the integral ( x ò 0 )
!    Output:  TJ --- Integration of J0(t) from 0 to x
!             TY --- Integration of Y0(t) from 0 to x
!    =======================================================

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: tj
REAL (dp), INTENT(OUT)  :: ty

REAL (dp)  :: a(18), a0, a1, af, bf, bg, eps, r, r2, rc, rs, ty1, ty2, x2, xp
INTEGER    :: k
REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp, el = .5772156649015329_dp

eps = 1.0D-12
IF (x == 0.0_dp) THEN
  tj = 0.0_dp
  ty = 0.0_dp
ELSE IF (x <= 20.0_dp) THEN
  x2 = x * x
  tj = x
  r = x
  DO  k = 1, 60
    r = -.25_dp * r * (2*k-1.0_dp) / (2*k+1.0_dp) / (k*k) * x2
    tj = tj + r
    IF (ABS(r) < ABS(tj)*eps) EXIT
  END DO
  ty1 = (el + LOG(x/2.0_dp)) * tj
  rs = 0.0_dp
  ty2 = 1.0_dp
  r = 1.0_dp
  DO  k = 1, 60
    r = -.25_dp * r * (2*k-1.0_dp) / (2*k+1.0_dp) / (k*k) * x2
    rs = rs + 1.0_dp / k
    r2 = r * (rs + 1.0_dp / (2*k+1))
    ty2 = ty2 + r2
    IF (ABS(r2) < ABS(ty2)*eps) EXIT
  END DO
  ty = (ty1 - x*ty2) * 2.0_dp / pi
ELSE
  a0 = 1.0_dp
  a1 = 5.0_dp / 8.0_dp
  a(1) = a1
  DO  k = 1, 16
    af = ((1.5_dp*(k + .5_dp)*(k + 5.0_dp/6.0_dp)*a1 - .5_dp*(k + .5_dp)*(k +  &
         .5_dp)*(k - .5_dp)*a0)) / (k+1)
    a(k+1) = af
    a0 = a1
    a1 = af
  END DO
  bf = 1.0_dp
  r = 1.0_dp
  DO  k = 1, 8
    r = -r / (x*x)
    bf = bf + a(2*k) * r
  END DO
  bg = a(1) / x
  r = 1.0_dp / x
  DO  k = 1, 8
    r = -r / (x*x)
    bg = bg + a(2*k+1) * r
  END DO
  xp = x + .25_dp * pi
  rc = SQRT(2.0_dp/(pi*x))
  tj = 1.0_dp - rc * (bf*COS(xp) + bg*SIN(xp))
  ty = rc * (bg*COS(xp) - bf*SIN(xp))
END IF
RETURN
END SUBROUTINE itjya
 
END MODULE itjya_func
 
 
 
PROGRAM mitjya
USE itjya_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:41

!       ===========================================================
!       Purpose: This program evaluates the integral of Bessel
!                functions J0(t) and Y0(t) with respect to t
!                from 0 to x using subroutine ITJYA
!       Input :  x  --- Upper limit of the integral ( x ò 0 )
!       Output:  TJ --- Integration of J0(t) from 0 to x
!                TY --- Integration of Y0(t) from 0 to x
!       Example:
!                   x         J0(t)dt          Y0(t)dt
!                ---------------------------------------
!                  5.0       .71531192       .19971938
!                 10.0      1.06701130       .24129032
!                 15.0      1.20516194       .00745772
!                 20.0      1.05837882      -.16821598
!                 25.0       .87101492      -.09360793
!                 30.0       .88424909       .08822971
!       ===========================================================

REAL (dp)  :: x, tj, ty

WRITE (*,*) 'Please enter x '
READ (*,*) x
WRITE (*,*) '   x         J0(t)dt          Y0(t)dt'
WRITE (*,*) '---------------------------------------'
CALL itjya(x, tj, ty)
WRITE (*,5000) x, tj, ty
STOP

5000 FORMAT (' ', f5.1, 2F16.8)
END PROGRAM mitjya
