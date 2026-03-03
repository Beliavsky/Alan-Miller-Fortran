MODULE elit_func
 
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
 

SUBROUTINE elit(hk, phi, fe, ee)

!       ==================================================
!       Purpose: Compute complete and incomplete elliptic
!                integrals F(k,phi) and E(k,phi)
!       Input  : HK  --- Modulus k ( 0 ó k ó 1 )
!                Phi --- Argument ( in degrees )
!       Output : FE  --- F(k,phi)
!                EE  --- E(k,phi)
!       ==================================================

REAL (dp), INTENT(IN)      :: hk
REAL (dp), INTENT(IN)      :: phi
REAL (dp), INTENT(OUT)     :: fe
REAL (dp), INTENT(OUT)     :: ee

REAL (dp), PARAMETER  :: pi = 3.14159265358979_dp
REAL (dp)  :: a, a0, b, b0, c, ce, ck, d, d0, fac, g, r
INTEGER    :: n

g = 0.0_dp
a0 = 1.0_dp
b0 = SQRT(1.0_dp - hk*hk)
d0 = (pi/180.0_dp) * phi
r = hk * hk
IF (hk == 1.0_dp .AND. phi == 90.0_dp) THEN
  fe = 1.0D+300
  ee = 1.0_dp
ELSE IF (hk == 1.0_dp) THEN
  fe = LOG((1.0_dp + SIN(d0))/COS(d0))
  ee = SIN(d0)
ELSE
  fac = 1.0_dp
  DO  n = 1, 40
    a = (a0+b0) / 2.0_dp
    b = SQRT(a0*b0)
    c = (a0-b0) / 2.0_dp
    fac = 2.0_dp * fac
    r = r + fac * c * c
    IF (phi /= 90.0_dp) THEN
      d = d0 + ATAN((b0/a0)*TAN(d0))
      g = g + c * SIN(d)
      d0 = d + pi * INT(d/pi + .5_dp)
    END IF
    a0 = a
    b0 = b
    IF (c < 1.0D-7) EXIT
  END DO

  ck = pi / (2.0_dp*a)
  ce = pi * (2.0_dp - r) / (4.0_dp*a)
  IF (phi == 90.0_dp) THEN
    fe = ck
    ee = ce
  ELSE
    fe = d / (fac*a)
    ee = fe * ce / ck + g
  END IF
END IF
RETURN
END SUBROUTINE elit
 
END MODULE elit_func
 
 
 
PROGRAM melit
USE elit_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:38

!       ==========================================================
!       Purpose: This program computes complete and incomplete
!                elliptic integrals F(k,phi) and E(k,phi) using
!                subroutine ELIT
!       Input  : HK  --- Modulus k ( 0 ó k ó 1 )
!                Phi --- Argument ( in degrees )
!       Output : FE  --- F(k,phi)
!                EE  --- E(k,phi)
!       Example:
!                k = .5

!                 phi     F(k,phi)       E(k,phi)
!                -----------------------------------
!                   0      .00000000      .00000000
!                  15      .26254249      .26106005
!                  30      .52942863      .51788193
!                  45      .80436610      .76719599
!                  60     1.08955067     1.00755556
!                  75     1.38457455     1.23988858
!                  90     1.68575035     1.46746221
!       ==========================================================

REAL (dp)  :: ee, fe, hk, phi

WRITE (*,*) 'Please enter k and phi (in degs.) '
READ (*,*) hk, phi
WRITE (*,*)
WRITE (*,*) '  phi     F(k,phi)       E(k,phi)'
WRITE (*,*) ' -----------------------------------'
CALL elit(hk, phi, fe, ee)
WRITE (*,5000) phi, fe, ee
STOP

5000 FORMAT (' ', f5.0, 2F15.8)
END PROGRAM melit
