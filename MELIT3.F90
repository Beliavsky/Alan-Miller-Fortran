MODULE elit3_func
 
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
 


SUBROUTINE elit3(phi,hk,c,el3)

!       =========================================================
!       Purpose: Compute the elliptic integral of the third kind
!                using Gauss-Legendre quadrature
!       Input :  Phi --- Argument ( in degrees )
!                 k  --- Modulus   ( 0 ó k ó 1.0 )
!                 c  --- Parameter ( 0 ó c ó 1.0 )
!       Output:  EL3 --- Value of the elliptic integral of the
!                        third kind
!       =========================================================


REAL (dp), INTENT(IN)      :: phi
REAL (dp), INTENT(IN)      :: hk
REAL (dp), INTENT(IN)      :: c
REAL (dp), INTENT(OUT)     :: el3

REAL (dp)  :: c0, c1, c2, f1, f2, t1, t2
LOGICAL    :: lb1, lb2
INTEGER    :: i
REAL (dp), PARAMETER  :: t(10) = (/ .9931285991850949_dp, .9639719272779138_dp,  &
    .9122344282513259_dp, .8391169718222188_dp, .7463319064601508_dp,  &
    .6360536807265150_dp, .5108670019508271_dp, .3737060887154195_dp,  &
    .2277858511416451_dp, .7652652113349734D-1 /)
REAL (dp), PARAMETER  :: w(10) = (/ .1761400713915212D-1, .4060142980038694D-1,  &
    .6267204833410907D-1, .8327674157670475D-1, .1019301198172404_dp,  &
    .1181945319615184_dp, .1316886384491766_dp, .1420961093183820_dp,  &
    .1491729864726037_dp, .1527533871307258_dp /)

lb1 = hk == 1.0_dp .AND. ABS(phi-90.0) <= 1.0D-8
lb2 = c == 1.0_dp .AND. ABS(phi-90.0) <= 1.0D-8
IF (lb1 .OR. lb2) THEN
  el3 = 1.0D+300
  RETURN
END IF
c1 = 0.87266462599716D-2 * phi
c2 = c1
el3 = 0.0_dp
DO  i = 1, 10
  c0 = c2 * t(i)
  t1 = c1 + c0
  t2 = c1 - c0
  f1 = 1.0_dp / ((1.0_dp - c*SIN(t1)*SIN(t1))*SQRT(1.0_dp - hk*hk*SIN(t1)*SIN(t1)))
  f2 = 1.0_dp / ((1.0_dp - c*SIN(t2)*SIN(t2))*SQRT(1.0_dp - hk*hk*SIN(t2)*SIN(t2)))
  el3 = el3 + w(i) * (f1+f2)
END DO
el3 = c1 * el3
RETURN
END SUBROUTINE elit3
 
END MODULE elit3_func
 
 
 
PROGRAM melit3
USE elit3_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:38

!       ==========================================================
!       Purpose: This program computes the elliptic integral of
!                the third kind using subroutine ELIT3
!       Input :  Phi --- Argument ( in degrees )
!                 k  --- Modulus   ( 0 ó k ó 1 )
!                 c  --- Parameter ( 0 ó c ó 1 )
!       Output:  EL3 ÄÄÄ Value of the elliptic integral of the
!                        third kind
!       ==========================================================

REAL (dp)  :: c, el3, hk, phi

WRITE (*,*) 'Please enter phi, k and c '
READ (*,*) phi, hk, c
CALL elit3(phi, hk, c, el3)
WRITE (*,5000) el3
STOP

5000 FORMAT (' EL3=', f12.8)
END PROGRAM melit3
