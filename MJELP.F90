MODULE jelp_func
 
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


SUBROUTINE jelp(u, hk, esn, ecn, edn, eph)

!    ========================================================
!    Purpose: Compute Jacobian elliptic functions sn u, cn u
!             and dn u
!    Input  : u   --- Argument of Jacobian elliptic fuctions
!             Hk  --- Modulus k ( 0 ó k ó 1 )
!    Output : ESN --- sn u
!             ECN --- cn u
!             EDN --- dn u
!             EPH --- phi ( in degrees )
!    ========================================================

REAL (dp), INTENT(IN)   :: u
REAL (dp), INTENT(IN)   :: hk
REAL (dp), INTENT(OUT)  :: esn
REAL (dp), INTENT(OUT)  :: ecn
REAL (dp), INTENT(OUT)  :: edn
REAL (dp), INTENT(OUT)  :: eph

REAL (dp)  :: a, a0, b, b0, c, d, dn, r(40), sa, t
INTEGER    :: j, n
REAL (dp), PARAMETER  :: pi = 3.14159265358979_dp

a0 = 1.0D0
b0 = SQRT(1.0D0 - hk*hk)
DO  n = 1, 40
  a = (a0+b0) / 2.0D0
  b = SQRT(a0*b0)
  c = (a0-b0) / 2.0D0
  r(n) = c / a
  IF (c < 1.0D-7) EXIT
  a0 = a
  b0 = b
END DO
dn = 2.0D0 ** n * a * u
DO  j = n, 1, -1
  t = r(j) * SIN(dn)
  sa = ATAN(t/SQRT(ABS(1.0D0 - t*t)))
  d = .5D0 * (dn+sa)
  dn = d
END DO
eph = d * 180.0D0 / pi
esn = SIN(d)
ecn = COS(d)
edn = SQRT(1.0D0 - hk*hk*esn*esn)
RETURN
END SUBROUTINE jelp
 
END MODULE jelp_func
 
 
 
PROGRAM mjelp
USE jelp_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:41

!    ============================================================
!    Purpose: This program computes Jacobian elliptic functions
!             sn u, cn u and dn u using subroutine JELP
!    Input  : u   --- Argument of Jacobian elliptic fuctions
!             Hk  --- Modulus k ( 0 ó k ó 1 )
!    Output : ESN --- sn u
!             ECN --- cn u
!             EDN --- dn u
!             EPH --- phi ( in degrees )
!    Example:
!             k = .5, ( K(k) = 1.68575035 ), and u = u0*K

!             u0       phi       sn u        cn u        dn u
!           ----------------------------------------------------
!            0.0      .0000    .0000000   1.0000000   1.0000000
!            0.5    47.0586    .7320508    .6812500    .9306049
!            1.0    90.0000   1.0000000    .0000000    .8660254
!            1.5   132.9414    .7320508   -.6812500    .9306049
!            2.0   180.0000    .0000000  -1.0000000   1.0000000
!    ============================================================

REAL (dp)  :: eph, esn, ecn, edn, hk, u

WRITE (*,*) 'Please enter k and u '
READ (*,*) hk, u
WRITE (*,*)
WRITE (*,*) '   k        u          phi        sn u        cn u        dn u'
WRITE (*,*) ' ----------------------------------------------------------------'
CALL jelp(u, hk, esn, ecn, edn, eph)
WRITE (*,5000) hk, u, eph, esn, ecn, edn
STOP

5000 FORMAT (' ', f5.3, f12.7, '  ', f9.5, 3F12.7)
END PROGRAM mjelp
