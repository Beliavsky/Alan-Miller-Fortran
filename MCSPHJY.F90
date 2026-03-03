MODULE csphjy_func
 
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


SUBROUTINE csphjy(n, z, nm, csj, cdj, csy, cdy)

!       ==========================================================
!       Purpose: Compute spherical Bessel functions jn(z) & yn(z)
!                and their derivatives for a complex argument
!       Input :  z --- Complex argument
!                n --- Order of jn(z) & yn(z) ( n = 0,1,2,... )
!       Output:  CSJ(n) --- jn(z)
!                CDJ(n) --- jn'(z)
!                CSY(n) --- yn(z)
!                CDY(n) --- yn'(z)
!                NM --- Highest order computed
!       Routines called:
!                MSTA1 and MSTA2 for computing the starting
!                point for backward recurrence
!       ==========================================================

INTEGER, INTENT(IN)        :: n
COMPLEX (dp), INTENT(IN)   :: z
INTEGER, INTENT(OUT)       :: nm
COMPLEX (dp), INTENT(OUT)  :: csj(0:n)
COMPLEX (dp), INTENT(OUT)  :: cdj(0:n)
COMPLEX (dp), INTENT(OUT)  :: csy(0:n)
COMPLEX (dp), INTENT(OUT)  :: cdy(0:n)

REAL (dp)     :: a0
COMPLEX (dp)  :: cf, cf0, cf1, cs, csa, csb
INTEGER       :: k, m

a0 = ABS(z)
nm = n
IF (a0 < 1.0D-60) THEN
  DO  k = 0, n
    csj(k) = 0.0_dp
    cdj(k) = 0.0_dp
    csy(k) = -1.0D+300
    cdy(k) = 1.0D+300
  END DO
  csj(0) = (1.0_dp,0.0_dp)
  cdj(1) = (.333333333333333_dp,0.0_dp)
  RETURN
END IF
csj(0) = SIN(z) / z
csj(1) = (csj(0) - COS(z)) / z
IF (n >= 2) THEN
  csa = csj(0)
  csb = csj(1)
  m = msta1(a0, 200)
  IF (m < n) THEN
    nm = m
  ELSE
    m = msta2(a0, n, 15)
  END IF
  cf0 = 0.0_dp
  cf1 = 1.0_dp - 100
  DO  k = m, 0, -1
    cf = (2*k+3) * cf1 / z - cf0
    IF (k <= nm) csj(k) = cf
    cf0 = cf1
    cf1 = cf
  END DO
  IF (ABS(csa) > ABS(csb)) cs = csa / cf
  IF (ABS(csa) <= ABS(csb)) cs = csb / cf0
  DO  k = 0, nm
    csj(k) = cs * csj(k)
  END DO
END IF
cdj(0) = (COS(z) - SIN(z)/z) / z
DO  k = 1, nm
  cdj(k) = csj(k-1) - (k+1) * csj(k) / z
END DO
csy(0) = -COS(z) / z
csy(1) = (csy(0) - SIN(z)) / z
cdy(0) = (SIN(z) + COS(z)/z) / z
cdy(1) = (2.0_dp*cdy(0) - COS(z)) / z
DO  k = 2, nm
  IF (ABS(csj(k-1)) > ABS(csj(k-2))) THEN
    csy(k) = (csj(k)*csy(k-1) - 1.0_dp/(z*z)) / csj(k-1)
  ELSE
    csy(k) = (csj(k)*csy(k-2) - (2*k-1)/z**3) / csj(k-2)
  END IF
END DO
DO  k = 2, nm
  cdy(k) = csy(k-1) - (k+1) * csy(k) / z
END DO
RETURN
END SUBROUTINE csphjy



FUNCTION msta1(x, mp) RESULT(fn_val)

!       ===================================================
!       Purpose: Determine the starting point for backward
!                recurrence such that the magnitude of
!                Jn(x) at that point is about 10^(-MP)
!       Input :  x     --- Argument of Jn(x)
!                MP    --- Value of magnitude
!       Output:  MSTA1 --- Starting point
!       ===================================================

REAL (dp), INTENT(IN)      :: x
INTEGER, INTENT(IN)        :: mp
INTEGER                    :: fn_val

REAL (dp)  :: a0, f, f0, f1
INTEGER    :: it, n0, n1, nn

a0 = ABS(x)
n0 = INT(1.1*a0) + 1
f0 = envj(n0,a0) - mp
n1 = n0 + 5
f1 = envj(n1,a0) - mp
DO  it = 1, 20
  nn = n1 - (n1-n0) / (1.0_dp - f0/f1)
  f = envj(nn,a0) - mp
  IF (ABS(nn-n1) < 1) EXIT
  n0 = n1
  f0 = f1
  n1 = nn
  f1 = f
END DO

fn_val = nn
RETURN
END FUNCTION msta1



FUNCTION msta2(x, n, mp) RESULT(fn_val)

!       ===================================================
!       Purpose: Determine the starting point for backward
!                recurrence such that all Jn(x) has MP
!                significant digits
!       Input :  x  --- Argument of Jn(x)
!                n  --- Order of Jn(x)
!                MP --- Significant digit
!       Output:  MSTA2 --- Starting point
!       ===================================================

REAL (dp), INTENT(IN)      :: x
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: mp
INTEGER                    :: fn_val

REAL (dp)  :: a0, ejn, f, f0, f1, hmp, obj
INTEGER    :: it, n0, n1, nn

a0 = ABS(x)
hmp = 0.5_dp * mp
ejn = envj(n, a0)
IF (ejn <= hmp) THEN
  obj = mp
  n0 = INT(1.1*a0)
ELSE
  obj = hmp + ejn
  n0 = n
END IF
f0 = envj(n0,a0) - obj
n1 = n0 + 5
f1 = envj(n1,a0) - obj
DO  it = 1, 20
  nn = n1 - (n1-n0) / (1.0_dp - f0/f1)
  f = envj(nn, a0) - obj
  IF (ABS(nn-n1) < 1) EXIT
  n0 = n1
  f0 = f1
  n1 = nn
  f1 = f
END DO

fn_val = nn + 10
RETURN
END FUNCTION msta2



FUNCTION envj(n, x) RESULT(fn_val)

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: x
REAL (dp)                  :: fn_val

fn_val = 0.5_dp * LOG10(6.28_dp*n) - n * LOG10(1.36_dp*x/n)
RETURN
END FUNCTION envj

END MODULE csphjy_func
 
 
 
PROGRAM mcsphjy
USE csphjy_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:37

!       ================================================================
!       Purpose: This program computes the spherical Bessel functions
!                jn(z), yn(z), and their derivatives for a complex
!                argument using subroutine CSPHJY
!       Input :  z --- Complex argument
!                n --- Order of jn(z) & yn(z) ( 0 ó n ó 250 )
!       Output:  CSJ(n) --- jn(z)
!                CDJ(n) --- jn'(z)
!                CSY(n) --- yn(z)
!                CDY(n) --- yn'(z)
!       Example: z = 4.0 +i 2.0

!     n     Re[jn(z)]       Im[jn(z)]       Re[jn'(z)]      Im[jn'(z)]
!   --------------------------------------------------------------------
!     0  -.80651523D+00  -.18941093D+00  -.37101203D-01   .75210758D+00
!     1   .37101203D-01  -.75210758D+00  -.67093420D+00   .11885235D+00
!     2   .60314368D+00  -.27298399D+00  -.24288981D+00  -.40737409D+00
!     3   .42955048D+00   .17755176D+00   .18848259D+00  -.24320520D+00
!     4   .12251323D+00   .22087111D+00   .19660170D+00   .17937264D-01
!     5  -.10242676D-01   .10975433D+00   .68951842D-01   .83020305D-01

!     n     Re[yn(z)]       Im[yn(z)]       Re[yn'(z)]      Im[yn'(z)]
!   --------------------------------------------------------------------
!     0   .21734534D+00  -.79487692D+00  -.77049661D+00  -.87010064D-02
!     1   .77049661D+00   .87010064D-02  -.92593503D-01  -.64425800D+00
!     2   .24756293D+00   .56894854D+00   .45127429D+00  -.25839924D+00
!     3  -.23845941D+00   .43646607D+00   .26374403D+00   .12439192D+00
!     4  -.27587985D+00   .20902555D+00  -.67092335D-01   .89500599D-01
!     5  -.70001327D-01   .18807178D+00  -.30472133D+00  -.58661384D-01
!       ================================================================

REAL (dp)     :: x, y
COMPLEX (dp)  :: csj(0:250), cdj(0:250), csy(0:250), cdy(0:250), z
INTEGER       :: k, n, nm, ns

WRITE (*,*) 'Please enter n,x,y (z=x+iy) '
READ (*,*) n, x, y
WRITE (*,5100) n, x, y
z = CMPLX(x, y, KIND=dp)
IF (n <= 8) THEN
  ns = 1
ELSE
  WRITE (*,*) 'Please enter order step Ns '
  READ (*,*) ns
END IF
CALL csphjy(n, z, nm, csj, cdj, csy, cdy)
WRITE(*,*)
WRITE(*,*) '  n      Re[jn(z)]        Im[jn(z)]        Re[jn''(Z)]       IM[JN''(Z)]'
WRITE(*,*) '------------------------------------------------------------------------'
DO  k = 0, nm, ns
  WRITE (*,5000) k, csj(k), cdj(k)
END DO
WRITE(*,*)
WRITE(*,*) '  n      Re[yn(z)]        Im[yn(z)]        Re[yn''(Z)]       IM[YN''(Z)]'
WRITE(*,*) '------------------------------------------------------------------------'
DO  k = 0, nm, ns
  WRITE (*,5000) k, csy(k), cdy(k)
END DO
STOP

5000 FORMAT (' ', i3, 4g17.8)
5100 FORMAT ('   Nmaz =', i3, ',     z = ', f8.1, '+ i', f8.1)
END PROGRAM mcsphjy
