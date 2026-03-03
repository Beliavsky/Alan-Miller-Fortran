MODULE cerror_func
 
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
 

SUBROUTINE cerror(z, cer)

!    ====================================================
!    Purpose: Compute error function erf(z) for a complex
!             argument (z=x+iy)
!    Input :  z   --- Complex argument
!    Output:  CER --- erf(z)
!    ====================================================

COMPLEX (dp), INTENT(IN)   :: z
COMPLEX (dp), INTENT(OUT)  :: cer

REAL (dp)     :: a0
COMPLEX (dp)  :: c0, cl, cr, cs, z1
INTEGER       :: k
REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp

a0 = ABS(z)
c0 = EXP(-z*z)
z1 = z
IF (REAL(z, KIND=dp) < 0.0_dp) z1 = -z
IF (a0 <= 5.8_dp) THEN
  cs = z1
  cr = z1
  DO  k = 1, 120
    cr = cr * z1 * z1 / (k + 0.5_dp)
    cs = cs + cr
    IF (ABS(cr/cs) < 1.0D-15) EXIT
  END DO
  cer = 2.0_dp * c0 * cs / SQRT(pi)
ELSE
  cl = 1.0_dp / z1
  cr = cl
  DO  k = 1, 13
    cr = -cr * (k - 0.5_dp) / (z1*z1)
    cl = cl + cr
    IF (ABS(cr/cl) < 1.0D-15) EXIT
  END DO
  cer = 1.0_dp - c0 * cl / SQRT(pi)
END IF
IF (REAL(z, KIND=dp) < 0.0_dp) cer = -cer
RETURN
END SUBROUTINE cerror

END MODULE cerror_func
 
 
 
PROGRAM mcerror
USE cerror_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:34

!    ============================================================
!    Purpose: This program computes the error function erf(z)
!             for a complex argument using subroutine CERROR
!    Input :  x   --- Real part of z
!             y   --- Imaginary part of z  ( y ó 3.0 )
!    Output:  ERR --- Real part of erf(z)
!             ERI --- Imaginary part of erf(z)
!    Example:
!                x       y       Re[erf(z)]      Im[erf(z)]
!              ---------------------------------------------
!               1.0     2.0      -.53664357     -5.04914370
!               2.0     2.0      1.15131087       .12729163
!               3.0     2.0       .99896328      -.00001155
!               4.0     2.0      1.00000057      -.00000051
!               5.0     2.0      1.00000000       .00000000
!    ============================================================

COMPLEX (dp)  :: cer, z
REAL (dp)     :: x, y

DO
  WRITE (*,*) 'X,Y=? '
  READ (*,*) x, y
  WRITE (*,*) '   x      y      Re[erf(z)]      Im[erf(z)]'
  WRITE (*,*) ' ---------------------------------------------'
  z = CMPLX(x, y, KIND=dp)
  CALL cerror(z, cer)
  WRITE (*,5000) z, cer
END DO
STOP

5000 FORMAT (' ', f5.1, '  ', f5.1, ' ', 2g16.8)
END PROGRAM mcerror
