MODULE e1z_func
 
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
 

SUBROUTINE e1z(z, ce1)

!       ====================================================
!       Purpose: Compute complex exponential integral E1(z)
!       Input :  z   --- Argument of E1(z)
!       Output:  CE1 --- E1(z)
!       ====================================================

COMPLEX (dp), INTENT(IN)   :: z
COMPLEX (dp), INTENT(OUT)  :: ce1

COMPLEX (dp)  :: cr, ct, ct0
REAL (dp)     :: a0, x
INTEGER       :: k
REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp, el = 0.5772156649015328_dp

x = REAL(z)
a0 = ABS(z)
IF (a0 == 0.0_dp) THEN
  ce1 = (1.0D+300,0.0_dp)
ELSE IF (a0 <= 10.0 .OR. x < 0.0 .AND. a0 < 20.0) THEN
  ce1 = (1.0_dp,0.0_dp)
  cr = (1.0_dp,0.0_dp)
  DO  k = 1, 150
    cr = -cr * k * z / (k+1) ** 2
    ce1 = ce1 + cr
    IF (ABS(cr) <= ABS(ce1)*1.0D-15) EXIT
  END DO
  ce1 = -el - LOG(z) + z * ce1
ELSE
  ct0 = (0.0_dp,0.0_dp)
  DO  k = 120, 1, -1
    ct0 = k / (1.0_dp + k/(z+ct0))
  END DO
  ct = 1.0_dp / (z+ct0)
  ce1 = EXP(-z) * ct
  IF (x <= 0.0 .AND. AIMAG(z) == 0.0) ce1 = ce1 - pi * (0.0_dp,1.0_dp)
END IF
RETURN
END SUBROUTINE e1z
 
END MODULE e1z_func
 
 
 
PROGRAM me1z
USE e1z_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:38

!       =========================================================
!       Purpose: This program computes the complex exponential
!                integral E1(z) using subroutine E1Z
!       Example:
!                     z            Re[E1(z)]       Im[E1(z)]
!                -----------------------------------------------
!                 3.0    2.0    -.90959209D-02   -.69001793D-02
!                 3.0   -2.0    -.90959209D-02    .69001793D-02
!                -3.0    2.0    -.28074891D+01    .59603353D+01
!                -3.0   -2.0    -.28074891D+01   -.59603353D+01
!                25.0   10.0    -.29302080D-12    .40391222D-12
!                25.0  -10.0    -.29302080D-12   -.40391222D-12
!               -25.0   10.0     .27279957D+10   -.49430610D+09
!               -25.0  -10.0     .27279957D+10    .49430610D+09
!       =========================================================

COMPLEX (dp)  :: z, ce1
REAL (dp)     :: x, y

WRITE (*,*) 'Please enter x and y ( z =x+iy ) '
READ (*,*) x, y
z = CMPLX(x, y, KIND=dp)
CALL e1z(z, ce1)
WRITE (*,*)
WRITE (*,*) '       z           Re[E1(z)]        Im[E1(z)]'
WRITE (*,*) ' -----------------------------------------------'
WRITE (*,5000) x, y, ce1
STOP

5000 FORMAT (' ', f5.1, '  ', f5.1, ' ', 2g17.8)
END PROGRAM me1z
