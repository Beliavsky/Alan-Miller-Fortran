MODULE cpsi_func
 
! From the book "Computation of Special Functions"
!      by Shanjie Zhang and Jianming Jin
!   Copyright 1996 by John Wiley & Sons, Inc.
! The authors state:
!   "However, we give permission to the reader who purchases this book
!    to incorporate any of these programs into his or her programs
!    provided that the copyright is acknowledged."

! Latest revision - 1 January 2002
! Corrections by Alan Miller (amiller @ bigpond.net.au)
! x1 & y1 undefined for positive x.

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
 
CONTAINS


SUBROUTINE cpsi(x, y, psr, psi)

!       =============================================
!       Purpose: Compute the psi function for a
!                complex argument
!       Input :  x   --- Real part of z
!                y   --- Imaginary part of z
!       Output:  PSR --- Real part of psi(z)
!                PSI --- Imaginary part of psi(z)
!       =============================================

REAL (dp), INTENT(IN OUT)  :: x
REAL (dp), INTENT(IN OUT)  :: y
REAL (dp), INTENT(OUT)     :: psr
REAL (dp), INTENT(OUT)     :: psi

REAL (dp), PARAMETER  :: a(8) = (/ -.8333333333333D-01, .83333333333333333D-02,  &
    -.39682539682539683D-02, .41666666666666667D-02, -.75757575757575758D-02, &
     .21092796092796093D-01, -.83333333333333333D-01, .4432598039215686_dp /)
REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp
REAL (dp)  :: ct2, ri, rr, th, tm, tn, x0, x1, y1, z0, z2
INTEGER    :: k, n

IF (y == 0.0_dp .AND. x == INT(x) .AND. x <= 0.0_dp) THEN
  psr = 1.0D+300
  psi = 0.0_dp
ELSE
  x1 = x
  y1 = y
  IF (x < 0.0_dp) THEN
    x = -x
    y = -y
  END IF
  x0 = x
  IF (x < 8.0_dp) THEN
    n = 8 - INT(x)
    x0 = x + n
  END IF
  IF (x0 == 0.0_dp .AND. y /= 0.0_dp) th = 0.5_dp * pi
  IF (x0 /= 0.0_dp) th = ATAN(y/x0)
  z2 = x0 * x0 + y * y
  z0 = SQRT(z2)
  psr = LOG(z0) - 0.5_dp * x0 / z2
  psi = th + 0.5_dp * y / z2
  DO  k = 1, 8
    psr = psr + a(k) * z2 ** (-k) * COS(2.0_dp*k*th)
    psi = psi - a(k) * z2 ** (-k) * SIN(2.0_dp*k*th)
  END DO
  IF (x < 8.0_dp) THEN
    rr = 0.0_dp
    ri = 0.0_dp
    DO  k = 1, n
      rr = rr + (x0-k) / ((x0-k)**2.0_dp + y*y)
      ri = ri + y / ((x0-k)**2.0_dp + y*y)
    END DO
    psr = psr - rr
    psi = psi + ri
  END IF
  IF (x1 < 0.0_dp) THEN
    tn = TAN(pi*x)
    tm = TANH(pi*y)
    ct2 = tn * tn + tm * tm
    psr = psr + x / (x*x + y*y) + pi * (tn - tn*tm*tm) / ct2
    psi = psi - y / (x*x + y*y) - pi * tm * (1.0_dp + tn*tn) / ct2
    x = x1
    y = y1
  END IF
END IF
RETURN
END SUBROUTINE cpsi
 
END MODULE cpsi_func
 
 
 
PROGRAM mcpsi
USE cpsi_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:37

!       =========================================================
!       Purpose: This program computes the psi function psi(z)
!                for a complex argument using subroutine CPSI
!       Input :  x   --- Real part of z
!                y   --- Imaginary part of z
!       Output:  PSR --- Real part of psi(z)
!                PSI --- Imaginary part of psi(z)
!       Examples:
!                   x       y      Re[psi(z)]     Im[psi(z)]
!                 -------------------------------------------
!                  3.0     2.0     1.16459152      .67080728
!                  3.0    -2.0     1.16459152     -.67080728
!                 -3.0     2.0     1.39536075     2.62465344
!                 -3.0    -2.0     1.39536075    -2.62465344
!       =========================================================

REAL (dp)  :: x, y, psr, psi

WRITE (*,*) 'Please enter x and y ( z=x+iy): '
READ (*,*) x, y
WRITE (*,5100) x, y
WRITE (*,*)
WRITE (*,*) '   x       y      Re[Psi(z)]      Im[Psi(z)]'
WRITE (*,*) ' ----------------------------------------------'
CALL cpsi(x, y, psr, psi)
WRITE (*,5000) x, y, psr, psi
STOP

5000 FORMAT (' ', f5.1, '   ', f5.1, '  ', 2g16.8)
5100 FORMAT (' x=', f6.2, '      y=', f6.2)
END PROGRAM mcpsi
