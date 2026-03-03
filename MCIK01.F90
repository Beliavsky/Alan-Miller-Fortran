MODULE cik01_func
 
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
 

SUBROUTINE cik01(z,cbi0,cdi0,cbi1,cdi1,cbk0,cdk0,cbk1,cdk1)

!       ==========================================================
!       Purpose: Compute modified Bessel functions I0(z), I1(z),
!                K0(z), K1(z), and their derivatives for a
!                complex argument
!       Input :  z --- Complex argument
!       Output:  CBI0 --- I0(z)
!                CDI0 --- I0'(z)
!                CBI1 --- I1(z)
!                CDI1 --- I1'(z)
!                CBK0 --- K0(z)
!                CDK0 --- K0'(z)
!                CBK1 --- K1(z)
!                CDK1 --- K1'(z)
!       ==========================================================

COMPLEX (dp), INTENT(IN)   :: z
COMPLEX (dp), INTENT(OUT)  :: cbi0
COMPLEX (dp), INTENT(OUT)  :: cdi0
COMPLEX (dp), INTENT(OUT)  :: cbi1
COMPLEX (dp), INTENT(OUT)  :: cdi1
COMPLEX (dp), INTENT(OUT)  :: cbk0
COMPLEX (dp), INTENT(OUT)  :: cdk0
COMPLEX (dp), INTENT(OUT)  :: cbk1
COMPLEX (dp), INTENT(OUT)  :: cdk1

REAL (dp)     :: a0, w0
COMPLEX (dp)  :: ca, cb, ci, cr, cs, ct, cw, z1, z2, zr, zr2
INTEGER       :: k, k0
REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp
REAL (dp), PARAMETER  :: a(12) = (/ 0.125_dp, 7.03125D-2, 7.32421875D-2,  &
      1.1215209960938D-1, 2.2710800170898D-1, 5.7250142097473D-1,  &
      1.7277275025845_dp, 6.0740420012735_dp, 2.4380529699556D01,  &
      1.1001714026925D02, 5.5133589612202D02, 3.0380905109224D03 /)
REAL (dp), PARAMETER  :: b(12) = (/ -0.375_dp, -1.171875D-1, -1.025390625D-1, &
     -1.4419555664063D-1, -2.7757644653320D-1, -6.7659258842468D-1,  &
     -1.9935317337513_dp, -6.8839142681099_dp, -2.7248827311269D01,  &
     -1.2159789187654D02, -6.0384407670507D02, -3.3022722944809D03 /)
REAL (dp), PARAMETER  :: a1(10) =  (/ 0.125_dp, 0.2109375_dp, 1.0986328125_dp,  &
      1.1775970458984D01, 2.1461706161499D02, 5.9511522710323D03,  &
      2.3347645606175D05, 1.2312234987631D07, 8.401390346421D08,  &
      7.2031420482627D10 /)

ci = (0.0_dp, 1.0_dp)
a0 = ABS(z)
z2 = z * z
z1 = z
IF (a0 == 0.0_dp) THEN
  cbi0 = (1.0_dp,0.0_dp)
  cbi1 = (0.0_dp,0.0_dp)
  cdi0 = (0.0_dp,0.0_dp)
  cdi1 = (0.5_dp,0.0_dp)
  cbk0 = (1.0D+300,0.0_dp)
  cbk1 = (1.0D+300,0.0_dp)
  cdk0 = -(1.0D+300,0.0_dp)
  cdk1 = -(1.0D+300,0.0_dp)
  RETURN
END IF
IF (REAL(z) < 0.0) z1 = -z
IF (a0 <= 18.0) THEN
  cbi0 = (1.0_dp,0.0_dp)
  cr = (1.0_dp,0.0_dp)
  DO  k = 1, 50
    cr = 0.25_dp * cr * z2 / (k*k)
    cbi0 = cbi0 + cr
    IF (ABS(cr/cbi0) < 1.0D-15) EXIT
  END DO

  cbi1 = (1.0_dp,0.0_dp)
  cr = (1.0_dp,0.0_dp)
  DO  k = 1, 50
    cr = 0.25_dp * cr * z2 / (k*(k+1))
    cbi1 = cbi1 + cr
    IF (ABS(cr/cbi1) < 1.0D-15) EXIT
  END DO
  cbi1 = 0.5_dp * z1 * cbi1
ELSE
  k0 = 12
  IF (a0 >= 35.0) k0 = 9
  IF (a0 >= 50.0) k0 = 7
  ca = EXP(z1) / SQRT(2.0_dp*pi*z1)
  cbi0 = (1.0_dp,0.0_dp)
  zr = 1.0_dp / z1
  DO  k = 1, k0
    cbi0 = cbi0 + a(k) * zr ** k
  END DO
  cbi0 = ca * cbi0
  cbi1 = (1.0_dp,0.0_dp)
  DO  k = 1, k0
    cbi1 = cbi1 + b(k) * zr ** k
  END DO
  cbi1 = ca * cbi1
END IF
IF (a0 <= 9.0) THEN
  cs = (0.0_dp,0.0_dp)
  ct = -LOG(0.5_dp*z1) - 0.5772156649015329_dp
  w0 = 0.0_dp
  cr = (1.0_dp,0.0_dp)
  DO  k = 1, 50
    w0 = w0 + 1.0_dp / k
    cr = 0.25_dp * cr / (k*k) * z2
    cs = cs + cr * (w0+ct)
    IF (ABS((cs-cw)/cs) < 1.0D-15) EXIT
    cw = cs
  END DO
  cbk0 = ct + cs
ELSE
  cb = 0.5_dp / z1
  zr2 = 1.0_dp / z2
  cbk0 = (1.0_dp,0.0_dp)
  DO  k = 1, 10
    cbk0 = cbk0 + a1(k) * zr2 ** k
  END DO
  cbk0 = cb * cbk0 / cbi0
END IF
cbk1 = (1.0_dp/z1 - cbi1*cbk0) / cbi0
IF (REAL(z) < 0.0) THEN
  IF (AIMAG(z) < 0.0) cbk0 = cbk0 + ci * pi * cbi0
  IF (AIMAG(z) > 0.0) cbk0 = cbk0 - ci * pi * cbi0
  IF (AIMAG(z) < 0.0) cbk1 = -cbk1 + ci * pi * cbi1
  IF (AIMAG(z) > 0.0) cbk1 = -cbk1 - ci * pi * cbi1
  cbi1 = -cbi1
END IF
cdi0 = cbi1
cdi1 = cbi0 - 1.0_dp / z * cbi1
cdk0 = -cbk1
cdk1 = -cbk0 - 1.0_dp / z * cbk1
RETURN
END SUBROUTINE cik01
 
END MODULE cik01_func
 
 
 
PROGRAM mcik01
USE cik01_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:35

!       =============================================================
!       Purpose: This program computes the modified Bessel functions
!                I0(z), I1(z), K0(z), K1(z), and their derivatives
!                for a complex argument using subroutine CIK01
!       Input :  z --- Complex argument
!       Output:  CBI0 --- I0(z)
!                CDI0 --- I0'(z)
!                CBI1 --- I1(z)
!                CDI1 --- I1'(z)
!                CBK0 --- K0(z)
!                CDK0 --- K0'(z)
!                CBK1 --- K1(z)
!                CDK1 --- K1'(z)
!       Example: z = 20.0 + i 10.0

!     n      Re[In(z)]      Im[In(z)]      Re[In'(z)]     Im[In'(z)]
!    -----------------------------------------------------------------
!     0   -.38773811D+08 -.13750292D+08 -.37852037D+08 -.13869150D+08
!     1   -.37852037D+08 -.13869150D+08 -.36982347D+08 -.13952566D+08

!     n      Re[Kn(z)]      Im[Kn(z)]      Re[Kn'(z)]     Im[Kn'(z)]
!    -----------------------------------------------------------------
!     0   -.37692389D-09  .39171613D-09  .38056380D-09 -.40319029D-09
!     1   -.38056380D-09  .40319029D-09  .38408264D-09 -.41545502D-09
!       =============================================================

COMPLEX (dp)  :: cbi0, cdi0, cbi1, cdi1, cbk0, cdk0, cbk1, cdk1, z
REAL (dp)     :: x, y

WRITE (*,*) '  Please enter x and y (z=x+iy) '
READ (*,*) x, y
z = CMPLX(x, y, KIND=dp)
WRITE (*,5200) x, y
CALL cik01(z, cbi0, cdi0, cbi1, cdi1, cbk0, cdk0, cbk1, cdk1)
WRITE(*,*)
WRITE(*,*) '  n      Re[In(z)]      Im[In(z)]      Re[In''(Z)]     IM[IN''(Z)]'
WRITE(*,*) ' -----------------------------------------------------------------'
WRITE(*,5000) cbi0, cdi0
WRITE(*,5100) cbi1, cdi1
WRITE(*,*)
WRITE(*,*) '  n      Re[Kn(z)]      Im[Kn(z)]      Re[Kn''(Z)]     IM[KN''(Z)]'
WRITE(*,*) ' -----------------------------------------------------------------'
WRITE(*,5000) cbk0, cdk0
WRITE(*,5100) cbk1, cdk1
STOP

5000 FORMAT (t4, '0', '  ', 4g15.7)
5100 FORMAT (t4, '1', '  ', 4g15.7)
5200 FORMAT (t4, 'z =', f7.2, ' + i', f7.2)
END PROGRAM mcik01
