MODULE cjy01_func
 
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
 

SUBROUTINE cjy01(z, cbj0, cdj0, cbj1, cdj1, cby0, cdy0, cby1, cdy1)

!       =======================================================
!       Purpose: Compute Bessel functions J0(z), J1(z), Y0(z),
!                Y1(z), and their derivatives for a complex
!                argument
!       Input :  z --- Complex argument
!       Output:  CBJ0 --- J0(z)
!                CDJ0 --- J0'(z)
!                CBJ1 --- J1(z)
!                CDJ1 --- J1'(z)
!                CBY0 --- Y0(z)
!                CDY0 --- Y0'(z)
!                CBY1 --- Y1(z)
!                CDY1 --- Y1'(z)
!       =======================================================

COMPLEX (dp), INTENT(IN)   :: z
COMPLEX (dp), INTENT(OUT)  :: cbj0
COMPLEX (dp), INTENT(OUT)  :: cdj0
COMPLEX (dp), INTENT(OUT)  :: cbj1
COMPLEX (dp), INTENT(OUT)  :: cdj1
COMPLEX (dp), INTENT(OUT)  :: cby0
COMPLEX (dp), INTENT(OUT)  :: cdy0
COMPLEX (dp), INTENT(OUT)  :: cby1
COMPLEX (dp), INTENT(OUT)  :: cdy1

REAL (dp), PARAMETER  :: a(12) = (/ -0.703125D-01, .112152099609375_dp,   &
     -.5725014209747314_dp, .6074042001273483D+01, -  &
      .1100171402692467D+03, .3038090510922384D+04, -  &
      .1188384262567832D+06, .6252951493434797D+07, -  &
      .4259392165047669D+09, .3646840080706556D+11, -  &
      .3833534661393944D+13, .4854014686852901D+15 /)
REAL (dp), PARAMETER  :: b(12) = (/.732421875D-01, -.2271080017089844_dp,  &
      .1727727502584457D+01, -.2438052969955606D+02,  &
      .5513358961220206D+03, -.1825775547429318D+05,  &
      .8328593040162893D+06, -.5006958953198893D+08,  &
      .3836255180230433D+10, -.3649010818849833D+12,  &
      .4218971570284096D+14, -.5827244631566907D+16 /)
REAL (dp), PARAMETER  :: a1(12) = (/ 0.1171875_dp, -.144195556640625_dp,  &
      .6765925884246826_dp, -.6883914268109947D+01,  &
      .1215978918765359D+03, -.3302272294480852D+04,  &
      .1276412726461746D+06, -.6656367718817688D+07,  &
      .4502786003050393D+09, -.3833857520742790D+11,  &
      .4011838599133198D+13, -.5060568503314727D+15 /)
REAL (dp), PARAMETER  :: b1(12) = (/ -0.1025390625_dp, .2775764465332031_dp,  &
     -.1993531733751297D+01, .2724882731126854D+02, -  &
      .6038440767050702D+03, .1971837591223663D+05, -  &
      .8902978767070678D+06, .5310411010968522D+08, -  &
      .4043620325107754D+10, .3827011346598605D+12, -  &
      .4406481417852278D+14, .6065091351222699D+16 /)

REAL (dp)     :: a0, rp2, w0, w1
INTEGER       :: k, k0
COMPLEX (dp)  :: ci, cp, cp0, cp1, cq0, cq1, cr, cs, ct1, ct2, cu, z1, z2
REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp, el = 0.5772156649015329_dp

rp2 = 2.0_dp / pi
ci = (0.0_dp,1.0_dp)
a0 = ABS(z)
z2 = z * z
z1 = z
IF (a0 == 0.0_dp) THEN
  cbj0 = (1.0_dp,0.0_dp)
  cbj1 = (0.0_dp,0.0_dp)
  cdj0 = (0.0_dp,0.0_dp)
  cdj1 = (0.5_dp,0.0_dp)
  cby0 = -(1.0D300,0.0_dp)
  cby1 = -(1.0D300,0.0_dp)
  cdy0 = (1.0D300,0.0_dp)
  cdy1 = (1.0D300,0.0_dp)
  RETURN
END IF
IF (REAL(z) < 0.0) z1 = -z
IF (a0 <= 12.0) THEN
  cbj0 = (1.0_dp,0.0_dp)
  cr = (1.0_dp,0.0_dp)
  DO  k = 1, 40
    cr = -0.25_dp * cr * z2 / (k*k)
    cbj0 = cbj0 + cr
    IF (ABS(cr) < ABS(cbj0)*1.0D-15) EXIT
  END DO
  cbj1 = (1.0_dp,0.0_dp)
  cr = (1.0_dp,0.0_dp)
  DO  k = 1, 40
    cr = -0.25_dp * cr * z2 / (k*(k+1))
    cbj1 = cbj1 + cr
    IF (ABS(cr) < ABS(cbj1)*1.0D-15) EXIT
  END DO
  cbj1 = 0.5_dp * z1 * cbj1
  w0 = 0.0_dp
  cr = (1.0_dp,0.0_dp)
  cs = (0.0_dp,0.0_dp)
  DO  k = 1, 40
    w0 = w0 + 1.0_dp / k
    cr = -0.25_dp * cr / (k*k) * z2
    cp = cr * w0
    cs = cs + cp
    IF (ABS(cp) < ABS(cs)*1.0D-15) EXIT
  END DO
  cby0 = rp2 * (LOG(z1/2.0_dp)+el) * cbj0 - rp2 * cs
  w1 = 0.0_dp
  cr = (1.0_dp,0.0_dp)
  cs = (1.0_dp,0.0_dp)
  DO  k = 1, 40
    w1 = w1 + 1.0_dp / k
    cr = -0.25_dp * cr / (k*(k+1)) * z2
    cp = cr * (2.0_dp*w1+1.0_dp/(k+1.0_dp))
    cs = cs + cp
    IF (ABS(cp) < ABS(cs)*1.0D-15) EXIT
  END DO
  cby1 = rp2 * ((LOG(z1/2.0_dp) + el)*cbj1 - 1.0_dp/z1 - .25_dp*z1*cs)
ELSE
  k0 = 12
  IF (a0 >= 35.0) k0 = 10
  IF (a0 >= 50.0) k0 = 8
  ct1 = z1 - .25_dp * pi
  cp0 = (1.0_dp,0.0_dp)
  DO  k = 1, k0
    cp0 = cp0 + a(k) * z1 ** (-2*k)
  END DO
  cq0 = -0.125_dp / z1
  DO  k = 1, k0
    cq0 = cq0 + b(k) * z1 ** (-2*k-1)
  END DO
  cu = SQRT(rp2/z1)
  cbj0 = cu * (cp0*COS(ct1) - cq0*SIN(ct1))
  cby0 = cu * (cp0*SIN(ct1) + cq0*COS(ct1))
  ct2 = z1 - .75_dp * pi
  cp1 = (1.0_dp,0.0_dp)
  DO  k = 1, k0
    cp1 = cp1 + a1(k) * z1 ** (-2*k)
  END DO
  cq1 = 0.375_dp / z1
  DO  k = 1, k0
    cq1 = cq1 + b1(k) * z1 ** (-2*k-1)
  END DO
  cbj1 = cu * (cp1*COS(ct2) - cq1*SIN(ct2))
  cby1 = cu * (cp1*SIN(ct2) + cq1*COS(ct2))
END IF
IF (REAL(z) < 0.0) THEN
  IF (AIMAG(z) < 0.0) cby0 = cby0 - 2.0_dp * ci * cbj0
  IF (AIMAG(z) > 0.0) cby0 = cby0 + 2.0_dp * ci * cbj0
  IF (AIMAG(z) < 0.0) cby1 = -(cby1 - 2.0_dp*ci*cbj1)
  IF (AIMAG(z) > 0.0) cby1 = -(cby1 + 2.0_dp*ci*cbj1)
  cbj1 = -cbj1
END IF
cdj0 = -cbj1
cdj1 = cbj0 - 1.0_dp / z * cbj1
cdy0 = -cby1
cdy1 = cby0 - 1.0_dp / z * cby1
RETURN
END SUBROUTINE cjy01
 
END MODULE cjy01_func
 
 
 
PROGRAM mcjy01
USE cjy01_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:36

!       ================================================================
!       Purpose: This program computes Bessel functions J0(z), J1(z),
!                Y0(z), Y1(z), and their derivatives for a complex
!                argument using subroutine CJY01
!       Input :  z --- Complex argument
!       Output:  CBJ0 --- J0(z)
!                CDJ0 --- J0'(z)
!                CBJ1 --- J1(z)
!                CDJ1 --- J1'(z)
!                CBY0 --- Y0(z)
!                CDY0 --- Y0'(z)
!                CBY1 --- Y1(z)
!                CDY1 --- Y1'(z)
!       Example: z =  4.0 + i  2.0

!     n     Re[Jn(z)]       Im[Jn(z)]       Re[Jn'(z)]      Im[Jn'(z)]
!   --------------------------------------------------------------------
!     0  -.13787022D+01   .39054236_dp   .50735255_dp   .12263041D+01
!     1  -.50735255_dp  -.12263041D+01  -.11546013D+01   .58506793_dp

!     n     Re[Yn(z)]       Im[Yn(z)]       Re[Yn'(z)]      Im[Yn'(z)]
!   --------------------------------------------------------------------
!     0  -.38145893_dp  -.13291649D+01  -.12793101D+01   .51220420_dp
!     1   .12793101D+01  -.51220420_dp  -.58610052_dp  -.10987930D+01
!       ================================================================

COMPLEX (dp)  :: z, cbj0, cdj0, cbj1, cdj1, cby0, cdy0, cby1, cdy1
REAL (dp)     :: x, y

WRITE (*,*) '  Please enter x,y (z=x+iy) '
READ (*,*) x, y
z = CMPLX(x, y, KIND=dp)
WRITE (*,5200) x, y
CALL cjy01(z, cbj0, cdj0, cbj1, cdj1, cby0, cdy0, cby1, cdy1)
WRITE(*,*)
WRITE(*,*) '  n      Re[Jn(z)]       Im[Jn(z)]       Re[Jn''(Z)]      IM[JN''(Z)]'
WRITE(*,*) ' --------------------------------------------------------------------'
WRITE(*,5000) cbj0, cdj0
WRITE(*,5100) cbj1, cdj1
WRITE(*,*)
WRITE(*,*) '  n      Re[Yn(z)]       Im[Yn(z)]       Re[Yn''(Z)]      IM[YN''(Z)]'
WRITE(*,*) ' --------------------------------------------------------------------'
WRITE(*,5000) cby0, cdy0
WRITE(*,5100) cby1, cdy1
STOP

5000 FORMAT ('   0', '  ', 4g16.8)
5100 FORMAT ('   1', '  ', 4g16.8)
5200 FORMAT ('   z =', f6.2, ' + i', f6.2)
END PROGRAM mcjy01
