MODULE ik01b_func
 
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
 

SUBROUTINE ik01b(x, bi0, di0, bi1, di1, bk0, dk0, bk1, dk1)

!       =========================================================
!       Purpose: Compute modified Bessel functions I0(x), I1(1),
!                K0(x) and K1(x), and their derivatives
!       Input :  x   --- Argument ( x ò 0 )
!       Output:  BI0 --- I0(x)
!                DI0 --- I0'(x)
!                BI1 --- I1(x)
!                DI1 --- I1'(x)
!                BK0 --- K0(x)
!                DK0 --- K0'(x)
!                BK1 --- K1(x)
!                DK1 --- K1'(x)
!       =========================================================

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: bi0
REAL (dp), INTENT(OUT)  :: di0
REAL (dp), INTENT(OUT)  :: bi1
REAL (dp), INTENT(OUT)  :: di1
REAL (dp), INTENT(OUT)  :: bk0
REAL (dp), INTENT(OUT)  :: dk0
REAL (dp), INTENT(OUT)  :: bk1
REAL (dp), INTENT(OUT)  :: dk1

REAL (dp)  :: t, t2

IF (x == 0.0_dp) THEN
  bi0 = 1.0_dp
  bi1 = 0.0_dp
  bk0 = 1.0D+300
  bk1 = 1.0D+300
  di0 = 0.0_dp
  di1 = 0.5_dp
  dk0 = -1.0D+300
  dk1 = -1.0D+300
  RETURN
ELSE IF (x <= 3.75_dp) THEN
  t = x / 3.75_dp
  t2 = t * t
  bi0 = (((((.0045813_dp*t2 + .0360768_dp)*t2 + .2659732_dp)*t2 +   &
      1.2067492_dp)*t2 + 3.0899424_dp)*t2 + 3.5156229_dp) * t2 + 1.0_dp
  bi1 = x * ((((((.00032411_dp*t2 + .00301532_dp)*t2 + .02658733_dp)*t2 +  &
      .15084934_dp)*t2 + .51498869_dp)*t2 + .87890594_dp)*t2 + .5_dp)
ELSE
  t = 3.75_dp / x
  bi0 = ((((((((.00392377_dp*t - .01647633_dp)*t + .02635537_dp)*t -   &
      .02057706_dp)*t + .916281D-2)*t - .157565D-2)*t + .225319D-2)*t +  &
      .01328592_dp)*t + .39894228_dp) * EXP(x) / SQRT(x)
  bi1 = ((((((((-.420059D-2*t + .01787654_dp)*t - .02895312_dp)*t +   &
      .02282967_dp)*t - .01031555_dp)*t + .163801D-2)*t - .00362018_dp)*t -  &
      .03988024_dp)*t + .39894228_dp) * EXP(x) / SQRT(x)
END IF
IF (x <= 2.0_dp) THEN
  t = x / 2.0_dp
  t2 = t * t
  bk0 = (((((.0000074_dp*t2 + .0001075_dp)*t2 + .00262698_dp)*t2 +  &
      .0348859_dp)*t2 + .23069756_dp)*t2 + .4227842_dp) * t2 -  &
      .57721566_dp - bi0 * LOG(t)
  bk1 = ((((((-.00004686_dp*t2 - .00110404_dp)*t2 - .01919402_dp)*t2 -  &
      .18156897_dp)*t2 - .67278579_dp)*t2 + .15443144_dp)*t2 + 1.0_dp) / x + bi1*LOG(t)
ELSE
  t = 2.0_dp / x
  t2 = t * t
  bk0 = ((((((.00053208_dp*t - .0025154_dp)*t + .00587872_dp)*t -   &
      .01062446_dp)*t + .02189568_dp)*t - .07832358_dp)*t + 1.25331414_dp) *  &
      EXP(-x) / SQRT(x)
  bk1 = ((((((-.00068245_dp*t + .00325614_dp)*t - .00780353_dp)*t +   &
      .01504268_dp)*t - .0365562_dp)*t + .23498619_dp)*t + 1.25331414_dp) *  &
      EXP(-x) / SQRT(x)
END IF
di0 = bi1
di1 = bi0 - bi1 / x
dk0 = -bk1
dk1 = -bk0 - bk1 / x
RETURN
END SUBROUTINE ik01b
 
END MODULE ik01b_func
 
 
 
PROGRAM mik01b
USE ik01b_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:40

!       =============================================================
!       Purpose: This program computes the modified Bessel functions
!                I0(x), I1(x), K0(x), K1(x), and their derivatives
!                using subroutine IK01B
!       Input :  x   --- Argument ( x ò 0 )
!       Output:  BI0 --- I0(x)
!                DI0 --- I0'(x)
!                BI1 --- I1(x)
!                DI1 --- I1'(x)
!                BK0 --- K0(x)
!                DK0 --- K0'(x)
!                BK1 --- K1(x)
!                DK1 --- K1'(x)
!       Example:

!         x      I0(x)         I0'(x)        I1(x)         I1'(x)
!       -------------------------------------------------------------
!        1.0   .126607D+01   .565159D+00   .565159D+00   .700907D+00
!       10.0   .281572D+04   .267099D+04   .267099D+04   .254862D+04
!       20.0   .435583D+08   .424550D+08   .424550D+08   .414355D+08
!       30.0   .781672D+12   .768532D+12   .768532D+12   .756055D+12
!       40.0   .148948D+17   .147074D+17   .147074D+17   .145271D+17
!       50.0   .293255D+21   .290308D+21   .290308D+21   .287449D+21

!         x      K0(x)         K0'(x)        K1(x)         K1'(x)
!       -------------------------------------------------------------
!        1.0   .421024D+00  -.601907D+00   .601907D+00  -.102293D+01
!       10.0   .177801D-04  -.186488D-04   .186488D-04  -.196449D-04
!       20.0   .574124D-09  -.588306D-09   .588306D-09  -.603539D-09
!       30.0   .213248D-13  -.216773D-13   .216773D-13  -.220474D-13
!       40.0   .839286D-18  -.849713D-18   .849713D-18  -.860529D-18
!       50.0   .341017D-22  -.344410D-22   .344410D-22  -.347905D-22
!       =============================================================

REAL (dp)  :: x, bi0, di0, bi1, di1, bk0, dk0, bk1, dk1

WRITE (*,*) 'Please enter x '
READ (*,*) x
WRITE (*,5000) x
WRITE (*,*) '  x       I0(x)          I0''(X)         I1(X)          I1''(X)'
WRITE (*,*) '-----------------------------------------------------------------'
CALL ik01b(x, bi0, di0, bi1, di1, bk0, dk0, bk1, dk1)
WRITE (*,5100) x, bi0, di0, bi1, di1
WRITE (*,*)
WRITE (*,*) '  x       K0(x)          K0''(X)         K1(X)          K1''(X)'
WRITE (*,*) '-----------------------------------------------------------------'
WRITE (*,5100) x, bk0, dk0, bk1, dk1
STOP

5000 FORMAT ('   x =', f5.1)
5100 FORMAT (' ', f4.1, 4g15.7)
END PROGRAM mik01b
