MODULE ittikb_func
 
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


SUBROUTINE ittikb(x,tti,ttk)

!    =========================================================
!    Purpose: Integrate [I0(t)-1]/t with respect to t from 0
!             to x, and K0(t)/t with respect to t from x to ì
!    Input :  x   --- Variable in the limits  ( x ò 0 )
!    Output:  TTI --- Integration of [I0(t)-1]/t from 0 to x
!             TTK --- Integration of K0(t)/t from x to ì
!    =========================================================

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: tti
REAL (dp), INTENT(OUT)  :: ttk

REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp, el = .5772156649015329_dp
REAL (dp)  :: e0, t, t1, x1

IF (x == 0.0D0) THEN
  tti = 0.0D0
ELSE IF (x <= 5.0D0) THEN
  x1 = x / 5.0D0
  t = x1 * x1
  tti = (((((((.1263D-3*t + .96442D-3)*t + .968217D-2)*t + .06615507D0)*t  &
        + .33116853D0)*t + 1.13027241D0)*t + 2.44140746D0)*t +  3.12499991D0)*t
ELSE
  t = 5.0D0 / x
  tti = (((((((((2.1945464D0*t - 3.5195009D0)*t - 11.9094395D0)*t +   &
        40.394734D0)*t - 48.0524115D0)*t + 28.1221478D0)*t - 8.6556013D0)*t  &
        + 1.4780044D0)*t - .0493843D0)*t + .1332055D0)*t + .3989314D0
  tti = tti * EXP(x) / (SQRT(x)*x)
END IF
IF (x == 0.0D0) THEN
  ttk = 1.0D+300
ELSE IF (x <= 2.0D0) THEN
  t1 = x / 2.0D0
  t = t1 * t1
  ttk = (((((.77D-6*t + .1544D-4)*t + .48077D-3)*t + .925821D-2)*t +   &
        .10937537D0)*t + .74999993D0)*t
  e0 = el + LOG(x/2.0D0)
  ttk = pi * pi / 24.0D0 + e0 * (.5D0*e0 + tti) - ttk
ELSE IF (x <= 4.0D0) THEN
  t = 2.0D0 / x
  ttk = (((.06084D0*t - .280367D0)*t + .590944D0)*t - .850013D0)*t + 1.234684D0
  ttk = ttk * EXP(-x) / (SQRT(x)*x)
ELSE
  t = 4.0D0 / x
  ttk = (((((.02724D0*t - .1110396D0)*t + .2060126D0)*t - .2621446D0)*t +   &
        .3219184D0)*t - .5091339D0)*t + 1.2533141D0
  ttk = ttk * EXP(-x) / (SQRT(x)*x)
END IF
RETURN
END SUBROUTINE ittikb
 
END MODULE ittikb_func
 
 
 
PROGRAM mittikb
USE ittikb_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:41

!    ============================================================
!    Purpose: This program computes the integral of [I0(t)-1]/t
!             with respect to t from 0 to x and K0(t)/t with
!             respect to t from x to ì using subroutine ITTIKB
!    Input :  x   --- Upper limit of the integral
!    Output:  TTI --- Integration of [I0(t)-1]/t from 0 to x
!             TTK --- Integration of K0(t)/t from x to ì
!    Example:
!                x     [1-I0(t)]/tdt      K0(t)/tdt
!             ---------------------------------------
!               5.0     .710478D+01     .586361D-03
!              10.0     .340811D+03     .156293D-05
!              15.0     .254373D+05     .598363D-08
!              20.0     .236735D+07     .267906D-10
!              25.0     .246534D+09     .131007D-12
!    ============================================================

REAL (dp)  :: x, tti, ttk

WRITE (*,*) 'Please enter x '
READ (*,*) x
WRITE (*,*) '   x     [1-I0(t)]/tdt      K0(t)/tdt'
WRITE (*,*) '---------------------------------------'
CALL ittikb(x, tti, ttk)
WRITE (*,5000) x, tti, ttk
STOP

5000 FORMAT (' ', f5.1, 2g16.6)
END PROGRAM mittikb
