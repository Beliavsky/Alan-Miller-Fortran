MODULE ittjyb_func
 
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
 

SUBROUTINE ittjyb(x, ttj, tty)

!       ==========================================================
!       Purpose: Integrate [1-J0(t)]/t with respect to t from 0
!                to x, and Y0(t)/t with respect to t from x to ì
!       Input :  x   --- Variable in the limits  ( x ò 0 )
!       Output:  TTJ --- Integration of [1-J0(t)]/t from 0 to x
!                TTY --- Integration of Y0(t)/t from x to ì
!       ==========================================================

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: ttj
REAL (dp), INTENT(OUT)  :: tty

REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp, el = .5772156649015329_dp
REAL (dp)  :: e0, f0, g0, t, t1, x1, xt

IF (x == 0.0D0) THEN
  ttj = 0.0D0
  tty = -1.0D+300
ELSE IF (x <= 4.0D0) THEN
  x1 = x / 4.0D0
  t = x1 * x1
  ttj = ((((((.35817D-4*t - .639765D-3)*t + .7092535D-2)*t -  &
        .055544803D0)*t + .296292677D0)*t - .999999326D0)*t +  1.999999936D0)*t
  tty = (((((((-.3546D-5*t + .76217D-4)*t - .1059499D-2)*t +   &
        .010787555D0)*t-.07810271D0)*t + .377255736D0)*t - 1.114084491D0)*t  &
        + 1.909859297D0)*t
  e0 = el + LOG(x/2.0D0)
  tty = pi / 6.0D0 + e0 / pi * (2.0D0*ttj - e0) - tty
ELSE IF (x <= 8.0D0) THEN
  xt = x + .25D0 * pi
  t1 = 4.0D0 / x
  t = t1 * t1
  f0 = (((((.0145369D0*t - .0666297D0)*t + .1341551D0)*t - .1647797D0)*t  &
       + .1608874D0)*t - .2021547D0)*t + .7977506D0
  g0 = ((((((.0160672D0*t - .0759339D0)*t + .1576116D0)*t - .1960154D0)*t  &
       + .1797457D0)*t - .1702778D0)*t + .3235819D0)*t1
  ttj = (f0*COS(xt) + g0*SIN(xt)) / (SQRT(x)*x)
  ttj = ttj + el + LOG(x/2.0D0)
  tty = (f0*SIN(xt) - g0*COS(xt)) / (SQRT(x)*x)
ELSE
  t = 8.0D0 / x
  xt = x + .25D0 * pi
  f0 = (((((.18118D-2*t - .91909D-2)*t + .017033D0)*t - .9394D-3)*t -  &
       .051445D0)*t - .11D-5)*t + .7978846D0
  g0 = (((((-.23731D-2*t + .59842D-2)*t + .24437D-2)*t - .0233178D0)*t +   &
       .595D-4)*t + .1620695D0)*t
  ttj = (f0*COS(xt) + g0*SIN(xt)) / (SQRT(x)*x) + el + LOG(x/2.0D0)
  tty = (f0*SIN(xt)-g0*COS(xt)) / (SQRT(x)*x)
END IF
RETURN
END SUBROUTINE ittjyb
 
END MODULE ittjyb_func
 
 
 
PROGRAM mittjyb
USE ittjyb_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:41

!    ===========================================================
!    Purpose: This program computes the integral of [1-J0(t)]/t
!             with respect to t from 0 to x and Y0(t)/t with
!             respect to t from x to ì using subroutine ITTJYB
!    Input :  x   --- Variable in the limits  ( x ò 0 )
!    Output:  TTJ --- Integration of [1-J0(t)]/t from 0 to x
!             TTY --- Integration of Y0(t)/t from x to ì
!    Example:
!               x      [1-J0(t)]/tdt       Y0(t)/tdt
!             ----------------------------------------
!              5.0     .1540347D+01    -.4632208D-01
!             10.0     .2177866D+01    -.2298791D-01
!             15.0     .2578551D+01     .3857453D-03
!             20.0     .2877311D+01     .8503154D-02
!             25.0     .3108231D+01     .3526339D-02
!    ===========================================================

REAL (dp)  :: ttj, tty, x

WRITE (*,*) 'Please enter x '
READ (*,*) x
WRITE (*,*) '   x     [1-J0(t)]/tdt      Y0(t)/tdt'
WRITE (*,*) '----------------------------------------'
CALL ittjyb(x, ttj, tty)
WRITE (*,5000) x, ttj, tty
STOP

5000 FORMAT (' ', f5.1, 2g17.7)
END PROGRAM mittjyb
