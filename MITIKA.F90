MODULE itika_func
 
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


SUBROUTINE itika(x, ti, tk)

!    =======================================================
!    Purpose: Integrate modified Bessel functions I0(t) and
!             K0(t) with respect to t from 0 to x
!    Input :  x  --- Upper limit of the integral  ( x ò 0 )
!    Output:  TI --- Integration of I0(t) from 0 to x
!             TK --- Integration of K0(t) from 0 to x
!    =======================================================

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: ti
REAL (dp), INTENT(OUT)  :: tk

REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp, el = .5772156649015329_dp
REAL (dp), PARAMETER  :: a(10) = (/ 0.625D0, 1.0078125D0, 2.5927734375D0, &
    9.1868591308594D0, 4.1567974090576D+1, 2.2919635891914D+2,  &
    1.491504060477D+3, 1.1192354495579D+4, 9.515939374212D+4,  &
    9.0412425769041D+5 /)
REAL (dp)  :: b1, b2, e0, r, rc1, rc2, rs, tw, x2
INTEGER    :: k

IF (x == 0.0D0) THEN
  ti = 0.0D0
  tk = 0.0D0
  RETURN
ELSE IF (x < 20.0D0) THEN
  x2 = x * x
  ti = 1.0D0
  r = 1.0D0
  DO  k = 1, 50
    r = .25D0 * r * (2*k-1.0D0) / (2*k+1.0D0) / (k*k) * x2
    ti = ti + r
    IF (ABS(r/ti) < 1.0D-12) EXIT
  END DO
  ti = ti * x
ELSE
  ti = 1.0D0
  r = 1.0D0
  DO  k = 1, 10
    r = r / x
    ti = ti + a(k) * r
  END DO
  rc1 = 1.0D0 / SQRT(2.0D0*pi*x)
  ti = rc1 * EXP(x) * ti
END IF
IF (x < 12.0D0) THEN
  e0 = el + LOG(x/2.0D0)
  b1 = 1.0D0 - e0
  b2 = 0.0D0
  rs = 0.0D0
  r = 1.0D0
  DO  k = 1, 50
    r = .25D0 * r * (2*k-1.0D0) / (2*k+1.0D0) / (k*k) * x2
    b1 = b1 + r * (1.0D0/(2*k+1) - e0)
    rs = rs + 1.0D0 / k
    b2 = b2 + r * rs
    tk = b1 + b2
    IF (ABS((tk-tw)/tk) < 1.0D-12) EXIT
    tw = tk
  END DO
  tk = tk * x
ELSE
  tk = 1.0D0
  r = 1.0D0
  DO  k = 1, 10
    r = -r / x
    tk = tk + a(k) * r
  END DO
  rc2 = SQRT(pi/(2.0D0*x))
  tk = pi / 2.0D0 - rc2 * tk * EXP(-x)
END IF
RETURN
END SUBROUTINE itika
 
END MODULE itika_func
 
 
 
PROGRAM mitika
USE itika_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:40

!    ============================================================
!    Purpose: This program evaluates the integral of modified
!             Bessel functions I0(t) and K0(t) with respect to t
!             from 0 to x using subroutine ITIKA
!    Input :  x  --- Upper limit of the integral  ( x ò 0 )
!    Output:  TI --- Integration of I0(t) from 0 to x
!             TK --- Integration of K0(t) from 0 to x
!    Example:
!                 x         I0(t)dt         K0(t)dt
!              --------------------------------------
!                5.0    .31848668D+02     1.56738739
!               10.0    .29930445D+04     1.57077931
!               15.0    .35262048D+06     1.57079623
!               20.0    .44758593D+08     1.57079633
!               25.0    .58991731D+10     1.57079633
!    ============================================================

REAL (dp)  :: ti, tk, x

WRITE (*,*) 'Please enter x '
READ (*,*) x
WRITE (*,*) '   x         I0(t)dt         K0(t)dt'
WRITE (*,*) ' --------------------------------------'
CALL itika(x, ti, tk)
WRITE (*,5000) x, ti, tk
STOP

5000 FORMAT (' ', f5.1, g17.8, f15.8)
END PROGRAM mitika
