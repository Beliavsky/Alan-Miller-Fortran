MODULE fcs_func
 
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


SUBROUTINE fcs(x, c, s)

!       =================================================
!       Purpose: Compute Fresnel integrals C(x) and S(x)
!       Input :  x --- Argument of C(x) and S(x)
!       Output:  C --- C(x)
!                S --- S(x)
!       =================================================

REAL (dp), INTENT(IN OUT)  :: x
REAL (dp), INTENT(OUT)     :: c
REAL (dp), INTENT(OUT)     :: s

REAL (dp), PARAMETER  :: pi = 3.141592653589793D0, eps = 1.0D-15
REAL (dp)  :: f, f0, f1, g, px, q, r, su, t, t0, t2, xa
INTEGER    :: k, m

xa = ABS(x)
px = pi * xa
t = .5D0 * px * xa
t2 = t * t
IF (xa == 0.0) THEN
  c = 0.0D0
  s = 0.0D0
ELSE IF (xa < 2.5D0) THEN
  r = xa
  c = r
  DO  k = 1, 50
    r = -.5D0 * r * (4.0D0*k-3.0D0) / k / (2.0D0*k-1.0D0) / (4.0D0*k+1.0D0) * t2
    c = c + r
    IF (ABS(r) < ABS(c)*eps) EXIT
  END DO
  s = xa * t / 3.0D0
  r = s
  DO  k = 1, 50
    r = -.5D0 * r * (4.0D0*k-1.0D0) / k / (2.0D0*k+1.0D0) / (4.0D0*k+3.0D0) * t2
    s = s + r
    IF (ABS(r) < ABS(s)*eps) GO TO 70
  END DO
ELSE IF (xa < 4.5D0) THEN
  m = INT(42.0+1.75*t)
  su = 0.0D0
  c = 0.0D0
  s = 0.0D0
  f1 = 0.0D0
  f0 = 1.0D-100
  DO  k = m, 0, -1
    f = (2*k+3) * f0 / t - f1
    IF (k == INT(k/2)*2) THEN
      c = c + f
    ELSE
      s = s + f
    END IF
    su = su + (2.0D0*k+1.0D0) * f * f
    f1 = f0
    f0 = f
  END DO
  q = SQRT(su)
  c = c * xa / q
  s = s * xa / q
ELSE
  r = 1.0D0
  f = 1.0D0
  DO  k = 1, 20
    r = -.25D0 * r * (4*k-1) * (4*k-3) / t2
    f = f + r
  END DO
  r = 1.0D0 / (px*xa)
  g = r
  DO  k = 1, 12
    r = -.25D0 * r * (4*k+1) * (4*k-1) / t2
    g = g + r
  END DO
  t0 = t - INT(t/(2.0D0*pi)) * 2.0D0 * pi
  c = .5D0 + (f*SIN(t0)-g*COS(t0)) / px
  s = .5D0 - (f*COS(t0)+g*SIN(t0)) / px
END IF

70 IF (x < 0.0D0) THEN
  c = -c
  s = -s
END IF
RETURN
END SUBROUTINE fcs
 
END MODULE fcs_func
 
 
 
PROGRAM mfcs
USE fcs_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:39

!       =====================================================
!       Purpose: This program computes the Fresnel integrals
!                C(x) and S(x) using subroutine FCS
!       Input :  x --- Argument of C(x) and S(x)
!       Output:  C --- C(x)
!                S --- S(x)
!       Example:
!                  x          C(x)          S(x)
!                -----------------------------------
!                 0.0      .00000000      .00000000
!                 0.5      .49234423      .06473243
!                 1.0      .77989340      .43825915
!                 1.5      .44526118      .69750496
!                 2.0      .48825341      .34341568
!                 2.5      .45741301      .61918176
!       =====================================================

REAL (dp)  :: c, s, x

WRITE (*,*) 'Please enter x '
READ (*,*) x
WRITE (*,*) '   x          C(x)          S(x)'
WRITE (*,*) ' -----------------------------------'
CALL fcs(x, c, s)
WRITE (*,5000) x, c, s
STOP

5000 FORMAT (' ', f5.1, 2F15.8)
END PROGRAM mfcs
