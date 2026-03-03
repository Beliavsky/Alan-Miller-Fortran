MODULE lgama_func
 
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
 

SUBROUTINE lgama(kf, x, gl)

!    ==================================================
!    Purpose: Compute gamma function â(x) or ln[â(x)]
!    Input:   x  --- Argument of â(x) ( x > 0 )
!             KF --- Function code
!                    KF=1 for â(x); KF=0 for ln[â(x)]
!    Output:  GL --- â(x) or ln[â(x)]
!    ==================================================

INTEGER, INTENT(IN)     :: kf
REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: gl

REAL (dp), PARAMETER  :: a(10) = (/ 8.333333333333333D-02,  &
   -2.777777777777778D-03, 7.936507936507937D-04, -5.952380952380952D-04,  &
    8.417508417508418D-04, -1.917526917526918D-03, 6.410256410256410D-03,  &
   -2.955065359477124D-02, 1.796443723688307D-01, -1.39243221690590D+00 /)
REAL (dp)  :: gl0, x0, x2, xp
INTEGER    :: k, n

x0 = x
IF (x == 1.0 .OR. x == 2.0) THEN
  gl = 0.0D0
  GO TO 30
ELSE IF (x <= 7.0) THEN
  n = 7-x
  x0 = x + n
END IF
x2 = 1.0D0 / (x0*x0)
xp = 6.283185307179586477_dp
gl0 = a(10)
DO  k = 9, 1, -1
  gl0 = gl0 * x2 + a(k)
END DO
gl = gl0 / x0 + 0.5D0 * LOG(xp) + (x0-.5D0) * LOG(x0) - x0
IF (x <= 7.0) THEN
  DO  k = 1, n
    gl = gl - LOG(x0-1.0D0)
    x0 = x0 - 1.0D0
  END DO
END IF

30 IF (kf == 1) gl = EXP(gl)
RETURN
END SUBROUTINE lgama
 
END MODULE lgama_func
 
 
 
PROGRAM mlgama
USE lgama_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:43

!    ===================================================
!    Purpose: This program computes the gamma function
!             â(x) for x > 0 using subroutine LGAMA
!    Examples:
!               x           â(x)
!             -------------------------
!              0.5     .1772453851D+01
!              2.5     .1329340388D+01
!              5.0     .2400000000D+02
!              7.5     .1871254306D+04
!             10.0     .3628800000D+06
!    ===================================================

REAL (dp)  :: gl, x
INTEGER    :: l

WRITE (*,*) '   x           â(x)'
WRITE (*,*) ' -------------------------'
DO  l = 0, 20, 5
  x = 0.5D0 * l
  IF (l == 0) x = 0.5
  CALL lgama(1, x, gl)
  WRITE (*,5000) x, gl
END DO
WRITE (*,*) 'Please enter x: '
READ (*,*) x
CALL lgama(1, x, gl)
WRITE (*,5000) x, gl
STOP

5000 FORMAT (' ', f5.1, g20.10)
END PROGRAM mlgama
